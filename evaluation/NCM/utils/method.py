import sys
import os
import subprocess
import random
import warnings
import pandas as pd
import numpy as np
from time import time

from scipy.interpolate import interp1d
from scipy.sparse.csr import csr_matrix

from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans, AgglomerativeClustering
import jenkspy

import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection

# optimizer
TRY_USE_NUMBA = True


def jit_if_available(func):
    # default "do nothing" decorator with the numba-like interface
    def decorated(*args, **kwargs):
        return func(*args, **kwargs)
    return decorated


if TRY_USE_NUMBA:
    try:
        from numba import jit  # as jit_if_available
        jit_if_available = jit()
    except ModuleNotFoundError:
        print("Install numba for a bit faster calculations")


def zscore(df):
    m = df.mean(axis=1)
    df = df.T - m
    df = df.T
    s = df.std(axis=1)
    df = df.T/s
    df = df.T
    # set to 0 not variable genes
    zero_var_genes = s[s == 0].index.values
    if len(zero_var_genes) > 0:
        print(len(zero_var_genes),
              "zero variance rows detected, assign zero z-scores ",
              file=sys.stderr)
    df.loc[zero_var_genes, :] = 0
    return df


def validate_input_matrix(exprs, tol=0.01, standradize=True, verbose=True):
    m = exprs.mean(axis=1)
    std = exprs.std(axis=1)
    # find zero variance rows
    zero_var = list(std[std == 0].index.values)
    if len(zero_var) > 0:
        print(len(zero_var), " zero variance rows will be dropped:",
              zero_var, file=sys.stderr)
        exprs = exprs.loc[std > 0]
        m = m[std > 0]
        std = std[std > 0]

    mean_passed = np.all(np.abs(m) < tol)
    std_passed = np.all(np.abs(std-1) < tol)
    if not (mean_passed and std_passed):
        print("Input is not standardized.", file=sys.stderr)
        if standradize:
            exprs = zscore(exprs)
            if not mean_passed:
                print("Centering mean to 0", file=sys.stderr)
            if not std_passed:
                print("Scaling std to 1", file=sys.stderr)
    else:
        if verbose:
            print("Input is standardized.", file=sys.stderr)
    if len(set(exprs.index.values)) < exprs.shape[0]:
        print("Row names are not unique.", file=sys.stderr)
    missing_values = exprs.isna().sum(axis=1)
    if missing_values[missing_values > 0].shape[0] > 0:
        print("Missing values detected in ",
              missing_values[missing_values > 0].index.values,
              file=sys.stderr)
    return exprs


def calc_snr_per_row(s, N, exprs, exprs_sums, exprs_sq_sums):
    bic = exprs[:, :s]
    bic_sums = bic.sum(axis=1)
    bic_sq_sums = np.square(bic).sum(axis=1)

    bg_counts = N - s
    bg_sums = exprs_sums - bic_sums
    bg_sq_sums = exprs_sq_sums - bic_sq_sums

    bic_mean, bic_std = calc_mean_std_by_powers((s, bic_sums, bic_sq_sums))
    bg_mean, bg_std = calc_mean_std_by_powers((bg_counts, bg_sums, bg_sq_sums))

    snr_dist = (bic_mean - bg_mean) / (bic_std + bg_std)

    return snr_dist


def calc_mean_std_by_powers(powers):
    count, val_sum, sum_sq = powers

    mean = val_sum / count  # what if count == 0?
    std = np.sqrt((sum_sq / count) - mean*mean)
    return mean, std


def calc_SNR(ar1, ar2):
    return (np.mean(ar1) - np.mean(ar2)) / (np.std(ar1) + np.std(ar2))


# Binarization #

def generate_null_dist(N, min_n_samples, n_permutations=10000, pval=0.001,
                       seed=42, verbose=True):
    # samples N values from standard normal distribution, and split them into
    # bicluster  background returns bicluster sizes tested, SNR thresholds for
    # each size, the distribution of SNR for each size

    n_permutations = max(n_permutations, int(1.0/pval*10))

    sizes = np.arange(min_n_samples, int(N/2)+1)
    if verbose:
        print("\tGenerate empirical distribuition of " +
              "SNR depending on the bicluster size ...")
        print(f"\t\ttotal samples: {N},\n\t\t" +
              f"number of samples in a bicluster: {sizes[0]} - " +
              f"{sizes[-1]},\n\t\tn_permutations: {n_permutations}")
        print("snr_pval threshold:", pval)

    # generate random expressions from st.normal
    exprs = np.zeros((n_permutations, N))
    np.random.seed(seed=seed)
    for i in range(n_permutations):
        exprs[i, ] = sorted(np.random.normal(size=N))

    exprs_sums = exprs.sum(axis=1)
    exprs_sq_sums = np.square(exprs).sum(axis=1)

    sizes = np.arange(min_n_samples, int(N/2)+1)
    null_distribution = np.zeros((sizes.shape[0], n_permutations))

    for s in sizes:
        null_distribution[s-min_n_samples, :] = -1 * \
            calc_snr_per_row(s, N, exprs, exprs_sums, exprs_sq_sums)

    return null_distribution


def get_trend(sizes, thresholds, plot=True):
    # smoothens the trend and retunrs a function min_SNR(size; p-val. cutoff)
    # rint("\t\t\tfraction:",round(min(0.1,15/len(sizes)),2))
    lowess = sm.nonparametric.lowess
    lowess_curve = lowess(thresholds, sizes, frac=min(
        0.1, 15/len(sizes)), return_sorted=True, is_sorted=False)
    # ,kind="nearest-up",fill_value="extrapolate")
    get_min_snr = interp1d(lowess_curve[:, 0], lowess_curve[:, 1])
    if plot:
        plt.plot(sizes, thresholds, "b--", lw=2)
        plt.plot(sizes, get_min_snr(sizes), "r-", lw=2)
        plt.xlabel("n_samples")
        plt.ylabel("SNR threshold")
        plt.ylim((0, 5))
        plt.show()
    return get_min_snr


def calc_e_pval(snr, size, min_n_samples, null_distribution):
    e_dist = null_distribution[size - min_n_samples, ]
    return (len(
        e_dist[e_dist >= abs(snr)])+1.0)/(null_distribution.shape[1]+1.0)


def jenks_binarization(exprs, null_distribution, min_n_samples, verbose=True,
                       plot=True, plot_SNR_thr=3.0, show_fits=[]):
    t0 = time()
    if verbose:
        print("\tJenks method is chosen ...")
    binarized_expressions = {}
    N = exprs.shape[1]
    n_bins = max(20, int(N/10))
    n_ambiguous = 0
    stats = {}
    for i, (gene, row) in enumerate(exprs.iterrows()):
        # TBD: add special treatment of zero expressions, try zero vs non-zero
        # before Jenks breaks

        values = row.values
        hist_range = values.min(), values.max()
        colors = ["grey", "grey"]

        neg_mask = values == hist_range[0]
        # if many min-values
        if neg_mask.astype(int).sum() >= min_n_samples:
            pos_mask = values > hist_range[0]

        else:
            breaks = jenkspy.jenks_breaks(values, nb_class=2)
            threshold = breaks[1]
            neg_mask = values < threshold
            pos_mask = values >= threshold

        down_group = values[neg_mask]
        up_group = values[pos_mask]

        snr, size, e_pval, direction = np.nan, np.nan, np.nan, None
        n_up = len(up_group)
        n_down = len(down_group)
        size = min(n_up, n_down)
        if size >= min_n_samples:
            # calculate SNR, size and empirical p-val
            snr = calc_SNR(up_group, down_group)
            # SNR p-val depends on biclister size
            e_pval = calc_e_pval(snr, size, min_n_samples, null_distribution)

            # -min_n_samples: # up-regulated group is bicluster
            if n_down-n_up >= 0:
                binarized_expressions[gene] = pos_mask.astype(int)
                colors[1] = 'red'
                direction = "UP"

            # -min_n_samples: # down-regulated group is bicluster
            if n_up-n_down > 0:
                # if not gene in binarized_expressions.keys():
                binarized_expressions[gene] = neg_mask.astype(int)
                # else:
                #    binarized_expressions[gene+"_DOWN"] = neg_mask.astype(int)
                colors[0] = "blue"
                direction = "DOWN"

            # in case of insignificant size difference
            # between up- and down-regulated groups
            # the bigger half is treated as signal too
            if abs(n_up-n_down) <= min_n_samples:
                colors = 'blue', 'red'

        stats[gene] = {"SNR": snr, "size": size,
                       "pval": e_pval, "direction": direction}

        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:", i)

        # plotting selected genes
        if gene in show_fits or abs(snr) > plot_SNR_thr:
            plot_binarized_feature(
                gene, down_group, up_group, colors, hist_range, snr, e_pval)

    stats = pd.DataFrame.from_dict(stats).T
    binarized_expressions = pd.DataFrame.from_dict(binarized_expressions)

    # logging
    if verbose:
        print("\tJenks binarization for " +
              "{} features completed in {:.2f} s".format(
                  len(exprs), time()-t0))

    return binarized_expressions, stats


def plot_binarized_feature(feature_name, down_group, up_group, colors,
                           hist_range, snr, e_pval):
    down_color, up_color = colors
    n_bins = int(max(20, (len(down_group)+len(up_group))/10))
    n_bins = min(n_bins, 200)
    tmp = plt.hist(down_group, bins=n_bins, alpha=0.5,
                   color=down_color, range=hist_range)
    tmp = plt.hist(up_group, bins=n_bins, alpha=0.5,
                   color=up_color, range=hist_range)
    tmp = plt.title("{}: SNR={:.2f}, neg={}, pos={}, p-val={:.2e}".format(
        feature_name, snr, len(down_group), len(up_group), e_pval))
    plt.show()


def select_pos_neg(row, null_distribution, min_n_samples, seed=42,
                   prob_cutoff=0.5, method="GMM"):
    """ find 'higher' (positive), and 'lower' (negative) signal in vals.
        vals are found with GM binarization
    """
    is_converged = None
    if method == "GMM":
        # this is to ignore convergence warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            row2d = row[:, np.newaxis]  # adding mock axis for GM interface
            # np.random.seed(seed=seed)
            model = GaussianMixture(n_components=2, init_params="kmeans",
                                    max_iter=len(row), n_init=1,
                                    covariance_type="spherical",
                                    random_state=seed).fit(row2d)
            is_convergend = model.converged_
            p0 = model.predict_proba(row2d)[:, 0]
            labels = np.zeros(len(row), dtype=bool)

            # let labels == True be always a smaller sample set
            if len(labels[p0 >= prob_cutoff]) >= len(
                    labels[p0 < 1-prob_cutoff]):
                labels[p0 < 1-prob_cutoff] = True
            else:
                labels[p0 >= prob_cutoff] = True

    elif method in ["kmeans", "ward"]:
        row2d = row[:, np.newaxis]  # adding mock axis
        if method == "kmeans":
            model = KMeans(n_clusters=2, max_iter=len(row), n_init=1,
                           random_state=seed)
        elif method == "ward":
            model = AgglomerativeClustering(n_clusters=2, linkage='ward')
        # elif method == "HC_ward":
        #    model = Ward(n_clusters=2)
        labels = np.zeros(len(row), dtype=bool)
        pred_labels = model.fit_predict(row2d)
        # let labels == True be always a smaller sample set
        if len(pred_labels[pred_labels == 1]) >= len(
                pred_labels[pred_labels == 0]):
            labels[pred_labels == 0] = True
        else:
            labels[pred_labels == 1] = True
    else:
        print("wrong method name", method,
              "must be ['GMM','kmeans','ward']", file=sys.stderr)

    # special treatment for cases when bic distribution is too wide and
    # overlaps bg distribution remove from bicluster samples with the sign
    # different from its median sign
    if len(row[labels]) > 0:
        if np.median(row[labels]) >= 0:
            labels[row < 0] = False
        else:
            labels[row > 0] = False

    n0 = len(labels[labels])
    n1 = len(labels[~labels])

    assert n0 + n1 == len(row)

    snr = np.nan
    e_pval = np.nan
    size = np.nan
    mask_pos = np.zeros(len(row), dtype=bool)
    mask_neg = np.zeros(len(row), dtype=bool)

    if n0 >= min_n_samples:
        snr = calc_SNR(row[labels], row[~labels])
        size = n0
        # SNR p-val depends on biclister size
        e_pval = calc_e_pval(snr, size, min_n_samples, null_distribution)
        if snr > 0:
            mask_pos = labels
            mask_neg = ~labels
        else:
            mask_neg = labels
            mask_pos = ~labels

    return mask_pos, mask_neg, abs(snr), size, e_pval, is_converged


def GM_binarization(exprs, null_distribution, min_n_samples, verbose=True,
                    plot=True, plot_SNR_thr=2, show_fits=[],
                    seed=1, prob_cutoff=0.5, method="GMM"):
    t0 = time()

    binarized_expressions = {}

    stats = {}
    for i, (gene, row) in enumerate(exprs.iterrows()):
        row = row.values
        pos_mask, neg_mask, snr, size, e_pval, is_converged = select_pos_neg(
            row, null_distribution, min_n_samples, seed=seed,
            prob_cutoff=prob_cutoff, method=method)

        direction = None

        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:", i)

        hist_range = row.min(), row.max()
        colors = ["grey", "grey"]
        up_group = row[pos_mask]
        down_group = row[neg_mask]
        n_up = len(up_group)
        n_down = len(down_group)

        # -min_n_samples: # up-regulated group is bicluster
        if n_down-n_up >= 0:
            binarized_expressions[gene] = pos_mask.astype(int)
            colors[1] = 'red'
            direction = "UP"

        # -min_n_samples: # down-regulated group is bicluster
        if n_up-n_down > 0:
            binarized_expressions[gene] = neg_mask.astype(int)
            colors[0] = "blue"
            direction = "DOWN"

        # in case of insignificant size difference
        # between up- and down-regulated groups
        # the bigger half is treated as signal too
        if abs(n_up-n_down) <= min_n_samples:
            colors = 'blue', 'red'

        # plotting selected genes
        if gene in show_fits or abs(snr) > plot_SNR_thr:
            plot_binarized_feature(
                gene, down_group, up_group, colors, hist_range, snr, e_pval)

        stats[gene] = {"pval": e_pval, "SNR": snr, "size": size,
                       "direction": direction, "convergence": is_converged}

    stats = pd.DataFrame.from_dict(stats).T
    binarized_expressions = pd.DataFrame.from_dict(binarized_expressions)

    # logging
    if verbose:
        print(
            "\tGMM binarization for {} features completed in {:.2f} s".format(
                len(exprs), time()-t0))

    return binarized_expressions, stats


def binarize(binarized_fname_prefix, exprs=None, method='GMM',
             save=True, load=False,
             min_n_samples=10, pval=0.001,
             plot_all=True, plot_SNR_thr=np.inf, show_fits=[],
             verbose=True, seed=random.randint(0, 100000), prob_cutoff=0.5,
             n_permutations=10000):
    '''
       binarized_fname_prefix is a basename of binarized data file;
       exprs is a dataframe with normalized features to be binarized.
    '''
    t0 = time()

    # a file with binarized gene expressions
    bin_exprs_fname = binarized_fname_prefix + ".seed=" + \
        str(seed)+".bin_method="+method+".binarized.tsv"
    # a file with statistics of binarization results
    bin_stats_fname = binarized_fname_prefix + ".seed=" + \
        str(seed) + ".bin_method=" + method + ".binarization_stats.tsv"
    # a file with background SNR distributions for each biclsuter size
    n_permutations = max(n_permutations, int(1.0/pval*10))
    bin_bg_fname = binarized_fname_prefix + ".seed=" + \
        str(seed)+".n="+str(n_permutations)+".background.tsv"

    load_failed = False
    if load:
        load_failed = False
        try:
            if verbose:
                print("Load binarized features from",
                      bin_exprs_fname, "\n", file=sys.stdout)
            # load binary expressions
            binarized_data = pd.read_csv(
                bin_exprs_fname, sep="\t", index_col=0)
        except ZeroDivisionError:
            print("file "+bin_exprs_fname +
                  " is not found and will be created", file=sys.stderr)
            load_failed = True
        try:
            # load stats
            stats = pd.read_csv(bin_stats_fname, sep="\t", index_col=0)
            if verbose:
                print("Load statistics from",
                      bin_stats_fname, "\n", file=sys.stdout)
        except ZeroDivisionError:
            print("file "+bin_stats_fname +
                  " is not found and will be created", file=sys.stderr)
            load_failed = True
        try:
            # load background distribution
            null_distribution = pd.read_csv(
                bin_bg_fname, sep="\t", index_col=0)
            null_distribution = null_distribution.values
            if verbose:
                print("Load background distribution from",
                      bin_bg_fname, "\n", file=sys.stdout)
        except ZeroDivisionError:
            print("file "+bin_bg_fname +
                  " is not found and will be created", file=sys.stderr)
            load_failed = True

    if not load or load_failed:
        if exprs is None:
            print("Provide either raw or binarized data.", file=sys.stderr)
            return None
        # binarize features
        start_time = time()
        if verbose:
            print("\nBinarization started ....\n")

        t0 = time()

        null_distribution = generate_null_dist(exprs.shape[1], min_n_samples,
                                               pval=pval,
                                               n_permutations=n_permutations,
                                               seed=seed, verbose=verbose)

        if verbose:
            print("\tSNR thresholds for individual features " +
                  "computed in {:.2f} seconds".format(time()-t0))

        if method == "Jenks":
            binarized_data, stats = jenks_binarization(
                exprs, null_distribution, min_n_samples,
                verbose=verbose,
                plot=plot_all,
                plot_SNR_thr=plot_SNR_thr, show_fits=show_fits)
        elif method in ["GMM", "kmeans", "ward"]:
            binarized_data, stats = GM_binarization(
                exprs, null_distribution, min_n_samples,
                plot=plot_all, plot_SNR_thr=plot_SNR_thr,
                prob_cutoff=prob_cutoff,
                show_fits=show_fits, verbose=verbose,
                seed=seed, method=method)
        else:
            print("Method must be 'GMM','kmeans','ward', or 'Jenks'.",
                  file=sys.stderr)
            return

        if verbose:
            print("\tBinarization runtime: {:.2f} s".format(
                time()-start_time), file=sys.stdout)

        if save:
            # save binarized data
            binarized_data.to_csv(bin_exprs_fname, sep="\t")
            if verbose:
                print("Binarized gene expressions are saved to",
                      bin_exprs_fname, file=sys.stdout)

            # save binarization statistics
            stats.to_csv(bin_stats_fname, sep="\t")
            if verbose:
                print("Statistics is saved to",
                      bin_stats_fname, file=sys.stdout)

            # save null distribution: null_distribution, size,threshold
            df = pd.DataFrame(null_distribution,
                              index=range(null_distribution.shape[0]),
                              columns=range(null_distribution.shape[1]))
            df.to_csv(bin_bg_fname, sep="\t")
            if verbose:
                print("Background sitribution is saved to",
                      bin_bg_fname, file=sys.stdout)

    # keep features passed binarization
    # find the threshold
    sizes = np.arange(min_n_samples, int(exprs.shape[1]/2)+1)
    thresholds = np.quantile(null_distribution, q=1-pval, axis=1)
    size_snr_trend = get_trend(sizes, thresholds, plot=False)
    stats["SNR_threshold"] = stats["size"].apply(lambda x: size_snr_trend(x))

    passed = stats.loc[stats["SNR"] > stats["SNR_threshold"], :]
    if verbose:
        print("\t\tUP-regulated features:\t%s" %
              (passed.loc[passed["direction"] == "UP", :].shape[0]),
              file=sys.stdout)
        print("\t\tDOWN-regulated features:\t%s" %
              (passed.loc[passed["direction"] == "DOWN", :].shape[0]),
              file=sys.stdout)

    # keep only binarized features
    binarized_data = binarized_data.loc[:, passed.index.values]
    # add sample names
    binarized_data.index = exprs.columns.values

    if plot_all:
        # plot null distribution of SNR(size) - for shuffled genes
        # n_permutation times for each size
        e_stats = []
        for i in range(null_distribution.shape[0]):
            # limit to 100 points per sample group size
            for j in range(max(100, null_distribution.shape[1])):
                e_stats.append(
                    {"size": min_n_samples+i, "SNR": null_distribution[i, j]})
        e_stats = pd.DataFrame.from_records(e_stats)
        # downsample to 1000 points in total
        if e_stats.shape[0] > 1000:
            e_stats = e_stats.sample(n=1000)
        if verbose:
            print(
                e_stats.shape[0], "points from null distribution plotted",
                file=sys.stderr)
        tmp = plt.figure(figsize=(20, 10))
        tmp = plt.scatter(
            e_stats["size"], e_stats["SNR"], alpha=1, color="grey")

        # plot binarization results for real genes
        passed = stats.loc[stats["SNR"] > stats["SNR_threshold"], :]
        failed = stats.loc[stats["SNR"] <= stats["SNR_threshold"], :]
        tmp = plt.scatter(failed["size"], failed["SNR"],
                          alpha=0.2, color="black")
        tmp = plt.scatter(passed["size"], passed["SNR"],
                          alpha=0.7, color="red")

        # plot cutoff
        sizes = sorted(list(set(stats["size"].values)))
        tmp = plt.plot(sizes, [size_snr_trend(x)
                       for x in sizes], c="yellow", lw=2)
        tmp = plt.xlabel("n_samples")
        tmp = plt.ylabel("SNR threshold")
        tmp = plt.ylim((0, 5))
        tmp = plt.show()

    return binarized_data, stats, null_distribution


# Cluster binarized genes ##

def run_WGCNA(binarized_expressions, fname,
              deepSplit=4, detectCutHeight=0.995,  # see WGCNA documentation
              verbose=False, rscr_path=False):
    t0 = time()

    deepSplit = int(deepSplit)
    if deepSplit not in [0, 1, 2, 3, 4]:
        print("deepSplit must be 1,2,3 or 4. See WGCNA documentation.",
              file=sys.stderr)
        return ([], [])
    if not 0 < detectCutHeight < 1:
        print(
            "detectCutHeight must be between 0-1. See WGCNA documentation.",
            file=sys.stderr)
        return ([], [])
    if verbose:
        print("\tRunning WGCNA for", fname, "...", file=sys.stdout)
    if not rscr_path:
        # assume run_WGCNA.R is in the same folder
        rscr_path = "/".join(os.path.realpath(__file__).split("/")
                             [:-1])+'/run_WGCNA.R'

    # save binarized expression to file
    binarized_expressions.to_csv(fname, sep="\t")

    # run Rscript
    process = subprocess.Popen(['Rscript', rscr_path, fname, str(deepSplit),
                                str(
        detectCutHeight)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    module_file = fname.replace(".tsv", ".modules.tsv")  # stdout.rstrip()
    try:
        modules_df = pd.read_csv(module_file, sep="\t", index_col=0)
    except ZeroDivisionError:
        stderr = stderr.decode('utf-8')
        print("WGCNA errors:", stderr, file=sys.stdout)
        modules_df = pd.DataFrame.from_dict({})
    if verbose:
        print("\tWGCNA runtime: modules detected in {:.2f} s.".format(
            time()-t0), file=sys.stdout)

    # read WGCNA output
    modules = []
    not_clustered = []
    module_dict = modules_df.T.to_dict()
    for i in module_dict.keys():
        genes = module_dict[i]["genes"].strip().split()
        if i == 0:
            not_clustered = genes
        else:
            modules.append(genes)

    # remove WGCNA input and output files
    try:
        os.remove(fname)
        os.remove(module_file)
    except ZeroDivisionError:
        pass

    if verbose:
        print("\tmodules: {}, not clustered features {} ".format(
            len(modules), len(not_clustered)), file=sys.stdout)
        if len(stderr) > 0:
            print(stderr, file=sys.stderr)

    return (modules, not_clustered)


def run_Louvain(similarity, verbose=True, similarity_cutoffs=1/3,
                m=False, plot=False):
    t0 = time()
    if similarity.shape[0] == 0:
        print("no features to cluster", file=sys.stderr)
        return [], [], None

    if verbose:
        print("\tRunning Louvain ...")
    from sknetwork.clustering import Louvain, modularity

    modularities = []
    feature_clusters = {}
    if m:
        # scan the whole range between 1/3 and 3/4 with step 0.01 and stop whem
        # modularity >= m
        similarity_cutoffs = np.arange(1/4, 4/5, 0.01)
    for cutoff in similarity_cutoffs:
        # scan the whole range of cutoffs
        sim_binary = similarity.copy()
        sim_binary[sim_binary < cutoff] = 0
        sim_binary[sim_binary > 0] = 1
        rsums = sim_binary.sum()
        non_zero_features = rsums[rsums > 0].index
        sim_binary = sim_binary.loc[non_zero_features, non_zero_features]
        gene_names = sim_binary.index.values
        sparse_matrix = csr_matrix(sim_binary)
        labels = Louvain(modularity='newman').fit_transform(
            sparse_matrix)  # modularity =  ['potts','dugue', 'newman']
        Q = modularity(sparse_matrix, labels)
        modularities.append(Q)
        feature_clusters[cutoff] = labels

    if len(similarity_cutoffs) == 1:
        # if one cutoff value is chosen
        best_cutoff = similarity_cutoffs[0]
        best_Q = Q
    elif m:
        best_cutoff = 0
        for i in range(len(modularities)):
            if modularities[i] >= m:
                best_cutoff = similarity_cutoffs[i]
                best_Q = modularities[i]
                labels = feature_clusters[best_cutoff]
                break

    else:
        # find the knee point
        from kneed import KneeLocator
        try:
            kn = KneeLocator(similarity_cutoffs, modularities,
                             curve='concave', direction='increasing',
                             online=True)
            best_cutoff = kn.knee
            best_Q = kn.knee_y
            labels = feature_clusters[best_cutoff]
        except ZeroDivisionError:

            print("Failed to identify similarity cutoff", file=sys.stderr)
            print("Cutoff:", best_cutoff, "set to 1/2", file=sys.stdout)
            best_cutoff = 1/2
            print("Modularity:", modularities, file=sys.stdout)
            plt.plot(similarity_cutoffs, modularities, 'bx-')
            plt.xlabel('similarity cutoff')
            plt.ylabel('modularity')
            plt.show()
            return [], [], None

    if plot and len(similarity_cutoffs) > 1:
        plt.plot(similarity_cutoffs, modularities, 'bx-')
        plt.vlines(best_cutoff, plt.ylim()[0], plt.ylim()[
                   1], linestyles='dashed', color="red")
        plt.xlabel('similarity cutoff')
        plt.ylabel('modularity')
        plt.show()

    modules = []
    not_clustered = []

    for label in set(labels):
        ndx = np.argwhere(labels == label).reshape(1, -1)[0]
        genes = gene_names[ndx]
        if len(genes) > 1:
            modules.append(genes)
        else:
            not_clustered.append(genes[0])
    if verbose:
        print("\tLouvain runtime: modules detected in {:.2f} s.".format(
            time()-t0), file=sys.stdout)
        print("\tmodules: {}, not clustered features {} ".format(
            len(modules), len(not_clustered)), file=sys.stdout)
        print("\t\tsimilarity cutoff: {:.2f} modularity: {:.3f}".format(
            best_cutoff, best_Q), file=sys.stdout)
    return modules, not_clustered, best_cutoff


def run_MCL(similarity, inflations=list(np.arange(11, 30)/10)+[3, 3.5, 4, 5],
            verbose=True):
    if similarity.shape[0] == 0:
        return [[], []]
    import markov_clustering as mc
    t0 = time()
    # MCL accepts positive similarity matrix
    gene_names = similarity.index.values
    sparse_matrix = csr_matrix(similarity.values)

    best_Q = -1
    best_result = None
    best_clusters = None
    best_inflation = None
    for inflation in inflations:
        # run MCL with default parameters
        result = mc.run_mcl(sparse_matrix, inflation=inflation)
        clusters = mc.get_clusters(result)    # get clusters
        # Newman modularity
        Q = mc.modularity(matrix=result, clusters=clusters)
        if Q > best_Q:
            best_Q = Q
            best_result = result
            best_clusters = clusters
            best_inflation = inflation
    modules = []
    not_clustered = []
    for cluster in set(best_clusters):
        genes = gene_names[np.array(cluster)]
        if len(genes) > 1:
            modules.append(sorted(genes))
        else:
            not_clustered.append(genes[0])

    if verbose:
        print("\tn_clusters:", len(modules),
              "not clustered:", len(not_clustered))
        print("\t\tinflation:", best_inflation, "modularity:",
              round(best_Q, 3), round(time()-t0), "s.")
    return [modules, not_clustered]


def get_similarity_jaccard_signif(df, qval=0.01, J=0.33):
    from fisher import pvalue
    # based on jaccard similarity
    t0 = time()
    genes = df.columns.values
    n_samples = df.shape[0]
    df = df.values
    results = []
    for i in range(df.shape[1]):
        for j in range(i+1, df.shape[1]):
            g1 = df[:, i]
            g2 = df[:, j]
            o = g1 * g2
            overlap = o.sum()
            u = g1 + g2
            union = u[u > 0].shape[0]
            jaccard = overlap/union
            g1_only = g1.sum() - overlap
            g2_only = g2.sum() - overlap
            p = pvalue(overlap, g1_only, g2_only, n_samples-union)
            results.append({"i": genes[i], "j": genes[j],
                            "J": jaccard, "p_val": p.right_tail,
                            "overlap": overlap, "union": union})
    results = pd.DataFrame.from_records(results)
    # correction form multiple testing
    bh_res, qvals = fdrcorrection(results["p_val"].values, alpha=0.05)
    results["q_val"] = qvals
    results = results.loc[results["q_val"] <= qval, :]
    results = results.loc[results["J"] >= J, :]
    similarity = {}
    for row in results.iterrows():
        g1 = row[1]["i"]
        g2 = row[1]["j"]
        jaccard = row[1]["J"]
        if g1 in similarity.keys():
            similarity[g1][g2] = jaccard
        else:
            similarity[g1] = {g2: jaccard}
        if g2 in similarity.keys():
            similarity[g2][g1] = jaccard
        else:
            similarity[g2] = {g1: jaccard}
    similarity = pd.DataFrame.from_dict(similarity)
    missed_genes = list(set(genes).difference(set(similarity)))
    for g in missed_genes:
        similarity[g] = 0
        similarity = similarity.T
        similarity[g] = 0
        similarity = similarity.T
    # make diagonal 1
    for g in genes:
        similarity.loc[g][g] = 1

    # fillna
    similarity = similarity.fillna(0)
    similarity = similarity.loc[sorted(
        similarity.index), sorted(similarity.index)]
    return similarity


def get_similarity_jaccard(binarized_data, verbose=True):  # ,J=0.5
    t0 = time()
    genes = binarized_data.columns.values
    n_samples = binarized_data.shape[0]
    size_threshold = int(min(0.45*n_samples, (n_samples)/2-10))
    n_genes = binarized_data.shape[1]
    df = np.array(binarized_data.T, dtype=bool)
    results = np.zeros((n_genes, n_genes))
    for i in range(0, n_genes):
        results[i, i] = 1
        g1 = df[i]

        for j in range(i+1, n_genes):
            g2 = df[j]
            o = g1 * g2
            u = g1 + g2
            jaccard = o.sum()/u.sum()
            # try matching complements
            if g1.sum() > size_threshold:
                g1_complement = ~g1
                o = g1_complement * g2
                u = g1_complement + g2
                jaccard_c = o.sum()/u.sum()
            elif g2.sum() > size_threshold:
                g2 = ~g2
                o = g1 * g2
                u = g1 + g2
                jaccard_c = o.sum()/u.sum()
            else:
                jaccard_c = 0
            jaccard = max(jaccard, jaccard_c)
            results[i, j] = jaccard
            results[j, i] = jaccard

    results = pd.DataFrame(data=results, index=genes, columns=genes)
    if verbose:
        print("\tJaccard similarities for " +
              "{} features computed in {:.2f} s.".format(
                  binarized_data.shape[1], time()-t0), file=sys.stdout)
    return results


def get_similarity_corr(df, verbose=True):
    t0 = time()
    corr = df.corr()  # .applymap(abs)
    corr = corr[corr > 0]  # to consider only direct correlations
    corr = corr.fillna(0)
    if verbose:
        print("\tPearson's r similarities for " +
              "{} features computed in {:.2f} s.".format(
                  df.shape[1], time()-t0), file=sys.stdout)
    return corr


def make_TOM(similarity):
    t0 = time()
    df = similarity.values
    fnames = similarity.index.values
    N = similarity.shape[0]
    TOM = np.zeros((N, N))
    for i in range(0, N):
        ai = df[i, ]
        sum_ai = ai.sum()
        for j in range(i, N):
            aj = df[j, ]
            sum_aj = aj.sum()
            aij = aj[i]
            f = (sum_aj+sum_ai)/2-aij
            a_product = np.dot(ai, aj)
            tom_ij = (aij + a_product - aij*aij)/(f+1-aij)
            TOM[i, j] = tom_ij
            TOM[j, i] = tom_ij
    TOM = pd.DataFrame(data=TOM, index=fnames, columns=fnames)
    print("\tTOM computed in {:.2f} s.".format(time()-t0), file=sys.stdout)
    return TOM


# Make biclusters #
def cluster_samples(data, min_n_samples=10, seed=0, method="kmeans"):
    # identify identify bicluster and backgound groups using 2-means
    if method == "kmeans" or method == "Jenks":
        labels = KMeans(n_clusters=2, random_state=seed,
                        n_init=5).fit(data).labels_
    elif method == "ward":
        labels = AgglomerativeClustering(
            n_clusters=2, linkage='ward').fit(data).labels_
    elif method == "GMM":
        labels = GaussianMixture(n_components=2, init_params="kmeans",
                                 max_iter=data.shape[0], n_init=5,
                                 covariance_type="spherical",
                                 random_state=seed).fit_predict(data)
    ndx0 = np.where(labels == 0)[0]
    ndx1 = np.where(labels == 1)[0]
    if min(len(ndx1), len(ndx0)) < min_n_samples:
        return {}
    if len(ndx1) > len(ndx0):
        samples = ndx0
    else:
        samples = ndx1

    bicluster = {"sample_indexes": set(samples), "n_samples": len(samples)}
    return bicluster


def modules2biclusters(modules, data_to_cluster,
                       method="kmeans",
                       min_n_samples=10,
                       min_n_genes=2, seed=0, verbose=True):
    '''Identifies optimal sample set for each module:
    splits samples into two sets in a subspace of each module
    '''
    t0 = time()
    biclusters = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0

    for mid in range(0, len(modules)):
        genes = modules[mid]
        if len(genes) >= min_n_genes:
            # cluster samples in a space of selected genes
            data = data_to_cluster.loc[genes, :].T
            bicluster = cluster_samples(
                data, min_n_samples=min_n_samples, seed=seed, method=method)
            if len(bicluster) > 0:
                bicluster["id"] = i
                bicluster["genes"] = set(genes)
                bicluster["n_genes"] = len(bicluster["genes"])
                biclusters[i] = bicluster
                i += 1

    if verbose:
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %
              (len(modules), round(time()-t0, 2)))
        print("Passed biclusters (>=%s genes, >= samples %s): %s" %
              (min_n_genes, min_n_samples, i-1), file=sys.stdout)
        print()

    return biclusters


def update_bicluster_data(bicluster, data):
    ''' distinguishes up- and down-regulated gene
    adds "samples" and "gene_indexes"
    calculates average z-score
    bicluster must contain "sample_indexes" and "genes"
    data must contain all features, not just binarized'''

    # add "samples" and "gene_indexes"
    sample_names = data.columns.values
    gene_names = data.index.values

    bic_samples = sample_names[list(bicluster["sample_indexes"])]
    bic_genes = list(bicluster["genes"])
    bg_samples = [x for x in sample_names if x not in bic_samples]
    bicluster["samples"] = set(bic_samples)
    bicluster["gene_indexes"] = set(
        [np.where(gene_names == x)[0][0] for x in bicluster["genes"]])

    # distinguish up- and down-regulated features
    m_bic = data.loc[bic_genes, bic_samples].mean(axis=1)
    m_bg = data.loc[bic_genes, bg_samples].mean(axis=1)
    genes_up = m_bic[m_bic >= m_bg].index.values
    genes_down = m_bic[m_bic < m_bg].index.values
    bicluster["genes_up"] = set(genes_up)
    bicluster["genes_down"] = set(genes_down)

    genes_up = m_bic[m_bic >= m_bg].index.values
    genes_down = m_bic[m_bic < m_bg].index.values

    # calculate average z-score for each sample
    if min(len(genes_up), len(genes_down)) > 0:  # take direction into account
        avg_zscore = (data.loc[genes_up, :].sum(
        )-data.loc[genes_down, :].sum())/bicluster["n_genes"]
    else:
        avg_zscore = data.loc[list(bicluster["genes"]), :].mean()

    # compute SNR for average z-score for this bicluster
    m = avg_zscore[bic_samples].mean() - avg_zscore[bg_samples].mean()
    s = avg_zscore[bic_samples].std() + avg_zscore[bg_samples].std()
    snr = np.abs(m)/s
    bicluster["SNR"] = snr
    return bicluster


def merge_biclusters(biclusters, data, J=0.8,
                     min_n_samples=5,
                     seed=42,
                     method="kmeans", verbose=True):
    #  bicluaters -> binary -> jaccard sim
    binary_representation = {}
    N = data.shape[1]
    for i in biclusters.keys():
        b = np.zeros(N)
        s_ndx = list(biclusters[i]["sample_indexes"])
        b[s_ndx] = 1
        binary_representation[i] = b
    binary_representation = pd.DataFrame.from_dict(binary_representation)
    binary_representation.index = data.columns.values
    bic_similarity = get_similarity_jaccard(
        binary_representation, verbose=verbose)
    # find groups of biclusters including the same sample sets
    merged, not_merged, similarity_cutoff = run_Louvain(
        bic_similarity, verbose=False, plot=False, similarity_cutoffs=[J])
    if len(merged) == 0 and verbose:
        print("No biclusters to merge", file=sys.stdout)
        return biclusters

    merged_biclusters = {}
    # add biclusters with no changes
    for bic_id in not_merged:
        merged_biclusters[bic_id] = biclusters[bic_id]

    # merge biclusters with overlapping sample sets
    for bic_group in merged:
        bic_group = sorted(bic_group)
        print("merged biclustres", bic_group, file=sys.stderr)
        new_bicluster = biclusters[bic_group[0]]
        # update genes
        for bic_id in bic_group[1:]:
            bic2 = biclusters[bic_id]
            new_bicluster["genes"] = new_bicluster["genes"] | bic2["genes"]
            new_bicluster["n_genes"] = len(new_bicluster["genes"])
        # update sample set for new bicluster
        # cluster samples in a space of selected genes
        new_bicluster.update(
            cluster_samples(data.loc[list(new_bicluster["genes"]), :].T,
                            min_n_samples=min_n_samples,
                            seed=seed, method=method))
        new_bicluster["n_samples"] = len(new_bicluster['sample_indexes'])
        merged_biclusters[bic_group[0]] = new_bicluster
    return merged_biclusters


def make_biclusters(feature_clusters,
                    binarized_data, data,
                    null_distribution,
                    merge=1.01,
                    min_n_samples=10, min_n_genes=2,
                    method="kmeans",
                    seed=42,
                    cluster_binary=False, verbose=True):
    sample_names = data.columns.values
    gene_names = data.index.values
    biclusters = []

    if cluster_binary:
        data_to_cluster = binarized_data.loc[:, :].T  # binarized expressions
    else:
        # z-scores
        data_to_cluster = data.loc[binarized_data.columns.values, :]

    if len(feature_clusters) == 0:
        print("No biclusters found.", file=sys.stderr)
    else:
        biclusters = modules2biclusters(feature_clusters, data_to_cluster,
                                        method=method,
                                        min_n_samples=min_n_samples,
                                        min_n_genes=2,
                                        verbose=False, seed=seed)

        # merge biclusters with highly similar sample sets
        if merge <= 1.0:
            biclusters = merge_biclusters(biclusters, data, method=method,
                                          J=merge,
                                          min_n_samples=min_n_samples,
                                          seed=seed, verbose=verbose)

        for i in list(biclusters.keys()):
            biclusters[i] = update_bicluster_data(biclusters[i], data)

    biclusters = pd.DataFrame.from_dict(biclusters).T
    # add direction
    biclusters["direction"] = "BOTH"
    biclusters.loc[biclusters["n_genes"] ==
                   biclusters["genes_up"].apply(len), "direction"] = "UP"
    biclusters.loc[biclusters["n_genes"] ==
                   biclusters["genes_down"].apply(len), "direction"] = "DOWN"

    # add p-value for bicluster SNR (computed for avg. zscores)
    # use the same distribution as for single features
    biclusters["e_pval"] = biclusters.apply(
        lambda row: calc_e_pval(row["SNR"],
                                row["n_samples"],
                                min_n_samples,
                                null_distribution), axis=1)

    # sort and reindex
    biclusters = biclusters.sort_values(
        by=["e_pval", "SNR"], ascending=[True, False])
    biclusters.index = range(0, biclusters.shape[0])

    biclusters = biclusters.loc[:, ["SNR", "e_pval", "n_genes", "n_samples",
                                    "genes", "samples", "direction",
                                    "genes_up", "genes_down", "gene_indexes",
                                    "sample_indexes"]]

    return biclusters


# reading and writing ##

def read_bic_table(file_name, parse_metadata=False):
    if not os.path.exists(file_name):
        return pd.DataFrame()
    biclusters = pd.read_csv(file_name, sep="\t", index_col=0, comment="#")
    if len(biclusters) == 0:
        return pd.DataFrame()
    else:
        biclusters.loc[:, ["genes_up", "genes_down"]] = biclusters.loc[:, [
            "genes_up", "genes_down"]].fillna("")
        biclusters["genes"] = biclusters["genes"].apply(
            lambda x: set(x.split(" ")))
        biclusters["genes_up"] = biclusters["genes_up"].apply(
            lambda x: set(x.split(" ")))
        biclusters["genes_down"] = biclusters["genes_down"].apply(
            lambda x: set(x.split(" ")))
        biclusters["samples"] = biclusters["samples"].apply(
            lambda x: set(x.split(" ")))
        biclusters["gene_indexes"] = biclusters["gene_indexes"].apply(
            lambda x: set(map(int, set(x.split(" ")))))
        biclusters["sample_indexes"] = biclusters["sample_indexes"].apply(
            lambda x: set(map(int, set(x.split(" ")))))

    # resulting_bics.set_index("id",inplace=True)
    if parse_metadata:
        f = open(file_name, 'r')
        metadata = f.readline()
        f.close()
        if metadata.startswith("#"):
            metadata = metadata.replace("#", "").rstrip()
            metadata = metadata.split("; ")
            metadata = dict([x.split("=") for x in metadata])
            return biclusters, metadata
    return biclusters


def write_bic_table(bics_dict_or_df, results_file_name, to_str=True,
                    add_metadata=False,
                    seed=None, min_n_samples=None,
                    bin_method=None, clust_method=None, pval=None,
                    alpha=None, beta_K=None, similarity_cutoff=None,
                    ds=None, dch=None, m=None):
    if add_metadata:
        metadata = "#seed="+str(seed)+"; "+"pval="+str(pval) + \
            "; "+"min_n_samples="+str(min_n_samples)+"; "
        metadata = metadata + "b="+bin_method+"; "
        metadata = metadata + "c="+clust_method+"; "
        if clust_method == "Louvain":
            metadata = metadata + "simiarity_cutoff=" + \
                str(similarity_cutoff)+"; modularity="+str(m)
        elif clust_method == "WGCNA":
            metadata = metadata + "ds="+str(ds)+"; dch="+str(dch)
        elif clust_method == "DESMOND":
            metadata = metadata + "alpha=" + \
                str(alpha)+"; " + "beta_K="+str(beta_K)
        else:
            print("Unknown 'clust_method'", clust_method, file=sys.stderr)
        with open(results_file_name, 'w') as f:
            f.write(metadata+"\n")
        write_mode = 'a'
    else:
        write_mode = 'w'

    bics = bics_dict_or_df.copy()
    if len(bics) == 0:
        print("No biclusters found", file=sys.stderr)
    else:
        if not bics.isinstance(pd.DataFrame):
            bics = pd.DataFrame.from_dict(bics)
        if to_str:
            bics["genes"] = bics["genes"].apply(
                lambda x: " ".join(map(str, sorted(x))))
            bics["genes_up"] = bics["genes_up"].apply(
                lambda x: " ".join(map(str, sorted(x))))
            bics["genes_down"] = bics["genes_down"].apply(
                lambda x: " ".join(map(str, sorted(x))))
            bics["samples"] = bics["samples"].apply(
                lambda x: " ".join(map(str, sorted(x))))
            bics["gene_indexes"] = bics["gene_indexes"].apply(
                lambda x: " ".join(map(str, sorted(x))))
            bics["sample_indexes"] = bics["sample_indexes"].apply(
                lambda x: " ".join(map(str, sorted(x))))
        bics.index = range(0, bics.shape[0])
        bics.index.name = "id"
        cols = bics.columns.values
    bics.to_csv(results_file_name, sep="\t", mode=write_mode)
