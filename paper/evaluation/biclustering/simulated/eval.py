import sys
import numpy as np
import pandas as pd
from fisher import pvalue
from statsmodels.stats.multitest import fdrcorrection

from utils.method import zscore

# from utils.eval import find_best_matches
from fisher import pvalue
from statsmodels.stats.multitest import fdrcorrection


def generate_exprs(data_sizes, g_size=5, frac_samples=[0.05, 0.1, 0.25, 0.5], m=2.0, std=1,
                   z=True,
                   outdir="./", outfile_basename="",
                   g_overlap=False, s_overlap=True,
                   seed=42, add_coexpressed=[]):
    n_genes, N = data_sizes
    biclusters = {}
    coexpressed_modules = []

    # generate background model
    np.random.seed(seed)
    exprs = pd.DataFrame(np.random.normal(loc=0, scale=1.0, size=(n_genes, N)))
    exprs.columns = ["s_" + str(x) for x in exprs.columns.values]
    exprs.index = ["g_" + str(x) for x in exprs.index.values]

    # implant bicluster
    N = exprs.shape[1]
    bic_g = []
    bic_s = []
    bg_g = set(exprs.index.values).difference(set(bic_g))
    bg_s = set(exprs.columns.values).difference(set(bic_s))
    bicluster_genes = []
    for s_frac in frac_samples:
        s_size = int(s_frac * N)
        # select random sets of samples and genes from the background
        bic_genes = list(np.random.choice(list(bg_g), size=g_size, replace=False))
        bic_samples = list(np.random.choice(list(bg_s), size=s_size, replace=False))
        bic_g += bic_genes
        bic_s += bic_samples
        # identify samples outside the bicluster
        if not g_overlap:
            bg_g = bg_g.difference(set(bic_g))
        if not s_overlap:
            bg_s = bg_s.difference(set(bic_s))
        # generate bicluster
        biclusters[s_frac] = {"genes": set(bic_genes), "samples": set(bic_samples), "frac": s_frac,
                              "n_genes": len(bic_genes), "n_samples:": len(bic_samples)}
        bic_exprs = np.random.normal(loc=m, scale=std, size=(g_size, s_size))
        # implant biclusters
        exprs.loc[bic_genes, bic_samples] += bic_exprs
        bicluster_genes += bic_genes

    # add modules of co-expressed genes 
    bg_g = set(exprs.index.values).difference(set(bicluster_genes))
    r = 0.5
    for module in add_coexpressed:
        module_genes = list(np.random.choice(list(bg_g), size=module, replace=False))
        n = exprs.loc[module_genes[0], :]
        for i in range(1, module):
            n_i = n * r + np.sqrt(1 - r ** 2) * exprs.loc[module_genes[i], :]
            exprs.loc[module_genes[i], :] = n_i
        print("\tco-exprs. module ", module, "r=",
              (exprs.loc[module_genes, :].T.corr().sum().sum() - module) / (module ** 2 / 2 - module))
        coexpressed_modules.append(module_genes)

    if z:
        # center to 0 and scale std to 1
        exprs = zscore(exprs)
    biclusters = pd.DataFrame.from_dict(biclusters).T
    # biclusters.set_index("frac",inplace = True,drop=True)
    biclusters_ = biclusters.copy()

    if outfile_basename:
        # overlap extension
        if s_overlap:
            if g_overlap:
                overlap_ext = ",overlap=yes"
            else:
                overlap_ext = ",overlap=s"
        elif g_overlap:
            overlap_ext = ",overlap=g"
        else:
            overlap_ext = ",overlap=no"
        # save expressions 
        exprs_file = outdir + "/" + outfile_basename + ".n_genes=" + str(g_size) + ",m=" + str(m) + ",std=" + str(std)
        exprs_file += overlap_ext + ".exprs_z.tsv"
        print("expressions:", exprs_file)
        exprs.to_csv(exprs_file, sep="\t")

        # save ground truth 
        biclusters["n_genes"] = biclusters["genes"].apply(lambda x: len(x))
        biclusters["n_samples"] = biclusters["samples"].apply(lambda x: len(x))
        biclusters["genes"] = biclusters["genes"].apply(lambda x: " ".join((map(str, sorted(x)))))
        biclusters["samples"] = biclusters["samples"].apply(lambda x: " ".join((map(str, sorted(x)))))

        biclusters_file = outdir + "/" + outfile_basename + ".n_genes=" + str(g_size) + ",m=" + str(m) + ",std=" + str(
            std)
        biclusters_file += overlap_ext + ".biclusters.tsv"
        biclusters.to_csv(biclusters_file, sep="\t")
        print("true bilusters:", biclusters_file)
        biclusters.to_csv(biclusters_file, sep="\t")

    return exprs, biclusters_, coexpressed_modules


def make_known_groups(annot, exprs, target_col="genefu_z", verbose=False):
    samples = set(exprs.columns.values).intersection(set(annot.index.values))
    if verbose:
        print("Total samples:", len(samples), file=sys.stdout)
    annot = annot.loc[list(samples), :]
    groups = set(annot.loc[:, target_col].values)

    known_groups = {}
    for group in groups:
        if group == group:
            group_samples = set(annot.loc[annot[target_col] == group, :].index.values)
            group_samples = group_samples.intersection(samples)
            if len(group_samples) > int(len(samples) / 2):
                print("take complement of ", group, file=sys.stderr)
                group_samples = samples.difference(group_samples)
            known_groups[group] = group_samples  # {"set":group_samples,"complement": samples.difference(group_samples)}
            if verbose:
                print(group, round(len(group_samples) / len(samples), 2),
                      len(group_samples), len(samples.difference(group_samples)))
    return known_groups


def apply_fdr(df_pval):
    df_fdr = {}
    for group in df_pval.columns.values:
        bh_res, adj_pval = fdrcorrection(df_pval[group].values, alpha=0.05)
        df_fdr[group] = adj_pval
    df_fdr = pd.DataFrame.from_dict(df_fdr)
    df_fdr.index = df_pval.index
    # df_fdr["associated"] = df_fdr.apply(lambda row: row[row<0.05].index.values,axis=1)
    return df_fdr


def evaluate_overlaps(biclusters, known_groups, all_elements, dimension="samples"):
    # compute exact Fisher's p-values and Jaccard overlaps for samples
    pvals = {}
    is_enriched = {}
    jaccards = {}
    N = len(all_elements)
    # sanity check - biclusters
    for i in biclusters.index.values:
        bic_members = biclusters.loc[i, dimension]
        if not bic_members.intersection(all_elements) == bic_members:
            print("bicluster {} elements {} are not in 'all_elements'".format(i, " ".join(
                bic_members.difference(all_elements))), file=sys.stderr)
            return
            # sanity check and sorting
    group_names = list(known_groups.keys())
    sorted_group_names = [group_names[0]]  # group names ordered by group size
    for group in group_names:
        group_members = known_groups[group]
        if not group_members.intersection(all_elements) == group_members:
            print(group, "elements are not in 'all_elements'", file=sys.stderr)
            return

        if group != group_names[0]:
            for i in range(len(sorted_group_names)):
                if len(group_members) < len(known_groups[sorted_group_names[i]]):
                    sorted_group_names = sorted_group_names[:i] + [group] + sorted_group_names[i:]
                    break
                elif i == len(sorted_group_names) - 1:
                    sorted_group_names = [group] + sorted_group_names
    # print(sorted_group_names)
    for group in sorted_group_names:
        group_members = known_groups[group]
        pvals[group] = {}
        is_enriched[group] = {}
        jaccards[group] = {}
        for i in biclusters.index.values:
            bic = biclusters.loc[i, :]
            bic_members = bic[dimension]
            shared = len(bic_members.intersection(group_members))
            bic_only = len(bic_members.difference(group_members))
            group_only = len(group_members.difference(bic_members))
            union = shared + bic_only + group_only
            pval = pvalue(shared, bic_only, group_only, N - union)
            if pval.right_tail < pval.left_tail:
                pvals[group][i] = pval.right_tail
                is_enriched[group][i] = True
                jaccards[group][i] = shared / union
            else:
                # save left-tail p-value and record that this is not enrichment
                pvals[group][i] = pval.left_tail
                is_enriched[group][i] = False
                # take complement for the biggest group
                if len(bic_members) > len(group_members):

                    # compute J for bicluster and group complement, (e.g. not_LumA instead of LumA)
                    shared_complement = len(group_members) - shared
                    union_complement = N - union + group_only
                else:
                    shared_complement = len(bic_members) - shared
                    union_complement = N - union + bic_only
                jaccards[group][i] = 0 if union_complement == 0 else shared_complement / union_complement

        # print(group,jaccards[group])

    pvals = pd.DataFrame.from_dict(pvals).loc[:, sorted_group_names]
    is_enriched = pd.DataFrame.from_dict(is_enriched).loc[:, sorted_group_names]
    jaccards = pd.DataFrame.from_dict(jaccards).loc[:, sorted_group_names]
    return pvals, is_enriched, jaccards


def find_best_matches(biclusters, known_groups, all_elements, FDR=0.05,
                      min_SNR=False, min_n_genes=False,
                      dimension="samples", verbose=False,
                      match_unique=True):
    # for each known group starting from the largest one,
    # identifies the best matching bicluster - a significant match with max.Jaccard
    # matched biclusters are removed from comparizon
    # if a bicluster is significantly under-represented in a known group, 
    # compares bicluster to group complement (e.g. not LumA instead of LumA) and sets is_enriched = False
    # returns all 
    results = {}
    if min_SNR:
        biclusters = biclusters.loc[biclusters["avgSNR"] >= min_SNR, :]
    if min_n_genes:
        biclusters = biclusters.loc[biclusters["n_genes"] >= min_n_genes, :]
    if biclusters is None or biclusters.shape[0] == 0:
        for group in sorted(list(known_groups.keys())):
            results[group] = {"J": 0, "J_weighted": 0}
        return results

    # calculate overlap p-vals
    try:
        df_pval, is_enriched, df_jaccard = evaluate_overlaps(biclusters, known_groups,
                                                             all_elements, dimension=dimension)
    except:
        # print("failed to calculate overlap p-values",file=sys.stderr)
        out = evaluate_overlaps(biclusters, known_groups, all_elements, dimension=dimension)
        return

    # BH-adjust for multiple testing
    df_fdr = apply_fdr(df_pval)

    not_matched_biclusters = df_fdr.index.values
    for group in df_fdr.columns.values:
        # choose biclusters with significant overlaps and not overlapping a bigger set
        passed_biclusters = set(df_fdr.loc[df_fdr[group] <= FDR, :].index.values).intersection(not_matched_biclusters)
        significant_matches_j = df_jaccard.loc[list(passed_biclusters), group]

        group_size = len(known_groups[group])
        if verbose:
            print("\t", group, "significant matches:", significant_matches_j.shape[0], file=sys.stdout)
        if significant_matches_j.shape[0] > 0:
            significant_matches_j = significant_matches_j.sort_values().tail(1)

            bm_id = significant_matches_j.index[0]

            # exclude best match
            if match_unique:
                not_matched_biclusters = [x for x in not_matched_biclusters if x != bm_id]

            bm = biclusters.loc[bm_id, :]
            j = df_jaccard.loc[bm_id, group]
            results[group] = {"group_size": group_size, "J": j,  # "J_weighted": j*len(known_groups[group])/N,
                              "is_enriched": is_enriched.loc[bm_id, group],
                              "best_match_id": bm_id}
            results[group].update(bm.to_dict())
        else:
            results[group] = {"group_size": group_size, "J": 0}
    results = pd.DataFrame.from_dict(results).T
    results.index.name = "known_group"
    total_bicluster_members = results["group_size"].sum()
    results["J_weighted"] = results["J"] * results["group_size"] / total_bicluster_members
    return results


def find_best_matching_biclusters(bics1, bics2, N, by="genes", adj_pval_thr=0.05):
    # takes two biluster dafaframes from read_bic_table
    # by = "genes" or "samples"
    # N - total number of objects (genes or samples)
    # finds best matches of bics1 biclusters among bics2 biclusters
    n_bics1 = bics1.shape[0]
    n_bics2 = bics2.shape[0]
    best_matches = {}  # OrderedDict({})
    for row1 in bics1.iterrows():
        bic1 = row1[1]
        i1 = row1[0]
        set1 = bic1[by]
        s1 = len(set1)

        best_matches[i1] = {}
        bm_J = 0
        bm_o = 0
        bm_adj_pval = 1
        bm_id = -1

        for row2 in bics2.iterrows():
            bic2 = row2[1]
            i2 = row2[0]
            set2 = bic2[by]
            o = len(set1.intersection(set2))
            J = 0
            adj_pval = 1
            if o > 0:
                s2 = len(set2)
                s2_ = s2 - o
                s1_ = s1 - o
                u = s1_ + s2_ + o
                bg = N - u
                if by == "genes":
                    pval = pvalue(o, s1_, s2_, bg).right_tail
                    J = o * 1.0 / u
                elif by == "samples":
                    # two-tailed pvalue is computed for samples
                    # because bicluster and background sets can be flipped 
                    pval = pvalue(o, s1_, s2_, bg)
                    if pval.right_tail <= pval.left_tail:
                        pval = pval.right_tail
                    elif max(s1,
                             s2) > 0.4 * N:  # try flipping the largest bicluster if it is close to 50% of the cohort
                        pval = pval.left_tail
                        if s1 > s2:  # flip s1
                            u = bg + s2
                            o = s2_
                        else:
                            u = bg + s1
                            o = s1_
                    J = o * 1.0 / u

                adj_pval = pval * n_bics2 * n_bics1 / 2
                if adj_pval < adj_pval_thr:

                    if J > bm_J or (J == bm_J and adj_pval < bm_adj_pval):
                        bm_J = J
                        bm_adj_pval = adj_pval
                        bm_id = i2
                        bm_o = o

        best_matches[i1][by] = set1
        best_matches[i1]["n_" + by] = s1
        if bm_id != -1:
            best_matches[i1]["bm_" + by] = bics2.loc[bm_id, by]
            best_matches[i1]["bm_n_" + by] = bics2.loc[bm_id, "n_" + by]
            best_matches[i1]["n_shared"] = bm_o
            best_matches[i1]["J"] = bm_J
            best_matches[i1]["adj_pval"] = bm_adj_pval
            best_matches[i1]["shared_" + by] = bics2.loc[bm_id, by].intersection(set1)
            # best_matches[i1]["bic_n_"+by] = bics2.loc[bm_id,"n_"+by]
            best_matches[i1]["bm_id"] = bm_id
    best_matches = pd.DataFrame.from_dict(best_matches).T
    return best_matches
