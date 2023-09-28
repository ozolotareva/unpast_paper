import sys
import os
import subprocess
import random
import warnings
import pandas as pd
import numpy as np
from time import time
import math

from scipy.interpolate import interp1d
from scipy.sparse.csr import csr_matrix
from scipy.stats import chi2_contingency

from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans,AgglomerativeClustering 
from fisher import pvalue

import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection

# optimizer
TRY_USE_NUMBA=True
def jit_if_available(func):
    # default "do nothing" decorator with the numba-like interface
    def decorated(*args, **kwargs):
        return func(*args, **kwargs)
    return decorated
if TRY_USE_NUMBA:
    try:
        from numba import jit  # as jit_if_available
        jit_if_available = jit()
    except:
        print("Numba is not available. Install numba for a bit faster calculations")

def zscore(df):
    m = df.mean(axis=1)
    df = df.T - m
    df = df.T
    s = df.std(axis=1)
    df = df.T/s
    df = df.T
    # set to 0 not variable genes 
    zero_var_genes = s[s==0].index.values
    if len(zero_var_genes)>0:
        print(len(zero_var_genes),"zero variance rows detected, assign zero z-scores ",file = sys.stderr)
    df.loc[zero_var_genes,:] = 0
    return df
        
def validate_input_matrix(exprs, tol=0.01,standradize=True,verbose = True):
    exprs.index = [str(x) for x in exprs.index.values]
    exprs.columns = [str(x) for x in exprs.columns.values]
    m = exprs.mean(axis=1)
    std = exprs.std(axis=1)
    # find zero variance rows
    zero_var = list(std[std==0].index.values)
    if len(zero_var)>0:
        print(len(zero_var)," zero variance rows will be dropped:",zero_var,file = sys.stderr)
        exprs = exprs.loc[std>0]
        m = m[std>0]
        std = std[std>0]
    
    mean_passed = np.all(np.abs(m)<tol)
    std_passed = np.all(np.abs(std-1)<tol)
    if not (mean_passed and std_passed):
        if verbose:
            print("Input is not standardized.",file = sys.stderr)
        if standradize:
            exprs = zscore(exprs)
            if not mean_passed:
                if verbose:
                    print("Centering mean to 0",file= sys.stderr)
            if not std_passed:
                if verbose:
                    print("Scaling std to 1",file= sys.stderr)
    else:
        if verbose:
            print("Input is standardized.",file = sys.stderr)
    if len(set(exprs.index.values)) < exprs.shape[0]:
        print("Row names are not unique.",file = sys.stderr)
    missing_values = exprs.isna().sum(axis=1)
    n_na = missing_values[missing_values>0].shape[0]
    if n_na>0:
        print("Missing values detected in ",
              missing_values[missing_values>0].index.values,file = sys.stderr)
        print("Features dropped:",n_na,file = sys.stderr)
        exprs = exprs.dropna()
    return exprs


def calc_snr_per_row(s, N, exprs, exprs_sums,exprs_sq_sums):
    bic = exprs[:,:s]
    bic_sums = bic.sum(axis=1)
    bic_sq_sums = np.square(bic).sum(axis=1)

    bg_counts = N - s
    bg_sums = exprs_sums - bic_sums
    bg_sq_sums = exprs_sq_sums - bic_sq_sums

    bic_mean, bic_std = calc_mean_std_by_powers((s,bic_sums,bic_sq_sums))
    bg_mean, bg_std = calc_mean_std_by_powers((bg_counts,bg_sums,bg_sq_sums))

    snr_dist = (bic_mean - bg_mean)/ (bic_std + bg_std)
    
    return snr_dist

def calc_mean_std_by_powers(powers):
    count, val_sum, sum_sq = powers

    mean = val_sum / count  # what if count == 0?
    std = np.sqrt((sum_sq / count) - mean*mean)
    return mean, std

def calc_SNR(ar1, ar2):
    return (np.mean(ar1) - np.mean(ar2)) / (np.std(ar1) + np.std(ar2))



######### Binarization #########
def generate_null_dist(N, sizes, n_permutations = 10000,pval = 0.001,
                       seed = 42,verbose = True):
    # samples 'N' values from standard normal distribution, and split them into bicluster and background groups
    # 'sizes' defines bicluster sizes to test
    # returns a dataframe with the distribution of SNR for each bicluster size (sizes x n_permutations )
    t0 = time()
    
    if verbose:
        print("\tGenerate background distribuition of SNR depending on the bicluster size ...",file = sys.stdout)
        print("\t\ttotal samples: %s,\n\t\tnumber of samples in a bicluster: %s - %s,\n\t\tn_permutations: %s"%(N,min(sizes),max(sizes),n_permutations),file = sys.stdout)
        print("\t\tsnr pval threshold:",pval,file = sys.stdout)
    
    exprs =  np.zeros((n_permutations,N)) # generate random expressions from st.normal
    #values = exprs.values.reshape(-1) # random samples from expression matrix
    #exprs = np.random.choice(values,size=exprs.shape[1])
    np.random.seed(seed=seed)
    for i in range(n_permutations):
        exprs[i,] = sorted(np.random.normal(size=N))
    
    exprs_sums = exprs.sum(axis=1)
    exprs_sq_sums = np.square(exprs).sum(axis=1)

    null_distribution = pd.DataFrame(np.zeros((sizes.shape[0],n_permutations)),index=sizes,columns = range(n_permutations))
    
    for s in sizes:
        null_distribution.loc[s,:] = -1*calc_snr_per_row(s, N, exprs, exprs_sums,exprs_sq_sums)
    
    if verbose:
        print("\tBackground ditribution generated in {:.2f} s".format(time()-t0),file = sys.stdout)
    return null_distribution


def get_trend(sizes, thresholds, plot= True, verbose = True):
    # smoothens the trend and retunrs a function min_SNR(size; p-val. cutoff)
    lowess = sm.nonparametric.lowess
    frac = max(5,min(math.floor(int(0.1*len(sizes))),15))/len(sizes)
    if verbose:
        print("\t\t\tLOWESS frac=",round(frac,2), file = sys.stdout)
    lowess_curve = lowess(thresholds,sizes,frac= frac,return_sorted=True,is_sorted=False)
    get_min_snr = interp1d(lowess_curve[:,0],lowess_curve[:,1])#,kind="nearest-up",fill_value="extrapolate")
    if plot:
        plt.plot(sizes, thresholds,"b--",lw=2)
        plt.plot(sizes,get_min_snr(sizes),"r-",lw=2)
        plt.xlabel("n_samples")
        plt.ylabel("SNR")
        plt.ylim((0,5))
        plt.show()
    return get_min_snr

def calc_e_pval(snr,size,null_distribution):
    e_dist = null_distribution.loc[int(size),:]
    return  (len(e_dist[e_dist>=abs(snr)])+1.0)/(null_distribution.shape[1]+1.0)

def jenks_binarization(exprs, min_n_samples,verbose = True,
                      plot=True, plot_SNR_thr= 3.0, show_fits = []):
    import jenkspy
    e_pval = -1
    t0= time()
    if verbose:
        print("\tJenks method is chosen ...")
    binarized_expressions = {}
    N = exprs.shape[1]
    n_bins = max(20,int(N/10))
    n_ambiguous = 0
    stats = {}
    for i, (gene, row) in enumerate(exprs.iterrows()):
        ### TBD: add special treatment of zero expressions, try zero vs non-zero before Jenks breaks
        
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
            neg_mask = values<threshold
            pos_mask = values>=threshold
        
        down_group = values[neg_mask]
        up_group = values[pos_mask]
        
        snr, size, e_pval,direction = np.nan,np.nan,np.nan, None
        n_up = len(up_group)
        n_down = len(down_group)
        size = min(n_up,n_down) 
        if size >= min_n_samples:
            # calculate SNR, size and empirical p-val
            snr = calc_SNR(up_group,down_group)
            
            if n_down-n_up >= 0: #-min_n_samples: # up-regulated group is bicluster
                binarized_expressions[gene] = pos_mask.astype(int)
                colors[1]='red'
                direction ="UP"

            if n_up-n_down > 0: #-min_n_samples: # down-regulated group is bicluster
                #if not gene in binarized_expressions.keys():
                binarized_expressions[gene] = neg_mask.astype(int)
                #else:
                #    binarized_expressions[gene+"_DOWN"] = neg_mask.astype(int)
                colors[0]="blue"
                direction ="DOWN"
                
            # in case of insignificant size difference 
            # between up- and down-regulated groups
            # the bigger half is treated as signal too
            if abs(n_up-n_down) <= min_n_samples:
                colors = 'blue', 'red'
                
                
                
        stats[gene] = {"SNR":snr,"size": size,"direction":direction}
   
        
        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tfeatures processed:",i)
        
        # plotting selected genes
        if gene in show_fits or abs(snr) > plot_SNR_thr:
            plot_binarized_feature(gene,down_group,up_group,colors,hist_range,snr,e_pval)
    
    stats = pd.DataFrame.from_dict(stats).T
    binarized_expressions = pd.DataFrame.from_dict(binarized_expressions)
    
    # logging
    if verbose:
        print("\tJenks binarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))
    
    return binarized_expressions, stats


def  plot_binarized_feature(feature_name, down_group,up_group,colors,hist_range,snr):
    down_color,up_color = colors
    n_bins = int(max(20,(len(down_group)+len(up_group))/10))
    n_bins = min(n_bins,200)
    fig, ax = plt.subplots()
    tmp = ax.hist(down_group, bins=n_bins, alpha=0.5, color=down_color,range=hist_range)
    tmp = ax.hist(up_group, bins=n_bins, alpha=0.5, color=up_color,range=hist_range)
    #tmp = plt.title("{}:    SNR={:.2f},    neg={}, pos={}".format(feature_name,snr,len(down_group),len(up_group)))
    n_samples = min(len(down_group),len(up_group))
    #tmp = ax.set_title("SNR={:.2f},   n_samples={}".format(snr,n_samples))
    ax.text(0.05,0.95,feature_name, ha='left', va='top',transform=ax.transAxes,  fontsize=24)
    ax.text(0.95,0.95,"SNR="+str(round(snr,2))+"\nn_samples="+str(n_samples), ha='right', va='top',transform=ax.transAxes,  fontsize=14)
    #tmp = plt.savefig("figs_binarization/"+feature_name+".hist.svg", bbox_inches='tight', transparent=True)
    plt.show()

def select_pos_neg(row, min_n_samples, seed=42, prob_cutoff=0.5, method = "GMM"):
    """ find 'higher' (positive), and 'lower' (negative) signal in vals. 
        vals are found with GM binarization
    """
    is_converged = None
    if method == "GMM":
        with warnings.catch_warnings(): # this is to ignore convergence warnings 
            warnings.simplefilter('ignore')
            row2d = row[:, np.newaxis]  # adding mock axis for GM interface
            #np.random.seed(seed=seed)
            model = GaussianMixture(n_components=2, init_params="kmeans",
                    max_iter=len(row), n_init = 1, 
                    covariance_type = "spherical",
                    random_state = seed).fit(row2d)
            is_convergend = model.converged_
            p0 = model.predict_proba(row2d)[:,0]
            labels = np.zeros(len(row),dtype=bool)

            # let labels == True be always a smaller sample set
            if len(labels[p0>=prob_cutoff]) >= len(labels[p0<1-prob_cutoff]): 
                labels[p0<1-prob_cutoff] = True
            else: 
                labels[p0>=prob_cutoff] = True
        
    elif method in ["kmeans","ward"]:
        row2d = row[:, np.newaxis]  # adding mock axis 
        if method == "kmeans":
            model = KMeans(n_clusters=2, max_iter=len(row), n_init = 1,
                        random_state = seed)
        elif method == "ward":
            model = AgglomerativeClustering(n_clusters=2, linkage='ward')
        #elif method == "HC_ward":
        #    model = Ward(n_clusters=2)
        labels = np.zeros(len(row),dtype=bool)
        pred_labels = model.fit_predict(row2d)
        # let labels == True be always a smaller sample set
        if len(pred_labels[pred_labels==1]) >= len(pred_labels[pred_labels==0]): 
            labels[pred_labels==0] = True
        else: 
            labels[pred_labels==1] = True
    else:
        print("wrong method name", method, "must be ['GMM','kmeans','ward']",file = sys.stderr )
    
    # special treatment for cases when bic distribution is too wide and overlaps bg distribution
    # remove from bicluster samples with the sign different from its median sign 
    if len(row[labels])>0:
        if np.median(row[labels])>=0:
            labels[row<0]=False
        else:
            labels[row>0]=False
    
    n0 = len(labels[labels])
    n1 = len(labels[~labels])
    
    assert n0 + n1 == len(row)

    snr = np.nan
    e_pval = np.nan
    size = np.nan
    mask_pos = np.zeros(len(row),dtype=bool)
    mask_neg = np.zeros(len(row),dtype=bool)
    
    if n0 >= min_n_samples:
        snr = calc_SNR(row[labels], row[~labels])
        size = n0
        
        if snr>0:
            mask_pos = labels
            mask_neg = ~labels
        else:
            mask_neg = labels
            mask_pos = ~labels

    return mask_pos, mask_neg, abs(snr), size, is_converged


def GM_binarization(exprs, min_n_samples, verbose=True, 
                    plot=True, plot_SNR_thr=2, show_fits=[],
                    seed=1,prob_cutoff=0.5, method = "GMM"):
    t0 = time()
    
    binarized_expressions = {}
    
    stats = {}
    for i, (gene, row) in enumerate(exprs.iterrows()):
        e_pval = -1
        row = row.values
        pos_mask, neg_mask, snr, size, is_converged = select_pos_neg(row, min_n_samples, seed=seed,prob_cutoff=prob_cutoff, method = method)
        
        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:",i)
        
        
        up_group = row[pos_mask]
        down_group = row[neg_mask]
        n_up = len(up_group)
        n_down = len(down_group)
        
        # if smaller sample group shows over- or under-expression
        if n_up <= n_down: # up-regulated group is bicluster
            binarized_expressions[gene] = pos_mask.astype(int)
            direction ="UP"
        else: # down-regulated group is bicluster
            binarized_expressions[gene] = neg_mask.astype(int)
            direction ="DOWN"
        
        stats[gene] = {"pval":0,"SNR":snr,"size":size,"direction":direction,"convergence":is_converged}
        
        if gene in show_fits or (abs(snr) > plot_SNR_thr and plot):
            hist_range = row.min(), row.max()
            
            # set colors to two sample groups
            # red - overexpression
            # blue - under-expression
            # grey - background (group size > 1/2 of all samples)
            colors = ["grey", "grey"]
            
            if n_down-n_up >= 0: # up-regulated group is bicluster
                colors[1]='red'

            if n_up-n_down > 0: # down-regulated group is bicluster
                colors[0]="blue"

            # in case of insignificant size difference 
            # between up- and down-regulated groups
            # the bigger half is treated as signal too
            if abs(n_up-n_down) <= min_n_samples:
                colors = 'blue', 'red'
                
            # plotting
            plot_binarized_feature(gene,down_group,up_group,colors,hist_range,snr)
        
    stats = pd.DataFrame.from_dict(stats).T
    
    binarized_expressions = pd.DataFrame.from_dict(binarized_expressions)

    # logging
    if verbose:
        print("\tBinarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))

    return binarized_expressions, stats


def binarize(binarized_fname_prefix, exprs=None, method='GMM',
             save = True, load = False,
             min_n_samples = 5, pval = 0.001,
             plot_all = True, plot_SNR_thr= np.inf,show_fits = [],
             verbose= True,seed=random.randint(0,100000),prob_cutoff=0.5,
             n_permutations=10000):
    '''
       binarized_fname_prefix is a basename of binarized data file;
       exprs is a dataframe with normalized features to be binarized.
    '''
    t0 = time()
    
    # a file with binarized gene expressions
    bin_exprs_fname = binarized_fname_prefix +".seed="+str(seed)+".bin_method="+method+".min_ns="+str(min_n_samples)+".binarized.tsv"
    # a file with statistics of binarization results
    bin_stats_fname = binarized_fname_prefix +".seed="+str(seed) + ".bin_method="+ method  + ".min_ns="+str(min_n_samples) + ".binarization_stats.tsv"
    # a file with background SNR distributions for each biclsuter size
    n_permutations = max(n_permutations,int(1.0/pval*10))
    bin_bg_fname = binarized_fname_prefix +".seed="+str(seed)+".n="+str(n_permutations)+".min_ns="+str(min_n_samples)+".background.tsv"
    
    if load:
        load_failed = False
        try:
            if verbose:
                print("Load binarized features from",bin_exprs_fname,"\n",file = sys.stdout)
            # load binary expressions
            binarized_data = pd.read_csv(bin_exprs_fname, sep ="\t",index_col=0)
        except:
            print("file "+bin_exprs_fname+" is not found and will be created",file = sys.stderr)
            load_failed = True
        try:
            # load stats 
            stats = pd.read_csv(bin_stats_fname, sep ="\t",index_col=0)
            if verbose:
                print("Load statistics from",bin_stats_fname,"\n",file = sys.stdout)
        except:
            print("file "+bin_stats_fname+" is not found and will be created",file = sys.stderr)
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

        if method=="Jenks":
            binarized_data, stats = jenks_binarization(exprs, min_n_samples,
                                                       verbose = verbose, plot=plot_all,
                                                       plot_SNR_thr=plot_SNR_thr, show_fits = show_fits)
        elif method in ["GMM","kmeans","ward"]:
            binarized_data, stats = GM_binarization(exprs, min_n_samples,
                                                    plot=plot_all, plot_SNR_thr= plot_SNR_thr,
                                                    prob_cutoff=prob_cutoff, 
                                                    show_fits = show_fits,verbose = verbose,
                                                    seed = seed, method = method)
        else:
            print("Method must be 'GMM','kmeans','ward', or 'Jenks'.",file=sys.stderr)
            return

    # load or generate empirical distributions for all bicluster sizes
    N = exprs.shape[1]
    # sizes of binarized features
    sizes1 = set([x for x in stats["size"].values if not np.isnan(x)])
    # no more than 100 of bicluster sizes are computed 
    step = max(int((N-min_n_samples)/100),1)
    sizes2 = set(map(int,np.arange(min_n_samples,int(N/2),step)))
    sizes = np.array(sorted(sizes1|sizes2))
    #sizes = np.arange(min_n_samples,int(exprs.shape[1]/2)+1)
    
    load_failed = False
    if load:
        try:
            #load background distribution
            null_distribution = pd.read_csv(bin_bg_fname, sep ="\t",index_col=0)
            null_distribution.columns = [int(x) for x in null_distribution.columns.values]
            if verbose:
                print("Loaded background distribution from",bin_bg_fname,"\n",file = sys.stdout)
            # check if any new sizes need to be precomputed
            precomputed_sizes = null_distribution.index.values
            add_sizes = np.array(sorted(set(sizes).difference(set(precomputed_sizes))))
            if len(add_sizes)>0:
                null_distribution2 = generate_null_dist(N, add_sizes,
                                               pval = pval,n_permutations=n_permutations,
                                               seed =seed,verbose = verbose)
                null_distribution2.columns = [int(x) for x in null_distribution2.columns.values]
                null_distribution = pd.concat([null_distribution,null_distribution2.loc[:,null_distribution.columns.values]],axis=0)
                if save:
                    null_distribution.loc[sorted(null_distribution.index.values),:].to_csv(bin_bg_fname, sep ="\t")
                    if verbose:
                        print("Background ditribution in %s is updated"%bin_bg_fname,file = sys.stdout)
                null_distribution = null_distribution.loc[sizes,:]
        except:
            print("file "+bin_bg_fname+" is not found and will be created",file = sys.stderr)
            load_failed = True
    if not load or load_failed:
        null_distribution = generate_null_dist(N, sizes,
                                               pval = pval,n_permutations=n_permutations,
                                               seed =seed,verbose = verbose)

    #if not load or load_failed:
    # add SNR p-val depends on bicluster size 
    stats = stats.dropna(subset=["size"])
    stats["pval"] = stats.apply(lambda row: calc_e_pval(row["SNR"],row["size"],null_distribution),axis=1)
    accepted, pval_adj = fdrcorrection(stats["pval"])
    stats["pval_BH"] = pval_adj

    # find SNR threshold 
    thresholds = np.quantile(null_distribution.loc[sizes,:].values,q=1-pval,axis=1)
    size_snr_trend = get_trend(sizes, thresholds, plot= False,verbose = verbose)
    stats["SNR_threshold"] = stats["size"].apply(lambda x: size_snr_trend(x))

    if save:        
        # save binarized data
        fpath= "/".join(bin_exprs_fname.split("/")[:-1])
        if not os.path.exists(fpath):
            os.makedirs(fpath)
            
        if not os.path.exists(bin_exprs_fname):
            binarized_data.to_csv(bin_exprs_fname, sep ="\t")
            if verbose:
                print("Binarized gene expressions are saved to",bin_exprs_fname,file = sys.stdout)

        # save binarization statistics
        if not os.path.exists(bin_stats_fname):
            stats.to_csv(bin_stats_fname, sep ="\t")
            if verbose:
                print("Statistics is saved to",bin_stats_fname,file = sys.stdout)

        # save null distribution: null_distribution, size,threshold 
        if not os.path.exists(bin_bg_fname):
            null_distribution.to_csv(bin_bg_fname, sep ="\t")
            if verbose:
                print("Background sitribution is saved to",bin_bg_fname,file = sys.stdout)
            
    
    ### keep features passed binarization
    passed = stats.loc[stats["SNR"]>stats["SNR_threshold"],:]
    #passed = stats.loc[stats["pval_BH"]<pval,:]
    
    if verbose:
        print("\t\tUP-regulated features:\t%s"%(passed.loc[passed["direction"]=="UP",:].shape[0]),file = sys.stdout)
        print("\t\tDOWN-regulated features:\t%s"%(passed.loc[passed["direction"]=="DOWN",:].shape[0]),file = sys.stdout)
        #print("\t\tambiguous features:\t%s"%(passed.loc[passed["direction"]=="UP,DOWN",:].shape[0]),file = sys.stdout)

    # keep only binarized features
    binarized_data = binarized_data.loc[:,list(passed.index.values)]
    # add sample names 
    binarized_data.index = exprs.columns.values
        
    if plot_all:
        # plot null distribution of SNR(size) - for shuffled genes n_permutation times for each size
        #e_stats = []
        #for s in sizes:
        #    for i in null_distribution.T.head(100).index.values:
        #        e_stats.append({"size":s,"SNR":null_distribution.loc[s,i]})
        #e_stats = pd.DataFrame.from_records(e_stats)
        # downsample to 1000 points in total
        #if e_stats.shape[0]>1000:
        #    e_stats = e_stats.sample(n=1000)
        #if verbose:
        #    print(e_stats.shape[0],"points from null distribution plotted",file = sys.stderr)
        
        #tmp  = plt.figure(figsize=(20,10))
        #tmp = plt.scatter(e_stats["size"],e_stats["SNR"],alpha = 0.2, color = "grey",label='background dist.')
        fig, ax  = plt.subplots(figsize = (10,4.5))
        
        # plot binarization results for real genes
        passed = stats.loc[stats["SNR"]> stats["SNR_threshold"],:]
        failed = stats.loc[stats["SNR"]<= stats["SNR_threshold"],:]
        tmp = ax.scatter(failed["size"],failed["SNR"],alpha = 1, color = "black",label='not passed')
        for i, txt in enumerate(failed["size"].index.values):
            ax.annotate(txt, (failed["size"][i], failed["SNR"][i]+0.1), fontsize=18)
        tmp = ax.scatter(passed["size"],passed["SNR"],alpha = 1, color = "red",label='passed')
        for i, txt in enumerate(passed["size"].index.values):
            ax.annotate(txt, (passed["size"][i], passed["SNR"][i]+0.1), fontsize=18)

        
        
        # plot cutoff 
        tmp = ax.plot(sizes,[size_snr_trend(x) for x in sizes], color="darkred",lw=2, ls = "--",
                       label="e.pval<"+str(pval))
        plt.gca().legend(('not passed','passed',"e.pval<"+str(pval)), fontsize=18)
        tmp = ax.set_xlabel("n_samples", fontsize=18)
        tmp = ax.yaxis.tick_right()
        tmp = ax.set_ylabel("SNR", fontsize=18)
        tmp = ax.set_ylim((0,4))
        #tmp = plt.savefig("figs_binarization/dist.svg", bbox_inches='tight', transparent=True)
        tmp = plt.show()
        

    return binarized_data, stats, null_distribution

    
#### Cluster binarized genes #####


def run_WGCNA_iterative(binarized_expressions,tmp_prefix="",
              deepSplit=0,detectCutHeight=0.995, nt = "signed_hybrid",# see WGCNA documentation
              max_power = 10,precluster=False,
              verbose = False,rscr_path=False, rpath = ""):
    
    t0 = time()
    
    not_clustered = binarized_expressions.columns.values
    binarized_expressions_ = binarized_expressions.loc[:,:].copy()
    stop_condition = False
    
    modules = []
    i=0
    while len(not_clustered)>=3 and not stop_condition:
        binarized_expressions_ = binarized_expressions_.loc[:,not_clustered]
        
        m,not_clustered = run_WGCNA(binarized_expressions_,tmp_prefix=tmp_prefix,
                  deepSplit=deepSplit,detectCutHeight=detectCutHeight, nt = nt,
                  max_power = max_power, precluster=precluster,
                  verbose = verbose,rscr_path=rscr_path, rpath = rpath)
        if verbose:
            print("\t\t\tWGCNA iteration %s, modules:%s, not clustered:%s"%(i,len(m),len(not_clustered)),file=sys.stdout)
        modules += m
        # stop when the number of not clustred samples does not change
        if len(m)==0:
            stop_condition = True
            if verbose:
                print("\t\t\tWGCNA iterations terminated at step ",i,file=sys.stdout)
        
        i+=1
    return (modules, not_clustered)

def run_WGCNA(binarized_expressions,tmp_prefix="",
              deepSplit=0,detectCutHeight=0.995, nt = "signed_hybrid",# see WGCNA documentation
              max_power = 10,precluster=False,
              verbose = False,rscr_path=False, rpath = ""):
    t0 = time()
    # create unique suffix for tmp files
    from datetime import datetime
    now = datetime.now()
    fname = "tmpWGCNA_" + now.strftime("%y.%m.%d_%H:%M:%S")+".tsv"
    if len(tmp_prefix)>0:
        fname = tmp_prefix+"."+fname 
    
    if verbose:
        print("\t\tWGCNA pre-clustering:",precluster,file=sys.stdout)
    if precluster:
        precluster = "T"
    else:
        precluster = "F"
        
    deepSplit = int(deepSplit)
    if not deepSplit in [0,1,2,3,4]:
        print("deepSplit must be 1,2,3 or 4. See WGCNA documentation.",file = sys.stderr)
        return ([],[])
    if not 0<detectCutHeight<1:
        print("detectCutHeight must be between 0 and 1. See WGCNA documentation.",file = sys.stderr)
        return ([],[])
    if verbose:
        print("\tRunning WGCNA for", fname, "...",file = sys.stdout)
    if not rscr_path:
        # assume run_WGCNA.R is in the same folder
        rscr_path = "/".join(os.path.realpath(__file__).split("/")[:-1])+'/run_WGCNA.R'
    
    binarized_expressions_ = binarized_expressions.loc[:,:].copy()
    
    # add suffixes to duplicated feature names
    feature_names = binarized_expressions.columns.values
    duplicated_feature_ndxs = np.arange(binarized_expressions.shape[1])[binarized_expressions.columns.duplicated()]
    
    if len(duplicated_feature_ndxs)>0:
        new_feature_names = []
        for i in range(len(feature_names)):
            fn = feature_names[i]
            if i in duplicated_feature_ndxs:
                fn = str(fn)+"*"+str(i)
            new_feature_names.append(fn)
        print("\t\t%s duplicated feature names detected."%len(duplicated_feature_ndxs),file=sys.stdout)
        dup_fn_mapping = dict(zip(new_feature_names,feature_names))
        binarized_expressions_.columns = new_feature_names
        
    # replace spaces in feature names
    # otherwise won't parse R output
    feature_names = (binarized_expressions.columns.values)
    feature_names_with_space = [x for x in feature_names if " " in x]
    if len(feature_names_with_space)>0:
        print("\t\t%s feature names containing spaces will be replaced."%len(feature_names_with_space),file=sys.stdout)
        fn_mapping = {}
        fn_mapping_back = {}
        for fn in feature_names:
            if " " in fn:
                fn_ = fn.replace(" ","_")
                fn_mapping[fn] = fn_
                fn_mapping_back[fn_] = fn 
        binarized_expressions_ =binarized_expressions.rename(fn_mapping,axis="columns")
    
    # save binarized expression to a file
    binarized_expressions_.to_csv(fname, sep ="\t")
    
    # run Rscript
    if len(rpath)>0:
        rpath = rpath+"/"
        
    if verbose:
        print("\tR command line:",file = sys.stdout)
        print("\t"+" ".join([rpath+'Rscript', rscr_path, fname, 
                             str(deepSplit), str(detectCutHeight), nt,str(max_power),precluster]),file = sys.stdout)
    
    process = subprocess.Popen([rpath+'Rscript', rscr_path, fname, 
                                str(deepSplit), str(detectCutHeight), nt,str(max_power),precluster], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    #stdout = stdout.decode('utf-8')
    module_file = fname.replace(".tsv",".modules.tsv") #stdout.rstrip()
    try:
        modules_df = pd.read_csv(module_file,sep = "\t",index_col=0)
    except:
        #print("WGCNA output:", stdout, file = sys.stdout)
        stderr = stderr.decode('utf-8')
        if verbose:
            print("WGCNA error:", stderr, file = sys.stdout)
        modules_df = pd.DataFrame.from_dict({})
    if verbose:
        print("\tWGCNA runtime: modules detected in {:.2f} s.".format(time()-t0),file = sys.stdout)
    
    # read WGCNA output
    modules = []
    not_clustered = []
    module_dict = modules_df.T.to_dict()
    for i in module_dict.keys():
        genes =  module_dict[i]["genes"].strip().split()
        # change feature names if they were modified
        
        # return spaces in feature names back if necessary
        if len(feature_names_with_space)>0:
            for j in range(len(genes)):
                if genes[j] in fn_mapping_back.keys():
                    genes[j] = fn_mapping_back[genes[j]]
        # remove suffixes from duplicated feature names
        if len(duplicated_feature_ndxs)>0:
            for j in range(len(genes)):
                if genes[j] in dup_fn_mapping.keys():
                    genes[j] = dup_fn_mapping[genes[j]]
        
        if i == 0:
            not_clustered = genes
        else:
            modules.append(genes)
    
    # remove WGCNA input and output files
    try:
        os.remove(fname) 
        os.remove(module_file)
    except:
        pass
    
    if verbose:
        print("\tmodules: {}, not clustered features {} ".format(len(modules),len(not_clustered)),file = sys.stdout)
        #print(stdout,file = sys.stdout)
        if len(stderr)>0:
            print(stderr,file = sys.stderr)
            
    
    return (modules,not_clustered)


def run_Louvain(similarity, similarity_cutoffs = np.arange(0.33,0.95,0.05), m=False,
                verbose = True,plot=False,modularity_measure = "newman"):
    t0 = time()
    if similarity.shape[0] == 0:
        print("no features to cluster",file =sys.stderr)
        return [], [], None
    
    
    if verbose:
        print("\tRunning Louvain ...")
        
    from sknetwork.clustering import Louvain, modularity
    
    modularities = []
    feature_clusters = {}
    best_Q = np.nan
    for cutoff in similarity_cutoffs:
        # scan the whole range of similarity cutoffs 
        # e.g. [1/4;9/10] with step 0.5
        sim_binary = similarity.copy()
        sim_binary[sim_binary<cutoff] = 0
        sim_binary[sim_binary!= 0] = 1
        rsums = sim_binary.sum()
        non_zero_features = rsums[rsums>0].index
        sim_binary = sim_binary.loc[non_zero_features,non_zero_features]
        gene_names = sim_binary.index.values
        sparse_matrix = csr_matrix(sim_binary)
        labels = Louvain(modularity =modularity_measure).fit_transform(sparse_matrix)
        Q = modularity(sparse_matrix, labels)
        modularities.append(Q)
        # if binary similarity matrix contains no zeroes
        # bugfix for Louvain()
        if sim_binary.min().min()==1: 
            labels = np.zeros(len(labels))
        feature_clusters[cutoff] = labels
    
    # if similarity_cutoffs contains only one value, choose it as best_cutoff
    if len(similarity_cutoffs)==1:
        best_cutoff = similarity_cutoffs[0]
        best_Q = Q
        
    # find best_cutoff automatically 
    else:
        # check if modularity(cutoff)=const 
        if len(set(modularities)) == 1:
            best_cutoff = similarity_cutoffs[-1]
            best_Q = modularities[-1]
            labels = feature_clusters[best_cutoff]
        
        #  if modularity!= const, scan the whole range of similarity cutoffs 
        #  e.g. [1/4;9/10] with step 0.05
        else:
            # find the knee point in the dependency modularity(similarity curoff)
            from kneed import KneeLocator
            # define the type of the curve
            curve_type = "increasing"
            if modularities[0] >= modularities[-1]:
                curve_type = "decreasing"
            if verbose:
                print("\tcurve type:",curve_type,file=sys.stdout)
            # detect knee and choose the one with the highest modularity
            try:
                kn = KneeLocator (similarity_cutoffs, modularities, curve='concave', direction=curve_type, online=True)
                best_cutoff = kn.knee
                best_Q = kn.knee_y
                labels = feature_clusters[best_cutoff]
            except:
                print("Failed to identify similarity cutoff",file=sys.stderr)
                print("Similarity cutoff: set to ",similarity_cutoffs[0], file = sys.stdout)
                best_cutoff = similarity_cutoffs[0]
                best_Q = np.nan
                print("Modularity:",modularities, file = sys.stdout)
                if plot:
                    plt.plot(similarity_cutoffs,modularities, 'bx-')
                    plt.xlabel('similarity cutoff')
                    plt.ylabel('modularity')
                    plt.show()
                    #return [], [], None
            if m:
                # if upper threshold for modularity m is specified 
                # chose the lowest similarity cutoff at which modularity reaches >= m
                best_cutoff_m = 1
                for i in range(len(modularities)):
                    if modularities[i] >= m:
                        best_cutoff_m = similarity_cutoffs[i]
                        best_Q_m = modularities[i]
                        labels_m = feature_clusters[best_cutoff]
                        break
                if best_cutoff_m<best_cutoff:
                    best_cutoff = best_cutoff_m
                    best_Q = best_Q_m
                    labels = labels_m
                    
    if plot and len(similarity_cutoffs)>1:
        plt.plot(similarity_cutoffs,modularities, 'bx-')
        plt.vlines(best_cutoff, plt.ylim()[0], plt.ylim()[1], linestyles='dashed',color= "red")
        plt.xlabel('similarity cutoff')
        plt.ylabel('modularity')
        plt.show()
    
    modules = []
    not_clustered = []
    
    for label in set(labels):
        ndx = np.argwhere(labels==label).reshape(1,-1)[0]
        genes = gene_names[ndx]
        if len(genes)>1:
            modules.append(genes) 
        else:
            not_clustered.append(genes[0])
    if verbose:
        print("\tLouvain runtime: modules detected in {:.2f} s.".format(time()-t0),file = sys.stdout)
        print("\tmodules: {}, not clustered features {} ".format(len(modules),len(not_clustered)),file = sys.stdout)
        print("\t\tsimilarity cutoff: {:.2f} modularity: {:.3f}".format(best_cutoff,best_Q),file=sys.stdout)
    return modules, not_clustered, best_cutoff


def run_MCL(similarity, inflations = list(np.arange(11,30)/10)+[3,3.5,4,5],verbose = True):
    if similarity.shape[0] == 0:
        return [[], []]
    import markov_clustering as mc
    t0 = time()
    # MCL accepts positive similarity matrix
    gene_names = similarity.index.values
    sparse_matrix= csr_matrix(similarity.values)
    
    best_Q = -1
    best_result = None
    best_clusters = None
    best_inflation = None
    for inflation in inflations:
        result = mc.run_mcl(sparse_matrix, inflation=inflation) # run MCL with default parameters
        clusters = mc.get_clusters(result)    # get clusters
        Q = mc.modularity(matrix=result, clusters=clusters) # Newman modularity
        #print(inflation, Q)
        if Q > best_Q:
            best_Q = Q
            best_result = result
            best_clusters = clusters
            best_inflation = inflation
    modules = []
    not_clustered = []
    for cluster in set(best_clusters):
        genes = gene_names[np.array(cluster)]
        if len(genes) >1:
            modules.append(sorted(genes))
        else:
            not_clustered.append(genes[0])
            
    if verbose:
        print("\tn_clusters:",len(modules),"not clustered:",len(not_clustered))
        print("\t\tinflation:", best_inflation, "modularity:", round(best_Q,3), round(time()-t0),"s.")
    return [modules, not_clustered]

def get_similarity_jaccard_signif(df,qval=0.01,J=0.33):
    # based on jaccard similarity
    t0 =  time()
    genes = df.columns.values
    n_samples = df.shape[0]
    df = df.values
    results = []
    for i in range(df.shape[1]):
        for j in range(i+1,df.shape[1]):
            g1= df[:,i]
            g2= df[:,j]
            o = g1 * g2
            overlap = o.sum()
            u = g1 + g2 
            union = u[u >0].shape[0]
            jaccard = overlap/union
            g1_only = g1.sum() - overlap
            g2_only = g2.sum() - overlap
            p =  pvalue(overlap,g1_only,g2_only,n_samples-union)
            results.append({"i":genes[i],"j":genes[j],"J":jaccard,"p_val":p.right_tail,
                            "overlap":overlap,"union":union})
    results = pd.DataFrame.from_records(results)
    # correction form multiple testing
    bh_res, qvals = fdrcorrection(results["p_val"].values,alpha=0.05)
    results["q_val"] = qvals
    results = results.loc[results["q_val"]<=qval,:]
    results = results.loc[results["J"]>=J,:]
    similarity ={}
    for row in results.iterrows():
        g1 = row[1]["i"]
        g2 = row[1]["j"]
        jaccard = row[1]["J"]
        if g1 in similarity.keys():
            similarity[g1][g2] = jaccard
        else:
            similarity[g1] = {g2:jaccard}
        if g2 in similarity.keys():
            similarity[g2][g1] = jaccard
        else:
            similarity[g2] = {g1:jaccard}
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
    similarity = similarity.loc[sorted(similarity.index),sorted(similarity.index)]
    return similarity

def get_similarity_jaccard(binarized_data,verbose = True): #,J=0.5
    t0 =  time()
    genes = binarized_data.columns.values
    n_samples = binarized_data.shape[0]
    size_threshold = int(min(0.45*n_samples,(n_samples)/2-10))
    #print("size threshold",size_threshold)
    n_genes = binarized_data.shape[1]
    df = np.array(binarized_data.T, dtype=bool)
    results = np.zeros((n_genes,n_genes))
    for i in range(0,n_genes):
        results[i,i] = 1
        g1= df[i]
        
        for j in range(i+1,n_genes):
            g2= df[j]
            o = g1 * g2 
            u = g1 + g2 
            jaccard = o.sum()/u.sum()
            # try matching complements
            if g1.sum()>size_threshold:
                g1_complement = ~g1
                o = g1_complement * g2 
                u = g1_complement + g2 
                jaccard_c = o.sum()/u.sum()
            elif g2.sum()>size_threshold:
                g2 = ~g2
                o = g1 * g2 
                u = g1 + g2 
                jaccard_c = o.sum()/u.sum()
            else:
                jaccard_c = 0
            jaccard = max(jaccard,jaccard_c)
            results[i,j] = jaccard
            results[j,i] = jaccard

    results = pd.DataFrame(data = results,index = genes, columns = genes)
    if verbose:
        print("\tJaccard similarities for {} features computed in {:.2f} s.".format(binarized_data.shape[1],time()-t0),file=sys.stdout)     
    return results


def get_similarity_corr(df,verbose = True): 
    t0 = time()
    corr = df.corr() #.applymap(abs)
    corr = corr[corr>0] # to consider only direct correlations 
    corr = corr.fillna(0)
    if verbose:
        print("\tPearson's r similarities for {} features computed in {:.2f} s.".format(df.shape[1],time()-t0),file=sys.stdout)
    return corr

def make_TOM(similarity):
    #from https://www.rdocumentation.org/packages/WGCNA/versions/1.71/topics/TOMsimilarityFromExpr
    t0 = time()
    df = similarity.values
    fnames = similarity.index.values
    N = similarity.shape[0]
    TOM = np.zeros((N,N))
    for i in range(0,N):
        ai = df[i,]
        sum_ai = ai.sum()
        for j in range(i,N):
            aj = df[j,]
            sum_aj = aj.sum()
            aij = aj[i]
            f = (sum_aj+sum_ai)/2-aij
            a_product = np.dot(ai,aj)
            tom_ij = (aij + a_product- aij*aij)/(f+1-aij)
            TOM[i,j] = tom_ij
            TOM[j,i] = tom_ij
    TOM= pd.DataFrame(data=TOM,index=fnames,columns =fnames)
    print("\tTOM computed in {:.2f} s.".format(time()-t0),file =sys.stdout)
    return TOM


######## Make biclusters #########

def cluster_samples(data,min_n_samples=5,seed=0,method="kmeans"):
    # identify identify bicluster and backgound groups using 2-means
    max_n_iter = max(max(data.shape),500)
    if method == "kmeans" or method == "Jenks":
        labels = KMeans(n_clusters=2,
                        random_state=seed,
                        init="random",
                        n_init=10,
                        max_iter=max_n_iter).fit(data).labels_
    elif method == "ward":
        labels = AgglomerativeClustering(n_clusters=2, linkage='ward').fit(data).labels_
    #elif method == "HC_ward":
    #        model = Ward(n_clusters=2).fit(data).labels_
    elif method == "GMM":
        labels = GaussianMixture(n_components=2, init_params="kmeans",
                                 max_iter=max_n_iter, n_init = 5, 
                                 covariance_type = "spherical",
                                 random_state = seed).fit_predict(data)
    ndx0 = np.where(labels == 0)[0]
    ndx1 = np.where(labels == 1)[0]
    if min(len(ndx1),len(ndx0))< min_n_samples:
        return {}
    if len(ndx1) > len(ndx0):
        samples = ndx0
    else: 
        samples = ndx1

    bicluster = {"sample_indexes":set(samples),"n_samples":len(samples)}
    return bicluster

def modules2biclusters_2(modules, data_to_cluster,
                       method = "kmeans",
                       min_n_samples=5,
                       min_n_genes=2,seed=0,verbose = True):
    '''Identifies optimal sample set for each module: 
    splits samples into two sets in a subspace of each module
    '''
    t0 = time()
    biclusters = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0
    
    for mid in range(0,len(modules)):
        genes = modules[mid]
        if len(genes) >= min_n_genes:
            # cluster samples in a space of selected genes
            data = data_to_cluster.loc[genes,:].T
            bicluster = cluster_samples(data,min_n_samples=min_n_samples,seed=seed,method = method)
            if len(bicluster)>0:
                bicluster["id"] = i
                bicluster["genes"] = set(genes)
                bicluster["n_genes"] = len(bicluster["genes"])
                
                bic = data_to_cluster.loc[genes,:].iloc[:,sorted(bicluster["sample_indexes"])]
                bg = data_to_cluster.loc[genes,:].iloc[:,[ndx for ndx in range(data_to_cluster.shape[1]) if not ndx in bicluster["sample_indexes"]]]
                bicluster["genes_up"] = bic[bic.mean(axis=1) >= bg.mean(axis=1)].index.values
                bicluster["genes_down"] = bic[bic.mean(axis=1) < bg.mean(axis=1)].index.values
                
                if len (bicluster["genes_up"])>0 and len(bicluster["genes_down"])>0:
                    bic_profile = data_to_cluster.loc[bicluster["genes_up"],:].mean(axis=0) - data_to_cluster.loc[sorted(bicluster["genes_down"]),:].mean(axis=0)
                else:
                    bic_profile = data_to_cluster.loc[bicluster["genes"],:].mean(axis=0)
                bicluster2 = cluster_samples(bic_profile.T.values.reshape(-1, 1),min_n_samples=min_n_samples,seed=seed,method = method)
                bicluster["sample_indexes"]= bicluster2["sample_indexes"]
                bicluster["n_samples"]= bicluster2["n_samples"]
                biclusters[i] = bicluster
                i+=1
      
    if verbose:
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %(len(modules),round(time()-t0,2)))
        print("Passed biclusters (>=%s genes, >= samples %s): %s"%(min_n_genes,min_n_samples,i-1), file = sys.stdout)
        print()
    
    return biclusters

def modules2biclusters(modules, data_to_cluster,
                       method = "kmeans",
                       min_n_samples=5,
                       min_n_genes=2,seed=0,verbose = True):
    '''Identifies optimal sample set for each module: 
    splits samples into two sets in a subspace of each module
    '''
    t0 = time()
    biclusters = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0
    
    for mid in range(0,len(modules)):
        genes = modules[mid]
        if len(genes) >= min_n_genes:
            # cluster samples in a space of selected genes
            data = data_to_cluster.loc[genes,:].T
            bicluster = cluster_samples(data,min_n_samples=min_n_samples,seed=seed,method = method)
            if len(bicluster)>0:
                bicluster["id"] = i
                bicluster["genes"] = set(genes)
                bicluster["n_genes"] = len(bicluster["genes"])
                biclusters[i] = bicluster
                i+=1
      
    if verbose:
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %(len(modules),round(time()-t0,2)))
        print("Passed biclusters (>=%s genes, >= samples %s): %s"%(min_n_genes,min_n_samples,i-1), file = sys.stdout)
        print()
    
    return biclusters

def update_bicluster_data(bicluster,data):
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
    bg_samples = [x for x in sample_names if not x in bic_samples]
    bicluster["samples"] = set(bic_samples)
    bicluster["gene_indexes"] = set([np.where(gene_names==x)[0][0] for x in bicluster["genes"]])

    # distinguish up- and down-regulated features 
    m_bic = data.loc[bic_genes,bic_samples].mean(axis=1)
    m_bg = data.loc[bic_genes,bg_samples].mean(axis=1)
    genes_up = m_bic[m_bic>=m_bg].index.values
    genes_down = m_bic[m_bic<m_bg].index.values
    bicluster["genes_up"] = set(genes_up)
    bicluster["genes_down"] = set(genes_down)

    genes_up = m_bic[m_bic>=m_bg].index.values
    genes_down = m_bic[m_bic<m_bg].index.values
    
    # calculate average z-score for each sample
    if min(len(genes_up),len(genes_down))>0: # take direction into account
        avg_zscore = (data.loc[genes_up,:].sum()-data.loc[genes_down,:].sum())/bicluster["n_genes"]
    else:
        avg_zscore = data.loc[list(bicluster["genes"]),:].mean()
    
    # compute SNR for average z-score for this bicluster
    m = avg_zscore[bic_samples].mean() - avg_zscore[bg_samples ].mean()
    s = avg_zscore[bic_samples].std() + avg_zscore[bg_samples].std()
    snr = np.abs(m)/s
    bicluster["SNR"] = snr
    return bicluster

def merge_biclusters(biclusters,data,J=0.8,
                     min_n_samples=5,
                     seed=42,
                     method = "kmeans",verbose=True):
    #  bicluaters -> binary -> jaccard sim
    binary_representation = {}
    N = data.shape[1]
    for i in biclusters.keys():
        b = np.zeros(N)
        s_ndx = list(biclusters[i]["sample_indexes"])
        b[s_ndx] = 1
        binary_representation[i] = b
    binary_representation = pd.DataFrame.from_dict(binary_representation)
    binary_representation.index=data.columns.values
    bic_similarity = get_similarity_jaccard(binary_representation, verbose = verbose)
    #bic_similarity[bic_similarity >= J] = 1
    #bic_similarity[bic_similarity < J] = 0
    # find groups of biclusters including the same sample sets
    merged, not_merged, similarity_cutoff = run_Louvain(bic_similarity,verbose = False,plot=False, similarity_cutoffs=[J])
    if len(merged)==0 and verbose:
        print("No biclusters to merge",file = sys.stdout)
        return biclusters 
    
    merged_biclusters = {}
    # add biclusters with no changes
    for bic_id in not_merged: 
        merged_biclusters[bic_id] = biclusters[bic_id]
    
    # merge biclusters with overlapping sample sets
    for bic_group in merged:
        bic_group = sorted(bic_group)
        if verbose:
            print("merged biclustres",bic_group,file = sys.stderr)
        new_bicluster = biclusters[bic_group[0]]
        # update genes
        for bic_id in bic_group[1:]:
            bic2 = biclusters[bic_id]
            new_bicluster["genes"] = new_bicluster["genes"] | bic2["genes"]
            new_bicluster["n_genes"] = len(new_bicluster["genes"])
        # update sample set for new bicluster
        # cluster samples in a space of selected genes
        new_bicluster.update(cluster_samples(data.loc[list(new_bicluster["genes"]),:].T,
                                    min_n_samples=min_n_samples,seed=seed,method = method))
        new_bicluster["n_samples"] = len(new_bicluster['sample_indexes'])
        merged_biclusters[bic_group[0]] = new_bicluster
    return merged_biclusters

def make_biclusters(feature_clusters,
                    binarized_data,data,
                    null_distribution,
                    merge = 1.01,
                    min_n_samples=10, min_n_genes=2,
                    method = "kmeans",
                    seed = 42,
                    cluster_binary = False,verbose= True):
    sample_names = data.columns.values
    gene_names = data.index.values
    biclusters = []

    if cluster_binary:
        data_to_cluster = binarized_data.loc[:,:].T # binarized expressions
    else:
        data_to_cluster = data.loc[binarized_data.columns.values,:] # z-scores

    if len(feature_clusters)==0:
        print("No biclusters found.", file = sys.stderr)
    else:
        biclusters = modules2biclusters(feature_clusters, data_to_cluster,method = method,
                                min_n_samples=min_n_samples, min_n_genes=2,
                                verbose = False,seed=seed)
        
        ### merge biclusters with highly similar sample sets
        if merge<=1.0:
            biclusters = merge_biclusters(biclusters,data, method = method, J=merge, 
                                          min_n_samples=min_n_samples,seed=seed,verbose = verbose)
        
        for i in list(biclusters.keys()):
            biclusters[i] = update_bicluster_data(biclusters[i],data)
        
    biclusters = pd.DataFrame.from_dict(biclusters).T
    # add direction
    biclusters["direction"] = "BOTH"
    biclusters.loc[biclusters["n_genes"] == biclusters["genes_up"].apply(len),"direction"] = "UP"
    biclusters.loc[biclusters["n_genes"] == biclusters["genes_down"].apply(len),"direction"] = "DOWN"
    
    # add p-value for bicluster SNR (computed for avg. zscores) 
    # use the same distribution as for single features
    #biclusters["e_pval"] = biclusters.apply(lambda row: calc_e_pval(row["SNR"], row["n_samples"], null_distribution),axis=1) 
    
    # sort and reindex
    biclusters = biclusters.sort_values(by=["SNR","n_genes"],ascending=[False,False])
    biclusters.index = range(0,biclusters.shape[0])
    
    biclusters = biclusters.loc[:,["SNR","n_genes","n_samples","genes","samples","direction",
                                   "genes_up","genes_down","gene_indexes","sample_indexes"]]
    
    
    return biclusters


#### consensus biclusters ####

def calc_bicluster_similarities(biclusters,exprs,
                                similarity="both",plot=True,
                                figsize=(17, 17),labels=False,
                                colorbar_off=True): 
    
    biclusters_dict = biclusters.T.to_dict()
    
    if similarity not in ["genes","samples","both"]:
        print("Similarity must be 'genes','samples','both'. Set to 'both'.",file=sys.stderr)
        similarity = "both"
        
    J_heatmap = {}
    s = set(exprs.columns.values)
    g = set(exprs.index.values)
    N_bics = len(biclusters.keys())
    ### plot pairwise comparisons:
    for bic_n1 in biclusters_dict.keys():
        b1 = biclusters_dict[bic_n1]
        g1 = b1["genes"]
        s1 = b1["samples"]
        J_heatmap[bic_n1] = {}
        for bic_n2 in biclusters_dict.keys():
            b2 = biclusters_dict[bic_n2]
            g2 = b2["genes"]
            s2 = b2["samples"]
            g_overlap = g1.intersection(g2)
            g_union = g1.union(g2)
            s_overlap = s1.intersection(s2)
            s_union = s1.union(s2)
            J_g = len(g_overlap)/len(g_union)
            J_s = len(s_overlap)/len(s_union)
            J_heatmap[bic_n1][bic_n2] = 0
                    
            # significance of gene overlap
            if similarity !="samples":
                if len(g_overlap)==0:
                    p_g = 1
                else:
                    if len(g)<5000:
                        # if not too many features, use Fisher's exact
                        p_g =  pvalue(len(g_overlap),len(g1.difference(g2)),len(g2.difference(g1)),len(g.difference(g1|g2)))
                        p_g = p_g.right_tail
                    else:
                        # otherwise replacing exact Fisher's with chi2 
                        chi2, p_g, dof, expected = chi2_contingency([[len(g_overlap), len(g1.difference(g2))],
                                                                     [len(g2.difference(g1)), len(g.difference(g1|g2))]]) 
            
            # significance of sample overlap
            if similarity !="genes":
                # skip if similarity==both and gene overlap is empty
                if similarity =="both" and len(g)==0: 
                    p_s = 1
                else:
                    if len(s)<5000:
                        p_s =  pvalue(len(s_overlap),len(s1.difference(s2)),len(s2.difference(s1)),len(s.difference(s1|s2)))
                        # for similarity==samples left-tail p-val matters too
                        # e.g. for biclusters constaining about a half of all samples
                        p_s = min(p_s.right_tail,p_s.left_tail)
                    else:
                        chi2, p_s, dof, expected = chi2_contingency([[len(s_overlap), len(s1.difference(s2))],
                                                                         [len(s2.difference(s1)), len(s.difference(s1|s2))]]) 
            
            
            if similarity =="genes":
                if p_g*(N_bics-1)*N_bics/2<0.05:
                    J_heatmap[bic_n1][bic_n2] = J_g
            elif similarity =="samples":
                if p_s*(N_bics-1)*N_bics/2<0.05:
                    J_heatmap[bic_n1][bic_n2] = J_s
            elif  similarity == "both": # both genes and samples considered
                # consider significant overlaps in g and s and save max. J
                if (p_g*(N_bics-1)*N_bics/2<0.05 and len(s_overlap)>0) or (p_s*(N_bics-1)*N_bics/2<0.05 and len(g_overlap)>0):
                        J_heatmap[bic_n1][bic_n2] = max(J_s,J_g)
                    
    J_heatmap = pd.DataFrame.from_dict(J_heatmap)
    
    if plot:
        import seaborn as sns
        g = sns.clustermap(J_heatmap,
                           yticklabels=labels, xticklabels=labels, 
                           linewidths=0,vmin=0,vmax=1,
                           figsize=figsize,center=0,annot=labels)
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
        if colorbar_off:
            g.cax.set_visible(False)
        plt.show()
        
    return J_heatmap
    
def make_consensus_biclusters(biclusters_list,exprs, frac_runs=0.5,
                              similarity = "both", # can be 'both','genes','samples' 
                              min_similarity = 0.33,
                              method="kmeans", modularity_measure = "potts",
                              min_n_genes =2, min_n_samples=5,
                              seed = -1, plot = False,
                              figsize=(17, 17),labels=False,colorbar_off=True):
    t0 = time()
    
    if seed == -1:
        seed = random.randint(0,1000000)
        print("Seed for sample clustering: %s"%(seed),file=sys.stderr)
    
    # list of biclusters from several runs
    n_runs = len(biclusters_list)
    biclusters = pd.concat(biclusters_list)
    if not "detected_n_times" in biclusters.columns.values:
        biclusters["detected_n_times"] = 1
    biclusters.index = range(biclusters.shape[0])
    
    # calculate pairwise bicluster similarity
    J_heatmap = calc_bicluster_similarities(biclusters,exprs,
                                            similarity = similarity,
                                            plot=plot,
                                            figsize=figsize,labels=labels,
                                            colorbar_off=colorbar_off)
    
    # if all biclusters are exactly the same
    if J_heatmap.min().min()==1:
        # return the first bicluster
        consensus_biclusters = biclusters.iloc[[0],:].copy()
        consensus_biclusters.index = [0]
        consensus_biclusters.loc[0,'detected_n_times'] = biclusters.shape[0]
        print("all biclusters are exactly the same",file = sys.stderr)
        return consensus_biclusters

    
    t1= time()
    print("%s s for similarity matrix"%round(t1-t0))
    
    # cluster biclusters by similarity
    matched, not_matched, cutoff = run_Louvain(J_heatmap,
                                               similarity_cutoffs = np.arange(min_similarity,0.9,0.05),m=False,
                                               verbose = True, plot=plot,
                                               modularity_measure = modularity_measure)
    t2 = time()
    #print(round(t2-t1),"s for Louvain ")
    
    
    # make consensus biclusters
    # for each group of matched biclusters, keep genes occuring at least n times
    # cluster samples again in a subspace of a new gene set
    consensus_biclusters = []

    # for each group of matched biclusters 
    for i in range(len(matched)):
        gsets = biclusters.loc[matched[i],"genes"].values
        detected_n_times = biclusters.loc[matched[i],"detected_n_times"].sum()
        # count gene occurencies
        gene_occurencies = {}
        for gene in set().union(*gsets):
            gene_occurencies[gene] = 0
            for gset in gsets:
                if gene in gset:
                    gene_occurencies[gene]+=1
        
        gene_occurencies = pd.Series(gene_occurencies).sort_values()
        passed_genes = sorted(gene_occurencies[gene_occurencies>=min(n_runs,len(matched))*frac_runs].index.values)
        not_passed_genes  = sorted(gene_occurencies[gene_occurencies<min(n_runs,len(matched))*frac_runs].index.values)

        
        if len(passed_genes)<min_n_genes:
            # move all biclusters to not matched
            not_matched += list(matched[i])
        else:
            # cluster samples in a new gene set
            bicluster = cluster_samples(exprs.loc[passed_genes,:].T,min_n_samples=min_n_samples,seed=seed,method=method)
            # if bicluster is not empty, add it to consenesus
            if "sample_indexes" in bicluster.keys():
                bicluster["genes"] = set(passed_genes)
                bicluster["n_genes"] = len(bicluster["genes"])
                bicluster = update_bicluster_data(bicluster,exprs)
                bicluster["detected_n_times"] = detected_n_times
                consensus_biclusters.append(bicluster)
    consensus_biclusters = pd.DataFrame.from_records(consensus_biclusters)
     
    # add not matched
    not_changed_biclusters = biclusters.loc[not_matched,:]
    #not_changed_biclusters["detected_n_times"] = 1
    consensus_biclusters = pd.concat([consensus_biclusters,not_changed_biclusters])

    # add direction
    consensus_biclusters["direction"] = "BOTH"
    consensus_biclusters.loc[consensus_biclusters["n_genes"] == consensus_biclusters["genes_up"].apply(len),"direction"] = "UP"
    consensus_biclusters.loc[consensus_biclusters["n_genes"] == consensus_biclusters["genes_down"].apply(len),"direction"] = "DOWN"

    # sort
    col_order = ['SNR','n_genes', 'n_samples', 'genes',  'samples',  'genes_up', 'genes_down',
                 'gene_indexes','sample_indexes',"direction",'detected_n_times']
    consensus_biclusters = consensus_biclusters.sort_values(by=["SNR","n_samples"],ascending=[False,True])
    consensus_biclusters = consensus_biclusters.loc[:,col_order]
    consensus_biclusters = consensus_biclusters.loc[consensus_biclusters["n_samples"]>=min_n_samples,:]
    
    consensus_biclusters.index = range(consensus_biclusters.shape[0])
    
    print(round(time()-t2),"s for making consensus biclusters from consensus gene sets")
    
    return consensus_biclusters


def make_consensus_biclusters2(biclusters_list,exprs,
                              similarity = "samples", # can be 'both','genes','samples' 
                              frac_runs=1/3,min_n_genes =2, min_n_samples=5,
                              p=0.05,
                              min_similarity = 0.33,
                              method="kmeans",
                              modularity_measure = "potts",
                              seed = -1, plot = False, verbose = True,
                              figsize=(10, 10),labels=False,colorbar_off=True):
    
    from utils.eval import find_best_matching_biclusters
    
    t0 = time()
    N = exprs.shape[0]
    
    if seed == -1:
        seed = random.randint(0,1000000)
        print("Seed for sample clustering: %s"%(seed),file=sys.stderr)
    
    # list of biclusters from several runs
    n_runs = len(biclusters_list)
    for i in range(n_runs):
        biclusters_list[i]["run"] = i
    biclusters = pd.concat(biclusters_list)
    biclusters.index = range(biclusters.shape[0])
    bic_ids = biclusters.index.values
    if not "detected_n_times" in biclusters.columns.values:
        biclusters["detected_n_times"] = 1
    
    n_bics = biclusters.shape[0]
    J_heatmap = pd.DataFrame(np.zeros((n_bics, n_bics)),index=bic_ids,columns=bic_ids)
    
    # add only best matches to Jaccard similarity matrix
    avg_J_sim = {}
    for i in range(n_runs):
        bics1 = biclusters.loc[biclusters["run"]==i,:]
        J_heatmap.loc[bics1.index.values,bics1.index.values] = np.identity(bics1.shape[0])
        avg_J_sim[i]={i:1}
        for j in range(n_runs):
            if i!=j:
                bics2 = biclusters.loc[biclusters["run"]==j,:]
                # find best matches between bics1 and bics2
                bm = find_best_matching_biclusters(bics1,bics2,exprs.shape,
                                                   by = similarity,min_g=min_n_genes,adj_pval_thr=p)
                bm = bm.dropna()
                if bm.shape[0]>0:
                    avg_J_sim[i][j] = np.mean(bm["J"])
                    df = bm.loc[:,["bm_id","J"]]
                    for row in df.iterrows():
                        J_heatmap.loc[row[0],row[1]["bm_id"]] += row[1]["J"]/2
                        J_heatmap.loc[row[1]["bm_id"],row[0]] += row[1]["J"]/2
                        
    
    # plot similarity heatmaps
    if plot: 
        import seaborn as sns
        
        avg_J_sim = pd.DataFrame.from_dict(avg_J_sim)
        avg_J_sim= avg_J_sim.fillna(0)
        if avg_J_sim.shape[0]>0:
            g = sns.clustermap(avg_J_sim,linewidths=0,vmin=0,vmax=1,figsize=(3,3),
                               center=0,annot=True)
            g.ax_cbar.set_visible(False)
            plt.show()
        
        labels = True
        if len(bic_ids)>20:
            labels = False
        g = sns.clustermap(J_heatmap,
                           yticklabels=labels, xticklabels=labels, 
                           linewidths=0,vmin=0,vmax=1,
                           figsize=figsize,center=0,annot=labels)
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
        if colorbar_off:
            g.cax.set_visible(False)
        plt.show()
        
    # if all biclusters are exactly the same
    if J_heatmap.min().min()==1:
        # return the first bicluster
        consensus_biclusters = biclusters.iloc[[0],:].copy()
        consensus_biclusters.index = [0]
        consensus_biclusters.loc[0,'detected_n_times'] = biclusters.shape[0]
        print("all biclusters are exactly the same",file = sys.stderr)
        return consensus_biclusters

    
    t1= time()
    print("%s s for similarity matrix"%round(t1-t0))
    
    # cluster biclusters by similarity
    matched, not_matched, cutoff = run_Louvain(J_heatmap,
                                               similarity_cutoffs = np.arange(min_similarity,0.9,0.05),m=False,
                                               verbose = True, plot=plot,
                                               modularity_measure = modularity_measure)
    t2 = time()
    
    # make consensus biclusters
    # for each group of matched biclusters, keep genes occuring at least n times
    # cluster samples again in a subspace of a new gene set
    consensus_biclusters = []

    # for each group of matched biclusters 
    for i in range(len(matched)):
        gsets = biclusters.loc[matched[i],"genes"].values
        detected_n_times = biclusters.loc[matched[i],"detected_n_times"].sum()
        # count gene occurencies
        gene_occurencies = {}
        for gene in set().union(*gsets):
            gene_occurencies[gene] = 0
            for gset in gsets:
                if gene in gset:
                    gene_occurencies[gene]+=1
        
        gene_occurencies = pd.Series(gene_occurencies).sort_values()
        passed_genes = sorted(gene_occurencies[gene_occurencies>=min(n_runs,len(matched))*frac_runs].index.values)
        not_passed_genes  = sorted(gene_occurencies[gene_occurencies<min(n_runs,len(matched))*frac_runs].index.values)

        
        if len(passed_genes)<min_n_genes:
            # move all biclusters to not matched
            not_matched += list(matched[i])
        else:
            # cluster samples in a new gene set
            bicluster = cluster_samples(exprs.loc[passed_genes,:].T,min_n_samples=min_n_samples,seed=seed,method=method)
            # if bicluster is not empty, add it to consenesus
            if "sample_indexes" in bicluster.keys():
                bicluster["genes"] = set(passed_genes)
                bicluster["n_genes"] = len(bicluster["genes"])
                bicluster = update_bicluster_data(bicluster,exprs)
                bicluster["detected_n_times"] = detected_n_times
                consensus_biclusters.append(bicluster)
    consensus_biclusters = pd.DataFrame.from_records(consensus_biclusters)
     
    # add not matched
    not_changed_biclusters = biclusters.loc[not_matched,:]
    #not_changed_biclusters["detected_n_times"] = 1
    consensus_biclusters = pd.concat([consensus_biclusters,not_changed_biclusters])

    # add direction
    consensus_biclusters["direction"] = "BOTH"
    consensus_biclusters.loc[consensus_biclusters["n_genes"] == consensus_biclusters["genes_up"].apply(len),"direction"] = "UP"
    consensus_biclusters.loc[consensus_biclusters["n_genes"] == consensus_biclusters["genes_down"].apply(len),"direction"] = "DOWN"

    # sort
    col_order = ['SNR','n_genes', 'n_samples', 'genes',  'samples',  'genes_up', 'genes_down',
                 'gene_indexes','sample_indexes',"direction",'detected_n_times']
    consensus_biclusters = consensus_biclusters.sort_values(by=["SNR","n_samples"],ascending=[False,True])
    consensus_biclusters = consensus_biclusters.loc[:,col_order]
    consensus_biclusters = consensus_biclusters.loc[consensus_biclusters["n_samples"]>=min_n_samples,:]
    
    consensus_biclusters.index = range(consensus_biclusters.shape[0])
    
    print(round(time()-t2),"s for making consensus biclusters from consensus gene sets")
    
    return consensus_biclusters
#### reading and writing #####

def read_bic_table(file_name, parse_metadata = False):
    if not os.path.exists(file_name):
        return pd.DataFrame()
    biclusters = pd.read_csv(file_name,sep = "\t",index_col=0,comment="#")
    if len(biclusters) ==0:
        return pd.DataFrame()
    else:
        biclusters.loc[:,["genes_up","genes_down"]] = biclusters.loc[:, ["genes_up","genes_down"]].fillna("")
        biclusters["genes"] = biclusters["genes"].apply(lambda x: set([g for g in x.split(" ") if not g =='']))
        biclusters["genes_up"] = biclusters["genes_up"].apply(lambda x: set([g for g in x.split(" ") if not g =='']))
        biclusters["genes_down"] = biclusters["genes_down"].apply(lambda x: set([g for g in x.split(" ") if not g =='']))
        biclusters["samples"] = biclusters["samples"].apply(lambda x: set(x.split(" ")))
        biclusters["gene_indexes"] = biclusters["gene_indexes"].apply(lambda x: set(map(int, set(x.split(" ")))))
        biclusters["sample_indexes"] = biclusters["sample_indexes"].apply(lambda x: set(map(int, set(x.split(" ")))))
    
    #resulting_bics.set_index("id",inplace=True)
    if parse_metadata:
        f = open(file_name, 'r')
        metadata = f.readline()
        f.close()
        if metadata.startswith("#"):
            metadata = metadata.replace("#","").rstrip()
            metadata = metadata.split("; ") 
            metadata = dict([x.split("=") for x in metadata])
            return biclusters, metadata 
    return biclusters

def write_bic_table(bics_dict_or_df, results_file_name,to_str=True,
                    add_metadata=False, 
                    seed = None, min_n_samples =None,
                    bin_method = None, clust_method = None, pval = None,
                    directions = [],
                    alpha=None, beta_K = None, similarity_cutoff = None,
                    ds = None, dch = None, m=None,
                   max_power=None, precluster=None, merge = None):
    if add_metadata:
        metadata = "#seed="+str(seed)+"; "+"pval="+str(pval)+"; "+"min_n_samples="+str(min_n_samples)+"; "
        metadata = metadata + "b="+bin_method+"; "
        metadata = metadata + "c="+clust_method+"; "
        if len(directions):
            metadata = metadata + "directions="+"-".join(directions)+"; "
        if clust_method == "Louvain":
            metadata = metadata + "simiarity_cutoff="+str(similarity_cutoff)+"; modularity="+str(m)
        elif "WGCNA" in clust_method:
            metadata = metadata + "ds="+str(ds)+"; dch="+str(dch)+"; max_power="+str(max_power)+"; precluster="+str(precluster)
        #elif clust_method == "DESMOND":
        #    metadata = metadata + "alpha="+str(alpha)+"; " + "beta_K="+str(beta_K)
        else:
            print("Unknown 'clust_method'",clust_method,file= sys.stderr)
        metadata = metadata + "; merge="+str(merge)
        with open(results_file_name, 'w') as f:
            f.write(metadata+"\n")
        write_mode = 'a'
    else:
        write_mode = 'w'
            
    bics = bics_dict_or_df.copy()
    if len(bics) ==0:
        print("No biclusters found",file=sys.stderr)
    else:
        if not type(bics) == type(pd.DataFrame()):
            bics = pd.DataFrame.from_dict(bics)
        if to_str:
            bics["genes"] = bics["genes"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["genes_up"] = bics["genes_up"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["genes_down"] = bics["genes_down"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["samples"] = bics["samples"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["gene_indexes"] = bics["gene_indexes"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["sample_indexes"] = bics["sample_indexes"].apply(lambda x:" ".join(map(str,sorted(x))))
        #bics.index = range(0,bics.shape[0])
        bics.index.name = "id"
        cols =  bics.columns.values
        #first_cols = ["SNR","e_pval","n_genes","n_samples","direction","genes","samples"]
        #bics = bics.loc[:,first_cols+sorted(list(set(cols).difference(first_cols)))]
    bics.to_csv(results_file_name ,sep = "\t", mode = write_mode)
