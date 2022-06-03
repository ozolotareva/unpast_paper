import sys
import os
import subprocess
import random
import warnings
import pandas as pd
import numpy as np
from time import time

from sklearn.mixture import GaussianMixture
from scipy.interpolate import interp1d
from scipy.sparse.csr import csr_matrix
from sknetwork.clustering import Louvain, modularity
import markov_clustering as mc
import jenkspy

from sklearn.cluster import KMeans
from method2 import calc_bic_SNR , identify_opt_sample_set

import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection


def validate_input_matrix(exprs, tol=0.01):
    m = exprs.mean(axis=1)
    std = exprs.std(axis=1)
    mean_passed = np.all(np.abs(m)<tol)
    std_passed = np.all(np.abs(std-1)<tol)
    if not (mean_passed and std_passed):
        print("Input is not standardized.",file = sys.stderr)
        return False
    if len(set(exprs.index.values)) < exprs.shape[0]:
        print("Row names are not unique.",file = sys.stderr)
        return False
    return True

def calc_SNR(ar1, ar2):
    return (np.mean(ar1) - np.mean(ar2)) / (np.std(ar1) + np.std(ar2))

######### Binarization #########

def random_splits(exprs, min_n_samples,n_permutations = 1000, seed = 42,verbose = True,pval = 0.001):
    # empirical distribuition of SNR depending on the bicluster size 
    # samples N values from expression matrix, and split them into bicluster  background
    # returns distributions
    n_permutations = max(1000,int(1.0/pval))*5
    N = exprs.shape[1]
    values = exprs.values.reshape(-1)
    
    sizes = np.arange(min_n_samples,int(N/2)+1)
    if verbose:
        print("\tGenerate empirical distribuition of SNR depending on the bicluster size ...")
        print("\t\ttotal samples: %s,\n\t\tnumber of samples in a bicluster: %s - %s,\n\t\tn_permutations: %s"%(N,sizes[0],sizes[-1],n_permutations))
        print("snr_pval threshold:",pval)
    empirical_snr = np.zeros((sizes.shape[0],n_permutations))
    #thresholds = np.zeros(sizes.shape[0])
    np.random.seed(seed=seed)
    for s in sizes:
        snrs = np.zeros(n_permutations)
        for i in range(0,n_permutations):
            x = np.random.normal(size = exprs.shape[1]) # np.random.choice(values,size=exprs.shape[1]) #
            x.sort()
            #x = (x - np.mean(x))/np.std(x)
            empirical_snr[s-min_n_samples,i] = calc_SNR(x[s:], x[:s])
    thresholds = np.quantile(empirical_snr,q=1-pval,axis=1)
    return sizes,thresholds, empirical_snr

def get_trend(sizes, thresholds, plot= True):
    # smoothens the trend and retunrs a function min_SNR(size; p-val. cutoff)
    print("\t\t\tfraction:",round(min(0.1,15/len(sizes)),2))
    lowess = sm.nonparametric.lowess
    lowess_curve = lowess(thresholds,sizes,frac=min(0.1,15/len(sizes)),return_sorted=True,is_sorted=False)
    get_min_snr = interp1d(lowess_curve[:,0],lowess_curve[:,1])#,kind="nearest-up",fill_value="extrapolate")
    if plot:
        plt.plot(sizes, thresholds,"b--",lw=2)
        plt.plot(sizes,get_min_snr(sizes),"r-",lw=2)
        plt.xlabel("n_samples")
        plt.ylabel("SNR threshold")
        plt.show()
    return get_min_snr

def jenks_binarization(exprs, empirical_snr,min_n_samples,verbose = True,
                      plot=True, plot_SNR_thr= 3.0, show_fits = []):
    t0= time()
    if verbose:
        print("\tJenks method is chosen ...")
    binarized_expressions = {"UP" : {}, "DOWN" : {}}
    N = exprs.shape[1]
    n_bins = max(20,int(N/10))
    n_ambiguous = 0
    stats = {}
    for i, (gene, row) in enumerate(exprs.iterrows()):
        values = row.values
        hist_range = values.min(), values.max()
        up_color, down_color = "grey", "grey"
        
        ### special treatment of zero expressions, try zero vs non-zero before Jenks breaks
        ### TBD
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
        
        SNR, size, e_pval = np.nan,np.nan,np.nan
        if min(len(down_group),len(up_group)) >= min_n_samples:
            # calculate SNR, size and empirical p-val
            SNR = calc_SNR(up_group,down_group)
            size = min(len(down_group),len(up_group)) 
            # SNR p-val depends on biclister size 
            e_snr_dist = empirical_snr[size - min_n_samples,]
            e_pval =  (len(e_snr_dist[e_snr_dist>=abs(SNR)])+1)/(empirical_snr.shape[1]+1)
            
            # in case of insignificant difference 
            # the bigger half is treated as signal too
            if abs(len(up_group)-len(down_group)) <= min_n_samples:
                n_ambiguous +=1 

            if len(down_group)-len(up_group) >= -min_n_samples:
                # up-regulated group is smaller than down-regulated
                binarized_expressions["UP"][gene] = pos_mask.astype(int)
                up_color='red'

            if len(up_group)-len(down_group) >= -min_n_samples:
                binarized_expressions["DOWN"][gene] = neg_mask.astype(int)
                down_color='blue'
        stats[gene] = {"SNR":SNR,"size": size,"pval":e_pval}
   
        
        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:",i)
                
        if gene in show_fits or SNR > plot_SNR_thr:
            print("Gene {}: SNR={:.2f}, pos={}, neg={}".format(gene, SNR, len(up_group), len(down_group)))
            plt.hist(down_group, bins=n_bins, alpha=0.5, color=down_color,range=hist_range)
            plt.hist(up_group, bins=n_bins, alpha=0.5, color=up_color,range=hist_range)
            plt.show()
    
    stats = pd.DataFrame.from_dict(stats).T
    stats = stats.dropna(subset = ["pval"])
    binarized_expressions["UP"] = pd.DataFrame.from_dict(binarized_expressions["UP"])
    #binarized_expressions["UP"] =  binarized_expressions["UP"].loc[:,stats.index]
    binarized_expressions["DOWN"] = pd.DataFrame.from_dict(binarized_expressions["DOWN"])
    #binarized_expressions["DOWN"] = binarized_expressions["DOWN"].loc[:,stats.index]
    
    
    # logging
    if verbose:
        print("\tJenks binarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))
        #print("\t\tup-regulated features:", binarized_expressions["UP"].shape[1])
        #print("\t\tdown-regulated features:", binarized_expressions["DOWN"].shape[1])
        #print("\t\tambiguous features:", n_ambiguous )
    
    return binarized_expressions, stats

def select_pos_neg(row, empirical_snr, min_n_samples, seed=42):
    """ find 'higher' (positive), and 'lower' (negative) signal in vals. 
        vals are found with GM binarization
    """

    with warnings.catch_warnings(): # this is to ignore convergence warnings 
        warnings.simplefilter('ignore')
        
        row2d = row[:, np.newaxis]  # adding mock axis for GM interface
        np.random.seed(seed=seed)
        model = GaussianMixture(
            n_components=2, init_params="kmeans",
            max_iter=len(row), n_init = 1, 
            covariance_type = "spherical",
            random_state = seed
        ).fit(row2d)
        
        labels = model.predict(row2d) 
        
        sig_snr = calc_SNR(row[labels==0], row[labels==1])
        
    
    # let labels == 1 be the bigger half
    if np.sum(labels == 0) > np.sum(labels == 1): 
        labels = 1 - labels
    n0 = np.sum(labels == 0)
    n1 = np.sum(labels == 1)
    assert n0 + n1 == len(row)

    signal_pretendents = []
    e_pval = np.nan
    size = np.nan
    if n0 >= min_n_samples:
        size = n0
        # SNR p-val depends on biclister size 
        e_snr_dist = empirical_snr[size - min_n_samples,]
        e_pval =  (len(e_snr_dist[e_snr_dist>=abs(sig_snr)])+1)/(empirical_snr.shape[1]+1)
        # signal (bicluster) should be big enough
        signal_pretendents.append(labels==0)
        if n1 - n0 < min_n_samples:
            # in case of insignificant difference 
            # the bigger half is treated as signal too
            signal_pretendents.append(labels==1) 
            #stat['n_inexplicit'] += 1 
            
    mask_pos = np.zeros_like(labels, bool)
    mask_neg = np.zeros_like(labels, bool)
    for mask in signal_pretendents:
        sig_snr = calc_SNR(row[mask], row[~mask])
        if sig_snr > 0:
            mask_pos |= mask
        else: 
            mask_neg |= mask

    return mask_pos, mask_neg, abs(sig_snr), size, e_pval, model.converged_


def GM_binarization(exprs,empirical_snr, min_n_samples, verbose=True, plot=True, plot_SNR_thr=2, show_fits=[],seed=1):
    t0 = time()
    if verbose:
        print("\tGMM method is chosen ...")
    N = exprs.shape[1]
    n_bins = max(20,int(N/10))
    stats = {}
    e_pvals = {}
    mask_pos, mask_neg = [], []
    for i, (gene, row) in enumerate(exprs.iterrows()):
        row = row.values
        row_mask_pos, row_mask_neg, snr, size, e_pval, is_converged = select_pos_neg(row, empirical_snr, min_n_samples, seed=seed)
        mask_pos.append(row_mask_pos.astype(int))
        mask_neg.append(row_mask_neg.astype(int))
        stats[gene] = {"pval":e_pval,"SNR":snr,"size":size,"convergence":is_converged}

        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:",i)
        s_pos = len(row[row_mask_pos])
        s_neg = len(row[row_mask_neg])
        if plot and ((abs(snr) > plot_SNR_thr and max(s_pos,s_neg)>0) or gene in show_fits):
            print("Gene %s: SNR=%s, pos=%s, neg=%s"%(gene, round(snr,2), s_pos, s_neg))
            row_mask_neutral = (~row_mask_pos) & (~row_mask_neg)
            if (e_pval<0.05 and max(s_pos,s_neg)>0):
                plt.hist(row[row_mask_neutral], bins=n_bins, alpha=0.5, color='lightgrey')
                plt.hist(row[row_mask_neg], bins=n_bins, alpha=0.5, color='blue')
                plt.hist(row[row_mask_pos], bins=n_bins, alpha=0.5, color='red')
            else:
                plt.hist(row, bins=n_bins, alpha=0.5, color='grey')
            plt.show()
 
    def _remove_empty_rows(df):
        # thx https://stackoverflow.com/a/22650162/7647325
        return df.loc[~(df==0).all(axis=1)]
    
    df_p = _remove_empty_rows(pd.DataFrame(mask_pos, index=exprs.index)).T
    df_n = _remove_empty_rows(pd.DataFrame(mask_neg, index=exprs.index)).T
    
    
    stats = pd.DataFrame.from_dict(stats).T
    stats = stats.loc[list(set(df_p.columns).union(set(df_n.columns))),:]

    # logging
    if verbose:
        print("\tGMM binarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))

    return {"UP":df_p, "DOWN":df_n}, stats


def binarize(binarized_fname_prefix, exprs=None, method='GMM',
             save = True, load = False,
             min_n_samples = 10, pval = 0.001,
             plot_all = True, plot_SNR_thr= np.inf,show_fits = [],
             verbose= True,seed=random.randint(0,100000)):
    '''
       binarized_fname_prefix is a basename of binarized data file;
       exprs is a dataframe with normalized features to be binarized.
    '''
    t0 = time()
    if load:
        # load from file binarized genes
        binarized_data = {}

        for d in ["UP","DOWN"]:
            fname = binarized_fname_prefix +"."+method+".binarized_"+d+".tsv"
            print("Load binarized features from",fname,file = sys.stdout)
            df = pd.read_csv(fname, sep ="\t",index_col=0)
            binarized_data[d] = df
        print("",file = sys.stdout)
        
        # load stats 
        fname = binarized_fname_prefix +"."+method+".binarization_stats.tsv"
        print("Load statistics from",fname,file = sys.stdout)
        stats = pd.read_csv(fname, sep ="\t",index_col=0)
        print("",file = sys.stdout)
        
    elif exprs is not None:
        # compute from expressions
        start_time = time()
        if verbose:
            print("\nBinarization started ....\n")

        t0 = time()
        sizes,thresholds,empirical_snr = random_splits(exprs, min_n_samples,
                                                       pval = pval,
                                                       seed =seed,verbose = verbose)
 
        size_snr_trend = get_trend(sizes, thresholds, plot= False)
        if verbose:
            print("\tSNR thresholds for individual features computed in {:.2f} seconds".format(time()-t0))

        if method=="Jenks":
            binarized_data, stats = jenks_binarization(exprs, empirical_snr,min_n_samples,
                                                       verbose = verbose, plot=plot_all,
                                                       plot_SNR_thr=plot_SNR_thr, show_fits = show_fits)
        elif method=="GMM":
            binarized_data, stats = GM_binarization(exprs, empirical_snr,min_n_samples,verbose = verbose,
                                             plot=plot_all, plot_SNR_thr= plot_SNR_thr, 
                                             show_fits = show_fits, seed = seed )
        else:
            print("Method must be either 'GMM' or 'Jenks'.",file=sys.stderr)
            return
        
        if verbose:
                print("\tbinarization runtime: {:.2f} s".format(time()-start_time ),file = sys.stdout)
        
        # keep only features with SNR above the thresholds
        stats["SNR_threshold"] = stats["size"].apply(lambda x: size_snr_trend(x))
        
        # z-scores for SNR values
        emean = empirical_snr.mean(axis=1)
        estd = empirical_snr.std(axis=1)
        
        def calc_SNR_zscore(row,m = emean, s = estd):
            SNR = row["SNR"]
            size = row["size"]
            ndx = np.where(sizes == size)[0][0]
            z = (SNR-emean[ndx])/estd[ndx]
            return z
        
        stats["z"] = stats.apply(lambda row: calc_SNR_zscore(row),axis=1)
        #binarized_data, stats = keep_binarized(binarized_data,stats,FDR=False,verbose=True)
        bh_res, adj_pval = fdrcorrection(stats["pval"].values,alpha=0.01)
        stats["qval"] = adj_pval
        stats = stats.sort_values(by="z",ascending = False)
        
        ### keep features passed binarization
        passed = set(stats.loc[stats["SNR"]>stats["SNR_threshold"],:].index)
        if verbose:
            print("\t%s features passed binarization "%len(passed),file = sys.stdout)
        for d in ["UP", "DOWN"]:
            df = binarized_data[d]
            df = df.loc[:,sorted(list(passed.intersection(df.columns)))]
            binarized_data[d] = df
            if verbose:
                print("\t\t%s-regulated features:\t%s"%(d,df.shape[1]),file = sys.stdout)
        ambiguous_features = set(binarized_data["UP"].columns).intersection(set(binarized_data["DOWN"].columns))
        if verbose:
            print("\t\tambiguous features:\t"+str(len(ambiguous_features)),file = sys.stdout)
        
        
        if save:
            # save to file binarized data
            sample_names = exprs.columns
            
            fname = binarized_fname_prefix +"."+method+".binarization_stats.tsv"
            print("Statistics is saved to",fname,file = sys.stdout)
            stats.to_csv(fname, sep ="\t")
            for d in ["UP","DOWN"]:
                df = binarized_data[d]
                df.index = sample_names
                fname = binarized_fname_prefix +"."+method+".binarized_"+d+".tsv"
                print("Binarized gene expressions are saved to",fname,file = sys.stdout)
                df.to_csv(fname, sep ="\t")
    else:
        print("Provide either raw or binarized data.", file=sys.stderr)
        return None
    
    if plot_all:
        # plot null distribution of SNR(size) - for shuffled genes n_permutation times for each size
        e_stats = []
        for i in range(empirical_snr.shape[0]):
            for j in range(empirical_snr.shape[1]):
                e_stats.append({"size":min_n_samples+i,"SNR":empirical_snr[i,j]})
        e_stats = pd.DataFrame.from_records(e_stats)
        tmp  = plt.figure(figsize=(20,10))
        tmp = plt.scatter(e_stats["size"],e_stats["SNR"],alpha = 1, color = "grey")
        
        # plot binarization results for real genes
        passed = stats.loc[stats["SNR"]> stats["SNR_threshold"],:]
        failed = stats.loc[stats["SNR"]<= stats["SNR_threshold"],:]
        tmp = plt.scatter(failed["size"],failed["SNR"],alpha = 0.5, color = "black")
        tmp = plt.scatter(passed["size"],passed["SNR"],alpha = 0.5, color = "red")
        
        # plot cutoff 
        sizes  = sorted(list(set(stats["size"].values)))
        tmp = plt.plot(sizes,[size_snr_trend(x) for x in sizes], c="yellow",lw=2)
        tmp = plt.xlabel("n_samples")
        tmp = plt.ylabel("SNR threshold")
        tmp = plt.show()
        

    return binarized_data, stats, empirical_snr

######## Clustering #########

def run_WGCNA(fname,p1=10,p2=10,verbose = False):
    # run Rscript
    if verbose:
        t0 = time()
        print("Running WGCNA for", fname, "...")
    process = subprocess.Popen(['Rscript','run_WGCNA.R', str(p1),str(p2),fname],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    module_file = fname.replace(".tsv",".modules.tsv")#stdout.rstrip()
    #print(module_file)
    modules_df = pd.read_csv(module_file,sep = "\t",index_col=0)
    if verbose:
        print("\tWGCNA runtime: modules detected in {:.2f} s.".format(time()-t0))
    
    # read WGCNA output
    modules = []
    not_clustered = []
    module_dict = modules_df.T.to_dict()
    for i in module_dict.keys():
        genes =  module_dict[i]["genes"].strip().split()
        if i == 0:
            not_clustered = genes
        else:
            modules.append(genes)
    if verbose:
        print("\t{} modules and {} not clustered genes".format(len(modules),len(not_clustered)))
    return (modules,not_clustered)

def make_biclsuter_jenks(exprs, bin_exprs,min_SNR =0,min_n_samples=-1,
                          verbose= True, plot=True):
    
    bicluster = {"n_genes":0}
    genes = exprs.index.to_list()
    ones_per_sample = bin_exprs.sum(axis=1)
    
    breaks = jenkspy.jenks_breaks(ones_per_sample, nb_class=2)
    border = breaks[1] 
    
    if plot and len(genes)>2:
        counts, pos, obj = plt.hist(ones_per_sample,bins=min(20,max(ones_per_sample)))
        tmp = plt.vlines(breaks[1],0,max(counts)*1.1,colors="red")
        tmp = plt.text(breaks[1],max(counts),str(breaks[1]))
        tmp = plt.show()
    
    samples = ones_per_sample[ones_per_sample>border].index.to_list()
    if len(samples) < min_n_samples:
        # not enough samples -> failed bicluster 
        return bicluster, genes   
    
    bg_samples = ones_per_sample[ones_per_sample<=border].index
    bg = exprs.loc[:, bg_samples]
    bic = exprs.loc[:,samples]
    SNR = (bic.mean(axis=1) - bg.mean(axis=1))/(bg.std(axis=1) + bic.std(axis=1))
    SNR = SNR.abs()
    #print(SNR)
    
    excluded_genes  = []
    if min_SNR > 0:
        excluded_genes = SNR[SNR<min_SNR].index.to_list()
        if len(excluded_genes)>0:
            SNR = SNR[SNR>=min_SNR]
            genes = SNR.index.to_list()
    
    if len(genes)<2:
        # not enough genes -> failed bicluster 
        return bicluster, genes + excluded_genes 
    
    avgSNR = SNR.mean()
    if verbose and len(genes)>2:
        print(SNR.shape[0],"x",len(samples),"avgSNR:",round(avgSNR,2))
        print("droped:(%s)"%len(excluded_genes)," ".join(excluded_genes))
        print("genes:(%s)"%len(genes)," ".join(genes))
        print("samples:(%s)"%len(samples)," ".join(samples))
    bicluster["samples"] = samples
    bicluster["n_samples"] = len(samples)
    bicluster["genes"] = genes
    bicluster["n_genes"] = len(genes)
    bicluster["avgSNR"] = avgSNR
    return bicluster, excluded_genes



def write_bic_table(bics_dict_or_df, results_file_name,to_str=True):
    bics = bics_dict_or_df.copy()
    if len(bics) ==0:
        print("No biclusters found",file=sys.stderr)
    else:
        if not type(bics) == type(pd.DataFrame()):
            bics = pd.DataFrame.from_dict(bics)
        if to_str:
            bics["genes"] = bics["genes"].apply(lambda x:" ".join(map(str,sorted(x))))
            bics["samples"] = bics["samples"].apply(lambda x:" ".join(map(str,sorted(x))))
        bics = bics.sort_values(by=["avgSNR","n_genes","n_samples"], ascending = False)
        bics.index = range(0,bics.shape[0])
        bics.index.name = "id"
        cols =  bics.columns.values
        first_cols = ["avgSNR","n_genes","n_samples","direction","genes","samples"]
        bics = bics.loc[:,first_cols+sorted(list(set(cols).difference(first_cols)))]
    bics.to_csv(results_file_name ,sep = "\t")
    
def modules2biclsuters_jenks(clustering_results,exprs,binarized_expressions,
                            min_SNR = 0.0,min_n_samples=20,
                            snr_pval=0.01,
                            bin_method = "GMM", clust_method = "WGCNA",
                            result_file_name = False,
                            directions=["UP","DOWN"],
                            plot=False,verbose = False):
    biclusters = {}
    not_clustered = {}
    i = 0
    for d in directions:
        biclusters[d] = {}
        exprs_bin = binarized_expressions[d]
        modules  = clustering_results[d][0].copy()
        not_clustered[d] = clustering_results[d][1].copy()

        for genes in modules:
            bicluster,excluded_genes = make_biclsuter_jenks(exprs.loc[genes,:], exprs_bin.loc[:,genes],
                                                                min_SNR = min_SNR,min_n_samples=min_n_samples,
                                                                plot=plot,verbose = verbose)
            not_clustered[d]+=excluded_genes

            # while any genes deleted, try to update samples of the bicluster
            while len(excluded_genes) >0 and bicluster["n_genes"]>1:
                genes = bicluster["genes"]
                bicluster,excluded_genes = make_biclsuter_jenks(exprs.loc[genes,:], exprs_bin.loc[:,genes],
                                                                min_SNR = min_SNR,min_n_samples=min_n_samples,
                                                                plot=plot,verbose = verbose)
                not_clustered[d]+=excluded_genes
                if bicluster["n_genes"] > 0:
                    genes = bicluster["genes"]

            if bicluster["n_genes"] > 1:
                bicluster["direction"] = d
                
                biclusters[d][i] = bicluster
                i+=1
            elif bicluster["n_genes"] == 1:
                not_clustered[d]+= bicluster["genes"]

        biclusters[d] = pd.DataFrame.from_dict(biclusters[d]).T
        biclusters[d] = biclusters[d].loc[:,["avgSNR","direction","n_genes","n_samples","genes","samples"]]
        biclusters[d] = biclusters[d].sort_values(by = ["avgSNR","n_genes","n_samples"],ascending = [False,False, False])
        biclusters[d].index = range(0,biclusters[d].shape[0])

        print("{}: {} features clustered into {} modules, {} - not clustered.".format(d,biclusters[d]["n_genes"].sum(),
                                                                                      biclusters[d].shape[0],
                                                                                      len(not_clustered[d])))
        
    
    biclusters = pd.concat([biclusters["UP"],biclusters["DOWN"]],axis =0)
    biclusters = biclusters.sort_values(by=["avgSNR"],ascending=[False])
    biclusters.index = range(0,biclusters.shape[0])
        
    if result_file_name:
        suffix  = ".bin="+bin_method+",clust="+clust_method
        write_bic_table(biclusters, result_file_name+suffix+".biclusters.tsv")
        
    amigouos_genes = set(not_clustered["UP"]).intersection(set(not_clustered["DOWN"]))
    df_up = binarized_expressions["UP"].T
    df_down = binarized_expressions["DOWN"].T
    nc_features = pd.concat([df_up.loc[not_clustered["UP"],:],df_down.loc[set(not_clustered["DOWN"]).difference(amigouos_genes),:] ])
    print("Unique not clustered: {}".format(nc_features.shape[0]),file=sys.stdout)
    
    if result_file_name:
        nc_features.to_csv(result_file_name+suffix+".not_clustered.tsv",sep = "\t")
        
    return biclusters, nc_features


#### K-means based biclustering

def run_2means(bic_genes,exprs,min_n_samples=10,seed=0):
    # identify identify bicluster and backgound groups using 2-means
    e = exprs[bic_genes,:].T
    
    labels = KMeans(n_clusters=2, random_state=seed,n_init=5).fit(e).labels_
    #labels = GaussianMixture(n_components=2, init_params="kmeans",
    #        max_iter=e.shape[0], n_init = 1, 
    #        covariance_type = "spherical",
    #        random_state = seed).fit_predict(e)
    ndx0 = np.where(labels == 0)[0]
    ndx1 = np.where(labels == 1)[0]
    if min(len(ndx1),len(ndx0))< min_n_samples:
        return {}
    if len(ndx1) > len(ndx0):
        samples = ndx0
    else: 
        samples = ndx1

    bic = {"gene_ids":set(bic_genes),"n_genes":len(bic_genes),
           "sample_ids":set(samples),"n_samples":len(samples)}
    return bic

def add_SNR_to_biclusters(bics, exprs, exprs_data):
    # calculates SNR for each bicluster in a list
    N, exprs_sums, exprs_sq_sums = exprs_data
    for i in range(len(bics)):
        gene_ids = list(bics[i]['gene_ids'])
        sample_ids = list(bics[i]['sample_ids'])
        # calcluate SNR 
        avgSNR = calc_bic_SNR(gene_ids, sample_ids, exprs, N, exprs_sums, exprs_sq_sums)
        bics[i]["avgSNR"] = abs(avgSNR) 
        # direction
        if avgSNR >0:
            bics[i]["direction"] = "UP"
        else:
            bics[i]["direction"] = "DOWN"
    return bics

def modules2biclusters(modules, genes2ids, exprs,
                       min_n_samples=10, min_n_genes=2,seed=0,verbose = True):
    # Identify optimal sample set for each module: split samples into two sets in a subspace of each module
    t0 = time()
    bics = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0
    
    for mid in range(0,len(modules)):
        gene_ids = modules[mid]
        if len(gene_ids) >= min_n_genes: 
            gene_ids = [genes2ids[x] for x in gene_ids]
            # define bicluster and background
            try:
                bic = run_2means(gene_ids,exprs,min_n_samples=min_n_samples,seed=seed)
            except:
                print(gene_ids)
                print(exprs)
                bic = run_2means(gene_ids,exprs,min_n_samples=min_n_samples,seed=seed)
            if len(bic)>0:
                bic["id"] = i
                bics[i] = bic
                i+=1
      
    if verbose:
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %(len(modules),round(time()-t0,2)))
        print("Passed biclusters (>=%s genes, >= samples %s): %s"%(min_n_genes,min_n_samples,i-1), file = sys.stdout)
        print()
    
    return bics

def make_biclusters(clustering_results,binarized_expressions,exprs,
                    min_n_samples=10, min_n_genes=2,
                    seed = 42,
                    save=False,out_dir="./",basename="",bin_method="",clust_method="",
                    cluster_binary = False):
    filtered_bics = []
    for d in ["UP","DOWN"]:
        genes = binarized_expressions[d].columns.values
        genes2ids = dict(zip(genes,range(0,len(genes))))

        if cluster_binary:
            exprs_np = binarized_expressions[d].T # binarized expressions
        else:
            exprs_np = exprs.loc[genes,:] # z-scoes
        ints2g_names = exprs_np.index.values
        ints2s_names = exprs_np.columns.values
        exprs_np = exprs_np.values    
        modules = clustering_results[d][0]

        bics = modules2biclusters(modules, genes2ids, exprs_np,
                                min_n_samples=min_n_samples, min_n_genes=2,
                                verbose = False,seed=seed)

        # make exprs_data for fast SNR calculations 
        exprs_np = exprs.loc[genes,:].values
        exprs_sums = exprs_np.sum(axis=1)
        exprs_sq_sums = np.square(exprs_np).sum(axis=1)
        N = exprs.shape[1]
        exprs_data = N, exprs_sums, exprs_sq_sums
        # add SNR and direction to biclsuters
        bics = add_SNR_to_biclusters(bics, exprs.loc[genes,:].values, exprs_data)
        # rename genes and samples
        for i in range(len(bics)):
            bics[i]["genes"] = set([ints2g_names[x] for x in bics[i]["gene_ids"]])
            bics[i]["samples"] = set([ints2s_names[x] for x in bics[i]["sample_ids"]])

        bics = pd.DataFrame.from_dict(bics).T
        bics.index = bics["id"]
        print(bics.shape)
        filtered_bics.append(bics)
    filtered_bics = pd.concat(filtered_bics,axis =0)
    filtered_bics = filtered_bics.sort_values(by=["avgSNR"],ascending=[False])
    filtered_bics.drop("id",inplace=True,axis=1)
    filtered_bics.index = range(0,filtered_bics.shape[0])

    ### TBD - merge similar up- and down-regulted

    if save:
        from method2 import write_bic_table
        suffix  = ".bin="+bin_method+",clust="+clust_method
        write_bic_table(filtered_bics, out_dir+basename+suffix+".biclusters.tsv")
        print(out_dir+basename+suffix+".biclusters.tsv")
    return filtered_bics   

def calc_bic_SNR(genes, samples, exprs, N, exprs_sums,exprs_sq_sums):
    bic = exprs[genes,:][:,samples]
    bic_sums = bic.sum(axis=1)
    bic_sq_sums = np.square(bic).sum(axis=1)

    bg_counts = N - len(samples)
    bg_sums = exprs_sums[genes]-bic_sums
    bg_sq_sums = exprs_sq_sums[genes]-bic_sq_sums
    
    bic_mean, bic_std = calc_mean_std_by_powers((len(samples),bic_sums,bic_sq_sums))
    bg_mean, bg_std = calc_mean_std_by_powers((bg_counts,bg_sums,bg_sq_sums))
    
    return  np.mean((bic_mean - bg_mean)/ (bic_std + bg_std))

def calc_mean_std_by_powers(powers):
    count, val_sum, sum_sq = powers

    mean = val_sum / count  # what if count == 0?
    std = np.sqrt((sum_sq / count) - mean*mean)
    return mean, std


def run_Louvain(similarity,verbose = True):
    t0 = time()
    gene_names = similarity.index.values
    louvain = Louvain(modularity = 'newman') # 'potts','dugue', 'newman'
    sparse_matrix = csr_matrix(similarity)
    labels = louvain.fit_transform(sparse_matrix)
    Q = modularity(sparse_matrix, labels)
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
        print("\t","modules:",len(modules),"not_clustered:",len(not_clustered))
        print("\t\tmodularity:",round(Q,3), "runtime:",round(time()-t0,2),"s")
    return [modules, not_clustered]


def run_MCL(similarity, inflations = list(np.arange(11,30)/10)+[3,3.5,4,5],verbose = True):
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

def get_similarity_jaccard(df,qval=0.01,J=0.33):
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

def get_similarity_corr(df,r=0):
    corr = df.corr()
    corr = corr[corr>r] 
    corr = corr.fillna(0)
    return corr