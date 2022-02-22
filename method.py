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
import jenkspy

from sklearn.cluster import KMeans
from method2 import calc_bic_SNR , identify_opt_sample_set



import matplotlib.pyplot as plt
import statsmodels.api as sm



def calc_SNR(ar1, ar2):
    return (np.mean(ar1) - np.mean(ar2)) / (np.std(ar1) + np.std(ar2))

######### Binarization #########

def rand_norm_splits(N, min_n_samples, snr_pval = 0.01,seed = 42,verbose = True):
    # empirical distribuition of SNR depending on the bicluster size 
    # generates normal samples of matching size and split them into bicluster and background
    min_n_perm = int(5*1.0/snr_pval)
    n_perm = max(min_n_perm,int(100000/N))
    sizes = np.arange(min_n_samples,int(N/2)+min_n_samples+1)#.shape
    if verbose:
        print("\tGenerate empirical distribuition of SNR depending on the bicluster size ...")
        print("\t\ttotal samples: %s, number of samples in bicluster: %s - %s, n_permutations: %s"%(N,sizes[0],sizes[-1],n_perm))
    snr_thresholds = np.zeros(sizes.shape)
    np.random.seed(seed=seed)
    for s in sizes:
        snrs = np.zeros(n_perm)
        for i in range(0,n_perm):
            x = np.random.normal(size = N)
            x.sort()
            snrs[i] = calc_SNR(x[s:], x[:s]) #(x[s:].mean()-x[:s].mean())/(x[s:].std()+x[:s].std())
        snr_thresholds[s-min_n_samples]=np.quantile(snrs,q=1-0.05)
    return sizes, snr_thresholds

def get_trend(sizes, thresholds, plot= True):
    # smoothens the trend and retunrs a function min_SNR(size; p-val. cutoff)
    lowess = sm.nonparametric.lowess
    lowess_curve = lowess(sizes, thresholds,frac=0.25,return_sorted=True,is_sorted=False)
    get_min_snr = interp1d(lowess_curve[:,1],lowess_curve[:,0],kind="nearest",fill_value="extrapolate")
    if plot:
        plt.plot(sizes, thresholds,"b--",lw=2)
        plt.plot(sizes,get_min_snr(sizes),"r-",lw=2)
        plt.xlabel("n_samples")
        plt.ylabel("SNR threshold")
        plt.show()
    return get_min_snr

def jenks_binarization(exprs, get_min_snr,min_n_samples,verbose = True,
                      plot=True, plot_SNR_thr= 3.0, show_fits = []):
    t0= time()
    if verbose:
        print("\tJenks method is chosen ...")
    binarized_expressions = {"UP" : {}, "DOWN" : {}}
    N = exprs.shape[1]
    n_bins = max(20,int(N/10))
    n_ambiguous = 0
    snrs = []
    sizes = []
    genes = []
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
        
        if min(len(down_group),len(up_group)) >= min_n_samples:

            # calculate SNR 
            SNR = calc_SNR(up_group,down_group)
            size = min(len(down_group),len(up_group))
            # define bicluster and background if SNR is signif. higer than random
            if SNR >= get_min_snr(size):
                snrs.append(SNR)
                sizes.append(size)
                genes.append(gene)

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
                
        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:",i)
                
        if gene in show_fits or SNR > plot_SNR_thr:
            print("Gene {}: SNR={:.2f}, pos={}, neg={}".format(gene, SNR, len(up_group), len(down_group)))
            plt.hist(down_group, bins=n_bins, alpha=0.5, color=down_color,range=hist_range)
            plt.hist(up_group, bins=n_bins, alpha=0.5, color=up_color,range=hist_range)
            plt.show()
    
    binarized_expressions["UP"] = pd.DataFrame.from_dict(binarized_expressions["UP"])
    binarized_expressions["DOWN"] = pd.DataFrame.from_dict(binarized_expressions["DOWN"])
    
    # logging
    if verbose:
        print("\tJenks binarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))
        print("\t\tup-regulated features:", binarized_expressions["UP"].shape[1])
        print("\t\tdown-regulated features:", binarized_expressions["DOWN"].shape[1])
        print("\t\tambiguous features:", n_ambiguous )
    stats = pd.DataFrame.from_records({"SNR":snrs,"size":sizes}, index = genes)
    return binarized_expressions, stats

def select_pos_neg(row, get_min_snr, min_n_samples, min_diff_samples, stat,seed=42):
    """ find 'higher' (positive), and 'lower' (negative) signal in vals. 
        vals are found with GM binarization
    """

    with warnings.catch_warnings(): # this is to ignore convergence warnings 
        warnings.simplefilter('ignore')
        
        row2d = row[:, np.newaxis]  # adding mock axis for GM interface
        np.random.seed(seed=seed)
        labels = GaussianMixture(
            n_components=2, init_params="kmeans",
            max_iter=300, n_init = 1, 
            covariance_type = "spherical",
            random_state = seed
        ).fit(row2d).predict(row2d) # Bayesian
        
        sig_snr = calc_SNR(row[labels==0], row[labels==1])
        stat['SNRs'].append(sig_snr) 
    
    # let labels == 1 be the bigger half
    if np.sum(labels == 0) > np.sum(labels == 1): 
        labels = 1 - labels
    n0 = np.sum(labels == 0)
    n1 = np.sum(labels == 1)
    assert n0 + n1 == len(row)

    # min_SNR threshold depends on biclister size
    min_SNR = get_min_snr(n0) 
    signal_pretendents = []
    if min_n_samples < n0: 
        # signal (bicluster) should be big enough
        signal_pretendents.append(labels==0)
        if n1 - n0 < min_diff_samples:
            # in case of insignificant difference 
            # the bigger half is treated as signal too
            signal_pretendents.append(labels==1) 
            if abs(sig_snr) > min_SNR: 
                stat['n_inexplicit'] += 1 
            
    mask_pos = np.zeros_like(labels, bool)
    mask_neg = np.zeros_like(labels, bool)
    for mask in signal_pretendents:
        sig_snr = calc_SNR(row[mask], row[~mask])
        if abs(sig_snr) > min_SNR:
            if sig_snr > 0:
                mask_pos |= mask
            else: 
                mask_neg |= mask

    return mask_pos, mask_neg

def GM_binarization(exprs, get_min_snr, min_n_samples, verbose=True, plot=True, plot_SNR_thr=2, show_fits=[],seed=1):
    t0 = time()
    N = exprs.shape[1]
    n_bins = max(20,int(N/10))
    stat = {
        'SNRs': [],
        'n_inexplicit': 0,
    }

    mask_pos, mask_neg = [], []
    for i, (gene, row) in enumerate(exprs.iterrows()):
        row = row.values
        row_mask_pos, row_mask_neg = select_pos_neg(row, get_min_snr, min_n_samples, min_n_samples, stat,seed=seed)
        mask_pos.append(row_mask_pos.astype(int))
        mask_neg.append(row_mask_neg.astype(int))

        # logging
        if verbose:
            if i % 1000 == 0:
                print("\t\tgenes processed:",i)
        SNR = stat['SNRs'][-1]
        s_pos = len(row[row_mask_pos])
        s_neg = len(row[row_mask_neg])
        if plot and ((abs(SNR) > plot_SNR_thr and max(s_pos,s_neg)>0) or gene in show_fits):
            print("Gene %s: SNR=%s, pos=%s, neg=%s"%(gene, round(SNR,2), s_pos, s_neg))
            row_mask_neutral = (~row_mask_pos) & (~row_mask_neg)

            plt.hist(row[row_mask_neutral], bins=n_bins, alpha=0.5, color='grey')
            plt.hist(row[row_mask_neg], bins=n_bins, alpha=0.5, color='blue')
            plt.hist(row[row_mask_pos], bins=n_bins, alpha=0.5, color='red')
            plt.show()
 
    def _remove_empty_rows(df):
        # thx https://stackoverflow.com/a/22650162/7647325
        return df.loc[~(df==0).all(axis=1)]
    df_p = _remove_empty_rows(pd.DataFrame(mask_pos, index=exprs.index)).T
    df_n = _remove_empty_rows(pd.DataFrame(mask_neg, index=exprs.index)).T

    # logging
    if verbose:
        print("\tGMM binarization for {} features completed in {:.2f} s".format(len(exprs),time()-t0))
        print("\tup-regulated features:", df_p.shape[1])
        print("\tdown-regulated features:", df_n.shape[1])
        print("\tambiguous features:", stat['n_inexplicit'])

    return {"UP":df_p, "DOWN":df_n}

def binarize(binarized_fname_prefix, exprs=None, method='GMM',
             save = True, load = False,
             min_n_samples = 10, snr_pval = 0.01,
             plot_all = True, plot_SNR_thr= np.inf,show_fits = [],
             verbose= True,seed=42):
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
        
    elif exprs is not None:
        # compute from expressions
        start_time = time()
        if verbose:
            print("\nBinarization started ....\n")

        t0 = time()
        sizes,thresholds = rand_norm_splits(exprs.shape[1],min_n_samples, snr_pval = snr_pval,seed=seed)
        get_min_snr = get_trend(sizes,thresholds, plot = plot_all)
        if verbose:
            print("\tSNR thresholds for individual features computed in {:.2f} seconds".format(time()-t0))

        if method=="Jenks":
            binarized_data, stats = jenks_binarization(exprs, get_min_snr,min_n_samples,verbose = verbose,
                                                  plot=plot_all, plot_SNR_thr=plot_SNR_thr, show_fits = show_fits)
        elif method=="GMM":
            binarized_data = GM_binarization(exprs,get_min_snr,min_n_samples,verbose = verbose,
                                                    plot=plot_all, plot_SNR_thr= plot_SNR_thr, show_fits = show_fits,
                                                    seed = seed)
        else:
            print("Method must be either 'GMM' or 'Jenks'.",file=sys.stderr)
            return
        
        if verbose:
                print("\tbinarization runtime: {:.2f} s".format(time()-start_time ),file = sys.stdout)
        if save:
            # save to file binarized data
            sample_names = exprs.columns
            for d in ["UP","DOWN"]:
                df = binarized_data[d]
                df.index = sample_names
                fname = binarized_fname_prefix +"."+method+".binarized_"+d+".tsv"
                print("Binarized gene expressions are saved to",fname,file = sys.stdout)
                df.to_csv(fname, sep ="\t")
    else:
        print("Provide either raw or binarized data.", file=sys.stderr)
        return None
    return binarized_data

######## Clustering #########

def run_WGCNA(fname,p1=10,p2=10,verbose = False):
    # run Rscript
    if verbose:
        t0 = time()
        print("Running WGCNA for", fname, "...")
    process = subprocess.Popen(['Rscript','run_WGCNA.R', '10','10',fname],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    module_file = stdout.rstrip()
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
    
    
from method import make_biclsuter_jenks, write_bic_table
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

def identify_opt_sample_set(bic_genes,exprs_bin,exprs_data,min_n_samples=8):
    # identify optimal samples set given gene set
    N, exprs_sums, exprs_sq_sums = exprs_data
    e = exprs_bin[bic_genes,:]
    
    labels = KMeans(n_clusters=2, random_state=0).fit(e.T).labels_
    ndx0 = np.where(labels == 0)[0]
    ndx1 = np.where(labels == 1)[0]
    if min(len(ndx1),len(ndx0))< min_n_samples:
        return {"avgSNR":-1}
    if len(ndx1) > len(ndx0):
        samples = ndx0
    else: 
        samples = ndx1
    
    avgSNR = calc_bic_SNR(bic_genes, samples, exprs_bin, N, exprs_sums, exprs_sq_sums)
    if avgSNR >0:
        direction = "UP"
    else:
        direction = "DOWN"

    if len(samples)>=min_n_samples: # len(samples)<N*0.5*1.1 - allow bicluster to be a little bit bigger than N/2
        bic = {"gene_ids":set(bic_genes),"n_genes":len(bic_genes),
               "sample_ids":set(samples),"n_samples":len(samples),
               "avgSNR":abs(avgSNR),"direction":direction}
        return bic
    else:
        return {"avgSNR":False}
    
def genesets2biclusters(modules, genes, exprs_np, exprs_data, ints2g_names, ints2s_names,
                        min_SNR = 0,min_n_samples=10, min_n_genes=2,
                        verbose = True):
    # Identify optimal sample set for each module: split samples into two sets in a subspace of each module
    # Filter out bad biclusters with too few genes or samples, or with low SNR
    t0 = time()
    filtered_bics = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0
    
    for mid in range(0,len(modules)):
        gene_ids = modules[mid]
        if len(gene_ids) >= min_n_genes: # makes sense to take 2+
            try:
                gene_ids = [genes[x] for x in gene_ids]
            except:
                print(mid,gene_ids)
                gene_ids = [genes[x] for x in gene_ids]
            bic = identify_opt_sample_set(gene_ids, exprs_np, exprs_data,
                                          min_n_samples=min_n_samples)
            avgSNR = bic["avgSNR"]
            if avgSNR == False:  # exclude biclusters with too few samples
                wrong_sample_number+=1
            elif avgSNR < min_SNR: # exclude biclusters with low avg. SNR 
                low_SNR += 1
            else:
                bic["id"] = i
                filtered_bics[i] = bic
                i+=1
      
    if verbose:
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %(len(modules),round(time()-t0,2)))
        print("Passed biclusters (>=%s genes, > %s SNR): %s"%(min_n_genes,min_SNR,i-1), file = sys.stdout)
        print("\tModules with not enough or too many samples:",wrong_sample_number, file = sys.stdout)      
        print("\tModules not passed avg. |SNR| threshold:", low_SNR, file = sys.stdout)
        print()
        
    # rename bicluster genes and samples 
    
    for i in filtered_bics.keys():
        bic = filtered_bics[i]
        bic["genes"] = set([ints2g_names[x] for x in bic["gene_ids"]])
        bic["samples"] = set([ints2s_names[x] for x in bic["sample_ids"]])
        if len(bic["genes"])>=min_n_genes and bic["avgSNR"]>min_SNR and verbose:
            print("\t".join(map(str,[str(bic["n_genes"])+"x"+str(bic["n_samples"]),
                                     round(bic["avgSNR"],3)," ".join(sorted(bic["genes"]))])),file = sys.stdout)
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