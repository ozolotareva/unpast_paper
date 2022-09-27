import sys
import copy
import random
import pandas as pd
import numpy as np
from time import time
import math
import itertools
import warnings
import os

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from fisher import pvalue
from scipy.interpolate import interp1d
import jenkspy

import matplotlib.pyplot as plt
import statsmodels.api as sm

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


def calc_SNR(ar1, ar2):
    return (np.mean(ar1) - np.mean(ar2)) / (np.std(ar1) + np.std(ar2))

def rand_norm_splits(N, min_n_samples, snr_pval = 0.05,seed = 1):
    # empirical distribuition of SNR depending on the bicluster size 
    # generates normal samples of matching size and split them into bicluster and background
    min_n_perm = int(5*1.0/snr_pval)
    n_perm = max(min_n_perm,int(100000/N))
    print("total samples: %s, min_n_samples: %s, n_permutations: %s"%(N,min_n_samples,n_perm))
    sizes = np.arange(min_n_samples,int(N/2+2))#.shape
    snr_thresholds = np.zeros(sizes.shape)
    np.random.seed(seed=seed)
    for s in sizes:
        snrs = np.zeros(n_perm)
        for i in range(0,n_perm):
            x = np.random.normal(size = N)
            x.sort()
            snrs[i] = calc_SNR(x[s:], x[:s]) #(x[s:].mean()-x[:s].mean())/(x[s:].std()+x[:s].std())
        snr_thresholds[s-min_n_samples]=np.quantile(snrs,q=1-snr_pval)
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
        ### TBG
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
                
        if (gene in show_fits or SNR > plot_SNR_thr) and plot:
            print("Gene %s: SNR=%s, pos=%s, neg=%s"%(gene, round(SNR,2), len(up_group), len(down_group)))
            plt.hist(down_group, bins=n_bins, alpha=0.5, color=down_color,range=hist_range)
            plt.hist(up_group, bins=n_bins, alpha=0.5, color=up_color,range=hist_range)
            plt.show()
    
    binarized_expressions["UP"] = pd.DataFrame.from_dict(binarized_expressions["UP"])
    binarized_expressions["DOWN"] = pd.DataFrame.from_dict(binarized_expressions["DOWN"])
    
    # logging
    if verbose:
        print("Total runtime",round(time()-t0,2), "s for ", len(exprs),"genes")
        #print("Genes passed SNR threshold of %s:"%round(min_SNR,2))
        print("\tup-regulated genes:", binarized_expressions["UP"].shape[1])
        print("\tdown-regulated genes:", binarized_expressions["DOWN"].shape[1])
        print("\tambiguous genes:", n_ambiguous )
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
        print("Total runtime",round(time()-t0,2), "s for ", len(exprs),"genes")
        #print("Genes passed SNR threshold of %s:"%round(min_SNR,2))
        print("\tup-regulated genes:", df_p.shape[1])
        print("\tdown-regulated genes:", df_n.shape[1])
        print("\tinexplicit genes:", stat['n_inexplicit'])

    return {"UP":df_p, "DOWN":df_n}

################## 2. Probabiliatic clustering #############
@jit_if_available
def calc_lp(gene_ndx,module_ndx,gene2Samples,
            nOnesPerSampleInModules,moduleSizes,
            moduleOneFreqs,p0,match_score,mismatch_score,bK_1,N,
            alpha,beta_K):
    
    # 1. Prepare vals: (n_ones_per_pat, m_size, alpha, gene_vector, beta_K)
    #      and return answer in special cases
    
    m_size = moduleSizes[module_ndx]
    if m_size == 0:
        # gene is removed, the module is empty
        return p0
    
    gene_vector = gene2Samples[gene_ndx,] 
    n_ones_per_pat = nOnesPerSampleInModules[module_ndx,]
    if gene_ndx == module_ndx: # remove the gene from its own module
        if m_size == 1:
            return p0
        m_size -=1
        n_ones_per_pat = n_ones_per_pat - gene_vector

    # if a module is composed of a single gene
    if m_size == 1:
        # just count number of matches and mismatches and
        # n_matches =  np.inner(n_ones_per_pat,gene_vector)
        n_matches =  np.sum(n_ones_per_pat[gene_vector==1])

        return n_matches*match_score + (N-n_matches)*mismatch_score + bK_1
    
    # 2. usual case if a module contains more than one gene
    return calc_lp_formula(n_ones_per_pat, m_size, alpha, gene_vector, beta_K)

@jit_if_available
def calc_lp_formula(n_ones_per_pat, m_size, alpha, gene_vector, beta_K): 
    beta_term = math.log(m_size+beta_K)
    
    # alpha_term
    # ones-matching
    oneRatios = (n_ones_per_pat+alpha/2)/(m_size+alpha)
    ones_matching_term = np.sum(np.log(oneRatios)[gene_vector == 1])

    # zero-matching
    # zeroRatios = (m_size-n_ones_per_pat+alpha/2)/(m_size+alpha)
    zeroRatios = 1 - oneRatios
    zeros_matching_term = np.sum(np.log(zeroRatios)[gene_vector == 0])

    return ones_matching_term + zeros_matching_term + beta_term

def calc_lp_column(module_ndx,gene2Samples,
            nOnesPerSampleInModules,moduleSizes,
            moduleOneFreqs,p0,match_score,mismatch_score,bK_1,N,
            alpha,beta_K):
    """Same as calc_lp, but for all the genes simultaneously"""
    m_size = moduleSizes[module_ndx]
    if m_size == 0:
        # all genes are removed, the module is empty, alfl answers are p0
        return p0

    n_ones_per_pat = nOnesPerSampleInModules[module_ndx,]
    if m_size == 1:
        # if a module is composed of a single gene
        n_matches = np.dot(gene2Samples, n_ones_per_pat)
        vals = n_matches*match_score + (N-n_matches)*mismatch_score + bK_1

    else: 
        # ones-matching
        oneRatios = (n_ones_per_pat+alpha/2)/(m_size+alpha)
        # ones_matching_term = np.dot((gene2Samples == 1), np.log(oneRatios))
        ones_matching_term = np_faster_dot(gene2Samples == 1, np.log(oneRatios))

        # zero-matching
        zeroRatios = 1 - oneRatios  # = (m_size-n_ones_per_pat+alpha/2)/(m_size+alpha)
        # zeros_matching_term = np.dot((gene2Samples == 0), np.log(zeroRatios))
        zeros_matching_term = np_faster_dot(gene2Samples == 0, np.log(zeroRatios))

        beta_term = math.log(m_size+beta_K)
        vals = ones_matching_term + zeros_matching_term + beta_term

    # calc LP[m,m] with the less optimized func
    vals[module_ndx] = calc_lp(module_ndx,module_ndx,gene2Samples,nOnesPerSampleInModules,moduleSizes,
                            moduleOneFreqs,p0,match_score,mismatch_score,bK_1,N,alpha,beta_K)
    return vals

@jit_if_available
def np_faster_dot(np_a, np_b): 
    # translating to float64 for jit compilation
    # may be a bit faster in some cases
    return np.dot(
        np_a.astype(np.float64),
        np_b.astype(np.float64)
    )

def set_initial_conditions(df, alpha,beta_K,verbose = True):
    t_0 = time()
    N = df.shape[0] # number of samples
    K = df.shape[1] # initial number of modules
    p0 = N*np.log(0.5)+np.log(beta_K)
    match_score = np.log((alpha*0.5+1)/(alpha))
    mismatch_score = np.log((alpha*0.5+0)/alpha)
    bK_1 = math.log(1+beta_K)
    #print("\t\tKxN=%sx%s"%(K,N))
    #print("\t\tp0=",p0)
    
    # p0, match_score, mismatch_score, bK_1
    genes = df.columns.values
    # 1. the number of genes inside each component, initially 1 for each gene
    moduleSizes=np.ones(K,dtype=np.int)
    
    # 2. a binary (int) matrix of size KxN that indicates the samples genes
    gene2Samples = df.T.values
    
    # 3. a binary matrix of size K by m that stores the total number of ones per sample in each module,
    # initially equal to 'gene2Samples'
    nOnesPerSampleInModules = copy.copy(gene2Samples)

    #4. initial module id
    gene2Module = list(range(0,K))
    gene2Module = np.array(gene2Module)
    
    #5. moduleOneFreqs
    moduleOneFreqs = []
    for g in range(0,K):
        moduleOneFreqs.append(float(sum(gene2Samples[g,]))/N)
    moduleOneFreqs = np.array(moduleOneFreqs)
    
    #6. setting initial LPs: i = gene, j=module
    LP = np.zeros((K,K),dtype=np.float)
    for i in range(0,K):
        #if (i+1)% 1000==0:
        #    if verbose:
        #        print("\t",i+1,"genes processed in ",round(time()- t_0,1),"s")
        for j in range(i,K):
            LP[i,j] = calc_lp(i,j,gene2Samples,
            nOnesPerSampleInModules,moduleSizes,
            moduleOneFreqs,p0,match_score,mismatch_score,bK_1,N,alpha, beta_K)
            LP[j,i] = LP[i,j]
    print("time:\tInitial state created in",round(time()-t_0, 1) , "s.", file = sys.stdout)
    return moduleSizes, gene2Samples, nOnesPerSampleInModules, gene2Module, moduleOneFreqs, LP

def adjust_lp(log_probs,n_exp_orders=7):
    # adjusting the log values before normalization to avoid under-flow
    max_p = max(log_probs)
    probs = []
    for lp in log_probs:
        # shift all probs to set max_prob less than log(max_np.float)  
        adj_lp = lp - max_p
        # set to minimal values all which less then 'n_orders' lesser than p_max to zeroes
        if adj_lp >= - n_exp_orders:
            probs.append(np.exp(adj_lp))
        else:
            probs.append(0)
    probs = probs/sum(probs)
    return probs

### functions for checking of convergence conditions ###

def count_skipping_genes(labels):
    n_skipping_genes = 0
    for gene in range(0,labels.shape[1]):
        states = labels[:,gene]
        unique = np.unique(states)
        if len(unique)> 1:
            n_skipping_genes +=1
    return n_skipping_genes


def check_convergence_conditions(n_skipping,n_min,n_max,tol=0.05, verbose = True):
    n_points = len(n_skipping)
    n_skipping = np.array(n_skipping,dtype=float)
    
    if n_max - n_min >0:
        # scale
        n_skipping = (n_skipping-n_min)/(n_max - n_min)*n_points
        # fit line
        A = np.vstack([range(0,n_points), np.ones(n_points)]).T
        k,b = np.linalg.lstsq(A, n_skipping, rcond=None)[0] # y = kx + b
    else:
        k = 0 
    
    if abs(k)<tol:
        convergence = True
    else:
        convergence = False
    if verbose:
        print("\tConverged:",convergence,"#skipping edges slope:",round(k,5))#, "RMS(Pn-Pn+1) slope:",round(k2,5))
    return convergence  

# @jit_if_available
def apply_changes(gene_ndx,new_module,curr_module,LP,gene2Samples, gene2Module, nOnesPerPatientInModules,moduleSizes,
                moduleOneFreqs, p0,match_score,mismatch_score, bK_1,alpha,beta_K,N,K, calc_LPs=True):
    """Moves the gene from the current module to the new one
        and updates gene2Module, nOnesPerPatientInModules and moduleSizes respectively.
        K - quantity of genes
        M (=K) - quantitiy of modules
        S - quantity of samples (patients)
        
        see also set_initial_conditions for more info

    Args: 
        gene_ndx(int): 
        new_module(int): 
        curr_module(int):
        LP(np.array KxM): matrix of transmission probabilities
        gene2Samples(np.array KxS): matrix of samples co-regulated with genes
        gene2Module(list): list of i'th module
        nOnesPerPatientInModules(np.array MxS): precalculated sums for optimization
            nOnesPerPatientInModules[m] = sum(gene2Samples[g] for g in genes if gene2Module[g] == m)
        moduleSizes(np.array M) - genes in module
        
        moduleOneFreqs - is not used
        N - is not used

        p0(float) - precalculated default probability
        match_score(float) - precalculated const
        mismatch_score(float) - precalculated const
        bK_1(float) - precalculated const
        
        alpha(float) - alg. global parameter
        beta_K(float) - alg. global parameter

        K(int) - count of Genes

        calc_LPs(bool) - should LPs be recalculated? 

    """
    # update the gene module membership
    gene2Module[gene_ndx] = new_module
    # for this edge no probabilities change
    
    # reduce curr_module size and nOnesPerPatientInModules
    gene_vector = gene2Samples[gene_ndx,]
    nOnesPerPatientInModules[curr_module,] = nOnesPerPatientInModules[curr_module,] - gene_vector
    moduleSizes[curr_module,]-=1
    
    # increase new_module
    nOnesPerPatientInModules[new_module,] = nOnesPerPatientInModules[new_module,] + gene_vector
    moduleSizes[new_module,]+=1
    
    # update LPs for all genes contacting curr and new modules
    if calc_LPs:
        for module in curr_module, new_module:
            ndx_old_val = LP[gene_ndx, module] # for this gene no probabilities changed
            LP[:, module] = calc_lp_column(module,gene2Samples,nOnesPerPatientInModules,moduleSizes,
                                                            moduleOneFreqs,p0,match_score,mismatch_score,bK_1,N,alpha,beta_K)
            LP[gene_ndx, module] = ndx_old_val






def sampling(LP,gene2Module, gene2Samples,nOnesPerPatientInModules,moduleSizes,
             moduleOneFreqs, p0, match_score, mismatch_score, bK_1, alpha, beta_K,
             max_n_steps=100,n_steps_averaged = 5,n_points_fit = 10,tol = 0.05,
             n_steps_for_convergence = 10,verbose=True):
    
    K = len(gene2Module)
    N  = gene2Samples.shape[1]
    t_ =  time()
    gene2Module_history = [copy.copy(gene2Module)]
    sampling_steps = n_steps_for_convergence 
    is_converged = False
    
    for step in range(1, max_n_steps):
        #if verbose:
        #    print("step", step,file = sys.stdout)
        not_changed_genes = 0
        t_0 = time()
        t_1=t_0
        i = 1
        for gene_ndx in range(0, K):
            # adjust LogP and sample a new module
            P_adj = adjust_lp(LP[gene_ndx,:], n_exp_orders=7)
            curr_module = gene2Module[gene_ndx]
            new_module = np.random.choice(range(0,K), p=P_adj) 

            # update network and matrices if necessary
            if new_module != curr_module:
                apply_changes(gene_ndx,new_module,curr_module,LP,gene2Samples, gene2Module, nOnesPerPatientInModules,moduleSizes,
                moduleOneFreqs, p0,match_score,mismatch_score, bK_1,alpha,beta_K, N, K)
                
            else:
                not_changed_genes +=1#
            i+=1
            #if i%1000 == 0:
            #    if verbose:
            #        print(i,"\t\tgenes processed in",round(time()- t_1, 1) , "s runtime...",file=sys.stdout)
            #    t_1 = time()
        if verbose:
            print("\t\tstep ",step,round(time() - t_0, 1) , "s", file = sys.stdout)
        
        gene2Module_history.append(copy.copy(gene2Module))
        
        if step == n_steps_averaged:
            is_converged = False
            n_times_cc_fulfilled = 0
            labels = np.asarray(gene2Module_history[step-n_steps_averaged:step])
            n_skipping = [count_skipping_genes(labels)]
            
        if step > n_steps_averaged:
            labels = np.asarray(gene2Module_history[step-n_steps_averaged:step])
            n_skipping.append(count_skipping_genes(labels))
            if verbose:
                print("Genes skipping in last %s steps: %s"%(n_steps_averaged, n_skipping[-1]),file=sys.stdout)
            
        
        if  step >= n_steps_averaged + n_points_fit:
            last_steps = n_skipping[-n_points_fit:]
            if len(np.unique(last_steps))==1: # stop when n_skipped is the same in last n_points_fit steps 
                if verbose:
                    print("The model converged after", step+1,"steps.", file = sys.stdout)
                    print("Consensus of last",sampling_steps,"states will be taken")
                    print("Sampling runtime",round(time()- t_ ,1) , "s", file = sys.stdout)
                return gene2Module_history, sampling_steps,n_skipping
            
            n_min, n_max = min(n_skipping),max(n_skipping)
            is_converged = check_convergence_conditions(last_steps,n_min, n_max,tol=tol,verbose = verbose)
                
            if is_converged:
                n_times_cc_fulfilled +=1
            else:
                n_times_cc_fulfilled = 0

            if n_times_cc_fulfilled > n_steps_for_convergence: # stop if convergence is True for the last n steps
                if verbose:
                    print("The model converged after", step+1,"steps.", file = sys.stdout)
                    print("Consensus of last",sampling_steps,"states will be taken")
                    print("Sampling runtime",round(time()- t_ ,1) , "s", file = sys.stdout)
                return gene2Module_history, sampling_steps,n_skipping
    
    if verbose:
        print("The model did NOT converge after", step+1,"steps.", file = sys.stdout)
        print("Consensus of last",sampling_steps,"states will be taken")
        print("Sampling runtime",round(time()- t_ ,1) , "s", file = sys.stdout)
        
    return gene2Module_history,sampling_steps,n_skipping


def plot_convergence(n_skipping_edges,thr_step,
                     alpha=1.0, beta_K=1.0,
                     n_steps_averaged=10,n_points_fit = 20, 
                     n_steps_for_convergence = 10,
                     save_plot = False):
    
    # plots numnber of oscilating genes
    steps = range(n_steps_averaged,n_steps_averaged+len(n_skipping_edges))
    fig, ax = plt.subplots(1, 1,figsize=(15,5))
    ax.set_title("Model convergence")
    ax.plot(steps, n_skipping_edges,'b.-')
    ax.axvline(thr_step,color="red",linestyle='--') 
    ax.set_ylabel("#genes oscilating on the last "+str(int(n_steps_averaged))+" steps")
    ax.set_xlabel('step')
    text = 'alpha=%s; beta/K=%s;\n'%(alpha,beta_K)
    text += "ns_a=%s; n_points_fit=%s; ns_c=%s"%(n_steps_averaged,n_points_fit,n_steps_for_convergence)
    ax.text(0.5, 0.90, text,
    verticalalignment='top', horizontalalignment='center',
    transform=ax.transAxes,
    color='black', fontsize=15)
    
    if save_plot:
        plt.savefig(save_plot, transparent=True)  


def get_consensus_modules(gene2module_history, f = 0.25, verbose = False):
    
    # g->m states from sampling phase
    labels = np.asarray(gene2module_history)
    
    # identify modules which genes ocsilate
    K = labels.shape[1] # number of genes and modules
    n_steps = labels.shape[0] # number of steps in sampling phase
    
    # consensus module -> genes 
    consensus = [[] for i in range(0,K)]
    
    # genes in many modules
    n_mult_genes = 0
    for gene_ndx in range(0,K):
        unique, counts = np.unique(labels[:,gene_ndx], return_counts=True)
        
        if len(unique) == 1:
            consensus[unique[0]].append(gene_ndx) 
            #print("Gene %s -> [%s]" %(gene_ndx, ",".join(map(str,unique))),file= sys.stdout)
        else:
            freqs = np.array(counts)/n_steps
            new_ndxs  = [unique[i] for i in range(len(freqs)) if freqs[i]>=f]
            for i in new_ndxs:
                consensus[i].append(gene_ndx) 
            if len(new_ndxs)>1 and verbose:
                n_mult_genes +=1 
                #print("Gene %s -> [%s]" %(gene_ndx, ",".join(map(str,new_ndxs))),file= sys.stdout)
                #print(gene_ndx ,":",unique, freqs)
            consensus.append(new_ndxs)
    print("Genes assigned to more than 1 module:",n_mult_genes, file =sys.stderr)
    if verbose:
        print()
        print("size(genes)\tn_modules")
        module_sizes = {}
        for m in consensus:
            s = len(m)
            if s in module_sizes.keys():
                module_sizes[s] +=1
            else:
                module_sizes[s] =1
        for s in sorted(module_sizes.keys()):
            print(s, "\t\t", module_sizes[s])
            
    return consensus
        
################################## 3. Post-processing ####################

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
    
def genesets2biclusters(consensus, exprs_np, exprs_data, ints2g_names, ints2s_names,
                        min_SNR = 0,min_n_samples=10, min_n_genes=2,
                        verbose = True):
    # Identify optimal sample set for each module: split samples into two sets in a subspace of each module
    # Filter out bad biclusters with too few genes or samples, or with low SNR
    t0 = time()
    filtered_bics = {}
    wrong_sample_number = 0
    low_SNR = 0
    i = 0
    
    for mid in range(0,len(consensus)):
        gene_ids = consensus[mid]
        if len(gene_ids) >= min_n_genes: # makes sense to take 2+
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
        print("time:\tIdentified optimal sample sets for %s modules in %s s." %(len(consensus),round(time()-t0,2)))
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

###### save and read modules #####
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

def read_bic_table(results_file_name):
    if not os.path.exists(results_file_name):
        return pd.DataFrame()
    resulting_bics = pd.read_csv(results_file_name,sep = "\t",index_col=0)
    if len(resulting_bics) ==0:
        return pd.DataFrame()
    else:
        resulting_bics["genes"] = resulting_bics["genes"].apply(lambda x: set(x.split(" ")))
        resulting_bics["samples"] = resulting_bics["samples"].apply(lambda x: set(x.split(" ")))
    #resulting_bics.set_index("id",inplace=True)
    
    return resulting_bics


#### permutation tests ####

def flip_direction(x):
    if x=="UP":
        return "DOWN"
    if x=="DOWN":
        return "UP"
    
def generate_null_dist(exprs,g_size=3,n_permutations=100):
    """Generates N biclusters of given size in genes"""
    
    null_dist_SNR = np.zeros(n_permutations)
    for i in range(0,n_permutations):
        e = exprs.sample(n = g_size).values
        labels = KMeans(n_clusters=2, random_state=0,n_init=1,max_iter=100).fit(e.T).labels_
        ndx0 = np.where(labels == 0)[0]
        ndx1 = np.where(labels == 1)[0]
        e1 = e[:,ndx1]
        e0 = e[:,ndx0]
        SNR = np.mean(abs(e0.mean(axis=1)-e1.mean(axis=1))/(e0.std(axis=1)+e1.std(axis=1)))
        null_dist_SNR[i] = SNR
    return null_dist_SNR


def calc_e_pv(row,snr_dist):
    SNR=row["avgSNR"]
    e_pval = (1+snr_dist[snr_dist >= SNR].shape[0])/snr_dist.shape[0]
    row["e_pval"] = e_pval
    return row

def do_permutations(basename,exprs,in_dir = ".",out_dir=".",alphas=[1.0], 
                betas=[1.0], p_vals =[0.005], runs=[1],FDR = 0.01):
    all_bics = {}
    all_bic_counts = []
    for pv in p_vals:
        for r in runs:
            for a in alphas:
                for b in betas:
                    bic_fname = basename +"."+str(r)+".alpha="+str(a)+",beta_K="+str(b)+",pv="+str(pv)+".biclusters.tsv"
                    print("Read",in_dir+bic_fname)
                    bics = read_bic_table(in_dir+bic_fname)
                    bic_fname = ".".join(bic_fname.split(".")[:-1])
                    bics["genes"] = bics["genes"].apply(lambda x:" ".join(sorted(x)))
                    bics["samples"] = bics["samples"].apply(lambda x:" ".join(sorted(x)))
                    bics = bics.sort_values(by="avgSNR",ascending = False)

                    # of biclusters of similar genes, keep those with higher SNR 
                    print("\tAll biclsuters:",bics.shape[0])
                    bics = bics.drop_duplicates(subset = "genes") # duplicates appear when alpha is large
                    print("\t\tdeduplicated:",bics.shape[0])
                    N = exprs.shape[1]
                    
                    # flip too large biclusters with > N/2 samples
                    if bics.loc[bics["n_samples"]>N/2,:].shape[0]>0:
                        print("\t\tflip %s biclusters"%bics.loc[bics["n_samples"]>N/2,:].shape[0])
                        bics.loc[bics["n_samples"]>N/2,"direction"] = bics.loc[bics["n_samples"]>N/2,"direction"].apply(flip_direction)
                        bics.loc[bics["n_samples"]>N/2,"n_samples"] = bics.loc[bics["n_samples"]>N/2,"n_samples"]*-1+N
                    
                    all_bics[bic_fname] = bics
                    bic_counts = bics.loc[:,["n_genes","genes"]].groupby("n_genes").agg("count")

                    bic_counts = bics.loc[:,["n_genes","genes"]].groupby("n_genes").agg("count")
                    bic_counts.columns = [bic_fname]
                    all_bic_counts.append(bic_counts)

    # define maximal number of biclusters of given size per run
    all_bic_counts= pd.concat(all_bic_counts,axis=1)
    all_bic_counts= all_bic_counts.fillna(0)
    g_sizes = all_bic_counts.index.values
    # number of permutations
    n_permutations = all_bic_counts.max(axis=1).apply(lambda x: max(int(x/FDR),100)).to_dict()

    # generate null distributions
    t0 = time()
    SNR_dists = {}
    for g_size in g_sizes:
        n_perm = int(n_permutations[g_size])
        print("size: %s - %s permutations"%(g_size,n_perm))
        t1 = time()
        snr_dist = generate_null_dist(exprs,g_size=g_size,n_permutations=n_perm)
        print("\tSNR=",round(np.quantile(snr_dist,0.95),2),round(time()-t1,2),"s")
        SNR_dists[g_size] = snr_dist
    print("Time to generate null distributions:",round(time()-t0,2),file =sys.stdout)
    
    for bic_fname in all_bics.keys():
        bics = all_bics[bic_fname]
        fdr_bics = []
        for g_size in g_sizes:
            df = bics.loc[bics["n_genes"]==g_size,:]
            snr_dist = SNR_dists[g_size]
            # calculate empirical p-val
            df = df.apply(lambda row: calc_e_pv(row,snr_dist),axis=1)
            if df.shape[0]>0:
                # apply BH-correction
                rejected, q_val = fdrcorrection(df["e_pval"], method='indep', is_sorted=False)
                df["q_val"] = q_val
                n_passed = df.loc[df["q_val"]<=FDR,:].shape[0]
                print("\t%s biclusters with %s genes"%(n_passed,g_size),file=sys.stdout)
                fdr_bics.append(df)
        fdr_bics = pd.concat(fdr_bics,axis=0)

        print(bic_fname,":\tTotal biclusters passed:",fdr_bics.loc[fdr_bics["q_val"]<=FDR,:].shape[0],file=sys.stdout)
        write_bic_table(fdr_bics, out_dir + bic_fname+".qval.tsv",to_str=False)
        fdr_bics = fdr_bics.loc[fdr_bics["q_val"]<=FDR,:]
        write_bic_table(fdr_bics,out_dir + bic_fname+".FDR_"+str(FDR)+".tsv",to_str=False)
        all_bics[bic_fname] = fdr_bics
    return all_bics