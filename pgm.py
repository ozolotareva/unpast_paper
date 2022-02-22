import sys
import copy
import random
import pandas as pd
import numpy as np
from time import time
import math
import os

import matplotlib.pyplot as plt

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
        print("\t\t\tconvergence condition:",convergence,"; curve slope: {:.2f}".format(k), file=sys.stdout)
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
        gene_ndxs = np.arange(K)
        np.random.shuffle(gene_ndxs)
        for gene_ndx in gene_ndxs:
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
                print("\t\t\tfeatures oscilating in the last {} steps: {}".format(n_steps_averaged, n_skipping[-1]),file=sys.stdout)
            
        
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
            consensus.append(new_ndxs)
            
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



def run_sampling(exprs_bin,alpha=1.0,beta_K=1.0,f=0.33,
                 max_n_steps=50, n_steps_averaged = 10,
                 n_points_fit = 10, tol = 0.1, n_steps_for_convergence = 10,
                 verbose =True,plot_all=True):

    moduleSizes, gene2Samples, nOnesPerSampleInModules, gene2Module, moduleOneFreqs, LP  = set_initial_conditions(exprs_bin,alpha,beta_K,verbose = verbose)
    K = len(moduleSizes)
    N = gene2Samples.shape[1]
    if verbose:
        print("\t\tLP memory, {:.2f} M".format(LP.nbytes/(1024*1024)),file = sys.stdout)
    
    # simplifying probability calculations
    max_log_float = np.log(np.finfo(np.float64).max)
    n_exp_orders = 7 # ~1000 times 
    p0 = N*np.log(0.5)+np.log(beta_K)
    match_score = np.log((alpha*0.5+1)/(alpha))
    mismatch_score = np.log((alpha*0.5+0)/alpha)
    bK_1 = math.log(1+beta_K)
    genes = exprs_bin.columns.values
    
    t0 = time()
    gene2Module_history,sampling_steps,n_skipping_genes = sampling(LP,gene2Module, gene2Samples, nOnesPerSampleInModules,moduleSizes,moduleOneFreqs,
                                                                   p0, match_score,mismatch_score, bK_1, alpha, beta_K, 
                                                                   max_n_steps=max_n_steps, n_steps_averaged = n_steps_averaged, 
                                                                   n_points_fit = n_points_fit, tol = tol, 
                                                                   n_steps_for_convergence = n_steps_for_convergence, verbose=verbose)

    print("time:\tSampling ({} steps) fininshed in {:.2f} s.".format(len(gene2Module_history),time()-t0), file = sys.stdout)

    # plot the number of oscilating genes per step
    if plot_all:
        #suffix = ".alpha="+str(alpha)+",beta_K="+str(beta_K)+",pv="+str(snr_pval)
        #save_plot = out_dir + basename +suffix + ".convergence.svg"
        save_plot = False
        plot_convergence(n_skipping_genes[0:], len(gene2Module_history)-sampling_steps-0,
                                 alpha=alpha, beta_K=beta_K,
                                 n_steps_averaged=n_steps_averaged,
                                 n_points_fit = n_points_fit, 
                                 n_steps_for_convergence = n_steps_for_convergence,
                                 save_plot=save_plot)
    
    # take the last (n_points_fit+n_steps_for_convergence) steps modules:
    # and get consensus feature-to-module membership
    consensus = get_consensus_modules(gene2Module_history[-sampling_steps:],verbose= verbose,f = f)
    
    # 
    modules = []
    not_clustered = []
    
    ints2g_names = exprs_bin.columns.values
    for mid in range(0,len(consensus)):
        gene_ids = consensus[mid]
        if len(gene_ids) > 1:
            modules.append(ints2g_names[gene_ids].tolist())
        elif len(gene_ids) == 1:
            not_clustered.append(ints2g_names[gene_ids[0]])
            
    if verbose:
        print("\t{} modules and {} not clustered genes".format(len(modules),len(not_clustered)),file = sys.stdout)
        print("time\tmodules detected in in {:.2f} s.".format(time()-t0),file = sys.stdout)
    return (modules,not_clustered)
