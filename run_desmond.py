#!/usr/bin/env python3

import argparse
import sys, os
import time, datetime
import random
import argparse
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from method2 import GM_binarization
from method2 import set_initial_conditions, sampling, get_consensus_modules
from method2 import genesets2biclusters,  write_bic_table


def run_DESMOND(
        exprs_file,
        out_dir = "./",
        basename = False,

        min_n_samples = -1,
        min_SNR = 1.5,
        
        alpha = 1.0,
        beta_K = 1.0,
            
        max_n_steps = 200,
        n_steps_averaged = 20,
        n_steps_for_convergence = 5,
        n_points_fit = 10,
    
        plot_all = True,
        verbose = True,
        seed = 0,
    ):
    random.seed(seed)
    np.random.seed(seed)

    biclusters = {} # UP and DOWN
    start_time = time.time()

    if basename:
        basename = basename
    else: 
        [date_h,mins] = str(datetime.datetime.today()).split(":")[:2]
        [date, hs] = date_h.split()
        basename = "results_"+hs+":"+mins+"_"+date 

    suffix  = ".alpha="+str(alpha)+",beta_K="+str(beta_K)+",minSNR="+str(min_SNR)
    if verbose:
        print("Will save output files to: %s*."% (out_dir + basename + suffix), file = sys.stdout)

    # read inputs
    exprs = pd.read_csv(exprs_file, sep="\t", index_col=0)

    # define minimal number of patients in a module
    if min_n_samples == -1:
        min_n_samples = int(max(10,0.05*exprs.shape[1])) # set to max(10, 5% of the cohort) 
    if verbose:
        print("Mininal number of samples in a module:",min_n_samples ,file=sys.stdout)
    ####1. Binarization ####
    sele_genes = ["SIDT1","CLIC6","BMPR1B", "FOXA1","GATA3","ESR1","ERBB2","GRB7"]
    binarized_expressions = GM_binarization(exprs,min_SNR,min_n_samples,verbose = verbose,
                                 plot=plot_all, plot_SNR_thr= 2.5,show_fits=sele_genes)    

    ### 2. Edge clustering ###
    
    # simplifying probability calculations
    N = exprs.shape[1]
    max_log_float = np.log(np.finfo(np.float64).max)
    n_exp_orders = 7 # ~1000 times 
    p0 = N*np.log(0.5)+np.log(beta_K)
    match_score = np.log((alpha*0.5+1)/(alpha))
    mismatch_score = np.log((alpha*0.5+0)/alpha)
    bK_1 = math.log(1+beta_K)
    
    i = 0
    for direction in ["UP","DOWN"]:
        t1 = time.time()
        exprs_bin = binarized_expressions[direction]
        genes = exprs_bin.columns.values
        # setting initial model state
        print("Searching for %s-regulated biclusters ..."%direction)
        moduleSizes, gene2Samples, nOnesPerSampleInModules, gene2Module, moduleOneFreqs, LP  = set_initial_conditions(exprs_bin, alpha,beta_K,verbose = verbose)
        K = len(moduleSizes)
        N = gene2Samples.shape[1]
        print("\t\tLP matrix memory usage: {:.2f}M".format(LP.nbytes/(1024*1024)),file = sys.stdout)

        # sampling
        t0 = time.time()
        gene2Module_history,n_final_steps,n_skipping_genes,P_diffs = sampling(LP,gene2Module, gene2Samples, nOnesPerSampleInModules,
                                                                              moduleSizes,moduleOneFreqs, p0, match_score,mismatch_score,
                                                                              bK_1, alpha, beta_K, max_n_steps=max_n_steps, tol = 0.1,
                                                                              n_steps_averaged = n_steps_averaged, n_points_fit = n_points_fit, 
                                                                              n_steps_for_convergence = n_steps_for_convergence, verbose=True)

        print("time:\tSampling (%s steps) fininshed in %s s." %(len(gene2Module_history),round(time.time()-t0,2)), file = sys.stdout)
        if plot_all:
            from method2 import plot_convergence
            plot_outfile = out_dir + basename +suffix+",ns_max=" + str(max_n_steps)+ ",ns_avg=" + str(n_steps_averaged) + ",ns_c="+str(n_steps_for_convergence) + ".convergence.svg"
            plot_convergence(n_skipping_genes, P_diffs,len(gene2Module_history)-n_final_steps,
                             n_steps_averaged, outfile=plot_outfile)

        # take the last (n_points_fit+n_steps_for_convergence) steps modules:
        # and get consensus edge-to-module membership
        consensus, nOnesPerSampleInModules, moduleSizes, moduleOneFreqs = get_consensus_modules(gene2Module_history[-n_final_steps:], 
                                                                                                LP, gene2Samples, gene2Module,
                                                                                                nOnesPerSampleInModules,moduleSizes, 
                                                                                                moduleOneFreqs, p0, 
                                                                                                match_score,mismatch_score,
                                                                                                bK_1,alpha,beta_K,N,K)

        print("\tEmpty modules:", len([x for x in moduleSizes if x == 0]),
              "\n\tNon-empty modules:",len([x for x in moduleSizes if x > 0]),file = sys.stdout)
        #### 3. Define biclusters and merge modules  ####
        exprs_np = exprs.loc[genes,:]
        ints2g_names = exprs_np.index.values
        ints2s_names = exprs_np.columns.values
        exprs_np = exprs_np.values
        exprs_sums = exprs_np.sum(axis=1)
        exprs_sq_sums = np.square(exprs_np).sum(axis=1)
        N = exprs.shape[1]
        exprs_data = N, exprs_sums, exprs_sq_sums
        # Identify optimal patient sets for each module: split patients into two sets in a subspace of each module
        # Filter out bad biclusters with too few genes or samples, or with low SNR
        filtered_bics = genesets2biclusters(exprs_np, exprs_data,moduleSizes,consensus,
                                min_SNR = 0.5,direction=direction,min_n_samples=min_n_samples,
                                verbose = verbose)
        # print info on found biclusters
        if verbose:
            for bic in filtered_bics:
                bic["id"] = i
                i+=1
                bic["genes"] = sorted([ints2g_names[x] for x in bic["genes"]])
                bic["samples"] = sorted([ints2s_names[x] for x in bic["samples"]])
                if len(bic["genes"])>2:
                    print("\t".join(map(str,[str(bic["n_genes"])+"x"+str(bic["n_samples"]),
                                             round(bic["avgSNR"],3)," ".join(bic["genes"])])),file = sys.stdout)
                        
        # save results 
        result_file_name = out_dir+basename+suffix+",direction="+direction
        write_bic_table(filtered_bics,result_file_name+".biclusters.tsv")
        if verbose:
            print("Runtime for %s-regulated biclusters: %s"%(direction, round(time.time()-t1,2)),file = sys.stdout)
        biclusters[direction] = filtered_bics
    
    biclusters = biclusters["UP"] + biclusters["DOWN"]
    result_file_name = out_dir+basename+suffix
    write_bic_table(biclusters,result_file_name+".biclusters.tsv")
    if verbose:
        print("Total runtime:",round(time.time()-start_time,2),file = sys.stdout)
        
    return biclusters

#exprs_file = "../../Expression/Harmonized_final/TCGA_micro_nolog2.z.13K_common.tsv"
#biclusters_TCGA_micro = run_DESMOND(exprs_file,basename = "TCGA_micro",min_SNR=2.0)

#exprs_file = "../../Expression/Harmonized_final/TCGA_micro_nolog2.z.13K_common.tsv"
#biclusters_TCGA_micro = run_DESMOND(exprs_file,basename = "TCGA_micro",min_SNR=2.0)


def parse_args():
    parser = argparse.ArgumentParser("DESMOND2 - is a tool for bicluster search")
    parser.add_argument('exprs_file', help="file with an input matrix")
    
    parser.add_argument('--out_dir', default="./")
    parser.add_argument('--basename', action='store_true')
    parser.add_argument('--min_n_samples', default=-1, type=int)
    parser.add_argument('--min_SNR', default=1.5, type=float)
    parser.add_argument('--alpha', default=1.0, type=float)
    parser.add_argument('--beta_K', default=1.0, type=float)
    parser.add_argument('--max_n_steps', default=200, type=int)
    parser.add_argument('--n_steps_averaged', default=20, type=int)
    parser.add_argument('--n_steps_for_convergence', default=5, type=int)
    parser.add_argument('--n_points_fit', default=10, type=int)
    parser.add_argument('--no_plot', dest='plot_all', action='store_false')
    parser.add_argument('--verbose', action='store_false')

    parser.add_argument('--seed', default=0, type=int, help="random seed")
    return parser.parse_args()
    

if __name__ == "__main__":
    args = parse_args()
    run_DESMOND(**vars(args))
