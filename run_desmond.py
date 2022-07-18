#!/usr/bin/env python3

import argparse
import random

def run_DESMOND(exprs_file, basename, out_dir="",  
                binarized_data = None, save=True, load = False,
                bin_method = "GMM", clust_method = "Louvain", cluster_binary=False, 
                min_n_samples = -1, show_fits = [],
                pval = 0.001,
                r = 1/3,
                alpha=1,beta_K = 1, 
                max_n_steps= 100, n_steps_for_convergence = 5,
                seed = -1,
                verbose = True, plot_all = False):
    
    import sys
    from time import time
    import pandas as pd
    
    start_time = time()
    
    if seed == -1:
        seed = random.randint(0,1000000)
        print("seed=",seed,file = sys.stdout)
        
    if basename:
        basename = basename
    else: 
        from datetime import datetime
        now = datetime.now()
        basename = "results_" + now.strftime("%y.%m.%d_%H:%M:%S")
        print("set output basename to", basename, file = sys.stdout)
        
    # read inputs
    exprs = pd.read_csv(exprs_file, sep="\t",index_col=0)
    if verbose:
        print("Read input from:",exprs_file ,file=sys.stdout)
        print("\t{} features x {} samples".format(exprs.shape[0],exprs.shape[1]) ,file=sys.stdout)

    #check if expressions are standardized (mean=0, std =1)
    from utils.method import validate_input_matrix
    exprs = validate_input_matrix(exprs)

    # define minimal number of samples
    if min_n_samples == -1:
        min_n_samples = round(min(0.5*exprs.shape[1],max(10,0.01*exprs.shape[1])))
    if verbose:
        print("Mininal number of samples in a bicluster:",min_n_samples ,file=sys.stdout)
    if min_n_samples < 10:
        print("min_n_samples is recommended to be >= 10", file= sys.stderr)

    ######### binarization #########
    from utils.method import binarize
    binarized_expressions,stats,empirical_snr = binarize(out_dir+"/"+basename, exprs=exprs,
                                 method=bin_method, save = save, load=load,
                                 min_n_samples = min_n_samples,pval=pval,
                                 plot_all = plot_all,show_fits = show_fits,
                                 verbose= verbose,seed=seed)

    ######### gene clustering #########
    if clust_method == "Louvain":
        from utils.method import run_Louvain, get_similarity_corr

        clustering_results = {}
        for d in ["UP","DOWN"]:
            similarity = get_similarity_corr(binarized_expressions[d],r = r)
            clustering_results[d] = run_Louvain(similarity, verbose = True)

    elif clust_method == "WGCNA":
        from utils.method import run_WGCNA

        clustering_results = {}
        for d in ["UP","DOWN"]:
            fname = out_dir+"/"+basename+ "."+bin_method+".binarized_"+d +".tsv"
            clustering_results[d] = run_WGCNA(fname, verbose = verbose)

    elif clust_method == "DESMOND":
        from utils.pgm import run_sampling

        # convergence
        n_steps_averaged = 10
        n_points_fit=10

        clustering_results ={}
        for d in ["UP","DOWN"]:
            exprs_bin = binarized_expressions[d]
            genes = exprs_bin.columns.values
            clustering_results[d] = run_sampling(exprs_bin,alpha=alpha,beta_K=beta_K,f=0.51,
                        max_n_steps=max_n_steps, n_steps_averaged = n_steps_averaged,
                        n_points_fit = n_points_fit, tol = 0.1,
                        n_steps_for_convergence = n_steps_for_convergence,
                        verbose =verbose,plot_all=plot_all)
    else:
        print("'clust_method' must be 'WGCNA' or 'DESMOND'.",file=sys.stderr)

    ######### making biclusters #########
    from utils.method import make_biclusters
    biclusters = make_biclusters(clustering_results,binarized_expressions,exprs,
                                min_n_samples=min_n_samples, min_n_genes=2,
                                seed = seed,cluster_binary=cluster_binary)

    from utils.method import write_bic_table
    suffix  = ".bin="+bin_method+",clust="+clust_method
    write_bic_table(biclusters, out_dir+basename+suffix+".biclusters.tsv",to_str=True,
                    add_metadata=True, seed = seed, min_n_samples = min_n_samples, pval = pval,
                    bin_method = bin_method, clust_method = clust_method, 
                    alpha=alpha, beta_K = beta_K,r = r)

    print(out_dir+basename+suffix+".biclusters.tsv")
    
    print("Total runtime: {:.2f} s".format(time()-start_time ),file = sys.stdout)
    
    return biclusters    

def parse_args():
    parser = argparse.ArgumentParser("DESMOND2 identifies differentially expressed biclusters in gene expression data.")
    parser.add_argument('--exprs', metavar="exprs.z.tsv", required=True, 
                        help=".tsv file with standardized gene expressions. The first column and row must contain unique gene and sample ids respectively.")
    parser.add_argument('--out_dir', metavar=".", default=".", help  = 'output folder')
    parser.add_argument('--basename', metavar="biclusters.tsv", default = False, type=str, help  = 'output files basename. If not specified, will be set to "results_"yy.mm.dd_HH:MM:SS""')
    parser.add_argument('-s','--min_n_samples', metavar=10, default=-1, type=int, help  = 'minimal number of samples in a bicluster.If not specified, will be automatically defined based on input sample size')
    parser.add_argument('-b','--binarization', metavar="GMM", default="GMM", type=str,
                        choices=['GMM', 'Jenks'], help='bianrization method')
    parser.add_argument('-p','--pval', metavar=0.001, default=0.001, type=float, help  = 'binarization p-value')
    parser.add_argument('-c','--clustering', metavar="Louvain", default="Louvain", type=str,
                        choices=['Louvain', 'WGCNA', 'DESMOND'], help='feature clustering method')
    parser.add_argument('-r', default=1/3, metavar="1/3", type=float, help='Pearsons correlation cutoff for similarity matrix (Louvain clustering)')
    parser.add_argument('--alpha', metavar=1.0, default=1.0, type=float, help = 'alpha parameter, positive float  (DESMOND clustering) ')
    parser.add_argument('--beta_K', metavar=1.0, default=1.0, type=float, help = 'beta/K parameter, positive float (DESMOND clustering)')
    parser.add_argument('--max_n_steps', metavar=200, default=200, type=int, help = 'maximal number of Gibbs sampling steps (DESMOND clustering)')
    parser.add_argument('--n_steps_for_convergence', metavar=5, default=5, type=int, help  = 'the number of Gibbs sampling steps when the convergence condition must hold (DESMOND clustering)')
    parser.add_argument('--plot', action='store_true', help = "show plots")
    parser.add_argument('--save_binary', action='store_true', help = "saves binarized expressions for up- and down-requlated genes to files named as <basename>.<binarization method>.binarized_[UP|DOWN].tsv. If WGCNA is clustering method, binarized expressions are always saved.")
    parser.add_argument('--verbose', action='store_true')
    seed = random.randint(0,1000000)
    parser.add_argument('--seed',metavar=seed, default=seed, type=int, help="random seed")
    
    return parser.parse_args()
    

if __name__ == "__main__":
    args = parse_args()
    save = args.save_binary
    if  args.clustering == "WGCNA":
        # WGCNA needs binarized expressions
        save = True
    run_DESMOND(args.exprs, args.basename, out_dir=args.out_dir,  
                binarized_data = None, save=save, load = False,
                bin_method = args.binarization, clust_method = args.clustering,
                pval = args.pval,
                cluster_binary=False, 
                min_n_samples = args.min_n_samples, show_fits = [],
                r = args.r, 
                alpha = args.alpha,beta_K = args.beta_K, # DESMOND
                max_n_steps = args.max_n_steps, 
                n_steps_for_convergence = args.n_steps_for_convergence,
                seed = args.seed,
                verbose = args.verbose, plot_all = args.plot)