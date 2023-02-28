#!/usr/bin/env python3
import argparse
import random
import numpy as np

def run(exprs_file, basename='', out_dir="./",
                save=True, load = False,
                ceiling = 3,
                bin_method = "kmeans", clust_method = "WGCNA", 
                min_n_samples = 5, 
                show_fits = [],
                pval = 0.01, # binarization p-value
                modularity=1/3, similarity_cutoffs = -1, # for Louvain
                ds = 0, dch = 0.995, # for WGCNA
                alpha=1,beta_K = 1, max_n_steps= 100, n_steps_for_convergence = 5, # for DESMOND
                cluster_binary=False, 
                merge = 1,
                seed = -1,
                verbose = True, plot_all = False):
    
    import sys
    from time import time
    import pandas as pd
    
    start_time = time()
    
    # make sure that out_dir has '/' suffix
    if out_dir[-1] != '/':
        out_dir += '/'
    
    if seed == -1:
        seed = random.randint(0,1000000)
        print("seed=",seed,file = sys.stdout)
        
    if not basename: 
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
    exprs = validate_input_matrix(exprs,verbose = verbose)
    
    #set extreme z-scores to -x and x, e.g. -3,3
    if ceiling:    
        exprs[exprs>ceiling] = ceiling
        exprs[exprs<-ceiling] = -ceiling
        exprs.fillna(-ceiling, inplace = True)
        if verbose:
            print("Standardized expressions will be limited to [-%s,%s]:"%(ceiling,ceiling),file=sys.stdout)
    
    if verbose:
        print("Mininal number of samples in a bicluster:",min_n_samples ,file=sys.stdout)
    if min_n_samples < 5:
        print("min_n_samples is recommended to be >= 5", file= sys.stderr)

    ######### binarization #########
    from utils.method import binarize
    
    binarized_expressions, stats, null_distribution  = binarize(out_dir+basename, exprs=exprs,
                                 method=bin_method, save = save, load=load,
                                 min_n_samples = min_n_samples,pval=pval,
                                 plot_all = plot_all,show_fits = show_fits,
                                 verbose= verbose,seed=seed,prob_cutoff=0.5)
    
    ######### gene clustering #########
    
    if verbose:
        print("Clustering features ...\n",file=sys.stdout)
    feature_clusters, not_clustered, used_similarity_cutoffs = [], [], []
    if clust_method == "Louvain":
        from utils.method import run_Louvain, make_TOM
        #from utils.method import get_similarity_corr 
        from utils.method import get_similarity_jaccard
        
        for d in ["DOWN","UP"]:
            df = binarized_expressions.loc[:,stats["direction"]==d]
            if df.shape[0]>1:
                
                similarity = get_similarity_jaccard(df,verbose = verbose)
                #similarity = get_similarity_corr(df,verbose = verbose)
                
                if similarity_cutoffs  == -1: # guess from the data
                    similarity_cutoffs = np.arange(0.3,0.9,0.01)
                # if similarity cuttofs is a single value turns it to a list
                try: 
                    similarity_cutoffs = [elem for elem in similarity_cutoffs]
                except:
                    similarity_cutoffs = [similarity_cutoffs]
                
                # if modularity m is defined, choses a similarity cutoff corresponding to this modularity
                # and rund Louvain clustering
                modules, single_features, similarity_cutoff = run_Louvain(similarity,
                                                                          similarity_cutoffs = similarity_cutoffs,
                                                                          m = modularity, 
                                                                          verbose = verbose)
                used_similarity_cutoffs.append(similarity_cutoff)
                feature_clusters+= modules
                not_clustered+= single_features
            elif df.shape[0]==1:
                not_clustered+= list(df.index.values)
                used_similarity_cutoffs.append(None)
        used_similarity_cutoffs = ",".join(map(str,used_similarity_cutoffs))
        
    elif clust_method == "WGCNA":
        from utils.method import run_WGCNA
        # create unique suffix  for tmp files
        from datetime import datetime
        now = datetime.now()
        suffix = ".tmp_" + now.strftime("%y.%m.%d_%H:%M:%S")
        for d in ["DOWN","UP"]:
            fname = out_dir+basename+ "."+bin_method+".pval="+str(pval)+".seed="+str(seed)+"."+d+suffix+".tsv"
            df = binarized_expressions.loc[:,stats["direction"]==d]
            if df.shape[0]>1:
                modules, single_features = run_WGCNA(df,fname,deepSplit=ds,detectCutHeight=dch, verbose = verbose)  
                feature_clusters+= modules
                not_clustered+= single_features

    elif clust_method == "DESMOND":
        from utils.pgm import run_sampling

        # convergence
        n_steps_averaged = 10
        n_points_fit=10

        clustering_results ={}
 
        exprs_bin = binarized_expressions
        genes = exprs_bin.columns.values
        feature_clusters, not_clustered  = run_sampling(exprs_bin,alpha=alpha,beta_K=beta_K,f=0.5,
                    max_n_steps=max_n_steps, n_steps_averaged = n_steps_averaged,
                    n_points_fit = n_points_fit, tol = 0.1,
                    n_steps_for_convergence = n_steps_for_convergence,
                    verbose =verbose,plot_all=plot_all)
    
    else:
        print("'clust_method' must be 'WGCNA' or 'Louvain', or 'DESMOND'.",file=sys.stderr)

    ######### making biclusters #########
    if len(feature_clusters)==0:
        print("No biclusters found",file = sys.stderr)
        return pd.DataFrame()
    
    from utils.method import make_biclusters
    biclusters = make_biclusters(feature_clusters,binarized_expressions,exprs,null_distribution,
                             method = bin_method, merge = merge,
                             min_n_samples=min_n_samples, min_n_genes=2,
                             seed = seed,cluster_binary=False,verbose = verbose)

    
    from utils.method import write_bic_table
    suffix  = ".seed="+str(seed)+".bin="+bin_method+",pval="+str(pval)+",clust="+clust_method
    if clust_method == "WGCNA":
        suffix2 = ",ds="+str(ds)+",dch="+str(dch)
        modularity, similarity_cutoff = None, None
    elif clust_method == "Louvain":
        suffix2 = ",m="+str(round(modularity,2))
        ds, dhs = None, None
    write_bic_table(biclusters, out_dir+basename+suffix+suffix2+".biclusters.tsv",to_str=True,
                    add_metadata=True, seed = seed, min_n_samples = min_n_samples, pval = pval,
                    bin_method = bin_method, clust_method = clust_method, 
                    alpha=alpha, beta_K = beta_K, similarity_cutoff = used_similarity_cutoffs,
                    m=modularity, ds = ds, dch = dch)

    if verbose:
        print(out_dir+basename+suffix+suffix2+".biclusters.tsv", file = sys.stdout)
        print("Total runtime: {:.2f} s".format(time()-start_time ), file = sys.stdout)
    
    return biclusters    

def parse_args():
    parser = argparse.ArgumentParser("UnPaSt identifies differentially expressed biclusters in gene expression data.")
    seed = random.randint(0,1000000)
    parser.add_argument('--seed',metavar=seed, default=seed, type=int, help="random seed")
    parser.add_argument('--exprs', metavar="exprs.z.tsv", required=True, 
                        help=".tsv file with standardized gene expressions. The first column and row must contain unique gene and sample ids, respectively.")
    parser.add_argument('--out_dir', metavar="./", default="./", help  = 'output folder')
    parser.add_argument('--basename', metavar="biclusters.tsv", default = False, type=str, help  = 'output files prefix. If not specified, will be set to "results_"yy.mm.dd_HH:MM:SS""')
    parser.add_argument('--ceiling', default=3, metavar="3",  type=float, required=False, 
                        help="Absolute threshold for z-scores. For example, when set to 3, z-scores greater than 3 are set to 3 and z-scores less than -3 are set to -3. No ceiling if set to 0.")
    parser.add_argument('-s','--min_n_samples', metavar=5, default=5, type=int, help  = 'minimal number of samples in a bicluster.')
    parser.add_argument('-b','--binarization', metavar="kmeans", default="kmeans", type=str,
                        choices=["kmeans","ward",'GMM', 'Jenks'], help='binarization method')
    parser.add_argument('-p','--pval', metavar=0.01, default=0.01, type=float, help  = 'binarization p-value')
    parser.add_argument('-c','--clustering', metavar="WGCNA", default="WGCNA", type=str,
                        choices=['Louvain', 'WGCNA'], help='feature clustering method')
    # Louvain parameters
    parser.add_argument('-m','--modularity', default=1/3, metavar="1/3", type=float, help='Modularity corresponding to a cutoff for similarity matrix (Louvain clustering)')
    parser.add_argument('-r','--similarity_cutoffs', default=-1, metavar="-1", type=float, help='A cutoff or a list of cuttofs for similarity matrix (Louvain clustering). If set to -1, will be chosen authomatically from [1/5,4/5] using elbow method')
    # WGCNA parameters 
    parser.add_argument('--ds', default=0, metavar="0", type=int,choices=[0,1,2,3,4], help='deepSplit parameter, see WGCNA documentation')
    parser.add_argument('--dch', default=0.995, metavar="0.995", type=float, help='dynamicTreeCut parameter, see WGCNA documentation')
    parser.add_argument('--merge', default=1, metavar="1", type=float,help = "Whether to merge biclustres similar in samples with Jaccard index not less then the specified.")
    parser.add_argument('--load_binary', action='store_true', help = "loads binarized features from <basename>.<bin_method>.seed="+str(seed)+".binarized.tsv, statistics from *.binarization_stats.tsv and the background SNR distribution from <basename>.<bin_method>.n=<n_permutations>.seed="+str(seed)+".background.tsv")
    parser.add_argument('--save_binary', action='store_true', help = "saves binarized features to a file named as <basename>.<bin_method>.seed="+str(seed)+".binarized.tsv. If WGCNA is clustering method, binarized expressions are always saved. Also, files *.binarization_stats.tsv and *.background.tsv with binarization statistincs and background SNR distributions respectively will be created")
    parser.add_argument('--verbose', action='store_true')
    #parser.add_argument('--plot', action='store_true', help = "show plots")
    
    
    return parser.parse_args()
    

if __name__ == "__main__":
    args = parse_args()
        
    biclusters = run(args.exprs, args.basename, out_dir=args.out_dir,  
                save = args.save_binary, load = args.load_binary,
                ceiling = args.ceiling,
                bin_method = args.binarization, 
                clust_method = args.clustering,
                pval = args.pval,
                min_n_samples = args.min_n_samples, 
                show_fits = [],
                modularity = args.modularity, similarity_cutoffs = args.similarity_cutoffs, # for Louvain
                ds = args.ds, dch = args.dch, # for WGCNA
                cluster_binary = False, 
                merge = args.merge,
                seed = args.seed,
                #plot_all = args.plot,
                verbose = args.verbose)
