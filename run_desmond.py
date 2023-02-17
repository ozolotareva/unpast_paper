from run_unpast import run


def run_DESMOND(*args, **kwargs):
    """ Wrapper of 'run_unpast' to support old function name """
    return run(*args, **kwargs)

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