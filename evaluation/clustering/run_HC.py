import sys, os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from scipy.cluster.hierarchy import linkage
# from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import cut_tree
# import matplotlib.pyplot as plt
from tqdm import tqdm
os.environ["OMP_NUM_THREADS"] = "1"
#  python3 run_kmeans.py ../datasets/DESMOND2_data_simulated/simulated/A/example.tsv  results


def reformat_cluster_results(clusters, input_df):
    res_reformat = []
    for labelno in range(max(clusters.label) + 1):
        # print(labelno)
        samples = set(clusters[clusters.label == labelno]['sample'].astype(str))
        n_samples = len(samples)
        frac = n_samples / input_df.shape[0]
        res_reformat.append([samples, n_samples, frac])
    res_reformat = pd.DataFrame(res_reformat, columns=['samples', 'n_samples', 'frac_samples'])
    return res_reformat


args = sys.argv
input_file = args[1]
result_file = args[2]
'''
input_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/A.n_genes=5,m=4,std=1,overlap=no.exprs_z.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/HC/clusters/A.n_genes=5,m=4,std=1,overlap=no_run1.tsv'
'''
#


df = pd.read_csv(input_file, sep='\t', index_col=0).T
distance_metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
# distance_metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulczynski1', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
distance_methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']

for method in distance_methods:
    for distance_metric in distance_metrics:
        print(method, distance_metric)
        if method in ['centroid', 'median', 'ward']:
            distance_metric = 'euclidean'
        mergings = linkage(df, method=method, metric=distance_metric)
        # dendrogram(mergings)
        # plt.show()
        for k in tqdm(range(1, 21)):
            cluster_labels = cut_tree(mergings, n_clusters=k).reshape(-1, )
            # print(cluster_labels)
            result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': cluster_labels}), input_df=df)
            result_k.to_csv(result_file.replace('.tsv', f'_method_{method}_distance_{distance_metric}_k_{k}.tsv'), sep='\t')