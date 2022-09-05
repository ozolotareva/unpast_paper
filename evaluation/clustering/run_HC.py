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
results_dir = args[2]
# results_dir = '/Users/fernando/Documents/Research/DESMOND2/methods/results/results_HC/'
# input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/simulated/A/example.tsv'
# input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/simulated/A/example_MEs.tsv'
df = pd.read_csv(input_file, sep='\t', index_col=0).T

for method in ['single', 'complete']:
    mergings = linkage(df, method=method, metric='euclidean')
    # dendrogram(mergings)
    # plt.show()
    for k in tqdm(range(1, 21)):
        cluster_labels = cut_tree(mergings, n_clusters=k).reshape(-1, )
        #print(cluster_labels)
        result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': cluster_labels}), input_df=df)
        result_k.to_csv(results_dir + input_file.split('/')[-1].replace('.tsv', f'_method_{method}_k_{k}.tsv'), sep='\t')