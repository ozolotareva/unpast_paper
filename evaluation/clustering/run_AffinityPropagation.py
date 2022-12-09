import sys, os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.cluster import AffinityPropagation
# import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
import numpy as np
import warnings
import os

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

seed_dict = dict()
seed_dict[1] = 57451
seed_dict[2] = 48699
seed_dict[3] = 22057
seed_dict[4] = 59467
seed_dict[5] = 43106

seed = seed_dict[int(result_file.split('_run')[-1].split('.')[0])]
np.random.seed(seed)
'''
input_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=5,m=4,std=1,overlap=no.exprs_z.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/HC/clustering/A.n_genes=5,m=4,std=1,overlap=no_run1.tsv'
'''

warnings.filterwarnings("ignore")
df = pd.read_csv(input_file, sep='\t', index_col=0).T

dampings = np.linspace(0.5, 1, 11)[:-1]
# max_iters = range(10, 1000, 100)
affinities = ['euclidean', 'precomputed']
# convergence_iters = range(1, 30, 5)
distance_metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean',
                    'hamming', 'jaccard', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
                    'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']

for mydamping in dampings:
    for myaffinity in affinities:
        for distance in distance_metrics:

            result_filename = result_file.replace('.tsv', f'_damping_{mydamping}_affinity_{myaffinity}_distance_{distance}_seed_{seed}.tsv')
            # print(result_filename)

            try:
                test = pairwise_distances(df, metric=distance)
                af = AffinityPropagation(damping=mydamping, preference=None, affinity='precomputed', random_state=seed)
                # af = AffinityPropagation(affinity='euclidean')
                clustering = af.fit(test)
                result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': clustering.labels_}), input_df=df)
                result_k.to_csv(result_filename, sep='\t')
            except:
                open(result_filename, 'a').close()
