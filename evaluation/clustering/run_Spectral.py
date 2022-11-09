import sys, os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.cluster import SpectralClustering
# import matplotlib.pyplot as plt
import numpy as np
import re

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
input_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/HC/clustering/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'
'''
#
df = pd.read_csv(input_file, sep='\t', index_col=0).T

n_clusters = range(1, 21)
eigen_solvers = ['arpack', 'lobpcg', 'amg']
affinities = ['nearest_neighbors', 'rbf']
N_neighbors = [2, 10, 20, 30, 40, 50]
label_assigners = ['kmeans', 'discretize', 'cluster_qr']

for n_cluster in n_clusters:
    for eigen_solver in eigen_solvers:
        for affinity in affinities:
            for label_assigner in label_assigners:
                if affinity == 'nearest_neighbors':
                    for nn in N_neighbors:
                        result_filename = result_file.replace('.tsv', f'_k_{n_cluster}_eignesolver_{eigen_solver}_affinity_{re.sub("_", "", affinity)}_assignlabels_{re.sub("_", "", label_assigner)}_nneighbors_{nn}_seed_{seed}.tsv')
                        if not os.path.exists(result_filename):
                            try:
                                sc = SpectralClustering(n_clusters=n_cluster, eigen_solver=eigen_solver, affinity=affinity, assign_labels=label_assigner, n_neighbors=nn).fit(df)
                                result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': sc.labels_}), input_df=df)
                                result_k.to_csv(result_filename, sep='\t')
                            except:
                                open(result_filename, 'a').close()
                else:
                    result_filename2 = result_file.replace('.tsv', f'_k_{n_cluster}_eignesolver_{eigen_solver}_affinity_{re.sub("_", "", affinity)}_assignlabels_{re.sub("_", "", label_assigner)}_seed_{seed}.tsv')
                    if not os.path.exists(result_filename2):
                        try:
                            sc = SpectralClustering(n_clusters=n_cluster, eigen_solver=eigen_solver, affinity=affinity, assign_labels=label_assigner).fit(df)
                            result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': sc.labels_}), input_df=df)
                            result_k.to_csv(result_filename2, sep='\t')
                        except:
                            open(result_filename2, 'a').close()

