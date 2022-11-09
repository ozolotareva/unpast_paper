import sys
import os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.cluster import BisectingKMeans
import re
import numpy as np
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


'''
input_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/clusters/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'
'''

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

df = pd.read_csv(input_file, sep='\t', index_col=0).T

methods = ["k-means++", "random"]
algorithms = ['lloyd', 'elkan']
bisecting_strategies = ['biggest_inertia', 'largest_cluster']

for k in tqdm(range(2, 21)):
    for algorithm in algorithms:
        for bisecting_strategy in bisecting_strategies:
            for initmethod in methods:
                # print(f'_initmethod_{re.sub("_", "", initmethod)}_algorithm_{re.sub("_", "", algorithm)}_bisectingstrategy_{re.sub("_", "", bisecting_strategy)}_k_{k}_seed_{seed}')
                result_filename = result_file.replace('.tsv', f'_initmethod_{re.sub("_", "", initmethod)}_algorithm_{re.sub("_", "", algorithm)}_bisectingstrategy_{re.sub("_", "", bisecting_strategy)}_k_{k}_seed_{seed}.tsv')
                if not os.path.exists(result_filename):
                    try:
                        bikmeans = BisectingKMeans(n_clusters=k, init=initmethod, bisecting_strategy=bisecting_strategy, algorithm=algorithm, random_state=seed)
                        result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': bikmeans.fit_predict(df)}), input_df=df)
                        # print(result_k)
                        result_k.to_csv(result_filename, sep='\t')

                    except:
                        open(result_filename, 'a').close()
