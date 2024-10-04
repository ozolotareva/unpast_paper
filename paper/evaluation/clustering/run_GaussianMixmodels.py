import sys
import os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.mixture import GaussianMixture
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
cov_types = ["spherical", "diag", "tied", "full"]
methods = ["kmeans", "random_from_data", "k-means++", "random"]

for k in tqdm(range(2, 21)):
    for cov_type in cov_types:
        for init_method in methods:
            # print(f'_covtype_{re.sub("_", "", cov_type)}_initmethod_{re.sub("_", "", init_method)}_k_{k}_seed_{seed}')
            result_filename = result_file.replace('.tsv', f'_covtype_{re.sub("_", "", cov_type)}_initmethod_{re.sub("_", "", init_method)}_k_{k}_seed_{seed}.tsv')
            if not os.path.exists(result_filename):
                try:
                    gmm = GaussianMixture(n_components=k, covariance_type=cov_type, init_params=init_method, random_state=seed)
                    result_k = reformat_cluster_results(
                        clusters=pd.DataFrame({'sample': df.index, 'label': gmm.fit_predict(df)}), input_df=df)
                    # print(result_k)
                    result_k.to_csv(result_filename, sep='\t')
                except:
                    # pass
                    open(result_filename, 'a').close()
