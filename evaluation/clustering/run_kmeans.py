import sys, os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.cluster import KMeans
from tqdm import tqdm
import numpy as np

# import matplotlib.pyplot as plt
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
# input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
df = pd.read_csv(input_file, sep='\t', index_col=0).T

cs = []
for k in tqdm(range(1, 21)):
    est_kmeans = KMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=20, random_state=0)
    est_kmeans.fit(df)
    result_k = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': est_kmeans.labels_}), input_df=df)
    result_k.to_csv(result_file.replace('.tsv', f'_k_{k}_seed_{seed}.tsv'), sep='\t')
    cs.append(est_kmeans.inertia_)

'''
#optional plotting
plt.plot(range(1, 21), cs, '-o')
plt.title('Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('CS')
plt.show()

kn = KneeLocator(range(len(cs)), cs, curve='convex', direction='decreasing')
if kn.knee:
    print('Found a knee')
    est_kmeans = KMeans(n_clusters=kn.knee, init='k-means++', max_iter=300, n_init=20, random_state=0)
    est_kmeans.fit(df)
    result_knee = reformat_cluster_results(clusters=pd.DataFrame({'sample': df.index, 'label': est_kmeans.labels_}), input_df=df)
    result_knee.to_csv(result_file.replace('.tsv', f'_knee_{kn.knee}.tsv'), sep='\t')
'''