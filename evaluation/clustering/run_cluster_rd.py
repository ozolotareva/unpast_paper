import os
import subprocess
import sys
# from evaluation.clustering.eval_cluster_methods import run_evalRD
# from eval_cluster_methods import run_evalRD
import glob
from tqdm import tqdm
import pandas as pd

args = sys.argv
tool_name = args[1]
script = args[2]
expr_file = args[3]
annot_file = args[4]
subtype_file = args[5]
result_file = args[6]
score_file = args[7]

'''
tool_name = 'kmeans'
script =  '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/run_kmeans.py'
expr_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/preprocessed_v6//METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv'
annot_file =   '/Users/fernando/Documents/Research/DESMOND2_data_simulated/preprocessed_v6/METABRIC_1904.annotation_v6.tsv'
subtype_file ='/Users/fernando/Documents/Research/DESMOND2_data_simulated/preprocessed_v6//METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv'
result_file ='/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/resultsRD/kmeans/clusters/METABRIC_run4.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/resultsRD/kmeans/kmeans_scoresRD.txt'
'''

if not os.path.exists(score_file):
    f = open(score_file, 'w+')
    print('File is created')


def get_command(tool_name, script_location, expr_file, out_file):
    command = []
    if tool_name in ['kmeans', 'WGCNAkmeans', 'WGCNAHC', 'HC', 'AffinityPropagation', 'MeanShift', 'Spectral', 'AgglomerativeClustering', 'DBSCAN', 'OPTICS', 'BIRCH', 'bikmeans', 'GMM']:
        command.append("python3")
        command.append(script_location)
        command.append(expr_file)
        command.append(out_file)
    return command


subprocess.check_call(get_command(tool_name, script, expr_file, result_file))

'''
listing = glob.glob(result_file[:-4] + '*')

subt_t = {}
# with open(score_file, 'a+') as fw:
for filename in tqdm(listing):
    # print(filename)
    # os.system(f'rm {filename}')
    # j_weighted = run_eval(expr_file=expr_file, result_file=filename, ground_truth_file=truth_file)
    j_s = run_evalRD(result_filepath=filename, subtypes_filepath=subtype_file, annotation_filepath=annot_file, expr_filepath=expr_file)
    keyname = filename.split('/')[-1]
    subt_t[keyname] = j_s['J']

results = pd.DataFrame.from_records(subt_t).T

if os.stat(score_file).st_size == 0:
    print('File is empty')
    results.to_csv(score_file, mode='a', index=True, header=True, sep='\t')

else:
    print('File is not empty')
    results.to_csv(score_file, mode='a', index=True, header=False, sep='\t')
'''
# print(f"J_weighted for {tool_name} and {filename}: {j_weighted}")
# fw.write(f"{filename.split('/')[-1]}\t{str(j_weighted)}\n")
