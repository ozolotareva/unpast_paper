import os
import subprocess
import sys
import glob
from tqdm import tqdm
sys.path.insert(1, '/home/bbb1417/DESMOND2_benchmarking/DESMOND2/utils/')
# from eval2 import run_eval
import re

args = sys.argv
tool_name = args[1]
script = args[2]
expr_file = args[3]
truth_file = args[4]
result_file = args[5]
score_file = args[6]

'''
os.chdir('/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering')
tool_name = 'WGCNAkmeans'
script = './run_WGCNAkmeans.py'
expr_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/WGCNAkmeans/clusters/B.n_genes=5,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/WGCNAkmeans/WGCNAkmeans_scores.txt'

script = './run_HC.py'
expr_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/A.n_genes=5,m=4,std=1,overlap=no.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/A.n_genes=5,m=4,std=1,overlap=no.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/HC/clusters/A.n_genes=5,m=4,std=1,overlap=no_run4.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/HC/HC_scores.txt'

tool_name = 'Spectral'
script = './run_Spectral.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/Spectral/clusters/B.n_genes=5,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/Spectral/Spectral_scores.txt'

tool_name = 'MeanShift'
script = './run_MeanShift.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/MeanShift/clusters/B.n_genes=5,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/MeanShift/Spectral_scores.txt'

tool_name = 'AffinityPropagation'
script = './run_AffinityPropagation.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/AffinityPropagation/clusters/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/AffinityPropagation/AffinityPropagation_scores.txt'


tool_name = 'AgglomerativeClustering'
script = './run_AgglomerativeClustering.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/AgglomerativeClustering/clusters/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/AgglomerativeClustering/AgglomerativeClustering_scores.txt'


tool_name = 'DBSCAN'
script = './run_DBSCAN.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/DBSCAN/clusters/B.n_genes=500,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/DBSCAN/DBSCAN_scores.txt'


tool_name = 'OPTICS'
script = './run_OPTICS.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/OPTICS/clusters/B.n_genes=500,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/OPTICS/OPTICS_scores.txt'

tool_name = 'BIRCH'
script = './run_BIRCH.py'
expr_file =  '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/BIRCH/clusters/B.n_genes=500,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/BIRCH/BIRCH_scores.txt'


'''


def get_command(tool_name, script_location, expr_file, out_file):
    command = []
    if tool_name in ['kmeans', 'WGCNAkmeans', 'WGCNAHC', 'HC', 'AffinityPropagation', 'MeanShift', 'Spectral', 'AgglomerativeClustering', 'DBSCAN', 'OPTICS', 'BIRCH', 'bikmeans', 'GMM']:
        command.append("python3")
        command.append(script_location)
        command.append(expr_file)
        command.append(out_file)
    return command


def pairwise(t):
    it = iter(t)
    return zip(it, it)


def parameters_to_string(parameters):
    pairs = list(pairwise(parameters))
    return ';'.join(['='.join(pair) for pair in pairs])


subprocess.check_call(get_command(tool_name, script, expr_file, result_file))

'''
listing = glob.glob(result_file[:-4] + '*')

with open(score_file, 'a+') as fw:
    for filename in listing:
        if 'k_1.' not in filename:
            # print(expression_files[keyval], bicluster_files[keyval], file)
            #print(filename)
            res = run_eval(expr_file=expr_file, ground_truth_file=truth_file, result_file=filename)
            # fw.write(f"{file}\t{sumjaccard}\n")
            if res is not None:
                scenario = filename.split('/')[-1].split('.')[0]
                gsize = re.search('.n_genes=(.*),m', filename).group(1)
                run = filename.split('_')[2][3:]
                parameters = filename[:-4].split('_')[3:]
                parameters = parameters_to_string(parameters)
                #print((f"{scenario}\t{gsize}\t{run}\t{method}\t{parameters}\t{res['performance_0.5']}\t{res['performance_0.25']}\t{res['performance_0.1']}\t{res['performance_0.05']}\t{res['overall_performance']}"))

                fw.write(f"{scenario}\t{gsize}\t{run}\t{tool_name}\t{parameters}\t{res['performance_0.5']}\t{res['performance_0.25']}\t{res['performance_0.1']}\t{res['performance_0.05']}\t{res['overall_performance']}\n")
                # b = [[file, sumjaccard]]
                # a = np.vstack((a,np.asarray(b,object)))
'''