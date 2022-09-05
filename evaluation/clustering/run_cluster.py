import os
import subprocess
import sys
from eval_cluster_methods import run_eval
import glob

args = sys.argv
tool_name = args[1]
script = args[2]
expr_file = args[3]
truth_file = args[4]
result_file = args[5]
score_file = args[6]

'''
os.chdir('/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/')
tool_name = 'WGCNA'
script = './run_WGCNAkmeans.py'
expr_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.biclusters.tsv'
result_file ='/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/kmeans/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/kmeans/scores.tsv'


tool_name = 'WGCNAkmeans'
script = './run_WGCNAkmeans.py'
expr_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv'
truth_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.biclusters.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/WGCNAkmeans/clusters/B.n_genes=5,m=4,std=1,overlap=yes_run5.tsv'
score_file = '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/WGCNAkmeans/WGCNAkmeans_scores.txt'
'''


def get_command(tool_name, script_location, expr_file, out_file):
    command = []
    if tool_name in ['kmeans', 'WGCNAkmeans', 'WGCNAHC', 'HC']:
        command.append("python3")
        command.append(script_location)
        command.append(expr_file)
        command.append(out_file)
    return command


subprocess.check_call(get_command(tool_name, script, expr_file, result_file))
listing = glob.glob(result_file[:-4] + '*')

with open(score_file, 'a+') as fw:
    for filename in listing:
        # os.system(f'rm {filename}')
        j_weighted = run_eval(expr_file=expr_file, result_file=filename, ground_truth_file=truth_file)
        print(f"J_weighted for {tool_name} and {filename}: {j_weighted}")
        fw.write(f"{filename.split('/')[-1]}\t{str(j_weighted)}\n")
