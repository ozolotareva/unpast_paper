import sys
import os
from pathlib import Path
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
# os.chdir('/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering')

args = sys.argv
script = args[0]
input_file = args[1]
results_file = args[2]

'''
input_file =  'input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv''
results_file =  '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/evaluation/clustering/results/WGCNAkmeans/clusters/B.n_genes=5,m=4,std=1,overlap=yes_run5.tsv'
'''
# results_dir = '/Users/fernando/Documents/Research/DESMOND2/methods/results/results_WGCNAkmeans/'
# input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/simulated/A/example.tsv'
dimreduced_input = input_file.replace('.tsv', '_MEs.tsv')

if not Path(dimreduced_input).exists():
    print('Running R WGCNA script')
    os.system(f"Rscript --vanilla run_WGCNA_dimreduce.R {input_file}")
os.system(f"python3 run_HC.py {dimreduced_input} {results_file}")
