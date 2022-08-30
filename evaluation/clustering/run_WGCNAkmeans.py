import sys, os
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
#os.chdir('/Users/fernando/Documents/Research/DESMOND2/methods')

args = sys.argv
input_file = args[1]
results_dir = args[2]

# results_dir = '/Users/fernando/Documents/Research/DESMOND2/methods/results/results_WGCNAkmeans/'
# input_file = '/Users/fernando/Documents/Research/DESMOND2/datasets/DESMOND2_data_simulated/simulated/A/example.tsv'
dimreduced_input = input_file.replace('.tsv', '_MEs.tsv')

print('Running R WGCNA script')
os.system(f"Rscript --vanilla run_WGCNA.R {input_file}")
os.system(f"python3 run_kmeans.py {dimreduced_input} {results_dir}")
