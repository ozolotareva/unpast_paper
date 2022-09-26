import subprocess
import pandas as pd
import os


def format_output(n_cluster):
    df_mofa2_cluster = pd.read_csv('MOFA2Cluster.csv')
    cluster = {}
    for i in range(1, n_cluster+1):
        df_sub = df_mofa2_cluster[df_mofa2_cluster['x']==i]
        cluster[i] = {
            'samples': set(df_sub.index),
            'n_samples': len(df_sub.index)
        }
    os.remove('MOFA2Cluster.csv')
    return pd.DataFrame(cluster).T

def read_runtime():
    runtime_string = open('MOFA2Runtime.txt', 'r').read().strip()
    runtime = float(runtime_string)
    os.remove('MOFA2Runtime.txt')
    return runtime

def run(exprs_file, n_factors, n_cluster, seed):
    # this saves the result to a file   
    # time is measured inside the R script
    subprocess.Popen(fr'Rscript ./methods/MOFA2.R {exprs_file} {n_factors} {n_cluster} {seed}', shell=True).wait()
    return format_output(n_cluster), read_runtime()