import subprocess
import pandas as pd
import os


def format_output(n_cluster):
    df_mocluster_cluster = pd.read_csv('moClusterCluster.csv')
    cluster = {}
    for i in range(1, n_cluster+1):
        df_sub = df_mocluster_cluster[df_mocluster_cluster['x']==i]
        cluster[i] = {
            'samples': set(df_sub.index),
            'n_samples': len(df_sub.index)
        }
    os.remove('moClusterCluster.csv')
    return pd.DataFrame(cluster).T

def read_runtime():
    runtime_string = open('moClusterRuntime.txt', 'r').read().strip()
    runtime = float(runtime_string)
    os.remove('moClusterRuntime.txt')
    return runtime

def run(exprs_file, n_dimensions, n_cluster):
    # this saves the result to a file   
    # time is measured inside the R script
    subprocess.Popen(fr'Rscript ./methods/moCluster.R {exprs_file} {n_dimensions} {n_cluster}', shell=True).wait()
    return format_output(n_cluster), read_runtime()