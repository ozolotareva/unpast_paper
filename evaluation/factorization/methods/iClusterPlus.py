# Run this with the environment "encore3". It has R 3.6 installed, such that it can run icluster

import itertools
import subprocess
import pandas as pd
import os
import numpy as np
from pathlib import Path
import pathlib

def get_parameters():
    files = [
        '/local/DESMOND2_data_simulated/simulated/A/A.n_genes=5,m=4,std=1,overlap=no.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/A/A.n_genes=50,m=4,std=1,overlap=no.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/B/B.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/B/B.n_genes=50,m=4,std=1,overlap=yes.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/B/B.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/C/C.n_genes=5,m=4,std=1,overlap=yes.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/C/C.n_genes=50,m=4,std=1,overlap=yes.exprs_z.tsv',
        '/local/DESMOND2_data_simulated/simulated/C/C.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv',
        
        
    ]
    n_cluster_list = list(range(4, 11, 2))
    type_list = ["gaussian","binomial","poisson","multinomial"]
    burnin_list = np.linspace(10, 500, 5)
    draw_list = np.linspace(10, 500, 5)
    lambda_list = np.linspace(1, 10, 5)
    lambda_list = np.append(lambda_list, 'null')
    maxiter_list = [10000]
    sdev_list = np.linspace(0.001, 0.2, 5)
    eps_list = np.linspace(1.0e-3, 1.0e-5, 5)

    combinations = []
    for file in files:
        for n_cluster in n_cluster_list:
            for type in type_list:
                for n_burnin in burnin_list:
                    for n_draw in draw_list:
                        for lambda_n in lambda_list:
                            for maxiter in maxiter_list:
                                for sdev in sdev_list:
                                    for eps in eps_list:
                                        combinations.append({
                                            'exprs_file':file,
                                            'lambda_n':lambda_n,
                                            'n_cluster':n_cluster,
                                            'lambda_scale':1, 
                                            'iter_max':maxiter,
                                            'eps':eps, 
                                            'type':type, 
                                            'burnin_n':n_burnin, 
                                            'draw_n':n_draw, 
                                            'sdev':sdev
                                        })
    return combinations

def format_output(n_cluster):
    df_cluster = pd.read_csv('iClusterPlusCluster.csv')
    cluster = {}
    for i in range(1, n_cluster+1):
        df_sub = df_cluster[df_cluster['x']==i]
        cluster[i] = {
            'samples': set(df_sub.index),
            'n_samples': len(df_sub.index)
        }
    os.remove('iClusterPlusCluster.csv')
    return pd.DataFrame(cluster).T


def read_runtime():
    runtime_string = open('iClusterPlusRuntime.txt', 'r').read().strip()
    runtime = float(runtime_string)
    os.remove('iClusterPlusRuntime.txt')
    return runtime

def run(exprs_file, lambda_n='null', n_cluster=4, lambda_scale=1, iter_max=20, eps=1.0e-4, type='gaussian', burnin_n=100, draw_n=200, sdev=0.05):
    # Given K, the number of cluster is K+1.
    n_cluster -= 1

    directory = "/iClusterPlusCluster/{}/{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
        exprs_file.split('/')[-1],
        lambda_n,
        n_cluster,
        lambda_scale,
        iter_max,
        eps,
        type,
        burnin_n,
        draw_n,
        sdev
    )
    directory = str(pathlib.Path().resolve()) + directory

    Path(directory).mkdir(parents=True, exist_ok=True)
    
    # this saves the result to a file
    # time is measured inside the R script
    print(f'Rscript ./methods/iClusterPlus.R {exprs_file} {lambda_n} {n_cluster} {lambda_scale} {iter_max} {eps} {type} {burnin_n} {draw_n} {sdev} {directory}')
    subprocess.Popen(fr'/usr/bin/Rscript ./methods/iClusterPlus.R {exprs_file} {lambda_n} {n_cluster} {lambda_scale} {iter_max} {eps} {type} {burnin_n} {draw_n} {sdev} {directory}', shell=True).wait()
    return #format_output(n_cluster), read_runtime()




