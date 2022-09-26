import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import os
from functools import partial
from itertools import repeat
from multiprocessing import Pool, freeze_support

from utils.eval import find_best_matches

from methods import NMF, PCA, sparse_PCA, moCluster, MOFA2

# let methods find different numbers of clusters
CLUSTER_RANGE = range(4, 11)

scenarios = ['A', 'B', 'C']
bicluster_sizes = [5, 50, 500]


def test_all(scenario, N):
    # expression file
    try:
        exprs_file = f"/local/DESMOND2_data_simulated/simulated/{scenario}/{scenario}.n_genes={N},m=4,std=1,overlap=yes.exprs_z.tsv"
        exprs = pd.read_csv(exprs_file,sep = "\t",index_col=0)
    except:
        exprs_file = f"/local/DESMOND2_data_simulated/simulated/{scenario}/{scenario}.n_genes={N},m=4,std=1,overlap=no.exprs_z.tsv"
        exprs = pd.read_csv(exprs_file,sep = "\t",index_col=0)
    

    # read ground truth from file
    try:
        ground_truth_file = f"/local/DESMOND2_data_simulated/simulated/{scenario}/{scenario}.n_genes={N},m=4,std=1,overlap=yes.biclusters.tsv"
        ground_truth = pd.read_csv(ground_truth_file,sep ="\t",index_col=0)
    except:
        ground_truth_file = f"/local/DESMOND2_data_simulated/simulated/{scenario}/{scenario}.n_genes={N},m=4,std=1,overlap=no.biclusters.tsv"
        ground_truth = pd.read_csv(ground_truth_file,sep ="\t",index_col=0)
    
    ground_truth["samples"] = ground_truth["samples"].apply(lambda x: set(x.split(" ")))
    if "genes" in ground_truth.columns.values:
        ground_truth["genes"] = ground_truth["genes"].apply(lambda x: set(x.split(" ")))

    # prepare a dict with sample groups corresponding to known bicluster
    known_groups = {}
    for group in ground_truth.index.values:
        known_groups[group] = ground_truth.loc[group,"samples"]

    SAMPLE_SIZES = len(exprs.columns)

    def test_method(method, args, known_groups):
        result, runtime = method(**args)
        try:
            best_matches = find_best_matches(result,known_groups, SAMPLE_SIZES, FDR=0.05)
            score = best_matches["J_weighted"].sum()
        except ZeroDivisionError:
            score = 0
        return score, result, runtime

    # collect all results from the methods
    evaluations = {}
    evaluations_all = {}
    
    # NMF - not transposed
    exprs_NMF = NMF.preprocess_exprs(exprs)
    scored_results = []
    for m in range(5):
        for k in CLUSTER_RANGE:
            print('Running NMF with k:', k, '...')
            args = {'exprs':exprs_NMF, 'k':k, 'random_state':m, 'transposed':False}
            test_result = test_method(NMF.run, args, known_groups)
            scored_results.append((k, *test_result, args))
            evaluations_all[f'NMF_H_randomState={m}_k={k}'] = (k, *test_result)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'NMF_H'] = best_result
    
    # NMF - transposed
    exprs_NMF = NMF.preprocess_exprs(exprs)
    scored_results = []
    for m in range(5):
        for k in CLUSTER_RANGE:
            print('Running NMF with k:', k, '...')
            args = {'exprs':exprs_NMF.T, 'k':k, 'random_state':m, 'transposed':True}
            test_result = test_method(NMF.run, args, known_groups)
            scored_results.append((k, *test_result, args))
            evaluations_all[f'NMF_W_randomState={m}_k={k}'] = (k, *test_result)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'NMF_W'] = best_result
    
    # PCA
    scored_results = []
    for n in CLUSTER_RANGE:
        print('Running PCA with n components:', n, '...')
        args = {'exprs':exprs.T, 'n':n}
        test_result = test_method(PCA.run, args, known_groups)
        scored_results.append((n, *test_result))
        evaluations_all[f'PCA_{n}'] = (n, *test_result, args)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'PCA'] = best_result
    
    # Sparse PCA
    scored_results = []
    for n in CLUSTER_RANGE:
        print('Running sparse_PCA with n components:', n, '...')
        args = {'exprs':exprs.T, 'n':n}
        test_result = test_method(sparse_PCA.run, args, known_groups)
        scored_results.append((n, *test_result, args))
        evaluations_all[f'sparse_PCA_{n}'] = (n, *test_result)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'sparse_PCA'] = best_result
    
    # moCluster
    scored_results = []
    for m in range(1, 11): # has no impact
        print('m:', m)
        for n in CLUSTER_RANGE:
            print('Running moCluster with n components:', n, '...')
            args = {'exprs_file':exprs_file, 'n_dimensions':m, 'n_cluster':n}
            test_result = test_method(moCluster.run, args, known_groups)
            scored_results.append((n, *test_result, args))
            evaluations_all[f'moCluster_m={m}_cluster={n}'] = (n, *test_result)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'moCluster'] = best_result
    
    # MOFA2
    scored_results = []
    for s in range(5):
        print('s:', s)
        for m in range(1, 11):
            print('m:', m)
            for n in CLUSTER_RANGE:
                print('Running MOFA2 with n components:', n, '...')
                args = {'exprs_file':exprs_file, 'n_factors':m, 'n_cluster':n, 'seed':s }
                test_result = test_method(MOFA2.run, args, known_groups)
                scored_results.append((n, *test_result, args))
                evaluations_all[f'MOFA2_factors={m}_cluster={n}_seed={s}'] = (n, *test_result)
    best_result = sorted(scored_results, key = lambda y: y[1], reverse=True)[0]
    evaluations[f'MOFA2'] = best_result
    
    # save data
    df_evaluations = pd.DataFrame(evaluations).T
    df_evaluations.columns = ['#cluster', 'Jaccard', 'cluster', 'time (s)', 'args']
    df_evaluations = df_evaluations.sort_values(by='Jaccard', ascending=False)
    df_evaluations.to_csv(f'simulated/evaluation_{scenario}_{N}.csv')
    
    df_evaluations_all = pd.DataFrame(evaluations_all).T
    df_evaluations_all.columns = ['#cluster', 'Jaccard', 'cluster', 'time (s)']
    df_evaluations_all = df_evaluations_all.sort_values(by='Jaccard', ascending=False)
    df_evaluations_all.to_csv(f'simulated/evaluation_all_{scenario}_{N}.csv')
        

if __name__ == '__main__':
    params = []
    for scenario in scenarios:
        for N in bicluster_sizes:
            params.append([scenario, N])
    print(params)
    freeze_support()
    with Pool(processes=4) as pool:
        pool.starmap(test_all, params)