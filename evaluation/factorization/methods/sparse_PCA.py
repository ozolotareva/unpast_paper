from sklearn.decomposition import SparsePCA
import sklearn
import time
import os
import pandas as pd
from .settings import S_PCA_ALPHA, S_PCA_MAX_ITER, S_PCA_METHOD, S_PCA_RIDGE_ALPHA, S_PCA_TOL, RANDOM_STATES, CLUSTER_RANGE
from .utils.miscellaneous import run_method
from .utils import interpret_results, resultsHandler
import json
import numpy as np


def generate_arg_list(exprs_file, output_folder, ground_truth_file, cluster_range=CLUSTER_RANGE):
    arguments = []
    cluster_range = list(cluster_range)
    cluster_range.append(None)
    for m in RANDOM_STATES:
        for n_components in cluster_range:
            for alpha in S_PCA_ALPHA:
                for max_iter in S_PCA_MAX_ITER:
                    for method in S_PCA_METHOD:
                        for ridge_alpha in S_PCA_RIDGE_ALPHA:
                            for tol in S_PCA_TOL:
                                output_path = os.path.join(output_folder, 
                                    f'n_components={n_components}',
                                    f'random_state={m}', 
                                    f'alpha={alpha}', 
                                    f'ridge_alpha={ridge_alpha}', 
                                    f'max_iter={max_iter}',
                                    f'method={method}',
                                    f'tol={tol}',
                                    )

                                args = {'exprs_file': exprs_file,
                                        'output_path': output_path,
                                        'ground_truth_file': ground_truth_file,
                                        'n_components': n_components,
                                        'random_state': m,
                                        'alpha': alpha,
                                        'ridge_alpha': ridge_alpha,
                                        'max_iter': max_iter,
                                        'method': method,
                                        'tol': tol,
                                    }
                                arguments.append(args)
    return arguments

def execute_algorithm(exprs, n_components, alpha, ridge_alpha, max_iter, method, tol, output_path, random_state=101, **_):
    start = time.time()
    with sklearn.config_context(assume_finite=True):
        transformer = SparsePCA(n_components=n_components, alpha=alpha, ridge_alpha=ridge_alpha, max_iter=max_iter, method=method, tol=tol, random_state=random_state)
        exprs_transformed = transformer.fit_transform(exprs)
    result = interpret_results.format_sklearn_output(exprs_transformed, n_components, exprs.index)
    end = time.time()
    runtime = end-start
    
    df = pd.DataFrame(transformer.components_, columns=exprs.columns)
    df.to_csv(os.path.join(output_path, 'pca_components.tsv'), sep='\t')
    gene_clusters = {}
    for index, row in df.iterrows():
        gene_cluster = list(row[row.map(lambda x : np.abs(x) > 0)].index)
        gene_clusters[index] = gene_cluster

    df_gene_clusters = pd.DataFrame([{'gene_cluster': key, 'genes': values} for key, values in gene_clusters.items()])
    df_gene_clusters.to_csv(os.path.join(output_path, 'gene_clusters.csv'))
    
    result['genes'] = df_gene_clusters['genes']
    result['n_genes'] = df_gene_clusters['genes'].map(len)
    return result, runtime

def run_simulated(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Skipping because result exists:', args["output_path"])
        return resultsHandler.read_result(args["output_path"])
    args['exprs'] = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
    result, runtime = run_method(execute_algorithm, args)
    resultsHandler.write_samples(args["output_path"], args['exprs'].index)
    del args['exprs']

    # save results
    resultsHandler.save(result, runtime, args["output_path"])
    return resultsHandler.read_result(args["output_path"])


def run_real(args, is_terminated=False):
    if is_terminated:
        try:
            return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
        except:
            return False, False
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Returning existing results:', args["output_path"])
    else:
        args['exprs'] = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
        result, runtime = run_method(execute_algorithm, args)
        resultsHandler.save(result, runtime, args["output_path"])
        del args['exprs']
    return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
