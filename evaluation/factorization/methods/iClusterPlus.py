# Run this with the environment "encore3". It has R 3.6 installed, such that it can run icluster

import subprocess
import pandas as pd
import os
from .settings import RANDOM_STATES, CLUSTER_RANGE, ICLUSTER_TYPES, ICLUSTER_BURNIN_NS, ICLUSTER_DRAW_NS, ICLUSTER_EPS_LIST, ICLUSTER_SDEVS, ICLUSTER_LAMBDA_NS, ICLUSTER_LAMBDA_SCALES, ICLUSTER_MAXITERS
from .utils.miscellaneous import run_method
from .utils import interpret_results, resultsHandler

def generate_arg_list(exprs_file, output_folder, ground_truth_file, cluster_range=CLUSTER_RANGE):
    arguments = []
    for seed in RANDOM_STATES:
        for n_cluster in cluster_range:
            for type in ICLUSTER_TYPES:
                for n_burnin in ICLUSTER_BURNIN_NS:
                    for n_draw in ICLUSTER_DRAW_NS:
                        for lambda_n in ICLUSTER_LAMBDA_NS:
                            for maxiter in ICLUSTER_MAXITERS:
                                for sdev in ICLUSTER_SDEVS:
                                    for eps in ICLUSTER_EPS_LIST:
                                        for lambda_scale in ICLUSTER_LAMBDA_SCALES:
                                            output_path = os.path.join(output_folder, 
                                                f'random_state={seed}', 
                                                f'n_cluster={n_cluster}',
                                                f'type={type}',
                                                f'n_burnin={n_burnin}',
                                                f'n_draw={n_draw}',
                                                f'lambda_n={lambda_n}',
                                                f'maxiter={maxiter}',
                                                f'sdev={sdev}',
                                                f'eps={eps}',
                                                f'lambda_scale={lambda_scale}',
                                            )

                                            arguments.append({
                                                'exprs_file': exprs_file,
                                                'output_path': output_path,
                                                'ground_truth_file': ground_truth_file,
                                                'lambda_n':lambda_n,
                                                'n_cluster':n_cluster,
                                                'lambda_scale':lambda_scale, 
                                                'iter_max':maxiter,
                                                'eps':eps, 
                                                'random_state':seed, 
                                                'type':type, 
                                                'burnin_n':n_burnin, 
                                                'draw_n':n_draw, 
                                                'sdev':sdev
                                            })
    return arguments

def format_output(output_path, n_cluster, labels):
    df_i_cluster = pd.read_csv(os.path.join(output_path, 'iclusterplus_result.csv'))
    cluster = {}
    for i in range(1, n_cluster+1):
        df_sub = df_i_cluster[df_i_cluster['x']==i]
        cluster[i] = {
            'samples': set([labels[x-1] for x in df_sub.index]),
            'n_samples': len(df_sub.index)
        }
    # os.remove(os.path.join(output_path, 'iclusterplus_result.csv'))
    return pd.DataFrame(cluster).T

def read_runtime(output_path):
    runtime_string = open(os.path.join(output_path, 'iclusterplus_runtime.txt'), 'r').read().strip()
    runtime = float(runtime_string)
    # os.remove(os.path.join(output_path, 'iclusterplus_runtime.txt'))
    return runtime

def execute_algorithm(exprs_file, lambda_n, n_cluster, lambda_scale, iter_max, eps, type, burnin_n, draw_n, sdev, output_path, random_state, **_):
    # Given K, the number of cluster is K+1.
    k = n_cluster - 1
    # this saves the result to a file
    # time is measured inside the R script
    subprocess.Popen(fr'/usr/bin/Rscript ./methods/iClusterPlus.R {exprs_file} {lambda_n} {k} {lambda_scale} {iter_max} {eps} {type} {burnin_n} {draw_n} {sdev} {output_path} {random_state}', shell=True).wait()
    return (output_path, n_cluster), read_runtime(output_path)

def run_simulated(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Skipping because result exists:', args["output_path"])
        return resultsHandler.read_result(args["output_path"])
    df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
    result, runtime = run_method(execute_algorithm, args)
    result = format_output(result[0], result[1], df_exprs.index)
    # save results
    resultsHandler.save(result, runtime, args["output_path"])
    resultsHandler.write_samples(args["output_path"], df_exprs.index)
    return resultsHandler.read_result(args["output_path"])
    
def run_real(args, is_terminated=False):
    if is_terminated:
        try:
            result = resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
            return result
        except:
            return False, False
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Returning existing results:', args["output_path"])
    else:
        print('Executing iClusterPlus')
        result, runtime = run_method(execute_algorithm, args)
        df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
        result = format_output(result[0], result[1], df_exprs.index)
        resultsHandler.save(result, runtime, args["output_path"])
    return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])




