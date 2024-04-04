import subprocess
import os
import pandas as pd
from .settings import RANDOM_STATES, CLUSTER_RANGE, MOFA2_FACTORS, MOFA2_ARD_WEIGHTS, MOFA2_SPIKESLAB_FACTORS, MOFA2_ARD_FACTORS, MOFA2_LIKELIHOODS, MOFA2_SPIKESLAB_WEIGHTS
from .utils.miscellaneous import run_method
from .utils import resultsHandler


def generate_arg_list(exprs_file, output_folder, ground_truth_file, cluster_range=CLUSTER_RANGE):
    arguments = []
    # MOFA2 - not transposed
    for m in RANDOM_STATES:
        for n_cluster in cluster_range:
            for n_factors in MOFA2_FACTORS:
                for ard_weights in MOFA2_ARD_WEIGHTS:
                    for ard_factors in MOFA2_ARD_FACTORS:
                        for likelihood in MOFA2_LIKELIHOODS:
                            for spikeslab_weights in MOFA2_SPIKESLAB_WEIGHTS:
                                for spikeslab_factors in MOFA2_SPIKESLAB_FACTORS:
                                    output_path = os.path.join(output_folder, 
                                        f'n_factors={n_factors}',
                                        f'n_cluster={n_cluster}',
                                        f'random_state={m}', 
                                        f'ard_weights={ard_weights}', 
                                        f'ard_factors={ard_factors}', 
                                        f'likelihood={likelihood}', 
                                        f'spikeslab_weights={spikeslab_weights}', 
                                        f'spikeslab_factors={spikeslab_factors}', 
                                        )

                                    args = {'exprs_file': exprs_file,
                                            'output_path': output_path,
                                            'ground_truth_file': ground_truth_file,
                                            'n_factors': n_factors,
                                            'n_cluster': n_cluster,
                                            'random_state': m,
                                            'ard_weights': ard_weights,
                                            'ard_factors': ard_factors,
                                            'likelihood': likelihood,
                                            'spikeslab_weights': spikeslab_weights,
                                            'spikeslab_factors': spikeslab_factors,
                                        }
                                    arguments.append(args)
    return arguments

def format_output(output_path, n_cluster):
    df_mofa2_cluster = pd.read_csv(os.path.join(output_path, 'mofa2_result.csv'))
    cluster = {}
    if n_cluster == 'all':
        n_cluster = len(df_mofa2_cluster['x'].unique())
    for i in range(1, n_cluster+1):
        df_sub = df_mofa2_cluster[df_mofa2_cluster['x']==i]
        cluster[i] = {
            'samples': set(df_sub.index),
            'n_samples': len(df_sub.index)
        }
    # os.remove(os.path.join(output_path, 'mofa2_result.csv'))
    # os.remove(os.path.join(output_path, 'mofa2_model.hdf5'))
    return pd.DataFrame(cluster).T

def read_runtime(output_path):
    runtime_string = open(os.path.join(output_path, 'mofa2_runtime.txt'), 'r').read().strip()
    runtime = float(runtime_string)
    os.remove(os.path.join(output_path, 'mofa2_runtime.txt'))
    return runtime

def execute_algorithm(exprs_file, n_factors, n_cluster, output_path, likelihood, spikeslab_factors, spikeslab_weights, ard_factors, ard_weights, random_state=101, **_):
    # this saves the result to a file
    # time is measured inside the R script

    subprocess.Popen(fr'/home/bba1401/anaconda3/envs/encore2/bin/Rscript ./methods/MOFA2.R {exprs_file} {n_factors} {n_cluster} {random_state} {output_path} {likelihood} {spikeslab_factors} {spikeslab_weights} {ard_factors} {ard_weights}', shell=True).wait()
    return format_output(output_path, n_cluster), read_runtime(output_path)

def run_simulated(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Skipping because result exists:', args["output_path"])
        return resultsHandler.read_result(args["output_path"])
    df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
    result, runtime = run_method(execute_algorithm, args)

    # save results
    resultsHandler.save(result, runtime, args["output_path"])
    resultsHandler.write_samples(args["output_path"], df_exprs.index)
    return resultsHandler.read_result(args["output_path"])

def run_real(args, is_terminated=False):
    if is_terminated:
        try:
            return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
        except Exception as e:
            print('Exception reading result file:', e)
            return False, False
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Returning existing results:', args["output_path"])
    else:
        result, runtime = run_method(execute_algorithm, args)
        resultsHandler.save(result, runtime, args["output_path"])
    return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])



    