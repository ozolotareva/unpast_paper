import subprocess
import os
import pandas as pd
from .settings import RANDOM_STATES, CLUSTER_RANGE, MOCLUSTER_RANGE
from .utils.miscellaneous import run_method
from .utils import interpret_results, resultsHandler


def generate_arg_list(exprs_file, output_folder, ground_truth_file, cluster_range=CLUSTER_RANGE):
    arguments = []
    for m in RANDOM_STATES:
        for n_cluster in cluster_range:
            for n_dimensions in MOCLUSTER_RANGE:
                output_path = os.path.join(output_folder, 
                    f'n_dimensions={n_dimensions}',
                    f'n_cluster={n_cluster}',
                    f'random_state={m}', 
                    )

                args = {'exprs_file': exprs_file,
                        'output_path': output_path,
                        'ground_truth_file': ground_truth_file,
                        'n_dimensions': n_dimensions,
                        'n_cluster': n_cluster,
                        'random_state': m,
                    }
                arguments.append(args)
    return arguments

def format_output(output_path, n_cluster):
    df_mocluster_cluster = pd.read_csv(os.path.join(output_path, 'mocluster_result.csv'))
    cluster = {}
    for i in range(1, n_cluster+1):
        df_sub = df_mocluster_cluster[df_mocluster_cluster['x']==i]
        cluster[i] = {
            'samples': set(df_sub.index),
            'n_samples': len(df_sub.index)
        }
    os.remove(os.path.join(output_path, 'mocluster_result.csv'))
    return pd.DataFrame(cluster).T

def read_runtime(output_path):
    runtime_string = open(os.path.join(output_path, 'mocluster_runtime.txt'), 'r').read().strip()
    runtime = float(runtime_string)
    os.remove(os.path.join(output_path, 'mocluster_runtime.txt'))
    return runtime

def execute_algorithm(exprs_file, n_dimensions, n_cluster, random_state, output_path, **_):
    # this saves the result to a file
    # time is measured inside the R script
    subprocess.Popen(fr'Rscript ./methods/moCluster.R {exprs_file} {n_dimensions} {n_cluster} {random_state} {output_path}', shell=True).wait()
    return format_output(output_path, n_cluster), read_runtime(output_path)

def run_simulated(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Skipping because result exists:', args["output_path"])
        return
    df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0).T
    result, runtime = run_method(execute_algorithm, args)

    # save results
    resultsHandler.save(result, runtime, args["output_path"])
    resultsHandler.write_samples(args["output_path"], df_exprs.index)

def run_real(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Returning existing results:', args["output_path"])
    else:
        result, runtime = run_method(execute_algorithm, args)
        resultsHandler.save(result, runtime, args["output_path"])
    return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])

