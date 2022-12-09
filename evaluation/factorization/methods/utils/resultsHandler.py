import os
import pandas as pd
from utils.eval import find_best_matches
from pathlib import Path
import numpy as np

OUTPUT_CLUSTER = 'result.csv'
OUTPUT_RUNTIME = 'runtime.txt'
OUTPUT_SAMPLES = 'samples.txt'

"""Something went wrong with some of the result files, maybe it was due to the scikit learn update."""
def _fix_result(df_result_broken):
    assert all(df_result_broken['samples'].map(lambda x: '{' in x or ('(' in x and 's' in x and 'e' in x and 't' in x and ')' in x)))
    df_result_broken['samples'] = df_result_broken['samples'].map(lambda x: set(x.replace(',', '').replace('{', '').replace('}', '').replace("'", '').split(' ')))
    df_result_broken['samples'] = df_result_broken['samples'].map(lambda x: set() if len(x) and list(x)[0] == 'set()' else x)
    return df_result_broken

def write_runtime(output_path, runtime):
    with open(os.path.join(output_path, OUTPUT_RUNTIME), "w") as text_file:
        text_file.write(f"{runtime}")

def write_samples(output_path, samples):
    samples = ','.join(list(samples))
    with open(os.path.join(output_path, OUTPUT_SAMPLES), "w") as text_file:
        text_file.write(f"{samples}")
        
def write_output(output_path, df_data):
    # df_data['samples'] = df_data['samples'].map(lambda x: ','.join(str(x)))
    df_data.to_csv(os.path.join(output_path, OUTPUT_CLUSTER))

def save(df_result, runtime, output_path):
    write_output(output_path, df_result)
    write_runtime(output_path, runtime)
    print(f'Saved {output_path}.')

def create_or_get_result_folder(output):
    Path(output).mkdir(parents=True, exist_ok=True)
    return os.path.isfile(os.path.join(output, OUTPUT_CLUSTER))

def read_result(output_path):
    df_data = pd.read_csv(os.path.join(output_path, OUTPUT_CLUSTER), index_col=0)
    # fucked up saving data, some contain the string 'set()' instead of NaN or an actual set {...} 
    try:
       df_data = _fix_result(df_data)
    except Exception:
        df_data['samples'] = df_data['samples'].map(lambda x: set(x.split(',')) if not isinstance(x, float) else set())
    return df_data

def read_runtime(output_path):
    with open(os.path.join(output_path, OUTPUT_RUNTIME), 'r') as f:
        return float(f.read())

def read_samples(output_path):
    with open(os.path.join(output_path, OUTPUT_SAMPLES), 'r') as f:
        samples = f.read().split(',')
        return samples
    
def calc_performance(found_clusters, all_samples, ground_truth_file,match_unique=True):
    
    ground_truth = pd.read_csv(ground_truth_file,sep ="\t",index_col=0)
    ground_truth["samples"] = ground_truth["samples"].apply(lambda x: set(x.split(" ")))
    if "genes" in ground_truth.columns.values:
        ground_truth["genes"] = ground_truth["genes"].apply(lambda x: set(x.split(" ")))
        
    # prepare a dict with sample groups corresponding to known bicluster
    known_groups = {}
    known_gsets = {}
    for group in ground_truth.index.values:
        known_groups[group] = ground_truth.loc[group,"samples"]
        known_gsets[group] = ground_truth.loc[group,"genes"]
    
    if found_clusters is None:
        performance= {}
        performance["total"] = 0
        best_matches = None
    else:
        best_matches = find_best_matches(found_clusters,known_groups,all_samples,
                                         FDR=0.05,verbose = False,match_unique=match_unique)
        J_total = best_matches["J_weighted"].sum()
        performance = best_matches.loc[:,["J"]].to_dict()["J"]
        # renaming subtype-specific performances to "performance_"+subt
        subtypes = list(performance.keys())
        for subt in subtypes:
            performance["performance_"+str(subt)] = performance.pop(subt)
        performance["overall_performance"] = J_total
    return performance, best_matches

def score_simulated_result(found_clusters, known_groups, all_samples):
    if not result['n_samples'].sum():
        return 0
    try:
        best_matches = find_best_matches(found_clusters, known_groups, all_samples, FDR=0.05)
        score = best_matches["J_weighted"].sum()
    except ZeroDivisionError:
        score = 0
    return score

def get_simulated_ground_truth(ground_truth_file):
    # read ground truth from file
    ground_truth = pd.read_csv(ground_truth_file,sep ="\t",index_col=0)
    ground_truth["samples"] = ground_truth["samples"].apply(lambda x: set(x.split(" ")))
    if "genes" in ground_truth.columns.values:
        ground_truth["genes"] = ground_truth["genes"].apply(lambda x: set(x.split(" ")))

    # prepare a dict with sample groups corresponding to known bicluster
    known_groups = {}
    for group in ground_truth.index.values:
        known_groups[group] = ground_truth.loc[group,"samples"]

    return known_groups

def evaluate_simulated(output_path, ground_truth_file, **_):
    result = read_result(output_path)
    samples = read_samples(output_path)
    # known_groups = get_simulated_ground_truth(ground_truth_file)
    # score = score_simulated_result(result, known_groups, samples)
    performance, best_matches = calc_performance(result, samples, ground_truth_file)
    runtime = read_runtime(output_path)
    return performance, runtime
