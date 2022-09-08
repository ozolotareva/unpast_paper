import os

import pandas as pd
from utils.eval import find_best_matches


def read_results(tool_name, result_file):
    if tool_name in ['isa2', 'fabia', 'qubic']:
        result = pd.read_csv(result_file, names=['samples', 'genes'], sep="\t")
        result["samples"] = result["samples"].apply(lambda x: set(x.strip().split(" ")))
        result["genes"] = result["genes"].apply(lambda x: set(x.strip().split(" ")))
        return result
    elif tool_name == 'debi':
        samples = []
        genes = []
        line_nr = -1
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                line_nr = (line_nr + 1) % 3
                if line_nr == 1:
                    genes.append(set(line.strip().split(" ")))
                elif line_nr == 2:
                    samples.append(set(line.strip().split(" ")))
        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'qubic2':
        samples = []
        genes = []
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith(" Genes"):
                    genes.append(set(line.strip().split(": ")[1].split(" ")))
                elif line.startswith(" Conds"):
                    samples.append(set(line.strip().split(": ")[1].split(" ")))

        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'dataframe':
        result = pd.read_csv(result_file, delimiter="\t", index_col=0)
        result["samples"] = result["samples"].apply(eval)
        result["genes"] = result["genes"].apply(eval)
        return result


def run_eval(tool_name, expr_file, ground_truth_file, result_file):
    # expression file
    if tool_name == 'qubic2':
        result_file = result_file+".blocks"

    exprs = pd.read_csv(expr_file, sep="\t", index_col=0, header=0)
    samples = list(exprs.columns)
    # read ground truth from file
    ground_truth = pd.read_csv(ground_truth_file, sep="\t", index_col=0)
    ground_truth["samples"] = ground_truth["samples"].apply(lambda x: set(x.split(" ")))
    if "genes" in ground_truth.columns.values:
        ground_truth["genes"] = ground_truth["genes"].apply(lambda x: set(x.split(" ")))

    # prepare a dict with sample groups corresponding to known bicluster
    known_groups = {}
    for group in ground_truth.index.values:
        known_groups[group] = ground_truth.loc[group, "samples"]

    result = read_results(tool_name, result_file)
    # print(samples)

    best_matches = find_best_matches(result, known_groups, samples, FDR=0.05)
    print(best_matches)
    # try:
    return (best_matches, result)
    # except:
    #     return (0.0,result)
#
#
# name = 'dataframe'
# expr_file = '/home/andim/Downloads/A.n_genes=50,m=4,std=1,overlap=no.exprs_z.tsv'
# truth = '/home/andim/Downloads/A.n_genes=50,m=4,std=1,overlap=no.biclusters.tsv'
# result = '/home/andim/Downloads/A.n_genes=50,m=4,std=1,overlap=no.exprs_z-biclusters_df.tsv'
#
#
# run_eval(name, expr_file, truth, result)