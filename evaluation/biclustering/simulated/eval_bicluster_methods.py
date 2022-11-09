import os

import pandas as pd
from utils.eval import find_best_matches


def read_results(tool_name, result_file):
    if tool_name in ['isa2', 'fabia', 'qubic']:
        result = pd.read_csv(result_file, names=['samples', 'genes'], sep="\t")
        result["samples"] = result["samples"].apply(lambda x: set(x.strip().split(" ")))
        result["genes"] = result["genes"].apply(lambda x: set(x.strip().split(" ")))
        return result
    elif tool_name == 'coalesce':
        samples = []
        genes = []
        gene_lines = []
        sample_lines = []
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith("Gene_names"):
                    gene_lines.append(line.strip())
                elif line.startswith("Condition_names"):
                    sample_lines.append(line.strip())
        for entry in range(0, len(gene_lines)):
            genes.append(set(sorted(gene_lines[entry].split(":\t")[1].strip().split("\t"))))
            samples.append(set(sorted(sample_lines[entry].split(":\t")[1].strip().split("\t"))))
        return pd.DataFrame({'samples': samples, 'genes': genes})
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
        gene_lines = []
        sample_lines = []
        with open(result_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith(" Genes"):
                    gene_lines.append(line.strip())
                elif line.startswith(" Conds"):
                    sample_lines.append(line.strip())
        for entry in range(0, len(gene_lines)):
            genes.append(
                set() if '[0]' in gene_lines[entry] else set(gene_lines[entry].split(": ")[1].strip().split(" ")))
            samples.append(
                set() if '[0]' in sample_lines[entry] else set(sample_lines[entry].split(": ")[1].strip().split(" ")))

        return pd.DataFrame({'samples': samples, 'genes': genes})
    elif tool_name == 'dataframe':
        result = pd.read_csv(result_file, delimiter="\t", index_col=0)
        result["samples"] = result["samples"].apply(eval)
        result["genes"] = result["genes"].apply(eval)
        return result


def run_eval(tool_name, expr_file, ground_truth_file, result_file):
    # expression file

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
    # print(result)
    # print(result)
    # result.to_csv("/home/andim/Downloads/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv.qubic2_df.tsv", sep="\t")
    # print(samples)
    try:
        best_matches = find_best_matches(result, known_groups, samples, FDR=0.05)
    except:
        print(f"Error when scoring {expr_file} for {tool_name} with result file {result_file}")
        best_matches = None
    # try:
    return (best_matches, result)
    # except:
    #     return (0.0,result)

#
#
# name = 'dataframe'
# expr_file = '/home/andim/Downloads/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
# truth = '/home/andim/Downloads/A.n_genes=500,m=4,std=1,overlap=no.biclusters.tsv'
# result = '/home/andim/Downloads/A.n_genes=500,m=4,std=1,overlap=no_prob_gene=0.8-pvalue_cond=0.05-pvalue_motif=0.05-zscore_cond=0.05-zscore_motif=0.05-pvalue_correl:=0.05-size_minimum:=25-size_maximum=1000-fraction_postprocess=0.3-random=42-biclusters_df.tsv'
# # result2 = '/home/andim/projects/ENCORE/DESMOND2/evaluation/biclustering/jbiclustge/profiles/coalesce/Biclustering_results5/coalesce_config_1/coalesce_qy4fh8i677sk_2022-10-21_14-28-56/JBiclustGE_csv.bicge'
# #
# (scores, result) = run_eval(name, expr_file, truth, result)
# print(scores)
#
# # print(scores["J_weighted"].sum())
#
# (scores2, result2) = run_eval(name, expr_file, truth, result2)
# print(f"Equal: {scores2.equals(scores)}")
# # print(scores['genes'].iloc[3])
# # print(result)
# print(scores["J_weighted"].sum())
# print(scores2["J_weighted"].sum())