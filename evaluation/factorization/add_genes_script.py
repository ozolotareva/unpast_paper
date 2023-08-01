# add genes to algoithms without genes

import os
import numpy as np
import pandas as pd
import sys,os
import random
import copy

import subprocess

import matplotlib.pyplot as plt
import seaborn as sns

from utils.eval import find_best_matches, generate_exprs

from methods import NMF, PCA, sparse_PCA, moCluster, MOFA2, iClusterPlus

gene_sets_are_defined = ['NMF', 'sparse_PCA']


classifications={"Intrinsic":["Luminal","Basal","Her2","Normal","Claudin-low"],
                "SCMOD2":["ER-/HER2-","ER+/HER2- Low Prolif","ER+/HER2- High Prolif","HER2+"],
                "IHC":["IHC_TNBC","IHC_ER","IHC_HER2","IHC_PR"]}

file_metabric_annotation = '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904.annotation_v6.tsv'
file_metabric_expression = '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904_17Kgenes.log2_exprs_v6.tsv'
file_metabric_subtypes = '/local/DESMOND2_data/v6/preprocessed_v6/METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv'
file_tcga_annotation = '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv'
file_tcga_expression = '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_v6.tsv'
file_tcga_subtypes = '/local/DESMOND2_data/v6/preprocessed_v6/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv'
file_gene_mapping = '/local/DESMOND2_data/v6/preprocessed_v6/gene_id_mapping.tsv'

out_dir = '/cosybio/project/hartung/unpast/unpast_real'

basename_t = "TCGA"
basename_m = "METABRIC" 
result_m = None
result_t = None

def add_genes(METHOD):
    global result_m, result_t
    
    method_name = METHOD.__name__.split('.')[-1]
    print('method_name:', method_name)

    #### Preparation
    # METABRIC
    file_path_m = file_metabric_expression
    output_path_m = os.path.join(out_dir, basename_m, method_name)
    ground_truth_file_m = file_metabric_annotation
    combinations_m = METHOD.generate_arg_list(file_path_m, output_path_m, ground_truth_file_m)
    # TCGA
    file_path_t = file_tcga_expression
    output_path_t = os.path.join(out_dir, basename_t, method_name)
    ground_truth_file_t = file_tcga_annotation
    combinations_t = METHOD.generate_arg_list(file_path_t, output_path_t, ground_truth_file_t)
    
    for _iteration, (comb_m, comb_t) in enumerate(zip(combinations_m, combinations_t)):
        
        if os.path.isfile(os.path.join(comb_m['output_path'], 'result.with_genes.tsv')) and os.path.isfile(os.path.join(comb_t['output_path'], 'result.with_genes.tsv')):
            print('Skipping because output files exist')
            continue
        
        result_m, _ = METHOD.run_real(comb_m, is_terminated=True)
        
        # save result as tsv for add genes input
        output_file = os.path.join(comb_m['output_path'], 'result.tsv')
        
        result_m['samples'] = result_m['samples'].map(list).map(lambda x: ' '.join(x))
        result_m.to_csv(output_file, sep='\t')
        
        print('output_file', output_file)
        
        # output_file = '/home/bba1401/Projects/unpast/DESMOND2/evaluation/factorization/METABRIC_1904_17Kgenes.Kmeans.clusters.tsv'

        process_m = subprocess.Popen([f"Rscript", "/home/bba1401/Projects/unpast/DESMOND2/utils/add_genes.R", output_file, file_metabric_expression, "FALSE"])
        
        
        
        result_t, _ = METHOD.run_real(comb_t, is_terminated=True)
        
        # save result as tsv for add genes input
        output_file = os.path.join(comb_t['output_path'], 'result.tsv')
        
        result_t['samples'] = result_t['samples'].map(list).map(lambda x: ' '.join(x))
        
        result_t.to_csv(output_file, sep='\t')
        print('output_file', output_file)

        process_t = subprocess.Popen([f"Rscript", "/home/bba1401/Projects/unpast/DESMOND2/utils/add_genes.R", output_file, file_tcga_expression, "TRUE"])
        
        process_m.wait()
        process_t.wait()
        
        
add_genes(MOFA2)