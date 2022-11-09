import os
import subprocess
import sys
import time

'''
real_data_path = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/preprocessed_v6/'
result_dir = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/resultsRD/'
script_folder = "/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/"

'''
real_data_path = '/home/bbb1417/DESMOND2_benchmarking/DESMOND2_data/preprocessed_v6/'
result_dir = '/home/bbb1417/DESMOND2_benchmarking/new_clusteringRD'
script_folder = "./"



if not os.path.exists(result_dir):
    # os.system(f"rm -rf {result_dir}")
    os.system(f"mkdir {result_dir}")

file_metabric_expression = f'{real_data_path}METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv'
file_metabric_annotation = f'{real_data_path}METABRIC_1904.annotation_v6.tsv'
file_metabric_subtypes = f'{real_data_path}METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv'

file_tcga_expression = f'{real_data_path}TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv'
file_tcga_annotation = f'{real_data_path}TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv'
file_tcga_subtypes = f'{real_data_path}TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv'

file_gene_mapping = f'{real_data_path}gene_id_mapping.tsv'

print(file_metabric_expression, file_metabric_annotation, file_metabric_subtypes)

print(file_tcga_expression, file_tcga_annotation, file_tcga_subtypes)

print(file_gene_mapping)

basename_t = "TCGA"
basename_m = "METABRIC"

tool_list = {
    #'kmeans': 'run_kmeans.py',
    #'WGCNAkmeans': 'run_WGCNAkmeans.py',
    #'HC': 'run_HC.py',
    #'WGCNAHC': 'run_WGCNAHC.py',
    #'AffinityPropagation': 'run_AffinityPropagation.py',
    #'Meanshift': 'run_MeanShift.py',
    #'Spectral': 'run_Spectral.py',
    #'AgglomerativeClustering': 'run_AgglomerativeClustering.py',
    #'DBSCAN': 'run_DBSCAN.py',
    #'OPTICS': 'run_OPTICS.py',
    #'BIRCH': 'run_BIRCH.py',
    'bikmeans': 'run_BisectingKmeans.py',
    'GMM': 'run_GaussianMixmodels.py'
}

commands = list()
running = list()

for tool_name in tool_list.keys():
    score_dir = os.path.join(result_dir, tool_name)

    if not os.path.exists(score_dir):
        os.system("mkdir " + score_dir)

    clusters_dir = os.path.join(score_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        os.system("mkdir " + clusters_dir)

    scores_file = os.path.join(score_dir,  f'{tool_name}_scoresRD.txt')

    for r in range(1, 6):
        commands.append(['python3', 'run_cluster_rd.py', tool_name, os.path.join(script_folder, tool_list[tool_name]), file_metabric_expression, file_metabric_annotation, file_metabric_subtypes, os.path.join(clusters_dir, f'METABRIC_run{r}.tsv'), scores_file])
        commands.append(['python3', 'run_cluster_rd.py', tool_name, os.path.join(script_folder, tool_list[tool_name]), file_tcga_expression, file_tcga_annotation, file_tcga_subtypes, os.path.join(clusters_dir, f'TCGABRCA_run{r}.tsv'), scores_file])


print(f"Commands running")
parallel_execs = int(sys.argv[1])
# parallel_execs = 1
while len(commands) > 0 or len(running) > 0:
    if len(running) < parallel_execs and len(commands) > 0:
        command = commands[0]
        print(f'starting command: {command}')
        commands = commands[1:]
        p = subprocess.Popen(command)
        running.append(p)
    done = []
    for i in range(0, len(running)):
        if running[i].poll() is not None:
            done.append(i)
    for i in done:
        del running[i]
    time.sleep(1)

print(f"All done, result scores are in their directories")
# collect_results.collect(result_dir)
