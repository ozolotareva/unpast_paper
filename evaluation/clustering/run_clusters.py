import os
import subprocess
import sys
import time

# os.chdir('/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering')
# test_case_folder = "/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated"
# result_dir = "/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results"
test_case_folder = "/home/bbb1417/DESMOND2_benchmarking/DESMOND2_data/simulated/"
script_folder = "./"

tool_list = {
    'kmeans': 'run_kmeans.py',
    # 'WGCNAkmeans': 'run_WGCNAkmeans.py',
    'HC': 'run_HC.py',
    # 'WGCNAHC': 'run_WGCNAHC.py',
    'AffinityPropagation': 'run_AffinityPropagation.py',
    'Meanshift': 'run_MeanShift.py',
    'Spectral': 'run_Spectral.py',
    'AgglomerativeClustering': 'run_AgglomerativeClustering.py',
    'DBSCAN': 'run_DBSCAN.py',
    # 'OPTICS': 'run_OPTICS.py',
    'BIRCH': 'run_BIRCH.py',
    'bikmeans': 'run_BisectingKmeans.py',
    'GMM': 'run_GaussianMixmodels.py'
}

WGCNA_expr_files = {}
expr_files = {}
bicluster_files = {}

for mode in os.listdir(test_case_folder):
    print(mode)
    mode_path = os.path.join(test_case_folder, mode)
    print(mode_path)
    for case_file in os.listdir(mode_path):
        print(case_file)
        file_path = os.path.join(mode_path, case_file)
        print(file_path)
        prefix = case_file.split(".")[0] + "." + case_file.split(".")[1]
        print(prefix)
        if "exprs" in case_file:
            if "_MEs" in file_path:
                WGCNA_expr_files[prefix] = file_path
            else:
                expr_files[prefix] = file_path
        elif "biclusters" in case_file:
            bicluster_files[prefix] = file_path

result_dir = "/home/bbb1417/DESMOND2_benchmarking/new_results_DESMOND2"

if not os.path.exists(result_dir):
    os.system(f"mkdir {result_dir}")

commands = list()
running = list()

for tool_name in tool_list.keys():
    score_dir = os.path.join(result_dir, tool_name)
    print(score_dir)
    if not os.path.exists(score_dir):
        os.system("mkdir " + score_dir)

    clusters_dir = os.path.join(score_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        os.system("mkdir " + clusters_dir)

    scores_file = os.path.join(score_dir, f'{tool_name}_scores.txt')
    for test_case in expr_files.keys():
        # print(test_case)
        expr_file = expr_files[test_case]
        true_file = bicluster_files[test_case]
        for r in range(1, 6):
            commands.append(
                ['python3', 'run_cluster.py', tool_name, os.path.join(script_folder, tool_list[tool_name]), expr_file,
                 true_file, os.path.join(clusters_dir, f'{test_case}_run{r}.tsv'), scores_file])

for tool_name in tool_list.keys():
    score_dir = os.path.join(result_dir, 'WGCNA_' + tool_name)
    print(score_dir)
    if not os.path.exists(score_dir):
        os.system("mkdir " + score_dir)

    clusters_dir = os.path.join(score_dir, 'clusters')
    if not os.path.exists(clusters_dir):
        os.system("mkdir " + clusters_dir)

    scores_file = os.path.join(score_dir, f'{tool_name}_scores.txt')
    for test_case in WGCNA_expr_files.keys():
        # print(test_case)
        expr_file = WGCNA_expr_files[test_case]
        true_file = bicluster_files[test_case]
        for r in range(1, 6):
            commands.append(
                ['python3', 'run_cluster.py', tool_name, os.path.join(script_folder, tool_list[tool_name]), expr_file,
                 true_file, os.path.join(clusters_dir, f'{test_case}_run{r}.tsv'), scores_file])

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
