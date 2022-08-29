import os
import subprocess

test_case_folder = "/local/DESMOND2_data_simulated/simulated/"
script_folder = "./"

tool_list = {
    'fabia': 'run_fabia.R',
    'isa': 'run_isa.R',
    'qubic': 'run_qubic.R'
}

expr_files = {}
bicluster_files = {}

for mode in os.listdir(test_case_folder):
    mode_path = os.path.join(test_case_folder, mode)
    for case_file in os.listdir(mode_path):
        file_path = os.path.join(mode_path, case_file)
        prefix = case_file.split(".")[0]
        if "exprs.tsv" in case_file:
            expr_files[prefix] = file_path
        elif "biclusters.tsv" in case_file:
            bicluster_files[prefix] = file_path


def get_output_file(tool_name, case_prefix):
    return os.path.join("/tmp/", f'{case_prefix}_{tool_name}-default.tsv')




for test_case in expr_files.keys():
    expr_file = expr_files[test_case]
    true_file = bicluster_files[test_case]
    for tool_name in tool_list.keys():
        out_file = get_output_file(tool_name, test_case)
        subprocess.check_call(['python3','run_bicluster.py',tool_name, os.path.join(script_folder,tool_list[tool_name]), expr_file, true_file, out_file])


