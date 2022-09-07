import os
import subprocess
import sys
import time
import collect_results

test_case_folder = "/local/DESMOND2_data_simulated/simulated/"
script_folder = "./"

tool_list = {
    'fabia': {'name': 'run_fabia.R', 'deterministic': False, 'precompute': False},
    'isa2': {'name': 'run_isa2.R', 'deterministic': False, 'precompute': False},
    'qubic': {'name': 'run_qubic.R', 'deterministic': True, 'precompute': False},
    'debi': {'name': './debi', 'deterministic': True, 'precompute': False},
    'qubic2': {'name': 'qubic2-master/qubic', 'deterministic': False, 'precompute': True}
}

expr_files = {}
bicluster_files = {}

for mode in os.listdir(test_case_folder):
    mode_path = os.path.join(test_case_folder, mode)
    for case_file in os.listdir(mode_path):
        file_path = os.path.join(mode_path, case_file)
        prefix = case_file.split(".")[0] + "." + case_file.split(".")[1]
        if "exprs" in case_file:
            expr_files[prefix] = file_path
        elif "biclusters" in case_file:
            bicluster_files[prefix] = file_path
    break


def get_output_file(tool_name, case_prefix):
    return os.path.join("/tmp/", f'{case_prefix}_{tool_name}-default.tsv')


result_dir = "/tmp/desmond2_bicluster_eval_results"
if os.path.exists(result_dir):
    os.system(f"rm -rf {result_dir}")
os.system(f"mkdir {result_dir}")

commands = list()
running = list()

for test_case in expr_files.keys():
    for tool_name in tool_list.keys():
        expr_file = expr_files[test_case]
        true_file = bicluster_files[test_case]
        score_dir = os.path.join(result_dir, tool_name)
        if not os.path.exists(score_dir):
            os.system("mkdir " + score_dir)
        score_dir = os.path.join(score_dir, 'default')
        if not os.path.exists(score_dir):
            os.system("mkdir " + score_dir)
        out_file = get_output_file(tool_name, test_case)

        if tool_list[tool_name]['precompute']:
            if tool_name == 'qubic2':
                disc_dir = score_dir + "_disc"
                if not os.path.exists(disc_dir):
                    os.system(f"mkdir {disc_dir}")
                discretization_input = os.path.join(disc_dir, os.path.split(expr_file)[1])
                if not os.path.exists(discretization_input + ".chars"):
                    os.system(f'cp {expr_file} {discretization_input}')

                    expr_file = discretization_input + ".chars"
                    out_file = os.path.join(score_dir, os.path.split(expr_file)[1] + ".blocks")
                    commands.append(
                            [os.path.join(script_folder, tool_list[tool_name]['name']), '-i', discretization_input,
                             '-F',
                             '-R'])
                    if not tool_list[tool_name]['deterministic']:
                        for r in range(1, 6):
                            commands.append(
                                ['cp', expr_file, expr_file.replace('.tsv.chars', '-run' + str(r)) + ".tsv.chars"])

        if tool_list[tool_name]['deterministic']:
            commands.append(
                ['python3', 'run_bicluster.py', tool_name, os.path.join(script_folder, tool_list[tool_name]['name']),
                 expr_file,
                 true_file, out_file, os.path.join(score_dir, f'{test_case}_default.score')])
        else:
            for r in range(1, 6):
                if tool_name == 'qubic2':
                    commands.append(
                    ['python3', 'run_bicluster.py', tool_name,
                     os.path.join(script_folder, tool_list[tool_name]['name']),
                     expr_file.replace('.tsv.chars', '-run' + str(r)) + ".tsv.chars",
                     true_file, out_file.replace('.tsv.chars.blocks', '-run' + str(r)) + ".tsv.chars.blocks", os.path.join(score_dir, f'{test_case}_default-run{r}.score')])
                else:
                    commands.append(
                        ['python3', 'run_bicluster.py', tool_name,
                         os.path.join(script_folder, tool_list[tool_name]['name']),
                         expr_file,
                         true_file, out_file,
                         os.path.join(score_dir, f'{test_case}_default-run{r}.score')])
parallel_execs = int(sys.argv[1])
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

print(f"All done, result scores are in {result_dir}")
collect_results.collect(result_dir)
