import copy
import os
import subprocess
import sys
import time
import collect_results

test_case_folder = "/local/DESMOND2_data_simulated/simulated/"
script_folder = "./"

tool_list = {
    'fabia': {'name': 'run_fabia.R', 'deterministic': False, 'precompute': False, 'params': {
        'alpha': [0.01, 0.02, 0.05], 'cyc': [250, 500, 1000], 'spl': [0.0, 0.5, 1.0, 2.0], 'spz': [0.5, 1.0, 2.0],
        'center': [2], 'lap': [0.75, 1.0, 2]
    }},
    'isa2': {'name': 'run_isa2.R', 'deterministic': False, 'precompute': False, 'params':{}},
    'qubic': {'name': 'run_qubic.R', 'deterministic': True, 'precompute': False, 'params': {
        'r': [1, 2, 5, 10, 25], 'q': [0.04, 0.06, 0.1],
        'c': [0.99, 0.95, 0.92], 'f': [0.5, 1, 5], 'P': False, 'C': False
    }},
    'debi': {'name': './debi', 'deterministic': True, 'precompute': False,
             'params': {'s': [0], 'o': [0, 0.5, 1], 'b': [-2, -1, 0, 1, 2], 'p': ['u']}},
    'qubic2': {'name': 'qubic2-master/qubic', 'deterministic': True, 'precompute': True,
               'params': {'c': True, 'n': True}}
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


def create_new_combos(combo, key, value):
    combos = []
    if type(value) is bool:
        combo1 = copy.copy(combo)
        combo1.update({key: value})
        combos.append(combo1)
        combo2 = copy.copy(combo)
        combo2.update({key: not value})
        combos.append(combo2)
    else:
        for v in value:
            new_combo = copy.copy(combo)
            new_combo.update({key: v})
            combos.append(new_combo)
    return combos


def create_param_combinations(tool_name):
    combos = [{}]
    params = tool_list[tool_name]['params']
    for k, v in params.items():
        new_combos = []
        for combo in combos:
            new_combos.extend(create_new_combos(combo, k, v))
        combos = new_combos
    return combos


def params_to_string(params):
    string = ""
    for k, v in params.items():
        if type(v) is bool:
            string += str(k) + "=" + ('T' if v else 'F')
        else:
            string += str(k) + "=" + str(v)
        string += "-"
    if len(string) > 0:
        return string[0:len(string) - 1]
    return string


def write_params(file, params):
    with open(file, 'w') as fw:
        for (k, v) in params.items():
            if type(v) is bool and v:
                fw.write(str(k) + "\t" + str(v) + "\n")
            else:
                fw.write(str(k) + "\t" + str(v) + "\n")


def get_output_file(tool_name, case_prefix):
    return os.path.join("/tmp/", f'{case_prefix}_{tool_name}-.tsv')


def run_tasklist(tasks, threads):
    while len(tasks) > 0 or len(running) > 0:
        if len(running) < threads and len(tasks) > 0:
            command = tasks[0]
            print(f'running command: {command}')
            tasks = tasks[1:]
            p = subprocess.Popen(command)
            running.append(p)
        done = []
        for i in range(0, len(running)):
            if running[i].poll() is not None:
                done.append(i)
        for i in done:
            del running[i]
        time.sleep(1)


result_dir = "/tmp/desmond2_bicluster_eval_results"
if os.path.exists(result_dir):
    os.system(f"rm -rf {result_dir}")
os.system(f"mkdir {result_dir}")

precomputing_multi = list()
commands_multi = list()
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
                disc_dir = "/tmp/qubic2_discretizations/"
                if not os.path.exists(disc_dir):
                    os.system(f"mkdir {disc_dir}")
                discretization_input = os.path.join(disc_dir, os.path.split(expr_file)[1])
                out_file = discretization_input + ".chars"
                os.system(f'cp -f {expr_file} {discretization_input}')
                if not os.path.exists(out_file):
                    precomputing_multi.append(
                        [os.path.join(script_folder, tool_list[tool_name]['name']), '-i', discretization_input,
                         '-F', '-R'])

                expr_file = discretization_input

        for params in create_param_combinations(tool_name):
            param_string = str(params_to_string(params))
            param_file = os.path.join(score_dir, os.path.split(expr_file)[1] + ".params")
            write_params(param_file, params)
            out_file_params = out_file.replace('.tsv', 'params=' + param_string + '.tsv')
            if tool_name == 'qubic2':
                precomputing_multi.append(['cp', '-f', expr_file + '.chars', out_file_params])
            if tool_list[tool_name]['deterministic']:
                command = ['python3', 'run_bicluster.py', tool_name,
                           os.path.join(script_folder, tool_list[tool_name]['name']),
                           expr_file, true_file, out_file_params,
                           os.path.join(score_dir, f'{test_case}_{param_string}.score'),
                           param_file]
                if tool_name in ['fabia']:
                    commands_multi.append(command)
                commands.append(command)
            else:
                for r in range(1, 6):
                    command = ['python3', 'run_bicluster.py', tool_name,
                               os.path.join(script_folder, tool_list[tool_name]['name']),
                               expr_file,
                               true_file, out_file_params,
                               os.path.join(score_dir, f'{test_case}_{param_string}-run{r}.score'), param_file]
                    if tool_name in ['fabia']:
                        commands_multi.append(command)
                    else:
                        commands.append(command)

allowed_threads = parallel_execs = int(sys.argv[1])
run_tasklist(precomputing_multi, 1)
run_tasklist(commands_multi, 1)
run_tasklist(commands, allowed_threads)

print(f"All done, result scores are in {result_dir}")
collect_results.collect(result_dir)
