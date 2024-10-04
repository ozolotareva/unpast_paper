import copy
import os
import subprocess
import sys
import eval_bicluster_methods

case_folder = "/local/DESMOND2_data/v6/preprocessed_v6/"
script_folder = "./"
rerun_evaluations = True

tool_list = {
    # deterministic False
    'fabia': {'name': 'run_fabia.R', 'deterministic': True, 'precompute': False, 'params': {
        'alpha': [0.001, 0.01, 0.05], 'spl': [0.0, 0.5],
        'spz': [0.0, 0.5, 1.0], 'cyc': [500],
        'center': [2]
    }},
    # deterministic False
    'isa2': {'name': 'run_isa2.R', 'deterministic': True, 'precompute': False,
             'params': {
                 "no_seeds": [1, 2,
                              3,
                              4, 5] + list(range(10, 110, 20)) + [100, 125, 150, 200]
             }},
    'qubic': {'name': 'run_qubic.R', 'deterministic': True, 'precompute': False, 'params': {
        'r': [1, 2
            , 5, 10, 25
              ], 'q': [0.04, 0.06, 0.1, 0.25],
        'c': [0.99, 0.95, 0.92
            , 0.85, 0.75, 0.51
              ], 'f': [0.5, 1, 5], 'type': [
            'default',
            'area']
    }},
    'coalesce': {
        'name': 'run_coalesce.py', 'deterministic': True, 'precompute': False,
        'params': {'prob_gene': [0.8, 0.95],
                   "pvalue_cond": [0.1, 0.05],
                   "pvalue_motif": [0.1, 0.05],
                   "zscore_cond": [0.05, 0.1],
                   "zscore_motif": [0.05, 0.1],
                   "pvalue_correl:": [0.05, 0.1],
                   "size_minimum:": [5, 10, 25, 50],
                   "size_maximum": [10, 100, 1000],
                   "fraction_postprocess": [0.3, 0.5, 0.7],
                   "random": [42]
                   }
    },
    # 'xmotifs': {
    #     'name': 'run_xmotifs.R', 'deterministic': True, 'precompute': False, 'params': {
    #         'ns': [5, 10, 25, 50, 75, 100], 'alpha': [0.001, 0.01, 0.05, 0.1, 0.15]
    #     }
    # },
    # 'debi': {'name': './debi', 'deterministic': False, 'precompute': False,
    #          'params': {
    #              's': [0, 1, 3], 'o': [0, 0.5, 1], 'b': [0, 0.5, 1, 2], 'p': ['u']
    #          }},
    'qubic2': {'name': 'qubic2-master/qubic', 'deterministic': True, 'precompute': True, 'discretization_files': [],
               'params': {
                   # 'C': True, 'N': True,
                   'k': [2, 3,
                         10] + list(range(10, 50, 20)),
                   'c': [0.51, 0.75, 0.92]
               },
               'params_pre': {'n': True, 'R': True,
                              'q': [0.06
                                    # , 0.2
                                    # , 0.35
                                    ], 'r': [1, 2
                                             # , 5
                                             ]}
               }
}

base_name_t = "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2"
base_name_m = "METABRIC_1904_17Kgenes.log2"
expr_file_t = os.path.join(case_folder, "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv")
expr_file_m = os.path.join(case_folder, "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv")
subtype_file_t = os.path.join(case_folder, "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv")
subtype_file_m = os.path.join(case_folder, "METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv")
annot_file_t = os.path.join(case_folder, "TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv")
annot_file_m = os.path.join(case_folder, "METABRIC_1904.annotation_v6.tsv")


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


def create_pre_param_combinations(tool_name):
    combos = [{}]
    params = tool_list[tool_name]['params_pre']
    for k, v in params.items():
        new_combos = []
        for combo in combos:
            new_combos.extend(create_new_combos(combo, k, v))
        combos = new_combos
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


def run_tasklist(tasks, threads, silenced=False):
    total = len(tasks)
    running = []
    while len(tasks) > 0 or len(running) > 0:
        if len(running) < threads and len(tasks) > 0:
            command = tasks[0]
            tasks = tasks[1:]
            if command[0] == 'silence' or silenced:
                if command[0] == 'silence':
                    command = command[1:]
                print(f'running command: {command}')
                p = subprocess.Popen(command, stdout=subprocess.DEVNULL)
            else:
                print(f'running command: {command}')
                p = subprocess.Popen(command)
            running.append(p)
        done = []
        for i in range(0, len(running)):
            if running[i].poll() is not None:
                done.append(i)
                finished = 1 + total - (len(tasks) + len(running))
                print(f"Current command batch: {finished}/{total} DONE!")
        for i in done:
            del running[i]
        # time.sleep(1)


result_dir = "/tmp/desmond2_bicluster_eval_results_bc"
# if os.path.exists(result_dir):
#     os.system(f"rm -rf {result_dir}")
os.system(f"mkdir {result_dir}")

precomputing_multi = list()
commands_multi = list()
commands_multi_late = list()
commands_late = list()
commands = list()
running = list()
debi_commands = list()

cases = [[base_name_t, expr_file_t, subtype_file_t, annot_file_t],
         [base_name_m, expr_file_m, subtype_file_m, annot_file_m]]

for case in cases:
    base_name = case[0]
    expr_file = case[1]
    subtype_file = case[2]
    annot_file = case[3]
    for tool_name in tool_list.keys():
        score_dir = os.path.join(result_dir, tool_name)
        if not os.path.exists(score_dir):
            os.system("mkdir " + score_dir)
        score_dir = os.path.join(score_dir, 'default')
        if not os.path.exists(score_dir):
            os.system("mkdir " + score_dir)
        out_file = get_output_file(tool_name, base_name)

        if tool_list[tool_name]['precompute']:
            if tool_name == 'qubic2':
                disc_dir = "/tmp/q2_disc_bc/"
                if not os.path.exists(disc_dir):
                    os.system(f"mkdir {disc_dir}")
                for params in create_pre_param_combinations(tool_name):
                    name = "pre=" + str(params_to_string(params))
                    disc_param_dir = os.path.join(disc_dir, name)
                    if not os.path.exists(disc_param_dir):
                        os.system(f"mkdir {disc_param_dir}")
                    discretization_input = os.path.join(disc_param_dir, os.path.split(expr_file)[1])
                    if not os.path.exists(discretization_input):
                        print(f"creating {discretization_input}")
                        os.system(f'cp -f {expr_file} {discretization_input}')
                    tool_list[tool_name]['discretization_files'].append((discretization_input, name))
                    if not os.path.exists(discretization_input + ".chars"):
                        command = ['silence', os.path.join(script_folder, tool_list[tool_name]['name']), '-i',
                                   discretization_input, '-F']
                        if 'n' in params and 'R' in params and params['n'] and params['R']:
                            continue
                        if 'n' in params and params['n']:
                            command.append('-n')
                        elif 'R' in params and params['R']:
                            command.append('-R')
                        if 'q' in params:
                            command.extend(['-q', str(params['q'])])
                        if 'r' in params:
                            command.extend(['-r', str(params['r'])])
                        precomputing_multi.append(command)

        for params in create_param_combinations(tool_name):
            param_string = str(params_to_string(params))
            param_file = os.path.join(score_dir, os.path.split(expr_file)[1] + "-params=" + param_string + ".params")
            write_params(param_file, params)
            out_file_params = out_file.replace('.tsv', 'params=' + param_string + '.tsv')
            if tool_name == 'qubic2':
                for (expr_file, name) in tool_list[tool_name]['discretization_files']:
                    score_file = os.path.join(score_dir, f'{base_name}_{name}_{param_string}.score')
                    already_done = os.path.exists(score_file)
                    if not already_done or rerun_evaluations:
                        out_file_params = (expr_file + ".chars").replace('.tsv',
                                                                         '-' + name + '-params=' + param_string + '.tsv')
                        if not os.path.exists(score_file.replace('.score', '-biclusters_df.tsv')):
                            precomputing_multi.append(['cp', '-f', expr_file + '.chars', out_file_params])
                        commands.append(['python3', 'run_bicluster.py', tool_name,
                                         os.path.join(script_folder, tool_list[tool_name]['name']),
                                         expr_file, subtype_file, annot_file, out_file_params,
                                         score_file,
                                         param_file, str(rerun_evaluations)])
            elif tool_list[tool_name]['deterministic']:
                score_file = os.path.join(score_dir, f'{base_name}_{param_string}.score')
                already_done = os.path.exists(score_file.replace('.score', '-biclusters_df.tsv'))
                if not already_done or rerun_evaluations:
                    command = ['python3', 'run_bicluster.py', tool_name,
                               os.path.join(script_folder, tool_list[tool_name]['name']),
                               expr_file, subtype_file, annot_file, out_file_params,
                               score_file,
                               param_file, str(rerun_evaluations)]
                    if tool_name in ['isa2', 'fabia'] and not already_done:
                        commands_multi_late.append(command)
                    elif tool_name in ['debi'] and not already_done:
                        debi_commands.append(command)
                    else:
                        commands.append(command)
            else:
                for r in range(0, 5):
                    score_file = os.path.join(score_dir, f'{base_name}_{param_string}-run{r}.score')
                    already_done = os.path.exists(score_file.replace('.score', '-biclusters_df.tsv'))
                    if not already_done or rerun_evaluations:
                        command = ['python3', 'run_bicluster.py', tool_name,
                                   os.path.join(script_folder, tool_list[tool_name]['name']),
                                   expr_file,
                                   subtype_file, annot_file, out_file_params,
                                   score_file, param_file, str(rerun_evaluations)]

                        if tool_name in ['isa2', 'fabia'] and not already_done:
                            commands_multi_late.append(command)
                        elif tool_name in ['debi'] and not already_done:
                            debi_commands.append(command)
                        else:
                            commands.append(command)
    # break
allowed_threads = parallel_execs = int(sys.argv[1])
run_tasklist(precomputing_multi, 1)
run_tasklist(commands_multi, 1)
run_tasklist(commands, allowed_threads)
# eval_bicluster_methods.evaluate(wd=result_dir, methods=['coalesce'])
if len(commands_multi_late) > 0 and len(commands) > 0:
    eval_bicluster_methods.evaluate(result_dir)
run_tasklist(commands_late, allowed_threads)
if len(commands_late) > 0:
    eval_bicluster_methods.evaluate(result_dir)
run_tasklist(commands_multi_late, 1)
if len(commands_multi_late) > 0:
    eval_bicluster_methods.evaluate(result_dir)
# run_tasklist(debi_commands, allowed_threads)
eval_bicluster_methods.evaluate(result_dir)
print(f"All done, result scores are in {result_dir}")
