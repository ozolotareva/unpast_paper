import subprocess
import sys, os
import eval_bicluster_methods

args = sys.argv

tool_name = args[1]
script = args[2]
expr_file = args[3]
truth_file = args[4]
result_file = args[5]
score_file = args[6]
config_file = args[7]


def read_config_file(file):
    params = dict()
    try:
        with open(file, 'r') as fh:
            for line in fh:
                vals = line.strip().split("\t")
                params[vals[0]] = vals[1]
    except:
        pass
    return params


def get_debi_command(script_location, expr_file, out_file, params):
    command = [script_location, expr_file, out_file]
    command.extend(["-b" + str(params.get('b', 2)), "-o" + str(params.get('o', 0.5)), "-p" + str(params.get('p', 'u')),
                    "-s" + str(params.get('s', 0))])
    return command


def get_qubic2_command(script_location, expr_file, disc_file, params):
    # disc_file = os.path.join('/tmp/qubic2_discretizations/', os.path.split(expr_file)[1])
    # os.system(f'cp -f {d_file} {disc_file}')

    command = [script_location, '-i', disc_file, '-o', '1000', '-d']
    if 'c' in params:
        command.append('-C')
    if 'n' in params:
        command.append('-N')
    return command


def get_command(tool_name, script_location, expr_file, out_file, param_file):
    command = []
    if tool_name in ['isa2', 'fabia', 'qubic']:
        command = ["Rscript", script_location, expr_file, out_file, param_file]
    else:
        params = read_config_file(param_file)
        if tool_name == 'debi':
            command = get_debi_command(script_location, expr_file, out_file, params)
        if tool_name == 'qubic2':
            command = get_qubic2_command(script_location, expr_file, out_file, params)
    print(command)
    return command


subprocess.check_call(get_command(tool_name, script, expr_file, result_file, config_file))
(scores, result) = eval_bicluster_methods.run_eval(tool_name=tool_name, expr_file=expr_file, result_file=result_file,
                                                   ground_truth_file=truth_file)

if ".chars" in result_file:
    result_file = result_file.replace('.chars', '')
try:
    result.to_csv(score_file.replace('.score', '-biclusters_df.tsv'), sep="\t")
except:
    os.system('touch ' + score_file.replace('.score', '-biclusters_df.tsv'))

j_weighted = 0.0
try:
    j_weighted = scores["J_weighted"].sum()
    scores.to_csv(score_file.replace('.score', '-scores_df.tsv'), sep="\t")
except:
    os.system(f'touch {score_file.replace(".score", "-scores_df.tsv")}')

print(f'move result file {result_file} to {score_file.replace(".score", "-raw.output")}')
os.system(f'mv {result_file} {score_file.replace(".score", "-raw.output")}')
print(f"J_weighted for {tool_name} and {expr_file}: {j_weighted}")
with open(score_file, 'w') as fw:
    fw.write(str(j_weighted))
