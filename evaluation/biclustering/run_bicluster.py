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


def get_command(tool_name, script_location, expr_file, out_file):
    command = []
    if tool_name == 'qubic2':
        disc_file = os.path.join(os.path.split(out_file)[0], os.path.split(expr_file)[1])
        os.system(f'cp {expr_file} {disc_file}')
        expr_file = expr_file.replace(".chars", "")
        command.extend([script_location, '-i', disc_file])
    if tool_name in ['isa2', 'fabia', 'qubic']:
        command.append("Rscript")

    if tool_name in ['isa2', 'fabia', 'qubic', 'debi']:
        command.append(script_location)
        command.append(expr_file)
        if tool_name in ['isa2', 'fabia', 'qubic', 'debi']:
            command.append(out_file)
    if tool_name == 'debi':
        command.append("-b2")
    if tool_name == 'qubic2':
        command.append('-d')
    return command


subprocess.check_call(get_command(tool_name, script, expr_file, result_file))
j_weighted = eval_bicluster_methods.run_eval(tool_name=tool_name, expr_file=expr_file, result_file=result_file,
                                             ground_truth_file=truth_file)
print(f'move result file {result_file} to {os.path.join(os.path.split(score_file)[0], os.path.split(result_file)[1])}')
os.system(f'mv {result_file} {os.path.join(os.path.split(score_file)[0], os.path.split(result_file)[1])}')
print(f"J_weighted for {tool_name} and {expr_file}: {j_weighted}")
with open(score_file, 'w') as fw:
    fw.write(str(j_weighted))
