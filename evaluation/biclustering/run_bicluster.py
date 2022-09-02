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
    if tool_name in ['isa2', 'fabia', 'qubic']:
        command.append("Rscript")

    if tool_name in ['isa2', 'fabia', 'qubic', 'debi']:
        command.append(script_location)
        command.append(expr_file)
        command.append(out_file)
    if tool_name == 'debi':
        command.append("-b2")
    return command


subprocess.check_call(get_command(tool_name, script, expr_file, result_file))
j_weighted = eval_bicluster_methods.run_eval(tool_name=tool_name ,expr_file=expr_file, result_file=result_file, ground_truth_file=truth_file)
os.system(f'rm {result_file}')
print(f"J_weighted for {tool_name} and {expr_file}: {j_weighted}")
with open(score_file, 'w') as fw:
    fw.write(str(j_weighted))
