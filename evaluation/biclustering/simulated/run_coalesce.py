import os
import subprocess
import sys
import uuid

in_file = sys.argv[1]
out_file = sys.argv[2]
param_file = sys.argv[3]

wd = f"/tmp/coalesce_{uuid.uuid4()}"
os.system(f"mkdir {wd}")
os.system(f"cp -r profiles/coalesce/* {wd}/")

params = dict()
with open(param_file, 'r') as fh:
    for line in fh:
        vals = line.strip().split("\t")
        params[vals[0]] = vals[1]

config_file = f"{wd}/algorithms/coalesce_config_1.conf"
out_lines = []
with open(config_file, 'r') as fr:
    for line in fr.readlines():
        sp = line.split("=")
        if sp[0] in params:
            out_lines.append(f"{sp[0]}={params[sp[0]]}")
        else:
            out_lines.append(line)

with open(config_file, 'w') as fw:
    fw.write("\n".join(out_lines))

os.system(f"cp {in_file} {wd}/dataset.tsv")
subprocess.check_call(["java", "-jar", "jbiclustge-cli.jar", "-run", wd])

result_file = ""
for dir in os.listdir(os.path.join(wd, "Biclustering_results1", "coalesce_config_1")):
    result_file = os.path.join(wd, "Biclustering_results1", "coalesce_config_1", dir, "JBiclustGE_csv.bicge")

os.system(f"cp {result_file} {out_file}")
os.system(f"rm -r {wd}")
