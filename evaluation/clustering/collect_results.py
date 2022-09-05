import os
from collections import defaultdict


def collect(out_dir):
    for tool in os.listdir(out_dir):
        tool_path = os.path.join(out_dir, tool)
        for config in os.listdir(tool_path):
            if '.txt' in config:
                continue
            config_path = os.path.join(tool_path, config)
            with open(os.path.join(tool_path, f'summary_{tool}_{config}.txt'), 'w') as fw:
                runs = defaultdict(dict)
                for result in os.listdir(config_path):
                    if 'summary' not in result:
                        run = result.split("-")[1].split(".")[0]
                        infile = result.split("-")[0]
                        with open(os.path.join(config_path, result), 'r') as fr:
                            for line in fr.readlines():
                                runs[infile][run] = float(line)
                for infile in runs.keys():
                    scores = list(runs[infile].values())
                    scores.sort()
                    score_string = ", ".join(map(str, scores))
                    fw.write(f'{infile}\t{(sum(scores)/len(scores))}\t{score_string}\n')
    print("Created summaries")