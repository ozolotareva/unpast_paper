import os
from collections import defaultdict
import pandas as pd


def get_overall_best_performer(summary_df, tool_path):
    out_file = os.path.join(tool_path, 'mean_overall_performance.tsv')
    try:
        summary_df.groupby('parameters')["overall_performance"].mean().sort_values(ascending=False).to_csv(out_file,
                                                                                                           sep='\t')
    except:
        os.system(f'touch {out_file}')


def collect(out_dir):
    tool_name_dict = {'qubic': 'QUBIC', 'qubic2': 'QUBIC2', 'isa2': 'ISA2', 'debi': 'DeBi', 'coalesce': 'COALESCE',
                      'xmotifs': 'xMotifs', 'fabia': 'FABIA'}
    print("Creating summaries...")
    for tool in os.listdir(out_dir):
        tool_path = os.path.join(out_dir, tool)
        if os.path.isdir(tool_path):
            tool_name = tool_name_dict[tool.lower()]
            df = []
            for config in os.listdir(tool_path):
                if not os.path.isdir(os.path.join(tool_path, config)):
                    continue
                config_path = os.path.join(tool_path, config)
                with open(os.path.join(tool_path, f'summary_{tool}_{config}.txt'), 'w') as fw:
                    runs = defaultdict(dict)
                    for result in os.listdir(config_path):
                        if result.endswith(".score"):
                            if "-run" in result:
                                run = int(result.split("-run")[1].split(".")[0]) - 1
                                infile = result.split("-run")[0]
                            else:
                                run = 0
                                infile = result.replace(".score", "")
                            with open(os.path.join(config_path, result), 'r') as fr:
                                for line in fr.readlines():
                                    runs[infile][run] = float(line)
                        if result.endswith("-scores_df.tsv"):
                            if "-run" in result:
                                run = int(result.split("-run")[1].split("-")[0]) - 1
                            else:
                                run = 0
                            parts = result.split('.')
                            underscore_splits = result.split('_')
                            scenario_params = result[len(parts[0]) + 1:len(underscore_splits[0]) + 1 + len(
                                underscore_splits[1])]
                            scenario_params_sep = scenario_params.split(',')
                            try:
                                param_part = result[len(parts[0]) + len(scenario_params) + 2:len(
                                    result) - 14] if '_pre' not in result else result[
                                                                               len(result.split("_pre")[0]) + 1:len(
                                                                                   result) - 14]
                                d = {"scenario": parts[0], "run": run,
                                     'gsize': int(scenario_params_sep[0].split('=')[1]),
                                     'parameters': param_part.replace('-', ';').split(';run')[0], "method": tool_name, 'seed': None}
                            except:
                                d = {"scenario": parts[0], "run": run,
                                     'gsize': int(scenario_params_sep[0].split('=')[1]),
                                     'parameters': 'default', "method": tool_name, 'seed': None}
                            J_total = 0.0
                            performance = dict()
                            try:
                                best_matches = pd.read_csv(os.path.join(config_path, result), sep="\t", index_col=0,
                                                           header=0)
                                J_total = best_matches["J_weighted"].sum()
                                performance = best_matches.loc[:, ["J"]].to_dict()["J"]
                                # renaming subtype-specific performances to "performance_"+subt
                                subtypes = list(performance.keys())
                                for subt in subtypes:
                                    performance["performance_" + str(subt)] = performance.pop(subt)
                            except:
                                pass
                            performance["overall_performance"] = J_total
                            d.update(performance)
                            df.append(d)
                    for infile in runs.keys():
                        scores = list(runs[infile].values())
                        scores.sort()
                        score_string = ", ".join(map(str, scores))
                        fw.write(f'{infile}\t{(sum(scores) / len(scores))}\t{score_string}\n')
            summary_df = pd.DataFrame.from_records(df)
            summary_df.to_csv(os.path.join(out_dir, f'{tool_name}_ABC.tsv'), sep="\t", index=False)
            get_overall_best_performer(summary_df, tool_path)
    print("Created summaries")


if __name__ == '__main__':
    import sys

    collect(sys.argv[1])
