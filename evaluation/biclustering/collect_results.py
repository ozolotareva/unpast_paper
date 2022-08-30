import sys, os

def collect(out_dir):

    for tool in os.listdir(out_dir):
        tool_path = os.path.join(out_dir, tool)
        for config in os.listdir(tool_path):
            config_path = os.path.join(tool_path,config)
            with open(os.path.join(tool_path,f'summary_{tool}_{config}.txt'), 'w') as fw:
                for result in os.listdir(config_path):
                    if 'summary' not in result:
                        with open(os.path.join(config_path,result), 'r') as fr:
                             for line in fr.readlines():
                                 fw.write(f'{result}\t{line}\n')
    print("Created summaries")
