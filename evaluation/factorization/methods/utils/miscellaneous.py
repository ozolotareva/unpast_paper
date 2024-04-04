from . import resultsHandler

def run_method(method, args):

    print(f'Running {args["output_path"]}...')
    # execute method
    result, runtime = method(**args)

    return result, runtime

def combination_to_string(d):
    del d['exprs_file']
    del d['output_path']
    if 'ground_truth_file' in d:
        del d['ground_truth_file']
                
    s = ''
    for key, value in d.items():
        s += f'{key}={value};'
    return s