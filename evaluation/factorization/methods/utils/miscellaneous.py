from . import resultsHandler

def run_method(method, args):

    print(f'Running {args["output_path"]}...')
    # execute method
    result, runtime = method(**args)

    return result, runtime

