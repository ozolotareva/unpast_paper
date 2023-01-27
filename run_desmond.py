from .run_unpast import run


def run_DESMOND(*args, **kwargs):
    """ Wrapper of 'run_unpast' to support old function name """
    return run(*args, **kwargs)