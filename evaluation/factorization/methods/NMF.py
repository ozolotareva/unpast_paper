from .utils import interpret_results
from sklearn.decomposition import NMF
import sklearn
import time


def preprocess_exprs(exprs):
    # remove negative values by shifiting by lowest value to positive
    # TODO this is very naive, think of a better way
    exprs_min = exprs.to_numpy().min()
    exprs = exprs - exprs_min
    return exprs

def run(exprs, k, random_state=101, transposed=False):
    start = time.time()
    with sklearn.config_context(assume_finite=True):
        model = NMF(n_components=k, init='random', random_state=random_state, max_iter=5000)
        W = model.fit_transform(exprs)
        H = model.components_
    if transposed:
        result = interpret_results.format_sklearn_output(W, k, exprs.index, transposed)
    else:
        result = interpret_results.format_sklearn_output(H, k, exprs.columns, transposed)
    end = time.time()
    runtime = end-start
    return result, runtime
