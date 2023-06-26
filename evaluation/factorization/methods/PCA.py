from sklearn.decomposition import PCA
import sklearn
from .utils import interpret_results
import time


def run(exprs, n, random_state=101):
    start = time.time()
    with sklearn.config_context(assume_finite=True):
        transformer = PCA(n_components=n, random_state=random_state)
        exprs_transformed = transformer.fit_transform(exprs)
    result = interpret_results.format_sklearn_output(exprs_transformed, n, exprs.index)
    end = time.time()
    runtime = end-start
    return result, runtime