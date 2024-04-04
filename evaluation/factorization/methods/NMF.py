from .utils import interpret_results, resultsHandler
from sklearn.decomposition import NMF
import sklearn
import time
import os
import pandas as pd
from .settings import RANDOM_STATES, NMF_MAXITER, CLUSTER_RANGE, NMF_ALPHAS, NMF_BETA_LOSS, NMF_SOLVER, NMF_INIT, NMF_SHUFFLE, NMF_TOL
from .utils.miscellaneous import run_method


def generate_arg_list(exprs_file, output_folder, ground_truth_file, cluster_range=CLUSTER_RANGE):
    nmf_arguments = []
    # NMF - not transposed
    cluster_range = list(cluster_range)
    cluster_range.append(None)
    for m in RANDOM_STATES:
        for k in cluster_range:
            for tol in NMF_TOL:
                for alpha_W in NMF_ALPHAS:
                    alpha_W = round(alpha_W, 2)
                    for alpha_H in NMF_ALPHAS:
                        alpha_H = round(alpha_H, 2)
                        for solver in NMF_SOLVER:
                            for shuffle in NMF_SHUFFLE:
                                for beta_loss in NMF_BETA_LOSS:
                                    for init in NMF_INIT:
                                        for max_iter in NMF_MAXITER:
                                            if solver == 'cd' and beta_loss == 'kullback-leibler':
                                                continue
                                            if solver == 'mu' and init == 'nndsvd':
                                                #The multiplicative update ('mu') solver cannot update zeros present in 
                                                #the initialization, and so leads to poorer results when used jointly with
                                                #init='nndsvd'. You may try init='nndsvda' or init='nndsvdar' instead.
                                                continue
                                            output_path = os.path.join(output_folder, 
                                                f'k={k}',
                                                f'init={init}',
                                                f'tol={tol}',
                                                f'random_state={m}',
                                                f'alpha_W={alpha_W}',
                                                f'alpha_H={alpha_H}',
                                                f'shuffle={shuffle}',
                                                f'solver={solver}',
                                                f'beta_loss={beta_loss}',
                                                f'max_iter={max_iter}'
                                                )

                                            args = {'exprs_file': exprs_file,
                                                    'output_path': output_path,
                                                    'ground_truth_file': ground_truth_file,
                                                    'k': k,
                                                    'init': init,
                                                    'tol': tol,
                                                    'random_state': m,
                                                    'transposed': False,
                                                    'alpha_W': alpha_W,
                                                    'alpha_H': alpha_H,
                                                    'shuffle': shuffle,
                                                    'solver': solver,
                                                    'beta_loss': beta_loss,
                                                    'max_iter': max_iter
                                                }
                                            nmf_arguments.append(args)
    return nmf_arguments

def preprocess_exprs(exprs):
    # remove negative values by shifiting by lowest value to positive
    exprs_min = exprs.to_numpy().min()
    exprs = exprs - exprs_min
    return exprs

def execute_algorithm(exprs, k, output_path, init='random', tol=1e-4, random_state=101, transposed=False, alpha_W=0, max_iter=1, alpha_H=0, shuffle=False, l1_ratio=0, solver='cd', beta_loss='frobenius', **_):
    start = time.time()
    with sklearn.config_context(assume_finite=True):
        model = NMF(n_components=k, init=init, tol=tol, random_state=random_state, max_iter=max_iter, alpha_W=alpha_W, alpha_H=alpha_H, shuffle=shuffle, l1_ratio=l1_ratio, solver=solver, beta_loss=beta_loss)
        W = model.fit_transform(exprs)
        H = model.components_
        pd.DataFrame(W).to_csv(os.path.join(output_path, 'W.csv'))
        pd.DataFrame(H).to_csv(os.path.join(output_path, 'H.csv'))
        
    # W we cluster based on columns, H based on rows!!
    if transposed:
        result = interpret_results.format_sklearn_output(W, k, exprs.index, transposed)
        result_genes = interpret_results.format_sklearn_output(H, k, exprs.columns, not transposed)
    else:
        result = interpret_results.format_sklearn_output(H, k, exprs.columns, transposed)
        result_genes = interpret_results.format_sklearn_output(W, k, exprs.index, not transposed)
    end = time.time()
    runtime = end-start
    return (result, result_genes), runtime

def run_simulated(args):
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        print('Skipping because result exists:', args["output_path"])
        return resultsHandler.read_result(args["output_path"])
    df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0)
    args['exprs'] = preprocess_exprs(df_exprs)
    result, runtime = run_method(execute_algorithm, args)
    result, result_genes = result
    resultsHandler.write_samples(args["output_path"], args['exprs'].columns)
    del args['exprs']
    result_genes.to_csv(os.path.join(args["output_path"], 'result_genes.csv'))
    # save results
    resultsHandler.save(result, runtime, args["output_path"])
    return resultsHandler.read_result(args["output_path"])

def run_real(args, is_terminated=False ):
    print('The scikit-learn version is {}.'.format(sklearn.__version__))
    if is_terminated:
        try:
            return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
        except:
            return False, False
    print('here1')
    if resultsHandler.create_or_get_result_folder(args["output_path"]):
        # print('Returning existing results:', args["output_path"])
        print('here2')
        pass
    else:
        print('here3')
        df_exprs = pd.read_csv(args['exprs_file'], sep='\t', index_col=0)
        args['exprs'] = preprocess_exprs(df_exprs)
        result, runtime = run_method(execute_algorithm, args)
        result, result_genes = result
        print('result', result)
        resultsHandler.save(result, runtime, args["output_path"])
        del args['exprs']
        result_genes.to_csv(os.path.join(args["output_path"], 'result_genes.csv'))
    return resultsHandler.read_result(args["output_path"]), resultsHandler.read_runtime(args["output_path"])
