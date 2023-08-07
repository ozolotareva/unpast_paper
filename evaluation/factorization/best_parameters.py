# METABRIC

mocluster = {
    'exprs_file': '',
    'output_path': '/home/bba1401/llaima_data/data/unpast_real/asthma/moCluster/',
    'ground_truth_file': '',
    'n_dimensions': 1,
    'n_cluster': 8,
    'random_state': 1,
    'center': 'TRUE',
    'method': 'globalScore',
    'option': 'lambda1',
    'scale': 'FALSE',
    'k': 0.1,
    'solver': 'fast',
    'maxiter': 1000
    }

sparse_pca = {
    'n_components': 10,
    'alpha': 5,
    'ridge_alpha': 0.001,
    'max_iter': 5000,
    'method': 'cd',
    'tol': 1e-08  
}

MOFA2 = {
    'exprs_file': '',
    'output_path': '',
    'ground_truth_file': '',
    'n_factors': 19,
    'n_cluster': 1,
    'random_state': 0,
    'likelihoods': 'gaussian',
    'spikeslab_factors': False,
    'spikeslab_weights': True,
    'ard_factors': False,
    'ard_weights': True,
    'likelihood': 'gaussian'
}

iclusterplus = {
    'K': 15,
    'type': "gaussian",
    'n.lambda': 10,
    'scale.lambda': 1,
    'maxiter': 20,
    'eps': 1e-4,
    'n.burnin': 200,
    'n.draw': 200,
    'sdev': 0.05
}

class AvgRank:
    NMF = {'k': 2, #8
    'init': 'nndsvda',
    'tol': 0.0001,
    'alpha_W': -0.1,
    'alpha_H': 0.0,
    'shuffle': False,
    'solver': 'cd',
    'beta_loss': 'frobenius',
    'max_iter': 1000}

    iClusterPlus = {
    'n_cluster': 2, #5
    'type': 'gaussian',
    'burnin_n': 200,
    'draw_n': 200,
    'lambda_n': 'null',
    'iter_max': 20,
    'sdev': 0.05,
    'eps': 0.0001,
    'lambda_scale': 1}

    moCluster = {'n_dimensions': 2,
    'n_cluster': 2, # 8
    'random_state': 1,
    'solver': 'fast',
    'center': True,
    'method': 'globalScore',
    'option': 'inertia',
    'scale': False,
    'k': 0.1}

    MOFA2 = {'n_factors': 18,
    'n_cluster': 2, #5
    'random_state': 3,
    'ard_weights': True,
    'ard_factors': False,
    'likelihood': 'gaussian',
    'spikeslab_weights': True,
    'spikeslab_factors': False}

    sparse_PCA = {'n_components': 2, # 20
    'random_state': 1,
    'alpha': 5,
    'ridge_alpha': 0.001,
    'max_iter': 1000,
    'method': 'cd',
    'tol': 1e-08}
    
class Default:
    NMF = {'k': 3,
    'init': 'nndsvda',
    'tol': 0.0001,
    'alpha_W': 0.0,
    'alpha_H': 0.0,
    'shuffle': False,
    'solver': 'cd',
    'beta_loss': 'frobenius',
    'max_iter': 200}

    iClusterPlus = {
    'n_cluster': 2,
    'type': 'gaussian',
    'burnin_n': 200,
    'draw_n': 200,
    'lambda_n': 'null',
    'iter_max': 20,
    'sdev': 0.05,
    'eps': 0.0001,
    'lambda_scale': 1}

    moCluster = {'n_dimensions': 2,
    'n_cluster': 8,
    'random_state': 1,
    'solver': 'fast',
    'center': True,
    'method': 'globalScore',
    'option': 'inertia',
    'scale': False,
    'k': 0.1}

    MOFA2 = {'n_factors': 15,
    'n_cluster': 2,
    'ard_weights': True,
    'ard_factors': False,
    'likelihood': 'gaussian',
    'spikeslab_weights': True,
    'spikeslab_factors': False}

    sparse_PCA = {'n_components': 2,
    'alpha': 1,
    'ridge_alpha': 0.01,
    'max_iter': 1000,
    'method': 'cd',
    'tol': 1e-08}