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