# METABRIC

class BestBrcaMETABRIC: # PAM50
    sparse_PCA = {'n_components': 8, 'alpha': 1, 'ridge_alpha': 0.1, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
    NMF = {'k': 3, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 200}
    moCluster = {'n_dimensions': 15, 'n_cluster': 20, 'solver': 'fast', 'center': True, 'method': 'globalScore', 'option': 'lambda1', 'scale': False, 'k': 0.1}
    MOFA2 = {'n_factors': 19, 'n_cluster': 11, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    iClusterPlus = {'lambda_n': 10, 'n_cluster': 13, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    
    
class BestBrcaTCGA: # PAM50
    sparse_PCA = {'n_components': 8, 'alpha': 5, 'ridge_alpha': 0.1, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
    NMF = {'k': 8, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 1000}
    moCluster = {'n_dimensions': 4, 'n_cluster': 10, 'solver': 'fast', 'center': True, 'method': 'globalScore', 'option': 'inertia', 'scale': False, 'k': 1}
    MOFA2 = {'n_factors': 2, 'n_cluster': 7, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    iClusterPlus = {'lambda_n': 5, 'n_cluster': 11, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    
    
class OptimizedBRCAForAsthma: # average ranked, on BRCA (TCGA & METABRIC)
    NMF = {'k': 2, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 1000}
    moCluster = {'n_dimensions': 5, 'n_cluster': 2, 'solver': 'svd', 'center': True, 'method': 'globalScore', 'option': 'uniform', 'scale': False, 'k': 1}
    iClusterPlus = {'lambda_n': 10, 'n_cluster': 2, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    MOFA2 = {'n_factors': 12, 'n_cluster': 2, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    sparse_PCA = {'n_components': 2, 'alpha': 1, 'ridge_alpha': 0.001, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
 
class OptimizedBRCAForAsthmaARI: 
    # average ranked, on BRCA (TCGA & METABRIC)
    NMF = {'k': 2, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 200}
    moCluster = {'n_dimensions': 5, 'n_cluster': 2, 'solver': 'svd', 'center': True, 'method': 'globalScore', 'option': 'lambda1', 'scale': False, 'k': 1}
    iClusterPlus = {'lambda_n': 25, 'n_cluster': 2, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    MOFA2 = {'n_factors': 2, 'n_cluster': 2, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    sparse_PCA = {'n_components': 2, 'alpha': 1, 'ridge_alpha': 0.1, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
 
class OptimizedBRCA: # average ranked, on BRCA (TCGA & METABRIC)
    NMF = {'k': 5, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 1000}
    moCluster = {'n_dimensions': 5, 'n_cluster': 5, 'solver': 'svd', 'center': True, 'method': 'globalScore', 'option': 'uniform', 'scale': False, 'k': 1}
    iClusterPlus = {'lambda_n': 10, 'n_cluster': 5, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    MOFA2 = {'n_factors': 12, 'n_cluster': 5, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    sparse_PCA = {'n_components': 5, 'alpha': 1, 'ridge_alpha': 0.001, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
   
   
class OptimizedBRCACluster5: # average ranked, on BRCA (TCGA & METABRIC)
    NMF = {'k': 5, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 1000}
    moCluster = {'n_dimensions': 5, 'n_cluster': 5, 'solver': 'svd', 'center': True, 'method': 'globalScore', 'option': 'uniform', 'scale': False, 'k': 1}
    iClusterPlus = {'lambda_n': 10, 'n_cluster': 5, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    MOFA2 = {'n_factors': 12, 'n_cluster': 5, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    sparse_PCA = {'n_components': 5, 'alpha': 1, 'ridge_alpha': 0.001, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
   
   
class OptimizedBRCA: # average ranked, on BRCA (TCGA & METABRIC)
    NMF = {'k': 8, 'init': 'nndsvda', 'tol': 0.0001, 'transposed': False, 'alpha_W': -0.1, 'alpha_H': 0.0, 'shuffle': False, 'solver': 'cd', 'beta_loss': 'frobenius', 'max_iter': 1000}
    moCluster = {'n_dimensions': 5, 'n_cluster': 13, 'solver': 'svd', 'center': True, 'method': 'globalScore', 'option': 'uniform', 'scale': False, 'k': 1}
    iClusterPlus = {'lambda_n': 10, 'n_cluster': 12, 'lambda_scale': 1, 'iter_max': 20, 'eps': 0.0001, 'type': 'gaussian', 'burnin_n': 200, 'draw_n': 200, 'sdev': 0.05}
    MOFA2 = {'n_factors': 12, 'n_cluster': 9, 'ard_weights': True, 'ard_factors': False, 'likelihood': 'gaussian', 'spikeslab_weights': True, 'spikeslab_factors': False}
    sparse_PCA = {'n_components': 9, 'alpha': 1, 'ridge_alpha': 0.001, 'max_iter': 1000, 'method': 'cd', 'tol': 1e-08}
    
class DefaultAsthma:
    NMF = {'k': 2,
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

    moCluster = {'n_dimensions': 10,
    'n_cluster': 2,
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
    
    
class DefaultBRCA:
    NMF = {'k': 5,
    'init': 'nndsvda',
    'tol': 0.0001,
    'alpha_W': 0.0,
    'alpha_H': 0.0,
    'shuffle': False,
    'solver': 'cd',
    'beta_loss': 'frobenius',
    'max_iter': 200}

    iClusterPlus = {
    'n_cluster': 5,
    'type': 'gaussian',
    'burnin_n': 200,
    'draw_n': 200,
    'lambda_n': 'null',
    'iter_max': 20,
    'sdev': 0.05,
    'eps': 0.0001,
    'lambda_scale': 1}

    moCluster = {'n_dimensions': 10,
    'n_cluster': 5,
    'random_state': 1,
    'solver': 'fast',
    'center': True,
    'method': 'globalScore',
    'option': 'inertia',
    'scale': False,
    'k': 0.1}

    MOFA2 = {'n_factors': 15,
    'n_cluster': 5,
    'ard_weights': True,
    'ard_factors': False,
    'likelihood': 'gaussian',
    'spikeslab_weights': True,
    'spikeslab_factors': False}

    sparse_PCA = {'n_components': 5,
    'alpha': 1,
    'ridge_alpha': 0.01,
    'max_iter': 1000,
    'method': 'cd',
    'tol': 1e-08}
    
        