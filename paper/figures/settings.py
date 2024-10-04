METHOD_MAIN = "UnPaSt"

# , 'WGCNA_AffinityPropagation', 'WGCNA_AgglomerativeClustering', 'WGCNA_bikmeans', 'WGCNA_BIRCH', 'WGCNA_DBSCAN', 'WGCNA_GMM', 'WGCNA_Meanshift', 'WGCNA_HC', 'WGCNA_kmeans', 'WGCNA_Spectral', 'MeanShift', 'Spectral', 
METHODS_CLUSTER = ["kmeans", 'BisectingKMeans', 'MiniBatchKMeans', 'HierarchicalClustering', 'AffinityPropagation', 'AgglomerativeClustering',  'BIRCH', 'DBSCAN', 'GMM', 'SpectralClustering', 'Optics', 'MeanShift', 'mclust']
METHODS_FACTORIZATION = ["NMF","sparse_PCA","iClusterPlus","moCluster","MOFA2"]
METHODS_BICLUSTER = ["QUBIC","QUBIC2","ISA2","FABIA","COALESCE", "BiCoN", "DESMOND","BiMax", "GrandForest", 'Plaid'] #"DeBi",

METHODS_NETWORK_CONSTRAINT = ["GrandForest", "DESMOND", "BiCoN"]

METHODS = [method for method_list in [[METHOD_MAIN], METHODS_CLUSTER, METHODS_FACTORIZATION, METHODS_BICLUSTER] for method in method_list]

METHOD_MAIN_COLOR = "#ff0000"
METHODS_CLUSTER_COLOR = "#0173b2"
METHODS_FACTORIZATION_COLOR = "#de8f05"
METHODS_BICLUSTER_COLOR = "#029e73"

METHOD_TYPE_PALETTE = {
    METHOD_MAIN: METHOD_MAIN_COLOR,
    'Clustering': METHODS_CLUSTER_COLOR,
    'Factorization': METHODS_FACTORIZATION_COLOR,
    'Biclustering': METHODS_BICLUSTER_COLOR
}

def find_method_palette(method):
    if method in METHODS_CLUSTER:
        return METHODS_CLUSTER_COLOR
    elif method in METHODS_FACTORIZATION:
        return METHODS_FACTORIZATION_COLOR
    elif method in METHODS_BICLUSTER:
        return METHODS_BICLUSTER_COLOR
    elif method == METHOD_MAIN:
        return METHOD_MAIN_COLOR
    return 'black'

METHOD_PALETTE = {method: find_method_palette(method) for method in METHODS}
METHOD_PALETTE['AP'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['AC'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['HC'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['MB-kmeans'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['Bi-kmeans'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['Spectral'] = METHODS_CLUSTER_COLOR
METHOD_PALETTE['sPCA'] = METHODS_FACTORIZATION_COLOR
_to_add = {}
for key, value in METHOD_PALETTE.items():
    if '_' in key:
        _to_add[key.replace('_', ' ')] = value
METHOD_PALETTE.update(_to_add)

METHOD_PALETTE_DATASETS = {'METABRIC': METHOD_PALETTE, 'TCGA': METHOD_PALETTE}

TABLE_KEYS = ['overall_performance', 'performance_0.5', 'performance_0.25', 'performance_0.1', 'performance_0.05']

SCENARIOS = ['A', 'B', 'C']
GENE_SIZES = [5, 50, 500]

REAL_DATASETS = ['METABRIC', 'TCGA']

REAL_COLUMNS = ['Normal', 'Her2', 'Basal', 'LumB', 'LumA', 'Luminal',
       'Claudin-low', 'HER2+', 'ER+/HER2- High Prolif',
       'ER+/HER2- Low Prolif', 'ER-/HER2-', 'IHC TNBC', 'IHC PR',
       'IHC ER', 'IHC HER2', 'overall performance Intrinsic',
       'overall performance SCMOD2', 'overall performance IHC']

METHOD_ABBREVIATIONS = {
    'AgglomerativeClustering': 'AC',
    'AffinityPropagation': 'AP',
    'sparse_PCA': 'sPCA',
    'MiniBatchKMeans': 'MB-kmeans',
    'BisectingKMeans': 'Bi-kmeans',
    'HierarchicalClustering': 'HC',
    'SpectralClustering': 'Spectral'
}

CANCER_GROUPS = [
    ('PAM50 + Luminal + Claudin-low', ['Basal', 'Her2', 'Normal', 'Claudin-low', 'LumA', 'LumB', 'Luminal',]), #   
#     ('Luminal', ['Luminal']),
#     ('Claudin-low', ['Claudin-low']),
#  ('SCMOD2', ['ER+/HER2- Low Prolif', 'ER-/HER2-', 'HER2+', 'ER+/HER2- High Prolif']),
    ('IHC', ['IHC_ER', 'IHC_HER2', 'IHC_PR', 'IHC_TNBC'])]

DEFAULT_PARAMETERS = {
    'kmeans': {
        'k': 5,
        'init': 'k-means++',
        'max_iter': 300,
        'tol': 0.0001,
    },
    'mclust': {
        'k': 'default'
    },
    'HierarchicalClustering': {
        'linkage': 'single',
#         'distance': 'euclidean',
        'k': 5
    },
    'AffinityPropagation': {
        'damping': 0.5,
#         'distance': 'euclidean',
        'max_iter': 200,
#         'affinity': 'precomputed',
        # 'preference': 'default',
        # 'convergence_iter': 15
    },
#     'AgglomerativeClustering': {
#         'k': 2,
#         'affinity': 'euclidean',
# #         'memory': None,
# #         'connectivity': None,
# #         'compute_full_tree': 'auto',
#         'linkage': 'ward',
# #         'distance_threshold': None,
# #         'compute_distances': False
#     },
    'BisectingKMeans': {
        'k': 5, #8
#         'initmethod': 'random',
#         'n_init': 1,
        'tol': 0.0001,
        'max_iter': 300,
#         'algorithm': 'lloyd',
        'bisecting_strategy': 'biggest_inertia'
    },
    'BIRCH': {
        'threshold': 0.5,
        'branching_factor': 50,
        'k': 5,
    },
    'DBSCAN': {
        'eps': 0.5,
        'min_samples': 5,
#         'metric': 'euclidean',
#         'algorithm': 'auto'
    },
    'GMM': {
        'covtype': 'full',
        'k': 5,
        'initmethod': 'kmeans'
    },
    'SpectralClustering': {
        'k': 5,
#         'eigen_solver': 'arpack',
#         'affinity': 'rbf',
        'assign_labels': 'kmeans',
        'n_neighbors': 10,
        'n_inits': 10
    },
    'NMF': {'k': 5, #8 'None'
    'init': 'nndsvda',
    'tol': 0.0001,
    'alpha_W': -0.1,
    'alpha_H': 0.0,
    'shuffle': False,
    'solver': 'cd',
    'beta_loss': 'frobenius',
    'max_iter': 200},

    'iClusterPlus': {
    'n_cluster': 5, #2
    'type': 'gaussian',
    'burnin_n': 200,
    'draw_n': 200,
    'lambda_n': 'null',
    'iter_max': 20,
    'sdev': 0.05,
    'eps': 0.0001,
    'lambda_scale': 1},

    'moCluster': {'n_dimensions': 10,
    'n_cluster': 5, # 8
    'solver': 'fast',
    'center': True,
    'method': 'globalScore',
    'option': 'lambda1',
    'scale': False,
    'k': 0.1}, # 'all' is not working, the documentation also suggests 0.1

    'MOFA2': {'n_factors': 15,
    'n_cluster': 5, #5
    'ard_weights': True,
    'ard_factors': False,
    'likelihood': 'gaussian',
    'spikeslab_weights': True,
    'spikeslab_factors': False},

    'sparse_PCA': {'n_components': 5, # 20 'None'
    'alpha': 1,
    'ridge_alpha': 0.01,
    'max_iter': 1000,
    'method': 'cd',
    'tol': 1e-08},
    
    'QUBIC': {
        'r': 1,
        'q': 0.06,
        'c': 0.95,
        'f': 1,
        'type': 'default'
    },
    
    'ISA2': {
        'no_seeds': 100
    },
    
    'FABIA': {
        'alpha': 0.01,
        'spl': 0.0,
        'spz': 0.5,
        'p': 200
    },
    
    'COALESCE': {
        'prob_gene': 0.95,
        'pvalue_cond': 0.05,
        'pvalue_correl': 0.05,
        'zscore_cond': 0.05
    },
    
    'BiMax': {
        'minr': 2,
        'minc': 2,
        'bin_thr': 0,
        'number': 100
    },
    
    'Plaid': {
        'row_release': 0.7,
        'col_release': 0.7,
        'shuffle': 3,
        'iter_startup': 5,
        'iter_layer': 10,
        'max_layers': 20
    },
    'MiniBatchKMeans': {
        'k':5,
        'max_iter':300,
        'batch_size':1024,
        'max_no_improvement':10,
        'reassignment_ratio':0.01
    },
    'MeanShift': {
        'max_iters':300,
        'cluster_all':'True',
        'min_bin_freq':1
    }
}
