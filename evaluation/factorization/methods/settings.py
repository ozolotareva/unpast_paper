import numpy as np

RANDOM_STATES = range(1, 6)
CLUSTER_RANGE = range(4, 21)
# CLUSTER_RANGE = [10]

# NMF
# NMF_ALPHAS = np.linspace(-0.4, 0.4, 9) # simulated
NMF_ALPHAS = np.linspace(-0.2, 0.2, 5) # np.linspace(-0.4, 0.4, 9)
NMF_BETA_LOSS = ['frobenius'] # , 'kullback-leibler' 
NMF_SOLVER = ['cd'] # , 'mu'
NMF_INIT = ['nndsvd', 'nndsvda', 'random']
NMF_SHUFFLE = [True, False]
NMF_TOL = [1e-4] 

# sparse PCA
S_PCA_ALPHA = [1, 5] # 1.1, 2.5, 0.5, 
S_PCA_RIDGE_ALPHA = [0.1, 0.01, 0.001]
S_PCA_MAX_ITER = [5000]
S_PCA_TOL = [1e-8] # 1e-9, 
S_PCA_METHOD = ['cd', 'lars']

# MOFA2
MOFA2_FACTORS = range(1, 21)

# MOCLUSTER
MOCLUSTER_RANGE = range(1, 21)
MOCLUSTER_RANDOM_STATES = [1]

# iCluster
ICLUSTER_TYPES = ["gaussian"]  # "binomial","poisson" throw errors, "multinomial" is too slow 
ICLUSTER_ALPHAS = [1]   # set to 1 as penalty for elasitcnet is not needed because elasticnet not used in this version of icluster
ICLUSTER_LAMBDA_NS = [5, 10, 25, 'null']   # null will be translated to "NULL" in R
ICLUSTER_LAMBDA_SCALES = [1]    # Value between 0 and 1
ICLUSTER_BURNIN_NS = [200]  # MCMC
ICLUSTER_DRAW_NS = [200]    # MCMC
ICLUSTER_MAXITERS = [20]
ICLUSTER_SDEVS = [0.05] # for random walk
ICLUSTER_EPS_LIST = [1e-4]
