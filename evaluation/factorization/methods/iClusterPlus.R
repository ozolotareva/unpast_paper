library(iClusterPlus)

require(data.table)

# check for arguments
args = commandArgs(trailingOnly=TRUE)
exprs_file <- args[1]
lambda_n <- if (args[2] != 'null') as.numeric(args[2]) else NULL
n_cluster <- as.numeric(args[3])
lambda_scale <- as.numeric(args[4])
iter_max <- as.numeric(args[5])
eps <- as.numeric(args[6])  # 1.0e-4
type <- args[7]  # "gaussian","binomial","poisson","multinomial"
burnin_n <- as.numeric(args[8])
draw_n <- as.numeric(args[9])
sdev <- as.numeric(args[10])
output_directory <- args[11]
seed <- as.numeric(args[12])

set.seed(seed)
#load Data:
exprs <- as.matrix(fread(exprs_file),rownames=1)

start.time <- Sys.time()

# lambda_n=NULL

#tune directly with one set number of cluster and one lambda
result = tune.iClusterPlus(cpus = 1,
      dt1 = exprs, type = c(type), K = n_cluster,
      n.lambda = lambda_n, scale.lambda = c(lambda_scale),
      maxiter = iter_max, eps = eps, n.burnin=burnin_n, n.draw=draw_n, sdev=sdev);

end.time <- Sys.time()
time.taken <- end.time - start.time

cluster_assignment <- result$fit[[1]]$clusters

lambda_n = if (is.null(lambda_n)) 'null'  else lambda_n
n_cluster = n_cluster + 1
file <- sprintf("%s/iclusterplus_result.csv", output_directory)

write.table(cluster_assignment,file = file, sep = ",")

# save runtime
runtime_file <- sprintf("%s/iclusterplus_runtime.txt", output_directory)
fileConn<-file(runtime_file)
writeLines(as.character(time.taken), fileConn)
close(fileConn)



