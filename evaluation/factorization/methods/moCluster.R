#!/usr/bin/Rscript

library(mogsa)
require(data.table)

# check for arguments
args = commandArgs(trailingOnly=TRUE)
exprs_file <- args[1]
dimensions <- args[2]
n_cluster <- args[3]
seed <- as.numeric(args[4])
output_folder <- args[5]
k <- as.numeric(args[6])
option <- args[7]
solver <- args[8]
center <- args[9] == 'True'
scale <- args[10] == 'True'
method <- args[11]


set.seed(seed)
exprs <- list(as.data.frame(as.matrix(fread(exprs_file),rownames=1)))

sapply(exprs, dim) # check dimensions of expression data

start.time <- Sys.time()
moas <- mbpca(exprs, ncomp = dimensions, k = k, method = method, option = option, 
              center=center, scale=scale, moa = TRUE, svd.solver = solver, maxiter = 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time

scrs <- moaScore(moas)

# samples <- moaCoef(moas)

# cluster the result scores
hcl <- hclust(dist(scrs))
cls <- cutree(hcl, k=n_cluster)

write.table(cls,file=file.path(output_folder, "mocluster_result.csv"), sep = ",")

fileConn<-file(file.path(output_folder, "mocluster_runtime.txt"))
writeLines(as.character(time.taken), fileConn)
close(fileConn)
