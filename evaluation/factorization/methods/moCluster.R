#!/usr/bin/Rscript

library(mogsa)
require(data.table)

# check for arguments
args = commandArgs(trailingOnly=TRUE)
exprs_file <- args[1]
dimensions <- args[2]
n_cluster <- args[3]

exprs <- list(as.data.frame(as.matrix(fread(exprs_file),rownames=1)))

sapply(exprs, dim) # check dimensions of expression data

start.time <- Sys.time()
moas <- mbpca(exprs, ncomp = dimensions, k = 0.1, method = "globalScore", option = "lambda1", 
              center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time

scrs <- moaScore(moas)

# samples <- moaCoef(moas)

# cluster the result scores
hcl <- hclust(dist(scrs))
cls <- cutree(hcl, k=n_cluster)

write.table(cls,file="moClusterCluster.csv", sep = ",")

fileConn<-file("moClusterRuntime.txt")
writeLines(as.character(time.taken), fileConn)
close(fileConn)