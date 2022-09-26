## ----message=FALSE------------------------------------------------------------
library(data.table)
library(MOFA2)
library(tidyr)

## -----------------------------------------------------------------------------
# check for arguments
args = commandArgs(trailingOnly=TRUE)
exprs_file <- args[1]
n_factors <- args[2]
n_cluster <- args[3]
seed <- args[4]
# /local/DESMOND2_data_simulated/simulated/C/C.n_genes=500,m=4,std=1,overlap=yes.exprs_z.tsv
data <- list(as.matrix(fread(exprs_file),rownames=1))
lapply(data,dim)

# create mofa object
MOFAobject <- create_mofa(data)

# create train options
data_opts <- get_default_data_options(MOFAobject)
# head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- as.numeric(n_factors)
# head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
train_opts$seed <- as.numeric(seed)
# head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(tempdir(),"mofa2_model.hdf5")

# train
start.time <- Sys.time()
MOFAobject.trained <- run_mofa(MOFAobject, outfile)
end.time <- Sys.time()
time.taken <- end.time - start.time

factors <- get_factors(MOFAobject.trained, factors = "all")


# cluster
# https://rdrr.io/github/bioFAM/MOFA2/man/cluster_samples.html
clusters <- cluster_samples(MOFAobject.trained, k=n_cluster, factors=1:n_factors)

write.table(clusters$cluster,file="MOFA2Cluster.csv", sep = ",")

fileConn<-file("MOFA2Runtime.txt")
writeLines(as.character(time.taken), fileConn)
close(fileConn)
