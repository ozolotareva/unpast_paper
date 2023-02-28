library(mclust)
library(glue)
library(stringr)

reformat_cluster_results <- function(mod_mclust){
  result <- data.frame(mod_mclust$classification)
  colnames(result) <- 'cluster'
  result$sample <- rownames(result)
  result_reformat <- data.frame()
  #colnames(result_reformat) <- c('samples', 'n_samples', 'frac')
  for (clusterno in unique(result$cluster)){
    samples <- result[result$cluster == clusterno, 'sample']
    samplesset <- paste0('{', paste(samples, collapse = ", "), '}')
    n_samples <- length(samples)
    frac <- n_samples /  dim(df)[1]
    result_reformat <- rbind(result_reformat, c(samplesset, n_samples, frac))
  }
  colnames(result_reformat) <- c('samples', 'n_samples', 'frac')
  return(result_reformat)
}


input_file = '/Users/fernando/Documents/Research/DESMOND2_data_simulated/simulated/A/A.n_genes=500,m=4,std=1,overlap=no.exprs_z.tsv'
result_file = '/Users/fernando/Documents/Research/DESMOND2/evaluation/clustering/results/DBSCAN/clusters/A.n_genes=500,m=4,std=1,overlap=no_run1.tsv'

seed_dict <- list(57451, 48699, 22057, 59467, 43106)
runno <-  as.integer(str_extract(result_file, "(?<=run)\\d+"))
seed <- seed_dict[[runno]]

#Reading input_file
df <- read.table(input_file)
df <- t(df)

set.seed(seed)
BIC <- mclustBIC(df, G=2:20)
#The mclustBIC function in the mclust R package is used for selecting the optimal number of clusters for model-based clustering of multivariate data. The function calculates the Bayesian Information Criterion (BIC) for a range of different Gaussian mixture models (GMMs) with varying numbers of clusters, and returns the GMM with the lowest BIC value, which is considered the optimal number of clusters for the given data.
# The BIC is a model selection criterion that balances the goodness of fit of the model to the data with the complexity of the model. It is defined as the negative log-likelihood of the data plus a penalty term that depends on the number of parameters in the model. The penalty term encourages simpler models with fewer parameters, thereby avoiding overfitting and improving the generalizability of the model.
#plot(BIC)
#summary(BIC)
set.seed(seed)
mod_mclust <- Mclust(df, x = BIC)
# summary(mod1, parameters = TRUE)
reformated_df <- reformat_cluster_results(mod_mclust = mod_mclust)
result_filename <- paste0(str_replace(result_file, '.tsv',''), glue('_seed_{seed}.tsv'))
write.table(reformated_df, result_filename, col.names = T, row.names = T, quote = F, sep = '\t')


# ICL <- mclustICL(df)
# LRT <- mclustBootstrapLRT(df, modelName = "VEI")

