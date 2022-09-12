library(mosbi)

args = commandArgs(trailingOnly = TRUE)

in_file = args[1]
out_file = args[2]
param_file = args[3]
param_file = '/home/andim/projects/ENCORE/DESMOND2/evaluation/biclustering/test.params'


no_seeds = 100

readParamsFile <- function(file) {
  out <- tryCatch({
    read.table(file, sep = "\t")
  }, error = function(cond) {
    return(as.table(matrix(c("", ""), ncol = 2)))
  },
    warning = function(cond) { return(as.table(matrix(c("", ""), ncol = 2))) }) }

params = readParamsFile(param_file)

for (idx in c(1:length(params[, 1]))) {
  if (params[idx,][[1]] == 'no_seeds') {
    no_seeds = as.numeric(params[idx,][[2]])
  }
}

data_matrix = as.matrix(read.table(in_file, sep = "\t", header = TRUE, row.names = 1))
BCqubic <- mosbi::run_isa(data_matrix, no.seeds = no_seeds)
results = sapply(BCqubic, function(x) { paste(paste(x@colname, collapse = " "), paste(x@rowname, collapse = " "), sep = "\t") })
writeLines(results, con = out_file, sep = "\n")
