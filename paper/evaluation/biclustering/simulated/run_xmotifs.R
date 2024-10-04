library(mosbi)

args = commandArgs(trailingOnly=TRUE)

in_file = args[1]
out_file = args[2]
param_file = args[3]


ns=10
alpha=0.05

params = read.table(param_file,sep="\t")
for(idx in c(1:length(params[,1]))){
  if(params[idx,][[1]] == 'ns'){
    ns = as.numeric(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'alpha'){
    alpha = as.double(params[idx,][[2]])
  }
}


data_matrix = as.matrix(read.table(in_file, sep="\t", header=TRUE, row.names=1))
BCqubic <- mosbi::run_xmotifs(data_matrix)
results = sapply(BCqubic, function(x){paste(paste(x@colname, collapse= " "),paste(x@rowname, collapse= " ") , sep= "\t")})
writeLines(results, con=out_file, sep="\n")
