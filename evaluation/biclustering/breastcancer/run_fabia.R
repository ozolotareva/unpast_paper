library(mosbi)

args = commandArgs(trailingOnly=TRUE)

in_file = args[1]
out_file = args[2]
param_file = args[3]

alpha=0.01
cyc=500
spl=0.0
spz=0.5
center=2
lap=1.0

readParamsFile <- function(file) {
  out <- tryCatch({
    read.table(file, sep = "\t")
  }, error = function(cond) {
    return(as.table(matrix(c("", ""), ncol = 2)))
  },
    warning = function(cond) { return(as.table(matrix(c("", ""), ncol = 2))) }) }

params = readParamsFile(param_file)

for(idx in c(1:length(params[,1]))){
  if(params[idx,][[1]] == 'alpha'){
    alpha = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'cyc'){
    cyc = as.numeric(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'spl'){
    spl = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'spz'){
    spz = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'center'){
    center = as.numeric(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'lap'){
    lap=as.double(params[idx,][[2]])
  }
}



data_matrix = as.matrix(read.table(in_file, sep="\t", header=TRUE, row.names=1))
BCqubic <- mosbi::run_fabia(data_matrix, p=200, cyc=cyc, alpha=alpha, spl=spl, spz=spz, center=center)
results = sapply(BCqubic, function(x){paste(paste(x@colname, collapse= " "),paste(x@rowname, collapse= " ") , sep= "\t")})
writeLines(results, con=out_file, sep="\n")
