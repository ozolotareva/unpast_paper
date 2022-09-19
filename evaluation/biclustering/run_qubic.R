library(mosbi)

args = commandArgs(trailingOnly=TRUE)

in_file = args[1]
out_file = args[2]
param_file = args[3]


r=1
q=0.06
c=0.95
f=1
P=FALSE
C=FALSE
type='default'

readParamsFile <- function(file) {
  out <- tryCatch({
    read.table(file, sep = "\t")
  }, error = function(cond) {
    return(as.table(matrix(c("", ""), ncol = 2)))
  },
    warning = function(cond) { return(as.table(matrix(c("", ""), ncol = 2))) }) }

params = readParamsFile(param_file)

for(idx in c(1:length(params[,1]))){
  if(params[idx,][[1]] == 'r'){
    r = as.numeric(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'q'){
    q = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'c'){
    c = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'f'){
    f = as.double(params[idx,][[2]])
  }
  if(params[idx,][[1]] == 'P'){
    P=TRUE
  }
  if(params[idx,][[1]] == 'C'){
    C=TRUE
  }
  if(params[idx,][[1]] == 'type'){
    type = params[idx,][[2]]
  }
}


data_matrix = as.matrix(read.table(in_file, sep="\t", header=TRUE, row.names=1))
BCqubic <- mosbi::run_qubic(data_matrix, o = 1000, r=r, q=q, c=c, f=f, P=P, C=C, type=type)
results = sapply(BCqubic, function(x){paste(paste(x@colname, collapse= " "),paste(x@rowname, collapse= " ") , sep= "\t")})
writeLines(results, con=out_file, sep="\n")
