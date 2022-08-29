library(mosbi)

args = commandArgs(trailingOnly=TRUE)

in_file = args[1]
out_file = args[2]

data_matrix = as.matrix(read.table(in_file, sep="\t", header=TRUE, row.names=1))
BCqubic <- mosbi::run_fabia(data_matrix)
results = sapply(BCqubic, function(x){paste(paste(x@colname, collapse= " "),paste(x@rowname, collapse= " ") , sep= "\t")})
writeLines(results, con=out_file, sep="\n")
