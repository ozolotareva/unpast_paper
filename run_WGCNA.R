# usage: Rscript run_WGCNA.R binarized_expressions.tsv power1 power2 
suppressPackageStartupMessages(library("WGCNA"))

args <- commandArgs(trailingOnly = TRUE)

p1 <- as.numeric(args[[1]]) # power for TOM, e.g. 10
p2 <- as.numeric(args[[2]]) # power for detectCutHeight, e.g. 10
fileBinExprs <- args[[3]]# e.g. "/home/olya/TUM/DESMOND/DESMOND2/tmp_results/TCGA-BRCA.pv=0.01,method=GMM,direction=DOWN.bin_exprs.tsv"

fileModules <- paste0(sub(".tsv","",fileBinExprs),".modules.tsv")

datExpr <- read.table(fileBinExprs,check.names=FALSE)
datExpr[] <- lapply(datExpr, as.numeric)

net = blockwiseModules(datExpr, power = p1,
TOMType = "signed",networkType = "signed hybrid", minModuleSize = 2,
mergeCutHeight = 0.05,
numericLabels = TRUE,
detectCutHeight = 1-10^(-1*p2),
verbose = 0)

moduleLabels = net$colors

sink(fileModules)
cat(paste0("module","\t","size","\t","genes","\n",sep=' '))
for (i in unique(moduleLabels)){
    #if (i != 0){
    cat(paste0(i,"\t",length(names(moduleLabels[moduleLabels==i])),"\t",paste(names(moduleLabels[moduleLabels==i]),collapse=" "),"\n"))
    #}
    #else{
    #    not_clustered = names(moduleLabels[moduleLabels==i])
    #}
}
sink()

cat(fileModules,"\n")