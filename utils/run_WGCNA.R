# usage: Rscript run_WGCNA.R binarized_expressions.tsv 
suppressPackageStartupMessages(library("WGCNA"))

args <- commandArgs(trailingOnly = TRUE)

fileBinExprs <- args[[1]]

fileModules <- paste0(sub(".tsv","",fileBinExprs),".modules.tsv")

datExpr <- read.table(fileBinExprs,check.names=FALSE)
datExpr[] <- lapply(datExpr, as.numeric)

#### finding power threshold #### 
powers = c(c(1:10), seq(from = 12, to=20, by=2))
rsqcuts <- seq(from = 0.9, to=0.05, by=-0.05)
i = 1
power <- NA
while (is.na(power)){
# Call the network topology analysis function
 sft <- pickSoftThreshold(datExpr, powerVector = powers,verbose = 0,
                          networkType = "signed hybrid",RsquaredCut =rsqcuts[[i]])
 power <- sft$powerEstimate
 i<-i+1
}
cat("\nR^2 threshold:",rsqcuts[[i]],"power:",power,"\n")

#### find modules ####
net = blockwiseModules(datExpr, power = power,
TOMType = "signed", networkType = "signed hybrid", minModuleSize = 2,
numericLabels = TRUE,
#detectCutHeight = 1-10^(-1*power),
#mergeCutHeight = 0.05,
verbose = 0)

moduleLabels = net$colors

sink(fileModules)
cat(paste0("module","\t","size","\t","genes","\n",sep=' '))
for (i in unique(moduleLabels)){
    cat(paste0(i,"\t",length(names(moduleLabels[moduleLabels==i])),"\t",paste(names(moduleLabels[moduleLabels==i]),collapse=" "),"\n"))
}
sink()

cat(fileModules,"\n")