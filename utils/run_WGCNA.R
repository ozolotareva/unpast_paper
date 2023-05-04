# usage: Rscript run_WGCNA.R binarized_expressions.tsv [deepSplit:0,1,2,3,4] [detectCutHeight:(0-1)]
suppressPackageStartupMessages(library("WGCNA"))

args <- commandArgs(trailingOnly = TRUE)

fileBinExprs <- args[[1]]
fileModules <- paste0(sub(".tsv","",fileBinExprs),".modules.tsv")

deepSplit <- as.integer(args[[2]])
detectCutHeight <- as.numeric(args[[3]])
nt <- "signed hybrid" # networkType = "unsigned", "signed hybrid"

datExpr <- read.csv(fileBinExprs,check.names=FALSE,sep = "\t",header = TRUE,row.names=1)
print(head(datExpr,3))
datExpr[] <- lapply(datExpr, as.numeric)

#### finding power threshold #### 
powers = c(c(1:10), seq(from = 12, to=20, by=2))
rsqcuts <- seq(from = 0.9, to=0.05, by=-0.05)
i = 1
power <- NA
while (is.na(power)){
# Call the network topology analysis function
 sft <- pickSoftThreshold(datExpr, powerVector = powers,verbose = 0,
                          networkType = nt,RsquaredCut =rsqcuts[[i]])  
 power <- sft$powerEstimate
 i<-i+1
}

#print(cat("\nR^2 threshold:",rsqcuts[[i]],"power:",power,"\n"))
print(cat("networkType:",nt,"\n"))
print(cat("maxBlockSize:",dim(datExpr)[[2]]+1,"\n"))

#### find modules ####
net = blockwiseModules(datExpr, power = power,
TOMType = "unsigned", 
networkType = nt, 
minModuleSize = 2,
numericLabels = TRUE,
maxBlockSize = 10000, #dim(datExpr)[[2]]+1,
detectCutHeight = detectCutHeight, #detectCutHeight = 0.995,
#mergeCutHeight = 0.05,
deepSplit = deepSplit,
verbose = 0)

moduleLabels = net$colors

sink(fileModules)
cat(paste0("module","\t","size","\t","genes","\n",sep=' '))
for (i in unique(moduleLabels)){
    cat(paste0(i,"\t",length(names(moduleLabels[moduleLabels==i])),"\t",paste(names(moduleLabels[moduleLabels==i]),collapse=" "),"\n"))
}
sink()

cat(fileModules,"\n")