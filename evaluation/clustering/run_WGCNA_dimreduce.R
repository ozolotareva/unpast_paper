# usage: Rscript run_WGCNA.R binarized_expressions.tsv
suppressPackageStartupMessages(library("WGCNA"))

args <- commandArgs(trailingOnly = TRUE)

# fileBinExprs <- '/Users/fernando/Documents/Research/DESMOND2/DESMOND2/data/simulated/A/example.tsv'

fileBinExprs <- args[[1]]

datExpr <- read.table(fileBinExprs, check.names=FALSE, sep = '\t', row.names = 1, header = T)
datExpr <- as.data.frame(t(datExpr))
datExpr[] <- lapply(datExpr, as.numeric)

#### finding power threshold ####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
rsqcuts <- seq(from = 0.9, to=0.05, by=-0.05)
i = 1
power <- NA

while (is.na(power)){
  # Call the network topology analysis function
  sft <- pickSoftThreshold(datExpr, powerVector = powers,verbose = 0, networkType = "signed hybrid", RsquaredCut =rsqcuts[[i]])
  power <- sft$powerEstimate
  i<-i+1
}

# modules are based on genes, not patients?
cat("\nR^2 threshold:", rsqcuts[[i]], "power:", power,"\n")

#### find modules ####
net <- blockwiseModules(datExpr, power = power,
                        TOMType = "signed", networkType = "signed hybrid", minModuleSize = 2,
                        numericLabels = TRUE,
                        #detectCutHeight = 1-10^(-1*power),
                        #mergeCutHeight = 0.05,
                       verbose = 0)


# moduleLabels = net$colors
# colmodules <- as.data.frame(colmodules)
# colmodules['gene'] <- rownames(colmodules)
# colmodules['label'] <- colmodules$colmodules - 1
# colmodules <- colmodules[,c('gene', 'label')]
# write.table(colmodules, paste0(sub(".tsv","",fileBinExprs),"_modules_genes.tsv"), sep='\t', quote = F, row.names = F, col.names = T)

write.table(as.data.frame(t(net$MEs)), sub(".tsv","_MEs.tsv", fileBinExprs), sep='\t', quote = F, row.names = T, col.names = T)

