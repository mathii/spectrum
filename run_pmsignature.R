## Run the analysis using the pmsignature package.
library(pmsignature)

n <- 2
rank <- 4
G <- readMPFile(paste0("~/spectrum/mpf/mpf.all.n",n, ".txt.gz"), numBases = 3)
BG_prob <- readBGFile(G)
Param <- getPMSignature(G, K = rank, BG = BG_prob)

Components <- getMembershipValue(Param)
write.table( Components, paste0("~/spectrum/plots/Components_pmsignature", ".n", n, ".r", rank, ".txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )
