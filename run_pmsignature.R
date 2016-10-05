## Run the analysis using the pmsignature package.
library(pmsignature)
source("~/spectrum/code/spectrumlib.R")

n <- 2
rank <- 4
G <- readMPFile(paste0("~/spectrum/mpf/mpf.all.n",n, ".txt.gz"), numBases = 3)
BG_prob <- readBGFile(G)
Param <- getPMSignature(G, K = rank+1, BG = BG_prob)

components <- getMembershipValue(Param)[,1:rank]
rownames(components) <- gsub("-", ".", rownames(components))

write.table( Components, paste0("~/spectrum/plots/Components_pmsignature", ".n", n, ".r", rank, ".txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )

pdf(paste0("~/spectrum/plots/","Components_pmsignature.n", n, ".r", rank, tag, ".pdf"), 12, 12)
plot.components(components, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(4,3,1,2), yploti=c(1,2,3,4))
dev.off()
