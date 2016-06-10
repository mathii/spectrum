source("~/spectrum/code/spectrumlib.R")

cosmic <- read.table("~/spectrum/cosmic data//signatures_probabilities.txt", as.is = TRUE, header = TRUE, sep="\t")

sigs <- c(1,11,8)

pdf(paste0("~/spectrum/cosmic data/cosmic.pdf"), width=12, height=3*length(sigs))
par(mfrow = c(2,1))
ms <- cosmic[,paste0("Signature.", sigs)]
rownames(ms) <- paste0(cosmic$Trinucleotide,".", substr(cosmic$Substitution.Type,3,3))
ms <- as.matrix(ms)

alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]

ms <- ms[alex[,2],]

plot.loadings(ms, n.loadings = length(sigs))
dev.off()