## Compare the decompositions from PCA, sparse PCA and Non-negative Matrix factorisation.
library("RColorBrewer")
## library("nsprcomp")
## library("NMF")
source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2

source("~/spectrum/code/spectrumlib.R")

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/spectrum_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )

freq2 <- freq2[rownames(freq2)!="ATA.C",]
## freq2 <- freq2-apply(freq2, 1, min)

## PCA
pca=prcomp(t(freq2), center=TRUE, scale.=TRUE)
## pca=princomp(t(freq2), center=TRUE, scale.=TRUE)


pdf(paste0("~/spectrum/plots/","Components_PCA.n", n, tag, ".pdf"), 12, 12)
plot.components(pca$x, name.map, src, cols=cols, n.components=5)
## plot.components(pca$scores, name.map, src, cols=cols, n.components=5)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_PCA.n", n, tag, ".pdf"), width=12, height=12)
plot.loadings(pca$rotation, rescale=FALSE)
## plot.loadings(loadings(pca))
dev.off()

## ## Sparse PCA
## spca=nsprcomp(t(freq2), center=TRUE, scale.=TRUE, k=10, ncomp=6)
## pdf(paste0("~/spectrum/plots/","Components_SPCA.n", n, tag, ".pdf"), 12, 12)
## plot.components(spca$x, name.map, src, cols=cols)
## dev.off()

## pdf(paste0("~/spectrum/plots/","Loadings_SPCA.n", n, tag, ".pdf"), width=12, height=12)
## plot.loadings(spca$rotation)
## dev.off()

## ## Non-negative Sparse PCA
## nspca=nsprcomp(t(freq2), center=TRUE, scale.=TRUE, k=10, ncomp=6, nneg=TRUE)
## pdf(paste0("~/spectrum/plots/","Components_NSPCA.n", n, tag, ".pdf"), 12, 12)
## plot.components(nspca$x, name.map, src, cols=cols)
## dev.off()

## pdf(paste0("~/spectrum/plots/","Loadings_NSPCA.n", n, tag, ".pdf"), width=12, height=12)
## plot.loadings(nspca$rotation)
## dev.off()

## Non-negative Matrix factorization
## rank=4
## nnegmf=nmf(as.matrix(freq2),rank=rank, nrun=20)
## pdf(paste0("~/spectrum/plots/","Components_NMF.n", n, tag, ".pdf"), 12, 12)
## plot.components(t(coef(nnegmf)), name.map, src, cols=cols, n.components=rank)
## dev.off()

## pdf(paste0("~/spectrum/plots/","Loadings_NMF.n", n, tag, ".pdf"), width=12, height=12)
## plot.loadings(basis(nnegmf), n.loadings=rank)
## dev.off()

