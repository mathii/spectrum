## Compare the decompositions from PCA, sparse PCA and Non-negative Matrix factorisation.
library("RColorBrewer")
library("nsprcomp")
library("NMF")
source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 1

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/spectrum_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
freq2 <- freq2[,!( colnames(freq2) %in% ancient)]

## Boilerplate - load info
info <- read.table("~/spectrum/data/location_info.txt", as.is=TRUE, header=TRUE, sep="\t")
info[info[,6]=="Genomic from saliva",6]<-"Genomic_from_saliva"
info[info[,6]=="?",6]<-"Unknown"
regions <- unique(info[,3])
cols <- brewer.pal(length(regions), "Set1")
cols[6] <- "darkgrey"
sources <- unique(info[,6])
source.cols <- brewer.pal(length(sources), "Set2")
names(cols) <- regions
names(source.cols) <- sources
name.map <- info[,3]
names(name.map) <- gsub("-", ".", info$ID, fixed=TRUE)
src <- info[,6]
names(src) <- gsub("-", ".", info$ID, fixed=TRUE)

## Rename and reorder to alexandrov format
alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]
rownames(freq2)<-amap[rownames(freq2)]
freq2<-freq2[amap,]
freq2 <- as.matrix(freq2)


freq2 <- freq2[rownames(freq2)!="ATA.C",]
freq2 <- freq2-apply(freq2, 1, min)

## PCA
pca=prcomp(t(freq2), center=TRUE, scale.=TRUE)

pdf(paste0("~/spectrum/plots/","Components_PCA.n", n, tag, ".pdf"), 12, 12)
plot.components(pca$x, name.map, src, cols=cols, n.components=5)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_PCA.n", n, tag, ".pdf"), width=12, height=12)
plot.loadings(pca$rotation)
dev.off()

## Sparse PCA
spca=nsprcomp(t(freq2), center=TRUE, scale.=TRUE, k=10, ncomp=6)
pdf(paste0("~/spectrum/plots/","Components_SPCA.n", n, tag, ".pdf"), 12, 12)
plot.components(spca$x, name.map, src, cols=cols)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_SPCA.n", n, tag, ".pdf"), width=12, height=12)
plot.loadings(spca$rotation)
dev.off()

## Non-negative Sparse PCA
nspca=nsprcomp(t(freq2), center=TRUE, scale.=TRUE, k=10, ncomp=6, nneg=TRUE)
pdf(paste0("~/spectrum/plots/","Components_NSPCA.n", n, tag, ".pdf"), 12, 12)
plot.components(nspca$x, name.map, src, cols=cols)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_NSPCA.n", n, tag, ".pdf"), width=12, height=12)
plot.loadings(nspca$rotation)
dev.off()

## Non-negative Matrix factorization
## rank=4
## nnegmf=nmf(as.matrix(freq2),rank=rank, nrun=20)
## pdf(paste0("~/spectrum/plots/","Components_NMF.n", n, tag, ".pdf"), 12, 12)
## plot.components(t(coef(nnegmf)), name.map, src, cols=cols, n.components=rank)
## dev.off()

## pdf(paste0("~/spectrum/plots/","Loadings_NMF.n", n, tag, ".pdf"), width=12, height=12)
## plot.loadings(basis(nnegmf), n.loadings=rank)
## dev.off()

