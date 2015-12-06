## Compare the decompositions from PCA, sparse PCA and Non-negative Matrix factorisation.
library("RColorBrewer")
library("NMF")
library("fastICA")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"


tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

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
reg <- info[,3]
names(reg) <- gsub("-", ".", info$ID, fixed=TRUE) 

## Rename and reorder to alexandrov format
alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]
rownames(freq2)<-amap[rownames(freq2)]
freq2<-freq2[amap,]
freq2 <- as.matrix(freq2)

## Subtract off minimal mutations. 
freq2 <- freq2-apply(freq2, 1, min)
freq2["ATA.C",] <- 0.0001         #Null row.

## Non-negative Matrix factorization
nnegmf=nmf(as.matrix(freq2),rank=rank, seed="ica", nrun=20)
pdf(paste0("~/spectrum/plots/","Components_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, ".pdf"), 12, 12)
plot.components(t(coef(nnegmf)), name.map, src, cols=cols, n.components=rank)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_", ifelse(spec=="spectrum", "", paste0(spec, "_")),"NMF.n", n, ".r", rank, tag, ".pdf"), width=12, height=12)
plot.loadings(basis(nnegmf), n.loadings=rank)
dev.off()

## outliers <- c("B_Crete.1", "B_Crete.2"   ,"B_Dai.4"    ,"B_Han.3", "B_Australian.4", "S_Miao.1", "S_Miao.2", "S_Russian.1", "S_Mongola.1")
outliers <- c("B_Mixe.1", "S_Pima.2")
freq2.out <- freq2[,!(colnames(freq2) %in% outliers)]
nnegmf.out=nmf(as.matrix(freq2.out),rank=rank, seed="ica", nrun=20)
pdf(paste0("~/spectrum/plots/","Components_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, "_out.pdf"), 12, 12)
plot.components(t(coef(nnegmf.out)), name.map, src, cols=cols, n.components=rank)
dev.off()

pdf(paste0("~/spectrum/plots/","Loadings_", ifelse(spec=="spectrum", "", paste0(spec, "_")),"NMF.n", n, ".r", rank, tag, "_out.pdf"), width=12, height=12)
plot.loadings(basis(nnegmf.out), n.loadings=rank)
dev.off()

## plot(coef(nnegmf.out)[1,], coef(nnegmf.out)[2, ])
## identify(coef(nnegmf.out)[1,], coef(nnegmf.out)[2, ], colnames(freq2.out))

## nnegmf.test <- nmf(as.matrix(freq2), 2:8, seed="random", nrun=20)
## pdf(paste0("~/spectrum/plots/","Diagnostics_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, tag, ".pdf"), 12, 12)
## plot(nnegmf.test)
## dev.off()
