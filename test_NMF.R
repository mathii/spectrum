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

write.table(t(coef(nnegmf)), paste0("~/spectrum/plots/","Components_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, ".txt"), row.names=T, col.names=F, quote=F) 
write.table(basis(nnegmf), paste0("~/spectrum/plots/","Loadings_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, ".txt"), row.names=T, col.names=F, quote=F) 

pdf(paste0("~/spectrum/plots/","Loadings_", ifelse(spec=="spectrum", "", paste0(spec, "_")),"NMF.n", n, ".r", rank, tag, "_out.pdf"), width=12, height=12)
plot.loadings(basis(nnegmf.out), n.loadings=rank)
dev.off()



## plot(coef(nnegmf)[3,], coef(nnegmf)[4, ])
## identify(coef(nnegmf)[3,], coef(nnegmf)[4, ], colnames(freq2.out))

## nnegmf.test <- nmf(as.matrix(freq2), 2:8, seed="random", nrun=20)
## pdf(paste0("~/spectrum/plots/","Diagnostics_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, tag, ".pdf"), 12, 12)
## plot(nnegmf.test)
## dev.off()
