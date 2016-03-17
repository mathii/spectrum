## Test overall rate for ancients

library("RColorBrewer")
## library("NMF")
## library("fastICA")
library("beeswarm")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
spec <- "spectrum"
what.anc <- "het_sgdp_f1"

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("Altai", "Denisova", "Loschbour", "LBK", "UstIshim")

ancientmap <- c("Neandertal", "Denisova", "Loschbour", "Stuttgart", "Ust' Ishim")
names(ancientmap) <- ancient

freq2 <- as.matrix(freq2)
nc <- NCOL(freq2)

for(i in 1:length(ancient)){
    adata <- read.table(paste0("~/spectrum/ancient/", ancient[i], ".", what.anc, ".counts.txt"), as.is=TRUE, header=TRUE)
    cmpnames <- sapply(adata[,1], reverse.complement)
    adata[,1] <- ifelse(adata[,1] %in% rownames(freq2), adata[,1], cmpnames)
    adata<-aggregate(adata[,2], by=list(adata[,1]), sum)
    
    nm <- adata[,1]
    adata <- adata[,2]/adata[adata[,1]=="ATA.C",2]
    names(adata) <- nm
    adata <- adata[rownames(freq2)]
    freq2 <- cbind(freq2, adata)
    colnames(freq2)[NCOL(freq2)] <- ancientmap[i]
}

freq2 <- freq2[rownames(freq2)!="ATA.C",]

pca=prcomp(t(freq2[,1:nc]), center=TRUE, scale.=TRUE)

acols <- "black"
names(acols) <- "Ancient (projected)"
asrc <- ancientmap
asrc[1:length(asrc)] <- "Ancient (projected)"
names(asrc) <- ancientmap

ppca<-predict(pca, t(freq2))

pdf(paste0("~/spectrum/plots/","Ancient_Components.",what.anc,".PCA.n", n, tag, ".pdf"), 12, 12)
plot.components(ppca, c(name.map, asrc) , c(src, asrc), cols=c(cols, acols) , n.components=5)
dev.off()

## pdf(paste0("~/spectrum/plots/","Ancient_Loadings_PCA.n", n, tag, ".pdf"), width=12, height=12)
## plot.loadings(ppca$rotation, rescale=FALSE)
## dev.off()
