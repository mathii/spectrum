## Find the line that separates cell lines from non cell lines on the PCA.

library("RColorBrewer")
source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/spectrum_matrix.n", n,tag, ".txt")

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


freq2 <- freq2[rownames(freq2)!="ATA.C",]
freq2 <- freq2-apply(freq2, 1, min)

cell.line <- ifelse(info$Source=="Genomic_from_cell_lines", "CellLine", "NonCellLine")
names(cell.line) <- gsub("-", ".", info$ID) 
cell.line <- cell.line[colnames(freq2)]
## PCA
pca=prcomp(t(freq2), center=TRUE, scale.=TRUE)

pdf(paste0("~/spectrum/plots/CellLine_PCA.n",n, ".pdf"), 12, 6)
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", col=ifelse(cell.line=="CellLine", "red", "blue"), pch=ifelse(cell.line=="CellLine", 13, 1))
legend("bottomright", c("Cell line sample", "Primary tissue sample"), col=c("red", "blue"), pch=c(13,1), bty="n")
dev.off()
