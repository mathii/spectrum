## TCT.T/TCA.T ratio for ancients..
library("RColorBrewer")
## library("NMF")
## library("fastICA")
library("beeswarm")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 3
spec <- "spectrum"
r1 <- "TCC.T"
r2 <- "TCA.T"

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
## freq2 <- freq2[,!( colnames(freq2) %in% ancient)]

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

ratio <- freq2[r1,]/freq2[r2,]

pdf(paste0("~/spectrum/plots/","Ratio_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".", r1, "-", r2, ".pdf"), 6, 6)
aa <- ancient
names(aa) <- ancient
regplot <- c(reg, aa)
ro <-  c(  "Africa", "Oceania", "EastAsia",  "CentralAsiaSiberia", "America", "SouthAsia", "WestEurasia", "Loschbour" , "LBK1b_leipzig_v2", "Ust_Ishim", "AltaiNeandertal", "Denisova")
regplot <- factor(regplot,ro)
cc <- c(cols[ro[1:7]], rep("black", 5))
beeswarm(log2(ratio)~regplot[names(ratio)], cex.axis=0.5, las=2, col=cc, pch=16, cex=0.5, xlab="", ylab="Log2 Ratio", bty="n", ylim=c(-1,1))
## bxplot(log2(ratio)~regplot[names(ratio)], add = TRUE, col=c("grey", "black", "grey"), cex=0.5)             
dev.off()
               
