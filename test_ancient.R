## TCT.T/TCA.T ratio for ancients..
library("RColorBrewer")
## library("NMF")
## library("fastICA")
library("beeswarm")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
spec <- "spectrum"
r1 <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
r2 <- c("TCA.T", "ACA.T", "TCA.T", "CCA.T")

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
## freq2 <- freq2[,!( colnames(freq2) %in% ancient)]

ancientmap <- c("Neandertal", "Denisova", "Loschbour", "Stuttgart", "Ust' Ishim")
names(ancientmap) <- ancient

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
logratio <- apply(log2(ratio), 2, sum)
lm <- mean(logratio[reg[names(logratio)]=="Africa"], na.rm=TRUE)
ls <- sd(logratio[reg[names(logratio)]=="Africa"], na.rm=TRUE)
logratio <- (logratio-lm)/ls

for(i in 1:length(ancientmap)){
    names(logratio)[names(logratio)==names(ancientmap)[i]] <- ancientmap[i]
}

pdf(paste0("~/spectrum/plots/","Ratio_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".", r1, "-", r2, ".pdf"), 6, 6)
aa <- ancientmap[ancient]
names(aa) <- ancientmap[ancient]
regplot <- c(reg, aa)
ro <-  c(  "Africa", "Oceania", "EastAsia",  "CentralAsiaSiberia", "America", "SouthAsia", "WestEurasia", "Loschbour" , "Stuttgart", "Ust' Ishim", "Neandertal", "Denisova")
regplot <- factor(regplot,ro)
cc <- c(cols[ro[1:7]], rep("black", 5))
beeswarm(logratio~regplot[names(logratio)], cex.axis=0.5, las=2, col=cc, pch=16, cex=0.5, xlab="", ylab=expression("M"[2]), bty="n", ylim=c(-5, 15))
## bxplot(logratio~regplot[names(logratio)], add = TRUE, col=c("grey", "black", "grey"), cex=0.5)             
dev.off()

### Does it correlate with ancestry estimtes               
lr.modern <- logratio[!(names(logratio) %in% ancientmap)]
pop.modern <-gsub('.[0-9]', "", gsub("^[A-Z]*_", "", names(lr.modern)))
laz.prop <- read.table("~/spectrum/data/lazaridis_proportions.txt", as.is=TRUE, header=TRUE)
rownames(laz.prop) <- laz.prop$POP
laz.prop <- laz.prop[,c(2:4)]

include <- pop.modern %in% rownames(laz.prop)
lr.modern.include <- lr.modern[include]
pop.modern.include <- pop.modern[include]

mod.eef <- lm(lr.modern.include~laz.prop[pop.modern.include, "EEF"])
mod.whg <- lm(lr.modern.include~laz.prop[pop.modern.include, "WHG"])
mod.ane <- lm(lr.modern.include~laz.prop[pop.modern.include, "ANE"])

#Does it correlate with any pca?
inc <-  which(reg[names(logratio)]=="WestEurasia")
we.logratio <- logratio[inc]
we.freq2 <- freq2[,inc]

we.freq2 <- we.freq2[rownames(we.freq2)!="ATA.C",]
