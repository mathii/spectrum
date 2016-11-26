## TCT.T/TCA.T ratio for ancients..
library("RColorBrewer")
## library("NMF")
## library("fastICA")
library("beeswarm")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
spec <- "spectrum"
what.anc <- "het_sgdp_f1"

r1 <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
r2 <- c("TCA.T", "ACA.T", "TCA.T", "CCA.T")
## r2 <- c("ATA.C","ATA.C","ATA.C","ATA.C")

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
## ancient <- c("Altai", "Denisova", "Loschbour", "LBK", "UstIshim", "Vindija")
## ancientmap <- c("Altai Neandertal", "Denisova", "Loschbour", "Stuttgart", "Ust' Ishim", "Vindija Neandertal")
ancient <- c("Altai", "Denisova", "Loschbour", "LBK", "UstIshim")
ancientmap <- c("Altai Neandertal", "Denisova", "Loschbour", "Stuttgart", "Ust' Ishim")

names(ancientmap) <- ancient

freq2 <- as.matrix(freq2)

ratio <- freq2[r1,,drop=FALSE]/freq2[r2,,drop=FALSE]
logratio <- apply(log2(ratio), 2, sum)
lm <- mean(logratio[reg[names(logratio)]=="Africa"], na.rm=TRUE)
ls <- sd(logratio[reg[names(logratio)]=="Africa"], na.rm=TRUE)
logratio <- (logratio-lm)/ls

ancient.sd <- rep(0, length(ancientmap))
ancient.mn <- rep(0, length(ancientmap))
ancient.lq <- rep(0, length(ancientmap))
ancient.uq <- rep(0, length(ancientmap))

names(ancient.sd) <- names(ancient.mn) <- ancientmap
names(ancient.lq) <- names(ancient.uq) <- ancientmap

for(i in 1:length(ancient)){
    adata <- read.table(paste0("~/spectrum/ancient/", ancient[i], ".", what.anc, ".counts.txt"), as.is=TRUE, header=TRUE)
    cmpnames <- sapply(adata[,1], reverse.complement)
    adata[,1] <- ifelse(adata[,1] %in% rownames(freq2), adata[,1], cmpnames)
    adata<-aggregate(adata[,2], by=list(adata[,1]), sum)

    nbs <- 1000
    bs <- rep(0,nbs)
    for(j in 1:nbs){
      bsa <- tabulate(sample(96, sum(adata[,2]), prob=adata[,2], replace=T), nbins=96)
      names(bsa) <- adata[,1]
      bs[j] <- (sum(log2(bsa[r1]/bsa[r2]))-lm)/ls
    }
    ancient.sd[i] <- sd(bs)
    ancient.mn[i] <- mean(bs)
    qq <- quantile(bs, c(0.05, 0.95))
    ancient.lq[i] <- qq[1]
    ancient.uq[i] <- qq[2]
    
    nm <- adata[,1]
    adata <- adata[,2]/adata[adata[,1]=="ATA.C",2]
    names(adata) <- nm
    adata <- adata[rownames(freq2)]
    logratio <- c(logratio, (sum(log2(adata[r1]/adata[r2]))-lm)/ls )
    names(logratio)[length(logratio)] <- ancientmap[i]
}


## plot(rowMeans(freq2)/sum(rowMeans(freq2)), adata)
## abline(0,1, col="red")
## identify(rowMeans(freq2)/sum(rowMeans(freq2)), adata, rownames(freq2))


pdf(paste0("~/spectrum/plots/","Corrected_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".", r1, "-", r2, ".", what.anc, ".pdf"), 6, 6)
aa <- ancientmap[ancient]
names(aa) <- ancientmap[ancient]
regplot <- c(reg, aa)
## anc.ro <- c( "Loschbour" , "Stuttgart", "Ust' Ishim", "Altai Neandertal", "Vindija Neandertal", "Denisova" )
anc.ro <- c( "Loschbour" , "Stuttgart", "Ust' Ishim", "Altai Neandertal", "Denisova" )

mod.ro <-  c(  "Africa", "Oceania", "EastAsia",  "CentralAsiaSiberia", "America", "SouthAsia", "WestEurasia")
ro <- c(mod.ro ,anc.ro)
regplot <- factor(regplot,ro)
cc <- c(cols[ro[1:7]], rep("black", length(ancient)))
beeswarm(logratio~regplot[names(logratio)], cex.axis=0.6, las=2, col=cc, pch=16, cex=0.7, xlab="", ylab=expression("M"[2]), bty="n")

for(i in 1:length(anc.ro)){
  x <- length(mod.ro) + i
  lq <- ancient.lq[anc.ro[i]]
  uq <- ancient.uq[anc.ro[i]]
  segments(x, lq, x, uq, lwd=3)
  points(x, logratio[anc.ro[i]], cex=1, pch=16)
}

## bxplot(logratio~regplot[names(logratio)], add = TRUE, col=c("grey", "black", "grey"), cex=0.5)             
dev.off()

### Does it correlate with ancestry estimtes               
## lr.modern <- logratio[!(names(logratio) %in% ancientmap)]
## pop.modern <-gsub('.[0-9]', "", gsub("^[A-Z]*_", "", names(lr.modern)))
## laz.prop <- read.table("~/spectrum/data/lazaridis_proportions.txt", as.is=TRUE, header=TRUE)
## rownames(laz.prop) <- laz.prop$POP
## laz.prop <- laz.prop[,c(2:4)]

## include <- pop.modern %in% rownames(laz.prop)
## lr.modern.include <- lr.modern[include]
## pop.modern.include <- pop.modern[include]

## mod.eef <- lm(lr.modern.include~laz.prop[pop.modern.include, "EEF"])
## mod.whg <- lm(lr.modern.include~laz.prop[pop.modern.include, "WHG"])
## mod.ane <- lm(lr.modern.include~laz.prop[pop.modern.include, "ANE"])

## #Does it correlate with any pca?
## inc <-  which(reg[names(logratio)]=="WestEurasia")
## we.logratio <- logratio[inc]
## we.freq2 <- freq2[,inc]

## we.freq2 <- we.freq2[rownames(we.freq2)!="ATA.C",]

prop2<-t(t(freq2)/colSums(freq2))
