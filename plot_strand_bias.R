## Plot the results for transcriptional strand bias.

source("~/spectrum/code/spectrumlib.R")

sig <- c("TCT.T", "TCC.T", "CCC.T", "ACC.T")
sig.rc <- c( "AGA.A", "GGA.A", "GGG.A", "GGT.A" )
sig.name <- "1"
in.pops <- "WestEurasia"
out.pops <- c("WestEurasia", "SouthAsia")

## sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
## sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )
## sig.name <- 2
## in.pops <- "America"
## out.pops <- "America"

n <- 2

sig.name.map <- c(4,1)
names(sig.name.map) <- c("1", "2")

## transcribed strand
plus <- read.table(paste0("~/spectrum/strand/chr1.n",n,".strand+.txt"), as.is=TRUE, header=TRUE)
row.names(plus) <- plus[,1]
plus <- plus[,3:NCOL(plus)]
for(chr in 2:22){
pp <- read.table(paste0("~/spectrum/strand/chr", chr, ".n",n,".strand+.txt"), as.is=TRUE, header=TRUE)
row.names(pp) <- pp[,1]
plus <- plus + pp[,3:NCOL(pp)]
}

## untranscribed strand
minus <- read.table(paste0("~/spectrum/strand/chr1.n",n,".strand-.txt"), as.is=TRUE, header=TRUE)
row.names(minus) <- minus[,1]
minus <- minus[,3:NCOL(minus)]
for(chr in 2:22){
mm <- read.table(paste0("~/spectrum/strand/chr", chr, ".n",n,".strand-.txt"), as.is=TRUE, header=TRUE)
row.names(mm) <- mm[,1]
minus <- minus + mm[,3:NCOL(mm)]
}

## Ratio
## ratio <- log((colSums(plus[sig,,drop=FALSE])/colSums(plus[sig.rc,,drop=FALSE]))/(colSums(minus[sig,,drop=FALSE])/colSums(minus[sig.rc,,drop=FALSE])))

ratio <- log((colSums(plus[sig,,drop=FALSE])+colSums(minus[sig.rc,,drop=FALSE]))/(colSums(plus[sig.rc,,drop=FALSE])+colSums(minus[sig,,drop=FALSE])))

## Effect strength
components<-read.table("~/spectrum/plots/Components_ica_NMF.n2.r4.txt", as.is=TRUE)
rownames(components) <- components[,1]
components <- components[2:NCOL(components)]
comp <- components[names(ratio),sig.name.map[sig.name]]

## plot(comp, ratio, col=cols[name.map[names(ratio)]])
## beeswarm(ratio~name.map[names(ratio)], col=cols[unique(name.map[names(ratio)])], pch=16)
mod<-lm(ratio~as.factor(name.map[names(ratio)]))
p <- anova(mod)["Pr(>F)"][[1]][1]

in.ratio <- ratio[name.map[names(ratio)] %in% in.pops]
out.ratio <- ratio[!(name.map[names(ratio)] %in% out.pops)]
pt <- t.test(in.ratio, out.ratio)$p.value

pdf(paste0("~/spectrum/plots/strand_bias.n", n, ".sig", sig.name, ".pdf"))
boxplot(ratio~name.map[names(ratio)], col=cols[unique(name.map[names(ratio)])], pch=16, cex=0.5, ylim=c(-0.3,0.3), bty="n", xlab="", xaxt="n", main=paste0("Signature ", sig.name, " (P=", format(p,digits=2, scientific=p<0.01), " , ", format(pt,digits=2, scientific=pt<0.01), ")"), ylab="Log-ratio strand bias", frame=F) 
mtext(unique(name.map[names(ratio)]), at=1:length(unique(name.map)), las=2, side=1, line=-4)
dev.off()

pdf(paste0("~/spectrum/plots/strand_bias.n", n, ".sig", sig.name, "_subplots.pdf"))
par(mfrow=c(2,2))
par(mar=c(2.1, 4.1, 4.1, 2.1))
for(i in 1:4){
    sg=sig[i]
    sg.rc=sig.rc[i]
    rat <- log((colSums(plus[sg,,drop=FALSE])+colSums(minus[sg.rc,,drop=FALSE]))/(colSums(plus[sg.rc,,drop=FALSE])+colSums(minus[sg,,drop=FALSE])))
    agg <- aggregate(rat, FUN=mean, by=list(name.map[names(rat)]))
    agg.sd <- aggregate(rat, FUN=sd, by=list(name.map[names(rat)]))
    agg.n <- aggregate(rat, FUN=length, by=list(name.map[names(rat)]))
    agg.sd$x <- agg.sd$x/sqrt(agg.n$x)

    mod<-lm(rat~as.factor(name.map[names(rat)]))
    p <- anova(mod)["Pr(>F)"][[1]][1]

    in.rat <- rat[name.map[names(rat)] %in% in.pops]
    out.rat <- rat[!(name.map[names(rat)] %in% out.pops)]
    pt <- t.test(in.rat, out.rat)$p.value
    
    plot(agg$x, col=cols[agg$Group.1], pch=16, cex=2, bty="n", xaxt="n", ylim=c(-0.3, 0.3), ylab="", xlab="", main=paste0(gsub( ".", ">", sg, fixed=TRUE), " (P=", format(p,digits=2, scientific=pt<0.01), " , ", format(pt,digits=2, scientific=p<0.01), ")"))
    segments(1:NROW(agg), agg$x-2*agg.sd$x, y1=agg$x+2*agg.sd$x,  col=cols[agg.sd$Group.1], lwd=1)
    segments(1:NROW(agg), agg$x-agg.sd$x, y1=agg$x+agg.sd$x,  col=cols[agg.sd$Group.1], lwd=2)
    
    abline(h=0)
}
dev.off()

