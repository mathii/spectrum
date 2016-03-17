## Analyse the 1kg data. doubletons.
library(beeswarm)

cols <- c("goldenrod", "peru", "navajowhite3", "lightgoldenrod3", "palegoldenrod", "brown3", "orangered", "sienna1", "brown", "tomato3", "orange", "green3", "greenyellow", "green", "green4", "limegreen", "blue", "turquoise3", "cyan", "cornflowerblue", "darkblue", "darkmagenta", "darkviolet", "maroon", "deeppink2", "magenta")
pops <- c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")
names(cols) <- pops
cols<-cols[order(names(cols))]

n <- 2

exclude <- c("HG01149", "NA20582", "NA12275", "NA19728", "NA20540")

## Include reverse complements here because we don't rc when generating. 
sig1 <- c("TCT.T", "TCC.T", "CCC.T", "ACC.T", "AGA.A", "GGA.A", "GGG.A", "GGT.A" )
sig2 <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T", "CGT.A", "CGG.A", "CGC.A", "CGA.A" )
    
chr <- 1
data <- read.table(paste0("~/spectrum/1kg/chr", chr, ".n", n), header=T)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]
for(chr in 2:22){
    d <- read.table(paste0("~/spectrum/1kg/chr", chr, ".n", n), header=T)
    rownames(d) <- d[,1]
    data <- data + d[,2:NCOL(d)]
}

for(n in 3){
    for(chr in 2:22){
        d <- read.table(paste0("~/spectrum/1kg/chr", chr, ".n", n), header=T)
        rownames(d) <- d[,1]
        data <- data + d[,2:NCOL(d)]
    }
}

data <- t(t(data)/colSums(data))

data <- data[,!(colnames(data)%in%exclude)]

panel <- read.table("~/spectrum/code/1kgP3_panel.txt", as.is=TRUE)
pmap <- panel[,2]
names(pmap) <- panel[,1]

s1 <- colSums(data[sig1,])

pop.means <- aggregate(s1, by=list(pmap[names(s1)]), mean)
pm1 <- pop.means[,2]
names(pm1) <- pop.means[,1]
pm1 <- sort(pm1)


s2 <- colSums(data[sig2,])
pop.means <- aggregate(s2, by=list(pmap[names(s2)]), mean)
pm2 <- pop.means[,2]
names(pm2) <- pop.means[,1]
pm2 <- sort(pm2)

pdf(paste0("~/spectrum/plots/1kg_sig1.","pdf"), width=12, height=4)
par(mar=c(5.1,4.1,2.1,2.1))
beeswarm(s1~pmap[names(s1)], pch=16, cex=0.3, col=cols, at=order(names(pm1)), xlab="Population", ylab="Signature 1", cex.axis=0.6)
dev.off()

pdf(paste0("~/spectrum/plots/1kg_sig2.", "pdf"), width=12, height=4)
par(mar=c(5.1,4.1,2.1,2.1))
beeswarm(s2~pmap[names(s2)], pch=16, cex=0.3, col=cols, at=order(names(pm2)), xlab="Population", ylab="Signature 2", cex.axis=0.6)
dev.off()

