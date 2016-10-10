## Plot the proportion of variants in each class, as a function of n, as fn increases

ns <- 2:30

source("~/spectrum/code/spectrumlib.R")

####################################################

## sig.name <- 1
## sig <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
## ylim=c(0.065, 0.105)
## spr<-0.5
## wts<-c(10,20,10,5,rep(1,length(ns)-4))

sig.name <- 2
sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
ylim=c(0.14, 0.22)
in.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")
spr <- 0.25
wts<-c(10,20,10,10,10,20,rep(1,length(ns)-6))

####################################################

sig.name.map <- c(4,1)
sig.map=c(2,4,3,1)

####################################################

regions <- sort(unique(info$Region))
ltys <- rep(1, length(cols))
names(ltys) <- names(cols)
if(sig.name==2){
    regions <- c("Africa", "America_Hi", "America_Lo", "CentralAsiaSiberia", "EastAsia", "Oceania", "SouthAsia", "WestEurasia")
    hi.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")
    cols <- c(cols, "America_Hi"="#984EA3", "America_Lo"="#984EA3")
    ltys <- rep(1, length(cols))
    names(ltys) <- names(cols)
    ltys["America_Lo"] <- 2
}

proportions <- counts <- as.data.frame(matrix(0, nrow=length(ns), ncol=length(regions)))
names(proportions) <- names(counts) <- regions
rownames(info) <- gsub("-", ".", info$ID)

if(sig.name==2){
    info[hi.ind,"Region"] <- "America_Hi"
    info[info[rownames(info),"Region"]=="America","Region"] <- "America_Lo"
}

for(i in 1:length(ns)){
    n <- ns[i]
    data <- read.table(paste0("~/spectrum/data/spectrum_matrix.n", n, ".txt"), as.is=TRUE, header=TRUE)
    freq <- t(t(data)/colSums(data))
    cnts <- read.table(paste0("~/spectrum/data/count_matrix.n", n, ".txt"), as.is=TRUE, header=TRUE)

    for(reg in regions){
        proportions[i,reg] <-  mean(colSums(freq[sig,info[colnames(freq),"Region"] %in% reg]))
        counts[i,reg] <-  sum(colSums(cnts[sig,info[colnames(cnts),"Region"] %in% reg]))
    }
}

pdf(paste0("~/spectrum/plots/fn_sig", sig.name, ".pdf"), width=8, height=5)
plot(ns, proportions[,regions[1]], pch=16, cex=0.5, col=cols[regions[1]], ylim=ylim, xlab="Allele count", ylab=bquote("Proportion of signature"~.(sig.name)~f[2]~"mutations"))
lines(smooth.spline(ns, proportions[,regions[1]], spar=spr, w=wts), col=cols[regions[1]], lty=ltys[regions[1]])
for(i in 2:length(regions)){
    points(ns, proportions[,regions[i]], pch=16, cex=0.5,  col=cols[regions[i]])
lines(smooth.spline(ns, proportions[,regions[i]], spar=spr, w=wts), col=cols[regions[i]], lty=ltys[regions[i]])
}
legend("topright", regions, cex=0.5, col=cols[regions], lty=ltys[regions], pch=16, bty="n", ncol=2)
dev.off()