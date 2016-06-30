## Plot the European fold-excess as a function of B statistic.

source("~/spectrum/code/spectrumlib.R")

n <- 2

bs <- 0:9
excess=-1+0*bs

sig <- c("TCT.T", "TCC.T", "CCC.T", "ACC.T")
sig.rc <- c( "AGA.A", "GGA.A", "GGG.A", "GGT.A" )
in.pops <- c("WestEurasia")
out.pops <- c("WestEurasia", "SouthAsia")

## sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
## sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )
## in.pops <- c("America")
## out.pops <- c("America")


for(bi in 1:10){
    data <- read.table(paste0("~/spectrum/Bstat/chr1.n",n,".b", bs[bi], ".txt"), as.is=TRUE, header=TRUE)
    row.names(data) <- data[,1]
    data <- data[,3:NCOL(data)]
    for(chr in 2:22){
        pp <- read.table(paste0("~/spectrum/Bstat/chr", chr ,".n",n,".b", bs[bi], ".txt"), as.is=TRUE, header=TRUE)
        row.names(pp) <- pp[,1]
        data <- data + pp[,3:NCOL(pp)]
    }

    proportion <- colSums(data[c(sig, sig.rc),])/colSums(data)
    excess[bi] <- mean(proportion[name.map[names(proportion)]%in% in.pops])-mean(proportion[!(name.map[names(proportion)]%in% out.pops)])
    
}

plot(excess, type="b", pch=16)
