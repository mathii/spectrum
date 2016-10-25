## Analyse the simulated spectrum
library(RcppCNPy)

sample.sizes <- c(100,100,100,40)                 #AFR, EUR, ASN, AMR
ns <- 1:20
weights <- c(80,20)                     #nonCpg / CpG weights
cols <- c("#377EB8" , "#4DAF4A" , "#E41A1C" , "#984EA3" )

## nonCpG <- read.fwf("~/spectrum/sims/nonCpG_spectrum.gz", rep(1, sum(sample.sizes)))
## CpG <- read.fwf("~/spectrum/sims/CpG_spectrum.gz", rep(1, sum(sample.sizes)))

cat("Loading\n")
nonCpG <- npyLoad("/Users/mathii/spectrum/sims/nonCpG_spectrum.npy")
CpG <- npyLoad("/Users/mathii/spectrum/sims/CpG_spectrum.npy")
cat("Analysing\n")

get.include.mutations <- function(muts, ss){
  inc <- as.data.frame(matrix(0, nrow=NROW(muts), ncol=length(ss)))
  brk <- c(0, cumsum(ss))
  for(i in 1:length(ss)){
    inc[,i] <- apply(muts[,(brk[i]+1):(brk[i+1])], 1, any)
  }
  return(inc)
}

inc.nonCpG <- get.include.mutations(nonCpG, sample.sizes )
inc.CpG <- get.include.mutations(CpG, sample.sizes )

proportions <- as.data.frame(matrix(0, nrow=length(ns), ncol=length(sample.sizes)))
for(i in 1:length(ns)){
  inc.n.cpg <- rowSums(CpG)==ns[i]
  inc.n.notCpG <- rowSums(nonCpG)==ns[i]
  for(j in 1:length(sample.sizes)){
    proportions[i,j] <- weights[2]*sum(inc.n.cpg & inc.CpG[,j])/(weights[2]*sum(inc.n.cpg & inc.CpG[,j])+weights[1]*sum(inc.n.notCpG & inc.nonCpG[,j]))
  }
}

plot(ns, proportions[,1], pch=16, cex=0.5, col=cols[1], ylim=c(0.1,0.3), xlab="Allele count", ylab="Proportion of signature mutations", xlim=c(0,10))
inc<-!is.na(proportions[,1])
lines(smooth.spline(ns[inc], proportions[inc,1], spar=0.25), col=cols[1])
for(i in 2:length(sample.sizes)){
    points(ns, proportions[,i], pch=16, cex=0.5,  col=cols[i])
    inc<-!is.na(proportions[,i])
    lines(smooth.spline(ns[inc], proportions[inc,i], spar=0.25), col=cols[i])
}
