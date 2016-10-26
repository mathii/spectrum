## Investigate what happens with repeat mutations.
## Type 1 mutations occur at rate a
## Type 2 at rate 1-a
## With probability b, a type 2 mutation is a double mutation.
## Assume that double mutations never happen on the same part of the tree
## That already mutated. 

source("~/rare-var/coalescent.r")

ntrees <- 2000
nsamples <- 100
b <- 0.5
a <- 0.86

spectrum.1 <- rep(0, nsamples)
spectrum.2 <- rep(0, nsamples)

for(i in 1:ntrees){
    t<-simulate.coalescent(100, plot.tree=FALSE, verbose=FALSE)
    freq <- colSums((t$data$seq))
    for(f in freq){spectrum.1[f] <- spectrum.1[f]+1}

    ## Signature 2 mutations
    t<-simulate.coalescent(100, plot.tree=FALSE, verbose=FALSE)
    seq <- t$data$seq
    i=1
    while(i < NCOL(seq)){
        if(runif(1)<b){
            seq[,i] <- as.numeric(seq[,i]|seq[,i+1])
            seq[,i+1] <- 0
            i=i+1
        }
        i=i+1
    }
    freq <- colSums(seq)
    freq <- freq[freq>0]
    for(f in freq){spectrum.2[f] <- spectrum.2[f]+1}
}

par(mfrow=c(1,2))
plot(a*spectrum.1/sum(spectrum.1), col="blue", xlim=c(1,20))
points((1-a)*spectrum.2/sum(spectrum.2), col="red")
lines(a*(1/(1:100))/sum(1/(1:100)), col="blue", lty=2)
lines((1-a)*(1/(1:100))/sum(1/(1:100)), col="red", lty=2)

plot((1-a)*spectrum.2/((1-a)*spectrum.2+a*spectrum.1), xlim=c(0,20), ylim=c(0, 0.4))
abline(h=1-a, lty=2)
