## Investigate what happens with repeat mutations.
## Type 1 mutations occur at rate a
## Type 2 at rate 1-a
## With probability b, a type 2 mutation is a double mutation.
## Assume that double mutations never happen on the same part of the tree
## That already mutated. 
set.seed(12345)

source("~/rare-var/kpop.R")

simulate.coalescent.norec<-function(sample=10, theta=10, sites=10, nrun=1, scaling=function(x){x},
	plot.tree=TRUE, return=TRUE, n.tree=5, all.seg=TRUE, inf.sites=TRUE, verbose=TRUE) {

	#No recombination
	nodes<-matrix(nrow=2*sample-1, ncol=5);
	colnames(nodes)<-c("Name", "Time", "D1", "D2", "Mutations");

	nodes[,1]<-c(1:nrow(nodes));
	nodes[,2:5]<-0;

	k<-sample;
	klist<-c(1:sample);
	time<-0;
	current.node<-k+1;

	while(k>1) {
		rate<-k*(k-1)/2;
		dt<--log(runif(1))/rate;
		time<-time+dt;
		l1<-ceiling(runif(1)*k);
		tmp<-klist[l1];
		klist[l1]<-klist[k];
		klist[k]<-tmp;
		l2<-ceiling(runif(1)*(k-1));
		nodes[current.node,2]<-time;
		nodes[current.node,3]<-klist[k];
		nodes[current.node,4]<-klist[l2];
		klist[l2]<-current.node;
		current.node<-current.node+1;
		k<-(k-1);
		klist<-klist[1:k];
	}

        ## Time scaling
        nodes[,2] <- scaling(nodes[,2])
        
        ## Add branch lengths to tree
        nodes <- annotate.tree(nodes)
        ## Add mutations to tree
        for(nn in 1:(2*sample-1)){
            nodes[nn,5] <- rpois(1, theta*nodes[nn,6]/2)
        }

	n.mtns<-sum(nodes[,5]);
	mut.list<-runif(n.mtns);
	seqs<-matrix(0, nrow=sample, ncol=n.mtns);
	cum.mut<-1;
	for (node in 1:(nrow(nodes)-1)) if (nodes[node,5]>0) {
		who.list<-find.descendants(node, nodes, c());
		for (mut in 1:nodes[node,5]) {
			seqs[who.list,cum.mut]=1;
			cum.mut<-cum.mut+1;
		}
	}
	seqs<-seqs[,order(mut.list)];
	mut.list<-mut.list[order(mut.list)];
	

	if (plot.tree == TRUE) plot.tree(nodes, add.mutations=T);

	if (return==TRUE) return(list(tree=nodes, 
		data=list(type="haplotype", n.seq=sample, l.seq=n.mtns, pos=mut.list, names=c(1:sample), seq=seqs)
		
		));

}

make.spectrum <- function(ntrees, nsamples, fun, double.prob=0.1){
for(i in 1:ntrees){
    t<-simulate.coalescent.norec(nsamples, plot.tree=FALSE, verbose=FALSE, scaling=fun)
    freq <- colSums((t$data$seq))
    for(f in freq){spectrum.1[f] <- spectrum.1[f]+1}

    ## Signature 2 mutations
    t<-simulate.coalescent.norec(nsamples, plot.tree=FALSE, verbose=FALSE, scaling=fun, theta=10*(1+double.prob))
    seq <- t$data$seq
    i=1
    while(i < NCOL(seq)){
        if(runif(1)<double.prob){
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
return(list("s1"=spectrum.1, "s2"=spectrum.2))
}


ntrees <- 200000
nsamples <- 50

spectrum.1 <- rep(0, nsamples)
spectrum.2 <- rep(0, nsamples)

identityfun <- approxfun(c(0,100),c(0,100), rule=1:2)
## bottleneckfun<-approxfun(c(0,0.01,0.02,100),c(0,0.01,0.0101,99.9901), rule=1:2)
## expfun<-approxfun(c(0,0.01,100),c(0, 0.1, 100.09), rule=1:2)
## expbottleneckfun<-approxfun(c(0, 0.01, 0.02, 100),c(0, 0.11, 0.1101,100.0901), rule=1:2)
expfun.10 <- function(x){
    t <- 0.01
    end.size <- 10
    g <- log(end.size)/t
    ## return(ifelse(x>t, log(g*t+1)/g+(x-t) ,log(g*x+1)/g))
    return(ifelse(x>t, (exp(g*t)-1)/g+(x-t) ,(exp(g*x)-1)/g))
}

expfun.100 <- function(x){
    t <- 0.01
    end.size <- 100
    g <- log(end.size)/t
    ## return(ifelse(x>t, log(g*t+1)/g+(x-t) ,log(g*x+1)/g))
    return(ifelse(x>t, (exp(g*t)-1)/g+(x-t) ,(exp(g*x)-1)/g))
}
    

expfun.1000 <- function(x){
    t <- 0.01
    end.size <- 1000
    g <- log(end.size)/t
    ## return(ifelse(x>t, log(g*t+1)/g+(x-t) ,log(g*x+1)/g))
    return(ifelse(x>t, (exp(g*t)-1)/g+(x-t) ,(exp(g*x)-1)/g))
}
    


s <- make.spectrum(ntrees, nsamples, identityfun, double.prob=1)
e <- make.spectrum(ntrees, nsamples, expfun.100, double.prob=1)

expect <- (1/(1:nsamples))/sum(1/(1:nsamples))
s.spec1 <- s$s1/sum(s$s1)
s.spec2 <- s$s2/sum(s$s2)
e.spec1 <- e$s1/sum(e$s1)
e.spec2 <- e$s2/sum(e$s2)


pdf("~/spectrum/plots/simulated_repeat_mutation_spectrum.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(s.spec2/expect, pch=16, type="o", col="#E41A1C", xlim=c(0,10), ylim=c(0,3), xlab="Allele count", ylab="Relative frequency spectrum (compared to ~1/n)",cex.axis=1.4, cex.lab=1.5)
points(e.spec1/expect, pch=16, type="o", col="#377EB8", xlim=c(1,10))
points(e.spec2/expect, pch=16, type="o", col="#4DAF4A", xlim=c(1,10))
legend("topright", c("Repeat mutations with constant size", "Single mutations with exponential growth", "Repeat mutations with exponential growth"), col=c("#E41A1C","#377EB8","#4DAF4A"), pch=16, lty=1, bty="n", cex=1.2)
abline(h=1, lty=2)
dev.off()

sc <- make.spectrum(ntrees, nsamples, identityfun, double.prob=0.1)
ec.10 <- make.spectrum(ntrees, nsamples, expfun.10, double.prob=0.1)
ec.100 <- make.spectrum(ntrees, nsamples, expfun.100, double.prob=0.1)
ec.1000 <- make.spectrum(ntrees, nsamples, expfun.1000, double.prob=0.1)

pdf("~/spectrum/plots/simulated_CpG_spectrum_proportion.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
a <- 0.85
plot((1-a)*sc$s2/((1-a)*sc$s2+a*sc$s1), xlim=c(0,20), ylim=c(0.13,0.23), type="o", pch=16, xlab="Allele count", ylab="Proportion of CpG mutations", lwd=2,cex.axis=1.4, cex.lab=1.5)
points((1-a)*ec.10$s2/((1-a)*ec.10$s2+a*ec.10$s1), xlim=c(0,20), pch=16, type="o", col="#66C2A5", lwd=2)
points((1-a)*ec.100$s2/((1-a)*ec.100$s2+a*ec.100$s1), xlim=c(0,20), pch=16, type="o", col="#FC8D62", lwd=2)
points((1-a)*ec.1000$s2/((1-a)*ec.1000$s2+a*ec.1000$s1), xlim=c(0,20), pch=16, type="o", col="#8DA0CB", lwd=2)
legend("topright", c("Constant size", "10-fold exponential growth", "100-fold exponential growth", "1000-fold exponential growth"), pch=16, lty=1, col=c("black", "#66C2A5", "#FC8D62", "#8DA0CB"), bty="n", lwd=2, cex=1.2)
dev.off()
    
save.image(file="~/spectrum/plots/repeat_mutation_data.Rdata")
