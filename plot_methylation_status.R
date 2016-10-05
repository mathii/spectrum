## Compare methylation results for signature 2.

source("~/spectrum/code/spectrumlib.R")

sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )
sig.name <- 2
in.pops <- "America"
out.pops <- "America"
in.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")
n <- 2

##all
all <- read.table(paste0("~/spectrum/data/count_matrix.n",n,".txt"), as.is=TRUE, header=TRUE)

##high methylation
plus <- read.table(paste0("~/spectrum/methyl/chr1.n",n,".methyl_gt50.txt"), as.is=TRUE, header=TRUE)
row.names(plus) <- plus[,1]
plus <- plus[,3:NCOL(plus)]
for(chr in 2:22){
pp <- read.table(paste0("~/spectrum/methyl/chr", chr, ".n",n,".methyl_gt50.txt"), as.is=TRUE, header=TRUE)
row.names(pp) <- pp[,1]
plus <- plus + pp[,3:NCOL(pp)]
}

##low methylation
minus <- read.table(paste0("~/spectrum/methyl/chr1.n",n,".methyl_lt50.txt"), as.is=TRUE, header=TRUE)
row.names(minus) <- minus[,1]
minus <- minus[,3:NCOL(minus)]
for(chr in 2:22){
mm <- read.table(paste0("~/spectrum/methyl/chr", chr, ".n",n,".methyl_lt50.txt"), as.is=TRUE, header=TRUE)
row.names(mm) <- mm[,1]
minus <- minus + mm[,3:NCOL(mm)]
}

all.count <- aggregate(colSums(all[sig,]), FUN=sum, by=list(name.map[colnames(plus)]))
hi.count <- aggregate(colSums(plus[sig,]), FUN=sum, by=list(name.map[colnames(plus)]))
lo.count <-  aggregate(colSums(minus[sig,]), FUN=sum, by=list(name.map[colnames(plus)]))

data_counts=cbind(hi.count$x, lo.count$x)
gm<-glm(data_counts ~ as.factor(hi.count[,1]), family="binomial")
print(summary(gm))

tab <- matrix( c( sum(plus[c(sig,sig.rc),!(colnames(plus)%in%in.ind)]), sum(minus[c(sig,sig.rc),!(colnames(minus)%in%in.ind)]), sum(plus[c(sig,sig.rc),in.ind]), sum(minus[c(sig,sig.rc),in.ind])),2,2)
print(fisher.test(t(tab)))
