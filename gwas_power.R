## Estimate sample size for gwas. 

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"

r1 <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
N.mutations <- 60                       #Per trio, say

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

we <- reg[colnames(norm2)]=="WestEurasia"

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )

alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]
rownames(freq2)<-amap[rownames(freq2)]
freq2<-freq2[amap,]
freq2 <- as.matrix(freq2)


norm2<-t(t(freq2)/colSums(freq2))

norm2 <- norm2[,!is.na(we)]
we <- reg[colnames(norm2)]=="WestEurasia"

## mutations in class
we.val <- mean(colSums(norm2[r1,we]))
nwe.val <- mean(colSums(norm2[r1,!we]))

## Proportional increase in the mutations in class
s.inc <- (we.val+we.val/nwe.val-1)/(1-we.val)

x0 <- N.mutations*nwe.val
x1 <- N.mutations*nwe.val*(1+s.inc)

## Need to get this many sds
sdev <- abs(qnorm(log(0.5e-7), log=TRUE))

Nfunc <- function(p){sdev*sdev*(x1/p/p+x0/(1-p)/(1-p))/(x1-x0)/(x1-x0)}
