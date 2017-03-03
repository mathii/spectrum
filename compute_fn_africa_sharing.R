## Calculation for a comment in the paper.
source("~/spectrum/code/spectrumlib.R")


n <- 2
sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )

pA <- paste0("~/spectrum/data/count_matrix.n", n, ".poly_Africa.txt")
npA <- paste0("~/spectrum/data/count_matrix.n", n, ".notpoly_Africa.txt")

pA <- read.table(pA, as.is=TRUE)
npA <- read.table(npA, as.is=TRUE)

pA.t <- colSums(pA[,])
npA.t <- colSums(npA[,])

pA.2 <- colSums(pA[sig,])
npA.2 <- colSums(npA[sig,])

pA.n2 <- colSums(pA[!(rownames(pA)%in%sig),])
npA.n2 <- colSums(npA[!(rownames(npA)%in%sig),])


pAstat <- sum(pA.2[reg[names(pA.2)]!="Africa"])/sum(pA.t[reg[names(pA.t)]!="Africa"])
npAstat <- sum(npA.2[reg[names(npA.2)]!="Africa"])/sum(npA.t[reg[names(npA.t)]!="Africa"])

c(pAstat, npAstat)

fstat <- sum(pA.2[reg[names(pA.2)]!="Africa"])/(sum(pA.2[reg[names(pA.2)]!="Africa"])+sum(npA.2[reg[names(npA.2)]!="Africa"])/2)

nfstat <- sum(pA.n2[reg[names(pA.n2)]!="Africa"])/(sum(pA.n2[reg[names(pA.n2)]!="Africa"])+sum(npA.n2[reg[names(npA.n2)]!="Africa"])/2)

c(fstat, nfstat, fstat-nfstat)
