## Estimate the increase in the total mutation rate.
set.seed(12345)

library("RColorBrewer")
library("NMF")
library("fastICA")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
spec <- "spectrum"

s1 <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
s2 <- c("ACG.T","CCG.T", "GCG.T", "TCG.T")

tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
prop2<-t(t(freq2)/colSums(freq2))

ps1 <- colSums(prop2[s1,])
ps2 <- colSums(prop2[s2,])

n.afr <- names(which(reg=="Africa" & names(reg) %in% colnames(prop2)))


incfunc <- function(k1, k2){return((k2-k1)/(1-k2))}
afr1 <- ps1[n.afr]
we1 <- ps1[names(which(reg=="WestEurasia"& names(reg) %in% colnames(prop2)))]



afr2 <- ps2[n.afr]
m2 <- mean(afr2)
amr2 <-ps2[c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")]
