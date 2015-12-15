## Estimate the increase in the total mutation rate.
set.seed(12345)

library("RColorBrewer")
library("NMF")
library("fastICA")

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"


tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
freq2 <- freq2[,!( colnames(freq2) %in% ancient)]

## Rename and reorder to alexandrov format
alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]
rownames(freq2)<-amap[rownames(freq2)]
freq2<-freq2[amap,]
freq2 <- as.matrix(freq2)

## Subtract off minimal mutations. 
freq2 <- freq2-apply(freq2, 1, min)
freq2["ATA.C",] <- 0.0001         #Null row.

## Non-negative Matrix factorization
nnegmf=nmf(as.matrix(freq2),rank=rank, seed="ica", nrun=20)

inname <- paste0("~/spectrum/data/count_matrix.n", n,tag, ".txt")
count2 <- read.table(inname, header=TRUE, as.is=TRUE )
count2 <- count2[,!( colnames(count2) %in% ancient)]
count2 <- t(t(count2)/colSums(count2))

ACAT <- count2["TAT.G",]

co <- coef(nnegmf)
ba <- basis(nnegmf)

plot.loadings(ba, n.loadings=rank)

s1 <- co[2,]*sum(ba[,2])*ACAT
reg <- info$Region
names(reg) <- gsub("-", ".", info$ID)

s1.WE <- median(s1[names(which(reg=="WestEurasia"))], na.rm=TRUE)
s1.NA <- median(s1[names(reg)[(!(reg %in% c("WestEurasia", "SouthAsia")))]], na.rm=TRUE)
dif <- s1.WE-s1.NA

s2 <- co[1,]*sum(ba[,1])*ACAT
afsa <- c("Chane", "Piapoco", "Quechua", "Mayan", "Nahua", "Zapotec", "Mixtex")
inc <- rep(FALSE, length(s2))
for(a in afsa){
    inc <- inc|grepl(a, names(s2))
}

dif2 <- median(s2[inc])-median(s2[!inc])
