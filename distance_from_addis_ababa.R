## Test signatures a function of heterozygosity (correcting for distance from Africa...)
source("~/spectrum/code/spectrumlib.R")
library(geosphere)

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
which <- 4                              #Which spectrum - 4=="Signature 1" 1=="Signature 2"
spec <- "spectrum"
what <- "ica_NMF"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
sig.map=c(2,4,3,1)

dataname <- paste0("~/spectrum/plots/","Components_", what,  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

het.rate <- read.table("~/spectrum/het_rate/het_rate_mlrho.txt", as.is=TRUE, header=TRUE)
hn <- gsub("-", ".", het.rate[,1])
het.rate <- het.rate[,2]
names(het.rate) <- hn
het.rate <- het.rate[rownames(data)]

rownames(data) <- gsub("-", ".", rownames(data))
rownames(info) <- gsub("-", ".", info$ID)
info <- info[rownames(data),]

lat <- info[rownames(data), "Latitude"]
lon <- info[rownames(data), "Longitude"]

include <- !(info[rownames(data),"Region"]=="Africa")
het.lm <- lm(data[include,which]~het.rate[include])

pdf(paste0("~/spectrum/hets/het_vs_sig.c", which , ".pdf"))
plot(het.rate[include], data[include,which], col=cols[info[include,"Region"]], pch=21, xlab="Heterozygosity", ylab=paste0("Signature ", sig.map[which], " loading"))
legend("topright", names(cols), col=cols, bty="n", pch=16)
dev.off()

## Compute distances as in Rmachandran (without Istanbul Waypoint though).
distance.from.cairo <-  c(distm(cbind(lon,lat), c(31.23, 30.04)))/1000
names(distance.from.cairo) <- rownames(data)
oceania <- info$Region=="Oceania"
distance.from.cairo[oceania] <- 7867.965+c(distm(cbind(lon,lat)[oceania,], c(104,11)))/1000
america <- info$Region=="America"
distance.from.cairo[america] <- 12173.961+c(distm(cbind(lon,lat)[america,], c(-130,54)))/1000

pdf(paste0("~/spectrum/hets/het_vs_dist.c", which , ".pdf"))
plot(distance.from.cairo[include], data[include,which], col=cols[info[include,"Region"]], pch=21, xlab="Distance from Cairo (km)", ylab=paste0("Signature ", sig.map[which], " loading"))
legend("topleft", names(cols), col=cols, bty="n", pch=16)
dev.off()


mod <- lm(data[include,which]~distance.from.cairo[include]+het.rate[include])
print(summary(mod))

## Also plot distance from addis
distance.from.addis.ababa <- distance.from.cairo
africa <- info$Region=="Africa"
distance.from.addis.ababa[africa] <- c(distm(cbind(lon,lat)[africa,], c(38,9)))/1000
distance.from.addis.ababa[!africa] <- distance.from.addis.ababa[!africa]+2608.397
labl <-  c(64, 139, 158, 183, 20,  47, 240)
pdf("~/spectrum/hets/distance_from_addis_ababa.pdf")
plot(distance.from.addis.ababa, het.rate, col=cols[info$Region], xlab="Distance from Addis Ababa (km)", ylab="Heterozygosity")
text(distance.from.addis.ababa[labl], het.rate[labl], rownames(data)[labl], cex=0.5, pos=4)
legend("topright", names(cols), col=cols, bty="n", pch=16)
dev.off()

## Now plot conditioned on San
conditioned<-read.table("het_conditioned_S_Khomani_San-1.txt", as.is=T)
conditioned.het <- conditioned[,3]/conditioned[,2]
names(conditioned.het) <- gsub("-", ".", conditioned[,1])
conditioned.het<-conditioned.het*het.rate["S_Mbuti.1"]/conditioned.het["S_Mbuti.1"]
conditioned.het <- conditioned.het[!grepl("Khomani_San", names(conditioned.het))]
conditioned.het <- conditioned.het[!grepl("Ju_hoan", names(conditioned.het))]
losch.dist <- c(distm( c(6.4, 49.81), c(31.23, 30.04))/1000+2608.397)
stut.dist <- c(distm( c(9.18, 48.78), c(31.23, 30.04))/1000+2608.397)
ust.dist <- c(distm( c(71.1, 57.7), c(31.23, 30.04))/1000+2608.397)
d.ex <- c(distance.from.addis.ababa, Loschbour=losch.dist, Stuttgart=stut.dist, Ust_Ishim=ust.dist)
d.ex <- d.ex[names(conditioned.het)]
nm <- names(conditioned.het)
dcol <- ifelse(nm %in% rownames(info), cols[info[nm,"Region"]], "black")
dpch <- ifelse(nm %in% c("Loschbour", "Stuttgart", "Ust_Ishim"), 16, 1)
pdf("~/spectrum/hets/distance_from_addis_ababa_conditional.pdf")
plot(d.ex, conditioned.het, col=dcol, xlab="Distance from Addis Ababa (km)", ylab="Heterozygosity", pch=dpch)
tt <- c("Loschbour", "Stuttgart", "Ust_Ishim")
text( d.ex[tt], conditioned.het[tt], tt, pos=c(1,1,4))
labl <-  c(64, 139, 158, 183, 20,  47)
text(distance.from.addis.ababa[labl], het.rate[labl], rownames(data)[labl], cex=0.5, pos=4)
dev.off()
