## Plot the components on a map.
source("~/spectrum/code/spectrumlib.R")
library("RColorBrewer")
library(mgcv)
library(maps)

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")

dataname <- paste0("~/spectrum/plots/","Components_ica_NMF",  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

rownames(info) <- gsub("-", ".", info$ID)
lat <- info[rownames(data), "Latitude"]
lon <- info[rownames(data), "Longitude"]

nxpix <- 361*2
nypix <- 301*2

for(i in 1:NCOL(data)){
    ## spln<-gam(data[,i]~s(lat, lon))
    ## grid <- data.frame(lat=rep(seq(-150,150, length.out=nypix), each=nxpix), lon=rep(seq(-180,180, length.out=nxpix), times=nypix))
    ## res<-matrix(predict(spln, grid), nrow=nxpix, ncol=nypix)

    pal <- brewer.pal(9, "YlOrRd")
    ## breaks=seq(min(min(data[,i]),min(res))*1.01, max(max(data[,i]), max(res))*1.01,length.out=10)
    breaks=seq(min(data[,i])*0.99, max(data[,i])*1.01,length.out=10)

    find.interval <- function(x, breaks){max(which(breaks<x))}

    oname <-  paste0("~/spectrum/plots/Map_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank,".c", i, tag, ".png")
    png(oname, width=12, height=6, res=200, units="in")
    ## image(x=seq(-180,180, length.out=nxpix), y=seq(-150,150, length.out=nypix), res, ylim=c(-60,80), xlab="Longitude", ylab="Lattitude", col=pal, breaks=breaks)
    ## maps::map(database="world", fill=FALSE, add=TRUE, interior=FALSE, ylim=c(-60,80))
    maps::map(database="world", fill=FALSE, add=FALSE, interior=FALSE, ylim=c(-60,80))
    points(jitter(lon,100), jitter(lat,100), pch=21, bg=paste0(pal[sapply(data[,i], find.interval, breaks=breaks)]), cex=1.5)
    dev.off()
}

## png(paste0("~/cteam/spectrum/Negative_map.png"), width=12, height=6, res=200, units="in")
## image(x=seq(-180,180, length.out=nxpix), y=seq(-150,150, length.out=nypix), 0*res, ylim=c(-60,80), xlab="Longitude", ylab="Lattitude", col=c("white"), breaks=c(-1,1))
## maps::map(database="world", fill=TRUE, col="black", add=TRUE, interior=FALSE, bg="white", ylim=c(-60,80))
## dev.off()
