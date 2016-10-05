## Test whether there is a correlation between the genetic/physical distance
## matrices and the signature loadings. 
library(geosphere)
source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"
what <- "ica_NMF"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")

dataname <- paste0("~/spectrum/plots/","Components_", what,  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

rownames(data) <- gsub("-", ".", rownames(data))
rownames(info) <- gsub("-", ".", info$ID)
lat <- info[rownames(data), "Latitude"]
lon <- info[rownames(data), "Longitude"]
geo.dist <- distm(cbind(lon, lat), cbind(lon, lat))/1000

gen.dist.tmp <- read.table("/Users/mathii/mapwarp/sgdp/dist.mat", as.is=TRUE)
gen.dist <- matrix(0, NROW(gen.dist.tmp), NCOL(gen.dist.tmp))
for(i in 1:NCOL(gen.dist.tmp)){
    gen.dist[,i] <- gen.dist.tmp[,i]
}

ids <- read.table("/Users/mathii/mapwarp/sgdp/dist.mat.id", as.is=TRUE)
rownames(gen.dist) <- colnames(gen.dist) <- gsub("-", ".", ids[,2])

gen.dist <- gen.dist[rownames(data),rownames(data)]

results <- data.frame("geo.dist"=rep(NA,rank), "gen.dist"=rep(NA,rank))

for(sig in 1:rank){
    dm <-  as.matrix(dist(data[,sig]))
    results[sig,] <- c(mantel.test(dm, geo.dist)$p, mantel.test(dm, gen.dist)$p)
}
