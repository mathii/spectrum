## Test whether factor loadings are correlated with principal components.

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")

dataname <- paste0("~/spectrum/plots/","Components_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

PCA <- read.table("~/spectrum/PCA/cteam_europe.evec", as.is=TRUE, header=FALSE)

data.we <- data[reg[rownames(data)]=="WestEurasia",]

lat<-info$Latitude
names(lat)<-gsub("-", ".", info$ID)
lat.we <- lat[rownames(data.we)]]
lon<-info$Longitude
names(lon)<-gsub("-", ".", info$ID)
lon.we <- lon[rownames(data.we)]]
