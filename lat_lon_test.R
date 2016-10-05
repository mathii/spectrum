source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
which <- 4                              #Which spectrum - 4=="Signature 1"
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

include <- info[rownames(data), "Region"] %in% c("WestEurasia")
lat.lm <- lm(data[include,which]~lat[include])
lon.lm <- lm(data[include,which]~lon[include])
