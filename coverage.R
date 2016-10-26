## Plot coverage against signature. 

source("~/spectrum/code/spectrumlib.R")

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
which <- 3                              #Which spectrum - 4=="Signature 1" 1=="Signature 2"
spec <- "spectrum"
what <- "ica_NMF"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")
sig.map=c(2,4,3,1)

dataname <- paste0("~/spectrum/plots/","Components_", what,  ifelse(spec=="spectrum", "", paste0(spec, "_")), ".n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

coverage <- read.table("~/spectrum/hets/coverage.txt", as.is=TRUE, header=TRUE)
hn <- gsub("-", ".", coverage[,1])
coverage <- coverage[,2]
names(coverage) <- hn
coverage <- coverage[rownames(data)]

rownames(data) <- gsub("-", ".", rownames(data))
rownames(info) <- gsub("-", ".", info$ID)
info <- info[rownames(data),]

pdf(paste0("~/spectrum/hets/cov_vs_dist.c", which , ".pdf"))
plot(coverage, data[,which], col=cols[info[,"Region"]], pch=21, xlab="Coverage", ylab=paste0("Signature ", sig.map[which], " loading"))
legend("topright", names(cols), col=cols, bty="n", pch=16)
dev.off()
