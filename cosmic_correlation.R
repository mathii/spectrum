## Compute the correlation between cosmic and our detected signatures.

exclude.cell.lines <- FALSE
n <- 2
rank <- 4
spec <- "spectrum"
tag <- ifelse(exclude.cell.lines, ".NoCellLines", "")

dataname <- paste0("~/spectrum/plots/","Loadings_ica_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, ".txt")
data <- read.table(dataname, as.is=TRUE)
rownames(data) <- data[,1]
data <- data[,2:NCOL(data)]

cosmic <- read.table("~/spectrum/cosmic/signatures_probabilities.txt", as.is=TRUE, header=TRUE, sep="\t")
ms <- cosmic[,paste0("Signature.", 1:30)]
rownames(ms) <- paste0(cosmic$Trinucleotide,".", substr(cosmic$Substitution.Type,3,3))
ms <- as.matrix(ms)

alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
ms <- ms[alex[,2],]

cc<-cor(ms, data)

for(i in 1:rank){
    cat(paste0("\nComponent ", i, "\n"))
    print(head(cc[order(-cc[,i]),i,drop=FALSE]))
}
