## Compare the decompositions from PCA, sparse PCA and Non-negative Matrix factorisation.

exclude.cell.lines <- FALSE
n <- 2
inname <- paste0("~/spectrum/data/spectrum_matrix.n", n, ifelse(exclude.cell.lines, ".NoCellLines", ""), ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
freq2 <- freq2[,!( colnames(freq2) %in% ancient)]

## Boilerplate - load info
info <- read.table("~/spectrum/data/location_info.txt", as.is=TRUE, header=TRUE, sep="\t")
info[info[,6]=="Genomic from saliva",6]<-"Genomic_from_saliva"
info[info[,6]=="?",6]<-"Unknown"
regions <- unique(info[,3])
cols <- brewer.pal(length(regions), "Set1")
cols[6] <- "darkgrey"
sources <- unique(info[,6])
source.cols <- brewer.pal(length(sources), "Set2")
names(cols) <- regions
names(source.cols) <- sources
name.map <- info[,3]
names(name.map) <- gsub("-", ".", info$ID, fixed=TRUE)
src <- info[,6]
names(src) <- gsub("-", ".", info$ID, fixed=TRUE)

plot.components <- function(components, name.map, src, cols="black", n.components=4){
    par(mfrow=c(floor(sqrt(n.components)),ceiling(n.components/floor(sqrt(n.components)))))
    for(i in 1:n.components){
        plot(components[,i], components[,i+1], col=cols[name.map[rownames(components)]],, xlab=paste0("PCA", i), ylab=paste0("PCA", i+1), pch=ifelse(src[rownames(components)]=="Genomic_from_cell_lines", 13, 1))
        legend("topleft", c("Cell lines", "Other"), pch=c(13,1), bty="n")
        legend("topright", names(cols), col=cols, bty="n", pch=16)
    }
}

plot.loadings <- function(loadings, n.loadings=4){
    
}

## PCA
pdf(paste0("~/spectrum/plots/","Spectrum_PCA", tag, ".n", n, ".pdf"), 12, 12)
pca=prcomp(freq2, center=TRUE, scale.=TRUE)
plot.components(pca$rotation, name.map, src, cols=cols)
dev.off()
