library("RColorBrewer")

## Boilerplate, load info
info <- read.table("~/spectrum/code/location_info.txt", as.is=TRUE, header=TRUE, sep="\t")
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
reg <- info[,3]
names(reg) <- gsub("-", ".", info$ID, fixed=TRUE) 

reverse.complement <- function(str){
    mm <- c("A", "C", "G", "T")
    names(mm) <- c("T", "G", "C", "A")
    bits <- strsplit(str, ".", fixed=TRUE)[[1]]
    for(i in 1:length(bits)){
        b=bits[i]
        bb=strsplit(b,"")[[1]]
        cc=mm[bb]
        bits[i]=paste0(rev(cc), collapse="")
    }
    return(paste0(bits, collapse="."))
}

## Plot principal (or other) components
plot.components <- function(components, name.map, src, cols="black", n.components=4, layout=NA, xploti=c(), yploti=c()){
    if(all(is.na(layout))){
        par(mfrow=c(floor(sqrt(n.components-1)),ceiling((n.components-1)/floor(sqrt(n.components-1)))))
    }else{
        par(mfrow=layout)
    }

    if(length(xploti)==0){
        xploti <- c(1:(n.components-1))
        yploti <- c(2:n.components)
    }
    
    for(i in 1:length(xploti)){
        plot(components[,xploti[i]], components[,yploti[i]], col=cols[name.map[rownames(components)]],, xlab=paste0("PCA", xploti[i]), ylab=paste0("PCA", yploti[i]), pch=ifelse(src[rownames(components)]=="Genomic_from_cell_lines", 13, 1))
        legend("topleft", c("Cell lines", "Other"), pch=c(13,1), bty="n")
        legend("topright", names(cols), col=cols, bty="n", pch=16)
    }
}

## Plot factor loadings
plot.loadings <- function(loadings, n.loadings=4){
    cols <- rep(c("#1EBFF0", "#000000", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5"), each=16)
    
    par(mfrow=c(n.loadings,1))
    for(i in 1:n.loadings){
        plot(loadings[,i], col=cols, xaxt="n", bty="n", xlab=paste0("Component ", i ), ylab="Loading", pch=16)
        segments(1:NROW(loadings), 0, 1:NROW(loadings), loadings[,i], col=cols, lwd=2)
        axis(1, at=1:NROW(loadings), labels=names(loadings[,i]), cex.axis=0.8, las=2)
    }
}

