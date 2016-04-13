## Test spectrum. 
library("RColorBrewer")

exclude.cell.lines <- FALSE
n <- 2
tag <- wut <- ""

cA <- commandArgs(TRUE)
if(length(cA)>0){
    n <- as.numeric(cA[1])
}
if(length(cA)>1){
    exclude.cell.lines <- as.logical(as.numeric(cA[2]))
}
if(length(cA)>2){
    wut <- paste0(".", cA[3])
}


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


exclude <- c()
ref <- c("HumanRef",  "panTro2",  "Ancestor")

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
names(name.map) <- info$ID
src <- info[,6]
names(src) <- info$ID

raw <- read.table(paste0("~/spectrum/counts/chr1.n", n , wut, ".txt"), as.is=TRUE, header=TRUE)
names(raw)<-gsub(".", "-", names(raw), fixed=TRUE)
data <- data.matrix(raw[,2:NCOL(raw)])
rownames(data) <- raw[,1]
data <- data[order(rownames(data)),]
for(chr in 2:22){
    raw <- read.table(paste0("~/spectrum/counts/chr", chr, ".n", n, wut, ".txt"), as.is=TRUE, header=TRUE)
    d <- data.matrix(raw[,2:NCOL(raw)])
    rownames(d) <- raw[,1]
    d <- d[order(rownames(d)),]
    data <- data+d
}

totals <- colSums(data)
min.total <- 1000
data <- data[,totals>min.total]
data <- data[,!(colnames(data) %in% exclude)]
data <- data[,!(colnames(data) %in% ref)]

if(exclude.cell.lines){
    data <- data[,(colnames(data) %in% names(src))]
    data <- data[,src[colnames(data)]!="Genomic_from_cell_lines"]
    tag <- "exCellLines"
}

primary.mut <- rownames(data)
for(i in 1:length(primary.mut)){
    primary.mut[i] <- sort(c(primary.mut[i], reverse.complement(primary.mut[i])))[2]
}

data2 <- aggregate(data, by=list(primary=primary.mut), sum)
prim.order <- data2$primary
data2 <- data2[,2:NCOL(data2)]
rownames(data2) <- prim.order

## Rename and reorder to alexandrov format
alex <- read.table("~/spectrum/code/alexandrovmap.txt", as.is=TRUE, header=FALSE)
amap <- alex[,2]
names(amap) <- alex[,1]
rownames(data2)<-amap[rownames(data2)]
data2<-data2[amap,]
freqs2 <- data2/c(data2["ATA.C",])

freqs2 <- as.matrix(freqs2)
data2 <- as.matrix(data2)

outname <- paste0("~/spectrum/data/spectrum_matrix.n", n, wut, ifelse(exclude.cell.lines, ".NoCellLines", ""), ".txt")
write.table(freqs2, outname, col.names=TRUE, row.names=TRUE, quote=F)

outname <- paste0("~/spectrum/data/count_matrix.n", n, wut, ifelse(exclude.cell.lines, ".NoCellLines", ""), ".txt")
write.table(data2, outname, col.names=TRUE, row.names=TRUE, quote=F)
