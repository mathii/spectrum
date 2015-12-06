## Test spectrum. 
library("RColorBrewer")

exclude.cell.lines <- FALSE
n <- 2
tag <- ""

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


ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
exclude <- c()
ref <- c("HumanRef",  "panTro2",  "Ancestor")

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
names(name.map) <- info$ID
src <- info[,6]
names(src) <- info$ID

raw <- read.table(paste0("~/spectrum/data/counts/chr1.n", n ,".txt"), as.is=TRUE, header=TRUE)
names(raw)<-gsub(".", "-", names(raw), fixed=TRUE)
data <- data.matrix(raw[,2:NCOL(raw)])
rownames(data) <- raw[,1]
data <- data[order(rownames(data)),]
for(chr in 2:22){
    raw <- read.table(paste0("~/spectrum/data/counts/chr", chr, ".n", n, ".txt"), as.is=TRUE, header=TRUE)
    d <- data.matrix(raw[,2:NCOL(raw)])
    rownames(d) <- raw[,1]
    d <- d[order(rownames(d)),]
    data <- data+d
}

totals <- colSums(data)
min.total <- 1000
data <- data[,totals>min.total]
## data <- data[,!(colnames(data) %in% ancient)]
data <- data[,!(colnames(data) %in% exclude)]
data <- data[,!(colnames(data) %in% ref)]
exclude.ABteam <- grep("[AB]\\.", colnames(data), value=TRUE)
data <- data[,!(colnames(data) %in% exclude.ABteam)]

## Exclude cell lines
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
## freqs2 <- t(t(data2)/colSums(data2))
## if(!(all(colSums(freqs2)==1))){stop("Columns do not sum to 1")}

freqs2 <- data2/c(data2["TAT.G",])

outname <- paste0("~/spectrum/data/spectrum_matrix.n", n, ifelse(exclude.cell.lines, ".NoCellLines", ""), ".txt")
write.table(freqs2, outname, col.names=TRUE, row.names=TRUE, quote=F)

outname <- paste0("~/spectrum/data/count_matrix.n", n, ifelse(exclude.cell.lines, ".NoCellLines", ""), ".txt")
write.table(data2, outname, col.names=TRUE, row.names=TRUE, quote=F)
