## Test for enrichment in the CCTT spectrum.

library("RColorBrewer")

source("~/spectrum/code/spectrumlib.R")
f=2

inname <- paste0("~/spectrum/data/count_matrix.n",f,".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
## Exclude ancient samples for this analysis
ancient <- c("AltaiNeandertal", "Denisova", "Loschbour", "LBK1b_leipzig_v2", "Ust_Ishim")
freq2 <- freq2[,!( colnames(freq2) %in% ancient)]


norm <- freq2["TAT.G",]

CCTT.tmp <- read.table("~/spectrum/CCTT/chr1.n1.txt", as.is=TRUE)
CCTT <- rep(0, NROW(CCTT.tmp))
names(CCTT) <- CCTT.tmp[,1]

for(chr in 1:22){
    CCTT.tmp <- read.table(paste0("~/spectrum/CCTT/chr",chr,".n",f ,".txt"), as.is=TRUE)
    CCTT <- CCTT+CCTT.tmp[,2]
}

freq <- unlist(CCTT[gsub(".", "-", names(norm), fixed=T)]/norm)

rownames(info) <- info$ID
reg <- info[gsub(".", "-", names(norm), fixed=T),"Region"]
src <- info[gsub(".", "-", names(norm), fixed=T), "Source"]

mod <- lm(freq~reg+src)
print(summary(mod))
