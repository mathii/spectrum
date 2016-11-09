## Plot the proportion of variants in each class, ALL vs f2

source("~/spectrum/code/spectrumlib.R")

####################################################

what <- ""
cA <- commandArgs(TRUE)
if(length(cA)>0){
  what <- paste0(".", cA[1])
}

####################################################

## sig.name <- 1
## sig <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
## ylim=c(0.065, 0.105)

sig.name <- 2
sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
ylim=c(0.13, 0.22)
in.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")

####################################################

sig.name.map <- c(4,1)
sig.map=c(2,4,3,1)

####################################################

regions <- sort(unique(info$Region))
ltys <- rep(1, length(cols))
names(ltys) <- names(cols)
if(sig.name==2){
    regions <- c("Africa", "America_Hi", "America_Lo", "CentralAsiaSiberia", "EastAsia", "Oceania", "SouthAsia", "WestEurasia")
    hi.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")
    ## cols <- c(cols, "America_Hi"="#984EA3", "America_Lo"="#984EA3")
    ## ltys <- rep(1, length(cols))
    ## names(ltys) <- names(cols)
    ## ltys["America_Lo"] <- 2
}

rownames(info) <- gsub("-", ".", info$ID)

## if(sig.name==2){
##     info[hi.ind,"Region"] <- "America_Hi"
##     info[info[rownames(info),"Region"]=="America","Region"] <- "America_Lo"
## }

f2.data <- read.table(paste0("~/spectrum/data/count_matrix.n2.txt"), as.is=TRUE, header=TRUE)
ALL.data <- read.table(paste0("~/spectrum/data/count_matrix.nALL.txt"), as.is=TRUE, header=TRUE)

f2.proportion <- colSums(f2.data[sig,])/colSums(f2.data)
ALL.proportion <- colSums(ALL.data[sig,])/colSums(ALL.data)

pdf(paste0("~/spectrum/plots/plot_all_against_f2_sig", sig.name, ".pdf"))
plot(f2.proportion, ALL.proportion, col=cols[reg[names(f2.proportion)]], ylim=c(0.1190, 0.1210), xlab="Proportion of f2 variants that are signature 2", ylab="Proportion of all variants, per genome, that are signature 2")
legend("bottomright", names(cols), col=cols, bty="n", pch=1)
dev.off()

t.test(ALL.proportion[reg[names(ALL.proportion)]%in%c("America", "CentralAsiaSiberia", "EastAsia", "Oceania", "SouthAsia", "WestEurasia")], ALL.proportion[hi.ind])
