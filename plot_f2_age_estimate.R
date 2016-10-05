## Plot the age estimates as a function of loadings

source("~/spectrum/code/spectrumlib.R")

####################################################

## sig <- "1"
## regions.to.include <- c("SouthAsia", "WestEurasia")

sig.name <- "2"
regions.to.include <- c("America", "EastAsia")
## regions.to.include <- c("EastAsia","WestEurasia","America","Oceania","SouthAsia","CentralAsiaSiberia")
sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )
####################################################

sig.name.map <- c(4,1)
names(sig.name.map) <- c("1", "2")

####################################################

all <- new.env()
load("~/f2_age/cteam/all/all_results.RData", envir=all)

all.ids <- unique(names(c(all$ID1.pop, all$ID2.pop)))

all$ID1 <- names(all$ID1.pop)
all$ID2 <- names(all$ID2.pop)
names(reg) <- info$ID
all.med <- rep(NA, length(all.ids))
names(all.med) <- all.ids
for(i in 1:length(all.ids)){
  all.med[i] <- median(all$t.hats[all$ID1==all.ids[i]|all$ID2==all.ids[i]])
}

sig.env <- new.env()
load(paste0("~/f2_age/cteam_sig", sig.name, "/all/all_results.RData"), envir=sig.env)

sig.env$ID1 <- names(sig.env$ID1.pop)
sig.env$ID2 <- names(sig.env$ID2.pop)
names(reg) <- info$ID
sig.med <- rep(NA, length(all.ids))
names(sig.med) <- all.ids
for(i in 1:length(all.ids)){
  sig.med[i] <- median(sig.env$t.hats[sig.env$ID1==all.ids[i]|sig.env$ID2==all.ids[i]])
}


## Quick hack for missing samples
all.med <- all.med[names(sig.med)]

components<-read.table("~/spectrum/plots/Components_ica_NMF.n2.r4.txt", as.is=TRUE)
rownames(components) <- gsub(".", "-", components[,1], fixed=TRUE)
components <- components[2:NCOL(components)]
comp <- components[names(sig.med),sig.name.map[sig]]

inc <- reg[names(sig.med)] %in% regions.to.include

plot(comp[inc], sig.med[inc], col=cols[reg[names(sig.med[inc])]], xlab=paste0("Signature ", sig.name, " loading"), ylab="Median f2 age (generations)")

spec="spectrum"
n <- 2
tag <- ""
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")
freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
freq <- t(t(freq2)/colSums(freq2))
ff <- colSums(freq[c(sig),])
names(ff) <- gsub(".", "-", names(ff), fixed=TRUE)
ff <- ff[names(sig.med)]

pdf(paste0("~/spectrum/plots/f2_age_sig", sig.name, ".pdf"))
plot(ff[inc], sig.med[inc], col=cols[reg[names(sig.med[inc])]], xlab=paste0("Proportion of signature ", sig.name, " f2 mutations"), ylab="Median age of signature 2 mutations (generations)", pch=16)
legend("bottomright", c("America", "East Asia"), col=c("#984EA3" , "#E41A1C" ), pch=16, bty="n") 
dev.off()
