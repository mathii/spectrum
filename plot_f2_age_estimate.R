## Plot the age estimates as a function of loadings

source("~/spectrum/code/spectrumlib.R")

####################################################

## sig.name <- "1"
## regions.to.include <- c( "WestEurasia", "EastAsia")
## regions.to.include.all <- c()
## sig <- c("TCC.T", "ACC.T", "TCT.T", "CCC.T")
## sig.rc <- c("GGA.A", "GGT.A", "AGA.A", "GGG.A" )

sig.name <- "2"
regions.to.include <- c("America", "EastAsia")
regions.to.include.all <- c("America", "EastAsia")

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
comp <- components[names(sig.med),sig.name.map[sig.name]]

inc <- reg[names(sig.med)] %in% regions.to.include

## plot(comp[inc], sig.med[inc], col=cols[reg[names(sig.med[inc])]], xlab=paste0("Signature ", sig.name, " loading"), ylab="Median f2 age (generations)")

spec="spectrum"
n <- 2
tag <- ""
inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n,tag, ".txt")
freq2 <- read.table(inname, header=TRUE, as.is=TRUE )
freq <- t(t(freq2)/colSums(freq2))
ff <- colSums(freq[c(sig),])
names(ff) <- gsub(".", "-", names(ff), fixed=TRUE)
ff <- ff[names(sig.med)]

lab <- c()
xlim <- c(0.14, 0.24)
ylim <- c(0,500)
if(sig.name=="2"){
    laby <- c(4,5,76)
    xlim <- c(0.14, 0.24)
    ylim <- c(0,500)
}
if(sig.name=="1"){
    laby <- c()    
    xlim <- c(0.06, 0.12)
    ylim <- c(0,500)
}

pdf(paste0("~/spectrum/plots/f2_age_sig", sig.name, ".pdf"))
plot(ff[inc], sig.med[inc], col=cols[reg[names(sig.med[inc])]], xlab=bquote("Proportion of signature"~.(sig.name)~f[2]~"mutations"), ylab=bquote("Median age of signature"~.(sig.name)~f[2]~"mutations (generations)"), pch=16, xlim=xlim, ylim=ylim, bty="n")
if(length(regions.to.include.all)){
    inc2 <- reg[names(sig.med)] %in% regions.to.include.all
    points(ff[inc2], all.med[inc2], col=cols[reg[names(sig.med[inc2])]], pch=1)
}
legend("bottomright", c(regions.to.include, paste0(regions.to.include.all, " (all f2 variants)")), col=c(cols[regions.to.include],cols[regions.to.include.all]), pch=c(rep(16, length(regions.to.include)), rep(1, length(regions.to.include.all))), bty="n") 
text(ff[inc][laby], sig.med[inc][laby], names(sig.med[inc])[laby], pos=4, cex=0.5)
dev.off()
