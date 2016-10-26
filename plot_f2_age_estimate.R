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
all.med <- all.med.private <- n.all <- n.all.private <- rep(NA, length(all.ids))
names(all.med) <- names(all.med.private) <- names(n.all) <- names(n.all.private) <- all.ids
for(i in 1:length(all.ids)){
    iall <-all$ID1==all.ids[i]|all$ID2==all.ids[i]
    iall.private <- (all$ID1==all.ids[i]|all$ID2==all.ids[i]) &(reg[all$ID1]== reg[all$ID2])
    all.med[i] <- median(all$t.hats[iall])
    all.med.private[i] <- median(all$t.hats[iall.private])
    n.all[i] <- sum(iall*all$haps$f2)
    n.all.private[i] <- sum(iall.private*all$haps$f2)

}


sig.env <- new.env()
load(paste0("~/f2_age/cteam_sig", sig.name, "/all/all_results.RData"), envir=sig.env)

sig.env$ID1 <- names(sig.env$ID1.pop)
sig.env$ID2 <- names(sig.env$ID2.pop)
names(reg) <- info$ID
sig.med <- sig.med.private <- n.sig <- n.sig.private <- rep(NA, length(all.ids))
names(sig.med) <- names(sig.med.private) <- names(n.sig) <- names(n.sig.private) <- all.ids
for(i in 1:length(all.ids)){
    isig <- sig.env$ID1==all.ids[i]|sig.env$ID2==all.ids[i]
    isig.private <- (sig.env$ID1==all.ids[i] | sig.env$ID2==all.ids[i])&( reg[sig.env$ID1]== reg[sig.env$ID2])
    sig.med[i] <- median(sig.env$t.hats[isig])
    sig.med.private[i] <- median(sig.env$t.hats[isig.private])
    n.sig[i] <- sum(isig*sig.env$haps$f2)
    n.sig.private[i] <- sum(isig.private*sig.env$haps$f2)
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

pdf(paste0("~/spectrum/plots/f2_age_private_sig", sig.name, ".pdf"))
plot(ff[inc], sig.med.private[inc], col=cols[reg[names(sig.med[inc])]], xlab=bquote("Proportion of signature"~.(sig.name)~f[2]~"mutations"), ylab=bquote("Median age of signature"~.(sig.name)~f[2]~"mutations (generations)"), pch=16, xlim=xlim, ylim=ylim, bty="n")
if(length(regions.to.include.all)){
    inc2 <- reg[names(sig.med)] %in% regions.to.include.all
    points(ff[inc2], all.med.private[inc2], col=cols[reg[names(sig.med[inc2])]], pch=1)
}
legend("bottomright", c(regions.to.include, paste0(regions.to.include.all, " (all f2 variants)")), col=c(cols[regions.to.include],cols[regions.to.include.all]), pch=c(rep(16, length(regions.to.include)), rep(1, length(regions.to.include.all))), bty="n") 
text(ff[inc][laby], sig.med.private[inc][laby], names(sig.med.private[inc])[laby], pos=4, cex=0.5)
dev.off()

prop.all <- n.sig/n.all
prop.all.private <- n.sig.private/n.all.private



## Count regional sharing matrix!
all.matrix <- matrix(0, nrow=length(regions), ncol=length(regions))
rownames(all.matrix) <- colnames(all.matrix) <- regions
sig.matrix <- matrix(0, nrow=length(regions), ncol=length(regions))
rownames(sig.matrix) <- colnames(sig.matrix) <- regions
for(reg1 in regions){
    for(reg2 in regions){
        total.all <- sum((reg[all$ID1]==reg1&reg[all$ID2]==reg2)|(reg[all$ID2]==reg1&reg[all$ID1]==reg2)*all$haps$f2)
        total.sig <- sum((reg[sig.env$ID1]==reg1&reg[sig.env$ID2]==reg2)|(reg[sig.env$ID2]==reg1&reg[sig.env$ID1]==reg2)*sig.env$haps$f2)
        all.matrix[reg1,reg2] <- all.matrix[reg1,reg2] <- total.all
        sig.matrix[reg1,reg2] <- sig.matrix[reg1,reg2] <- total.sig
    }
}

sig.100.matrix <- sig.matrix*0
all.100.matrix <- all.matrix*0
for(reg1 in regions){
    for(reg2 in regions){
        total.all.100 <- sum((reg[all$ID1]==reg1&reg[all$ID2]==reg2)|(reg[all$ID2]==reg1&reg[all$ID1]==reg2)*(all$t.hats<100)*all$haps$f2)
        total.sig.100 <- sum((reg[sig.env$ID1]==reg1&reg[sig.env$ID2]==reg2)|(reg[sig.env$ID2]==reg1&reg[sig.env$ID1]==reg2)*(sig.env$t.hats<100)*sig.env$haps$f2)
        all.100.matrix[reg1,reg2] <- all.100.matrix[reg1,reg2] <- total.all.100
        sig.100.matrix[reg1,reg2] <- sig.100.matrix[reg1,reg2] <- total.sig.100
    }
}


plot((all.matrix/colSums(all.matrix))["America",], (sig.matrix/colSums(sig.matrix))["America",], xlab="Proportion of all f2 variants that America shares with each region", ylab=bquote("Proportion of signature"~.(sig.name)~f[2]~"mutations that America shares with each region"))
text((all.matrix/colSums(all.matrix))["America",], (sig.matrix/colSums(sig.matrix))["America",], regions, pos=4, cex=0.5)
abline(0,1,col="red")

for(regs in regions){
    plot((all.matrix/colSums(all.matrix))[regs,], (sig.matrix/colSums(sig.matrix))[regs,], xlab="Proportion of all f2 variants shared with each region", ylab=bquote("Proportion of signature"~.(sig.name)~f[2]~"mutations shared with each region"), main=regs)
    text((all.matrix/colSums(all.matrix))[regs,], (sig.matrix/colSums(sig.matrix))[regs,], regions, pos=4, cex=0.5)
    abline(0,1,col="red")
}
