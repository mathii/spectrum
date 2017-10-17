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
## ylim=c(0.06262, 0.06282)

sig.name <- 2
sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
ylim=c(0.1190, 0.1210)
hi.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")

####################################################

sig.name.map <- c(4,1)
sig.map=c(2,4,3,1)

####################################################

regions <- sort(unique(info$Region))
ltys <- rep(1, length(cols))
names(ltys) <- names(cols)

rownames(info) <- gsub("-", ".", info$ID)

## if(sig.name==2){
##     info[hi.ind,"Region"] <- "America_Hi"
##     info[info[rownames(info),"Region"]=="America","Region"] <- "America_Lo"
## }

f2.data <- read.table(paste0("~/spectrum/data/count_matrix.n2", what, ".txt"), as.is=TRUE, header=TRUE)
ALL.data <- read.table(paste0("~/spectrum/data/count_matrix.nALL", what, ".txt"), as.is=TRUE, header=TRUE)

f2.data <- f2.data[,!(colnames(f2.data)=="S_Daur.1")]
ALL.data <- ALL.data[,!(colnames(ALL.data)=="S_Daur.1")]                   

f2.proportion <- colSums(f2.data[sig,])/colSums(f2.data)
ALL.proportion <- colSums(ALL.data[sig,])/colSums(ALL.data)

f2.sig.total <- colSums(f2.data[sig,])
f2.variants.total <- colSums(f2.data)
ALL.sig.total <- colSums(ALL.data[sig,])
ALL.variants.total <- colSums(ALL.data)
ALL.notsig.total <- colSums(ALL.data[!(rownames(ALL.data) %in% sig),])

pdf(paste0("~/spectrum/plots/plot_all_against_f2_sig", sig.name, what, ".pdf"))
plot(f2.proportion, ALL.proportion, col=cols[reg[names(f2.proportion)]], xlab=paste0("Proportion of f2 variants that are signature ", sig.name), ylab=paste0("Proportion of all variants, per genome, that are signature ", sig.name), pch=ifelse(grepl("^B", names(f2.proportion)), 2, 1))
legend("bottomright", c(names(cols), "S Panel", "B Panel"), col=c(cols, "black", "black"), bty="n", pch=c(rep(1, 8), 2))
dev.off()

## t.test(ALL.proportion[reg[names(ALL.proportion)]%in%c("America", "CentralAsiaSiberia", "EastAsia", "Oceania", "SouthAsia", "WestEurasia")], ALL.proportion[hi.ind])

pdf(paste0("~/spectrum/plots/plot_all_against_f2_sig_total", sig.name, what, ".pdf"))
plot(ALL.notsig.total, ALL.sig.total, col=cols[reg[names(ALL.sig.total)]] )
dev.off()

 AP<-ALL.proportion[!grepl("^B", names(ALL.proportion))]

##  t.test(AP[reg[names(AP)]=="Africa"], AP[!(reg[names(AP)] %in% c("Africa"))])
## Sig2 Africa vs non -africa P=4.6e-11 0.1083719 0.1084019
## t.test(AP[names(AP) %in% hi.ind], AP[(reg[names(AP)] %in% c("EastAsia"))])
## Sig2 Hi vs other p=0.0035

GC.AT<-substr(rownames(ALL.data),2,2) %in% c("G", "C")&substr(rownames(ALL.data),5,5)%in% c("A", "T")
AT.GC <- substr(rownames(ALL.data),2,2) %in% c("A", "T")&substr(rownames(ALL.data),5,5)%in% c("G", "C")
Other <- !(GC.AT|AT.GC)
AD<-ALL.data[!grepl("^B", colnames(ALL.data))]
Biased.GC <- data.frame()
d1 <- colSums(AD[GC.AT,])
dd <- data.frame("Type"="GC.AT", "ID"=names(d1), "Value"=c(t(d1))/mean(c(t(d1))), stringsAsFactors=FALSE)
Biased.GC <- rbind(Biased.GC, dd)
d1 <- colSums(AD[AT.GC,])
dd <- data.frame("Type"="AT.GC", "ID"=names(d1), "Value"=c(t(d1))/mean(c(t(d1))), stringsAsFactors=FALSE)
Biased.GC <- rbind(Biased.GC, dd)
d1 <- colSums(AD[Other,])
dd <- data.frame("Type"="Other", "ID"=names(d1), "Value"=c(t(d1))/mean(c(t(d1))), stringsAsFactors=FALSE)
Biased.GC <- rbind(Biased.GC, dd)
Biased.GC$pop <- reg[Biased.GC$ID]
Biased.GC$Africa <- ifelse(Biased.GC$pop=="Africa", "Africa", "Non-Africa")
boxplot(Value~Africa+Type, data=Biased.GC)

Ts<-(substr(rownames(ALL.data),2,2)=="C"&substr(rownames(ALL.data),5,5)=="T")|(substr(rownames(ALL.data),2,2)=="T"&substr(rownames(ALL.data),5,5)=="C")|(substr(rownames(ALL.data),2,2)=="A"&substr(rownames(ALL.data),5,5)=="G")|(substr(rownames(ALL.data),2,2)=="G"&substr(rownames(ALL.data),5,5)=="A")
Tv <- !Ts
Ts.Tv <- data.frame()
d1 <- colSums(AD[Ts,])
dd <- data.frame("Type"="Ts", "ID"=names(d1), "Total"=d1, "Relative"=c(t(d1))/mean(c(t(d1))), stringsAsFactors=FALSE)
Ts.Tv <- rbind(Ts.Tv, dd)
d1 <- colSums(AD[Tv,])
dd <- data.frame("Type"="Tv", "ID"=names(d1), "Total"=d1, "Relative"=c(t(d1))/mean(c(t(d1))), stringsAsFactors=FALSE)
Ts.Tv <- rbind(Ts.Tv, dd)
Ts.Tv$pop <- reg[Ts.Tv$ID]
Ts.Tv$Africa <- ifelse(Ts.Tv$pop=="Africa", "Africa", "Non-Africa")
boxplot(Value~Africa+Type, data=Ts.Tv)
