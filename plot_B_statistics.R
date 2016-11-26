## Plot the European fold-excess as a function of B statistic.
## also used to plot recombination rate. 

source("~/spectrum/code/spectrumlib.R")

n <- 2
what <- "Bstat"                         
## what <- "rrate.1kb"
bs <- 0:9
excess=-1+0*bs

sig <- c("TCT.T", "TCC.T", "CCC.T", "ACC.T")
sig.rc <- c( "AGA.A", "GGA.A", "GGG.A", "GGT.A" )
sig.name <- "1"
in.ind <- NA
in.pops <- c("WestEurasia")
out.pops <- c("WestEurasia", "SouthAsia")

## sig <- c("ACG.T", "CCG.T", "GCG.T", "TCG.T")
## sig.rc <- c("CGT.A", "CGG.A", "CGC.A", "CGA.A" )
## in.ind <- c("S_Chane.1", "S_Piapoco.2", "S_Quechua.3", "S_Mayan.1", "S_Mayan.2", "S_Quechua.1", "S_Nahua.1", "S_Quechua.2", "S_Nahua.2", "S_Zapotec.1", "S_Mixtec.1")
## in.pops <- c("America")
## out.pops <- c("America")
## sig.name <- "2"


data <- data.frame()

for(bi in 1:10){
    data <- read.table(paste0("~/spectrum/", what, "/chr1.n",n,".b", bs[bi], ".txt"), as.is=TRUE, header=TRUE)
    row.names(data) <- data[,1]
    data <- data[,3:NCOL(data)]
    for(chr in 2:22){
        pp <- read.table(paste0("~/spectrum/", what, "/chr", chr ,".n",n,".b", bs[bi], ".txt"), as.is=TRUE, header=TRUE)
        row.names(pp) <- pp[,1]
        data <- data + pp[,3:NCOL(pp)]
    }

    proportion <- colSums(data[c(sig, sig.rc),])/colSums(data)

    if(!all(is.na(in.ind))){
        in.prop <- proportion[names(proportion)%in% in.ind]
        out.prop <- proportion[!(names(proportion)%in% in.ind)]        
    }else{
        in.prop <- proportion[name.map[names(proportion)]%in% in.pops]
        out.prop <- proportion[!(name.map[names(proportion)]%in% out.pops)]
    }
    if(bi==1){
        props <- data.frame(prop=in.prop, in.group=rep("in", length(in.prop)), bq=rep(bs[bi], length(in.prop)))
    } else{
        props <- rbind(props, data.frame(prop=in.prop, in.group=rep("in", length(in.prop)), bq=rep(bs[bi], length(in.prop))))
    }
    props <- rbind(props, data.frame(prop=out.prop, in.group=rep("out", length(out.prop)), bq=rep(bs[bi], length(out.prop))))
    
    excess[bi] <- mean(proportion[name.map[names(proportion)]%in% in.pops])-mean(proportion[!(name.map[names(proportion)]%in% out.pops)])
    
}

xlab=""
if(what=="Bstat"){
    xlab <- "B statistic decile"
}else if(grepl("^rrate", what)){
    xlab <- "Recombination rate decile"
}

pdf(paste0("~/spectrum/plots/", what, "_sig", sig.name, ".pdf")) 
boxplot(prop~in.group+bq, data=props, col=rep(c( cols[in.pops], "white"), 10), frame=F, xaxt="n", xlab=xlab, ylab=paste0("Proportion of signature ", sig.name ," mutations"))
mtext(paste0((0:9)*10, "-", (1:10)*10, "%"), side=1, at=2*(1:10)-0.5, cex=0.75, line=-1)

## if(grepl("^rrate", what)){
##     print("rrate")
## }else if(sig.name=="1"){
##     mod<-lm(prop~as.factor(in.group)+bq, data=props)
##     ss=summary(mod)
##     a=ss$coefficients[1,1]
##     b=ss$coefficients[3,1]
##     c=ss$coefficients[2,1]
##     abline(a-b/2,b/2, col="blue", lty=2, lwd=2)
##     abline(a+c-b,b/2, col="blue", lty=2, lwd=2)
## }else if(sig.name=="2"){
##     props$bq2 <- props$bq*props$bq
##     props$bq3 <- props$bq2*props$bq

##     mod<-lm(prop~as.factor(in.group)+bq+bq2, data=props)
##     ss=summary(mod)
##     a=ss$coefficients[1,1]
##     b=ss$coefficients[2,1]
##     c=ss$coefficients[3,1]
##     d=ss$coefficients[4,1]

##     xs=seq(-1,10,length.out=100)
##     v1 <- a+c*xs+d*xs*xs
##     v2 <-  a+b+c*xs+d*xs*xs
##     lines(2*xs+1, v1, col="blue", lty=2, lwd=2) 
##     lines(2*xs+2, v2,  col="blue", lty=2, lwd=2)
## }
dev.off()

print(summary(mod))
