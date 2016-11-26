## Compare the decompositions from PCA, sparse PCA and Non-negative Matrix factorisation.
set.seed(12345)

library("RColorBrewer")
library("NMF")
library("fastICA")

source("~/spectrum/code/spectrumlib.R")

spec <- "spectrum"

n <- 2
rank <- 4
tag <- ""
method <- "ica"
diagnostics <- FALSE
subtract <- FALSE

cA <- commandArgs(TRUE)
if(length(cA)>0){
  n <- cA[1]
}
if(length(cA)>1){
    rank <- as.numeric(cA[2])
}
if(length(cA)>2 & !(cA[3] %in% c("", "."))){
    tag <- paste0(".", cA[3])
}
if(length(cA)>3){
    diagnostics <- as.logical(as.numeric(cA[4]))
}
if(length(cA)>4){
    method <- cA[5]
}
if(length(cA)>5){
    subtract <- as.logical(as.numeric(cA[6]))
}
if(length(cA)>6){
    spec <- cA[7]
}


inname <- paste0("~/spectrum/data/", spec ,"_matrix.n", n ,tag, ".txt")

freq2 <- read.table(inname, header=TRUE, as.is=TRUE )

## Subtract off minimal mutations. 

if(subtract){
    freq2 <- freq2-apply(freq2, 1, min)
}
if(subtract & spec=="spectrum"){
    freq2["ATA.C",] <- 0.00001         #Null row.
}
if(spec=="totalnorm"){
  freq2 <- t(t(freq2)/colSums(freq2))
}
## Non-negative Matrix factorization
seed <- ifelse(method=="random", rep(123456, 6), method)

outtag <-  paste0(ifelse(spec=="spectrum", "", paste0("_", spec)), ifelse(subtract, "_subtract", ""), "_", method)

logfile <- paste0("~/spectrum/plots/","Info",  outtag, "_NMF.n", n, ".r", rank, tag, ".log")

if(!diagnostics){
cat(paste0("ARGS: ", cA, "\n"), file=logfile)
nnegmf <- nmf(as.matrix(freq2),rank=rank, nrun=200, seed = seed )

coeff <- t(coef(nnegmf))

if(spec=="count"){
    coeff <- coeff/rowSums(coeff)
}

pdf(paste0("~/spectrum/plots/","Components",  outtag, "_NMF.n", n, ".r", rank, tag, ".pdf"), 12, 12)
if(n=="2" & rank==4 & method=="ica" & spec=="totalnorm"){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(2,3,3,3), yploti=c(4,1,4,2))
}else if(n=="2" & rank==4 & method=="ica"){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(3,4,1,2), yploti=c(1,2,4,3))
}else if(n=="2" & rank==4 & method=="random"){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(1,2,3,4), yploti=c(4,3,2,1))
}else if(n=="3" & rank==3){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(2,2), yploti=c(1,3))
}else if(n=="3" & rank==4){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2), xploti=c(4,2,3,4), yploti=c(3,1,1,2))
}else if(rank==2){
    plot.components(coeff, name.map, src, cols=cols, n.components=rank)
}else{
    plot.components(coeff, name.map, src, cols=cols, n.components=rank, layout=c(2,2))
}
dev.off()

bas <- basis(nnegmf)
nbas <- substr(rownames(bas),1,3)

ref <- read.table("~/spectrum/ref/trinucleotide_context.txt", as.is=TRUE, header=FALSE)
cmpnames <- sapply(ref[,1], reverse.complement)
ref[,1] <- ifelse(ref[,1] %in% nbas, ref[,1], cmpnames)
ref<-aggregate(ref[,2], by=list(ref[,1]), sum)
ref[,2] <- ref[,2]/sum(as.numeric(ref[,2]))
refmap <- ref[,2]
names(refmap) <-ref[,1]
scale <- refmap[nbas]

## if(n==3 & rank==4){bas <- bas[,c(2,4,3,1)]}
## if(n==3 & rank==3){bas <- bas[,c(2,1,3)]}
pdf(paste0("~/spectrum/plots/","Loadings", outtag,"_NMF.n", n, ".r", rank, tag, ".pdf"), width=12, height=rank*4)
plot.loadings(bas/scale, n.loadings=rank)
dev.off()

write.table(coeff, paste0("~/spectrum/plots/","Components",  outtag, "_NMF.n", n, ".r", rank, tag, ".txt"), row.names=T, col.names=F, quote=F) 
write.table(bas/scale, paste0("~/spectrum/plots/","Loadings",  outtag, "_NMF.n", n, ".r", rank, tag, ".txt"), row.names=T, col.names=F, quote=F) 

features <- lapply(extractFeatures(nnegmf), function(x){rownames(freq2)[x]})
residuals <- residuals(nnegmf)
cat(paste0("RES: ", residuals, "\n"), file=logfile, append=TRUE)
for(i in 1:length(features)){
    cat(paste0("F", i, ": ", features[[i]], "\n", collapse=" "), file=logfile, append=TRUE)
}
## Print errors for best fit. 
abs.errors <- as.matrix((abs((fitted(nnegmf)-freq2))))
sq.errors <- abs.errors*abs.errors
log.errors <- as.matrix((abs(log(fitted(nnegmf)/freq2))))
cat(paste0("MAXAD: ", max(abs.errors), "\n"), file=logfile, append=TRUE)
cat(paste0("MNAD: ", mean(abs.errors), "\n"), file=logfile, append=TRUE)
cat(paste0("MAXSD: ", max(sq.errors), "\n"), file=logfile, append=TRUE)
cat(paste0("MNSD: ", mean(sq.errors), "\n"), file=logfile, append=TRUE)
cat(paste0("MAXLD: ", max(log.errors), "\n"), file=logfile, append=TRUE)
cat(paste0("MNLD: ", mean(log.errors), "\n"), file=logfile, append=TRUE)
}

## outliers <- c("B_Crete.1", "B_Crete.2"   ,"B_Dai.4"    ,"B_Han.3", "B_Australian.4", "S_Miao.1", "S_Miao.2", "S_Russian.1", "S_Mongola.1")
## outliers <- c("B_Mixe.1", "S_Pima.2")
## freq2.out <- freq2[,!(colnames(freq2) %in% outliers)]
## nnegmf.out=nmf(as.matrix(freq2.out),rank=rank, seed="ica", nrun=20)
## pdf(paste0("~/spectrum/plots/","Components_",  ifelse(spec=="spectrum", "", paste0(spec, "_")), "NMF.n", n, ".r", rank, tag, "_out.pdf"), 12, floor(sqrt(rank))*6)
## plot.components(t(coef(nnegmf.out)), name.map, src, cols=cols, n.components=rank)
## dev.off()

## pdf(paste0("~/spectrum/plots/","Loadings_", ifelse(spec=="spectrum", "", paste0(spec, "_")),"NMF.n", n, ".r", rank, tag, "_out.pdf"), width=12, height=12)
## plot.loadings(basis(nnegmf.out), n.loadings=rank)
## dev.off()

## plot(coef(nnegmf)[3,], coef(nnegmf)[4, ])
## identify(coef(nnegmf)[3,], coef(nnegmf)[4, ], colnames(freq2.out))

if(diagnostics){
    nnegmf.test <- nmf(as.matrix(freq2), 2:8, seed=seed, nrun=200)

    pdf(paste0("~/spectrum/plots/","Diagnostics",  outtag, "NMF.n", n, tag, ".pdf"), 12, 12)
    ## Yikes!
    plot(plot(nnegmf.test))
    dev.off()

    pdf(paste0("~/spectrum/plots/","DiagnosticsSmall",  outtag, "NMF.n", n, tag, ".pdf"), 12, 4)
    plot(plot(nnegmf.test, what=c("dispersion", "rss", "silhouette")))
    dev.off()
}
