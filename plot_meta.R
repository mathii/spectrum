## Manhattan and QQ plots of the GWAS meta-analyis results.

source("~/selection/code/lib/mh_plot_lib.R")

what="meta"
tagstr<-"P_value"
cat("1")
results <- read.table(paste0("~/spectrum/gonl/results/", what, ".chr1.test.txt.gz"), as.is=TRUE, header=TRUE, sep=" ", fill=TRUE)
results<-results[,c("chr", "pos", tagstr)]
for(chr in 2:22){
    cat(paste0("\r", chr))
    rs <- read.table(paste0("~/spectrum/gonl/results/", what, ".chr", chr,".test.txt.gz"), as.is=TRUE, header=TRUE, sep=" ", fill=TRUE)
    results<-rbind(results, rs[,c("chr", "pos", tagstr)])
}
cat("\n")

results<-results[!is.na(results[,3]),]
colnames(results) <- c("CHR", "POS", "PVAL")

png(paste0("~/spectrum/gonl/results/", what, ".result.png"), width=1800, height=600)
MH.plot(results, thin.thin=1000)
dev.off()

png(paste0("~/spectrum/gonl/results/", what, ".qq.png"), width=600, height=600)
qqPlotOfPValues(results$PVAL, thin.thin=1000)
dev.off()
