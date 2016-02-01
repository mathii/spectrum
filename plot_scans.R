## Manhattan and QQ plots of the GWAS results.

source("~/selection/code/lib/mh_plot_lib.R")
plots <- c("child", "maternal", "paternal")
type <- c("add", "dom", "rec")

for( what in plots ){
    cat(paste0(what, "\n"))
    for( typ in type){
        cat(paste0(typ, "\n"))
        tagstr<-paste0("frequentist_",typ,"_pvalue")
        cat("1")
        results <- read.table(paste0("~/spectrum/gonl/results/", what, ".chr1.test.txt.gz"), as.is=TRUE, header=TRUE, sep=" ", fill=TRUE)
        results<-results[,c("chromosome", "position", tagstr)]
        for(chr in 2:22){
            cat(paste0("\r", chr))
            rs <- read.table(paste0("~/spectrum/gonl/results/", what, ".chr", chr,".test.txt.gz"), as.is=TRUE, header=TRUE, sep=" ", fill=TRUE)
            results<-rbind(results, rs[,c("chromosome", "position", tagstr)])
        }
        cat("\n")
  ## a=qchisq(median(results$frequentist_add_pvalue), df=1, lower.tail=TRUE)
  ## b=qchisq(0.5, df=1, lower.tail=TRUE)
  ## lambda=a/b

        tagstr<-paste0("frequentist_",typ,"_pvalue")
        results<-results[,c("chromosome", "position", tagstr)]
        results<-results[!is.na(results[,3]),]
        
        colnames(results) <- c("CHR", "POS", "PVAL")
        png(paste0("~/spectrum/gonl/results/", what,"_", typ, ".result.png"), width=1800, height=600)
        MH.plot(results, thin.thin=1000)
        dev.off()
        
        png(paste0("~/spectrum/gonl/results/", what, "_", typ, ".qq.png"), width=600, height=600)
        qqPlotOfPValues(results$PVAL, thin.thin=1000)
        dev.off()
    }
}
