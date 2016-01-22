## Manhattan and QQ plots of the GWAS results.

source("~/selection/code/lib/mh_plot_lib.R")
plots <- c("child", "maternal", "paternal")

for( what in plots ){
  results <- read.table(paste0("~/spectrum/gonl/results/", what, ".results.txt.gz"), as.is=TRUE, header=TRUE)
  results<-results[!is.na(results[,4]),2:4]

  ## a=qchisq(median(results$frequentist_add_pvalue), df=1, lower.tail=TRUE)
  ## b=qchisq(0.5, df=1, lower.tail=TRUE)
  ## lambda=a/b

  colnames(results) <- c("CHR", "POS", "PVAL")
  png(paste0("~/spectrum/gonl/results/", what, ".result.png"), width=1800, height=600)
  MH.plot(results, thin.thin=1000)
  dev.off()

  png(paste0("~/spectrum/gonl/results/", what, ".qq.png"), width=600, height=600)
  qqPlotOfPValues(results$PVAL, thin.thin=1000)
  dev.off()

}
