## Replot the random initialization loadings to check whether they are consistent with cpg loadings remover

loadings <- read.table("~/spectrum/plots/Loadings_random_NMF.n2.r4.txt", as.is=TRUE)
rownames(loadings) <- loadings[,1]
loadings <- loadings[,2:5]

s2 <- c("ACG.T","CCG.T", "GCG.T", "TCG.T")

ll <- t(t(loadings)/colSums(loadings))
ll[s2,] <- ll[s2,]-apply(ll[s2,], 1, min)

pdf(paste0("~/spectrum/plots/","Loadings_NMF.n2.r4.random_replotted.pdf"), width=12, height=4*4)
plot.loadings(ll, n.loadings=4)
dev.off()
