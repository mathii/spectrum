## Output EMu and run it.
source("~/spectrum/code/spectrumlib.R")


n <- 2

inname <- paste0("~/spectrum/data/count_matrix.n", n,".txt")
count <- read.table(inname, header=TRUE, as.is=TRUE )

tcount <- t(count)
write.table(tcount, paste0("~/spectrum/EMu/tcount_matrix.n", n, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

res <- system(paste0("~/Packages/EMu/build/EMu --mut ~/spectrum/EMu/tcount_matrix.n",2,".txt --opp human-genome --pre ~/spectrum/EMu/EMu_results_n",n), intern=TRUE)
rg=res[grep("According to BIC", res)]
pos <- regexpr("[0-9] spectra", rg)[1]
rank <- as.numeric(substr(rg, pos, pos))

components <- read.table(paste0("~/spectrum/EMu/EMu_results_n",n,"_",rank,"_assigned.txt"), as.is=TRUE)
rownames(components) <- rownames(tcount)
loadings <- read.table(paste0("~/spectrum/EMu/EMu_results_n",n,"_",rank,"_ml_spectra.txt"), as.is=TRUE)
loadings <- t(loadings)
rownames(loadings) <- colnames(tcount)

components <- components/rowSums(components)

write.table( components, paste0("~/spectrum/plots/Components_EMu", ".n", n, ".r", rank, ".txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )

pdf(paste0("~/spectrum/plots/","Components_EMu.n", n, ".r", rank, ".pdf"), 12, 12)
plot.components(components, name.map, src, cols=cols, n.components=rank, layout=c(2,2))
dev.off()

write.table(loadings, paste0("~/spectrum/plots/","Loadings_EMu.n", n, ".r", rank, ".txt"), row.names=T, col.names=F, quote=F) 

pdf(paste0("~/spectrum/plots/","Loadings_EMu.n", n, ".r", rank, ".pdf"), width=12, height=rank*4)
plot.loadings(loadings, n.loadings=rank)
dev.off()
