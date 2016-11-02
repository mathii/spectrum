## Plot the psmc rate (CpG vs nonCpg mutations).

data <- read.table("~/spectrum/psmc/output/S_Yoruba-1.0.txt", as.is=T)
data.CpG <- read.table("~/spectrum/psmc/output/S_Yoruba-1.CpG.0.txt", as.is=T)
data.nonCpG <- read.table("~/spectrum/psmc/output/S_Yoruba-1.nonCpG.0.txt", as.is=T)

plot(data[data[,1]>1000,1], data[data[,1]>1000,2], type="s", col="red", log="x", ylim=c(0,3), ylab="Ne", xlab="Time (years)")

CpG.scale=0.1
lines(data.CpG[,1]/CpG.scale, data.CpG[,2]/CpG.scale, type="s", col="blue")

nonCpG.scale=1-CpG.scale
lines(data.nonCpG[data.nonCpG[,1]>1000,1]/nonCpG.scale, data.nonCpG[data.nonCpG[,1]>1000,2]/nonCpG.scale, type="s", col="green")
