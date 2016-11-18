## Analyse the SGDP total mutation rate from ext table 1 of the sgdp paper.

NonWEurNonAfr <- c("Oceania", "America", "SouthAsia", "CentralAsiaSiberia", "EastAsia")
NonAfr <- c("WestEurasia", "Oceania", "America", "SouthAsia", "CentralAsiaSiberia", "EastAsia")
NonSan <- c("WestEurasia", "Oceania", "America", "SouthAsia", "CentralAsiaSiberia", "EastAsia", "Africa")

data <- read.table("~/spectrum/code/sgdp_mtuatation_accumulation_table.txt", header=TRUE, as.is=TRUE)
data$Khoesan<-ifelse(data$Pop1=="Khoesan", 1, ifelse(data$Pop2=="Khoesan", -1, 0))
data$Africa<-ifelse(data$Pop1=="Africa", 1, ifelse(data$Pop2=="Africa", -1, 0))
data$WestEurasia<-ifelse(data$Pop1=="WestEurasia", 1, ifelse(data$Pop2=="WestEurasia", -1, 0))
data$NonSan<-ifelse(data$Pop1%in%NonSan, 1, ifelse(data$Pop2%in%NonSan, -1, 0))
data$NonAfr<-ifelse(data$Pop1%in%NonAfr, 1, ifelse(data$Pop2%in%NonAfr, -1, 0))
data$NonWEurNonAfr<-ifelse(data$Pop1%in%NonWEurNonAfr, 1, ifelse(data$Pop2%in%NonWEurNonAfr, -1, 0))
