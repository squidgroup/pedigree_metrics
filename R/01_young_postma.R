dd<-read.csv("Data/Raw/young_postma_2023_data.csv")

names(dd)

tab <- table(dd$method.2,dd$year)

length(unique(subset(dd,dd$method.2 %in% c("animal model (MCMC)","animal model (REML)"))$study.ID))

am <- colSums(tab[2:3,])
other <- colSums(tab[4:8,])
total <- colSums(tab)

plot(am, type="l")
lines(other)

plot(am/total~as.numeric(names(am)), type="l")

am_papers<- aggregate(estimate.ID~study.ID+	year+	reference+	species.english+	species.latin+	taxonomic.class+	population+	country,subset(dd,dd$method.2 %in% c("animal model (MCMC)","animal model (REML)")),length)
write.csv(am_papers,row.names=FALSE,file="to_extract.csv")

max(am_papers$year)



dd2 <- subset(dd,dd$method.2 %in% c("animal model (MCMC)","animal model (REML)"))
quantile(dd2$n,c(0,0.1,0.25,0.5,0.75,0.90,0.95,1),na.rm=TRUE)
hist((dd2$n), xlim=c(0,4000), breaks=seq(0,40000,100))