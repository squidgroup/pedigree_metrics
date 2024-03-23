dd<-read.csv("/Users/joelpick/github/maternal_effects/Data/Raw/young_postma_2023_data.csv")

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
