rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))

studies_data <- read.csv(paste0(wd,"Data/Raw/studies_results.csv"))


si_metric <- studies_info[match(studies_data$Study.ID,studies_info$id),c("id","species_english","year","taxonomic.class")]

sd_metrics<-strsplit(studies_data$Metrics.reported,", ")
sd_metrics2<-do.call(rbind,lapply(1:length(sd_metrics),function(i){
  if(length(sd_metrics[[i]])>0) data.frame(id=studies_data[i,"Study.ID"],year=si_metric[i,"year"],metric=sd_metrics[[i]])
}))

# table(sd_metrics2[,"metric"])
par(mar=c(4,15,1,1))
barplot(table(sd_metrics2[,"metric"])[c("relatedness distribution","average relatedness","inbreeding","mean pedigree depth","max pedigree depth","sibship size","number of founders","grandparents","siblings","maternities/paternities","Number of individuals in the pedigree", "None")],las=2, horiz = TRUE, col=palette.colors()[c(2,2,3,4,4,5,6,6,6,6,6,7)+1], xlab="Number of Studies")


par(mar=c(4,8,1,1))
barplot(cbind(pruning=table(studies_data$Is.pruning.mentioned.),
simulation=table(studies_data$Are.pre.analysis.simulations.mentioned.),
data=c(NA,NA)
), col=palette.colors()[2:3], las=2, horiz = TRUE, xlab="Number of Studies", ylim=c(0,5))
barplot(cbind(c(NA,NA,NA),c(NA,NA,NA),table(ifelse(studies_data$A.or.invA.matrix.provided.=="Yes", "A", ifelse(studies_data$Pedigree.provided.=="Yes", "ped", "1"))))
,add=TRUE, horiz = TRUE, col=palette.colors()[c(2,4,5)], xaxt="n")
legend("top", c("No","Yes","A matrix","Pedigree"),pch=16, col=palette.colors()[2:5])





barplot(table(studies_info$year))

barplot(table(studies_data$Pedigree.provided.,si_metric$year),add=TRUE, col=c("cornflowerblue","blue"))
legend("topleft",c("no pedigree available","pedigree available","not extracted"), fill=c("cornflowerblue","blue","grey"))








metric_year <- table(sd_metrics2$year,sd_metrics2$metric)
all_n <-table(si_metric$year)
year <-as.numeric(rownames(metric_year))

plot(NA,ylim=c(0,max(metric_year)),xlim=c(1992,2020))
lines(metric_year[,"siblings"]~as.numeric(rownames(metric_year)))
lines(metric_year[,"None"]~as.numeric(rownames(metric_year)))

plot(NA,ylim=c(0,1),xlim=c(1992,2020))
lines(metric_year[,"siblings"]/all_n~as.numeric(rownames(metric_year)), col="blue")
lines(metric_year[,"None"]/all_n~as.numeric(rownames(metric_year)), col="red")
lines(metric_year[,"Number of individuals in the pedigree"]/all_n~as.numeric(rownames(metric_year)))

plot(as.vector(metric_year[,"relatedness distribution"]/all_n)~as.numeric(rownames(metric_year)), type="b")

plot(as.vector(metric_year[,"maternities/paternities"]/all_n)~as.numeric(rownames(metric_year)), type="b")

par(mfrow=c(4,4),mar=c(4,4,1,1))
for( i in 1:ncol(metric_year)){
  metric <- as.vector(metric_year[,i]/all_n)
  plot(metric~year, ylab=colnames(metric_year)[i], type="l", col="grey")
  points(metric~year, col=c("red","blue")[as.factor(year>2010)], pch=19)
}
