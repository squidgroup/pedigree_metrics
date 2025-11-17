rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))

studies_data <- read.csv(paste0(wd,"Data/Raw/studies_results.csv"))

nrow(studies_info)
length(unique(studies_data$Study.ID))
nrow(studies_data)
# missinf 147, 167


si_metric <- studies_info[match(studies_data$Study.ID,studies_info$id),c("id","species_english","year","taxonomic.class")]

sd_metrics<-strsplit(studies_data$Metrics.reported,", ")
sd_metrics2<-do.call(rbind,lapply(1:length(sd_metrics),function(i){
  if(length(sd_metrics[[i]])>0) data.frame(id=studies_data[i,"Study.ID"],year=si_metric[i,"year"],metric=sd_metrics[[i]])
}))
sd_metrics2$metric <- ifelse(sd_metrics2$metric=="Number of individuals in the pedigree","records",sd_metrics2$metric)

table(sd_metrics2$metric)

metric_year <- table(sd_metrics2$year,sd_metrics2$metric)
year <-as.numeric(rownames(metric_year))
all_n <-table(si_metric$year)


time_plot <- function(y, ...){
  plot(y~year, type="l", ylim=c(0,1),...)
  points(y~year, pch=19, col=c("red","blue")[as.factor(year>2010)])  
}




setEPS()
pdf(paste0(wd,"Figures/Fig_metrics.pdf"), height=6, width=12)

{
layout(matrix(c(1,1,2,4,3,5),nrow=2,ncol=3), width=c(3,2,2))
# table(sd_metrics2[,"metric"])
par(mar=c(4,10,1,1))
barplot(table(sd_metrics2[,"metric"])[c("relatedness distribution","average relatedness","inbreeding","mean pedigree depth","max pedigree depth","sibship size","number of founders","grandparents","siblings","maternities/paternities","records", "None")],las=2, horiz = TRUE, col=palette.colors()[c(2,2,3,4,4,5,6,6,6,6,6,7)+1], xlab="Number of Studies", xlim=c(0,nrow(studies_data)))


par(mar=c(4,4,1,1))

barplot(table(studies_data$Pedigree.provided.,si_metric$year), col=c("cornflowerblue","blue"), xlab="Number of Studies")
legend("topleft",c("no pedigree available","pedigree available"), fill=c("cornflowerblue","blue","grey"))

ped_avail <- table(studies_data$Pedigree.provided.,si_metric$year)
ainv_avail <- table(studies_data$A.or.invA.matrix.provided.,si_metric$year)

plot(ped_avail["Yes",]/colSums(ped_avail)~year, ylab="Proportion with data available", type="b", pch=19, ylim=c(0,1), xlab="Year")
points((ped_avail["Yes",] + ainv_avail["Yes",])/colSums(ped_avail)~year, pch=21)
legend("topleft", c("ainv or pedigree","pedigree"), pch=c(21,21),pt.bg=c(0,1))


sims<- table(studies_data$Are.pre.analysis.simulations.mentioned.,si_metric$year)
pruning<- table(studies_data$Is.pruning.mentioned.,si_metric$year)

time_plot(sims["Yes",]/colSums(sims),ylab="Proportion mentioning simulations", xlab="Year")
legend("topleft", c("2010 or earlier","post 2010"), pch=c(19,19),col=c("red","blue"))
time_plot(pruning["Yes",]/colSums(pruning),ylab="Proportion mentioning pruning", xlab="Year")

}
dev.off()



metric_names<-c("relatedness distribution","average relatedness","inbreeding","mean pedigree depth","max pedigree depth","sibship size","number of founders","grandparents","siblings","maternities/paternities","records", "None")


setEPS()
pdf(paste0(wd,"Figures/Fig_metrics_all.pdf"), height=8, width=10)

par(mfrow=c(4,3),mar=c(4,4,1,1))
for( i in metric_names){
  time_plot(as.vector(metric_year[,i]/all_n), ylab=i)
}
dev.off()

year_bin <- as.factor(year>2010)

summary(glm(cbind(metric_year[,"average relatedness"],all_n-metric_year[,"average relatedness"])~year + year_bin, family="binomial"))


colnames(metric_year)

table(metric_year$metric)