rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))



# null<- function()
# sims<-c(parallel::mclapply(1:200, function(i){
# 	data<-data.frame(p=rnorm(nrow(ped),0,1),animal=as.factor(ped[,1]))
# 	assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 
# 	mod_1 <- asreml::asreml(
# 			fixed= p~1
# 			, random= ~vm(animal,ped.ainv)
# 		  , data= data, trace=FALSE)
# 	rm("ped.ainv", envir = .GlobalEnv)
# 	cat(i," ")
# 	summary(mod_1)$varcomp["vm(animal, ped.ainv)",1]
# },mc.cores=6), recursive=TRUE)
# hist(sims)



files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
type <- sapply(strsplit(files,"\\."),function(x) x[2])
peds<-lapply(1:length(files), function(i){
  ped <- read.table(paste0(wd,"Data/Raw/pedigrees/",files[i]),header=TRUE,sep=if(type[i]=="txt"){""}else{","}, na.strings=c("NA",0))
  ped2 <- data.frame(animal=ped[,1],dam=ped[,"dam"],sire=ped[,"sire"])
  pedantics::fixPedigree(ped2)
})
names(peds) <- files
lapply(peds,head)
lapply(peds,tail)
lapply(peds,nrow)

all_stat <- do.call(cbind,lapply(peds,ped_stat))

all_m <- do.call(c,parallel::mclapply(peds,samp_var_h2,0.3,mc.cores=6))
all_sim <- sapply(peds,samp_var_sim,0.3)

all_mbias_sims <- sapply(peds,samp_var_sim,h2=0.1,m2=0.2)

save(all_stat, all_m, all_sim,all_mbias_sims,file=paste0(wd,"Data/Processed/ped_stats.Rdata"))


all_mbias_sims_prec <- sapply(peds,samp_var_sim_mat,h2=0.3,m2=0.2)


# cor(t(rbind(do.call(cbind,all_stat),samp_var=all_m)))

# pairs(t(rbind(do.call(cbind,all_stat),samp_var=all_m)))

all_stat==0

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))


par(mfrow=c(5,5),mar=c(2,2,3,1))
for(i in 1:nrow(all_stat)){
 plot(log(all_stat[i,]),1/sqrt(all_sim[2,]),pch=19,main=rownames(all_stat)[i])
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}
# apply(all_stat,1, function(y) plot(1/sqrt(all_m),log(y)))

plot(all_sim[2,],all_m,pch=19, log="xy")
abline(0,1)


par(mfrow=c(5,5),mar=c(2,2,3,1))
for(i in 1:nrow(all_stat)){
 plot(log(all_stat[i,]),all_mbias_sims[1,],pch=19,main=rownames(all_stat)[i])
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}
plot(all_mbias_sims[1,]~all_m,pch=19, log="xy")
abline(0,1)


par(mfrow=c(5,5),mar=c(2,2,3,1))
for(i in 1:nrow(all_stat)){
 plot(log(all_stat[i,]),1/sqrt(all_mbias_sims_prec[2,]),pch=19,main=rownames(all_stat)[i])
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}
plot(all_mbias_sims_prec[2,],all_m,pch=19, log="xy")
abline(0,1)

prec_va_vm <- 1/sqrt(all_mbias_sims_prec[2,])
prec_va <- 1/sqrt(all_sim)

par(mfrow=c(5,5),mar=c(2,2,3,1))
for(i in 1:nrow(all_stat)){
plot(log(all_stat[i,]),
	(prec_va_vm - prec_va)/prec_va,pch=19,main=rownames(all_stat)[i])
abline(h=0)
}


studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))

studies_data <- read.csv(paste0(wd,"Data/Raw/studies_results.csv"))

id <- gsub("\\.|[a-z]","",files) #sapply(strsplit(files,"\\."),function(x) x[1])




si <- studies_info[match(id,studies_info$id),c("id","species_english","taxonomic.class")]

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))

si$n <- all_stat[1,]
boxplot(n~taxonomic.class,si)


si_metric <- studies_info[match(studies_data$Study.ID,studies_info$id),c("id","species_english","year","taxonomic.class")]

sd_metrics<-strsplit(studies_data$Metrics.reported,", ")
sd_metrics2<-do.call(rbind,lapply(1:length(sd_metrics),function(i){
	if(length(sd_metrics[[i]])>0) data.frame(id=studies_data[i,"Study.ID"],year=si_metric[i,"year"],metric=sd_metrics[[i]])
}))

table(sd_metrics2[,"metric"])
par(mar=c(2,10,1,1))
barplot(table(sd_metrics2[,"metric"])[c("pedigree error","relatedness distribution","average relatedness","inbreeding","mean pedigree depth","max pedigree depth","sibship size","number of founders","grandparents","siblings","maternities/paternities","Number of individuals in the pedigree", "None")],las=2, horiz = TRUE, col=c(1,2,2,3,4,4,5,6,6,6,6,6,7))


metric_year <- table(sd_metrics2$year,sd_metrics2$metric)
all_n <-table(si_metric$year)

plot(NA,ylim=c(0,max(metric_year)),xlim=c(1992,2020))
lines(metric_year[,"siblings"]~as.numeric(rownames(metric_year)))
lines(metric_year[,"None"]~as.numeric(rownames(metric_year)))

plot(NA,ylim=c(0,1),xlim=c(1992,2020))
lines(metric_year[,"siblings"]/all_n~as.numeric(rownames(metric_year)))
lines(metric_year[,"None"]/all_n~as.numeric(rownames(metric_year)))
lines(metric_year[,"Number of individuals in the pedigree"]/all_n~as.numeric(rownames(metric_year)))


for(i in 1:ncol(metric_year)) lines(metric_year[,i]~as.numeric(rownames(metric_year)))

barplot(cbind(pruning=table(studies_data$Is.pruning.mentioned.),
simulation=table(studies_data$Are.pre.analysis.simulations.mentioned.),
pedigree=table(studies_data$Pedigree.provided.),
 Amatrix=table(studies_data$A.or.invA.matrix.provided.)
))

table(studies_data$Pedigree.provided.,si_metric$year)
barplot(table(studies_data$Pedigree.provided.,si_metric$year))


###############

pairs(t(all_stat))
pairs(t(all_stat)[,1:6])
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
   par(usr = c(0, 1, 0, 1))
   r <- abs(cor(x, y,use="pairwise.complete.obs"))
   txt <- format(c(r, 0.123456789), digits = digits)[1]
   txt <- paste0(prefix, txt)
   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
   text(0.5, 0.5, txt, cex = cex.cor * r)
}
all_stat[all_stat==0]<-NA
pairs(log(t(all_stat)[,1:12]), upper.panel = panel.cor,gap=0, pch=19)
pairs(log(t(all_stat)[,c(1,13:23)]), upper.panel = panel.cor,gap=0, pch=19)


cbind(rownames(all_stat),eigen(cor(t(all_stat)))$vectors[,2])

princomp(t(all_stat),cor=TRUE)