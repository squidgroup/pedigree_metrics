rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))
load(paste0(wd,"Data/Processed/ped_mde.Rdata"))
load(paste0(wd,"Data/Processed/ped_mbias.Rdata"))

all_sim <- all_sim[,!colnames(all_sim) %in% c("D33.csv","D38.txt")]

lambdas <- lambdas[!names(lambdas) %in% c("D33.csv","D38.txt")]
all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.3)

all_mbias_sims <- all_mbias_sims[,!colnames(all_mbias_sims) %in% c("D33.csv","D38.txt")]



files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
id <- gsub("\\.|[a-z]","",files)

studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))
dryad_ped_info <- read.csv(paste0(wd,"Data/Raw/dryad_ped_info.csv"))
all_ped_info <- rbind(studies_info[,c("id","species_english","taxonomic.class")], dryad_ped_info[,c("id","species_english","taxonomic.class")])

si <- all_ped_info[match(id,all_ped_info$id),c("id","species_english","taxonomic.class")]


par(mfrow=c(3,8),mar=c(2,2,4,1))
for(i in 1:23){
 plot(all_stat[i,],all_mbias_sims[1,],pch=19,main=rownames(all_stat)[i], log="x", col=c(4,3,2,1)[as.factor(si$taxonomic.class)])
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}
plot(all_mbias_sims[1,]~all_sv,pch=19, main="Visscher & Goddard 2015", log="x", col=c(4,3,2,1)[as.factor(si$taxonomic.class)])
