rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))
load(paste0(wd,"Data/Processed/ped_mde.Rdata"))

all_sim <- all_sim[,!colnames(all_sim) %in% c("D33.csv","D38.txt")]

lambdas <- lambdas[!names(lambdas) %in% c("D33.csv","D38.txt")]
all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.3)

all_mde <- all_mde[,!colnames(all_mde) %in% c("D33.csv","D38.txt")]

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
id <- gsub("\\.|[a-z]","",files)

studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))
dryad_ped_info <- read.csv(paste0(wd,"Data/Raw/dryad_ped_info.csv"))
all_ped_info <- rbind(studies_info[,c("id","species_english","taxonomic.class")], dryad_ped_info[,c("id","species_english","taxonomic.class")])

si <- all_ped_info[match(id,all_ped_info$id),c("id","species_english","taxonomic.class")]


par(mfrow=c(3,8),mar=c(2,2,4,1))#mar=c(2,2,3,1))
for(i in 1:23){
 plot(all_stat[i,],sqrt(all_sim[2,]),pch=19,main=rownames(all_stat)[i], log="xy", col=c(4,3,2,1)[as.factor(si$taxonomic.class)])
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}

# plot(all_sim[2,],all_m,pch=19, log="xy", main="Visscher & Goddard 2015")
plot(sqrt(all_sim[2,]),sqrt(all_sv),pch=19, log="xy", main="Visscher & Goddard 2015", col=c(4,3,2,1)[as.factor(si$taxonomic.class)])
abline(0,1)



###############
## MDE
###############


  par(mfrow=c(1,2), mar=c(5,5,1,1))
plot(all_stat[2,],sqrt(all_sim[2,]),pch=19, ylab="Standard error when h2=0.3",xlab="Number of maternities", xlim=c(0,11000))
plot(all_stat[2,],all_mde[3,],pch=19, ylab="Minimum detectable effect",xlab="Number of maternities", xlim=c(0,11000))
