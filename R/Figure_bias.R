rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))
load(paste0(wd,"Data/Processed/ped_mde.Rdata"))
load(paste0(wd,"Data/Processed/ped_mbias.Rdata"))

all_stat <- all_stat[,colnames(all_stat) %in% colnames(all_sim)]

all_sim <- all_sim[,colnames(all_sim) %in% colnames(all_stat)]

lambdas <- lambdas[names(lambdas) %in% colnames(all_stat)]
all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.3)

all_mbias_sims <- all_mbias_sims[,colnames(all_mbias_sims) %in% colnames(all_stat)]



files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
id <- gsub("\\.|[a-z]","",files)

studies_info <- read.csv(paste0(wd,"Data/Raw/studies_info.csv"))
dryad_ped_info <- read.csv(paste0(wd,"Data/Raw/dryad_ped_info.csv"))
all_ped_info <- rbind(studies_info[,c("id","species_english","taxonomic.class")], dryad_ped_info[,c("id","species_english","taxonomic.class")])

si <- all_ped_info[match(id,all_ped_info$id),c("id","species_english","taxonomic.class")]


includes0<-apply(all_stat,1,function(x)any(x==0))
all_stat2<-all_stat
all_stat2[includes0,]<-all_stat2[includes0,]+0.5


setEPS()
pdf(paste0(wd,"Figures/Fig_bias.pdf"), height=8, width=14)
{

par(mfrow=c(4,6),mar=c(3.5,3,2,1),mgp=c(2,0.5,0))#mar=c(2,2,3,1))
for(i in 1:23){
  r2 <- summary(lm(log(all_mbias_sims[1,])~log(all_stat2[i,])))$adj.r.squared

  plot(all_stat2[i,],all_mbias_sims[1,],pch=19, col= scales::alpha(1,0.5),xlab=rownames(all_stat)[i], ylab="Bias in Va", log="xy", main=paste0("r2=",round(r2,3)))
  }


r2 <- summary(lm(log(all_mbias_sims[1,])~log(sqrt(all_sv))))$adj.r.squared

# plot(all_sim[2,],all_m,pch=19, log="xy", main="Visscher & Goddard 2015")
plot(sqrt(all_sv),all_mbias_sims[1,],pch=19, log="xy", xlab="V-G SE", ylab="Bias in Va", col= scales::alpha(1,0.5), main=paste0("r2=",round(r2,3)))

}
dev.off()
