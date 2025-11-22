rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

# save(all_sim_low,file=paste0(wd,"Data/Processed/ped_sims_low.Rdata"))

load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims_low.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))

all_sim <- all_sim[,colnames(all_sim) %in% names(lambdas)]
all_sim_low <- all_sim_low[,colnames(all_sim_low) %in% names(lambdas)]

all_sv <- sapply(lambdas,lambda_to_samp_var,h2=c(0.1,0.3))

par(mfrow=c(1,2))
plot(sqrt(all_sv[1,]),sqrt(all_sim_low[2,]),log="xy")
abline(0,1)
plot(sqrt(all_sv[2,]),sqrt(all_sim[2,]),log="xy")
abline(0,1)



length(all_sim_low[2,])
length(all_sv)

summary(lm(sqrt(all_sim_low[2,])~sqrt(all_sv[1,])))
summary(lm(sqrt(all_sim[2,])~sqrt(all_sv[2,])))


all_sim <- all_sim[,colnames(all_sim) %in% names(lambdas)]

all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.1)
