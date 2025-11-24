rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))
load(paste0(wd,"Data/Processed/ped_h2_se_simple.Rdata"))


lambdas <- lambdas[names(lambdas) %in% colnames(all_stat)]
all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.3)
length(all_sv)
all_sim <- all_sim[2,colnames(all_sim) %in% names(all_sv)]

all_h2_se_simple <- all_h2_se_simple[names(all_h2_se_simple) %in% names(all_sv)]


r2 <- summary(lm(log(sqrt(all_sim))~log(sqrt(all_sv))))$adj.r.squared
r2_simple <- summary(lm(log(sqrt(all_sim))~log(sqrt(all_h2_se_simple))))$adj.r.squared



# plot(all_sim[2,],all_m,pch=19, log="xy", main="Visscher & Goddard 2015")

par(mfrow=c(1,2))
plot(sqrt(all_sv),sqrt(all_sim),pch=19, log="xy", xlab="V-G SE", ylab="Expected SE", col= scales::alpha(1,0.5), main=paste0("r2=",round(r2,3)))
abline(0,1)

plot(sqrt(all_h2_se_simple),sqrt(all_sim),pch=19, log="xy", xlab="simple SE", ylab="Expected SE",ylim=range(c(sqrt(all_h2_se_simple),sqrt(all_sim))), col= scales::alpha(1,0.5), main=paste0("r2=",round(r2_simple,3)))
abline(0,1)
