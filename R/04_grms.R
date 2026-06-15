rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

library(gremlin)
library(Matrix)


####
stickleback <- as.matrix(read.csv(paste0(wd,"/Data/Raw/grms/stickleback_GRM.csv"),row.names=1))

stickleback[upper.tri(stickleback)]<-t(stickleback)[upper.tri(stickleback)]
# chol(as(stickleback,"dgCMatrix"))
# stickleback[1:10,1:10]

####

taz_devil <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))

# nrow(taz_devil)
# taz_devil[1:10,1:10]
isSymmetric(taz_devil)

td_grm<-forceSymmetric(taz_devil)

####
make_grm <- function(df){
	grm <- matrix(NA,ncol=length(unique(df[,1])),nrow=length(unique(df[,1])))
	rownames(grm)<-colnames(grm)<-unique(df[,1])
	for(i in 1:nrow(df)) grm[df[i,1],df[i,2]] <- df[i,3]
	grm
}

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)
dragon_grm <- make_grm(dragon)
# dragon_grm[1:10,1:10]

isSymmetric(dragon_grm)


####

koala <- read.csv(paste0(wd,"/Data/Raw/grms/koala_GRM_2021.txt"))
koala2 <- sapply(koala[,1], function(x){
	y<-strsplit(x," ")[[1]]
	as.numeric(y[length(y)])
	})
# hist(koala2, breaks=100)

koala_grm<-matrix(koala2,nrow=sqrt(length(koala2)),ncol=sqrt(length(koala2)))
colnames(koala_grm)<-rownames(koala_grm)<-1:nrow(koala_grm)
# koala_grm[1:10,1:10]
isSymmetric(koala_grm)

####
library(Matrix)
soays_grm <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
# class(soays)
soays_grm <- forceSymmetric(soays_grm)
# soays_grm[1:10,1:10]
# save(soays_grm,file=paste0(wd,"/Data/Raw/grms/soays_grm.Rdata"))


####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm<-as.matrix(b)
# roe_grm[1:10,1:10]

isSymmetric(roe_grm)

####


# save(roe_lambda,soays_lambda,koala_lambda,dragon_lambda,taz_devil_lambda,s_lambda,file=paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))

set.seed(202606111)
s_lambda <- eigen(stickleback, symmetric=TRUE,only.values=TRUE)$values
s_se <- sqrt(lambda_to_samp_var(s_lambda,h2=0.3))
s_sim <- samp_var_sim_GRM(stickleback, h2=0.3, n_sim=1000)

set.seed(202606112)
td_lambda <- eigen(td_grm, symmetric=TRUE,only.values=TRUE)$values
td_se <-sqrt(lambda_to_samp_var(td_lambda,h2=0.3))
td_sim <- samp_var_sim_GRM(td_grm, h2=0.3, n_sim=1000)

set.seed(202606113)
dragon_lambda <- eigen(dragon_grm, symmetric=TRUE,only.values=TRUE)$values
dragon_se <- sqrt(lambda_to_samp_var(dragon_lambda,h2=0.3))
dragon_sim <- samp_var_sim_GRM(dragon_grm, h2=0.3, n_sim=1000)

set.seed(202606114)
koala_lambda <- eigen(koala_grm, symmetric=TRUE,only.values=TRUE)$values
koala_se <- sqrt(lambda_to_samp_var(koala_lambda,h2=0.3))
koala_sim <- samp_var_sim_GRM(koala_grm, h2=0.3, n_sim=1000)

set.seed(202606115)
roe_lambda <- eigen(roe_grm, symmetric=TRUE,only.values=TRUE)$values
roe_se <-sqrt(lambda_to_samp_var(roe_lambda,h2=0.3))
roe_sim <- samp_var_sim_GRM(roe_grm, h2=0.3, n_sim=1000)


## this one takes ages to run!! 
set.seed(202606116)
soays_lambda <- eigen(soays_grm, symmetric=TRUE,only.values=TRUE)$values
soays_se <-sqrt(lambda_to_samp_var(soays_lambda,h2=0.3))
soays_sim <- samp_var_sim_GRM(soays_grm, h2=0.3, n_sim=1000)


# grm_sim <- rbind(
# 	roe=roe_sim,
# 	soays = soays_sim,
# 	koala = koala_sim,
# 	dragon = dragon_sim,
# 	taz_devil = taz_devil_sim,
# 	stickleback = stickleback_sim

# 	)

# save(grm_sim,file=paste0(wd,"/Data/Processed/grm_sims.Rdata"))

grm_dat <- rbind(
	roe = c(roe_se, length(roe_lambda),roe_sim),
soays = c(soays_se, length(soays_lambda),soays_sim),
koala = c(koala_se, length(koala_lambda),koala_sim),
dragon = c(dragon_se, length(dragon_lambda),dragon_sim),
taz_devil = c(td_se, length(td_lambda),td_sim),
stickleback = c(s_se, length(s_lambda),s_sim)
)
colnames(grm_dat)[1:2] <- c("se","n")

save(grm_dat,file=paste0(wd,"/Data/Processed/grm_dat.Rdata"))

####















