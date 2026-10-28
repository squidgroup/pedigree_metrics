rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

library(gremlin)
library(Matrix)

convert_grm <- function(grm){
	new_grm <- Matrix::nearPD(grm,base.matrix = TRUE,ensureSymmetry=TRUE)$mat
	colnames(new_grm)<-rownames(new_grm)<-colnames(grm)
	new_grm
}


####
s_grm_raw <- as.matrix(read.csv(paste0(wd,"/Data/Raw/grms/stickleback_GRM.csv"),row.names=1))

s_grm_raw[upper.tri(s_grm_raw)]<-t(s_grm_raw)[upper.tri(s_grm_raw)]
# s_grm_raw[1:10,1:10]

s_grm<-convert_grm(s_grm_raw)

set.seed(202606111)
s_lambda <- eigen(s_grm, symmetric=TRUE,only.values=TRUE)$values
s_se <- sqrt(lambda_to_samp_var(s_lambda,h2=0.3))
s_sim <- samp_var_sim_GRM(s_grm, h2=0.3, n_sim=1000)

stickleback_dat = c(se=s_se, n=length(s_lambda),s_sim)

save(stickleback_dat,file=paste0(wd,"/Data/Processed/grm_stickleback_dat.Rdata"))


####

td_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))
# taz_devil[1:10,1:10]

td_grm<-convert_grm(td_grm_raw)

set.seed(202606112)
td_lambda <- eigen(td_grm, symmetric=TRUE,only.values=TRUE)$values
td_se <-sqrt(lambda_to_samp_var(td_lambda,h2=0.3))
td_sim <- samp_var_sim_GRM(td_grm, h2=0.3, n_sim=1000)

td_dat = c(se=td_se, n=length(td_lambda),td_sim)

save(td_dat,file=paste0(wd,"/Data/Processed/grm_td_dat.Rdata"))

####

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)

dragon_grm_raw<-matrix(dragon[,3],nrow=sqrt(nrow(dragon)),ncol=sqrt(nrow(dragon)))
colnames(dragon_grm_raw)<-rownames(dragon_grm_raw)<-1:nrow(dragon_grm_raw)

dragon_grm <- convert_grm(dragon_grm_raw)
# dragon_grm[1:10,1:10]

set.seed(202606113)
dragon_lambda <- eigen(dragon_grm, symmetric=TRUE,only.values=TRUE)$values
dragon_se <- sqrt(lambda_to_samp_var(dragon_lambda,h2=0.3))
dragon_sim <- samp_var_sim_GRM(dragon_grm, h2=0.3, n_sim=1000)

dragon_dat = c(se=dragon_se, n=length(dragon_lambda),dragon_sim)

save(dragon_dat,file=paste0(wd,"/Data/Processed/grm_dragon_dat.Rdata"))

####

koala <- read.csv(paste0(wd,"/Data/Raw/grms/koala_GRM_2021.txt"))
koala2 <- sapply(koala[,1], function(x){
	y<-strsplit(x," ")[[1]]
	as.numeric(y[length(y)])
	})

koala_grm_raw<-matrix(koala2,nrow=sqrt(length(koala2)),ncol=sqrt(length(koala2)))
colnames(koala_grm_raw)<-rownames(koala_grm_raw)<-1:nrow(koala_grm_raw)
koala_grm <- convert_grm(koala_grm_raw)
# koala_grm[1:10,1:10]

set.seed(202606114)
koala_lambda <- eigen(koala_grm, symmetric=TRUE,only.values=TRUE)$values
koala_se <- sqrt(lambda_to_samp_var(koala_lambda,h2=0.3))
koala_sim <- samp_var_sim_GRM(koala_grm, h2=0.3, n_sim=1000)

koala_dat = c(se=koala_se, n=length(koala_lambda),koala_sim)

save(koala_dat,file=paste0(wd,"/Data/Processed/grm_koala_dat.Rdata"))

####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm_raw<-as.matrix(b)
roe_grm <- convert_grm(roe_grm_raw)

# roe_grm[1:10,1:10]

set.seed(202606115)
roe_lambda <- eigen(roe_grm, symmetric=TRUE,only.values=TRUE)$values
roe_se <-sqrt(lambda_to_samp_var(roe_lambda,h2=0.3))
roe_sim <- samp_var_sim_GRM(roe_grm, h2=0.3, n_sim=1000)

roe_dat = c(se=roe_se, n=length(roe_lambda),roe_sim)

save(roe_dat,file=paste0(wd,"/Data/Processed/grm_roe_dat.Rdata"))

####

soays_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
# class(soays)
soays_grm <- convert_grm(soays_grm_raw)
# soays_grm[1:10,1:10]
# save(soays_grm,file=paste0(wd,"/Data/Raw/grms/soays_grm.Rdata"))

## this one takes days to run!! 
set.seed(202606116)
soays_lambda <- eigen(soays_grm, symmetric=TRUE,only.values=TRUE)$values
soays_se <-sqrt(lambda_to_samp_var(soays_lambda,h2=0.3))
soays_sim <- samp_var_sim_GRM(soays_grm, h2=0.3, n_sim=1000)

soays_dat = c(se=soays_se, n=length(soays_lambda),soays_sim)

save(soays_dat,file=paste0(wd,"/Data/Processed/grm_soays_dat.Rdata"))
#start 10:35 saturday
#30 by 13:55 saturday
####


