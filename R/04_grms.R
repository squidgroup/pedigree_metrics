rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

library(gremlin)
library(Matrix)

convert_grm <- function(grm) Matrix::nearPD(grm,base.matrix = TRUE,ensureSymmetry=TRUE)$mat


####
s_grm_raw <- as.matrix(read.csv(paste0(wd,"/Data/Raw/grms/stickleback_GRM.csv"),row.names=1))

s_grm_raw[upper.tri(s_grm_raw)]<-t(s_grm_raw)[upper.tri(s_grm_raw)]
# s_grm_raw[1:10,1:10]

s_grm<-convert_grm(s_grm_raw)


####

td_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))
# taz_devil[1:10,1:10]

td_grm<-convert_grm(taz_devil)



####

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)

dragon_grm_raw<-matrix(dragon[,3],nrow=sqrt(nrow(dragon)),ncol=sqrt(nrow(dragon)))
colnames(dragon_grm_raw)<-rownames(dragon_grm_raw)<-1:nrow(dragon_grm_raw)

dragon_grm <- convert_grm(dragon_grm_raw)
# dragon_grm[1:10,1:10]


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

####

soays_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
# class(soays)
soays_grm <- convert_grm(soays_grm_raw)
# soays_grm[1:10,1:10]
# save(soays_grm,file=paste0(wd,"/Data/Raw/grms/soays_grm.Rdata"))


####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm_raw<-as.matrix(b)
roe_grm <- convert_grm(roe_grm_raw)

# roe_grm[1:10,1:10]

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


## this one takes days to run!! 
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

grm_dat2 <- rbind(
	roe = c(roe_se, length(roe_lambda)),
soays = c(soays_se, length(soays_lambda)),
koala = c(koala_se, length(koala_lambda)),
dragon = c(dragon_se, length(dragon_lambda)),
taz_devil = c(td_se, length(td_lambda)),
stickleback = c(s_se, length(s_lambda))
)
colnames(grm_dat2)[1:2] <- c("se","n")


grm_dat[,1:2]<-grm_dat2

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


grm_dat <- rbind(grm_dat,koala2 =c(koala_se2, length(koala_lambda2),grm_dat["koala",3:4]) )
####





s_lambda_raw <- eigen(s_grm_raw, symmetric=TRUE,only.values=TRUE)$values
s_se_raw <- sqrt(lambda_to_samp_var(s_lambda_raw,h2=0.3))

td_lambda_raw <- eigen(td_grm_raw, symmetric=TRUE,only.values=TRUE)$values
td_se_raw <- sqrt(lambda_to_samp_var(td_lambda_raw,h2=0.3))

dragon_lambda_raw <- eigen(dragon_grm_raw, symmetric=TRUE,only.values=TRUE)$values
dragon_se_raw <- sqrt(lambda_to_samp_var(dragon_lambda_raw,h2=0.3))

roe_lambda_raw <- eigen(roe_grm_raw, symmetric=TRUE,only.values=TRUE)$values
roe_se_raw <- sqrt(lambda_to_samp_var(roe_lambda_raw,h2=0.3))

koala_lambda_raw <- eigen(koala_grm_raw, symmetric=TRUE,only.values=TRUE)$values
koala_se_raw <- sqrt(lambda_to_samp_var(koala_lambda_raw,h2=0.3))

soays_lambda_raw <- eigen(soays_grm_raw, symmetric=TRUE,only.values=TRUE)$values
soays_se_raw <- sqrt(lambda_to_samp_var(soays_lambda_raw,h2=0.3))


se_h2_comp <- cbind(
	adjusted=c(s_se,td_se,dragon_se,koala_se,roe_se,soays_se_raw),
	raw=c(s_se_raw,td_se_raw,dragon_se_raw,koala_se_raw,roe_se_raw,soays_se_raw))
save(se_h2_comp,file=paste0(wd,"/Data/Processed/grm_dat_compare.Rdata"))





