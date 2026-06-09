rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))



####
# s_chol <- chol(stickleback)
# chol(taz_devil)
# chol(dragon_grm)
# chol(koala_grm)
# chol(roe_grm)
# MASS::mvrnorm(n=1, mu=rep(0,555), Sigma=s_chol, tol=1e-6)

# stickleback <- read.table(paste0(wd,"/Data/Raw/grms/stickleback_GRM.txt"),header=TRUE)
stickleback <- as.matrix(read.csv(paste0(wd,"/Data/Raw/grms/stickleback_GRM.csv"),row.names=1))

stickleback[upper.tri(stickleback)]<-t(stickleback)[upper.tri(stickleback)]
# chol(as(stickleback,"dgCMatrix"))
# stickleback[1:10,1:10]

s_lambda <- eigen(stickleback, symmetric=TRUE,only.values=TRUE)$values
s_se <- sqrt(lambda_to_samp_var(s_lambda,h2=0.3))



####

taz_devil <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))

# nrow(taz_devil)
# taz_devil[1:10,1:10]

taz_devil_lambda <- eigen(taz_devil, symmetric=TRUE,only.values=TRUE)$values
td_se <-sqrt(lambda_to_samp_var(taz_devil_lambda,h2=0.3))


####

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)

dragon_grm<-matrix(dragon[,3],nrow=length(unique(dragon[,1])),ncol=length(unique(dragon[,1])))
# hist(dragon_grm, breaks=100)
# dragon_grm[1:10,1:10]


dragon_lambda <- eigen(dragon_grm, symmetric=TRUE,only.values=TRUE)$values
dragon_se <- sqrt(lambda_to_samp_var(dragon_lambda,h2=0.3))


####

koala <- read.csv(paste0(wd,"/Data/Raw/grms/koala_GRM_2021.txt"))
koala2 <- sapply(koala[,1], function(x){
	y<-strsplit(x," ")[[1]]
	as.numeric(y[length(y)])
	})
# hist(koala2, breaks=100)

koala_grm<-matrix(koala2,nrow=sqrt(length(koala2)),ncol=sqrt(length(koala2)))

# koala_grm[1:10,1:10]

koala_lambda <- eigen(koala_grm, symmetric=TRUE,only.values=TRUE)$values
koala_se <- sqrt(lambda_to_samp_var(koala_lambda,h2=0.3))

####

soays_grm <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
# class(soays)
soays_lambda <- eigen(soays_grm, symmetric=TRUE,only.values=TRUE)$values

soays_se <-sqrt(lambda_to_samp_var(soays_lambda,h2=0.3))

# soays[1:10,1:10]


####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm<-b
# roe_grm[1:10,1:10]


roe_lambda <- eigen(roe_grm, symmetric=TRUE,only.values=TRUE)$values
roe_se <-sqrt(lambda_to_samp_var(roe_lambda,h2=0.3))


####

save(roe_lambda,soays_lambda,koala_lambda,dragon_lambda,taz_devil_lambda,s_lambda,file=paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))

grm_dat <- rbind(
	roe = c(roe_se, length(roe_lambda)),
soays = c(soays_se, length(soays_lambda)),
koala = c(koala_se, length(koala_lambda)),
dragon = c(dragon_se, length(dragon_lambda)),
taz_devil = c(td_se, length(taz_devil_lambda)),
stickleback = c(s_se, length(s_lambda))
)
colnames(grm_dat) <- c("se","n")

save(grm_dat,file=paste0(wd,"/Data/Processed/grm_dat.Rdata"))

####

