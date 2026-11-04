rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


####

load(file=paste0(wd,"/Data/Processed/grm_roe_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_soays_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_koala_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_dragon_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_td_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_stickleback_dat.Rdata"))


####
s_grm_raw <- as.matrix(read.csv(paste0(wd,"/Data/Raw/grms/stickleback_GRM.csv"),row.names=1))

s_grm_raw[upper.tri(s_grm_raw)]<-t(s_grm_raw)[upper.tri(s_grm_raw)]
# s_grm_raw[1:10,1:10]

s_lambda_raw <- eigen(s_grm_raw, symmetric=TRUE,only.values=TRUE)$values
s_se_raw <- sqrt(lambda_to_samp_var(s_lambda_raw,h2=0.3))
####

td_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))

td_lambda_raw <- eigen(td_grm_raw, symmetric=TRUE,only.values=TRUE)$values
td_se_raw <- sqrt(lambda_to_samp_var(td_lambda_raw,h2=0.3))


####

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)

dragon_grm_raw<-matrix(dragon[,3],nrow=sqrt(nrow(dragon)),ncol=sqrt(nrow(dragon)))
colnames(dragon_grm_raw)<-rownames(dragon_grm_raw)<-1:nrow(dragon_grm_raw)

dragon_lambda_raw <- eigen(dragon_grm_raw, symmetric=TRUE,only.values=TRUE)$values
dragon_se_raw <- sqrt(lambda_to_samp_var(dragon_lambda_raw,h2=0.3))


####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm_raw<-as.matrix(b)

roe_lambda_raw <- eigen(roe_grm_raw, symmetric=TRUE,only.values=TRUE)$values
roe_se_raw <- sqrt(lambda_to_samp_var(roe_lambda_raw,h2=0.3))


####

koala <- read.csv(paste0(wd,"/Data/Raw/grms/koala_GRM_2021.txt"))
koala2 <- sapply(koala[,1], function(x){
	y<-strsplit(x," ")[[1]]
	as.numeric(y[length(y)])
	})

koala_grm_raw<-matrix(koala2,nrow=sqrt(length(koala2)),ncol=sqrt(length(koala2)))
colnames(koala_grm_raw)<-rownames(koala_grm_raw)<-1:nrow(koala_grm_raw)

koala_lambda_raw <- eigen(koala_grm_raw, symmetric=TRUE,only.values=TRUE)$values
koala_se_raw <- sqrt(lambda_to_samp_var(koala_lambda_raw,h2=0.3))


####

soays_grm_raw <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
soays_lambda_raw <- eigen(soays_grm_raw, symmetric=TRUE,only.values=TRUE)$values
soays_se_raw <- sqrt(lambda_to_samp_var(soays_lambda_raw,h2=0.3))


se_h2_comp <- cbind(
	adjusted=c(
	roe_dat["se"],
	soays_dat["se"],
	koala_dat["se"],
	dragon_dat["se"],
	td_dat["se"],
	stickleback_dat["se"]
	),
	raw=c(roe_se_raw,
		soays_se_raw,
		koala_se_raw,
		dragon_se_raw,
		td_se_raw,
		s_se_raw)
)

save(se_h2_comp,file=paste0(wd,"/Data/Processed/grm_dat_compare.Rdata"))





### Supplementary figure comparing SE approximation on raw GRMs, compared with ones adjusted to be positive definite

setEPS()
pdf(paste0(wd,"Figures/FigSM_grm_compare.pdf"), height=5, width=5)
{
	par(mar=c(5,5,1,1))
plot(se_h2_comp, log="xy", cex=1.5, pch=19,col=scales::alpha(1,0.5), xlab=expression(Adjusted~italic(SE[hat(h^2)] ^{"  "~(2)} )), ylab=expression(Raw~italic(SE[hat(h^2)] ^{"  "~(2)} )))
abline(0,1)

}
dev.off()
