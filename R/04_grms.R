rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


####

load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))

lambdas <- lambdas[!names(lambdas) %in% c("D33.csv","D38.txt")]
all_sv <- sqrt(sapply(lambdas,lambda_to_samp_var,h2=0.3))
n <- sapply(lambdas,length)


####

stickleback <- read.table(paste0(wd,"/Data/Raw/grms/stickleback_GRM.txt"),header=TRUE)

tail(stickleback)
table(stickleback[,2])
length(unique(stickleback[,1]))
s_grm<-matrix(stickleback[,3],nrow=length(unique(stickleback[,1])),ncol=length(unique(stickleback[,1])))
diag(s_grm)<-1
hist(s_grm, breaks=100)

hist(stickleback[stickleback[,3]>0.1,3])

s_lambda <- eigen(s_grm, symmetric=TRUE,only.values=TRUE)$values
s_se <- sqrt(lambda_to_samp_var(s_lambda,h2=0.3))



####

taz_devil <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_deveil_grm.csv"), sep=",", header=TRUE, row.names=1))

nrow(taz_devil)

taz_devil_lambda <- eigen(taz_devil, symmetric=TRUE,only.values=TRUE)$values
td_se <-sqrt(lambda_to_samp_var(taz_devil_lambda,h2=0.3))


####

dragon <- read.table(paste0(wd,"/Data/Raw/grms/dragon_GRM_2020.txt"),header=TRUE)

dragon_grm<-matrix(dragon[,3],nrow=length(unique(dragon[,1])),ncol=length(unique(dragon[,1])))
hist(dragon_grm, breaks=100)


dragon_lambda <- eigen(dragon_grm, symmetric=TRUE,only.values=TRUE)$values
dragon_se <- sqrt(lambda_to_samp_var(dragon_lambda,h2=0.3))


####

koala <- read.csv(paste0(wd,"/Data/Raw/grms/koala_GRM_2021.txt"))
koala2 <- sapply(koala[,1], function(x){
	y<-strsplit(x," ")[[1]]
	as.numeric(y[length(y)])
	})
hist(koala2, breaks=100)

koala_grm<-matrix(koala2,nrow=sqrt(length(koala2)),ncol=sqrt(length(koala2)))
koala_lambda <- eigen(koala_grm, symmetric=TRUE,only.values=TRUE)$values
koala_se <- sqrt(lambda_to_samp_var(koala_lambda,h2=0.3))

####

setEPS()
pdf(paste0(wd,"Figures/Fig_grm.pdf"), height=8, width=8)
{

plot(all_sv~n, log="xy", ylim=range(c(all_sv,s_se)), col=scales::alpha(1,0.5), pch=19, ylab="Estimated Standard Error")
points(td_se~length(taz_devil_lambda),pch=19, col="red")
points(s_se~length(s_lambda),pch=19, col="blue")
points(dragon_se~length(dragon_lambda),pch=19, col="green")
points(koala_se~length(koala_lambda),pch=19, col="orange")

}
dev.off()

