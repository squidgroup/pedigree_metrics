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

taz_devil <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/taz_devil_grm.csv"), sep=",", header=TRUE, row.names=1))

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

soays_grm <- as.matrix(read.table(paste0(wd,"/Data/Raw/grms/soays_GRM.txt"),header=TRUE))
class(soays)
soays_lambda <- eigen(soays_grm, symmetric=TRUE,only.values=TRUE)$values

soays_se <-sqrt(lambda_to_samp_var(soays_lambda,h2=0.3))

soays[1:10,1:10]


####
load(paste0(wd,"/Data/Raw/grms/roe_deer_grm.Rdata"))

roe_grm<-b
roe_lambda <- eigen(roe_grm, symmetric=TRUE,only.values=TRUE)$values
roe_se <-sqrt(lambda_to_samp_var(roe_lambda,h2=0.3))


####

save(roe_lambda,soays_lambda,koala_lambda,dragon_lambda,taz_devil_lambda,s_lambda,file=paste0(wd,"/Data/Processed/GRM_lambdas.Rdata"))

####

setEPS()
pdf(paste0(wd,"Figures/Fig_grm.pdf"), height=8, width=8)
{

plot(all_sv~n, log="xy", ylim=range(c(all_sv,s_se)), col=scales::alpha(1,0.5), pch=19, ylab="Estimated Standard Error")
points(td_se~length(taz_devil_lambda),pch=21, cex=1.5, col="red",bg="red")
points(s_se~length(s_lambda),pch=22, cex=1.5, col="blue",bg="blue")
points(dragon_se~length(dragon_lambda),pch=23, cex=1.5, col="green",bg="green")
points(koala_se~length(koala_lambda),pch=24, cex=1.5, col="orange",bg="orange")

points(roe_se~length(roe_lambda),pch=25, cex=1.5, col="purple",bg="purple")

points(soays_se~length(soays_lambda),pch=21, cex=1.5, col="pink",bg="pink")

legend(
	"topright",
	c("pedigrees","Taz. Devil","Stickleback","Dragon","Koala","Roe Deer","Soay Sheep"),
	pch=c(19,21:25,21),
	pt.bg=c(scales::alpha(1,0.5),"red","blue","green","orange","purple","pink"),
	col=c(scales::alpha(1,0.5),"red","blue","green","orange","purple","pink")

	)

}
dev.off()

