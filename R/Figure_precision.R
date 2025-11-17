rm(list=ls())


wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_stats.Rdata"))
load(paste0(wd,"Data/Processed/ped_sims.Rdata"))
load(paste0(wd,"Data/Processed/ped_lambdas.Rdata"))
load(paste0(wd,"Data/Processed/ped_mde.Rdata"))

all_stat <- all_stat[,colnames(all_stat) %in% colnames(all_sim)]

all_sim <- all_sim[,colnames(all_sim) %in% colnames(all_stat)]

lambdas <- lambdas[names(lambdas) %in% colnames(all_stat)]
all_sv <- sapply(lambdas,lambda_to_samp_var,h2=0.3)

# all_mde <- all_mde[,!colnames(all_mde) %in% c("D29.csv","D33.csv","D38.txt","D37.csv")]

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
pdf(paste0(wd,"Figures/Fig_precision.pdf"), height=8, width=14)

{
par(mfrow=c(4,6),mar=c(3.5,3,2,1),mgp=c(2,0.5,0))#mar=c(2,2,3,1))
for(i in 1:23){

 r2 <- summary(lm(log(sqrt(all_sim[2,]))~log(all_stat2[i,])))$adj.r.squared

 plot(all_stat2[i,],sqrt(all_sim[2,]),pch=19, col= scales::alpha(1,0.5),xlab=rownames(all_stat)[i], ylab="Expected SE", log="xy", main=paste0("r2=",round(r2,3)))# , col=c(4,3,2,1)[as.factor(si$taxonomic.class)]
 # text(25, min(log(all_stat[i,]))+(max(log(all_stat[i,]))-min(log(all_stat[i,])))/10,rownames(all_stat)[i])
}

r2 <- summary(lm(log(sqrt(all_sim[2,]))~log(sqrt(all_sv))))$adj.r.squared

# plot(all_sim[2,],all_m,pch=19, log="xy", main="Visscher & Goddard 2015")
plot(sqrt(all_sv),sqrt(all_sim[2,]),pch=19, log="xy", xlab="V-G SE", ylab="Expected SE", col= scales::alpha(1,0.5), main=paste0("r2=",round(r2,3)))
abline(0,1)

}
dev.off()

# summary(lm(log(sqrt(all_sim[2,]))~log(sqrt(all_sv))))$adj.r.squared


###############
## MDE
###############

setEPS()
pdf(paste0(wd,"Figures/Fig_precision_mat.pdf"), height=5, width=7)
{
  par(mfrow=c(1,1), mar=c(5,5,1,1))

plot(all_stat[2,],sqrt(all_sim[2,]),pch=19, ylab="Expected SE when h2=0.3",xlab="Number of maternities", xlim=c(0,11000), col= scales::alpha(1,0.5), ylim=c(0,0.2))

x<-seq(min(log(all_stat2["maternities",])),max(log(all_stat2["maternities",])),0.1)
preds <- as.vector(exp(model.matrix(~ x) %*% coef(lm(log(sqrt(all_sim[2,]))~log(all_stat2["maternities",])))))

lines(preds~exp(x), col="red")

}
dev.off()


# plot(all_stat[2,],all_mde[3,],pch=19, ylab="Minimum detectable effect",xlab="Number of maternities", xlim=c(0,11000))
