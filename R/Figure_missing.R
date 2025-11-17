rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_missing.Rdata"))


total_n<-sapply(peds,nrow)
pheno_n<-sapply(peds,function(x)sum(!is.na(x$dam)))

missing_n <- (total_n - pheno_n)/total_n
range(missing_n)
hist(missing_n)

all_sv <- sapply(lambdas,function(x) c(lambda_to_samp_var(x[["lambda_full"]],h2=0.3),lambda_to_samp_var(x[["lambda_missing"]],h2=0.3)))


setEPS()
pdf(paste0(wd,"Figures/Fig_missing.pdf"), height=4, width=10)

{
  par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(total_n, sqrt(all_missing["var_all",]), log="xy", pch=19, col=scales::alpha(1,0.4), xlim=range(c(pheno_n,total_n)), ylim=range(sqrt(all_missing),na.rm=TRUE), xlab="Records", ylab="Expected standard error")
points(pheno_n, sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("green",0.6))

points(total_n, sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("blue",0.5))


plot(sqrt(all_sv[1,]), sqrt(all_missing["var_all",]), pch=19,col=scales::alpha(1,0.4), xlim=range(sqrt(all_sv)), ylim=range(sqrt(all_missing),na.rm=TRUE), log="xy", xlab="V-G standard error", ylab="Expected standard error")
points(sqrt(all_sv[2,]), sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("green",0.6))

points(sqrt(all_sv[1,]), sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("blue",0.5))


abline(0,1)
}
dev.off()
