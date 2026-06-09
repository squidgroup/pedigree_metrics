rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

load(paste0(wd,"Data/Processed/ped_missing.Rdata"))

total_n<-sapply(peds,nrow)
pheno_n<-sapply(peds,function(x)sum(!is.na(x$dam)))

missing_n <- (total_n - pheno_n)/total_n
range(missing_n)
# hist(missing_n)



setEPS()
pdf(paste0(wd,"Figures/Fig_missing.pdf"), height=5, width=10)
  
{
  par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.5)
  plot(total_n, sqrt(all_missing["var_all",]), log="xy", pch=19, col=scales::alpha(1,0.4), xlim=range(c(pheno_n,total_n)), ylim=range(sqrt(all_missing),na.rm=TRUE), xlab="Records", ylab="Expected standard error")
  points(pheno_n, sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("green",0.6))
  points(total_n, sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("blue",0.5))
  mtext("a)", side=3, adj=-0.2,padj=1, cex=2)

  plot(sqrt(all_sv_missing[,1]), sqrt(all_missing["var_all",]), pch=19,col=scales::alpha(1,0.4), xlim=range(sqrt(all_sv_missing)), ylim=range(sqrt(all_missing),na.rm=TRUE), log="xy", xlab=expression(italic(SE[h^2]^{(1)} )), ylab="Expected standard error")
  points(sqrt(all_sv_missing[,2]), sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("green",0.6))
  points(sqrt(all_sv_missing[,1]), sqrt(all_missing["var_missing",]), pch=19, col=scales::alpha("blue",0.5))
  abline(0,1)
  mtext("b)", side=3, adj=-0.2,padj=1, cex=2)


}
dev.off()
