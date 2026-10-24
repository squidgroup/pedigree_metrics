rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

####

load(file=paste0(wd,"/Data/Processed/grm_roe_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_soay_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_koala_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_dragon_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_td_dat.Rdata"))
load(file=paste0(wd,"/Data/Processed/grm_stickleback_dat.Rdata"))



grm_dat <- rbind(
	roe= roe_dat,
	soays = soays_dat,
	koala = koala_dat,
	dragon = dragon_dat,
	taz_devil = taz_devil_dat,
	stickleback = stickleback_dat

	)


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





### Supplementary figure comparing SE approximation on raw GRMs, compared with ones adjusted to be positive definite

setEPS()
pdf(paste0(wd,"Figures/FigSM_grm_compare.pdf"), height=5, width=5)
{
plot(se_h2_comp, log="xy", cex=1.5, pch=19,col=scales::alpha(1,0.5), xlab=expression(Adjusted~italic(SE[h^2]^{(2)} )), ylab=expression(Raw~italic(SE[h^2]^{(2)} )))
abline(0,1)

}
dev.off()
