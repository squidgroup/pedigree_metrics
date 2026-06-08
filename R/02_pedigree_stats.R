rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))
library(gremlin)
library(pedtricks)

####
# Read in Data
####

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
type <- sapply(strsplit(files,"\\."),function(x) x[2])
peds<-lapply(1:length(files), function(i){
  ped <- read.table(paste0(wd,"Data/Raw/pedigrees/",files[i]),header=TRUE,sep=if(type[i]=="txt"){""}else{","}, na.strings=c("NA",0))
  ped2 <- data.frame(animal=ped[,1],dam=ped[,"dam"],sire=ped[,"sire"])
  pedtricks::fix_ped(ped2)
})
names(peds) <- files

# pedigree lengths
range(sapply(peds,nrow))




####
# Calculate stats for each pedigrees
####

# summary(ped_stats(peds[[1]], includeA=FALSE))
# ped_stat(peds[[1]])
all_stat <- do.call(cbind,parallel::mclapply(peds,function(ped) summary(ped_stats(ped, includeA=FALSE)), mc.cores=6))
colnames(all_stat) <- files
save(all_stat,file=paste0(wd,"Data/Processed/ped_stats.Rdata"))


## correlation between number of materities and number of individuals with parents
mat_vs_parent <- do.call(rbind,parallel::mclapply(peds,function(ped) {
  c(
    maternities = sum(!is.na(ped[,2])),
    any_parent=sum(!(is.na(ped[,2])&is.na(ped[,3]) ) )
    )
}, mc.cores=6))
cor(log(mat_vs_parent))
plot(mat_vs_parent, log="xy")


####
# Calculate Visscher and Goddard 2015 sampling variance for each pedigree
####
 
all_sv <- do.call(rbind,parallel::mclapply(peds, function(ped){
  sv <- pedtricks::h2_samp_var(ped,c(0.1,0.3,0.5))
  cat("I ")
  sv
  } ,mc.cores=6)
)
save(all_sv,file=paste0(wd,"Data/Processed/ped_h2_se.Rdata"))



####
# Calculate Visscher et al. 2014 sampling variance for each pedigrees
####

samp_var_h2_simple <- function(ped){
  sp<-Matrix::summary(nadiv::makeA(ped[,1:3]))
  sp2<-sp[sp[,1]!=sp[,2],]
  total_links <- (nrow(ped) * (nrow(ped) - 1) / 2 )
  var_a <-  var(c(sp2[,3], rep(0,total_links-nrow(sp2))))
  N <- nrow(ped)
  2/(N^2 * var_a)
}

all_h2_se_simple <- sapply(peds,function(ped){
  samp_var_h2_simple(ped)
})
save(all_h2_se_simple,file=paste0(wd,"Data/Processed/ped_h2_se_simple.Rdata"))



####
# Calculate sampling variance for each pedigrees via simulation across 3 heritabilities
####
set.seed(20241105)
all_sim <- sapply(peds,samp_var_sim,h2=0.3,n_sim=1000)

save(all_sim,file=paste0(wd,"Data/Processed/ped_sims.Rdata"))


set.seed(20251119)
all_sim_low <- sapply(peds,samp_var_sim,h2=0.1,n_sim=1000)

save(all_sim_low,file=paste0(wd,"Data/Processed/ped_sims_low.Rdata"))


set.seed(20251120)
all_sim_high <- sapply(peds,samp_var_sim,h2=0.5,n_sim=1000)

save(all_sim_high,file=paste0(wd,"Data/Processed/ped_sims_high.Rdata"))




####
# Calculate minimum detectable effect for each pedigrees via simulation
####
set.seed(20241106)
all_mde <- sapply(peds,MDE_sim,1000)
# MDE_sim(peds[[1]],10)

save(all_mde,file=paste0(wd,"Data/Processed/ped_mde.Rdata"))



####
# Calculate bias from unmodelled maternal effects for each pedigrees via simulation
####
set.seed(20241107)
all_mbias_sims <- sapply(peds,samp_var_sim,h2=0.1,m2=0.2, n_sim=500)
save(all_mbias_sims,file=paste0(wd,"Data/Processed/ped_mbias.Rdata"))

