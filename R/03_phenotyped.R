rm(list=ls())
library(gremlin)

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
type <- sapply(strsplit(files,"\\."),function(x) x[2])

## create reduced pedigrees pruned assuming that only individuals with a mother have a phenotype.
peds<-parallel::mclapply(1:length(files), function(i){
  ped <- read.table(paste0(wd,"Data/Raw/pedigrees/",files[i]),header=TRUE,sep=if(type[i]=="txt"){""}else{","}, na.strings=c("NA",0))
  ped2 <- data.frame(animal=ped[,1],dam=ped[,"dam"],sire=ped[,"sire"])
  ped3 <- pedtricks::fix_ped(ped2)
  MCMCglmm::prunePed(ped3, keep=ped3[!is.na(ped3$dam),1], make.base=TRUE)
}, mc.cores=6)
names(peds) <- files


all_sv_missing <- do.call(rbind,parallel::mclapply(peds,function(ped){
  subset_ids <- ped[!is.na(ped$dam),1]
  A<-nadiv::makeA(ped)
  A2<-A[subset_ids,subset_ids]

  lambda1 <- eigen(A, symmetric=TRUE,only.values=TRUE)$values
  lambda2 <- eigen(A2, symmetric=TRUE,only.values=TRUE)$values
  cat("I ")
  c(full=lambda_to_samp_var(lambda1, h2=0.3),
  missing=lambda_to_samp_var(lambda2, h2=0.3))

}, mc.cores=6))




## take pruned pedigree.
## simulations h2=0.3
## model all individuals 
## model only individuals with mothers
## VG of 

# ped<-peds[[2]]
# h2=0.3
# n_sim=200
samp_var_sim_missing <- function(ped, h2, n_sim=200){
  
  data<-squidSim::get_population_data(
    squidSim::simulate_population(
      data_structure=ped,
      parameters=list(
        id=list(vcov=h2),
        residual=list(vcov=1-h2)
      ),
      pedigree=list(id=ped),
      n_pop=n_sim
    ),
  list=TRUE)

  assign("ped.ainv", nadiv::makeAinv(ped)$Ainv, envir = .GlobalEnv)

# x<-data[[8]]
  sims<-do.call(rbind,parallel::mclapply(data, function(x){

    x$id <- as.factor(x$id)

    data_missing <- x[!is.na(ped[,"dam"]),]

    mod_1 <- gremlin::gremlin(
            y~1
            , random= ~id
            , data= x
            , ginverse = list(id = ped.ainv)
            ,v=0
            )

    va_mod_1 <- deltaSE(list(h2 ~ V1 / (V1 + V2)), mod_1)["Estimate"]
    va_mod_2 <- tryCatch({
     
      mod_2 <- gremlin::gremlin(
        y~1
        , random= ~id
        , data= data_missing
        , ginverse = list(id = ped.ainv)
        ,v=0
        )
      deltaSE(list(h2 ~ V1 / (V1 + V2)), mod_1)["Estimate"]
      },
    error=function(e)0)
    
    c(va_mod_1,va_mod_2)
      
  },mc.cores=6))  
  rm("ped.ainv", envir = .GlobalEnv)
  cat("I ")
  c(var_all=var(sims[,1]),var_missing=var(sims[,2]))
}

# samp_var_sim_missing(peds[[2]],0.3, n_sim=1000)

all_missing <- sapply(peds,samp_var_sim_missing,0.3, n_sim=1000)
save(peds,all_missing,all_sv_missing,file=paste0(wd,"Data/Processed/ped_missing.Rdata"))
#save(peds,all_missing,file=paste0(wd,"Data/Processed/ped_missing2.Rdata"))

