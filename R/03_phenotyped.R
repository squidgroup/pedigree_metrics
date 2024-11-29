rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))

files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
type <- sapply(strsplit(files,"\\."),function(x) x[2])
peds<-parallel::mclapply(1:length(files), function(i){
  ped <- read.table(paste0(wd,"Data/Raw/pedigrees/",files[i]),header=TRUE,sep=if(type[i]=="txt"){""}else{","}, na.strings=c("NA",0))
  ped2 <- data.frame(animal=ped[,1],dam=ped[,"dam"],sire=ped[,"sire"])
  ped3 <- pedantics::fixPedigree(ped2)
  MCMCglmm::prunePed(ped3, keep=ped3[!is.na(ped3$dam),1], make.base=TRUE)
}, mc.cores=6)



# load(paste0(wd,"Data/Processed/phenotyped.Rdata"))
# plot(sim_all[2,],sim_juv[2,], log="xy")


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

  assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv)

# x<-data[[8]]
  sims<-do.call(rbind,parallel::mclapply(data, function(x){

    x$id <- as.factor(x$id)

    data_missing <- x[!is.na(ped[,"dam"]),]

    mod_1 <- asreml::asreml(
        fixed= y~1
        , random= ~vm(id,ped.ainv)
        , data= x, trace=FALSE)
    va_mod_2<- tryCatch({
      mod_2 <- asreml::asreml(
        fixed= y~1
        , random= ~vm(id,ped.ainv)
        , data= data_missing, trace=FALSE)
      summary(mod_2)$varcomp["vm(id, ped.ainv)",1]
      },
    error=function(e)0)
    
    c(summary(mod_1)$varcomp["vm(id, ped.ainv)",1],va_mod_2)
      

  },mc.cores=6))  
  rm("ped.ainv", envir = .GlobalEnv)
  c(var_all=var(sims[,1]),var_missing=var(sims[,2]))
}

# samp_var_sim_missing(peds[[2]],0.3, n_sim=1000)

all_missing <- sapply(peds,samp_var_sim_missing,0.3, n_sim=1000)
save(all_missing,file=paste0(wd,"Data/Processed/ped_missing.Rdata"))



n<-sapply(peds, function(x){
  c(nrow(x),sum(!is.na(x$dam)))
})

{
par(mfrow=c(1,2),mar=c(4,4,1,1))

plot(sqrt(t(sim_all)),log="xy")
abline(0,1)

plot(t(n))
abline(0,1)
}

{par(mfrow=c(1,3),mar=c(4,4,1,1))

plot(n[1,],sqrt(sim_all[1,]),log="xy")

plot(n[1,],sqrt(sim_all[2,]),log="xy")

plot(n[2,],sqrt(sim_all[2,]),log="xy")
}

save(sim_all,file=paste0(wd,"Data/Processed/phenotyped.Rdata"))