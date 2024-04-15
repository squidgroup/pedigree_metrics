rm(list=ls())

wd <- "/Users/joelpick/github/pedigree_metrics/"

source(paste0(wd,"R/00_functions.R"))


files<-list.files(paste0(wd,"Data/Raw/pedigrees"))
type <- sapply(strsplit(files,"\\."),function(x) x[2])
peds<-lapply(1:length(files), function(i){
  ped <- read.table(paste0(wd,"Data/Raw/pedigrees/",files[i]),header=TRUE,sep=if(type[i]=="txt"){""}else{","}, na.strings=c("NA",0))
  ped2 <- data.frame(animal=ped[,1],dam=ped[,"dam"],sire=ped[,"sire"])
  pedantics::fixPedigree(ped2)
})




plot(sim_all[2,],sim_juv[2,], log="xy")







# ped<-peds[[2]]
# h2=0.3;m2=0
samp_var_sim <- function(ped, h2, m2=0){
  ## give all individuals without a mother a dummy mother, so get a phenotype
  ped_ds <- ped[,1:2]
  ped_ds[,"dam"]<-as.character(ped_ds[,"dam"])
  ped_ds[is.na(ped_ds[,"dam"]),"dam"] <- paste0("dd",1:sum(is.na(ped_ds[,"dam"])))
  data<-squidSim::get_population_data(
    squidSim::simulate_population(
      data_structure=ped_ds,
      parameters=list(
        id=list(vcov=h2),
        dam=list(vcov=m2),
        residual=list(vcov=1-h2-m2)
      ),
      pedigree=list(id=ped),
      n_pop=200
    ),
  list=TRUE)

  assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv)

# x<-data[[8]]
  sims<-do.call(rbind,parallel::mclapply(data, function(x){

    x$id <- as.factor(x$id)

    x2 <- x[!is.na(ped[,"dam"]),]

    mod_1 <- asreml::asreml(
        fixed= y~1
        , random= ~vm(id,ped.ainv)
        , data= x, trace=FALSE)
    va_mod_2<- tryCatch({
      mod_2 <- asreml::asreml(
        fixed= y~1
        , random= ~vm(id,ped.ainv)
        , data= x2, trace=FALSE)
      summary(mod_2)$varcomp["vm(id, ped.ainv)",1]
      },
    error=function(e)0)
    
# ,ai.sing=TRUE, fail="soft"
    c(summary(mod_1)$varcomp["vm(id, ped.ainv)",1],va_mod_2)
      

  },mc.cores=6))  
  rm("ped.ainv", envir = .GlobalEnv)
  c(var_all=var(sims[,1]),var_p=var(sims[,2]))
}

sim_all <- sapply(peds,samp_var_sim,0.3)

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