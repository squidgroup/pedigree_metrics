
rm(list=ls())

library(pedtricks)
wd <- "/Users/joelpick/github/maternal_effects/"


ped_BT<-read.csv(paste0(wd,"Data/Raw/BT_Pick.csv"))[,1:3]
ped_RD<-read.table(paste0(wd,"Data/Raw/RD_Gauzere.txt"),header=TRUE)[,c(1,3,2)]

nrow(ped_BT)

ggpedigree(ped_BT)
ggpedigree(ped_RD)