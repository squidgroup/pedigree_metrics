
ped_stat <- function(ped, phenotyped=NULL){	

	colnames(ped) <- c("animal","dam","sire")

	# dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	# dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	# ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	# ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))


	## parent pairs
	ped$pair<-paste(ped$dam,ped$sire)
	ped$pair<-ifelse(ped$pair== "NA NA", NA, ped$pair)


	## FOUNDERS
	founders <- ped[is.na(ped$pair),1]

	pd<-do.call(rbind,parallel::mclapply(founders, function(parent) {
		# loop though ped for each founder
		for(j in 1:100){
			parent <- subset(ped,(dam %in% parent | sire %in% parent ))[,1]
			if(length(parent)==0) break
		}
		j
	},mc.cores=8))
# table(pd)



  ## grandparents
	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	## grandparent pairs
	ped$MGP <- paste(ped$MGM,ped$MGF)
	ped$PGP <- paste(ped$PGM,ped$PGF)
	# ped$MGP<-ifelse(ped$MGP== "NA NA", NA, ped$MGP)
	# ped$PGP<-ifelse(ped$PGP== "NA NA", NA, ped$PGP)
	ped$MGP<-ifelse(grepl("NA",ped$MGP), NA, ped$MGP)
	ped$PGP<-ifelse(grepl("NA",ped$PGP), NA, ped$PGP)

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]
	## subset to phenotyped individuals
	ped2 <- subset(ped,animal %in% phenotyped)	


####################
# ---- parents and grandparents
####################

	gp <- apply(ped[ped[,1] %in% phenotyped,c("dam","sire","MGM","MGF","PGM","PGF")],2,function(x) sum(x%in% phenotyped) )


####################
# ---- sibs
####################

	mat_sib <- sum(table(ped2$dam) * (table(ped2$dam)-1)/2)
	pat_sib <- sum(table(ped2$sire) * (table(ped2$sire)-1)/2)
	ped2$pair <- paste(ped2$dam,ped2$sire)
	n_pair <-table(ped2$pair[!grepl("NA",ped2$pair)])
	FS <- sum(n_pair * (n_pair-1)/2)
	
## relatedness

	sp<-Matrix::summary(nadiv::makeA(ped2[,1:3]))
# sum(sp[,3])/ (nrow(ped) * (nrow(ped) - 1) / 2 + nrow(ped))

	sp2<-sp[sp[,1]!=sp[,2],]

	
####################
# ---- put together
####################

# [1] "records "                    "maternities"                
#  [3] "paternities"                 "full sibs"                  
#  [5] "maternal sibs"               "maternal half sibs"         
#  [7] "paternal sibs"               "paternal half sibs"         
#  [9] "maternal grandmothers"       "maternal grandfathers"      
# [11] "paternal grandmothers"       "paternal grandfathers"      
# [13] "maximum pedigree depth"      "founders"                   
# [15] "mean maternal sibsip size"   "mean paternal sibsip size"  
# [17] "non-zero F"                  "F > 0.125"                  
# [19] "mean pairwise relatedness"   "pairwise relatedness>=0.125"
# [21] "pairwise relatedness>=0.25"  "pairwise relatedness>=0.5"

	c(
		records = length(phenotyped),
		maternities = sum(!is.na(ped[,2])),
		paternities = sum(!is.na(ped[,3])),
		FS = FS, 
		MS = mat_sib, 
		MHS = mat_sib-FS, 
		PS = pat_sib, 
		PHS = pat_sib-FS,
		MGM = as.vector(gp["MGM"]),
		MGF = as.vector(gp["MGF"]), 
		PGM = as.vector(gp["PGM"]), 
		PGF = as.vector(gp["PGF"]),
		"mean pedigree depth" = mean(pd),
		"maximum pedigree depth" = max(pd),
		## take the founders and loop over to see how far done can trace decendants?
		founders= sum(is.na(ped$pair)),
		"mean maternal sibsip size" = mean(aggregate(animal~dam,ped,length)$animal),
    "mean paternal sibsip size"= mean(aggregate(animal~sire,ped,length)$animal),

		"non-zero F" = sum(MCMCglmm::inverseA(ped[,1:3])$inbreeding>0),
		"F > 0.125" = sum(MCMCglmm::inverseA(ped[,1:3])$inbreeding>0.125),

		"mean pairwise relatedness" = sum(sp2[,3])/ (nrow(ped) * (nrow(ped) - 1) / 2 ),

		"pairwise relatedness>=0.125" = sum(sp2[,3]>=0.125),
    "pairwise relatedness>=0.25" = sum(sp2[,3]>=0.25),
    "pairwise relatedness>=0.5" = sum(sp2[,3]>=0.5)
	)
	
}

# cbind(ped_stat(ped), pedantics::pedStatSummary(pedantics::pedigreeStats(ped[,1:3],includeA=FALSE,lowMem=FALSE,graphicalReport=FALSE))
# )


# pedantics::pedStatSummary(pedantics::pedigreeStats(ped[,1:3],includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))
 





samp_var_h2 <- function(ped, h2){
	A<-nadiv::makeA(ped[,1:3])
	lambda <- eigen(A, symmetric=TRUE,only.values=TRUE)$values
	N <- nrow(ped)
	sapply(h2, function(x){
		a <- sum( (lambda-1)^2 / (1 + x*(lambda-1))^2 )
		b <- sum( (lambda-1) / (1 + x*(lambda-1)) )
		2/(a-b^2/N)
	})	
}

samp_var_h2_ext <- function(ped, h2){
	A<-nadiv::makeA(ped[,1:3])
	lambda <- eigen(A, symmetric=TRUE,only.values=TRUE)$values
	N <- nrow(ped)
	sapply(h2, function(x){
		a <- sum( (lambda-1)^2 / (1 + x*(lambda-1))^2 )
		b <- sum( (lambda-1) / (1 + x*(lambda-1)) )
		c(a=a,b=b,N=N,sv=2/(a-b^2/N))
	})	
}
# extract a,b, N and samp var
# minimum detectable effect



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

	sims<-c(parallel::mclapply(data, function(x){
		x$id <- as.factor(x$id)
		mod_1 <- asreml::asreml(
				fixed= y~1
				, random= ~vm(id,ped.ainv)
			  , data= x, trace=FALSE)
		summary(mod_1)$varcomp["vm(id, ped.ainv)",1]
	},mc.cores=6), recursive=TRUE)	
	rm("ped.ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims))
}

MDE_sim <- function(ped,n_sims=1000){

	assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv)

	sims<-do.call(rbind,parallel::mclapply(1:n_sims, function(x){
		dat <- data.frame(y=rnorm(nrow(ped)), id=as.factor(ped[,1]))
		mod_1 <- asreml::asreml(
				fixed= y~1
				, random= ~vm(id,ped.ainv)
			  , data= dat, trace=FALSE)
		summary(mod_1)$varcomp["vm(id, ped.ainv)",1]
	},mc.cores=6))	
	rm("ped.ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims), MDE=quantile(sims,0.95))
}


# ped<-peds[[2]]
samp_var_sim_mat <- function(ped, h2, m2=0){
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

#x<-data[[1]]
	sims<-c(parallel::mclapply(data, function(x){
		x <- x[!is.na(ped[,"dam"]),] #subset(x,!is.na(dam))
		x$id <- as.factor(x$id)
		x$dam <- as.factor(x$dam)
		mod_1 <- asreml::asreml(
				fixed= y~1
				, random= ~vm(id,ped.ainv) + dam
			  , data= x, trace=FALSE)
		summary(mod_1)$varcomp["vm(id, ped.ainv)",1]
	},mc.cores=6), recursive=TRUE)	
	rm("ped.ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims))
}
