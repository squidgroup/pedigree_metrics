

samp_var_sim <- function(ped, h2, m2=0, n_sim=200){
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
			n_pop=n_sim
		),
	list=TRUE)

	assign("ped.ainv", nadiv::makeAinv(ped)$Ainv, envir = .GlobalEnv)

	sims<-c(parallel::mclapply(data, function(x){
		x$id <- as.factor(x$id)
			mod_1 <- gremlin::gremlin(
				y~1
				, random= ~id
			  , data= x
				, ginverse = list(id = ped.ainv)
				,v=0
			  )
	  summary(mod_1)$varcompSummary["G.id","Estimate"]

		# mod_1 <- asreml::asreml(
		# 		fixed= y~1
		# 		, random= ~vm(id,ped.ainv)
		# 	  , data= x, trace=FALSE)
		# summary(mod_1)$varcomp["vm(id, ped.ainv)",1]
	},mc.cores=6), recursive=TRUE)	
	rm("ped.ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims))
}


MDE_sim <- function(ped,n_sims=1000){

	assign("ped.ainv", nadiv::makeAinv(ped)$Ainv, envir = .GlobalEnv)

	sims<-do.call(rbind,parallel::mclapply(1:n_sims, function(x){

		dat <- data.frame(y=rnorm(nrow(ped)), id=as.factor(ped[,1]))
		
		mod_1 <- gremlin::gremlin(
			y~1
			, random= ~id
		  , data= dat
			, ginverse = list(id = ped.ainv)
			,v=0
		  )

		summary(mod_1)$varcompSummary["G.id","Estimate"]

	},mc.cores=6))	
	rm("ped.ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims), MDE=quantile(sims,0.95))
}

# grm<-soays_grm
# h2=0.3
#  n_sim=10
samp_var_sim_GRM <- function(grm, h2=0.3, n_sim=200){

	grm<-as(grm,"dgCMatrix")

	data<-squidSim::get_population_data(
		squidSim::simulate_population(
			data_structure=data.frame(id=rownames(grm)),
			parameters=list(
				id=list(vcov=h2),
				residual=list(vcov=1-h2)
			),
			cov_str=list(id=grm),
			n_pop=n_sim
		),
	list=TRUE)
# x<- data.frame(id=rownames(grm), y=rnorm(nrow(grm)))
	ainv <- solve(grm)
	# assign("ped.ainv", nadiv::makeAinv(ped)$Ainv, envir = .GlobalEnv)

	sims<-c(parallel::mclapply(data, function(x){
		x$id <- as.factor(x$id)
			mod_1 <- gremlin::gremlin(
				y~1
				, random= ~id
			  , data= x
				, ginverse = list(id = ainv)
				,v=0
			)
	  summary(mod_1)$varcompSummary["G.id","Estimate"]

		# mod_1 <- asreml::asreml(
		# 		fixed= y~1
		# 		, random= ~vm(id,ped.ainv)
		# 	  , data= x, trace=FALSE)
		# summary(mod_1)$varcomp["vm(id, ped.ainv)",1]
	},mc.cores=6), recursive=TRUE)	
	# rm("ainv", envir = .GlobalEnv)
	c(mean=mean(sims),var=var(sims))
}



## function to go from eigenvalue of A matrix to approximation of sampling variance for h2
lambda_to_samp_var <- function(lambda,h2){
    N <- length(lambda)
  sapply(h2, function(x){
    a <- sum( (lambda-1)^2 / (1 + x*(lambda-1))^2 )
    b <- sum( (lambda-1) / (1 + x*(lambda-1)) )
    2/(a-b^2/N)
  })  
}
