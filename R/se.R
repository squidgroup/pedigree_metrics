
rm(list=ls())


library(asreml)
library(parallel)
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")


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


se_ped_sims <- function(n_sims=100,years=5, n_females=100, fecundity=4, p_sire=0.5, juv_surv = 0.5){

	ped <- simulate_pedigree(
					years = years,
					n_females = n_females,
					fecundity = fecundity,
					p_sire = p_sire,
					p_polyandry=1, 					# all females have same probability or EPP
					juv_surv = juv_surv,
					adult_surv = 0,					# discrete generations
					immigration = 0, 				# closed population
					constant_pop = TRUE     # constant population size
					)$pedigree

	## expected sampling variance
	samp_var <- sqrt(samp_var_h2(ped,c(0.1,0.3,0.5,0.7,0.9)))

	assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 

	out <- as.data.frame(do.call(rbind,parallel::mclapply(c(0.1,0.3,0.5,0.7,0.9),function(h2){
		t(replicate(n_sims,{
			
			# generate data
			data <- data.frame(
				p = MCMCglmm::rbv(ped[,1:3],h2) 
						+ rnorm(nrow(ped),0,sqrt(1-h2)),
				animal=as.factor(ped$animal)
				)
			
			## standardise response to unit variance 
			data$pS <- scale(data$p)

			# run model
			mod_1 <- asreml(
				fixed= pS~1
				, random= ~vm(animal,ped.ainv)
			  , data= data, trace=FALSE)
			c(sim_h2 =h2, summary(mod_1)$varcomp[1,1:2],recursive=TRUE)
		}))
	}, mc.cores=5)))

	list(samp_var=samp_var,sims=out)

}

plot_sim <- function(sims){
	sims_sum<-aggregate(component~sim_h2,sims$sims,sd)

	plot(component~sim_h2,sims_sum, ylim=range(c(sims$samp_var,sims_sum$component)), pch=19, ylab="sampling sd")
	points(sims$samp_var~c(0.1,0.3,0.5,0.7,0.9), pch=19,col="red")
	legend("topleft", c("observed","expected"), pch=19, col=c(1,2))
}
	

sims1 <- se_ped_sims(n_sims=100,
				years = 5,
				n_females = 50,
				fecundity = 8,
				p_sire = 0.5,
				juv_surv = 0.25)
aggregate(std.error~sim_h2,sims1$sims,mean)

plot_sim(sims1)

sims2 <- se_ped_sims(n_sims=100,
				years = 5,
				n_females = 100,
				fecundity = 4,
				p_sire = 0,
				juv_surv = 0.5)

plot_sim(sims2)

sims3 <- se_ped_sims(n_sims=100,
				years = 5,
				n_females = 100,
				fecundity = 4,
				p_sire = 1,
				juv_surv = 0.5)

plot_sim(sims3)


sims4 <- se_ped_sims(n_sims=100,
				years = 2,
				n_females = 200,
				fecundity = 4,
				p_sire = 1,
				juv_surv = 0.5)

plot_sim(sims4)

(rbind(sims1$samp_var,sims2$samp_var,sims3$samp_var,sims4$samp_var))^2

(rbind(aggregate(component~sim_h2,sims1$sims,sd)$component,
aggregate(component~sim_h2,sims2$sims,sd)$component,
aggregate(component~sim_h2,sims3$sims,sd)$component,
aggregate(component~sim_h2,sims4$sims,sd)$component))^2

(rbind(aggregate(std.error~sim_h2,sims1$sims,mean)$std.error,
aggregate(std.error~sim_h2,sims2$sims,mean)$std.error,
aggregate(std.error~sim_h2,sims3$sims,mean)$std.error,
aggregate(std.error~sim_h2,sims4$sims,mean)$std.error))^2



ped <- simulate_pedigree(
				years = 2,
				n_females = 100,
				fecundity = 4,
				p_sire = 0.5,
				p_polyandry=1,
				juv_surv = 0.5,
				adult_surv = 0,					# discrete generations
				immigration = 0, 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
nrow(ped)

# A<-nadiv::makeA(ped[,1:3])

# samp_var <- 2/(var(log(eigen(2*A, symmetric=TRUE,only.values=TRUE)$values))*nrow(ped))

# samp_var <- 2/(var(log(eigen(2*A, symmetric=TRUE,only.values=TRUE)$values))*nrow(ped))

# ped.ainv<- asreml::ainverse(ped)

# out <- as.data.frame(do.call(rbind,parallel::mclapply(c(0.1,0.3,0.5,0.7,0.9),function(h2){
# 	t(replicate(100,{
# 			data <- data.frame(
# 				p = MCMCglmm::rbv(ped[,1:3],h2) 
# 						+ rnorm(nrow(ped),0,sqrt(1-h2)),
# 				animal=as.factor(ped$animal)
# 				)
# 			data$pS <- scale(data$p)
# 			mod_1 <- asreml(
# 				fixed= pS~1
# 				, random= ~vm(animal,ped.ainv)
# 			  , data= data, trace=FALSE)
# 			c(sim_h2 =h2, summary(mod_1)$varcomp[1,1:2],recursive=TRUE)
# 			}))
# }, mc.cores=5)))

# aggregate(std.error~sim_h2,out,mean)
# aggregate(component~sim_h2,out,sd)
# aggregate(component~sim_h2,out,mean)

# plot(component~sim_h2,out)


# ped2 <- simulate_pedigree(
# 				years = 5,
# 				n_females = 50,
# 				fecundity = 8,
# 				p_sire = 0.5,
# 				p_polyandry=1,
# 				juv_surv = 0.25,
# 				adult_surv = 0,					# discrete generations
# 				immigration = 0, 				# closed population
# 				constant_pop = TRUE     # constant population size
# 				)$pedigree


# A2<-nadiv::makeA(ped2[,1:3])

# samp_var2 <- 2/(var(log(eigen(2*A2, symmetric=TRUE,only.values=TRUE)$values))*nrow(ped2))





# ped2.ainv<- asreml::ainverse(ped2)

# out2 <- as.data.frame(do.call(rbind,parallel::mclapply(c(0.1,0.3,0.5,0.7,0.9),function(h2){
# 	t(replicate(100,{
# 			data <- data.frame(
# 				p = MCMCglmm::rbv(ped2[,1:3],h2) 
# 						+ rnorm(nrow(ped2),0,sqrt(1-h2)),
# 				animal=as.factor(ped2$animal)
# 				)
	
# 			mod_1 <- asreml(
# 				fixed= p~1
# 				, random= ~vm(animal,ped2.ainv)
# 			  , data= data, trace=FALSE)
# 			c(sim_h2 =h2, summary(mod_1)$varcomp[1,1:2],recursive=TRUE)
# 			}))
# }, mc.cores=5)))

# samp_var2



# aggregate(std.error~sim_h2,out2,mean)
# aggregate(component~sim_h2,out2,sd)
# aggregate(component~sim_h2,out2,mean)


# plot(std.error~sim_h2,out2)


