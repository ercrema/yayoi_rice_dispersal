# Load Library and Data ----
library(here)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
load(here("data","c14rice.RData"))

# General Setup ----
# Data
d <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error,constraint_trapezium=rep(1,constants$N.areas))
# Inits
buffer <- 100
theta.init <- DateInfo$median.dates
delta.init <- SiteInfo$Diff + buffer
alpha.init <- SiteInfo$Earliest + buffer/2
# Constants
constants$constraint_trapezium  <- rep(1,constants$N.areas)
constants$constraint_dispersal  <- 1


# MCMC RunScript (Uniform Model0) ----
unif.model0  <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dunif(max=a[id.area[k]],min=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>b[j])
		}
	})

	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	init.a = init.b =  numeric(length=constants$N.areas)
	for (i in 1:constants$N.areas)
	{
		tmp.alpha <- alpha.init[which(constants$id.area==i)]
		tmp.delta <- delta.init[which(constants$id.area==i)]
		tmp.beta <- tmp.alpha - tmp.delta
		init.a[i] <- runif(1,max(tmp.alpha),5000)
		init.b[i]  <- runif(1,500,min(tmp.beta))
	}

	inits$a  <- init.a
	inits$b  <- init.b

	# Compile and Run model	
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
}

# MCMC RunScript (Uniform Model1) ----

unif.model1 <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		# Model Start/End Dates at Individual Sites	
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dunif(max=a[id.area[k]],min=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>b[j])
		}

		# Define Dispersal Constraint
		constraint_dispersal ~ dconstraint(a[1]>a[2] & a[1]>a[3] & a[3]>a[4] & a[4]>a[5] & a[5]>a[6] & a[6]>a[7] & a[7]>a[8])  
	})
	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	check = TRUE
	while(check)
	{
		init.a = init.b = rep(NA,length=constants$N.areas)
		for (i in 1:constants$N.areas)
		{
			tmp.alpha <- alpha.init[which(constants$id.area==i)]
			tmp.delta <- delta.init[which(constants$id.area==i)]
			tmp.beta <- tmp.alpha - tmp.delta
			if (i == 1) {init.a[i] <- runif(1,max(alpha.init),5000)}
			if (i > 1) {
				if(max(tmp.alpha)<init.a[i-1])
				{
					init.a[i] <- runif(1,max(tmp.alpha),init.a[i-1])
				}
				if(max(tmp.alpha)>=init.a[i-1])
				{
					break()
				}
			}
			
			init.b[i]  <- runif(1,500,min(tmp.beta))
		}
		if (all(!is.na(init.a))){check=FALSE}
	}

	inits$a  <- init.a
	inits$b  <- init.b
	
	#Add constraint
	d$constraint_dispersal  <- 1

	#Compile and Run
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
	return(results)
}

# MCMC RunScript (Uniform Model2) ----

unif.model2 <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		# Model Start/End Dates at Individual Sites	
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dunif(max=a[id.area[k]],min=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>b[j])
		}

		# Define Dispersal Constraint
		constraint_dispersal ~ dconstraint(a[1]>a[2] & a[1]>a[3] & a[3]>a[4] & a[4]>a[5] & a[4]>a[6] & a[4]>a[7] & a[4]>a[8])  
	})
	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	check = TRUE
	while(check)
	{
		init.a = init.b = rep(NA,length=constants$N.areas)
		for (i in 1:constants$N.areas)
		{
			tmp.alpha <- alpha.init[which(constants$id.area==i)]
			tmp.delta <- delta.init[which(constants$id.area==i)]
			tmp.beta <- tmp.alpha - tmp.delta
			if (i == 1) {init.a[i] <- runif(1,max(alpha.init),5000)}
			if (i > 1) {
				if(max(tmp.alpha)<init.a[i-1])
				{
					init.a[i] <- runif(1,max(tmp.alpha),init.a[i-1])
				}
				if(max(tmp.alpha)>=init.a[i-1])
				{
					break()
				}
			}
			
			init.b[i]  <- runif(1,500,min(tmp.beta))
		}
		if (all(!is.na(init.a))){check=FALSE}
	}

	inits$a  <- init.a
	inits$b  <- init.b
	
	#Add constraint
	d$constraint_dispersal  <- 1

	#Compile and Run
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
	return(results)
}

# MCMC RunScript (Trapezoidal Model0) ----

trap.model0 <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		# Model Start/End Dates at Individual Sites	
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dTrapezoidal(a=a[id.area[k]],m1=m1[id.area[k]],m2=m2[id.area[k]],b=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			m1[j] ~dunif(500,5000);
			m2[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>m1[j] & m1[j]>m2[j] & m2[j]>b[j])
		}
	})

	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	init.a = init.b  = init.m1 = init.m2  = numeric(length=constants$N.areas)
	for (i in 1:constants$N.areas)
	{
		tmp.alpha <- alpha.init[which(constants$id.area==i)]
		tmp.delta <- delta.init[which(constants$id.area==i)]
		tmp.beta <- tmp.alpha - tmp.delta
		init.a[i] <- runif(1,max(tmp.alpha),5000)
		init.b[i]  <- runif(1,500,min(tmp.beta))
		tmp.m  <- sort(runif(2,init.b[i],init.a[i]),decreasing=TRUE)
		init.m1[i] <- tmp.m[1]
		init.m2[i]  <- tmp.m[2]
	}

	inits$a  <- init.a
	inits$b  <- init.b
	inits$m1  <- init.m1
	inits$m2  <- init.m2

        # Compile and Run model	
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
	return(results)
}


# MCMC RunScript (Trapezoidal Model1) ----

trap.model1 <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		# Model Start/End Dates at Individual Sites	
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dTrapezoidal(a=a[id.area[k]],m1=m1[id.area[k]],m2=m2[id.area[k]],b=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			m1[j] ~dunif(500,5000);
			m2[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>m1[j] & m1[j]>m2[j] & m2[j]>b[j])
		}

		# Define Dispersal Constraint
		constraint_dispersal ~ dconstraint(a[1]>a[2] & a[1]>a[3] & a[3]>a[4] & a[4]>a[5] & a[5]>a[6] & a[6]>a[7] & a[7]>a[8])  
	})
	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	check = TRUE
	while(check)
	{
		init.a = init.b = init.m1 = init.m2 = rep(NA,length=constants$N.areas)
		for (i in 1:constants$N.areas)
		{
			tmp.alpha <- alpha.init[which(constants$id.area==i)]
			tmp.delta <- delta.init[which(constants$id.area==i)]
			tmp.beta <- tmp.alpha - tmp.delta
			if (i == 1) {init.a[i] <- runif(1,max(alpha.init),5000)}
			if (i > 1) {
				if(max(tmp.alpha)<init.a[i-1])
				{
					init.a[i] <- runif(1,max(tmp.alpha),init.a[i-1])
				}
				if(max(tmp.alpha)>=init.a[i-1])
				{
					break()
				}
			}
			
			init.b[i]  <- runif(1,500,min(tmp.beta))
			tmp.m  <- sort(runif(2,init.b[i],init.a[i]),decreasing=TRUE)
			init.m1[i] <- tmp.m[1]
			init.m2[i]  <- tmp.m[2]
		}
		if (all(!is.na(init.a))){check=FALSE}
	}

	inits$a  <- init.a
	inits$b  <- init.b
	inits$m1  <- init.m1
	inits$m2  <- init.m2
	
	#Add constraint
	d$constraint_dispersal  <- 1

	#Compile and Run
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
	return(results)
}

# MCMC RunScript (Trapezoidal Model2) ----

trap.model2 <- function(seed, d, theta.init, alpha.init, delta.init, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		# Model Start/End Dates at Individual Sites	
		for (k in 1:N.sites)
		{
			delta[k] ~ dexp(gamma)
			alpha[k] ~ dTrapezoidal(a=a[id.area[k]],m1=m1[id.area[k]],m2=m2[id.area[k]],b=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Duration
		gamma ~ dunif(0.002,1)

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(500,5000);
			m1[j] ~dunif(500,5000);
			m2[j] ~ dunif(500,5000);
			b[j] ~ dunif(500,5000);
			constraint_trapezium[j] ~ dconstraint(a[j]>m1[j] & m1[j]>m2[j] & m2[j]>b[j])
		}

		# Define Dispersal Constraint
		constraint_dispersal ~ dconstraint(a[1]>a[2] & a[1]>a[3] & a[3]>a[4] & a[4]>a[5] & a[4]>a[6] & a[4]>a[7] & a[4]>a[8])  
	})
	# Define Inits
	set.seed(seed)
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init)
	inits$gamma <- runif(1,0.002,1)
	check = TRUE
	while(check)
	{
		init.a = init.b = init.m1 = init.m2 = rep(NA,length=constants$N.areas)
		for (i in 1:constants$N.areas)
		{
			tmp.alpha <- alpha.init[which(constants$id.area==i)]
			tmp.delta <- delta.init[which(constants$id.area==i)]
			tmp.beta <- tmp.alpha - tmp.delta
			if (i == 1) {init.a[i] <- runif(1,max(alpha.init),5000)}
			if (i > 1) {
				if(max(tmp.alpha)<init.a[i-1])
				{
					init.a[i] <- runif(1,max(tmp.alpha),init.a[i-1])
				}
				if(max(tmp.alpha)>=init.a[i-1])
				{
					break()
				}
			}
			
			init.b[i]  <- runif(1,500,min(tmp.beta))
			tmp.m  <- sort(runif(2,init.b[i],init.a[i]),decreasing=TRUE)
			init.m1[i] <- tmp.m[1]
			init.m2[i]  <- tmp.m[2]
		}
		if (all(!is.na(init.a))){check=FALSE}
	}

	inits$a  <- init.a
	inits$b  <- init.b
	inits$m1  <- init.m1
	inits$m2  <- init.m2
	
	#Add constraint
	d$constraint_dispersal  <- 1

	#Compile and Run
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('theta')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
	return(results)
}


# Run MCMCs ----

# MCMC Setup
ncores  <-  3
cl <- makeCluster(ncores)
seeds <- c(123,456,789)
niter  <- 1000000
nburnin  <- 500000
thin  <- 50

out.unif.model0  <-  parLapply(cl = cl, X = seeds, fun = unif.model0, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)
out.unif.model1  <-  parLapply(cl = cl, X = seeds, fun = unif.model1, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)
out.unif.model2  <-  parLapply(cl = cl, X = seeds, fun = unif.model2, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)

out.trap.model0  <-  parLapply(cl = cl, X = seeds, fun = trap.model0, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)
out.trap.model1  <-  parLapply(cl = cl, X = seeds, fun = trap.model1, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)
out.trap.model2  <-  parLapply(cl = cl, X = seeds, fun = trap.model2, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  niter = niter, nburnin = nburnin,thin = thin)


out.unif.model0 <- mcmc.list(out.unif.model0)
out.unif.model1 <- mcmc.list(out.unif.model1)
out.unif.model2 <- mcmc.list(out.unif.model2)
out.trap.model0 <- mcmc.list(out.trap.model0)
out.trap.model1 <- mcmc.list(out.trap.model1)
out.trap.model2 <- mcmc.list(out.trap.model2)


# Diagnostics ----

rhat.unif.model0 <- gelman.diag(out.unif.model0,multivariate = FALSE)
rhat.unif.model1 <- gelman.diag(out.unif.model1,multivariate = FALSE)
rhat.unif.model2 <- gelman.diag(out.unif.model2,multivariate = FALSE)
rhat.trap.model0 <- gelman.diag(out.trap.model0,multivariate = FALSE)
rhat.trap.model1 <- gelman.diag(out.trap.model1,multivariate = FALSE)
rhat.trap.model2 <- gelman.diag(out.trap.model2,multivariate = FALSE)

ess.unif.model0 <- effectiveSize(out.unif.model0)
ess.unif.model1 <- effectiveSize(out.unif.model1)
ess.unif.model2 <- effectiveSize(out.unif.model2)
ess.trap.model0 <- effectiveSize(out.trap.model0)
ess.trap.model1 <- effectiveSize(out.trap.model1)
ess.trap.model2 <- effectiveSize(out.trap.model2)

a.unif.model0 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.unif.model0[[1]][,grep("theta",colnames(out.unif.model0[[1]]))],verbose = F)
a.unif.model1 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.unif.model1[[1]][,grep("theta",colnames(out.unif.model1[[1]]))],verbose = F)
a.unif.model2 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.unif.model2[[1]][,grep("theta",colnames(out.unif.model2[[1]]))],verbose = F)
a.trap.model0 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.trap.model0[[1]][,grep("theta",colnames(out.trap.model0[[1]]))],verbose = F)
a.trap.model1 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.trap.model1[[1]][,grep("theta",colnames(out.trap.model1[[1]]))],verbose = F)
a.trap.model2 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.trap.model2[[1]][,grep("theta",colnames(out.trap.model2[[1]]))],verbose = F)


# Save output ----
save(out.unif.model0,rhat.unif.model0,ess.unif.model0,a.unif.model0,file=here("results","unif_model0.RData"))
save(out.unif.model1,rhat.unif.model1,ess.unif.model1,a.unif.model1,file=here("results","unif_model1.RData"))
save(out.unif.model2,rhat.unif.model2,ess.unif.model2,a.unif.model2,file=here("results","unif_model2.RData"))
save(out.trap.model0,rhat.trap.model0,ess.trap.model0,a.trap.model0,file=here("results","trap_model0.RData"))
save(out.trap.model1,rhat.trap.model1,ess.trap.model1,a.trap.model1,file=here("results","trap_model1.RData"))
save(out.trap.model2,rhat.trap.model2,ess.trap.model2,a.trap.model2,file=here("results","trap_model2.RData"))



