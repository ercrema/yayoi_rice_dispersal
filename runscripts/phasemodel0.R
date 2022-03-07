# Load Library and Data ----
library(here)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
load(here("data","c14rice.RData"))

# General Setup ----
# Data
d <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error,constraint_uniform=rep(1,constants$N.areas))
# Constraint for ignoring inference outside calibration range
d$cra.constraint = rep(1,constants$N.dates)
# Inits
buffer <- 100
theta.init <- DateInfo$median.dates
delta.init <- SiteInfo$Diff + buffer
alpha.init <- SiteInfo$Earliest + buffer/2

# Dummy extension of the calibration curve
constants$calBP <- c(1000000,constants$calBP,-1000000)
constants$C14BP <- c(1000000,constants$C14BP,-1000000)
constants$C14err <- c(1000,constants$C14err,1000)

# Regional Init
init.a  <- aggregate(Earliest~Area,FUN=max,data=SiteInfo)[,2] + 100
init.b  <- aggregate(Latest~Area,FUN=min,data=SiteInfo)[,2] - 100




# MCMC RunScript (Uniform Model0) ----
unif.model0  <- function(seed, d, theta.init, alpha.init, delta.init, init.a, init.b, constants, nburnin, thin, niter)
{
	#Load Library
	library(nimbleCarbon)
	#Define Core Model
	model <- nimbleCode({
		for (k in 1:N.sites)
		{
			delta[k] ~ dgamma(gamma1,(gamma1-1)/gamma2)
			alpha[k] ~ dunif(max=a[id.area[k]],min=b[id.area[k]]);
		}

		for (i in 1:N.dates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			cra.constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		# Set Prior for Each Region
		for (j in 1:N.areas){
			a[j] ~ dunif(50,5000);
			b[j] ~ dunif(50,5000);
			constraint_uniform[j] ~ dconstraint(a[j]>b[j])
		}
		# Hyperprior for duration
		gamma1 ~ dunif(1,20)
		gamma2 ~ T(dnorm(mean=200,sd=100),1,500)
	})

	# Define Inits
	inits <- list(theta=theta.init, alpha=alpha.init, delta=delta.init, a=init.a, b=init.b)
	inits$gamma1  <- 10
	inits$gamma2  <- 200

	# Compile and Run model	
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,control=list(adaptInterval=20000,adaptFactorExponent=0.1))
	conf$addMonitors(c('theta','delta','alpha'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
}

# MCMC RunScript (Uniform Model1) ----
# Run MCMCs ----

# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 8000000
nburnin  <- 4000000
thin  <- 200

out.unif.model0  <-  parLapply(cl = cl, X = seeds, fun = unif.model0, d = d,constants = constants, theta.init = theta.init, alpha.init = alpha.init, delta.init = delta.init,  init.a = init.a, init.b = init.b, niter = niter, nburnin = nburnin,thin = thin)

out.unif.model0 <- mcmc.list(out.unif.model0)

# Diagnostics ----

rhat.unif.model0 <- gelman.diag(out.unif.model0,multivariate = FALSE)
ess.unif.model0 <- effectiveSize(out.unif.model0)
a.unif.model0 <- agreementIndex(d$cra,d$cra_error,calCurve='intcal20',theta=out.unif.model0[[1]][,grep("theta",colnames(out.unif.model0[[1]]))],verbose = F)

# Save output ----
save(out.unif.model0,rhat.unif.model0,ess.unif.model0,a.unif.model0,file=here("results","phase_model0.RData"))
