# Load Libraries ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(quantreg)
library(parallel)

# Load and prepare data ----
load(here('data','c14rice.RData'))

# Quantile Regression based on median earliest date at each site ----
# Add distance from origin to SiteInfo
SiteInfo$dist_org <- constants$dist_org
# Compute Quantile Regression
fit.median  <- rq(Earliest ~ dist_org, tau = 0.99, data=SiteInfo,alpha=0.95)
# Derive estimated rate of dispersal
-1/summary(fit.median)$coefficients[2,]
## between 0.62 to 11.4 km/year, with estimate at 4.08 km/year

# Bayesian Quantile Regression -----

## Setup Data, Constants, and Inits ----

# Data
## Consider only the earliest date at each site
subset.DateInfo  <- subset(DateInfo, EarliestAtSite == TRUE)
subset.DateInfo  <- subset.DateInfo[order(subset.DateInfo$SiteID, decreasing = F),]
## Generate list of observed data
dat  <- list(cra = subset.DateInfo$cra, cra_error = subset.DateInfo$cra_error)

# Constants
# Remove unused slots
constants$N.sites <- NULL
constants$id.sites  <- NULL
constants$dist_mat  <- NULL
# Update number of dates
constants$N.dates  <- nrow(subset.DateInfo)
# Define Quantile
constants$tau <- 0.99

# Init (theta)
theta.init  <- subset.DateInfo$median.dates

## Main Scipt ----
runFun <- function(seed, dat, theta.init, constants, nburnin, thin, niter)
{
	options(warn=-1)
	library(nimbleCarbon)
	model <- nimbleCode({
		for (i in 1:N.dates){
			# Model
			mu[i] <- alpha - beta*dist_org[i]
			theta[i] ~ dAsymLaplace(mu=mu[i],sigma=sigma,tau=tau)

			c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
		}
		#priors
		alpha ~ dnorm(3000,sd=200);
		beta ~ dexp(1)
		sigma ~ dexp(0.01)
	}) 
	set.seed(seed)
	inits  <- list(alpha=rnorm(1,3000,200),beta=rexp(1,1),sigma=rexp(1,0.01),theta=theta.init)
	model.asymlap <- nimbleModel(model,constants = constants,data=dat,inits=inits)

	#The while loop ensures that the none of the log probabilities of theta are -Inf. This issues 
	#a warning message that crashes the parallel processing script. Notice that if done sequentially there is no problem.
	while(any(model.asymlap$logProb_theta==-Inf))
	{
	inits  <- list(alpha=rnorm(1,3000,200),beta=rexp(1,1),sigma=rexp(1,0.01),theta=theta.init)
	model.asymlap <- nimbleModel(model,constants = constants,data=dat,inits=inits)	
	}
	cModel.asymlap <- compileNimble(model.asymlap)
	conf.asymlap <- configureMCMC(model.asymlap)
	conf.asymlap$addMonitors('theta')
	MCMC.asymlap <- buildMCMC(conf.asymlap)
	cMCMC.asymlap <- compileNimble(MCMC.asymlap)
	results <- runMCMC(cMCMC.asymlap, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(results)
}       

# Setup and Execution of MCMC in Parallel ----
ncores  <-  3
cl <- makeCluster(ncores)
seeds  <-  c(12,45,78)
niter  <- 1000000
nburnin  <- 500000
thin  <- 50

chain_output = parLapply(cl = cl, X = seeds, fun = runFun, d = dat, constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
        
# Convert into a mcmc.list object for diagnostic (see below)
quantreg_sample <- coda::mcmc.list(chain_output)
rhat <- coda::gelman.diag(quantreg_sample)
## Store Output ----
save(quantreg_sample,file=here('results','quantreg_res.RData'))




