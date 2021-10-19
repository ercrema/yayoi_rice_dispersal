# Load Libraries ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(quantreg)

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
subset.DateInfo  <- subset.DateInfo[order(subset.DateInfo$SiteID, decreasing = F),[
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
runFun <- function(seed, d, theta.init, constants, nburnin, thin, niter)
{
	library(nimbleCarbon)
	model <- nimbleCode({
		for (i in 1:N){
			# Model
			mu[i] <- alpha + beta*d_org[i]
			theta[i] ~ dAsymLaplace(mu=mu[i],sigma=sigma,tau=tau)

			c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
		}
		#priors
		alpha ~ dnorm(3000,sd=200);
		beta ~ dnorm(0,sd=3)
		sigma ~ dexp(0.05)
	}) 
	set.seed(seed)
	inits  <- list(alpha=rnorm(1,3000,200),beta=rnorm(1,0,3),sigma=rexp(1,0.05),theta=theta.init)
	model.asymlap <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel.asymlap <- compileNimble(model.asymlap)
	conf.asymlap <- configureMCMC(model.asymlap)
	conf.asymlap$addMonitors('theta')
	MCMC.asymlap <- buildMCMC(conf.asymlap)
	cMCMC.asymlap <- compileNimble(MCMC.asymlap)
	results <- runMCMC(cMCMC.asymlap, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(results)
}       

# Setup and Execution of MCMC in Parallel ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds  <-  c(12,34,56,78)
niter  <- 100000
nburnin  <- 50000
thin  <- 5

chain_output = parLapply(cl = cl, X = seeds, fun = runFun, d = d,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
# Convert into a mcmc.list object for diagnostic (see below)
quantreg_sample=coda::mcmc.list(chain_output)

## Store Output ----
save(quantreg_sample,file=here('results','quantreg_res.RData'))




