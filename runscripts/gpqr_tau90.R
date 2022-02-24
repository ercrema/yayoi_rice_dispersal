# Load Libraries and Data ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(parallel)

## Data Setup ----
# Read 14C dates
load(here('data','c14rice.RData'))

# Dummy extension of the calibration curve
constants$calBP <- c(1000000,constants$calBP,-1000000)
constants$C14BP <- c(1000000,constants$C14BP,-1000000)
constants$C14err <- c(1000,constants$C14err,1000)


## Define Run Function ----
runFun  <- function(seed,dat,theta.init,constants,niter,nburnin,thin)
{
	library(nimbleCarbon)
	library(truncnorm)
	library(cascsim)

	# Variance Covariance Function	
	cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0),sigmasq = double(0)) 
				   {
					   returnType(double(2))
					   n <- dim(dists)[1]
					   result <- matrix(nrow = n, ncol = n, init = FALSE)
					   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
					   diag(deltaij) <- 1
					   for(i in 1:n)
						   for(j in 1:n)
							   result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
					   return(result)
				   })
	Ccov_ExpQ <- compileNimble(cov_ExpQ)
	assign('cov_ExpQ',cov_ExpQ,envir=.GlobalEnv)

	# Handle constraints on dipsersal rate
	dat$lim  <- rep(1,constants$N.sites)
	# Core model
	model <- nimbleCode({
		for (i in 1:N.sites){
			# Model
			rate[i] <- -1/(s[i]-beta1)
			lim[i] ~ dconstraint(rate[i]>0)
			mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
			theta[i] ~ dAsymLaplace(mu=mu[i],sigma=sigma,tau=tau)
			mu.date[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu.date[i],sd=sd[i]);
		}
		#priors
		beta0 ~ dnorm(3000,sd=200);
		beta1 ~ dexp(1)
		sigma ~ dexp(0.01)
		etasq ~ dexp(20);
		rho ~ T(dgamma(10,(10-1)/150),1,1350); #mode 150
		mu_s[1:N.sites] <- 0;
		cov_s[1:N.sites, 1:N.sites] <- cov_ExpQ(dist_mat[1:N.sites, 1:N.sites], rho, etasq, 0.000001)
		s[1:N.sites] ~ dmnorm(mu_s[1:N.sites], cov = cov_s[1:N.sites, 1:N.sites])
	}) 

	# MCMC initialisation
	set.seed(seed)
	inits  <-  list()
	inits$theta  <- theta.init
	inits$beta0 <- rnorm(1,3000,200)
	inits$beta1 <- rexp(1,1)
	inits$sigma  <- rexp(1,0.01)
	inits$rho  <- rtgamma(1,shape=10,scale=(10-1)/200,min=1,max=1350)
	inits$etasq  <- rexp(1,20)
	inits$s  <- rep(0,constants$N.sites)
	inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
	inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$N.sites)
	inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

	# Model Compilation
	model.gpqr <- nimbleModel(model,constants = constants,data=dat,inits=inits)
	cModel.gpqr <- compileNimble(model.gpqr)

	# MCMC configuration
	conf.gpqr <- configureMCMC(model.gpqr)
	conf.gpqr$addMonitors('s')
	conf.gpqr$addMonitors('rho')
	conf.gpqr$addMonitors('etasq')
	conf.gpqr$removeSamplers('s[1:206]')
	conf.gpqr$removeSamplers('beta1')
	conf.gpqr$addSampler(c('beta1','s[1:206]'), type='AF_slice') 
	MCMC.gpqr <- buildMCMC(conf.gpqr)
	cMCMC.gpqr <- compileNimble(MCMC.gpqr)

	# MCMC execution
	results <- runMCMC(cMCMC.gpqr, nchain=1,niter = niter, thin=thin, nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed) 
	return(results)
}


## Define Constants and Data ----
DateInfo  <- subset(DateInfo,EarliestAtSite==TRUE)
DateInfo  <- DateInfo[order(DateInfo$SiteID,decreasing = FALSE),]

# Data
dat  <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error)

# Constants
constants$tau <- 0.90

# Fixed Inits
theta.init <- medCal(calibrate(dat$cra,dat$cra_error))

## Setup and Execution of MCMC in Parallel ----
ncores <- 4
cl <- makeCluster(ncores)
# Run the model in parallel:
seeds <- c(12,45,67,89)
niter = 8000000
nburnin = 4000000
thin = 400
chain_output <- parLapply(cl = cl, X = seeds, fun = runFun, dat = dat,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
# Convert into a mcmc.list object for diagnostic (see below)
gpqr_tau90 <- coda::mcmc.list(chain_output)
rhat.t90 <- coda::gelman.diag(gpqr_tau90,multivariate = F)
# range(rhat.t90$psrf[,1])
# ii  <- which(rhat.t90$psrf[,1]>1.01)  
ess.t90  <- coda::effectiveSize(gpqr_tau90)
# range(ess.t90)

## Store Output ----
save(rhat.t90,ess.t90,gpqr_tau90,file=here('results','gpqr_tau90.RData'))
