# Load Libraries and Data ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(parallel)
library(sf)
library(maptools)
library(coda)

## Data Setup ----
# Load Simulated Data
load(here('data','tactical_sim_gpqr.RData'))

# Constants
data(intcal20)
constants  <- list()
constants$N <- true.param$n
constants$dist_mat  <- spDists(as_Spatial(sim.sites),longlat=TRUE)
constants$dist_org  <- spDistsN1(as_Spatial(sim.sites),true.param$origin.point,longlat=TRUE)
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$tau  <- 0.90

# Data
dat  <- list()
dat$cra  <- sim.sites$cra
dat$cra_error  <- sim.sites$cra.error

# Theta Init
theta.init  <-  sim.sites$med.date

# MCMC Function ----
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
	dat$lim  <- rep(1,constants$N)
	# Core model
	model <- nimbleCode({
		for (i in 1:N){
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
		mu_s[1:N] <- 0;
		cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
		s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
	}) 

	# MCMC initialisation
	set.seed(seed)
	inits  <-  list()
	inits$theta  <- theta.init
	inits$beta0 <- rnorm(1,3000,200)
	inits$beta1 <- rexp(1,0.5)
	inits$sigma  <- rexp(1,0.01)
	inits$rho  <- rtgamma(1,shape=10,scale=(10-1)/200,min=1,max=1350)
	inits$etasq  <- rexp(1,20)
	inits$s  <- rep(0,constants$N)
	inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
	inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$N)
	inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

	# Model Compilation
	model.gpqr <- nimbleModel(model,constants = constants,data=dat,inits=inits)
	cModel.gpqr <- compileNimble(model.gpqr)

	# MCMC configuration
	conf.gpqr <- configureMCMC(model.gpqr)
	conf.gpqr$addMonitors('s')
	conf.gpqr$addMonitors('rho')
	conf.gpqr$addMonitors('etasq')
	conf.gpqr$removeSamplers('s[1:150]')
	conf.gpqr$removeSamplers('beta1')
	conf.gpqr$addSampler(c('beta1','s[1:150]'), type='AF_slice') 
	MCMC.gpqr <- buildMCMC(conf.gpqr)
	cMCMC.gpqr <- compileNimble(MCMC.gpqr)

	# MCMC execution
	results <- runMCMC(cMCMC.gpqr, nchain=1,niter = niter, thin=thin, nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed) 
	return(results)
}


## Setup and Execution of MCMC in Parallel ----
ncores <- 4
cl <- makeCluster(ncores)
# Run the model in parallel:
seeds <- c(12,45,67,89)
niter = 1000000
nburnin = 500000
thin = 50
chain_output <- parLapply(cl = cl, X = seeds, fun = runFun, dat = dat,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
# Convert into a mcmc.list object for diagnostic (see below)
gpqr_tactsim <- coda::mcmc.list(chain_output)
rhat <- coda::gelman.diag(gpqr_tactsim,multivariate = F)
range(rhat$psrf[,1])
ii  <- which(rhat$psrf[,1]>1.01)  
ess  <- coda::effectiveSize(gpqr_tactsim)
range(ess)

## Store Output ----
save(gpqr_tactsim,file=here('results','gpqr_tactsim.RData'))
