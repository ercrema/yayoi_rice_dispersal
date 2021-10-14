# Load Libraries and Data ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(parallel)
library(coda)

## Data Setup ----
# Load Simulated Data
load(here('results','tactical_sim_res.RData'))

# Data
dat  <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error)

# Constants
data(intcal20)
Ndates <- nrow(DateInfo)
Nsites  <- nrow(SiteInfo)
constants  <- list(Nsites=Nsites,Ndates=Ndates)
constants$id.sites  <- DateInfo$SiteID
constants$dmat  <- sim.constants$dist_mat
constants$d  <- sim.constants$dist_org
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$tau  <- 0.99

# Fixed Inits
buffer <- 100
theta.init <- DateInfo$med
delta.init <- SiteInfo$Diff + buffer
alpha.init <- SiteInfo$Earliest + buffer/2


# MCMC Function ----
runFun  <- function(seed,dat,constants,theta.init,delta.init,alpha.init,niter,nburnin,thin)
{
	library(nimbleCarbon) 
	library(truncnorm)
	# Variance Covariance Function	
	cov_GPL2 <- nimbleFunction(
				   run = function(dists = double(2), rhosq = double(0), etasq = double(0), sigmasq = double(0)) {
					   returnType(double(2))
					   n <- dim(dists)[1]
					   result <- matrix(nrow = n, ncol = n, init = FALSE)
					   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
					   diag(deltaij) <- 1
					   for(i in 1:n)
						   for(j in 1:n)
							   result[i, j] <- etasq*exp(-rhosq*dists[i,j]^2)+sigmasq*deltaij[i,j]
					   return(result)
				   })
	Ccov_GPL2 <- compileNimble(cov_GPL2)
	assign('cov_GPL2',cov_GPL2,envir=.GlobalEnv)

	# Core model
	gpqr <- nimbleCode({
		for (k in 1:Nsites){
			# Model
			mu[k] <- beta0 + (s[k]+beta1)*d[k]
			alpha[k] ~ dAsymLaplace(mu=mu[k],sigma=sigma,tau=tau)
			delta[k] ~ dgamma(shape=omega,rate=phi)
		}
		for (i in 1:Ndates){
			theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
			mu.date[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu.date[i],sd=sd[i]);
		}
		#priors
		beta0 ~ dnorm(3000,sd=200);
		beta1 ~ dnorm(0,sd=3)
		sigma ~ dexp(0.01)
		etasq ~ dexp(10);
		rhosq ~ dexp(1000);
		omega ~ T(dnorm(mean=3,sd=1),0,)
		phi ~ T(dnorm(mean=0.01,sd=0.01),0,)
		mu_s[1:Nsites] <- 0;
		cov_s[1:Nsites, 1:Nsites] <- cov_GPL2(dmat[1:Nsites, 1:Nsites], rhosq, etasq, 0.00001)
		s[1:Nsites] ~ dmnorm(mu_s[1:Nsites], cov = cov_s[1:Nsites, 1:Nsites])
	}) 

	# MCMC initialisation
	set.seed(seed)
	inits  <-  list(delta=delta.init,alpha=alpha.init)
	inits$theta  <- theta.init
	inits$beta0 <- rnorm(1,3000,200)
	inits$beta1 <- -rexp(1,1)
	inits$sigma  <- rexp(1,0.01)
	inits$rhosq  <- rexp(1,1000)
	inits$etasq  <- rexp(1,10)
	inits$omega  <- rtruncnorm(1,mean=3,sd=1,a=0)
	inits$phi  <- rtruncnorm(1,mean=0.01,sd=0.01,a=0)
	inits$s  <- rep(0,constants$Nsites)
	inits$cov_s <- Ccov_GPL2(constants$dmat, inits$rhosq, inits$etasq, 0.0001)
	inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$Nsites)
	inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix

	# Model Compilation
	model.gpqr <- nimbleModel(gpqr, constants = constants,data=dat,inits=inits)
	cModel.gpqr <- compileNimble(model.gpqr)

	# MCMC configuration
	conf.gpqr <- configureMCMC(model.gpqr)
	conf.gpqr$addMonitors('s')
	conf.gpqr$addMonitors('rhosq')
	conf.gpqr$addMonitors('etasq')
	conf.gpqr$addMonitors('alpha')
	conf.gpqr$addMonitors('delta')
	conf.gpqr$addMonitors('omega')
	conf.gpqr$addMonitors('phi')
	conf.gpqr$removeSamplers('s[1:150]')
	conf.gpqr$removeSamplers('beta1')
	conf.gpqr$addSampler(c('beta1','s[1:150]'), type='AF_slice') 
	MCMC.gpqr <- buildMCMC(conf.gpqr)
	cMCMC.gpqr <- compileNimble(MCMC.gpqr)

	# MCMC execution
	results <- runMCMC(cMCMC.gpqr, nchain=1,niter = niter, thin=thin, nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed) 
	return(results)
}

# Setup and Execute MCMC in Parallel  ----
ncores <- 3
cl <- makeCluster(ncores)
# Run the model in parallel:
seeds <- c(12, 34, 56)
niter = 1000000
nburnin = 500000
thin = 50
chain_output <- parLapply(cl = cl, X = seeds, fun = runFun, dat = dat,
constants = constants, theta = theta.init, alpha.init = alpha.init, delta.init = delta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
# Convert into a mcmc.list object for diagnostic (see below)
res <- mcmc.list(chain_output)
rhat <- gelman.diag(res,multivariate = F)
ess <- effectiveSize(res) 

## Store Output ----
save(res,rhat,ess,file=here('results','res_tactical.RData'))












