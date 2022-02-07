#Load Libraries and Data ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
load(here('data','tactical_sim_phase.RData'))

# Model assuming independence of samples ----
model1 <- nimbleCode({
	for (i in 1:Ndates)
	{
		theta[i] ~ dunif(b,a);
		# Calibration
		mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		error[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=mu[i],sd=error[i]);
	}
	a ~ dunif(0,10000);
	b ~ dunif(0,10000);
	unif.const ~ dconstraint(b<a);
})

data("intcal20") 
constants1 <- list(Ndates=sim.constants$Ndates,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
d1 <- list(cra=d.sim$cra,cra_error=d.sim$cra_error,unif.const=1)
theta.init = medCal(calibrate(d1$cra,d1$cra_error,verbose = FALSE))
inits1 <- list(a=5000,b=500,theta=theta.init)

#Run MCMC
mcmc.samples1<- nimbleMCMC(code = model1,constants = constants1,data = d1,niter = 2000000, nchains = 3, thin=100, nburnin = 1000000,monitors=c('a','b','theta'), inits=inits1, samplesAsCodaMCMC=TRUE)

#Diagnostics
rhat1  <- gelman.diag(mcmc.samples1,multivariate = FALSE)
ess1  <- effectiveSize(mcmc.samples1)




# Model integrating sample interdependence ----
model2 <- nimbleCode({
	for (k in 1:Nsites)
	{
		delta[k] ~ dgamma(gamma1,(gamma1-1)/gamma2)
		alpha[k] ~ dunif(max=a,min=b);
	}

	for (i in 1:Ndates){
		theta[i] ~ dunif(alpha[id.sites[i]]-delta[id.sites[i]]+1,alpha[id.sites[i]]);
		mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		error[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=mu[i],sd=error[i]);
	}

	# Set Prior for Duration
	a ~ dunif(0,10000);
	b ~ dunif(0,10000);	unif.const ~ dconstraint(b<a);
	gamma1 ~ dunif(1,20);
	gamma2 ~ T(dnorm(mean=200,sd=100),1,500);
})

data("intcal20") 
constants2 <- list(Ndates=sim.constants$Ndates,
		   Nsites=sim.constants$Nsites,
		   calBP=intcal20$CalBP,
		   C14BP=intcal20$C14Age,
		   C14err=intcal20$C14Age.sigma,
		   id.sites=sim.constants$id.sites)
d2 <- list(cra=d.sim$cra,cra_error=d.sim$cra_error,unif.const=1)
theta.init  <-  medCal(calibrate(d1$cra,d1$cra_error,verbose = FALSE))
dd  <-  data.frame(theta=theta.init,id=constants2$id.sites)
buffer  <- 100
earliest  <- aggregate(theta~id,max,data=dd)
latest  <- aggregate(theta~id,min,data=dd)
diff.age  <- earliest$theta - latest$theta
delta.init  <- diff.age + buffer
alpha.init  <- earliest$theta + buffer/2


inits2 <- list(a=5000,b=500,gamma1=5,gamma2=200,theta=theta.init,alpha=alpha.init,delta=delta.init)

#Run MCMC
mcmc.samples2<- nimbleMCMC(code = model2,constants = constants2,data = d2,niter = 2000000, nchains = 3, thin=100, nburnin = 1000000,monitors=c('a','b','theta','gamma1','gamma2'), inits=inits2, samplesAsCodaMCMC=TRUE)

#Diagnostics
rhat2  <- gelman.diag(mcmc.samples2,multivariate = FALSE)
ess2  <- effectiveSize(mcmc.samples2)

# Save output ----
save(mcmc.samples1,rhat1,ess1,mcmc.samples2,rhat2,ess2,file=here('results','phasemodel_tactsim.RData'))
