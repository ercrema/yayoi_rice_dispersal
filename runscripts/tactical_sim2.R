library(rcarbon)
library(nimbleCarbon)
library(here)

set.seed(123)
Nsites <- 10
Ndates <- 30
id.sites <- c(1:Nsites,sample(1:Nsites,size=Ndates-Nsites,replace=TRUE,prob=dexp(1:Nsites,rate=1)/sum(dexp(1:Nsites,rate=1))))

sim.model <- nimbleCode({
	for (k in 1:Nsites)
	{
		delta[k] ~ dgamma(5,(5-1)/200);
		alpha[k] ~ dunif(max=a,min=b);
		beta[k] <- alpha[k] - delta[k] + 1;
	}

	for (i in 1:Ndates){
		theta[i] ~ dunif(beta[id.sites[i]],alpha[id.sites[i]]);
	}
})

sim.constants <- list()
sim.constants$Nsites <- Nsites
sim.constants$Ndates  <- Ndates
sim.constants$id.sites  <- id.sites
sim.constants$a <- 3500
sim.constants$b <- 3000

set.seed(123)
simModel <- nimbleModel(code = sim.model,constants = sim.constants)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine and store output
cra = uncalibrate(round(simModel$theta))$rCRA
cra_error = rep(20,length(cra))
d.sim <- list(cra=cra,cra_error=cra_error,id.sites=id.sites)
save(d.sim,sim.constants,file=here('results','tactical_sim_res2.RData'))
