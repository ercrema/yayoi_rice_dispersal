# Load Libraries ----
library(here)
library(nimbleCarbon)
library(rnaturalearth)
library(sp)
library(maptools)
library(rgeos)
library(rcarbon)

# Generate Spatial Window for Analyses ----
sf_japan <- ne_states(country = "japan") |> subset(!name_vi %in%  c("Okinawa","Hokkaido"))
sampling.win <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan)))
## Exclude small islands
sampling.win <- disaggregate(sampling.win) 
sampling.win  <- sampling.win[order(raster::area(sampling.win),decreasing=TRUE)[1:3]]

# Tactical Simulation  ----

## Simulate Site Locations and Distances ----
set.seed(12345)
Nsites = 150 #number of simulated sites
sim.sites <- spsample(sampling.win, n = Nsites, type = 'random') 

dist_mat  <- spDists(sim.sites,longlat=TRUE) #inter-site distance matrix
origin_point  <- c(129.95798382,33.44851467) #Nabatake Site 
dist_org  <-  spDistsN1(sim.sites,origin_point,longlat=TRUE) #distance from nabatake site

# Associate dates to sites
Ndates = 400 #total number of dates
# ensure that each site has at least one date
id.sites <- c(1:Nsites,sample(1:Nsites,size=Ndates-Nsites,replace=TRUE))

## Setup Dispersal Model in Nimble ----

# Source and compile Variance-Covariance function
source(here('src','cov_GPL2.R'))
Ccov_GPL2 <- compileNimble(cov_GPL2)

# Core Model
dispersalmodel <- nimbleCode({
	for (k in 1:Nsites){
		# Model
		mu[k] <- beta0 + (s[k]+beta1)*dist_org[k]
		alpha[k] ~ dnorm(mean=mu[k],sd=sigma)
		delta[k] ~ dgamma(shape=omega,rate=phi)
		alpha2[k] <- alpha[k] - delta[k] - 1 
	}
	for (i in 1:Ndates){
		Y[i] ~ dunif(min=alpha2[id.sites[i]],max=alpha[id.sites[i]]);
	}
	mu_s[1:Nsites] <- 0
	cov_s[1:Nsites, 1:Nsites] <- cov_GPL2(dist_mat[1:Nsites, 1:Nsites], rhosq, etasq, 0.000001)
	s[1:Nsites] ~ dmnorm(mu_s[1:Nsites], cov = cov_s[1:Nsites, 1:Nsites])
}) 

## Define Parameters ----
sim.constants  <- list()
sim.constants$Ndates <- Ndates
sim.constants$Nsites <- Nsites
sim.constants$dist_org <- dist_org
sim.constants$dist_mat <- dist_mat
sim.constants$id.sites <- id.sites
# Model parameters
sim.constants$beta0 <- 3000
sim.constants$beta1 <- -0.7
sim.constants$sigma <- 100
sim.constants$etasq <- 0.05
sim.constants$rhosq <- 0.00005
sim.constants$omega <- 2
sim.constants$phi <- 0.005




## Simulate Parameters ----
set.seed(123)
simModel <- nimbleModel(code = dispersalmodel,constants = sim.constants)
                       
simModel$simulate('s')
simModel$simulate('mu')
simModel$simulate('alpha')
simModel$simulate('delta')
simModel$simulate('alpha2')
simModel$simulate('Y')

### Combine Results ----
coords  <-  as.data.frame(coordinates(sim.sites))
dates <- data.frame(cra = round(uncalibrate(round(simModel$Y))$ccCRA), cra_error = 20)
dates$med  <-  medCal(calibrate(dates$cra,dates$cra_error))

sim.DateInfo = data.frame(SiteID=id.sites,cra=dates$cra,cra_error=dates$cra_error,med=dates$med)
sim.SiteInfo = data.frame(SiteID=1:nrow(coords),Longitude=coords$x,Latitude=coords$y,Earliest=NA,Latest=NA,Diff=NA,alpha=simModel$alpha,s=simModel$s)

for (i in 1:Nsites)
{
	tmp.dates <- subset(sim.DateInfo,SiteID==sim.SiteInfo$SiteID[i])
	sim.SiteInfo$Earliest[i] = max(tmp.dates$med)
	sim.SiteInfo$Latest[i] = min(tmp.dates$med)
	sim.SiteInfo$Diff[i] = sim.SiteInfo$Earliest[i] - sim.SiteInfo$Latest[i]
}


sim.sites <- as(sim.sites,'SpatialPointsDataFrame')
sim.sites@data = data.frame(s=simModel$s)
sim.sites@data$rate=-1/(sim.constants$beta1 + simModel$s)
sim.sites@data$arrival = BPtoBCAD(sim.constants$beta0 + (sim.constants$beta1 + simModel$s) * sim.constants$dist_org)
sim.sites <- as(sim.sites,'sf')

# range(sim.sites$s)
# range(sim.sites$rate)
# range(sim.sites$arrival)

# Save simulation output ----
save(sim.constants,sim.SiteInfo,sim.DateInfo,sim.sites,file=here('results','tactical_sim_res.RData'))
