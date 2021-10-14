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
win <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan)))
## Exclude small islands
win <- disaggregate(win) 
win  <- win[order(raster::area(win),decreasing=TRUE)[1:3]]
win.sf  <-  as(win,'sf')

# Tactical Simulation  ----
## Simulate Site Locations and Distances ----
set.seed(1323)
Nsites = 150 #number of simulated sites
sites <- spsample(win, n = Nsites, type = 'random') 

dist_mat  <- spDists(sites,longlat=TRUE) #inter-site distance matrix
origin_point  <- c(129.95798382,33.44851467) #Nabatake Site 
dist_org  <-  spDistsN1(sites,origin_point,longlat=TRUE) #distance from nabatake site

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
	cov_s[1:Nsites, 1:Nsites] <- cov_GPL2(dist_mat[1:Nsites, 1:Nsites], rhosq, etasq, 0.00001)
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
sim.constants$etasq <- 0.1
sim.constants$rhosq <- 0.002
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
coords  <-  as.data.frame(coordinates(sites))
dates <- data.frame(cra = round(uncalibrate(round(simModel$Y))$ccCRA), cra_error = 20)
dates$med  <-  medCal(calibrate(dates$cra,dates$cra_error))

DateInfo = data.frame(SiteID=id.sites,cra=dates$cra,cra_error=dates$cra_error,med=dates$med)
SiteInfo = data.frame(SiteID=1:nrow(coords),Easting=coords$x,Northing=coords$y,Earliest=NA,Latest=NA,Diff=NA,alpha=simModel$alpha,s=simModel$s)

for (i in 1:Nsites)
{
	tmp.dates <- subset(DateInfo,SiteID==SiteInfo$SiteID[i])
	SiteInfo$Earliest[i] = max(tmp.dates$med)
	SiteInfo$Latest[i] = min(tmp.dates$med)
	SiteInfo$Diff[i] = SiteInfo$Earliest[i] - SiteInfo$Latest[i]
}


# Save simulation output ----
save(sim.constants,SiteInfo,DateInfo,win,file=here('results','tactical_sim_res.RData'))

# plot(sim.constants$dist_org,BPtoBCAD(SiteInfo$Earliest),pch=20)
# ggplot(win.sf) + geom_sf() + geom_point(data=SiteInfo,aes(x=Easting,y=Northing,col=s)) + scale_color_gradient2()
# ggplot(win.sf) + geom_sf() + geom_point(data=SiteInfo,aes(x=Easting,y=Northing,col=alpha)) + scale_color_gradient2(low = 'blue',mid ='white',high = 'red',midpoint = 2500)


