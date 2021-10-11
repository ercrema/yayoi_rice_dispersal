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
win <- as(sf_japan, "SpatialPolygons") |>
 unionSpatialPolygons(IDs = rep(1, nrow(sf_japan))) |>
  spTransform(CRSobj = CRS("+proj=utm +zone=54 ellps=WGS84"))
## Exclude small islands
win <- disaggregate(win) 
win  <- win[order(raster::area(win),decreasing=TRUE)[1:3]]
win.sf  <-  as(win,'sf')

# Tactical Simulation  ----
## Simulate Site Locations and Distances ----
set.seed(1323)
Nsites = 150 #number of simulated sites
sites <- spsample(win, n = Nsites, type = 'random') 

# Compute Distance from Northern Kyushu Site
origin_point = c(-520847.3,3751061)
distance.tmp <- rbind(origin_point, coordinates(sites)) |> 
 dist() |> as.matrix()
d <- distance.tmp[1, -1] / 1000
dmat <- distance.tmp[-1,-1] / 1000

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
		mu[k] <- beta0 + (s[k]+beta1)*d[k]
		alpha[k] ~ dnorm(mean=mu[k],sd=sigma)
		delta[k] ~ dgamma(shape=omega,rate=phi)
		alpha2[k] <- alpha[k] - delta[k] - 1 
	}
	for (i in 1:Ndates){
		Y[i] ~ dunif(min=alpha2[id.sites[i]],max=alpha[id.sites[i]]);
	}
	#fixed parameterts
	beta0 <- 3000
	omega <- 3
	phi <- 0.01
	beta1 <- -0.7
	sigma <- 100
	etasq <- 0.2
	rhosq <- 0.0001
	mu_s[1:Nsites] <- 0
	cov_s[1:Nsites, 1:Nsites] <- cov_GPL2(dmat[1:Nsites, 1:Nsites], rhosq, etasq, 0.00001)
	s[1:Nsites] ~ dmnorm(mu_s[1:Nsites], cov = cov_s[1:Nsites, 1:Nsites])
}) 


## Simulate Parameters ----
set.seed(123)
simModel <- nimbleModel(code = dispersalmodel,constants = list(Ndates = Ndates,d = d, dmat = dmat,Nsites = Nsites, id.sites = id.sites))
                       
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
save(SiteInfo,DateInfo,win.sf,dmat,d,file=here('results','tactical_sim_res.RData'))

# ggplot(win.sf) + geom_sf() + geom_point(data=SiteInfo,aes(x=Easting,y=Northing,col=s)) + scale_color_gradient2()
# ggplot(win.sf) + geom_sf() + geom_point(data=SiteInfo,aes(x=Easting,y=Northing,col=alpha)) + scale_color_gradient2(low = 'blue',mid ='white',high = 'red',midpoint = 2500)


