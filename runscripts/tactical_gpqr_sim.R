# Load Libraries ----
library(here)
library(nimbleCarbon)
library(rnaturalearth)
library(sp)
library(maptools)
library(sf)
library(ggplot2)
library(viridis)
library(rgeos)
library(rcarbon)
source(here('src','gpqrSim.R'))

# Generate Spatial Window for Analyses ----
sf_japan <- ne_states(country = "japan") |> subset(!name_vi %in%  c("Okinawa","Hokkaido"))
sampling.win <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan)))
## Exclude small islands
sampling.win <- disaggregate(sampling.win) 
sampling.win  <- sampling.win[order(raster::area(sampling.win),decreasing=TRUE)[1:3]]

# Target Parameters ----
true.param  <- list()
true.param$n  <- 150 #number of sites & dates
true.param$origin.point <- c(129.95,33.45) #dispersal origin point
true.param$beta0 <- 3100 #mean date at origin point
true.param$beta1 <- 0.6 #reciprocal of dispersal rate
true.param$sigma <- 100 #variability of arrival date
true.param$etasq <- 0.06 #variability of dispersal rate
true.param$rho <- 150 #range of spatial autocorrelation
true.param$seed <- 1233 #random seed

# Simulate Data ----
sim.sites  <- gpqrSim(win=sampling.win,n=true.param$n,beta0=true.param$beta0,beta1=true.param$beta1,sigma=true.param$sigma,origin.point=true.param$origin.point,etasq=true.param$etasq,rho=true.param$rho,seed=true.param$seed)

# Save simulation output ----
save(true.param,sim.sites,file=here('data','tactical_sim_gpqr.RData'))
