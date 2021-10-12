# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(dplyr)
library(here)

# Read charred rice data  ----
dat <- read.csv(here("data", "R14CDB_v1_4_5.csv")) |> subset(Rice == "TRUE" & UseForAnalyses == 'TRUE')
#**NOTE** update raw CSV data for the final version so that subsetting is not required

# Subset data ----
# Only Dates between 3000 and 1000 C14age, excluding Hokkaido and Okinawa
dat <- subset(dat, !Region %in% c("Hokkaido", "Okinawa") & C14Age <= 3000 & C14Age>1000 )

# Assign Site ID ----
dat$SiteID  <- as.numeric(factor(dat$SiteName))

# Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates
dat$median.dates = medCal(calibrate(dat$C14Age,dat$C14Error,calCurve='intcal20'))

# Collect site level information
earliest_dates <- aggregate(median.dates~SiteID,data=dat,FUN=max) #Earliest medCal Date for Each Site
latest_dates <- aggregate(median.dates~SiteID,data=dat,FUN=min) #Latest medCal Date for Each Site
n_dates <- aggregate(median.dates~SiteID,data=dat,FUN=length) #Number of medCal Date for Each Site
SiteInfo <- data.frame(SiteID = earliest_dates$SiteID,
		       Earliest = earliest_dates$median.dates,
		       Latest = latest_dates$median.dates,
		       Diff = earliest_dates$median.dates - latest_dates$median.dates,
		       N_dates = n_dates$median.dates) |> unique()
SiteInfo <- left_join(SiteInfo,unique(select(dat,Latitude,Longitude,SiteID,SiteName,SiteName_En)))

# Collect date level information
DateInfo <- unique(select(dat,ID,LabCode,SiteID,cra=C14Age,cra_error=C14Error,median.dates=median.dates)) |> arrange(ID) 
DateInfo$EarliestAtSite  <- FALSE

for (i in unique(SiteInfo$SiteID))
{
	tmp.index  <-  which(DateInfo$SiteID==i)
	ii  <- tmp.index[which.max(DateInfo$median.dates[tmp.index])]
	DateInfo$EarliestAtSite[ii]  <- TRUE
}

# Compute Great-Arc Distances in km
sites <- SiteInfo
coordinates(sites) <- c('Longitude','Latitude')
proj4string(sites)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
dist_mat  <- spDists(sites,longlat=TRUE) #inter-site distance matrix
origin_point  <- c(129.95798382,33.44851467) #Nabatake Site 
dist_org  <-  spDistsN1(sites,origin_point,longlat=TRUE) #distance from nabatake site

# Create list with constants and data

## Data
d  <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error)

## Constants
data(intcal20)
constants <- list()
constants$N.sites <- nrow(SiteInfo)
constants$N.dates  <- nrow(DateInfo)
constants$id.sites <- DateInfo$SiteID
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$calBP <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma

# Save everything on a R image file ----
save(constants,d,SiteInfo,DateInfo,file=here('data','c14rice.RData'))


