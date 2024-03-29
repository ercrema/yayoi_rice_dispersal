# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(dplyr)
library(here)

# Read charred rice data  ----
dat <- read.csv(here("data", "R14CDB.csv")) |> subset(UseForAnalyses == 'TRUE'& C14Age > 1000)

# Pool Samples TKA-23237 amd TKA-23238
poolLabCodes  <- c('TKA-23237','TKA-23238')
i  <- which(dat$LabCode%in%poolLabCodes)
pooledDates  <- poolDates(dat$C14Age[i],dat$C14Error[i])[2:3]
dat  <- rbind.data.frame(dat,dat[i[1],])
dat$ID[nrow(dat)] <- max(dat$ID) + 1
dat$LabCode[nrow(dat)] <- 'Combined_TKA23237_TKA23238'
dat$C14Age[nrow(dat)] <- as.numeric(pooledDates[1])
dat$C14Error[nrow(dat)] <- as.numeric(pooledDates[2])
dat$Context[nrow(dat)] <- 'Charred rice grain embedded in sherd fabric (Combined Dates)'
dat  <- dat[-i,]


# Assign Site ID ----
dat$SiteID  <- as.numeric(factor(dat$SiteName_jp))

# Assign Regions & Areas ----
regions <- read.csv(here("data","prefecture_region_match.csv"))
dat <- left_join(dat,regions)

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
SiteInfo <- left_join(SiteInfo,unique(select(dat,Latitude,Longitude,SiteID,SiteName_jp,SiteName_en,Prefecture,Region,Area)))
SiteInfo$area.id <- as.numeric(as.factor(SiteInfo$Area))

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
origin_point  <- c(130.0087,33.411) #Ukikunden Shell Midden 
dist_org  <-  spDistsN1(sites,origin_point,longlat=TRUE) #distance from nabatake site

# Create list with constants and data

## Data
dat  <- list(cra=DateInfo$cra,cra_error=DateInfo$cra_error)

## Constants
data(intcal20)
constants <- list()
constants$N.sites <- nrow(SiteInfo)
constants$N.dates  <- nrow(DateInfo)
constants$N.areas  <- length(unique(SiteInfo$Area))
constants$id.sites <- DateInfo$SiteID
constants$id.area  <- SiteInfo$area.id
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$calBP <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma

# Save everything on a R image file ----
save(sites,constants,dat,SiteInfo,DateInfo,file=here('data','c14rice.RData'))
