# Load Libraries and spatial data ----
library(here)
library(rcarbon) 
library(dplyr)
library(coda)
# Data File S1 Rice Radiocarbon Database ----
datafile.S1  <- read.csv(here('data','R14CDB.csv'))
ii  <- which(datafile.S1$C14Age>0)
x  <- calibrate(datafile.S1$C14Age[ii],datafile.S1$C14Error[ii],calCurve='intcal20')
datafile.S1$medianDate  <- datafile.S1$sigmaRange  <-  NA
datafile.S1$medianDate[ii] <- paste(abs(BPtoBCAD(medCal(x))),ifelse(BPtoBCAD(medCal(x))<0,'BC','AD'))
datafile.S1$sigmaRange[ii] <- hpdi(x,calendar='BCAD',credMass=0.9545,asList=FALSE)
datafile.S1$used = datafile.S1$UseForAnalyses
datafile.S1$used[which(datafile.S1$Prefecture%in%c('Hokkaido','Okinawa'))]  <- FALSE
datafile.S1 <- select(datafile.S1,ID,LabCode,C14Age,C14Error,medianDate,sigmaRange,Method,SiteName_jp,SiteName_en,Latitude,Longitude,Context,Prefecture,Region,PublicationYear,Reference,Notes,used)

write.table(datafile.S1,file=here('manuscript','supplementary_tables','datafile_S1.csv'),col.names=c('ID','LabCode','C14Age','C14AgeError','MedianDate','Two-Sigma Range','Method','SiteName_jp','SiteName_en','Laitude','Longitude','Context','Prefecture','Region','PublicationYear','Reference','Notes','USED'),sep=',',row.names=FALSE)

# Table S1 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.9) ----
load(here('results','gpqr_tau90.RData'))
gpqr.tau90.comb  <- do.call(rbind,gpqr_tau90)
params = c('beta0','beta1','rho','etasq')
gpqr.tau90.comb[,'beta1'] = 1/gpqr.tau90.comb[,'beta1']
meds = apply(gpqr.tau90.comb[,params],2,median) |> round(3)
lo90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]}) |> round(3)
hi90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]}) |> round(3)
rhats = gelman.diag(gpqr_tau90)$psrf[params,1]
ess = effectiveSize(gpqr_tau90)[params]
params[2]='1/beta1'
table.S1 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S1,file=here('manuscript','supplementary_tables','table_S1.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)


# Table S2 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.99) ----
load(here('results','gpqr_tau99.RData'))
gpqr.tau99.comb  <- do.call(rbind,gpqr_tau99)
params = c('beta0','beta1','rho','etasq')
gpqr.tau90.comb[,'beta1'] = 1/gpqr.tau90.comb[,'beta1']
meds = apply(gpqr.tau99.comb[,params],2,median) |> round(3)
lo90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]}) |> round(3)
hi90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]}) |> round(3)
rhats = gelman.diag(gpqr_tau99)$psrf[params,1]
ess = effectiveSize(gpqr_tau99)[params]
params[2]='1/beta1'
table.S2 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S2,file=here('manuscript','supplementary_tables','table_S2.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)

# Table S3 Summary Data Per Region ----
load(here('data','c14rice.RData'))
pref.regions  <- read.csv(here('data','prefecture_region_match.csv'))
region.names  <- as.character(as.roman(1:8))
region.prefs  <- character(8)
for (i in 1:8)
{
	region.prefs[i] = paste0(subset(pref.regions,Area==paste0('Area',i))$Prefecture,collapse="|")
}
n.sites <- as.numeric(aggregate(SiteID~Area,data=SiteInfo,length)$SiteID)
n.dates <- as.numeric(aggregate(N_dates~Area,data=SiteInfo,sum)$N_dates)
table.S3  <- data.frame(Area=region.names,Prefectures=region.prefs,n_dates=n.dates,n_sites=n.sites)

write.table(table.S3,file=here('manuscript','supplementary_tables','table_S3.csv'),col.names=c('Area','Prefectures','Number of Dates','Number of Sites'),sep=',',row.names=FALSE)



# Table S4 posterior estimates for nu ----
load(here("results","phase_model_a.RData"))
load(here("results","phase_model_b.RData"))
out.comb.unif.modela  <- do.call(rbind,out.unif.model_a)
out.comb.unif.modelb  <- do.call(rbind,out.unif.model_b)
post.nu.modela  <- out.comb.unif.modela[,paste0('a[',1:8,']')] |> round()
post.nu.modelb  <- out.comb.unif.modelb[,paste0('a[',1:8,']')] |> round()
hpdi.modela  <- apply(post.nu.modela,2,function(x){HPDinterval(as.mcmc(x),prob = .90)}) 
hpdi.modelb  <- apply(post.nu.modelb,2,function(x){HPDinterval(as.mcmc(x),prob = .90)}) 
med.modela  <- apply(post.nu.modela,2,median)
med.modelb  <- apply(post.nu.modelb,2,median)


foo  <- function(x)
{
	x = BPtoBCAD(x)
	ifelse(x<0,paste(abs(x),'BC'),paste(x,'AD'))
}
models  <- rep(c('Model a','Model b'),each=8) 
area  <- rep(as.character(as.roman(1:8)),2)
meds  <- c(foo(med.modela),foo(med.modelb))
hi90  <- c(foo(hpdi.modela[1,]),foo(hpdi.modelb[1,]))
lo90  <- c(foo(hpdi.modela[2,]),foo(hpdi.modelb[2,]))
rhat  <- c(rhat.unif.model_a$psrf[1:8,1],rhat.unif.model_b$psrf[1:8,1]) |> round(digits=3)
ess  <- c(ess.unif.model_a[1:8],ess.unif.model_b[1:8]) |> round()
table.S4  = data.frame(models,area,meds,lo90,hi90,rhat,ess)
write.table(table.S4,file=here('manuscript','supplementary_tables','table_S4.csv'),col.names=c('Model','Area','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)
