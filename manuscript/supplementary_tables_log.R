# Load Libraries and spatial data ----
library(here)
library(rcarbon) 
library(dplyr)
library(coda)
# Table S1 Rice Radiocarbon Database ----
table.S1  <- read.csv(here('data','R14CDB.csv'))
ii  <- which(table.S1$C14Age>0)
x  <- calibrate(table.S1$C14Age[ii],table.S1$C14Error[ii],calCurve='intcal20')
table.S1$medianDate  <- table.S1$sigmaRange  <-  NA
table.S1$medianDate[ii] <- paste(abs(BPtoBCAD(medCal(x))),ifelse(BPtoBCAD(medCal(x))<0,'BC','AD'))
table.S1$sigmaRange[ii] <- hpdi(x,calendar='BCAD',credMass=0.9545,asList=FALSE)
table.S1$used = table.S1$UseForAnalyses
table.S1$used[which(table.S1$Prefecture%in%c('Hokkaido','Okinawa'))]  <- FALSE
table.S1 <- select(table.S1,ID,LabCode,C14Age,C14Error,medianDate,sigmaRange,Method,SiteName_jp,SiteName_en,Latitude,Longitude,Context,Prefecture,Region,PublicationYear,Reference,Notes,used)

write.table(table.S1,file=here('manuscript','supplementary_tables','table_S1.csv'),col.names=c('ID','LabCode','C14Age','C14AgeError','MedianDate','Two-Sigma Range','Method','SiteName_jp','SiteName_en','Laitude','Longitude','Context','Prefecture','Region','PublicationYear','Reference','Notes','USED'),sep=',',row.names=FALSE)

# Table S2 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.9) ----
load(here('results','gpqr_tau90.RData'))
gpqr.tau90.comb  <- do.call(rbind,gpqr_tau90)
params = c('beta0','beta1','rho','etasq')
meds = apply(gpqr.tau90.comb[,params],2,median)
lo90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]})
hi90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]})
rhats = gelman.diag(gpqr_tau90)$psrf[params,1]
ess = effectiveSize(gpqr_tau90)[params]
table.S2 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S2,file=here('manuscript','supplementary_tables','table_S1.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)


# Table S3 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.99) ----
load(here('results','gpqr_tau99.RData'))
gpqr.tau99.comb  <- do.call(rbind,gpqr_tau99)
params = c('beta0','beta1','rho','etasq')
meds = apply(gpqr.tau99.comb[,params],2,median)
lo90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]})
hi90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]})
rhats = gelman.diag(gpqr_tau99)$psrf[params,1]
ess = effectiveSize(gpqr_tau99)[params]
table.S3 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S3,file=here('manuscript','supplementary_tables','table_S2.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)

# Table S4 Summary Data Per Region ----
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
table.S4  <- data.frame(Area=region.names,Prefectures=region.prefs,n_dates=n.dates,n_sites=n.sites)

write.table(table.S4,file=here('manuscript','supplementary_tables','table_S3.csv'),col.names=c('Area','Prefectures','Number of Dates','Number of Sites'),sep=',',row.names=FALSE)



# Table S5 posterior estimates for nu and upsilon ----
load(here("results","phase_model0.RData"))
load(here("results","phase_model1.RData"))
load(here("results","phase_model2.RData"))
out.comb.unif.model0  <- do.call(rbind,out.unif.model0)
out.comb.unif.model1  <- do.call(rbind,out.unif.model1)
out.comb.unif.model2  <- do.call(rbind,out.unif.model2)
post.nu.model0  <- out.comb.unif.model0[,paste0('a[',1:8,']')] |> round()
post.nu.model1  <- out.comb.unif.model1[,paste0('a[',1:8,']')] |> round()
post.nu.model2  <- out.comb.unif.model2[,paste0('a[',1:8,']')] |> round()
hpdi.model0  <- apply(post.nu.model0,2,function(x){HPDinterval(as.mcmc(x),prob = .90)}) 
hpdi.model1  <- apply(post.nu.model1,2,function(x){HPDinterval(as.mcmc(x),prob = .90)}) 
hpdi.model2  <- apply(post.nu.model2,2,function(x){HPDinterval(as.mcmc(x),prob = .90)}) 
med.model0  <- apply(post.nu.model0,2,median)
med.model1  <- apply(post.nu.model1,2,median)
med.model2  <- apply(post.nu.model2,2,median)


foo  <- function(x)
{
	x = BPtoBCAD(x)
	ifelse(x<0,paste(abs(x),'BC'),paste(x,'AD'))
}
models  <- rep(c('Model 0','Model 1','Model 2'),each=8) 
area  <- rep(as.character(as.roman(1:8)),3)
meds  <- c(foo(med.model0),foo(med.model1),foo(med.model2))
hi90  <- c(foo(hpdi.model0[1,]),foo(hpdi.model1[1,]),foo(hpdi.model2[1,]))
lo90  <- c(foo(hpdi.model0[2,]),foo(hpdi.model1[2,]),foo(hpdi.model2[2,]))
rhat  <- c(rhat.unif.model0$psrf[1:8,1],rhat.unif.model1$psrf[1:8,1],rhat.unif.model2$psrf[1:8,1]) |> round(digits=3)
ess  <- c(ess.unif.model0[1:8],ess.unif.model1[1:8],ess.unif.model2[1:8]) |> round()
table.S5  = data.frame(models,area,meds,lo90,hi90,rhat,ess)
write.table(table.S5,file=here('manuscript','supplementary_tables','table_S4.csv'),col.names=c('Model','Area','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)
