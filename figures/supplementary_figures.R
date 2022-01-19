# Load Libraries and spatial data ----
library(here)
library(truncnorm)
library(cascsim)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(maptools)
library(sf)
library(rgeos)
library(viridis)
library(latex2exp)
library(gridExtra)
library(diagram)
library(quantreg)
source(here('src','orderPPlot.R'))
source(here('src','diffplot.R'))
source(here('src','gpqrSim.R'))

# Figure S1 (Impact of Hallstatt Plateau) ----
# Load Observed Data
load(here('data','c14rice.RData'))
# Load quantile regression results
load(here('results','quantreg_res.RData'))
## Compute Fitted Model Confidence Intervals:

# rq and median calibrated date
rq.ci <- predict.rq(fit.rq,newdata=data.frame(dist_org=0:1300),interval='confidence')

# Bayesian model 
qr.ch1 <- quantreg_sample[[1]]
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg,beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta',colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg,2,median)
post.quantreg <- apply(post.alpha.beta.quantreg,1,function(x){x[1]-x[2]*0:1300})
post.ci <- t(apply(post.quantreg,1,quantile,c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('figures','figureS1.pdf'),width=8.5,height=7)
plot(NULL,xlim=c(0,1300),ylim=c(3200,900),axes=F,xlab='Distance from Nabatake Site (in km)',ylab='Cal BP')
rect(xleft=-100,xright=1400,ybottom=2720,ytop=2350,col=col.alpha('grey',0.2),border=NA)
abline(h=2720,lty=4)
abline(h=2350,lty=4)
axis(1) 
axis(2)
axis(4,at=BCADtoBP(c(-1000,-600,-200,200,600,1000)),labels=c('1000BC','600BC','200BC','200AD','600AD','1000AD'))
points(constants$dist_org,SiteInfo$Earliest)
points(constants$dist_org,post.theta.med,pch=20)
for (i in 1:nrow(SiteInfo))
{
   lines(rep(constants$dist_org[i],2),c(SiteInfo$Earliest[i],post.theta.med[i]),lty=2)
}
lines(0:1300,rq.ci[,1],lty=1,lwd=2,col='blue')
polygon(x=c(0:1300,1300:0),c(rq.ci[,2],rev(rq.ci[,3])),col=col.alpha('lightblue',0.4),border=NA)

lines(0:1300,post.ci[,2],lty=1,lwd=2,col='indianred')
polygon(x=c(0:1300,1300:0),c(post.ci[,1],rev(post.ci[,3])),col=col.alpha('indianred',0.2),border=NA)

text(x=245,y=2400,labels='Hallstat Plateau')
legend('bottomright',legend=c('Median Calibrated Date',TeX('Median Posterior $\\theta$'),'Quantile Regression on Median Dates','Bayesian Quantile Regression with Measurement Error'),pch=c(1,20,NA,NA),lwd=c(NA,NA,2,2),col=c(1,1,'blue','indianred'))
box()
dev.off()

# Figure S2 (Simulated dispersal rate and local deviations) ----
# Load Spatial Data 
win.sf  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')
# Load Tactical Simulation
load(here('results','tactical_sim_res.RData'))
# Load Tactical Simulation GPQR
load(here('results','res_tactical.RData'))
# Extract and Compute Posteriors
tactical_fitted  <- rbind(res[[1]],res[[2]],res[[3]])
fitted_s  <- tactical_fitted[,grep('s\\[',colnames(tactical_fitted))]
fitted_beta1  <- tactical_fitted[,'beta1']
fitted_slope <- fitted_s - fitted_beta1
fitted_rate <- -1/fitted_slope
# Assign to sim.sites posterior mean of each relevant parameters
sim.sites$fitted_s <- apply(fitted_s,2,mean)
sim.sites$fitted_rate <- apply(fitted_rate,2,mean)

# Plot true s
fS2a  <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=s),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='a',fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

# Plot predicted s
fS2b  <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=fitted_s),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='b',fill='Predicted s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

# Plot true local rate
fS2c <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=rate),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='c',fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="magma",limits=c(1,3.5)) +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

# Plot predicted local rare
fS2d <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=fitted_rate),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='c',fill='Predicted \n Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="magma",limits=c(1,3.5)) +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


pdf(file=here('figures','figureS2.pdf'),width=8,height=8)
grid.arrange(fS2a,fS2b,fS2c,fS2d,ncol=2,nrow=2)
dev.off()

# Figure S3 (Posterior vs True values of s) ----

# Load Tactical Simulation GPQR
load(here('results','res_tactical.RData'))
# Load Tactical Simulation GPQR
# Extract and Compute Posteriors
tactical_fitted  <- rbind(res[[1]],res[[2]],res[[3]])
fitted_s  <- tactical_fitted[,grep('s\\[',colnames(tactical_fitted))]
fitted_alpha  <- tactical_fitted[,grep('alpha\\[',colnames(tactical_fitted))]
fitted_beta1  <- tactical_fitted[,'beta1']
fitted_slope <- fitted_s - fitted_beta1
fitted_rate <- -1/fitted_slope

pdf(file=here('figures','figureS3.pdf'),width=8,height=8)
par(mfrow=c(2,2))

# panel a : s
lo_s <- apply(fitted_s,2,quantile,0.025)
hi_s  <- apply(fitted_s,2,quantile,0.975)
rr.y = c(min(lo_s),max(hi_s))
rr.x = range(sim.SiteInfo$s)
col = rep('darkgrey',sim.constants$Nsites)
col[which(sim.SiteInfo$s>hi_s | sim.SiteInfo$s <lo_s)] = 'darkorange'
plot(sim.SiteInfo$s,apply(fitted_s,2,median),pch=20,xlim=rr.x,ylim=rr.y,xlab=TeX('Simulated $s_k$'),ylab=TeX('Predicted $s_k$'),main='a',col=col)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.SiteInfo$s[i],2),c(lo_s[i],hi_s[i]),lwd=0.5,col=col[i])
}
abline(a=0,b=1,lty=2,col=1,lwd=2)

# panel b : alpha
lo_alpha <- apply(fitted_alpha,2,quantile,0.025)
hi_alpha  <- apply(fitted_alpha,2,quantile,0.975)
rr.y = c(min(lo_alpha),max(hi_alpha))
rr.x = range(sim.SiteInfo$alpha)
col = rep('darkgrey',sim.constants$Nsites)
col[which(sim.SiteInfo$alpha>hi_alpha | sim.SiteInfo$alpha <lo_alpha)] = 'darkorange'

plot(sim.SiteInfo$alpha,apply(fitted_alpha,2,median),pch=20,xlim=rev(rr.x),ylim=rev(rr.y),xlab=TeX('Simulated $\\alpha_k$'),ylab=TeX('Predicted $\\alpha_k$'),main='b',col=col)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.SiteInfo$alpha[i],2),c(lo_alpha[i],hi_alpha[i]),lwd=0.5,col=col[i])
}
abline(a=0,b=1,lty=2,col=1,lwd=2)


#panel c: slope (s+beta1)
sim.slope = sim.SiteInfo$s+sim.constants$beta1
lo_slope <- apply(fitted_slope,2,quantile,0.025)
hi_slope  <- apply(fitted_slope,2,quantile,0.975)
rr.y = c(min(lo_slope),max(hi_slope))
rr.x = range(sim.slope)
col = rep('darkgrey',sim.constants$Nsites)
col[which(sim.slope>hi_slope | sim.slope <lo_slope)] = 'darkorange'

plot(sim.slope,apply(fitted_slope,2,median),pch=20,xlim=rr.x,ylim=rr.y,xlab='Simulated Slope',ylab='Predicted Slope',main='c',col=col)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.slope[i],2),c(lo_slope[i],hi_slope[i]),lwd=0.5,col=col[i])
}
abline(a=0,b=1,lty=2,col=1,lwd=2)


#panel d: rate (-1/slope)
sim.rate = -1/sim.slope
lo_rate <- apply(fitted_rate,2,quantile,0.025)
hi_rate  <- apply(fitted_rate,2,quantile,0.975)
rr.y = c(min(lo_rate),max(hi_rate))
rr.x = range(sim.rate)
col = rep('darkgrey',sim.constants$Nsites)
col[which(sim.rate>hi_rate | sim.rate <lo_rate)] = 'darkorange'

plot(sim.rate,apply(fitted_rate,2,median),pch=20,xlim=rr.x,ylim=rr.y,xlab='Simulated Dispersal Rate',ylab='Predicted Dispersal Rate',main='d',col=col)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.rate[i],2),c(lo_rate[i],hi_rate[i]),lwd=0.5,col=col[i])
}
abline(a=0,b=1,lty=2,col=1,lwd=2)


dev.off()

# Figure S4 (Posterior vs True values of beta0,beta1,rho,etasq) ----
# Load Tactical Simulation
load(here('results','tactical_sim_res.RData'))
# Load GPQR on tactical simulation
load(here('results','res_tactical.RData'))

tactical_fitted  <- rbind(res[[1]],res[[2]],res[[3]])
fitted_beta0  <- tactical_fitted[,'beta0']
fitted_beta1  <- tactical_fitted[,'beta1']
fitted_etasq  <- tactical_fitted[,'etasq']
fitted_rhosq  <- tactical_fitted[,'rhosq']
fitted_omega  <- tactical_fitted[,'omega']
fitted_phi  <- tactical_fitted[,'phi']
pdf(file=here('figures','figureS4.pdf'),width=6,height=9)

par(mfrow=c(3,2))
postHPDplot(fitted_beta0,main=TeX('Posterior $\\beta_0$'),HPD = 0.95)
abline(v=qnorm(0.90,mean=sim.constants$beta0,sd=sim.constants$sigma),lty=2)

postHPDplot(1/fitted_beta1,main=TeX('Posterior $\\beta_1$'),HPD = 0.95)
abline(v=-1/sim.constants$beta1,lty=2)

postHPDplot(fitted_etasq,main=TeX('Posterior $\\eta^2$'),HPD = 0.95)
abline(v=sim.constants$etasq,lty=2)

postHPDplot(fitted_rhosq,main=TeX('Posterior $\\rho^2$'),HPD = 0.95,rnd = 5)
abline(v=sim.constants$rhosq,lty=2)

postHPDplot(fitted_omega,main=TeX('Posterior $\\omega$'),HPD = 0.95)
abline(v=sim.constants$omega,lty=2)

postHPDplot(fitted_phi,main=TeX('Posterior $\\phi$'),HPD = 0.95)
abline(v=sim.constants$phi,lty=2)

dev.off()

# Figure S5 (Prior Predictive Check beta0 and beta1) ----
nsim <- 500
beta0.prior <- rnorm(nsim,mean=3000,sd=200)
# beta1.prior  <- rnorm(nsim,mean=0,sd=3)
beta1.prior  <- -rexp(nsim,rate=1)


pdf(file=here('figures','figureS5.pdf'),width=6,height=5)
plot(NULL,xlim=c(0,1300),ylim=c(3400,1500),type='n',xlab='Distance (km)',ylab='Cal BP',axes=F)
axis(1)
axis(2,at=seq(3400,1600,-400))
axis(4,at=BCADtoBP(c(-1400,-1000,-600,-200,200,600)),labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'))
box()
for(i in 1:nsim){abline(a=beta0.prior[i],b=beta1.prior[i],col=rgb(0,0,0,0.1))}
dev.off()

# Figure S6 (Intepretation of etasq and rho) ----
# spatial window
sf_japan <- ne_states(country = "japan") |> subset(!name_vi %in%  c("Okinawa","Hokkaido"))
sampling.win <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan)))
sampling.win <- disaggregate(sampling.win) 
sampling.win  <- sampling.win[order(raster::area(sampling.win),decreasing=TRUE)[1:3]]
win.sf  <- as(sampling.win,'sf')
win = sampling.win
# fixed params
n = 300
origin.point = c(129.959,33.4485)
beta0 = 3000
beta1 = 1
sigma = 100
seed = 144231
# sweep params
etasq = c(0.01,0.05,0.15)
rho = c(10,100,500)
cov.param = expand.grid(etasq=etasq,rho=rho)
tmp  <- vector('list',length=nrow(cov.param))
out  <- vector('list',length=nrow(cov.param))


for (i in 1:nrow(cov.param))
{
	tmp[[i]]  <- gpqrSim(win=win,n=n,beta0=beta0,beta1=beta1,sigma=sigma,origin.point=origin.point,etasq=cov.param$etasq[i],rho=cov.param$rho[i],seed=seed)

	out[[i]] <- ggplot() +
		geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
		geom_sf(data=tmp[[i]],mapping = aes(fill=rate),pch=21,col='darkgrey',size=1.5) + 
		xlim(129,143) + 
		ylim(31,42) +
# 		labs(title=paste0('etasq=',cov.param$etasq[i],' rho=',cov.param$rho[i]),fill='Dispersal Rate \n (km/yr)') + 
 		labs(title=TeX(sprintf("$\\eta^2 = %g \\, \\rho = %g$",cov.param$etasq[i],cov.param$rho[i])),fill='Dispersal Rate \n (km/yr)') + 
 		scale_fill_viridis(option = 'turbo',limits = c(0.7, 3), oob = scales::squish) +
#  		scale_fill_viridis() +
		theme(plot.title = element_text(hjust = 0.5,size=11), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=7),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

}


pdf(file=here('figures','figureS6.pdf'),width=8,height=8)
grid.arrange(grobs=out,nrow=3,ncol=3)
dev.off()

# Figure S7 (Prior Predictive Check etasq and rhosq) ----
nsim  <- 500
etasq.prior  <- rexp(nsim,10)
rho.prior  <- rgamma(nsim,10,(10-1)/200)

pdf(file=here('figures','figureS7.pdf'),width=6,height=5)
plot(NULL,xlab='Distance (km)',ylab='Covariance',xlim=c(0,1000),ylim=c(0,1))
for (i in 1:nsim)
{
	curve(etasq.prior[i]*exp(-0.5*(x/rho.prior[i])^2),add=TRUE,from=0,to=1000,col=rgb(0,0,0,0.1))
}
dev.off()

# Figure S8 (Traceplot of beta0, beta1, rhosq, and etasq for tau=0.9) ----

load(here('results','gpqr_tau90.RData'))

pdf(file=here('figures','figureS8.pdf'),width=8,height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau90[,'beta0'],main=TeX('$\\beta_0$'),smooth=TRUE)
traceplot(gpqr_tau90[,'beta1'],main=TeX('$\\beta_1$'),smooth=TRUE)
traceplot(gpqr_tau90[,'rho'],main=TeX('$\\rho$'),smooth=TRUE)
traceplot(gpqr_tau90[,'etasq'],main=TeX('$\\eta^2$'),smooth=TRUE)
dev.off()

# Figure S9 (Traceplot of beta0, beta1, rhosq, and etasq for tau=0.99) ----

load(here('results','gpqr_tau99.RData'))

pdf(file=here('figures','figureS9.pdf'),width=8,height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau99[,'beta0'],main=TeX('$\\beta_0$'),smooth=TRUE)
traceplot(gpqr_tau99[,'beta1'],main=TeX('$\\beta_1$'),smooth=TRUE)
traceplot(gpqr_tau99[,'rho'],main=TeX('$\\rho$'),smooth=TRUE)
traceplot(gpqr_tau99[,'etasq'],main=TeX('$\\eta^2$'),smooth=TRUE)
dev.off()

# Table S1 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.9) ----

gpqr.tau90.comb  <- do.call(rbind,gpqr_tau90)
params = c('beta0','beta1','rho','etasq')
meds = apply(gpqr.tau90.comb[,params],2,median)
lo90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]})
hi90 = apply(gpqr.tau90.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]})
rhats = gelman.diag(gpqr_tau90)$psrf[params,1]
ess = effectiveSize(gpqr_tau90)[params]
table.S1 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S1,file=here('tables','table_S1.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)


# Table S2 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.99) ----

gpqr.tau99.comb  <- do.call(rbind,gpqr_tau99)
params = c('beta0','beta1','rho','etasq')
meds = apply(gpqr.tau99.comb[,params],2,median)
lo90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]})
hi90 = apply(gpqr.tau99.comb[,params],2,function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]})
rhats = gelman.diag(gpqr_tau99)$psrf[params,1]
ess = effectiveSize(gpqr_tau99)[params]
table.S2 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S2,file=here('tables','table_S1.csv'),col.names=c('Parameter','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)

# Figure S10 (Marginal and Joint posteriors of beta0, beta1, rhosq, etasq for tau = 0.9) ----

pdf(file=here('figures','figureS10.pdf'),width=8,height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau90.comb[,'beta0'],main=TeX('$\\beta_0$'),xlab='Cal BP',ylab='')
postHPDplot(gpqr.tau90.comb[,'beta1'],main=TeX('$\\beta_1$'),xlab='',ylab='')
postHPDplot(gpqr.tau90.comb[,'rho'],main=TeX('$\\rho$'),xlab='km',ylab='')
postHPDplot(gpqr.tau90.comb[,'etasq'],main=TeX('$\\eta^2$'),xlab='',ylab='')
dev.off()

# Figure S11 (Marginal and Joint posteriors of beta0, beta1, rhosq, etasq for tau = 0.99) ----

pdf(file=here('figures','figureS11.pdf'),width=8,height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau99.comb[,'beta0'],main=TeX('$\\beta_0$'),xlab='Cal BP',ylab='')
postHPDplot(gpqr.tau99.comb[,'beta1'],main=TeX('$\\beta_1$'),xlab='',ylab='')
postHPDplot(gpqr.tau99.comb[,'rho'],main=TeX('$\\rho$'),xlab='km',ylab='')
postHPDplot(gpqr.tau99.comb[,'etasq'],main=TeX('$\\eta^2$'),xlab='',ylab='')
dev.off()

# Figure S12 Tactical Simulation Posterior Predictive Check for nu and upsilon ----



# Table S3 posterior estimates for nu and upsilon ----
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
table.S3  = data.frame(models,area,meds,lo90,hi90,rhat,ess)
write.table(table.S3,file=here('tables','table_S3.csv'),col.names=c('Model','Area','Median Posterior','90% HPDI (low)','90% HPDI (high)','Rhat','ESS'),sep=',',row.names=FALSE)

# Figure S13 Marginal Posterior Distribution of nu, model 0 ----
model0.long  <- data.frame(value=as.numeric(post.nu.model0),Area = rep(as.character(as.roman(1:8)),each=nrow(post.nu.model0)))

pdf(file=here('figures','figureS13.pdf'),height=10,width=7)
ggplot(model0.long, aes(x = value, y = Area,fill='lighblue')) + 
	geom_density_ridges() +
	scale_x_reverse(limits=c(3300,1800),breaks=BCADtoBP(c(-1200,-1000,-800,-600,-400,-200,1)),labels=c(1200,1000,800,600,400,200,1)) +
	scale_fill_manual(values='lightblue') +
	theme(legend.position = "none") +
	xlab('BC')
dev.off()

# Figure S14 Probability Matrix of nu, model 0 ----
pdf(file=here('figures','figureS14.pdf'),width=7,height=7.5)
orderPPlot(post.nu.model0,name.vec=paste("Area",as.character(as.roman(1:8))))
dev.off()

# Figure S15 Difference Matrix plot of nu, model 0 ----
pdf(file=here('figures','figureS15.pdf'),width=16,height=11)
mat <- cbind(c(1,9:15),c(37,2,16:21),c(rep(37,2),3,22:26),c(rep(37,3),4,27:30),c(rep(37,4),5,31:33),c(rep(37,5),6,34:35),c(rep(37,6),7,36),c(rep(37,7),8))
layout(mat)
par(mar=c(0,0,0,0))
for (i in 1:8){plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F);text(0.5,0.5,paste('Area',as.character(as.roman(i))),cex=3)}
par(mar=c(3,0,0,1))
for (i in 1:8){
	for (j in 1:8){
		if (i < j)
		{
			diffDens(post.nu.model0[,i],post.nu.model0[,j],xlim=c(-1200,1200))
		}
	}
}
dev.off()



# Figure S16 Marginal Posterior Distribution of nu and upsilon, model 1 ----
model1.long  <- data.frame(value=as.numeric(post.nu.model1),Area = rep(as.character(as.roman(1:8)),each=nrow(post.nu.model1)))

pdf(file=here('figures','figureS16.pdf'),height=10,width=7)
ggplot(model1.long, aes(x = value, y = Area,fill='lighblue')) + 
	geom_density_ridges() +
	scale_x_reverse(limits=c(3300,1800),breaks=BCADtoBP(c(-1200,-1000,-800,-600,-400,-200,1)),labels=c(1200,1000,800,600,400,200,1)) +
	scale_fill_manual(values='lightblue') +
	theme(legend.position = "none") +
	xlab('BC')
dev.off()

# Figure S17 Marginal Posterior Distribution of nu and upsilon, model 2 ----
model2.long  <- data.frame(value=as.numeric(post.nu.model2),Area = rep(as.character(as.roman(1:8)),each=nrow(post.nu.model2)))

pdf(file=here('figures','figureS17.pdf'),height=10,width=7)
ggplot(model2.long, aes(x = value, y = Area,fill='lighblue')) + 
	geom_density_ridges() +
	scale_x_reverse(limits=c(3300,1800),breaks=BCADtoBP(c(-1200,-1000,-800,-600,-400,-200,1)),labels=c(1200,1000,800,600,400,200,1)) +
	scale_fill_manual(values='lightblue') +
	theme(legend.position = "none") +
	xlab('BC')
dev.off()




