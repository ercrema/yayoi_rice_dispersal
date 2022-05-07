# Load Libraries and spatial data ----
library(here)
library(truncnorm)
library(cascsim)
library(corrplot)
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
library(coda)
source(here('src','orderPPlot.R'))
source(here('src','diffplot.R'))
source(here('src','gpqrSim.R'))


# Figure S1 (Impact of Measurement Error and Calibration Plateau)----
d.sim = c(0,200,400,600)
d.dates = BCADtoBP(c(-760,-690,-620,-550))
pdf(here('manuscript','supplementary_figures','figureS1.pdf'),width=5,height=4)
par(mar=c(5,4,2,4))
plot(d.sim,d.dates,pch=20,ylim=c(2800,2300),xlab='Distance',ylab='',type='b')
axis(4,at=BCADtoBP(c(-800,-700,-600,-500,-400)),labels=c(800,700,600,500,400))
mtext(side=2,line=3,'cal BP')
mtext(side=4,line=3,'BC')
d.caldates = medCal(calibrate(uncalibrate(d.dates)$ccCRA,rep(30,4)))
lines(d.sim,d.caldates,pch=20,col='red',type='b')
legend('topleft',legend=c('True Date','Median Calibrated Date'),pch=c(20,20),col=c(1,2))
dev.off()


# Figure S2 (Bayesian Quantile Regression Model with Measurement Error) ----
# Load Observed Data
load(here('data','c14rice.RData'))
# Load quantile regression results
load(here('results','quantreg_res.RData'))
## Compute Fitted Model Confidence Intervals:

# rq and median calibrated date
rq.ci <- predict.rq(fit.rq,newdata=data.frame(dist_org=0:1300),interval='confidence')

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
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

pdf(file=here('manuscript','supplementary_figures','figureS2.pdf'),width=8.5,height=7)
plot(NULL,xlim=c(0,1300),ylim=c(3200,900),axes=F,xlab='Distance from Ukikunden Site (in km)',ylab='Cal BP')
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

text(x=245,y=2400,labels='Calibration Plateau')
legend('bottomright',legend=c('Median Calibrated Date',TeX('Median Posterior $\\theta$'),'Quantile Regression on Median Dates','Bayesian Quantile Regression with Measurement Error'),pch=c(1,20,NA,NA),lwd=c(NA,NA,2,2),col=c(1,1,'blue','indianred'),cex=0.8)
box()
dev.off()


# Figure S3 (Posterior Dispersal Rate of non-spatial quantile regression) ----
pdf(here('manuscript','supplementary_figures','figureS3.pdf'),height=5,width=5.5)
postHPDplot(1/post.beta.quantreg,xlab='km/year',ylab='Probability Density',prob=.90,main=TeX('Posterior of $1/\\beta_1$'))
dev.off()


# Figure S4 (Variance Covariance Function Parameters) ----
pdf(here('manuscript','supplementary_figures','figureS4.pdf'),width=5,height=5)
etasq = 0.05
rho = 150
curve(etasq*exp(-0.5*(x/rho)^2),from=0,to=500,xlab=TeX('$d_{i,j}$'),ylab=TeX('$cov(i,j)$'))
abline(h=0.05,lty=2)
lines(x=c(rho,rho),y=c(-0.01,etasq*exp(-0.5)),lty=2)
lines(x=c(-20,rho),y=c(etasq*exp(-0.5),etasq*exp(-0.5)),lty=2)
text(400,0.046,label=TeX('$\\eta^2 = 0.05$'))
text(75,etasq*exp(-0.5)-0.004,label=TeX('$\\eta^2e^{-0.5}$'))
text(200,0.01,label=TeX('$\\rho = 150$'))
dev.off()

# Figure S5 (Impact of rho and etasq on variability in dispersal rate) ----
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


pdf(file=here('manuscript','supplementary_figures','figureS5.pdf'),width=8,height=8)
grid.arrange(grobs=out,nrow=3,ncol=3)
dev.off()


# Figure S6 (Relationship between slope and rate of dispersal)  ----
pdf(here('manuscript','supplementary_figures','figureS6.pdf'),width=5,height=5)
par(mar=c(5,6,2,2))
curve(-1/x,from=-2.5,to=-0.05,xlab=TeX('$Slope: \\, s - \\beta_1$'),ylab=TeX('$Speed: \\, \\frac{-1}{s-\\beta_1}$'),axes=TRUE)
dev.off()


# Figure S7 (Prior Predictive Check beta0, beta1, s) ----
nsim <- 5000
beta0.prior <- rnorm(nsim,mean=3000,sd=200)
beta1.prior  <- rexp(nsim,rate=1)
s.prior  <- rnorm(nsim,mean=0,sd=sqrt(rexp(nsim,rate=20)))
slope  <- s.prior - beta1.prior
beta0.prior  <- beta0.prior[which((-1/slope)>0)]
slope  <- slope[which((-1/slope)>0)]
nsim2  <- length(slope)
dists  <- -100:1400
slope.mat = matrix(NA,nrow=nsim2,ncol=length(dists))
for (i in 1:nsim2)
{
	slope.mat[i,] <- beta0.prior[i] + slope[i]*c(dists)	
}

pdf(file=here('manuscript','supplementary_figures','figureS7.pdf'),width=6,height=6)
plot(NULL,xlim=c(0,1300),ylim=c(3400,1500),type='n',xlab='Distance (km)',ylab='Cal BP',axes=F)
axis(1)
axis(2,at=seq(3400,1600,-400))
axis(4,at=BCADtoBP(c(-1400,-1000,-600,-200,200,600)),labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'))
box()
polygon(x=c(dists,rev(dists)),y=c(apply(slope.mat,2,quantile,prob=0.025),rev(apply(slope.mat,2,quantile,prob=0.975))),border=NA,col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists,rev(dists)),y=c(apply(slope.mat,2,quantile,prob=0.25),rev(apply(slope.mat,2,quantile,prob=0.75))),border=NA,col=rgb(0.25,0.41,0.88,0.5))
abline(a=3000,b=-1,lty=2)
text(x=1200,y=1900,label='1km/yr')
abline(a=3000,b=-1/0.5,lty=2)
text(x=700,y=1800,label='0.5km/yr')
abline(a=3000,b=-1/5,lty=2)
text(x=1250,y=2800,label='5km/yr')
legend('bottomright',legend=c('50% percentile range','95% percentile range'),fill=c(rgb(0.67,0.84,0.9,0.5),rgb(0.25,0.41,0.88,0.5)))
dev.off()


# Figure S8 (Prior Predictive Check etasq and rho) ----
nsim  <- 1000
etasq.prior  <- rexp(nsim,20)
rho.prior  <- rtgamma(nsim,10,(10-1)/150,1,1350)
cov.mat = matrix(NA,nrow=nsim,ncol=length(0:1000))
for (i in 1:nsim)
{
 cov.mat[i,] = etasq.prior[i]*exp(-0.5*(0:1000/rho.prior[i])^2)
}

pdf(file=here('manuscript','supplementary_figures','figureS8.pdf'),width=6,height=5)
plot(NULL,xlab='Distance (km)',ylab='Covariance',xlim=c(0,1000),ylim=c(0,0.2))
polygon(c(0:1000,1000:0),c(apply(cov.mat,2,quantile,0.025),rev(apply(cov.mat,2,quantile,0.975))),border=NA,col=rgb(0.67,0.84,0.9,0.5))
polygon(c(0:1000,1000:0),c(apply(cov.mat,2,quantile,0.5),rev(apply(cov.mat,2,quantile,0.75))),border=NA,col=rgb(0.25,0.41,0.88,0.5))
legend('topright',legend=c('50% percentile range','95% percentile range'),fill=c(rgb(0.67,0.84,0.9,0.5),rgb(0.25,0.41,0.88,0.5)))
dev.off()


# Figure S9 (Tactical Simulation) ----
load(here('data','tactical_sim_gpqr.RData'))
load(here('results','gpqr_tactsim.RData'))
gpqr_tactsim_post  <- do.call(rbind,gpqr_tactsim)
tactsim_post_s  <- gpqr_tactsim_post[,paste0('s[',1:nrow(sim.sites),']')]
tactsim_post_beta1  <- gpqr_tactsim_post[,'beta1']
tactsim_post_rate  <-  -1/(tactsim_post_s-tactsim_post_beta1)
sim.sites$pred.rate  <- apply(tactsim_post_rate,2,median)


s9a  <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=rate),pch=21,col='darkgrey',size=3) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="turbo",limits=c(0,3.5),oob = scales::squish) +
	annotate("text", x = 140, y = 33, label = TeX('$\\beta_0 = 3000$')) +
	annotate("text", x = 140, y = 32.5, label = TeX('$\\beta_1 = 0.6$')) +
	annotate("text", x = 140, y = 32, label = TeX('$\\eta^2 = 0.08$')) +
	annotate("text", x = 140, y = 31.5, label = TeX('$\\rho = 150$')) +
	ggtitle('Simulated Dispersal Rates') +
	theme(legend.position = c(0.2, 0.8),legend.background=element_rect(fill = alpha("white",0.5)),axis.title.x=element_blank(),axis.title.y=element_blank())


s9b  <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=pred.rate),pch=21,col='darkgrey',size=3) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="turbo",limits=c(0,3.5),oob = scales::squish) +
	ggtitle('Predicted Dispersal Rates') +
	theme(legend.position = c(0.2, 0.8),legend.background=element_rect(fill = alpha("white",0.5)),axis.title.x=element_blank(),axis.title.y=element_blank())


pdf(here('manuscript','supplementary_figures','figureS9.pdf'),width=10,height=7)
grid.arrange(s9a,s9b,ncol=2)
dev.off()




# Figure S10 (Posterior vs True values of s for Tactical Simulation) ----
load(here('results','gpqr_tactsim.RData'))
gpqr_tactsim_post  <- do.call(rbind,gpqr_tactsim)
tactsim_post_s  <- gpqr_tactsim_post[,paste0('s[',1:nrow(sim.sites),']')]
tactsim_post_s_med  <- apply(tactsim_post_s,2,median)
tactsim_post_s_lo  <- apply(tactsim_post_s,2,quantile,0.025)
tactsim_post_s_hi  <- apply(tactsim_post_s,2,quantile,0.975)
rr = c(min(c(tactsim_post_s_lo,sim.sites$s)),max(c(tactsim_post_s_hi,sim.sites$s)))

pdf(here('manuscript','supplementary_figures','figureS10.pdf'),height=6,width=6)
plot(NULL,xlim=rr,ylim=rr,xlab='Simulated s',ylab='Predicted s')
points(sim.sites$s,tactsim_post_s_med,pch=20)
for (i in 1:nrow(sim.sites))
{
	lines(x=c(sim.sites$s[i],sim.sites$s[i]),y=c(tactsim_post_s_lo[i],tactsim_post_s_hi[i]))
}
abline(a=0,b=1,lty=2,col='red')
dev.off()

# Figure S11 (Posterior vs True values of beta0,beta1,rho,etasq for Tactical Simulation) ----
tactsim_post_beta0  <- gpqr_tactsim_post[,'beta0']
tactsim_post_beta1 <- gpqr_tactsim_post[,'beta1']
tactsim_post_etasq  <- gpqr_tactsim_post[,'etasq']
tactsim_post_rho  <- gpqr_tactsim_post[,'rho']
true_beta0_with_tau09  <- qnorm(0.9,true.param$beta0,true.param$sigma) 


pdf(here('manuscript','supplementary_figures','figureS11.pdf'),height=8,width=8)
par(mfrow=c(2,2))
postHPDplot(tactsim_post_beta0,xlab='Cal BP',ylab='Posterior Probability',main=TeX('$\\beta_0$'),prob = 0.95)
abline(v=true_beta0_with_tau09,lty=2)
postHPDplot(tactsim_post_beta1,xlab='',ylab='Posterior Probability',main=TeX('$\\beta_1$'),prob=0.95)
abline(v=true.param$beta1,lty=2)
postHPDplot(tactsim_post_etasq,xlab='',ylab='Posterior Probability',main=TeX('$\\eta^2$'),prob=0.95)
abline(v=true.param$etasq,lty=2)
postHPDplot(tactsim_post_rho,xlab='km',ylab='Posterior Probability',main=TeX('$\\rho$'),prob=0.95)
abline(v=true.param$rho,lty=2)
dev.off()

# Figure S12 (Traceplot of beta0, beta1, rhosq, and etasq for tau=0.9) ----

load(here('results','gpqr_tau90.RData'))

pdf(file=here('manuscript','supplementary_figures','figureS12.pdf'),width=8,height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau90[,'beta0'],main=TeX('$\\beta_0$'),smooth=TRUE)
traceplot(gpqr_tau90[,'beta1'],main=TeX('$\\beta_1$'),smooth=TRUE)
traceplot(gpqr_tau90[,'rho'],main=TeX('$\\rho$'),smooth=TRUE)
traceplot(gpqr_tau90[,'etasq'],main=TeX('$\\eta^2$'),smooth=TRUE)
dev.off()

# Figure S13 (Traceplot of beta0, beta1, rhosq, and etasq for tau=0.99) ----

load(here('results','gpqr_tau99.RData'))

pdf(file=here('manuscript','supplementary_figures','figureS13.pdf'),width=8,height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau99[,'beta0'],main=TeX('$\\beta_0$'),smooth=TRUE)
traceplot(gpqr_tau99[,'beta1'],main=TeX('$\\beta_1$'),smooth=TRUE)
traceplot(gpqr_tau99[,'rho'],main=TeX('$\\rho$'),smooth=TRUE)
traceplot(gpqr_tau99[,'etasq'],main=TeX('$\\eta^2$'),smooth=TRUE)
dev.off()

# Figure S14 (Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.9) ----
gpqr.tau90.comb  <- do.call(rbind,gpqr_tau90)
pdf(file=here('manuscript','supplementary_figures','figureS14.pdf'),width=8,height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau90.comb[,'beta0'],main=TeX('$\\beta_0$'),xlab='Cal BP',ylab='')
postHPDplot(gpqr.tau90.comb[,'beta1'],main=TeX('$\\beta_1$'),xlab='',ylab='')
postHPDplot(gpqr.tau90.comb[,'rho'],main=TeX('$\\rho$'),xlab='km',ylab='')
postHPDplot(gpqr.tau90.comb[,'etasq'],main=TeX('$\\eta^2$'),xlab='',ylab='')
dev.off()

# Figure S15 (Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.99) ----
gpqr.tau99.comb  <- do.call(rbind,gpqr_tau99)
pdf(file=here('manuscript','supplementary_figures','figureS15.pdf'),width=8,height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau99.comb[,'beta0'],main=TeX('$\\beta_0$'),xlab='Cal BP',ylab='')
postHPDplot(gpqr.tau99.comb[,'beta1'],main=TeX('$\\beta_1$'),xlab='',ylab='')
postHPDplot(gpqr.tau99.comb[,'rho'],main=TeX('$\\rho$'),xlab='km',ylab='')
postHPDplot(gpqr.tau99.comb[,'etasq'],main=TeX('$\\eta^2$'),xlab='',ylab='')
dev.off()

# Figure S16 Tactical Simulation Posterior Predictive Check for nu and upsilon ----
load(here("results","phasemodel_tactsim.RData"))
post.model.a  <- do.call(rbind,mcmc.samples1)[,1:2]
post.model.b  <- do.call(rbind,mcmc.samples2)[,1:2]
dens.a.nu  <- density(post.model.a[,1],bw = 5)
dens.a.upsilon  <- density(post.model.a[,2],bw=5)
dens.b.nu  <- density(post.model.b[,1],bw=5)
dens.b.upsilon  <- density(post.model.b[,2],bw=5)

pdf(file=here('manuscript','supplementary_figures','figureS16.pdf'),width=8,height=8)
plot(NULL,xlim=c(3900,2500),ylim=c(0,0.022),xlab='Cal BP',ylab='Posterior Probability') 
polygon(c(dens.a.nu$x,rev(dens.a.nu$x)),c(rep(0,length(dens.a.nu$x)),rev(dens.a.nu$y)),border=NA,col=rgb(0,0.4,0,0.5))
polygon(c(dens.a.upsilon$x,rev(dens.a.upsilon$x)),c(rep(0,length(dens.a.upsilon$x)),rev(dens.a.upsilon$y)),border=NA,col=rgb(0,0.4,0,0.5))
polygon(c(dens.b.nu$x,rev(dens.b.nu$x)),c(rep(0,length(dens.b.nu$x)),rev(dens.b.nu$y)),border=NA,col=rgb(1,0.55,0,0.5))
polygon(c(dens.b.upsilon$x,rev(dens.b.upsilon$x)),c(rep(0,length(dens.b.upsilon$x)),rev(dens.b.upsilon$y)),border=NA,col=rgb(1,0.55,0,0.5))
abline(v=c(3500,3000),lty=2)
axis(3,at=c(3500,3000),labels=c(TeX('$\\nu$'),TeX('$\\upsilon$')))
legend('topright',legend=c('Non hierarchichal','Hierarchichal'),fill=c('darkgreen','darkorange'))
dev.off()

# Figure S17 Prior Predictive check for delta ----
nsim  <- 5000
set.seed(123)
gamma1  <- runif(nsim,1,20)
gamma2  <- rtruncnorm(nsim,mean=200,sd=100,1,500)
delta.mat = matrix(NA,ncol=1000,nrow=nsim)
for (i in 1:nsim) {delta.mat[i,] = dgamma(1:1000,gamma1[i],(gamma1[i]-1)/gamma2[i])}

pdf(file=here('manuscript','supplementary_figures','figureS17.pdf'),height=6,width=6)
plot(NULL,xlab=TeX('$\\delta$'),ylab='Probability Density',xlim=c(1,1000),ylim=c(0,0.02))
polygon(x=c(1:1000,1000:1),y=c(apply(delta.mat,2,quantile,prob=0.025),rev(apply(delta.mat,2,quantile,prob=0.975))),border=NA,col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:1000,1000:1),y=c(apply(delta.mat,2,quantile,prob=0.25),rev(apply(delta.mat,2,quantile,prob=0.75))),border=NA,col=rgb(0.25,0.41,0.88,0.5))
legend('topright',legend=c('50% percentile range','95% percentile range'),fill=c(rgb(0.67,0.84,0.9,0.5),rgb(0.25,0.41,0.88,0.5)))
dev.off()

# Figure S18 Marginal Posterior Distribution of nu, model a ----
load(here("results","phase_model_a.RData"))
out.comb.unif.model.a  <- do.call(rbind,out.unif.model_a)
post.nu.model.a  <- out.comb.unif.model.a[,paste0('a[',1:8,']')] |> round()

model.a.long  <- data.frame(value=as.numeric(post.nu.model.a),Area = rep(as.character(as.roman(1:8)),each=nrow(post.nu.model.a)))

pdf(file=here('manuscript','supplementary_figures','figureS18.pdf'),height=10,width=7)
ggplot(model.a.long, aes(x = value, y = Area,fill='lighblue')) + 
	geom_density_ridges() +
	scale_x_reverse(limits=c(3300,1800),breaks=BCADtoBP(c(-1200,-1000,-800,-600,-400,-200,1)),labels=c(1200,1000,800,600,400,200,1)) +
	scale_fill_manual(values='lightblue') +
	theme(legend.position = "none") +
	xlab('BC')
dev.off()

# Figure S19 Marginal Posterior Distribution of nu and upsilon, model b ----
load(here("results","phase_model_b.RData"))
out.comb.unif.model.b  <- do.call(rbind,out.unif.model_b)
post.nu.model.b  <- out.comb.unif.model.b[,paste0('a[',1:8,']')] |> round()
model.b.long  <- data.frame(value=as.numeric(post.nu.model.b),Area = rep(as.character(as.roman(1:8)),each=nrow(post.nu.model.b)))

pdf(file=here('manuscript','supplementary_figures','figureS19.pdf'),height=10,width=7)
ggplot(model.b.long, aes(x = value, y = Area,fill='lightblue')) + 
	geom_density_ridges() +
	scale_x_reverse(limits=c(3300,1800),breaks=BCADtoBP(c(-1200,-1000,-800,-600,-400,-200,1)),labels=c(1200,1000,800,600,400,200,1)) +
	scale_fill_manual(values='lightblue') +
	theme(legend.position = "none") +
	xlab('BC')
dev.off()

# Figure S20 Probability Matrix of nu, model a ----
pdf(file=here('manuscript','supplementary_figures','figureS20.pdf'),width=7,height=7.5)
orderPPlot(post.nu.model.a,name.vec=paste("Area",as.character(as.roman(1:8))))
dev.off()

# Figure S21 Difference Matrix plot of nu, model a ----
pdf(file=here('manuscript','supplementary_figures','figureS21.pdf'),width=16,height=11)
mat <- cbind(c(1,9:15),c(37,2,16:21),c(rep(37,2),3,22:26),c(rep(37,3),4,27:30),c(rep(37,4),5,31:33),c(rep(37,5),6,34:35),c(rep(37,6),7,36),c(rep(37,7),8))
layout(mat)
par(mar=c(0,0,0,0))
for (i in 1:8){plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F);text(0.5,0.5,paste('Area',as.character(as.roman(i))),cex=3)}
par(mar=c(3,0,0,1))
for (i in 1:8){
	for (j in 1:8){
		if (i < j)
		{
			diffDens(post.nu.model.a[,i],post.nu.model.a[,j],xlim=c(-1200,1200),prob=0.9)
		}
	}
}
dev.off()



