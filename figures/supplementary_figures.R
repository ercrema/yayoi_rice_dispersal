# Load Libraries and spatial data ----
library(here)
library(truncnorm)
library(ggplot2)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(maptools)
library(viridis)
library(latex2exp)
library(gridExtra)
library(quantreg)

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

# Figure S2 (Simulated dispersal rate, local deviations, and arrival dates) ----
# Load Spatial Data 
win.sf  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')
# Load Tactical Simulation
load(here('results','tactical_sim_res.RData'))

fS2a  <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=s),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='a',fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

fS2b <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=rate),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='b',fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="magma") +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

fS2c <- ggplot() +
	geom_sf(data=win.sf,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sim.sites,mapping = aes(fill=arrival),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='c',fill='1st Perc.') + 
	scale_fill_continuous(type = "viridis", breaks = c(-1000,-600,-200,200), labels = c('1000BC','600BC','200BC','200AD'),limits=c(-1000,250)) +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

pdf(file=here('figures','figureS2.pdf'),width=10,height=4)
grid.arrange(fS2a,fS2b,fS2c,ncol=3)
dev.off()

# Figure S3 (Posterior vs True values of s and alpha) ----

# Load GPQR on tactical simulation
load(here('results','res_tactical.RData'))

tactical_fitted  <- res[[1]] #Extract chain #1
fitted_s  <- tactical_fitted[,grep('s\\[',colnames(tactical_fitted))]
fitted_alpha  <- tactical_fitted[,grep('alpha\\[',colnames(tactical_fitted))]

pdf(file=here('figures','figureS3.pdf'),width=8,height=4.5)
par(mfrow=c(1,2))
plot(sim.SiteInfo$s,apply(fitted_s,2,mean),pch=20,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5),xlab=TeX('Simulated $s_k$'),ylab=TeX('Predicted $s_k$'),main='a')
lo_s <- apply(fitted_s,2,quantile,0.025)
hi_s  <- apply(fitted_s,2,quantile,0.975)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.SiteInfo$s[i],2),c(lo_s[i],hi_s[i]),lwd=0.5)
}
abline(a=0,b=1,lty=2,col=2,lwd=2)

plot(sim.SiteInfo$alpha,apply(fitted_alpha,2,mean),pch=20,xlim=c(1600,3200),ylim=c(1600,3200),xlab=TeX('Simulated $\\alpha_k$'),ylab=TeX('Predicted $\\alpha_k$'),main='b')
lo_alpha <- apply(fitted_alpha,2,quantile,0.025)
hi_alpha  <- apply(fitted_alpha,2,quantile,0.975)

for (i in 1:nrow(sim.SiteInfo))
{
lines(rep(sim.SiteInfo$alpha[i],2),c(lo_alpha[i],hi_alpha[i]),lwd=0.5)
}
abline(a=0,b=1,lty=2,col=2,lwd=2)

dev.off()

# Figure S4 (Posterior vs True values of beta0,beta1,rhosq,etasq,omega, and phi) ----

fitted_beta0  <- tactical_fitted[,'beta0']
fitted_beta1  <- tactical_fitted[,'beta1']
fitted_etasq  <- tactical_fitted[,'etasq']
fitted_rhosq  <- tactical_fitted[,'rhosq']
fitted_omega  <- tactical_fitted[,'omega']
fitted_phi  <- tactical_fitted[,'phi']

pdf(file=here('figures','figureS4.pdf'),width=6,height=9)

par(mfrow=c(3,2))
postHPDplot(fitted_beta0,main=TeX('Posterior $\\beta_0$'),HPD = 0.95)
abline(v=sim.constants$beta0,lty=2)

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

# Figure S6 (Prior Predictive Check etasq and rhosq) ----
nsim  <- 500
etasq.prior  <- rexp(nsim,10)
rhosq.prior  <- rexp(nsim,1000)

pdf(file=here('figures','figureS6.pdf'),width=6,height=5)
plot(NULL,xlab='Distance (km)',ylab='Covariance',xlim=c(0,500),ylim=c(0,0.6))
for (i in 1:nsim)
{
	curve(etasq.prior[i]*exp(-rhosq.prior[i]*x^2),add=TRUE,from=0,to=500,col=rgb(0,0,0,0.1))
}
dev.off()

# Figure S7 (Prior Predictive Check of omega and phi) ----
nsim  <- 500
omega.prior  <- rtruncnorm(nsim,a=0,mean=3,sd=1)
phi.prior  <- rtruncnorm(nsim,a=0,mean=0.01,sd=0.01)


pdf(file=here('figures','figureS7.pdf'),width=6,height=5)
plot(NULL,xlab='Duration',ylab='Probability Density',xlim=c(0,1500),ylim=c(0,0.02))
for (i in 1:nsim)
{
	curve(dgamma(x,omega.prior[i],phi.prior[i]),add=TRUE,from=0,to=1500,col=rgb(0,0,0,0.1))
}
dev.off()

# Figure S8 (Traceplot of beta0, beta1, rhosq, etasq, omega, and phi) ----

# Figure S9 (Marginal posterior of beta0, beta1, rhosq, etasq, omega, and phi) ----

# Figure S10 (Joint Posterior beta0 & beta1) ----

# Figure S11 (Joint Posterior rhosq and etasq) ----

# Figure S12 (Joint Posterior omega and phi) ----
