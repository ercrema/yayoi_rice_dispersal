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
library(diagram)
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

# Figure S3 (Posterior vs True values of s and alpha) ----

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

# Figure S4 (Posterior vs True values of beta0,beta1,rhosq,etasq,omega, and phi) ----
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

# Figure S13 (Posterior Arrival Date Trapezoidal Model ----
load(here('results','trap_model0.RData'))
load(here('results','trap_model1.RData'))
load(here('results','trap_model2.RData'))

extract <- function(x)
{
	tmp = do.call(rbind,x)
	tmp2 = tmp[,grep('^a\\[',colnames(tmp))]
	qta = apply(tmp2,2,quantile,prob=c(0,0.025,0.25,0.5,0.75,0.975,1))
	return(qta)
}

post.bar <- function(x,i,h,col)
{
	lines(c(x[1],x[7]),c(i,i),col=col)
	rect(xleft=x[2],xright=x[6],ybottom=i-h/5,ytop=i+h/5,border=NA,col=col)
	rect(xleft=x[3],xright=x[5],ybottom=i-h/3,ytop=i+h/3,border=NA,col=col)
	lines(c(x[4],x[4]),c(i-h/2,i+h/2),lwd=2,col='grey44')
}


main.col <- c(rgb(51,34,136,maxColorValue=255),rgb(136,204,238,maxColorValue = 255),rgb(68,170,153,maxColorValue = 255),rgb(17,119,51,maxColorValue = 255),rgb(153,153,51,maxColorValue = 255),rgb(221,204,119,maxColorValue = 255),rgb(204,102,119,maxColorValue = 255),rgb(136,34,85,maxColorValue = 255))

pdf(file=here('figures','figureS13.pdf'),width=5.5,height=7,pointsize=9)
par(mar=c(3,3,3,1))
# Map
# Posterior Arrival Times
plot(NULL,xlim=c(4000,2000),ylim=c(0.5,31.5),xlab='',ylab='',axes=F)
tmp0 = extract(out.trap.model0)
tmp1 = extract(out.trap.model1)
tmp2 = extract(out.trap.model2)

iseq0 = seq(1,by=4,length.out=8)
iseq1 = seq(2,by=4,length.out=8)
iseq2 = seq(3,by=4,length.out=8)
abline(h=seq(4,by=4,length.out=7),col='darkgrey',lty=2)
for (i in 1:8)
{
	post.bar(tmp0[,i],i=iseq0[i],h=0.7,col=main.col[i])
	post.bar(tmp1[,i],i=iseq1[i],h=0.7,col=main.col[i])
	post.bar(tmp2[,i],i=iseq2[i],h=0.7,col=main.col[i])
}
axis(2,at=iseq1,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1,at=BCADtoBP(c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1)),labels=c('2000BC','1740BC','1500BC','1250BC','1000BC','750BC','500BC','250BC','1AD'),tck=-0.01)
axis(3,at=seq(3900,1800,-300),labels=paste0(seq(3900,1800,-300),'BP'),tck=-0.01)
box()
text(x=c(2600,2600,2600),y=c(iseq0[1],iseq1[1],iseq2[1]),labels=c('Model0','Model1','Model2'))
dev.off()

