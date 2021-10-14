# Load R libraries ----
library(rnaturalearth)
library(maptools)
library(rgeos)
library(ggplot2)
library(here)
library(viridis)
library(rcarbon)
library(gridExtra)
library(ggspatial)

# Load Data & Results ----
load(here('data','c14rice.RData'))
load(here('results','gpqr_res.RData'))

# Figure 1 (Site Distribution and SPD of charred rice dates) ----

# Obtain Background  Map
win  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')

# Convert sites into sf 
sites.sf  <- as(sites,'sf')

# Distribution map
f1a  <- ggplot() + 
	geom_sf(data=win,aes(),fill='cornsilk3',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,size=.8) + 
	xlim(129,143) + 
	ylim(31,42) +
	annotation_scale(location='br') +
	labs(title='a') + 
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2))
# SPD/cKDE
caldates.rice  <- calibrate(d$cra,d$cra_error,normalised=FALSE)
spd.rice  <- spd(caldates.rice,timeRange=c(3000,500),spdnormalised=TRUE)
boot.rice.dates  <- sampleDates(caldates.rice,boot=TRUE,nsim=1000)
boot.ckde  <- ckde(boot.rice.dates,timeRange=c(3000,500),bw=50)

df.combo  <-  data.frame(BCAD = BPtoBCAD(spd.rice$grid$calBP),
		    SPD = spd.rice$grid$PrDens,
		    mCKDE = apply(boot.ckde$res.matrix,1,quantile,0.5,na.rm=TRUE), 
		    loBoot = apply(boot.ckde$res.matrix,1,quantile,0.025,na.rm=TRUE),
		    hiBoot = apply(boot.ckde$res.matrix,1,quantile,0.975,na.rm=TRUE))
df.spd <- reshape2::melt(df.combo, id.var = "BCAD") |> subset(variable%in%c('SPD','mCKDE'))
df.boot <- reshape2::melt(df.combo, id.var = "BCAD") |> subset(variable%in%c('hiBoot','loBoot'))

f1b  <- ggplot(data=df.combo) + 
	geom_line(aes(x=BCAD,y=mCKDE,color='cKDE (Median)',linetype='cKDE (Median)',size='cKDE (Median)')) +
	geom_ribbon(aes(ymin=loBoot, ymax=hiBoot, x=BCAD, fill = "cKDE Bootstrap 95% CI"), alpha = 0.3) +
	geom_line(aes(x=BCAD,y=SPD,color='SPD',linetype='SPD',size='SPD')) +
	geom_point(data=DateInfo,aes(x=BPtoBCAD(median.dates),shape='Median Calibrated Dates'),y=0.0015,alpha=0.5) + 
	scale_colour_manual("",values=c("black","black")) +
	scale_linetype_manual("",values=c('solid','dashed')) +
	scale_size_manual("",values=c(0.5,0.2)) +
	scale_fill_manual("",values="cadetblue") +
	scale_shape_manual("",values=3) +
	scale_x_continuous('', breaks=c(-1000,-600,-200,200,600,1000),labels=c('1000BC','600BC','200BC','200AD','600AD','1000AD'),limits=c(-1100,1000)) + 
	labs(title="b", y = "Normalised Summed Probability") +
	theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.22,0.8),legend.spacing.y = unit(-11, "pt"),legend.background = element_rect(fill = "transparent"), legend.text=element_text(size=6),axis.title=element_text(size=9))

pdf(file=here('figures','figure1.pdf'),width = 7.2, height = 3.7,pointsize=10)
grid.arrange(f1a,f1b,ncol=2)
dev.off()

# Figure 3 (Posterior Mean of dispersal rate deviations and arrival time) ----
post.gpqr <- gpqr_uniform_sample[[1]] #Extract posterior from chain #1
nmcmc  <- nrow(post.gpqr) #number of MCMC samples
post.s  <- post.gpqr[,grep('s\\[',colnames(post.gpqr))]
post.arrival <- matrix(NA,nmcmc,constants$N.sites)
post.rate <- matrix(NA,nmcmc,constants$N.sites)
for (i in 1:nmcmc)
{
	post.rate[i,] = - 1 / (post.gpqr[i,'betaD'] + post.s[i,])
	post.arrival[i,]  <- BPtoBCAD(post.gpqr[i,'intercept'] + (post.gpqr[i,'betaD'] + post.s[i,]) * constants$dist_org)
}

sites@data$s.m <- apply(post.s,2,median)
sites@data$arrival.m  <- apply(post.arrival,2,median)
sites@data$rate.m  <- apply(post.rate,2,median)

sites.sf <- as(sites,'sf')

f2a <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=s.m),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='a',fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


f2b <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=rate.m),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='b',fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="magma") +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))



f2c <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=arrival.m),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='c',fill='1th Perc.') + 
	scale_fill_continuous(type = "viridis", breaks = c(-1000,-600,-200,200), labels = c('1000BC','600BC','200BC','200AD'),limits=c(-1000,250)) +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.2,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.1, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


pdf(file=here('figures','figure2.pdf'),width=2.9,height=7)
grid.arrange(f2a,f2b,f2c,nrow=3,padding=0)
dev.off()
