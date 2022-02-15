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
library(dplyr)
library(rgeos)
library(elevatr)
library(raster)
library(diagram)
library(RColorBrewer)
# Load Data & Results ----
load(here('data','c14rice.RData'))
load(here('results','gpqr_tau90.RData'))
# load(here('results','gpqr_tau99.RData'))
load(here('results','phase_model0.RData'))
load(here('results','phase_model1.RData'))
load(here('results','phase_model2.RData'))

# Figure 1 (Sampling Site Distribution) ----

# Obtain Background  Map
win  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')

japan <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
df.pref.reg = data.frame(Pref = c("Aichi","Gifu","Akita","Ishikawa","Aomori","Chiba","Ehime","Saga","Fukuoka","Miyazaki","Kochi","Oita","Nagasaki","Kumamoto","Hiroshima","Fukui","Gunma","Hokkaido","Fukushima","Kanagawa","Hyogo","Ibaraki","Iwate","Kagoshima","Kagawa","Kyoto","Mie","Miyagi","Niigata","Nagano","Nara","Tokyo","Okayama","Okinawa","Osaka","Toyama","Tottori","Saitama","Shiga","Shimane","Shizuoka","Tochigi","Tokushima","Yamagata","Yamanashi","Wakayama","Yamaguchi"), Regions = c("Chubu","Chubu","Tohoku","Chubu","Tohoku","Kanto","Shikoku","Kyushu","Kyushu","Kyushu","Shikoku","Kyushu","Kyushu","Kyushu","Chugoku","Chubu","Kanto","Hokkaido","Tohoku","Kanto","Kansai","Kanto","Tohoku","Kyushu","Shikoku","Kansai","Kansai","Tohoku","Chubu","Chubu","Kansai","Kanto","Chugoku","Kyushu","Kansai","Chubu","Chugoku","Kanto","Kansai","Chugoku","Chubu","Kanto","Shikoku","Tohoku","Chubu","Kansai","Chugoku"))


japan@data <- left_join(japan@data,df.pref.reg,by=c('name'='Pref'))
japan <- gUnaryUnion(japan,id=japan@data$Regions)
japan.sf <- as(japan,'sf')
elevation <- get_elev_raster(locations = japan.sf, z = 4,clip = "locations",src='aws')
cropped_elev <- crop(elevation,japan.sf)
elevate <- as.data.frame(cropped_elev,xy = TRUE)
colnames(elevate)[3] = "elevation_value"
elevate <- elevate[complete.cases(elevate),]


# Convert sites into sf 
sites.sf  <- as(sites,'sf')

# Distribution map
myPalette <- colorRampPalette((brewer.pal(9, "YlOrBr")))
f1 <-   ggplot() +
	geom_sf(data=win,aes(),fill='grey37',show.legend=FALSE,lwd=0) +
	geom_raster(data = elevate , aes(x = x, y = y,fill = elevation_value )) +
	geom_sf(data=japan.sf,alpha=0,lwd=0.8,col='grey68') +
	geom_sf(data=sites.sf,size=1.5) + 
	scale_fill_gradientn(colours = myPalette(100), limits=c(0, 3500)) +
	xlim(129,143) + 
	ylim(31,42) +
	annotate(geom='text',x=138.7, y=39.5, label="Tohoku",size=4) +
	annotate(geom='text',x=141.7, y=36, label="Kanto",size=4) +
	annotate(geom='text',x=136.7, y=38, label="Chubu",size=4) +
	annotate(geom='text',x=136.3, y=33, label="Kinki",size=4) +
	annotate(geom='text',x=134.2, y=32.8, label="Shikoku",size=4) +
	annotate(geom='text',x=133, y=31.8, label="Kyushu",size=4) +
	annotate(geom='text',x=131.1, y=35.6, label="Chugoku",size=4) +
	labs(fill = "Elevation (meters)",x='',y='') +
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.3,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.2, 'in'),legend.key.size = unit(0.2, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),axis.title=element_blank())

pdf(file=here('manuscript','main_figures','figure1.pdf'),width=5,height=5)
f1
dev.off()


# Figure 2 (SPD of charred rice dates) -----
# SPD/cKDE
caldates.rice  <- calibrate(DateInfo$cra,DateInfo$cra_error,normalised=FALSE)
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

f2  <- ggplot(data=df.combo) + 
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
	labs(y = "Normalised Summed Probability") +
	theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.22,0.8),legend.spacing.y = unit(-11, "pt"),legend.background = element_rect(fill = "transparent"), legend.text=element_text(size=6),axis.title=element_text(size=9),plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"))

pdf(file=here('manuscript','main_figures','figure2.pdf'),width = 5, height = 4)
f2
dev.off()



# Figure 3 (Posterior Mean of dispersal rate deviations) ----
post.gpqr  <- do.call(rbind,gpqr_tau90)
nmcmc  <- nrow(post.gpqr) #number of MCMC samples
post.s  <- post.gpqr[,grep('s\\[',colnames(post.gpqr))]
# post.arrival <- matrix(NA,nmcmc,constants$N.sites)
post.rate <- matrix(NA,nmcmc,constants$N.sites)
for (i in 1:nmcmc)
{
	post.rate[i,] = -1 / (post.s[i,]-post.gpqr[i,'beta1'])
# 	post.arrival[i,]  <- BPtoBCAD(post.gpqr[i,'beta0'] + (post.gpqr[i,'beta1'] + post.s[i,]) * constants$dist_org)
}

sites@data$s.m <- apply(post.s,2,median)
# sites@data$arrival.m  <- apply(post.arrival,2,median)
sites@data$rate.m  <- apply(post.rate,2,median)

sites.sf <- as(sites,'sf')

f3a <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=s.m),pch=21,col='darkgrey',size=1.2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='a',fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5,size=6), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=4),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=4.5),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


f3b <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=rate.m),pch=21,col='darkgrey',size=1.2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title='b',fill='Dispersal Rate \n (km/yr)') + 
	scale_fill_viridis(option="turbo") +
	theme(plot.title = element_text(hjust = 0.5,size=6), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=4),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=4.5),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

pdf(file=here('manuscript','main_figures','figure3.pdf'),width=2.9,height=4.7)
grid.arrange(f3a,f3b,ncol=2,padding=0)
dev.off()

#Figure 4 Estimated Arrival Date ----
# Setup Functions and Variables 
extract <- function(x)
{
	tmp = do.call(rbind,x)
	tmp2 = tmp[,grep('^a\\[',colnames(tmp))]
	qta = apply(tmp2,2,quantile,prob=c(0,0.025,0.25,0.5,0.75,0.975,1))
	return(qta)
}

post.bar <- function(x,i,h,col)
{
# 	lines(c(x[1],x[7]),c(i,i),col=col)
	rect(xleft=x[2],xright=x[6],ybottom=i-h/5,ytop=i+h/5,border=NA,col=col)
	rect(xleft=x[3],xright=x[5],ybottom=i-h/3,ytop=i+h/3,border=NA,col=col)
	lines(c(x[4],x[4]),c(i-h/2,i+h/2),lwd=2,col='grey44')
}


main.col <- c(rgb(51,34,136,maxColorValue=255),rgb(136,204,238,maxColorValue = 255),rgb(68,170,153,maxColorValue = 255),rgb(17,119,51,maxColorValue = 255),rgb(153,153,51,maxColorValue = 255),rgb(221,204,119,maxColorValue = 255),rgb(204,102,119,maxColorValue = 255),rgb(136,34,85,maxColorValue = 255))

win  <- ne_countries(continent = 'asia',scale=10)
japan <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
df.pref.reg <- read.csv(here('data','prefecture_region_match.csv'))
japan@data <- left_join(japan@data,df.pref.reg,by=c('name'='Prefecture'))
japan <- gUnaryUnion(japan,id=japan@data$Area)
japan.sf <- as(japan,'sf')
win <- gUnaryUnion(win,id=win@data[,1])

pdf(file=here('manuscript','main_figures','figure4.pdf'),width=7,height=3.9,pointsize=4)
par(mfrow=c(1,2),mar=c(3,3,3,1))
# Map
plot(win,xlim=c(127,143),ylim=c(31,43),col='lightgrey')
legend('bottomright',legend='a',bty = 'n',cex=2)
plot(japan.sf,col=main.col,border=NA,add=TRUE)
axis(1,at=seq(129,143,2),tck=-0.01,padj=-0.8)
axis(2,at=seq(30,42,2),tck=-0.01,padj=0.8)
mtext(side=1,line =1.4,'Longitude')
mtext(side=2,line=1.4,'Latitude')
box()
rect(xleft=127.5,xright=136,ybottom=37,ytop=43,col='white')
text(130,34,labels='I',cex=1.5)
text(132,32,labels='II',cex=1.5)
text(134,32.5,labels='III',cex=1.5)
text(137,33.2,labels='IV',cex=1.5)
text(138.5,33.8,labels='V',cex=1.5)
text(141.5,36,labels='VI',cex=1.5)
text(139,39,labels='VII',cex=1.5)
text(142.3,41,labels='VIII',cex=1.5)

# Model Diagrams
aw = 0.05
al = 0.1
##Model 0
text(x=129,y=42.5,labels='Model 0',cex=1)
text(x=seq(128,by=1,length.out=8),y=rep(41.9,8),labels=c('I','II','III','IV','V','VI','VII','VIII'),cex=1)
##Model 1
text(x=129,y=41,labels='Model 1',cex=1)
text(x=seq(128,by=1,length.out=8),y=rep(40.4,8),labels=c('I','II','III','IV','V','VI','VII','VIII'),cex=1)
pos2 = cbind(seq(128,by=1,length.out=8),rep(40.4,8))
straightarrow(from=pos2[1,],to=pos2[2,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd = 1)
curvedarrow(from = pos2[1, ], to = pos2[3, ]-c(0,0.25),curve = 0.3, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)
for (i in 3:7)
{
  straightarrow(from = pos2[i, ], to = pos2[i+1, ],arr.type='triangle',segment=c(0.3,0.8),endhead=T,arr.width=aw,arr.length=al,lwd=1)
}
##Model 2
text(x=129,y=39.2,labels='Model 2',cex=1)
text(x=seq(128,by=1,length.out=8),y=rep(38.6,8),labels=c('I','II','III','IV','V','VI','VII','VIII'),cex=1)
pos3 = cbind(seq(128,by=1,length.out=8),rep(38.6,8))
straightarrow(from=pos3[1,],to=pos3[2,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[1, ], to = pos3[3, ]-c(0,0.25),curve = 0.4, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1) 
straightarrow(from=pos3[3,],to=pos3[4,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
straightarrow(from=pos3[4,],to=pos3[5,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[6, ]-c(0,0.25),curve = 0.4, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[7, ]-c(0,0.25),curve = 0.42, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[8, ]+c(0,0.25),curve = -0.3, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)

# Posterior Arrival Times
plot(NULL,xlim=c(3600,1900),ylim=c(0.5,31.5),xlab='',ylab='',axes=F)
legend('bottomright',legend='b',bty = 'n',cex=2)
tmp0 = extract(out.unif.model0)
tmp1 = extract(out.unif.model1)
tmp2 = extract(out.unif.model2)

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
axis(1,at=BCADtoBP(c(-1500,-1250,-1000,-750,-500,-250,1)),labels=c('1500BC','1250BC','1000BC','750BC','500BC','250BC','1AD'),tck=-0.01)
axis(3,at=seq(3400,1800,-300),labels=paste0(seq(3400,1800,-300),'BP'),tck=-0.01)
box()
text(x=BCADtoBP(c(-400,-400,-400)),y=c(iseq0[1],iseq1[1],iseq2[1]),labels=c('Model0','Model1','Model2'))
dev.off()

