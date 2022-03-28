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
library(pals)
library(latex2exp)
# Load Data & Results ----
load(here('data','c14rice.RData'))
load(here('results','gpqr_tau90.RData'))
load(here('results','gpqr_tau99.RData'))
load(here('results','phase_model_a.RData'))
load(here('results','phase_model_b.RData'))

# Figure 1 (Sampling Site Distribution) ----

# Obtain Background  Map
win  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')

japan <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
df.pref.reg = data.frame(Pref = c("Aichi","Gifu","Akita","Ishikawa","Aomori","Chiba","Ehime","Saga","Fukuoka","Miyazaki","Kochi","Oita","Nagasaki","Kumamoto","Hiroshima","Fukui","Gunma","Hokkaido","Fukushima","Kanagawa","Hyogo","Ibaraki","Iwate","Kagoshima","Kagawa","Kyoto","Mie","Miyagi","Niigata","Nagano","Nara","Tokyo","Okayama","Okinawa","Osaka","Toyama","Tottori","Saitama","Shiga","Shimane","Shizuoka","Tochigi","Tokushima","Yamagata","Yamanashi","Wakayama","Yamaguchi"), Regions = c("Chubu","Chubu","Tohoku","Chubu","Tohoku","Kanto","Shikoku","Kyushu","Kyushu","Kyushu","Shikoku","Kyushu","Kyushu","Kyushu","Chugoku","Chubu","Kanto","Hokkaido","Tohoku","Kanto","Kansai","Kanto","Tohoku","Kyushu","Shikoku","Kansai","Kansai","Tohoku","Chubu","Chubu","Kansai","Kanto","Chugoku","Kyushu","Kansai","Chubu","Chugoku","Kanto","Kansai","Chugoku","Chubu","Kanto","Shikoku","Tohoku","Chubu","Kansai","Chugoku"))


japan@data <- left_join(japan@data,df.pref.reg,by=c('name'='Pref'))
japan <- gUnaryUnion(japan,id=japan@data$Regions)
japan.sf <- as(japan,'sf')
elevation <- get_elev_raster(locations = japan.sf, z = 4,clip = "locations",src='aws')
slope  <- terrain(elevation,opt='slope')
aspect  <- terrain(elevation,opt='aspect')
hill  <- hillShade(slope,aspect,40,270)
cropped_elev <- crop(elevation,japan.sf)
cropped_hill  <- crop(hill,japan.sf)
elevate <- as.data.frame(cropped_elev,xy = TRUE)
hs  <- as.data.frame(cropped_hill,xy = TRUE)
colnames(elevate)[3] = "elevation_value"
colnames(hs)[3] = "hs_value"
elevate <- elevate[complete.cases(elevate),]
hs <- hs[complete.cases(hs),]

# Convert sites into sf 
sites.sf  <- as(sites,'sf')

# Key Sites
key.sites  <- read.csv(here('data','key_sites_for_map.csv'))
coordinates(key.sites)  <- c('lon','lat')
proj4string(key.sites)  <- proj4string(sites)
key.sites.sf  <- as(key.sites,'sf')

# Distribution map
shade <- ggplot(hs, aes(x, y)) +
	geom_raster(aes(x=x,y=y,fill = hs_value), alpha = 0.5) +
	scale_fill_gradient2(low = "white", high = "white", mid = "black", midpoint = 0.6, guide = "none") +
	xlim(129,143) + 
	ylim(31,42)

grob.shade <- ggplotGrob(shade)
grob.shade <- grob.shade$grobs[[6]]$children[[3]]

f1 <-   ggplot() +
	geom_sf(data=win,aes(),fill='grey65',show.legend=FALSE,lwd=0) +
# 	geom_raster(data = elevate , aes(x = x, y = y,fill = elevation_value )) +
	annotation_custom(grob = grob.shade) +
	geom_sf(data=japan.sf,alpha=0,lwd=0.2,lty=1,col='white') +
	geom_sf(data=sites.sf,size=1.5,col='black',pch=21,fill='darkorange') + 
	geom_sf(data=key.sites.sf,size=2.5,col='black',fill='black',pch=24) +
 	scale_fill_gradientn(colours = kovesi.linear_gow_65_90_c35(100), limits=c(0, 3500)) +
	xlim(129,143) + 
	ylim(31,42) +
	annotate(geom='text',x=138.7, y=39.5, label="Tohoku",size=4) +
	annotate(geom='text',x=141.7, y=36, label="Kanto",size=4) +
	annotate(geom='text',x=135.5, y=37, label="Chubu",size=4) +
	annotate(geom='text',x=137.2, y=33.8, label="Kinki",size=4) +
	annotate(geom='text',x=134.2, y=32.8, label="Shikoku",size=4) +
	annotate(geom='text',x=133, y=31.8, label="Kyushu",size=4) +
	annotate(geom='text',x=131.5, y=36.2, label="Chugoku",size=4) +
	annotate("point",x=139,y=31.5,pch=21,size=1.5,col='black',fill='darkorange') +
	annotate("point",x=139,y=31,pch=24,size=2.5,col='black',fill='black')  +
	annotate("text", x=139.6, y=31.5, label='Sampled sites',size=3,hjust=0) +
	annotate("text", x=139.6, y=31, label='Key sites',size=3,hjust=0) +
	labs(fill = "Elevation (meters)",x='',y='') +
# 	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position=c(0.3,0.8),legend.text = element_text(size=7),legend.key.width= unit(0.2, 'in'),legend.key.size = unit(0.2, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=8),plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),axis.title=element_blank())
	theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.2),legend.position='none',plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),axis.title=element_blank())


f1 = f1 +
	annotate(geom='segment',x=c(137.99,140.13,129.99),y = c(36.63,40.81,33.71),xend = c(137.75,139.15,130.20),yend = c(37.47,41.12, 34.27)) +
	geom_text(data = data.frame(x = c(130.50,138.77, 130.27, 139.31, 132.44, 137.68),
				    y = c(34.02,41.23, 34.55, 35.00, 35.50, 37.72),
				    label = c("3","7,8", "1,2", "6", "4", "5")),mapping = aes(x = x, y = y, label = label),inherit.aes=FALSE)

yseq = seq(41.5,by=-0.35,length.out=nrow(key.sites.sf))
for (i in 1:nrow(key.sites.sf))
{
f1 = f1 + annotate(geom='text',x=130.5,y=yseq[i],label=paste0(i,'. ',key.sites.sf$name[i]),size=3,hjust=0)
}

pdf(file=here('manuscript','main_figures','figure1.pdf'),width=5,height=5)
f1
dev.off()


# Figure 2 (Posterior Mean of dispersal rate deviations) ----
post.gpqr.tau90  <- do.call(rbind,gpqr_tau90)
post.gpqr.tau99  <- do.call(rbind,gpqr_tau99)
nmcmc  <- nrow(post.gpqr.tau90) #number of MCMC samples (same for tau90 and tau99)
post.s.tau90  <- post.gpqr.tau90[,grep('s\\[',colnames(post.gpqr.tau90))]
post.s.tau99  <- post.gpqr.tau99[,grep('s\\[',colnames(post.gpqr.tau99))]

# post.arrival <- matrix(NA,nmcmc,constants$N.sites)
post.rate.tau90 <- post.rate.tau99  <-  matrix(NA,nmcmc,constants$N.sites)
for (i in 1:nmcmc)
{
	post.rate.tau90[i,] = -1 / (post.s.tau90[i,]-post.gpqr.tau90[i,'beta1'])
	post.rate.tau99[i,] = -1 / (post.s.tau99[i,]-post.gpqr.tau99[i,'beta1'])
}

sites@data$s.m.tau90 <- apply(post.s.tau90,2,median)
sites@data$s.m.tau99 <- apply(post.s.tau99,2,median)
sites@data$rate.m.tau90  <- apply(post.rate.tau90,2,median)
sites@data$rate.m.tau99  <- apply(post.rate.tau99,2,median)

sites.sf <- as(sites,'sf')

f2a <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=s.m.tau90),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title=TeX(r"(Posterior Median of s with $\tau = 0.90$)"),fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5,vjust = -1.5, size=10), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=6),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=7),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


f2b <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=s.m.tau99),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title=TeX(r"(Posterior Median of s with $\tau = 0.99$)"),fill='s') + 
	scale_fill_gradient2(low='blue',high='red',mid='white') +
	theme(plot.title = element_text(hjust = 0.5,vjust = -1.5, size=10), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=6),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=7),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))


f2c <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=rate.m.tau90),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.90$)"),fill='Dispersal Rate (km/year)') + 
	scale_fill_viridis(option="turbo",limits=c(0,4)) +
	theme(plot.title = element_text(hjust = 0.5,vjust = -1.5, size=10), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=6),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=7),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

f2d <- ggplot() +
	geom_sf(data=win,aes(),fill='grey66',show.legend=FALSE,lwd=0) +
	geom_sf(data=sites.sf,mapping = aes(fill=rate.m.tau99),pch=21,col='darkgrey',size=2) + 
	xlim(129,143) + 
	ylim(31,42) +
	labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.99$)"),fill='Dispersal Rate (km/year)') + 
	scale_fill_viridis(option="turbo",limits=c(0,4)) +
	theme(plot.title = element_text(hjust = 0.5,vjust=-1.5,size=10), panel.background = element_rect(fill='lightblue'),panel.grid.major = element_line(size = 0.1),legend.position=c(0.2,0.8),legend.text = element_text(size=6),legend.key.width= unit(0.1, 'in'),legend.key.size = unit(0.08, "in"),legend.background=element_rect(fill = alpha("white", 0.5)),legend.title=element_text(size=7),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin = unit(c(0,0,0,0), "in"))

pdf(file=here('manuscript','main_figures','figure2.pdf'),width=7,height=7)
grid.arrange(f2a,f2b,f2c,f2d,ncol=2,padding=0)
dev.off()

#Figure 3 Estimated Arrival Date ----
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

pdf(file=here('manuscript','main_figures','figure3.pdf'),width=7,height=3.9,pointsize=4)
par(mfrow=c(1,2),mar=c(3,3,3,1))
# Map
plot(win,xlim=c(127,143),ylim=c(31,43),col='lightgrey')
plot(japan.sf,col=main.col,border=NA,add=TRUE)
axis(1,at=seq(129,143,2),tck=-0.01,padj=-0.8)
axis(2,at=seq(30,42,2),tck=-0.01,padj=0.8)
mtext(side=1,line =1.4,'Longitude')
mtext(side=2,line=1.4,'Latitude')
box()
rect(xleft=127.5,xright=136,ybottom=38.2,ytop=43,col='white')
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
##Model a
text(x=129,y=42.5,labels='Model a',cex=1)
text(x=seq(128,by=1,length.out=8),y=rep(41.9,8),labels=c('I','II','III','IV','V','VI','VII','VIII'),cex=1)
##Model b
text(x=129,y=40.6,labels='Model b',cex=1)
text(x=seq(128,by=1,length.out=8),y=rep(40,8),labels=c('I','II','III','IV','V','VI','VII','VIII'),cex=1)
pos3 = cbind(seq(128,by=1,length.out=8),rep(40,8))

straightarrow(from=pos3[1,],to=pos3[2,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[1, ], to = pos3[3, ]-c(0,0.25),curve = 0.4, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1) 
straightarrow(from=pos3[3,],to=pos3[4,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
straightarrow(from=pos3[4,],to=pos3[5,],arr.type='triangle',segment=c(0.3,0.8),endhead=TRUE,arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[6, ]-c(0,0.25),curve = 0.4, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[7, ]-c(0,0.25),curve = 0.42, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)
curvedarrow(from = pos3[4, ], to = pos3[8, ]+c(0,0.25),curve = -0.3, arr.pos = 1,endhead = TRUE,segment=c(0.1,0.5),arr.width=aw,arr.length=al,lwd=1)

# Posterior Arrival Times
plot(NULL,xlim=c(3600,1900),ylim=c(0.5,23.5),xlab='',ylab='',axes=F)
tmp.a = extract(out.unif.model_a)
tmp.b = extract(out.unif.model_b)
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	post.bar(tmp.a[,i],i=iseq.a[i],h=0.9,col=main.col[i])
	post.bar(tmp.b[,i],i=iseq.b[i],h=0.9,col=main.col[i])
}


axis(2,at=iseq.a+0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1,at=BCADtoBP(c(-1500,-1250,-1000,-750,-500,-250,1)),labels=c('1500BC','1250BC','1000BC','750BC','500BC','250BC','1AD'),tck=-0.01)
axis(3,at=seq(3400,1800,-300),labels=paste0(seq(3400,1800,-300),'BP'),tck=-0.01)
box()

post.bar(c(2600,2500,2400,2300,2200,2100,2000),i=1.5,h=0.9,col='lightgrey')
arrows(x0=2500,x1=2100,y0=0.3,y1=0.3,angle = 90,code = 3,length = 0.01)
arrows(x0=2400,x1=2200,y0=0.8,y1=0.8,angle = 90,code = 3,length = 0.01)
text(x=2080,y=0.8,"50% HPDI",cex=0.8)
text(x=1980,y=0.3,"95% HPDI",cex=0.8)
text(x=2070,y=2.3,"Median Posterior",cex=0.8)
lines(x=c(2300,2230),y=c(2.1,2.3))
text(x=3490,y=2,'Model a',cex=1.1)
text(x=3490,y=1,'Model b',cex=1.1)

dev.off()

