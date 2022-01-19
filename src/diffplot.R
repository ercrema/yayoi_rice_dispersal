diffDens  <- function(x,y,prob=0.9,...)
{
	require(coda)
	z  <- x - y 
	nsample = length(z)
	left  <- c(HPDinterval(mcmc(z),prob=prob)[1],0)
	right  <- c(0,HPDinterval(mcmc(z),prob=prob)[2])
	plotRight=plotLeft=TRUE
	if (any(right<0)){left[2]=right[2];plotRight=FALSE}
	if (any(left>0)){right[1]=left[1];plotLeft=FALSE}
	dens = density(z)
	hpdi.left.x = dens$x[which(dens$x>=left[1]&dens$x<=left[2])]
	hpdi.left.y = dens$y[which(dens$x>=left[1]&dens$x<=left[2])]
	hpdi.right.x = dens$x[which(dens$x>=right[1]&dens$x<=right[2])]
	hpdi.right.y = dens$y[which(dens$x>=right[1]&dens$x<=right[2])]
	plot(dens$x,dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,...)
	if(plotLeft){polygon(x=c(hpdi.left.x,rev(hpdi.left.x)),y=c(hpdi.left.y,rep(0,length(hpdi.left.y))),border=NA,col='lightpink')}
	if(plotRight){polygon(x=c(hpdi.right.x,rev(hpdi.right.x)),y=c(hpdi.right.y,rep(0,length(hpdi.right.y))),border=NA,col='lightblue')}
	lines(dens)
	abline(v=0,lty=2,lwd=1.5)
	axis(1,cex.axis=0.85,padj=-0.5)
	xlim = range(axTicks(1))
	axis(1,at=seq(xlim[1],xlim[2],100),labels=NA,tck=-.01)
# 	axis(2)
# 	box()
}
