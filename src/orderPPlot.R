orderPPlot  <- function(x,name.vec)
{
	order.prob <- matrix(NA,nrow=ncol(x),ncol=ncol(x),dimnames=list(name.vec,name.vec))
	for (i in 1:ncol(x)){
		for (j in 1:ncol(x)){
			if (i>j)
			{
				order.prob[i,j] = round(sum(x[,i]<x[,j])/nrow(x),2)
			}
		}
	}
	corrplot(order.prob,is.corr=FALSE,diag=FALSE,col=COL2('RdBu'),addCoef.col = 'white',method='color',type='lower',tl.col="black")
}
