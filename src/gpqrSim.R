gpqrSim  <- function(win,n=150,seed=123,beta0=3000,beta1=0.7,sigma=100,etasq=0.05,rho=100,origin.point=c(129.959,33.4485))
{
	require(nimbleCarbon)
	require(rcarbon)
	require(sf)
	require(sp)
	require(maptools)
	set.seed(seed)
	sites <- spsample(win, n = n, type = 'random') 
	dist_mat  <- spDists(sites,longlat=TRUE)
	dist_org  <-  spDistsN1(sites,origin.point,longlat=TRUE)

	cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0),sigmasq = double(0)) 
				   {
					   returnType(double(2))
					   n <- dim(dists)[1]
					   result <- matrix(nrow = n, ncol = n, init = FALSE)
					   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
					   diag(deltaij) <- 1
					   for(i in 1:n)
						   for(j in 1:n)
							   result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
					   return(result)
				   })
	Ccov_ExpQ <- compileNimble(cov_ExpQ)
	assign('cov_ExpQ',cov_ExpQ,envir=.GlobalEnv)

	dispersalmodel <- nimbleCode({
		for (i in 1:N){
			# Model
			mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
			theta[i] ~ dnorm(mean=mu[i],sd=sigma)
		}
		mu_s[1:N] <- 0
		cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
		s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
	})

	#Define Parameters
	constants  <- list()
	constants$N <- n
	constants$dist_org <- dist_org
	constants$dist_mat <- dist_mat
	constants$beta0  <- beta0
	constants$beta1  <- beta1
	constants$sigma  <- sigma
	constants$etasq  <- etasq
	constants$rho  <- rho

	#Simulate
	set.seed(seed)
	simModel  <- nimbleModel(code=dispersalmodel,constants=constants)
	simModel$simulate('s')
	simModel$simulate('mu')
	simModel$simulate('theta')

	#Combine Results
	out.df  <- data.frame(ID=1:n,theta=simModel$theta)
	out.df$cra  <- round(uncalibrate(round(out.df$theta))$ccCRA)
	out.df$cra.error  <- 20
	out.df$med.date  <- medCal(calibrate(out.df$cra,out.df$cra.error,verbose=F))
	out.df$s  <- simModel$s
	out.df$rate  <- -1/(out.df$s-beta1)
	out.df$mu  <- beta0 + (out.df$s-beta1)*dist_org
	out  <- as(sites,'SpatialPointsDataFrame')
	out@data  <- out.df
	out.sf  <- as(out,'sf')

	#Store Output
	return(out.sf)
}
