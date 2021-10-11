# yayoi_rice_dispersal
Modelling Dispersal of Rice in Prehistoric Japan



## Model Parameters and Variables
### Variables
 - d_i ... distance from source
 - dmat ... distance matrix
 - cra_i ... observed 14C age i
 - cra.error_i ... observed 14C age error

### Parameters
 - beta0 ... intercept
 - beta1 ... average dispersal rate
 - s_i ... local deviation in the dispersal rate at site k
 - sigma ... regression error
 - tau ... regression quantile
 - alpha_k ... start date at site k
 - delta_k ... duration site k
 - omega ... duration gamma distribution shape 
 - phi ... duration gamma distribution rate
 - mu_omega ...hyperprior mean for omega
 - sigma_omega ... hyperprios sd for omega
 - mu_phi ...hyperprior mean for phi
 - sigma_phi ... hyperprios sd for phi
 - theta
 - eta ... covariance function parameter
 - rho ... covariance function parameter

