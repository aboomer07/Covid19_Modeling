# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/26/2021

library(zoo)

gen_distribution <- function(study_len, mean, variance, type, delta) {

  r <- 0:(study_len*delta)
  r <- r/delta

  if (type == 'norm') {
    a <- mean
    b <- variance
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dnorm(r, a, sqrt(b))
  }

  if (type == 'lnorm') {
    a <- log(mean^2 / (sqrt(mean^2 + variance)))
    b <- log(1 + variance / mean^2)
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dlnorm(r,a,sqrt(b))
  }

  if (type == 'gamma') {
    a <- (mean^2)/variance
    b <- mean/variance
    # b <- 1/b
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dgamma(r, a, b)
  }

  if (type == 'weibull') {
    a <- as.numeric(weibullpar(mean, sqrt(variance))[1])
    b <- as.numeric(weibullpar(mean, sqrt(variance))[2])
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dweibull(r, a, b)
  }

  distribution <- list(a = a, b = b, omega = omega, range = r)
  return(distribution)
}


#-----------------------------------------------------------------#
#--------------------------- Some checks -------------------------#
#-----------------------------------------------------------------#
delta = 24
omega <- gen_distribution(20, 7, 2, "weibull", delta)
plot(omega$omega~omega$range)
#Check if integral is 1 
area <- c()
for (i in 1:(length(omega$range)-1)){
  area[i] <- (omega$range[i+1]-omega$range[i])*omega$omega[i]
}
sum(area)  #should be 1

#Check if mean and variance are the same as the inputted ones
prod <- c()
for (i in 1:length(omega$range)){
  prod[i] <- omega$range[i]*omega$omega[i]
}
sum(prod)/delta #should be 7

squares <- c()
for (i in 1:length(omega$range)){
  squares[i] <- omega$range[i]^2*omega$omega[i]
}
(sum(squares)/delta)-(sum(prod)/delta)^2 #should be 2


###############################################################
########### Simulate Serial Interval Data #####################
###############################################################

samp_pois <- function(params, increasingR = TRUE) {

  R_start <- params[['R_start']]; R_end <- params[['R_end']]
  R_val <- params[['R_val']]; study_len <- params[['study_len']]
  num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]

  #Determine Rt
  if (increasingR == TRUE){
    Rt <- seq(R_start, R_end, length.out = study_len)
  }

  if(increasingR == FALSE){
    Rt <- rep(R_val, study_len)
  }

  sims <- num_people

  # Generate Distribution for serial interval
  omega <- gen_distribution(study_len, sim_mu, sim_var, sim_type, delta)
  
  #Generate daily Lambda
  omega_sum <- zoo::rollapply(omega$omega, delta, sum, by = delta)
  Lambda <- c() 
  for (d in 1:study_len){
  Lambda[d] <- Rt * zoo::rollapply(omega$omega,  / delta
  }

  func <- function(t) rep(t, sum(rpois(sims, Lambda[t])))
  samplesDaily <- do.call(c, sapply(1:study_len, func))

  # Get output that includes true distribution and simulated secondary cases
  serinfect <- list(daily = samplesDaily, dist = omega$omega)
  return(serinfect)

}

test <- samp_pois(params, increasingR = FALSE)
plot(test$dist, type = "l")
serial_ests(test$daily)
