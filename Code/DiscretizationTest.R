# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/26/2021

library(zoo)

gen_distribution <- function(study_len, mean, variance, type, delta) {

  k <- seq(0, study_len, length.out = study_len * delta)

  if (type == 'norm') {
    a <- mean
    b <- variance
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dnorm(k, a, b)
  }

  if (type == 'lnorm') {
    a <- log(mean^2 / (sqrt(mean^2 + variance)))
    b <- log(1 + variance / mean^2)
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dlnorm(k,a,b)
  }

  if (type == 'gamma') {
    a <- mean^2 / variance
    b <- (mean / variance)
    # b <- 1/b
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dgamma(k, a, b)
  }

  if (type == 'weibull') {
    a <- as.numeric(weibullpar(mean, variance)[1])
    b <- as.numeric(weibullpar(mean, variance)[2])
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dweibull(k, a, b)
  }

  distribution <- list(a = a, b = b, omega = omega)
  return(distribution)
}

omega <- gen_distribution(20, 7, 2, "weibull", 24)


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
  Lambda[d] <- Rt[d] * omega_sum[d] / delta
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
