# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/26/2021

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


###############################################################
########### Simulate Serial Interval Data #####################
###############################################################

samp_pois <- function(params) {

  R_val <- params[['R_val']]; study_len <- params[['study_len']]
  num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]

  Rt <- rep(R_val, each = (study_len*delta/length(R_val))+1)
  sims <- num_people

  # Generate Distribution for serial interval
  omega <- gen_distribution(study_len, sim_mu, sim_var, sim_type, delta)
  #Generate lambda parameter for the poisson draw
  Lambda <- Rt * sum(omega$omega) / delta




  range <- 1:(study_len*delta)

  func <- function(t) rep(t, sum(rpois(sims, Lambda[t])))
  samplescont <- do.call(c, sapply(range, func))

  # Make infections daily
  daily <- floor(samplescont / delta) + 1   # round discretized infections in samplescont to day (day 1 = 0-1)


  # Get output that includes true distribution and simulated secondary cases
  serinfect <- list(samplescont = samplescont, daily = daily, dist = omega$omega)
  return(serinfect)
}