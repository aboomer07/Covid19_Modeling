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

# note: R(t) in serial interval study needs to be sufficient large to see something
#       since it reflects number of how many persons an infected is expected to infect
#       (around 4 is realistic)

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

  # Generate daily Lambda
  omega_sum <- zoo::rollapply(omega$omega, width = delta, FUN = sum, by = delta)
  Lambda <- Rt * omega_sum /delta

  # repeat the respective day of infection as many times as there are infections
  func <- function(t) rep(t, sum(rpois(sims, Lambda[t])))
  samplesDaily <- do.call(c, sapply(1:study_len, func))

  # Get output that includes true distribution and simulated secondary cases
  serinfect <- list(daily = samplesDaily, dist = omega$omega)
  return(serinfect)

}





nour_sim_data <- function(params) {

  R <- params[['R_val']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; tau_m <- params[['tau_m']]
  days <- params[['n_days']]; init_infec <- params$init_infec

  # Setup
  data <- data.frame(matrix(nrow = days * delta, ncol = 4))
  colnames(data) <- c('index', 't', 'R_t', 'I_dot')
  data$index <- 1:dim(data)[1]
  data$t <- rep(1:days, each = delta)
  data$R_t <- rep(R, each = (delta * days / length(R))) #TODO

  # Generate cases for Burn-in
  data$infected[1:(tau_m * delta)] <- init_infec

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta)$omega)

  # simulate outbreak
  # start of outbreak
  start <- ((tau_m) * delta) + 1

  for (n in start:dim(data)[1]) {

    I_vec <- data[data$index %in% (n - tau_m*delta):(n - 1),"infected"]
    total_infec <- sum(I_vec * omega) / delta

    data$I_dot <- R[t] * total_infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(t) %>%
    summarise(infected = round(sum(I_dot)), R_val = mean(R_t))

  return(daily_infec)
}
