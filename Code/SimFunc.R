# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(fitdistrplus)
library(RColorBrewer)
library(data.table)
library(np)
library(KernSmooth)
library(mixdist)



###############################################################
############ Generate 'true' distribution #####################
###############################################################

gen_distribution <- function(k, mean, variance, type, delta) {

  k <- 1:(k*delta)

  if (type == 'norm') {
    a <- mean * delta
    b <- variance * delta
    print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dnorm(k, a, b)
  }

  if (type == 'lnorm') {
    a <- log(mean^2 / (sqrt(mean^2 + variance))) + log(delta)
    b <- log(1 + variance / mean^2)
    print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dlnorm(k, a,b)
  }

  if (type == 'gamma') {
    a <- mean^2 / variance
    b <- (mean / variance) / delta
    print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dgamma(k, a, b)    
  }

  if (type == 'weibull') {
    a <- as.numeric(weibullpar(mean, variance)[1])
    b <- as.numeric(weibullpar(mean, variance)[2]) * delta
    print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dweibull(k, a, b)
  }

 return(omega)
}


###############################################################
########### Simulate Serial Interval Data #####################
###############################################################

samp_pois <- function(R_val, study_len, num_people, sim_mu, sim_sig, sim_type, delta) {

  data <- data.frame(rep(R_val, each = (study_len*delta / length(R_val))))
  names(data) <- "Rt"
  data$Lambda <- rep(0, study_len*delta)

  sims <- num_people
  samplescont <- c()

  range <- 1:nrow(data)

  # Generate Distribution for serial interval
  dist <- gen_distribution(study_len, sim_mu, sim_sig, sim_type, delta)
  for (t in range) {

    # get mean infectivity for respective day of study
    data[t,]$Lambda <- data[t,]$Rt * dist[t]

    #Add up all the random poisson infections
    for (sim in 1:sims) {
      samplescont <- c(samplescont, rep(t, rpois(1, data[t,]$Lambda)))
    }
  }

  # Make infections daily
  daily <- c()
  day <- seq(1, study_len*delta, delta)
  
  for (d in 1:length(day)) {
    for (i in 1:length(samplescont)) {
      if ((samplescont[i] >= day[d]) & (samplescont[i] < day[d+1])) {
        daily <- c(daily, d)
      }
    }
  }

  # Get output that includes true distribution and simulated secondary cases
  serinfect <- list(samplescont = samplescont, daily = daily, dist = dist)
  return(serinfect)
}


###############################################################
################ Estimate Serial Interval #####################
###############################################################
serial_ests <- function(samps) {
  num_vars <- 5
  params <- list()

  #fit  distribution with a gamma
  fit <- fitdist(samps, "gamma")

  a_est <- as.numeric(fit$estimate[1])
  b_est <- as.numeric(fit$estimate[2])

  # a_sd <- as.numeric(fit$sd[1])
  # b_sd <- as.numeric(fit$sd[2])
  
  #Transform parameters into mean and variance
  a_norm <- as.numeric(a_est/b_est)
  b_norm <- as.numeric(a_est/(b_est^2))
  
  # mean <- a_est / b_est
  # var <- a_est * (b_est^2)
  # mean <- mean - 1
  # b_est <- mean / sqrt(var)
  # a_est <- mean / b_est

  #lower_a <- a_est - (1.96 * a_sd)
  #upper_a <- a_est + (1.96 * a_sd)

  #lower_b <- b_est - (1.96 * b_sd)
  #upper_b <- b_est + (1.96 * b_sd)

  #vals <- expand.grid(seq(lower_a, upper_a, ((upper_a - lower_a) / num_vars)),
  #seq(lower_b, upper_b, ((upper_b - lower_b) / num_vars)))

  #names(vals) <- c("Param_a", 'Param_b')
  #de <- data.frame(a_est, b_est)
  #names(de) <- c("Param_a", 'Param_b')
  #vals <- rbind(vals, de)
  #vals$True_a <- a_est
  #vals$True_b <- b_est

  vals <- data.frame(a_est, b_est, a_norm, b_norm)
  names(vals) <- c("shape", "rate", "mean", "variance")
  return(vals)
}

serial_ests_nonpara <- function(samps, correction, range) {
  h_nsr <- 1.059 * sd(samps) * length(samps)^(-1 / 5)
  kernel_est <- bkde(samps, bandwidth = h_nsr + correction, range.x = range, gridsize = max(range))
}

###############################################################
################ Simulate Incidence Data ######################
###############################################################
nour_sim_data <- function(sim_mu, sim_var, sim_type, delta, days = n_days, tau_m = window, R = R_val, N = n) {

  data <- data.frame(matrix(nrow = days * delta, ncol = 4))
  colnames(data) <- c('t', 'days', 'R_t', 'infective')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$infective[1:(tau_m * delta)] <- round(seq(10, 100, length.out = tau_m * delta))
  data$R_t <- rep(R, each = (delta * days / length(R)))

  dist <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta))

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {
    I_vec <- data[data$t %in% (t - N):(t - 1),]$infective
    R_mean <- mean(data[which(data$t == t),]$R_t)
    total_infec <- sum(I_vec * dist)

    infec <- R_mean * total_infec
    data[which(data$t == t),]$infective <- infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infective_day = round(mean(infective)), R_val = mean(R_t))

  return(daily_infec)
}

###############################################################
#################### Estimate Rt ##############################
###############################################################
Rt_est <- function(df, vals, type) {
  start <- 10
  data <- data.frame(matrix(nrow = n_days * nrow(vals), ncol = 5))
  names(data) <- c('Date', 'est_a', 'est_b', 'Rt', 'Est_Rt')

  data$Date <- rep(1:n_days, nrow(vals))
  data$est_a <- rep(vals$shape, each = n_days)
  data$est_b <- rep(vals$rate, each = n_days)
  data$Rt <- rep(df['R_val'][[1]], nrow(vals))
  #data$True_est_a <- rep(vals[, 3], each = n_days)
  #data$True_est_b <- rep(vals[, 4], each = n_days)

  for (i in start:nrow(data)) {
    if (data[i,]$Date <= start) {
      data[i,]$Est_Rt <- NA
    }
    else {
      a <- data[i,]$est_a
      b <- data[i,]$est_b
      mean <- a/b
      var <- a/(b^2)
      t <- data[i,]$Date

      dist <- rev(gen_distribution(t - 1, mean, var, type, 1))
      I <- df[which(df$days == t),]$infective_day
      I_window <- df[df$days %in% 1:(t - 1),]$infective_day
      data[i,]$Est_Rt <- (I) / (sum(I_window * dist))
    }
  }

  return(data)
}

Rt_est_nonpara <- function(df, samps, corrections) {
  start <- 10
  data <- data.frame(matrix(nrow = n_days * length(corrections), ncol = 4))
  names(data) <- c('Date', 'Cor_Par', 'Rt', 'Est_Rt')

  data$Date <- rep(1:n_days, length(corrections))
  data$Cor_Par <- rep(corrections, each = n_days) # instead of distribution parameters we can correct/opt bandwidth
  data$Rt <- rep(df['R_val'][[1]], length(corrections))

  for (i in start:nrow(data)) {
    if (data[i,]$Date <= start) {
      data[i,]$Est_Rt <- NA
    }
    else {
      correction <- data[i,]$Cor_Par
      t <- data[i,]$Date

      dist <- rev(serial_ests_nonpara(samps, correction, range = c(1, (t - 1)))$y)
      I <- df[which(df$days == t),]$infective_day
      I_window <- df[df$days %in% 1:(t - 1),]$infective_day
      data[i,]$Est_Rt <- (I) / (sum(I_window * dist))
    }
  }
  return(data)
}

MSE_est <- function(df) {

  df$SSE <- (df$Rt - df$Est_Rt)^2
  df <- df %>%
    group_by(est_a, est_b) %>%
    summarize(MSE = mean(SSE, na.rm = TRUE), True_est_a = mean(True_est_a), True_est_b = mean(True_est_b))

  return(df)
}