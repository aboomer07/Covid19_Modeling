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

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

###############################################################
############ Generate 'true' distribution #####################
###############################################################

gen_params <- function(mean, variance, type, delta) {

  if (type == 'norm') {
    a <- mean * delta
    b <- variance * delta
    return(c(a, b))
  }

  if (type == 'lnorm') {
    a <- log(mean^2 / (sqrt(mean^2 + variance))) + log(delta)
    b <- log(1 + variance / mean^2)
    return(c(a, b))
  }

  if (type == 'gamma') {
    a <- mean^2 / variance
    b <- (mean / variance) / delta
    return(c(a, b))
  }

  if (type == 'weibull') {
    a <- as.numeric(weibullpar(mean, variance)[1])
    b <- as.numeric(weibullpar(mean, variance)[2]) * delta
    return(c(a, b))
  }

}

gen_distribution <- function(k, mean, var, type, delta) {
  k <- 1:(k * (1 / delta))

  params <- gen_params(mean, var, type, delta)
  a <- params[1]
  b <- params[2]

  if (type == 'gamma') {
    cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

    omega <- k * cdf_gamma(k, a, b) +
      (k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
    omega <- omega + a *
      b *
      (2 * cdf_gamma(k - 1, a + 1, b) -
        cdf_gamma(k - 2, a + 1, b) -
        cdf_gamma(k, a + 1, b))
    omega <- sapply(omega, function(e) max(0, e))

    return(omega)
  }

  if (type == 'weibull') {
    omega <- dweibull(k, a, b)
    return(omega)
  }

  if (type == 'norm') {
    omega <- dnorm(k, a, b)
    return(omega)
  }

  if (type == 'lnorm') {
    omega <- dlnorm(k, a, b)
    return(omega)
  }

}

###############################################################
########### Simulate Serial Interval Data #####################
###############################################################
samp_pois <- function(R_val, study_len, num_people, sim_mu, sim_sig, sim_type, delta) {

  data <- data.frame(rep(R_val, (study_len / length(R_val))))
  names(data) <- "Rt"
  data$sim_mu <- rep(sim_mu, study_len)
  data$sim_sig <- rep(sim_sig, study_len)
  data$Lambda <- rep(0, study_len)

  sims <- num_people
  samples <- c()

  range <- 1:nrow(data)

  #Generate Mean Incidence based on gamma prior
  for (t in range) {

    dist <- gen_distribution(t, sim_mu, sim_sig, sim_type, delta)
    data[t,]$Lambda <- data[t,]$Rt * dist[t]

    #Add up all the random poisson infections
    for (sim in 1:sims) {
      samples <- c(samples, rep(t, rpois(1, data[t,]$Lambda)))
    }
  }
  return(samples)
}


###############################################################
################ Estimate Serial Interval #####################
###############################################################
serial_ests <- function(samps) {
	num_vars <- 10
	params <- list()

	#fit  distribution with a gamma
    fit <- fitdist(samps, "gamma")

    a_est <- fit$estimate[1]
    b_est <- fit$estimate[2]

    a_sd <- fit$sd[1]
    b_sd <- fit$sd[2]

    mean <- a_est / b_est
    var <- a_est * (b_est^2)
    mean <- mean - 1
    b_est <- mean / sqrt(var)
    a_est <- mean / b_est

    lower_a <- a_est - (1.96 * a_sd)
    upper_a <- a_est + (1.96 * a_sd)

    lower_b <- b_est - (1.96 * b_sd)
    upper_b <- b_est + (1.96 * b_sd)

    vals <- expand.grid(seq(lower_a, upper_a, ((upper_a - lower_a) / num_vars)),
                        seq(lower_b, upper_b, ((upper_b - lower_b) / num_vars)))

    names(vals) <- c("Param_a", 'Param_b')
    de <- data.frame(a_est, b_est)
    names(de) <- c("Param_a", 'Param_b')
    vals <- rbind(vals, de)
    vals$True_a <- a_est
    vals$True_b <- b_est

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
  data <- data.frame(matrix(nrow = n_days * nrow(vals), ncol = 6))
  names(data) <- c('Date', 'est_a', 'est_b', 'est_type', 'Rt', 'Est_Rt')

  data$Date <- rep(1:n_days, nrow(vals))
  data$est_a <- rep(vals[, 1], each = n_days)
  data$est_b <- rep(vals[, 2], each = n_days)
  data$est_type <- rep(vals[, 3], each = n_days)
  data$Rt <- rep(df['R_val'][[1]], nrow(vals))
  data$True_est_a <- rep(vals[, 4], each = n_days)
  data$True_est_b <- rep(vals[, 5], each = n_days)

  for (i in start:nrow(data)) {
    if (data[i,]$Date <= start) {
      data[i,]$Est_Rt <- NA
    }
    else {
      a <- data[i,]$est_a
      b <- data[i,]$est_b
      type <- data[i,]$est_type
      t <- data[i,]$Date

      dist <- rev(gen_distribution(t - 1, a, b, type))
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
    group_by(est_a, est_b, est_type) %>%
    summarize(MSE = mean(SSE, na.rm = TRUE), True_est_a = mean(True_est_a), True_est_b = mean(True_est_b))

  return(df)
}