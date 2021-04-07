# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(tidyverse)
library(fitdistrplus)
library(RColorBrewer)
library(data.table)
library(np)
library(KernSmooth)
library(mixdist)
library(msm)
library(ggridges)
imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

# source(paste0(getwd(), "/Code/Params.R")

###############################################################
############ Generate 'true' distribution #####################
###############################################################

# helper function for gamma dist
# cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, rate = b)

gen_distribution <- function(k, mean, variance, type, delta) {

  k <- 1:(k*delta)

  if (type == 'norm') {
    a <- mean * delta
    b <- variance * delta
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dnorm(k, a, b)
  }

  if (type == 'lnorm') {
    a <- log(mean^2 / (sqrt(mean^2 + variance))) + log(delta)
    b <- log(1 + variance / mean^2)
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dlnorm(k, a,b)
  }

  if (type == 'gamma') {
    a <- mean^2 / variance
    b <- (mean / variance) / delta
    # b <- 1/b
    #print(paste("Distribution:", type, "Serial interval parameters: a =",a, "b =", b))
    omega <- dgamma(k, a, b)
  }

  if (type == 'weibull') {
    a <- as.numeric(weibullpar(mean, variance)[1])
    b <- as.numeric(weibullpar(mean, variance)[2])*delta
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

  Rt <- rep(R_val, each = (study_len*delta/length(R_val)))
  sims <- num_people

  # Generate Distribution for serial interval
  omega <- gen_distribution(study_len, sim_mu, sim_var, sim_type, delta)
  #Generate lambda parameter for the poisson draw
	Lambda <- Rt * omega$omega

  range <- 1:(study_len*delta)

  func <- function(t) rep(t, sum(rpois(sims, Lambda[t])))
  samplescont <- do.call(c, sapply(range, func))

  # Make infections daily
  daily <- c()
  day <- seq(1, study_len*delta+2, delta)
  
  for (d in 1:length(day)) {
    for (i in 1:length(samplescont)) {
      if ((samplescont[i] >= day[d]) & (samplescont[i] < day[d+1])) {
        daily <- c(daily, d)
      }
    }
  }

  # Get output that includes true distribution and simulated secondary cases
  serinfect <- list(samplescont = samplescont, daily = daily, dist = omega$omega)
  return(serinfect)
}
 
###############################################################
################ Estimate Serial Interval #####################
###############################################################
serial_ests <- function(samps) {
  num_vars <- 5
  params <- list()

  #fit  distribution with a gamma
  fit <- fitdist(samps, "gamma", method = "mle")

  a_est <- as.numeric(fit$estimate[1])
  b_est <- as.numeric(fit$estimate[2])
  a_sd <- as.numeric(fit$sd[1])
  b_sd <- as.numeric(fit$sd[2])

  a_cov <- as.numeric(fit$vcov[1,2])
  a_norm <- as.numeric(a_est/b_est)
  b_norm <- as.numeric(a_est/(b_est^2))
 
  vals <- list(shape = a_est, rate = b_est, shape_sd = a_sd, rate_sd = b_sd, cov = a_cov, meanhat = a_norm, varhat = b_norm)
  return(vals)
}

serial_ests_nonpara <- function(samps, range, bandwidth) {
  if (bandwidth == 'nsr'){
    h <- 1.059 * sd(samps) * length(samps)^(-1 / 5)
  }
  if (bandwidth == 'iqr'){
    h <- 0.9 * length(samps)^(-1 / 5) * (IQR(samps)/1.34)
  }
  kernel_est <- bkde(samps, bandwidth = h, range.x = range, gridsize = max(range))
}

###############################################################
################ Simulate Incidence Data ######################
###############################################################
nour_sim_data <- function(params) {

  R <- params[['R_val']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; tau_m <- params[['tau_m']]
  days <- params[['n_days']]

  data <- data.frame(matrix(nrow = days * delta, ncol = 4))
  colnames(data) <- c('t', 'days', 'R_t', 'infected')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$infected <- NA
  data$infected[1:(tau_m * delta)] <- round(seq(10, 100, length.out = tau_m * delta))
  data$R_t <- rep(R, each = (delta * days / length(R)))

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta))

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {
    I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$infected
    R_mean <- mean(data[which(data$t == t),]$R_t)
    total_infec <- sum(I_vec * omega)

    infec <- R_mean * total_infec
    data[which(data$t == t),]$infected <- infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day = round(mean(infected)), R_val = mean(R_t))

  return(daily_infec)
}


####################################################################################
#################### Simulate several times to get distribution ####################
####################################################################################

params_distribution <- function(params) {

  R_val <- params[['R_val']]; study_len <- params[['study_len']]
  num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; simulations <- params[['simulations']]

  estimates <- data.frame(ncol = 4, nrow = simulations)

  for (i in 1:simulations){
    serinfect <- samp_pois(params)
    vals <- serial_ests(serinfect$daily) # here we obtain the params for Rt_est
    print(i)
    estimates[i, 1] <- vals$shape
    estimates[i, 2] <- vals$rate
    estimates[i, 3] <- vals$shape_sd
    estimates[i, 4] <- vals$rate_sd
    estimates[i, 5] <- vals$cov
    estimates[i, 6] <- vals$meanhat
    estimates[i, 7] <- vals$varhat
    estimates[i, 8] <- length(serinfect$daily)
    names(estimates) <- c("shape_hat", "rate_hat", "shape_sd", "rate_sd",
               "cov_hat", "mean_hat", "var_hat", "n")
    }

  #Delta method to get standard errors of the mean and the variance
  deltag <- matrix(ncol = 2, nrow = simulations)

  for (i in 1:simulations) {
    gradient <- matrix(ncol = 2, nrow = 2)
    gradient[1,1] <- 1/estimates$rate_hat[i]
    gradient[1,2] <- -estimates$shape_hat[i]/estimates$rate_hat[i]^2
    gradient[2,1] <- 1/estimates$rate_hat[i]^2
    gradient[2,2] <- -2*estimates$shape_hat[i]/estimates$rate_hat[i]^3

    varcov <- matrix(ncol = 2, nrow =2)
    varcov[1,2] <- estimates$cov[i]
    varcov[2,1] <- estimates$cov[i]
    varcov[1,1] <- estimates$shape_hat[i]^2
    varcov[2,2] <- estimates$rate_hat[i]^2

    temp.matrix <- t(gradient) %*% varcov %*% gradient
    deltag[i, 1] <- sqrt(temp.matrix[1,1]/estimates$n[i])
    deltag[i, 2] <- sqrt(temp.matrix[2,2]/estimates$n[i])

  }

  #Store parameters of interest
  avg_params <- data.frame()[1, ]
  avg_params$avg_shape_hat <- mean(estimates$shape_hat)
  avg_params$avg_shape_sd <- mean(estimates$shape_sd)
  avg_params$avg_rate_hat <- mean(estimates$rate_hat)
  avg_params$avg_rate_sd <- mean(estimates$rate_sd)
  avg_params$avg_mean_hat <- mean(estimates$mean_hat)
  avg_params$avg_mean_sd <- mean(deltag[1])
  avg_params$avg_var_hat <- mean(estimates$var_hat)
  avg_params$avg_var_sd <- mean(deltag[2])

  estimates <- list(distribution = estimates, avg_params = avg_params)

}


nonpara_eval <- function(params, bw){

  R_val <- params[['R_val']]; study_len <- params[['study_len']]
  num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; sims <- params[['simulations']]

  simulations <- list()
  for (i in 1:sims){
    simulation <- data.frame(matrix(ncol = 4, nrow = study_len))
    serinfect <- samp_pois(params)

    out <- serial_ests_nonpara(serinfect$daily, range = c(1, study_len), bandwidth = bw)

    simulation[ ,1] <- out$x
    simulation[, 2] <- out$y
    simulation[, 3] <- i
    simulation[, 4] <- bw
    simulations[[i]] <- simulation
  }
  full <- do.call('rbind', simulations)
  #omega <- gen_distribution(study_len, sim_mu, sim_var, sim_type, 1)$omega
  #full$Omega <- rep(omega, study_len*sims)
  colnames(full) <- c('X', 'Y', 'Sim', 'Bandwidth')
  return(full)
}

plot_nonpara_distplot <- function(est_dist, plot_type, params){
  omega <- data.frame(matrix(nrow = params[['study_len']], ncol = 2))
  colnames(omega) <- c('X', 'Y')
  omega$Y <- gen_distribution(params[['study_len']], params[['sim_mu']], params[['sim_var']], params[['sim_type']], 1)$omega
  omega$X <- 1:params[['study_len']]
  est_dist$X <- as.factor(est_dist$X)
  if (plot_type == 'ridge'){
    ggplot() +
      geom_density_ridges_gradient(data = est_dist, aes(Y, X),
                                   scale = 3, rel_min_height = 0.01, size = 0.3) +
      theme_ridges() + theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
      geom_point(data = omega, aes(Y, X), color = 'green', shape = 5) +
      xlim(0, 0.35) +
      coord_flip() +
      ylab('Days') + xlab('Estimated Y-values') +
      ggsave(paste0(outpath, 'SerialEst_nonpara_ridgeplot_', length(est_dist$X), '.png'), width = 10, height = 5)
  }
  if (plot_type == 'violin'){
    ggplot(est_dist, aes(X, Y)) +
      geom_violin() + theme_minimal() +
      theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
      ylab('') + xlab('Days')
      ggsave(paste0(outpath, 'SerialEst_nonpara_violinplot_', length(est_dist$X), '.png'), width = 10, height = 5)
  }
}


plot_nonpara_eval <- function(true_dist, est_dist) {
  conf <- est_dist %>% group_by(X) %>%
    summarise(upper = max(Y),
              lower = min(Y))
  png(paste0(outpath, 'SerialEst_nonpara_', length(est_dist$X), '.png'))
  plot(true_dist, type = 'l', ylim=c(0,max(conf$upper)))
  lines(conf$lower, lty = "dotted", col = "blue")
  lines(conf$upper, lty = "dotted", col = "blue")
  dev.off()
}


###############################################################
#################### Estimate Rt ##############################
###############################################################
Rt_est <- function(df, vals, params, deterministic = F, correct_bias = F, variant = F) {
  n_days <- params[['n_days']]
  start <- params[['study_len']]
  type <- params[['sim_type']]
  data <- data.frame(matrix(nrow = n_days, ncol = 5))
  names(data) <- c('Date', 'est_a', 'est_b', 'Rt', 'Est_Rt')

  if (variant){
    data$Rt1 <- rep(df$R1_val)
    data$Rt2 <- rep(df$R2_val)
  }
  else{
     data$Rt <- rep(df$R_val)
  }

  data$Date <- seq(1, n_days)

  if (deterministic){
    mean <- params[['sim_mu']]
    var <- params[['sim_var']]

    for (i in start:nrow(data)) {
      if (data[i,]$Date <= start) {
        data[i,]$Est_Rt <- NA
      }
      else {
        t <- data[i,]$Date

        omega <- rev(gen_distribution(t - 1, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% 1:(t - 1),]$infected_day
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
      }
    }
  }

  else{
    data$est_a <- rep(vals$shape, n_days)
    data$est_b <- rep(vals$rate, n_days)

    for (i in start:nrow(data)) {
      if (data[i,]$Date <= start) {
        data[i,]$Est_Rt <- NA
      }
      else {
        #a <- data[i,]$est_a
        #b <- data[i,]$est_b
        mean <- vals$meanhat
        var <- vals$varhat
        t <- data[i,]$Date

        omega <- rev(gen_distribution(t - 1, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% 1:(t - 1),]$infected_day
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
      }
    }
  }
  if (correct_bias){
    data$Est_Rt <- data$Est_Rt * 1/df$S_pct
  }
  return(data)
}

Rt_est_nonpara <- function(df, samps, bw, params, correct_bias = F, variant = F) {
  start <- params[['study_len']]
  n_days <- params[['n_days']]
  data <- data.frame(matrix(nrow = n_days, ncol = 3))
  names(data) <- c('Date', 'Rt', 'Est_Rt')

  data$Date <- rep(1:n_days)

  if (variant){
    data$Rt1 <- rep(df$R1_val)
    data$Rt2 <- rep(df$R2_val)
  }

  else{
    data$Rt <- rep(df['R_val'][[1]])
  }

  for (i in start:nrow(data)) {
    if (data[i,]$Date <= start) {
      data[i,]$Est_Rt <- NA
    }
    else {
      t <- data[i,]$Date

      dist <- rev(serial_ests_nonpara(samps, range = c(1, (t - 1)), bandwidth = bw)$y)
      I <- df[which(df$days == t),]$infected_day
      I_window <- df[df$days %in% 1:(t - 1),]$infected_day
      data[i,]$Est_Rt <- (I) / (sum(I_window * dist))
    }
  }
  if (correct_bias){
    data$Est_Rt <- data$Est_Rt * 1/df$S_pct
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


MSE_est2 <- function(df) {
  
  df$SSE <- (df$true - df$est)^2
  MSE<-mean(df$SSE)
  
  return(MSE)
}