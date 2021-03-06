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
library(zoo)
library(msm)
library(ggridges)
#setwd("..")
imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

# source(paste0(getwd(), "/Code/Params.R"))

###############################################################
############ Generate 'true' distribution #####################
###############################################################

# helper function for gamma dist
# cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, rate = b)

gen_distribution <- function(study_len, mean, variance, type, delta) {

  r <- 0:(study_len*delta-1)
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

###############################################################
#################### Generate 'true' Rt #######################
###############################################################

Rt_gen <- function(R_type, n_days){

  if (R_type == "constant"){
    Rt <- rep(1.8, n_days)
  }

  else if (R_type == "increasing"){
    Rt <- seq(1.5, 3, length.out = n_days)
  }

  else if (R_type == "decreasing"){
    Rt <- seq(3, 1.5, length.out = n_days)
  }

  else if (R_type == "panic"){
    panic_func <- function(x){1.079069 + 0.227532*x - 0.006662227*x^2 + 0.00006154452*x^3 - 1.795452e-7*x^4}
    Rt <- panic_func(1:n_days)
  }

  else if (R_type == "cave"){
    cave_func <- function(x){0.5909007 + 0.1099206*x - 0.0008213363*x^2}
    Rt <- cave_func(1:n_days)
  }

  #Make sure there are no negative values
  for (t in seq_along(Rt)){
      ifelse(Rt[t] < 0, Rt[t] <- 0, Rt[t] <- Rt[t])
  }

  return(Rt)
}


###############################################################
########### Simulate Serial Interval Data #####################
###############################################################

samp_pois <- function(params){

  R_val <- params[['R_val']]; study_len <- params[['study_len']]
  num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; R_type <- params[['R_type']]
  n_days <- params[["n_days"]]

  Rt <- Rt_gen(R_type, n_days)
  Rt <- rep(mean(Rt[study_len:(2*study_len)]), study_len)

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
  kernel_est <- bkde(samps, bandwidth = h, gridsize = max(range), range.x = range)
}

###############################################################
################ Simulate Incidence Data ######################
###############################################################

nour_sim_data <- function(params) {

  R <- params[['R_val']]; sim_mu <- params[['sim_mu']]
  sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  delta <- params[['delta']]; tau_m <- params[['tau_m']]
  n_days <- params[['n_days']]; init_infec <- params[['init_infec']]
  R_type <- params[['R_type']]

  # Setup
  data <- data.frame(matrix(nrow = n_days * delta, ncol = 4))
  colnames(data) <- c('index', 'days', 'R_t', 'infected')
  data$index <- 1:dim(data)[1]
  data$days <- rep(1:n_days, each = delta)
  data$R_t <- rep(Rt_gen(R_type, n_days), each = delta)

  # Generate cases for Burn-in
  data$infected[1:(tau_m * delta)] <- init_infec

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta)$omega)

  # simulate outbreak
  # start of outbreak
  start <- ((tau_m) * delta) + 1

  pb <- txtProgressBar(min = start, max = dim(data)[1], style = 3)
  for (n in start:dim(data)[1]) {

    I_vec <- data[data$index %in% (n - tau_m*delta):(n - 1),"infected"]
    total_infec <- sum(I_vec * omega) / delta
    I_dot <- data$R_t[n] * total_infec

    data$infected[n] <- I_dot

    #setTxtProgressBar(pb, n)
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day = round(sum(infected)), R_val = mean(R_t), .groups = 'drop')

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

  estimates <- list(distribution = estimates, avg_params = avg_params, Vg = deltag)
  return(estimates)

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
                                   scale = 3, rel_min_height = 0.01, size = 0.1) +
      theme_ridges() + theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
      geom_point(data = omega, aes(Y, X), color = 'darkred', shape = 18, size = 3) +
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
  plot(true_dist, type = 'l', ylim=c(0,max(conf$upper)), ann = F)
  lines(conf$lower, lty = "dotted", col = "blue")
  lines(conf$upper, lty = "dotted", col = "blue")
}


###############################################################
#################### Estimate Rt ##############################
###############################################################

Rt_est <- function(df, vals, params, deterministic = F, correct_bias = F, variant = F, sep_Rt = F, sep_S = F) {
  
  n_days <- params[['n_days']]; start <- params[['study_len']] + 1 # should start at start or start +1?
  type <- params[['sim_type']]; start_variant <- params[["start_variant"]] +1
  tau_m <- params[['tau_m']]; R_type <- params[['R_type']]; delta <- params[['delta']]

  # note: df contains daily infections, nrow df is number of days

  # different Rts?
  if (sep_Rt) {
    data <- data.frame(matrix(nrow = n_days, ncol = 7))
    names(data) <- c('Date', 'est_a', 'est_b', 'Rt', 'Est_Rt', 'Est_Rt1', 'Est_Rt2')
  } else{
    data <- data.frame(matrix(nrow = n_days, ncol = 5))
    names(data) <- c('Date', 'est_a', 'est_b', 'Rt', 'Est_Rt')
  }

  # R values
  if (variant){
    data$Rt1 <- rep(df$R1_val)  # TODO: should be through RT_gen??
    data$Rt2 <- rep(df$R2_val)
  }
  else{
     data$Rt <- Rt_gen(R_type, n_days)
  }

  data$Date <- seq(1, n_days)
  data$Est_Rt <- 1
  data[data$Date < start, ]$Est_Rt <- NA
  if (sep_Rt) {
    data$Est_Rt1 <- 1
    data$Est_Rt2 <- 1
    data[data$Date < start, ]$Est_Rt1 <- NA
    data[data$Date < (start_variant + tau_m), ]$Est_Rt2 <- NA
  }

  if (deterministic){
    mean <- params[['sim_mu']]
    var <- params[['sim_var']]

    for (i in start:nrow(data)) {
      if (data[i,]$Date >= start) {
        t <- data[i,]$Date

        omega <- rev(gen_distribution(tau_m, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% (t - tau_m):(t - 1),"infected_day"]
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
        if (sep_Rt) {
          I1 <- df[which(df$days == t),]$I1_daily
          I2 <- df[which(df$days == t),]$I2_daily
          I1_vec <- df[df$days %in% (t - tau_m):(t - 1),]$I1_daily
          I2_vec <- df[df$days %in% (t - tau_m):(t - 1),]$I2_daily
          data$Est_Rt1[i] <- (I1) / (sum(I1_vec * omega))
          if (t > (start_variant + tau_m)) {
            data$Est_Rt2[i] <- (I2) / (sum(I2_vec * omega))
          }
        }
      }
    }
  }

  else{
    data$est_a <- rep(vals$shape, n_days)
    data$est_b <- rep(vals$rate, n_days)

    for (i in start:nrow(data)) {
      if (data[i,]$Date < start) {
        data[i,]$Est_Rt <- NA
        }
      if (sep_Rt) {
        if (data[i,]$Date < start) {
          data$Est_Rt1[i] <- NA
        }
        if (data$Date[i] < (start_variant + tau_m)) {
          data$Est_Rt2[i] <- NA
        }
      }
      else {
        #a <- data[i,]$est_a
        #b <- data[i,]$est_b
        mean <- vals$meanhat
        var <- vals$varhat
        t <- data[i,]$Date

        omega <- rev(gen_distribution(tau_m, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% (t - tau_m):(t - 1),"infected_day"]
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
        if (sep_Rt) {
          I1 <- df[which(df$days == t),]$I1_daily
          I2 <- df[which(df$days == t),]$I2_daily
          I1_vec <- df[df$days %in% (t - tau_m):(t - 1),]$I1_daily
          I2_vec <- df[df$days %in% (t - tau_m):(t - 1),]$I2_daily
          data$Est_Rt1[i] <- (I1) / (sum(I1_vec * omega))
          if (t > (start_variant + tau_m)) {
            data$Est_Rt2[i] <- (I2) / (sum(I2_vec * omega))
          }
        }
      }
    }
  }

  if (sep_Rt) {
    df$I1_Prop <- df$I1_daily / df$infected_day
    df$I2_Prop <- df$I2_daily / df$infected_day
    df$I2_Prop[is.na(df$I2_Prop)] <- 0
    data$Rt_Avg <- (df$R1_val * df$I1_Prop) + (df$R2_val * df$I2_Prop)
  }

  if (correct_bias){
    if (sep_Rt) {
      if (sep_S) {
        data$Est_Rt1 <- data$Est_Rt1 * 1/df$S1_pct
        data$Est_Rt2 <- data$Est_Rt2 * 1/df$S2_pct
      }
      else {
        data$Est_Rt1 <- data$Est_Rt1 * 1/df$S_pct
        data$Est_Rt2 <- data$Est_Rt2 * 1/df$S_pct
        data$Est_Rt <- data$Est_Rt * 1/df$S_pct
      }
    }
    else {
      data$Est_Rt <- data$Est_Rt * 1/df$S_pct
    }
  }
  return(data)
}

Rt_est_nonpara <- function(df, samps, bw, params, correct_bias = F, variant = F) {
  start <- params$study_len +1
  n_days <- params$n_days
  study_len <- params$study_len
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
    if (data$Date[i] < start) {
      data$Est_Rt[i] <- NA
    }
    else {
      t <- data$Date[i]

      dist <- rev(serial_ests_nonpara(samps, range = c(1, study_len), bandwidth = bw)$y)
      I <- df[which(df$days == t),]$infected_day
      I_window <- df[df$days %in% (t - study_len):(t - 1),"infected_day"]
      data$Est_Rt[i] <- (I) / (sum(I_window * dist))
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