# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/5/2021

source(paste0(getwd(), "/Code/Params.R"))
source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))

si_sim <- function(params) {
  # unpack variables
  days <- params[['n_days']]; delta <- params[['delta']]
  R <- params[['R_val']]; N <- params[['pop']]; tau_m <- params[['tau_m']]
  sim_mu <- params[['sim_mu']]; sim_var <- params[['sim_var']]
  sim_type <- params[['sim_type']]


  # Setup environment
  data <- data.frame(matrix(nrow = days * delta, ncol = 6))
  colnames(data) <- c('t', 'days', 'R_t', 'infected', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$R_t <- rep(R, each = (delta * days / length(R)))
  data$N <- N

  # start of epidemic (burn-in)
  # Assume 1 initial ifected and that new ifections happen every 4 days depending on specified R (for burn-in)
  data$infected[1] <- 1
  data$S[1] <- data$N[1] - data$infected[1]
  data$I <- data$infected

  # Simulate initial infected spreading disease during burn-in period
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$infected[t-1]
    # data$infected[t] <- sum(data$infected[1:t-1]) * R / (delta*2)
    data$infected[t] <- 10/delta
    data$I[t] <- data$I[t-1] + data$infected[t]
  }

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta)$omega)

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {

    # update susceptible population
    data$S[t] <- data$S[t-1] - data$infected[t-1]

    # compute I(u):
    ## compute I(p) first:
    ## get omega(n*delta) I(t-n*delta) => total_infec (= All infected for a TSI from 1 until tau_m at all deltas)
    I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$infected
    R_mean <- data$R_t[t] #mean(data$R_t[t])
    I_p <- R_mean * omega * I_vec * data$S[t] / data$N[t]
    ## get I(u)
    infec <- sum(I_p)
    data$infected[t] <- infec

    # update total infected
    data$I[t] <- data$I[t-1] + infec
  }

  # plot(data$infected, type='l')
  # abline(v=tau_m*delta, col='red')

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day = sum(infected), R_val = mean(R_t),
              S_daily = last(S), N = last(N))

  daily_infec$I_cum <- cumsum(daily_infec$infected_day)
  daily_infec$S_pct <- daily_infec$S_daily / daily_infec$N
  daily_infec$I_pct <- daily_infec$I_cum / daily_infec$N

  return(daily_infec)
}

## Plot
si_plot <- function (model){
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2, 
    ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I_pct, type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"), 
    col = c("blue", "orange"), pch = 16, bty = "n")
}

si_plot_detail <- function (model){
  layout(matrix(c(rep(1,6), rep(2,3), rep(3,3)), nrow=4, ncol=3, byrow = T))
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2,
    ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab="")
  lines(x = model$days, y = model$I_pct, type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"),
    col = c("blue", "orange"), pch = 20, bty = "n")
  plot(x = model$days, y = model$infected_day, type="l", lwd=2, col = "green",
       ylab = "Infected", xlab = "")
  legend("topright", legend = "New Infections", col="green", pch = 16, bty="n")
  plot(x = model$days, y = model$R_val, type="l", lwd=2, col = "red",
       ylab = "R(t)", xlab = "Days")
  legend("topright", legend = "R(t)", col="red", pch = 16, bty="n")
}


# test
params[['delta']] <- 24
params[['n_days']] <- 300
params[['sim_type']] <- 'weibull'
si_model <- si_sim(params)
si_plot(si_model)
si_plot_detail(si_model)

############################## SII MODEL #######################################
# do the same but add two different infected

sii_sim <- function(params) {

  days <- params[['n_days']]; delta <- params[['delta']]
  R1 <- params[['R_val']]; N <- params[['pop']]; tau_m <- params[['tau_m']]
  sim_mu1 <- params[['sim_mu']]; sim_var1 <- params[['sim_var']]
  sim_type1 <- params[['sim_type']]; sim_mu2 <- params[['sim_mu_variant']]
  sim_var2 <- params[['sim_var_variant']]; R2 <- params[['R_val_variant']]
  start_variant <- params[['start_variant']] + tau_m
  sim_type2 <- params[['sim_type_variant']];

  # setup environment
  # start_variant <- start_variant*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 8))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'I1', 'I2', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$N <- N
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$I1[1:(tau_m * delta)] <- 5
  data[1:(start_variant * delta),]$I2 <- 0
  data$I2[(start_variant * delta+1):(start_variant*delta + (tau_m * delta))] <- 5

  data$S[1] <- data$N[1] - data$I1[1] - data$I2[1]
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]
  }

  # get serial intervals
  omega1 <- rev(gen_distribution(tau_m, sim_mu1, sim_var1, sim_type1, delta)$omega)
  omega2 <- rev(gen_distribution(tau_m, sim_mu2, sim_var2, sim_type2, delta)$omega)

  # TODO Finish model

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {
    if (t < (start_variant * delta + tau_m * delta + 1)){
      data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]

      I1_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I1
      data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S[t] / data$N[t]
    }

    else{
      data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]

      I1_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I1
      data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S[t] / data$N[t]
      I2_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I2
      data$I2[t] <- data$R_t2[t] * sum(I2_vec * omega2) * data$S[t] / data$N[t]
    }
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(I1_daily = sum(I1),
              I2_daily = sum(I2),
              R1_val = mean(R_t1), R2_val = mean(R_t2),
              S_daily = last(S), N = last(N))

  daily_infec$I1_cum <- cumsum(daily_infec$I1_daily)
  daily_infec$I2_cum <- cumsum(daily_infec$I2_daily)
  daily_infec$S_pct <- daily_infec$S_daily / daily_infec$N
  daily_infec$I1_pct <- daily_infec$I1_cum / daily_infec$N
  daily_infec$I2_pct <- daily_infec$I2_cum / daily_infec$N
  daily_infec$infected_day <- daily_infec$I1_daily + daily_infec$I2_daily

  return(daily_infec)
}

## Plot
si_plot <- function (model){
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2,
    ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I1_pct, type = "l", col = "orange", lwd=2)
  lines(x = model$days, y = model$I2_pct, type = "l", col = "green", lwd=2)
  legend("topright",
    legend = c("Susceptible1", "Infected1", "Infected2"),
    col = c("blue", 'red', "orange", 'green'), pch = 16, bty = "n")
}

# test
params[['R_val']] <- 1.5
params[['R_val_variant']] <- 1.9
params[['sim_type']] <- 'gamma'
params[['sim_type_variant']] <- 'gamma'
params[['start_variant']] <- 50
params[['n_days']] <- 150

si_model <- sii_sim(params)
si_plot(si_model)

############################## SSII MODEL ######################################
# do the same but add two different infected + different susceptibility groups

ssii_sim <- function(params) {

  days <- params[['n_days']]; delta <- params[['delta']]
  R1 <- params[['R_val']]; N <- params[['pop']]; tau_m <- params[['tau_m']]
  sim_mu1 <- params[['sim_mu']]; sim_var1 <- params[['sim_var']]
  sim_type1 <- params[['sim_type']]; sim_mu2 <- params[['sim_mu_variant']]
  sim_var2 <- params[['sim_var_variant']]; R2 <- params[['R_val_variant']]
  start_variant <- params[['start_variant']] + tau_m
  sim_type2 <- params[['sim_type_variant']];
  cross <- params[['sii_cross']]

  # setup environment
  # start_variant <- start_variant*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 9))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'I1', 'I2', 'S1', 'S2', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$N <- N
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$I1[1:(tau_m * delta)] <- 5
  data[1:(start_variant * delta),]$I2 <- 0
  data$I2[((start_variant * delta)+1):((start_variant + tau_m) * delta)] <- 5

  data$S1[1] <- data$N[1]
  data$S2[1] <- data$N[1]
  for (t in 2:(tau_m * delta)) {
    data$S1[t] <- data$S1[t-1] - data$I1[t-1] - data$I2[t-1]
    data$S2[t] <- data$S2[t-1] - data$I2[t-1] - data$I1[t-1]
  }

  # get serial intervals
  omega1 <- rev(gen_distribution(tau_m, sim_mu1, sim_var1, sim_type1, delta)$omega)
  omega2 <- rev(gen_distribution(tau_m, sim_mu2, sim_var2, sim_type2, delta)$omega)

  # TODO Finish model

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {

    # update susceptible
    data$S1[t] <- data$S1[t-1] - data$I1[t-1] - data$I2[t-1]
    data$S2[t] <- data$S2[t-1] - data$I2[t-1] - data$I1[t-1]
    I1_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I1
    data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S1[t] / data$N[t]

    if (t >= (((start_variant + tau_m) * delta) + 1)){
      # update susceptible, substract both normal virus and variant
      data$S2[t] <- data$S2[t-1] - data$I2[t-1] - data$I1[t-1] + (cross * data$I1[t-start+1])
      data$S1[t] <- data$S1[t] - data$I2[t-1] + (cross * data$I2[t-start+1])
      I2_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I2
      data$I2[t] <- data$R_t2[t] * sum(I2_vec * omega2) * data$S2[t] / data$N[t]
    }
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(I1_daily = sum(I1),
              I2_daily = sum(I2),
              R1_val = mean(R_t1), R2_val = mean(R_t2),
              S1_daily = last(S1), S2_daily = last(S2), N = last(N))

  daily_infec$I1_cum <- cumsum(daily_infec$I1_daily)
  daily_infec$I2_cum <- cumsum(daily_infec$I2_daily)
  daily_infec$S1_pct <- daily_infec$S1_daily / daily_infec$N
  daily_infec$S2_pct <- daily_infec$S2_daily / daily_infec$N
  daily_infec$I1_pct <- daily_infec$I1_cum / daily_infec$N
  daily_infec$I2_pct <- daily_infec$I2_cum / daily_infec$N
  daily_infec$infected_day <- daily_infec$I1_daily + daily_infec$I2_daily

  return(daily_infec)
}

## Plot
sii_plot <- function (model){
  plot(x = model$days, y = model$S1_pct, type="l", col = "blue", lwd=2,
    ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$S2_pct, type = "l", col = "red", lwd=2)
  lines(x = model$days, y = model$I1_pct, type = "l", col = "orange", lwd=2)
  lines(x = model$days, y = model$I2_pct, type = "l", col = "green", lwd=2)
  legend("topright",
    legend = c("Susceptible1", 'Susceptible2', "Infected1", "Infected2"),
    col = c("blue", 'red', "orange", 'green'), pch = 16, bty = "n")
}



