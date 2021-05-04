# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/5/2021

source(paste0(getwd(), "/Code/Params.R"))
source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))

si_sim <- function(params) {
  # unpack variables
  days <- params$n_days
  delta <- params$delta
  R <- params$R_val
  tau_m <- params$tau_m
  sim_mu <- params$sim_mu
  sim_var <- params$sim_var
  sim_type <- params$sim_type
  init_infec <- params$init_infec

  start <- ((tau_m) * delta) + 1

  # Setup environment
  data <- data.frame(matrix(nrow = days * delta, ncol = 5))
  colnames(data) <- c('t', 'days', 'R_t', 'infected', 'S')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$R_t <- rep(R, each = (delta * days / length(R)))

  # start of epidemic (burn-in)
  # Assume 1 initial ifected and that new ifections happen every 4 days depending on specified R (for burn-in)
  data$infected[1:(tau_m * delta)] <- init_infec
  data$S[1] <- params$pop - data$infected[1]
  data$I <- data$infected

  # Simulate initial infected spreading disease during burn-in period
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$infected[t-1]
    # data$infected[t] <- sum(data$infected[1:t-1]) * R / (delta*2)
    # data$infected[t] <- init_infec
    data$I[t] <- data$I[t-1] + data$infected[t]
  }

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta)$omega)

  # simulate outbreak
  for (t in start:dim(data)[1]) {

    # update susceptible population
    data$S[t] <- data$S[t-1] - data$infected[t-1]

    # compute I(u):
    ## compute I(p) first:
    ## get omega(n*delta) I(t-n*delta) => total_infec (= All infected for a TSI from 1 until tau_m at all deltas)
    I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$infected
    I_p <- data$R_t[t] * omega * I_vec * data$S[t] / params$pop

    ## get I(u)
    infec <- sum(I_p)
    data$infected[t] <- infec

    # update total infected
    data$I[t] <- data$I[t-1] + infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day=sum(infected), R_val=mean(R_t), S_daily=last(S))

  daily_infec$I_cum <- cumsum(daily_infec$infected_day)
  daily_infec$S_pct <- daily_infec$S_daily / params$pop
  daily_infec$I_pct <- daily_infec$I_cum / params$pop

  return(daily_infec)
}

## Plot
si_plot <- function (model, params){
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2, 
    ylim = c(0, 1), 
    ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I_pct, type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"), 
    col = c("blue", "orange"), pch = 16, bty = "n")
}

si_plot_detail <- function (model){
  layout(matrix(c(rep(1,6), rep(2,3), rep(3,3)), nrow=4, ncol=3, byrow = T))
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2,
    ylim = c(0, 1), 
    ylab = "Susceptible and Infected Population", xlab="")
  lines(x = model$days, y = model$I_pct, type = "l", col = "yellow", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"),
    col = c("blue", "orange"), pch = 20, bty = "n")
  plot(x = model$days, y = model$infected_day, type="l", lwd=2, col = "green",
       ylim=c(0, 1), ylab = "Infected", xlab = "")
  legend("topright", legend = "New Infections", col="green", pch = 16, bty="n")
  plot(x = model$days, y = model$R_val, type="l", lwd=2, col = "red",
       ylim=c(0, 3), ylab = "R(t)", xlab = "Days")
  legend("topright", legend = "R(t)", col="red", pch = 16, bty="n")
  layout(matrix(1, nrow=1, ncol=1, byrow = T))
}


# test

############################## SII MODEL #######################################
# do the same but add two different infected


sii_sim <- function(params) {

  days <- params$n_days
  delta <- params$delta
  R1 <- params$R_val
  tau_m <- params$tau_m
  sim_mu1 <- params$sim_mu
  sim_var1 <- params$sim_var
  sim_type1 <- params$sim_type
  sim_mu2 <- params$sim_mu_variant
  sim_var2 <- params$sim_var_variant
  R2 <- params$R_val_variant
  start_variant <- params$start_variant
  sim_type2 <- params$sim_type_variant

  start <- ((tau_m) * delta) + 1

  # setup environment
  # start_variant <- start_variant*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 7))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'I1', 'I2', 'S')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$I1[1:(tau_m * delta)] <- 5
  data[1:(start_variant * delta),]$I2 <- 0
  data$I2[(start_variant * delta+1):((start_variant + tau_m) * delta)] <- 5
  data$infected_day <- data$I1
  data$S[1] <- params$pop
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]
  }

  # get serial intervals
  omega1 <- rev(gen_distribution(tau_m, sim_mu1, sim_var1, sim_type1, delta)$omega)
  omega2 <- rev(gen_distribution(tau_m, sim_mu2, sim_var2, sim_type2, delta)$omega)

  # simulate outbreak
  for (t in start:dim(data)[1]) {
    if (t < ((start_variant + tau_m) * delta + 1)){
      data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]

      I1_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I1
      data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S[t] / params$pop
      data$infected_day[t] <- data$I1[t]
    }

    else{
      data$S[t] <- data$S[t-1] - data$I1[t-1] - data$I2[t-1]

      I1_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I1
      data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S[t] / params$pop
      I2_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I2
      data$I2[t] <- data$R_t2[t] * sum(I2_vec * omega2) * data$S[t] / params$pop
      data$infected_day[t] <- data$I1[t] + data$I2[t]
    }
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(I1_daily = sum(I1), I2_daily = sum(I2),
              infected_day = sum(infected_day),
              R1_val = mean(R_t1), R2_val = mean(R_t2),
              S_daily = last(S))

  daily_infec$I1_cum <- cumsum(daily_infec$I1_daily)
  daily_infec$I2_cum <- cumsum(daily_infec$I2_daily)
  daily_infec$I_cum <- daily_infec$I1_cum + daily_infec$I2_cum
  daily_infec$S_pct <- daily_infec$S_daily / params$pop
  daily_infec$I1_pct <- daily_infec$I1_cum / params$pop
  daily_infec$I2_pct <- daily_infec$I2_cum / params$pop
  daily_infec$I_pct <- daily_infec$I_cum / params$pop

  return(daily_infec)
}

sii_plot <- function (model, Rt){
  layout_mat <- matrix(c(c(1, 3), c(2, 3)), nrow=2, ncol=2, byrow = T)
  layout(layout_mat)
  # layout(mat = matrix(c(2, 3, 0, 1), nrow = 2, ncol = 2),
  #      heights = c(2, 2), widths = c(2, 2))
  plot(x = model$days, y = model$S_pct, type="l", col = "blue", lwd=2,
    ylim = c(0, 1), 
    ylab = "Susceptible", xlab = "Days")
  lines(x = model$days, y = model$I1_pct, col = "yellow", lwd=2, lty='dotted')
  # lines(x = model$days, y = model$I2_pct, col = "green", lwd=2)
  lines(x = model$days, y = model$I_pct, col='yellow', lwd=2)
  legend("left",
    legend = c("Susceptible", 'Infections1', 'Total Infected'),
    # legend = c("Susceptible", "Infected1", "Infected2"),
    col = c("blue", 'yellow', "yellow"), pch = 16, bty = "n",
    lty=c('solid', 'dotted', 'solid'))

  plot(model$days, y=model$I1_daily, type='l', col="brown", lwd=2,
    ylab='Daily Infected', xlab='Days', ylim=c(0, max(max(model$I1_daily), max(model$I2_daily))))
  lines(model$days, y=model$I2_daily, col='green', lwd=2)
  legend('topleft',
    legend=c('Daily Cases', 'Daily Cases Variant'),
    col=c('brown', 'green'), pch=16, bty='n')

  plot(x = Rt$Date, y = Rt$Rt1, type="l", col="red", lwd=2,
    ylim = c(0, 3), ylab = "Rt", xlab = "Days")
  lines(x = Rt$Date, y = Rt$Rt2, col="orange", lwd=2)
  lines(x = Rt$Date, y = Rt$Est_Rt1, col="black", lwd=2, lty='dotted')
  lines(x = Rt$Date, y = Rt$Est_Rt2, col='black', lwd=2, lty='dotted')
  lines(x = Rt$Date, y = Rt$Est_Rt, col="black", lwd=2)
  # lines(x = Rt$Date, y = Rt$Rt_Avg, col='green', lwd=2)
  legend("bottomleft",
    legend = c("True Rt1", "True Rt2", 'Est Rt1', 'Est Rt2', "Est Rt Overall"),
    col = c("red", "orange", 'black', 'black', 'black'), pch = 16, bty = "n",
    lty=c('solid', 'solid', 'dotted','dotted', 'solid', 'solid'))
 }


############################## SSII MODEL ######################################
# do the same but add two different infected + different susceptibility groups

ssii_sim <- function(params) {

  days <- params$n_days
  delta <- params$delta
  R1 <- params$R_val
  tau_m <- params$tau_m
  sim_mu1 <- params$sim_mu
  sim_var1 <- params$sim_var
  sim_type1 <- params$sim_type
  sim_mu2 <- params$sim_mu_variant
  sim_var2 <- params$sim_var_variant
  R2 <- params$R_val_variant
  start_variant <- params$start_variant + tau_m
  sim_type2 <- params$sim_type_variant
  cross <- params$sii_cross

  # setup environment
  # start_variant <- start_variant*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 8))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'I1', 'I2', 'S1', 'S2')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$I1[1:(tau_m * delta)] <- 5
  data[1:(start_variant * delta),]$I2 <- 0
  data$I2[((start_variant * delta)+1):((start_variant + tau_m) * delta)] <- 5

  data$S1[1] <- params$pop
  data$S2[1] <- params$pop
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
    data$I1[t] <- data$R_t1[t] * sum(I1_vec * omega1) * data$S1[t] / params$pop

    if (t >= (((start_variant + tau_m) * delta) + 1)){
      # update susceptible, substract both normal virus and variant
      data$S2[t] <- data$S2[t-1] - data$I2[t-1] - data$I1[t-1] + (cross * data$I1[t-start+1])
      data$S1[t] <- data$S1[t] - data$I2[t-1] + (cross * data$I2[t-start+1])
      I2_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$I2
      data$I2[t] <- data$R_t2[t] * sum(I2_vec * omega2) * data$S2[t] / params$pop
    }
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(I1_daily = sum(I1),
              I2_daily = sum(I2),
              R1_val = mean(R_t1), R2_val = mean(R_t2),
              S1_daily = last(S1), S2_daily = last(S2))

  daily_infec$I1_cum <- cumsum(daily_infec$I1_daily)
  daily_infec$I2_cum <- cumsum(daily_infec$I2_daily)
  daily_infec$I_cum <- daily_infec$I1_cum + daily_infec$I2_cum
  daily_infec$infected_day <- daily_infec$I1_daily + daily_infec$I2_daily
  daily_infec$S1_pct <- daily_infec$S1_daily / params$pop
  daily_infec$S2_pct <- daily_infec$S2_daily / params$pop
  daily_infec$I1_pct <- daily_infec$I1_cum / params$pop
  daily_infec$I2_pct <- daily_infec$I2_cum / params$pop
  daily_infec$I_pct <- daily_infec$I_cum / params$pop

  return(daily_infec)
}

## Plot
ssii_plot <- function (model, Rt){
  layout(nrow=2, ncol=1)
  plot(x = model$days, y = model$S1_pct, type="l", col = "black", lwd=2,
    ylim = c(0, 1), 
    ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$S2_pct, col = "orange", lwd=2)
  lines(x = model$days, y = model$I1_pct, col = "blue", lwd=2)
  lines(x = model$days, y = model$I2_pct, col = "green", lwd=2)
  lines(x = model$days, y = model$I_pct, col='red', lwd=2)
  legend("topright",
    legend = c("Susceptible1", "Susceptible2", "Infected1", "Infected2", "Total Infected"),
    col = c("black", "orange", 'blue', 'green', 'red'), pch = 16, bty = "n")

  plot(x = Rt$Date, y = Rt$R_t1, type="l", col = "black", lwd=2,
    ylim = c(0, 3), ylab = "Rt", xlab = "Days")
  lines(x = Rt$Date, y = Rt$R_t2, col = "orange", lwd=2)
  lines(x = Rt$Date, y = Rt$Est_Rt, col = "blue", lwd=2)
  lines(x = Rt$Date, y = Rt$Est_Rt1, col = "green", lwd=2)
  lines(x = Rt$Date, y = Rt$Est_Rt2, col='red', lwd=2)
  legend("topright",
    legend = c("True Rt1", "True Rt2", "Est Rt Overall", 'Est Rt1', 'Est Rt2'),
    col = c("black", "orange", 'blue', 'green', 'red'), pch = 16, bty = "n")
}



