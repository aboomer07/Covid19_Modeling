# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/5/2021

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))


si_sim <- function(sim_mu, sim_var, sim_type, delta, days = n_days,
                          tau_m = tau_m, R = R_val, N = population) {
  # Setup environment
  data <- data.frame(matrix(nrow = days * delta, ncol = 6))
  colnames(data) <- c('t', 'days', 'R_t', 'infected', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$R_t <- rep(R, each = (delta * days / length(R)))
  data$N <- N

  # start of epidemic (burn-in)
  # Assume 1 initial ifected and that new ifections happen every 4 days depending on specified R (for bun-in)
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

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta))

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
    I_p <- R_mean * omega * I_vec * data$S[t] / data$N[t] #* (1/delta)
    ## get I(u)
    infec <- sum(I_p)
    data$infected[t] <- infec

    # update total infected
    # data$I[t] <- data$I[t] + infec
  }

  plot(data$infected, type='l')
  abline(v=tau_m*delta, col='red')

  #which(data$infected == max(data$infected), arr.ind=T)

  # return(data)

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
  plot(x = model$days, y = model$S_pct,
       type="l", col = "blue", lwd=2,
       ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I_pct,
        type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"), col = c("blue", "orange"), pch = 16, bty = "n")
}


# test
delta <- 24
n_days <- 300
population <- 6e6
sim_type <- "weibull"
si_model <- si_sim(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5, N = population)
si_plot(si_model)

si_model <- si_sim(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5, N = population)
si_plot(si_model)

si_plot(daily_infec)


test <- nour_sim_data(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5)








#################################### SII MODEL  #######################################
# do the same but add two different infected

sii_sim <- function(sim_mu1, sim_var1, sim_type1, sim_mu2, sim_var2, sim_type2,
                    delta, days = n_days, tau_m = tau_m, R1, R2, N, start_mutation) {

  # setup environment
  # start_mutation <- start_mutation*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 8))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'new_infected1', 'new_infected2', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$N <- N
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$new_infected1[1:(tau_m * delta)] <- 5
  data[1:(start_mutation * delta),]$new_infected2 <- 0
  data$new_infected2[(start_mutation * delta+1):(start_mutation*delta + (tau_m * delta))] <- 5

  data$S[1] <- data$N[1] - data$new_infected1[1] - data$new_infected2[1]
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$new_infected1[t-1] - data$new_infected2[t-1]
  }

  # get serial intervals
  omega1 <- rev(gen_distribution(tau_m, sim_mu1, sim_var1, sim_type1, delta))
  omega2 <- rev(gen_distribution(tau_m, sim_mu2, sim_var2, sim_type2, delta))

# TODO Finish model

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {

    if (t < (start_mutation * delta + tau_m * delta + 1)){
      # update susceptible
      data[which(data$t == t),]$S <- data[which(data$t == t)-1,]$S - data[which(data$t == t)-1,]$new_infected1

      I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$new_infected1
      R1_mean <- mean(data[which(data$t == t),]$R_t1)
      total_infec <- sum(I_vec * omega1)
      # multiply TSI with S/N
      infec <- R1_mean * total_infec * data[which(data$t == t),]$S / data[which(data$t == t),]$N
      data[which(data$t == t),]$new_infected1 <- infec
    }
    else{
      # update susceptible, substract both normal virus and variant
      data[which(data$t == t),]$S <- data[which(data$t == t)-1,]$S - (data[which(data$t == t)-1,]$new_infected1 + data[which(data$t == t)-1,]$new_infected2)

      I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$new_infected1
      R1_mean <- mean(data[which(data$t == t),]$R_t1)
      total_infec <- sum(I_vec * omega1)
      # multiply TSI with S/N
      infec <- R1_mean * total_infec * data[which(data$t == t),]$S / data[which(data$t == t),]$N
      data[which(data$t == t),]$new_infected1 <- infec

      I_vec_var <- data[data$t %in% (t - tau_m*delta):(t - 1),]$new_infected2
      R2_mean <- mean(data[which(data$t == t),]$R_t2)
      total_infec_var <- sum(I_vec_var * omega2)
      # multiply TSI with S/N
      infec_var <- R2_mean * total_infec_var * data[which(data$t == t),]$S / data[which(data$t == t),]$N
      data[which(data$t == t),]$new_infected2 <- infec_var
    }
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day1 = sum(new_infected1),
              infected_day2 = sum(new_infected2),
              R1_val = mean(R_t1),
              R2_val = mean(R_t2),
              S_daily = last(S), N = last(N))

  daily_infec$I1_cum <- cumsum(daily_infec$infected_day1)
  daily_infec$I2_cum <- cumsum(daily_infec$infected_day2)
  daily_infec$S_pct <- daily_infec$S_daily / daily_infec$N
  daily_infec$I1_pct <- daily_infec$I1_cum / daily_infec$N
  daily_infec$I2_pct <- daily_infec$I2_cum / daily_infec$N

  return(daily_infec)
}



## Plot
si_plot <- function (model){
  plot(x = model$days, y = model$S_pct,
       type="l", col = "blue", lwd=2,
       ylim = c(0,1), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I1_pct,
        type = "l", col = "orange", lwd=2)
  lines(x = model$days, y = model$I2_pct,
      type = "l", col = "green", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected1", "Infected2"), col = c("blue", "orange", 'green'), pch = 16, bty = "n")
}


# test
R1 <- 1.5; R2 <- 1.9; sim_mu1 <- 7; sim_var1 <- 2; sim_mu2 <- 7; sim_var2 <- 2
sim_type1 <- "gamma"; sim_type2 <- "gamma"; N <- 60000000; start_mutation <- 75

si_model <- sii_sim(sim_mu1, sim_var1, sim_type1,
                    sim_mu2, sim_var2, sim_type2,
                    delta, n_days, tau_m, R1, R2, N, start_mutation)
si_plot(si_model)