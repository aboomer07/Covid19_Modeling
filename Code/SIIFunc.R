# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 4/5/2021

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))


si_sim <- function(sim_mu, sim_var, sim_type, delta, days = n_days,
                          tau_m = tau_m, R = R_val, N = population) {

  data <- data.frame(matrix(nrow = days * delta, ncol = 6))
  colnames(data) <- c('t', 'days', 'R_t', 'infected', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$infected[1:(tau_m * delta)] <- round(seq(10, 10, length.out = tau_m * delta))/delta
  data$R_t <- rep(R, each = (delta * days / length(R)))
  data$N <- N
  data$S[1] <- data$N[1] - data$infected[1]
  # remember that infected are new daily infections, not total
  for (t in 2:(tau_m * delta)){
    data$S[t] <- data$S[t-1] - data$infected[t-1]
  }

  omega <- rev(gen_distribution(tau_m, sim_mu, sim_var, sim_type, delta))

  # simulate outbreak
  start <- ((tau_m) * delta) + 1

  for (t in start:dim(data)[1]) {

    # update susceptible
    data[which(data$t == t),]$S <- data[which(data$t == t)-1,]$S - data[which(data$t == t)-1,]$infected

    I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$infected
    R_mean <- mean(data[which(data$t == t),]$R_t)
    total_infec <- sum(I_vec * omega)
    # multiply TSI with S/N
    infec <- R_mean * total_infec * data[which(data$t == t),]$S / data[which(data$t == t),]$N
    data[which(data$t == t),]$infected <- infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day = round(sum(infected)), R_val = mean(R_t),
              S_daily = mean(S), N = mean(N)) %>%
    mutate(I_daily = round(N - S_daily))

  return(daily_infec)
}



## Plot
si_plot <- function (model){
  plot(x = model$days, y = model$S_daily,
       type="l", col = "blue", lwd=2,
       ylim = c(0,population), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I_daily,
        type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"), col = c("blue", "orange"), pch = 16, bty = "n")
}


# test
delta <- 24
n_days <- 300
population <- 1e6
si_model <- si_sim(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5, N = population)
si_plot(si_model)


test <- nour_sim_data(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5)




#################################### SII MODEL  #######################################
# do the same but add two different infected
R1 <- 1.2; R2 <- 1.7; sim_mu1 <- 7; sim_var1 <- 2; sim_mu2 <- 3; sim_var2 <- 1
sim_type <- "gamma"; sim_type2 <- "gamma"; N <- 60000000; start_mutation <- 20



sii_sim <- function(sim_mu1, sim_var1, sim_type1, sim_mu2, sim_var2, sim_type2,
                    delta, days = n_days, tau_m = tau_m, R1, R2, N, start_mutation) {

  # setup environment
  start_mutation <- start_mutation*delta
  data <- data.frame(matrix(nrow = days * delta, ncol = 8))
  colnames(data) <- c('t', 'days', 'R_t1', 'R_t2', 'new_infected1', 'new_infected2', 'S', 'N')
  data$t <- 1:dim(data)[1]
  data$days <- rep(1:days, each = delta)
  data$N <- N
  data$R_t1 <- rep(R1, each = (delta * days / length(R1)))
  data$R_t2 <- rep(R2, each = (delta * days / length(R2)))

  # start of the two pandemics (1 assumed to start in beginning, 2 is variable)
  data$new_infected1[1:(tau_m * delta)] <- round(seq(10, 100, length.out = tau_m * delta))
  data[1:(start_mutation - delta),]$new_infected2 <- 0
  data$new_infected2[(start_mutation - delta+1):(start_mutation+(tau_m * delta))] <- round(seq(10, 100, length.out = tau_m * delta))

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

    # update susceptible
    data[which(data$t == t),]$S <- data[which(data$t == t)-1,]$S - data[which(data$t == t)-1,]$infected

    I_vec <- data[data$t %in% (t - tau_m*delta):(t - 1),]$infected
    R_mean <- mean(data[which(data$t == t),]$R_t)
    total_infec <- sum(I_vec * omega)
    # multiply TSI with S/N
    infec <- R_mean * total_infec * data[which(data$t == t),]$S / data[which(data$t == t),]$N
    data[which(data$t == t),]$infected <- infec
  }

  # aggregate to daily (avg = /delta)
  daily_infec <- data %>%
    group_by(days) %>%
    summarise(infected_day = round(mean(infected)), R_val = mean(R_t),
              S_daily = mean(S), N = mean(N)) %>%
    mutate(I_daily = round(N - S_daily))

  return(daily_infec)
}



## Plot
si_plot <- function (model){
  plot(x = model$days, y = model$S_daily,
       type="l", col = "blue", lwd=2,
       ylim = c(0,population), ylab = "Susceptible and Infected Population", xlab = "Days")
  lines(x = model$days, y = model$I_daily,
        type = "l", col = "orange", lwd=2)
  legend("topright", legend = c("Susceptible", "Infected"), col = c("blue", "orange"), pch = 16, bty = "n")
}


# test
si_model <- si_sim(sim_mu, sim_var, sim_type, delta, n_days, tau_m, R = 1.5, N = population)
si_plot(si_model)
