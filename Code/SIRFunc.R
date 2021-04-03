# Objective: script that contains function to simulate cases from SIR model that can then be used in master file
# Created by: jacobpichelmann
# Created on: 03.04.21

# import helper functions
source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))

#import libraries
library(deSolve)


sir_sim <- function(N, I, R = 0, Rt, gamma, length){
  ## Create an SIR function
  sir <- function(times, state, parameters) {

    with(as.list(c(state, parameters)), {

      dS <- -beta * S * I/N
      dI <-  beta * S * I/N - gamma * I
      dR <-                 gamma * I

      return(list(c(dS, dI, dR)))
    })
}

  # initialize population compartments
  init <- c(S = N - I - R, I = I, R = R)
  # set parameters
  parameters <- c(beta = Rt * gamma, gamma = gamma)
  # set time span of outbreak
  times <- seq(1, length, by = 1)

  # solve SIR model using ode (General Solver for Ordinary Differential Equations)
  out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
  out$I <- round(out$I) # round infected

  # get in correct format so we can use R_t function
  colnames(out)[3] <- 'infected_day'
  colnames(out)[1] <- 'days'
  out$R_val <- Rt

  return(out)
}

# showcase how it works
n_days <- 100
cases <- sir_sim(N = 1e6, I = 10, Rt = 1.4, gamma = 1/18, length = n_days)

samps <- samp_pois(1.4, 20, num_people = 100, sim_mu = 6.6, sim_sig = 1.1, sim_type = 'gamma',
                   delta = 24)

Rt_est_nonpara(cases, samps$daily, bw = 'iqr')