# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(accelerometry)
library(tidyverse)

# Import simulated data (OutbreakSim)
imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

sim_const <- read.csv2(paste0(imppath, 'SimulatedSIR.csv'))[-1, ]
sim_dyn <- read.csv2(paste0(imppath, 'SimulatedSIR_dyn.csv'))[-1, ]
# sim_dyn <- sim_dyn[20:dim(sim_dyn)[1], names(sim_dyn)]

###################################
########## Estimating R0 ##########
###################################

# Austrian, US and French method

params <- list(list(4.46, 2.63), list(5.2, 5.1), list(4, 4))
names(params) <- c('Austria', 'USA', 'France')

full_eval <- list()

for (i in 1:length(params)) {
  mean <- params[[i]][[1]]
  std <- params[[i]][[2]]
  est_r0 <- estimate_R(incid=as.numeric(sim_dyn$Delta), method='parametric_si', config=make_config(list(mean_si=mean, std_si=std)))[[1]]
  est_r0 <- est_r0[c('Mean(R)', 'Std(R)', 't_end')]
  names(est_r0) <- c("R_Mean", 'R_Std', 'Days')
  eval <- left_join(sim_dyn, est_r0, by = 'Days')
  eval['lower'] <- eval['R_Mean'] - (1.96 * eval['R_Std'])
  eval['upper'] <- eval['R_Mean'] + (1.96 * eval['R_Std'])
  eval$Country <- rep(names(params)[[i]], dim(eval)[[1]])
  full_eval[[i]] <- eval
}

df <- do.call(rbind, full_eval)

# German method

# R0 with serial time of 4 days
int_list <- list()
intervals <- c(4)
for (inter in c(1)) {
  interval <- intervals[inter]
  r0_ger <- rep(NA, nrow(sim_dyn))
  for (t in (3 + interval):nrow(sim_dyn)) {
   r0_ger[t] <- sum(as.numeric(sim_dyn$Delta)[t-0:3]) / sum(as.numeric(sim_dyn$Delta)[t-interval:(3 + interval)])
  }
  r0_ger <- data.frame(r0_ger)
  # bind to evaluation dataframe
  r0_ger <- cbind(sim_dyn, r0_ger)

  long <- reshape2::melt(r0_ger, id.vars = 'Days')
  long['Interval'] <- rep(interval, dim(long)[1])

  int_list[[inter]] <- long
}

#ger <- do.call(rbind, int_list)
#ggplot(data=ger) +
#  geom_line(aes(x = Days, y=value, color = variable)) +
#  facet_wrap(~Interval, nrow=2, ncol=2)

# Italian method
#R0 with rolling window 5 
#Notice that window = 5 is for the MA calculation, 
#while window_interval should 
#represent the generation time, which we take to be 4

r0_itafun <- function(ts, omega, tau) {
  data1 <- accelerometry::movingaves(ts, window=tau)
  names(data1) <- names(ts)[3:(length(ts)-2)]
  res <- sapply((1+omega+tau):length(data1), function(t) {
    data1[t]/data1[t-omega]
  })
  return(res)
}

intervals <- c(4)
italist <- list()
for (inter in c(1)) {
  interval <- intervals[inter]
  r0_ita <- r0_itafun(as.numeric(sim_dyn$Delta), interval, 1)
  r0_ita <- data.frame(r0_ita) #Length of 92 since res only starts at 5th observation and MA window is +-2
}
r0_ita$Days = 1:dim(r0_ita)[1]

evalita <- left_join(sim_dyn, r0_ita[c('r0_ita', 'Days')], by = 'Days')

#ita <- do.call(rbind, italist)
#ggplot(data = ita) +
#  geom_line(aes(x = dates, y = value, color = variable)) +
#  facet_wrap(~Interval, nrow = 2, ncol = 2) +
#  ggsave("ItalyRobustness.png")


# names(r0_ger) <- c("dates", "R_sim", "R_Mean")
colnames(r0_ger)[colnames(r0_ger) == 'r0_ger'] <- 'R_Mean'
r0_ger['R_Std'] <- rep(0, dim(r0_ger)[1])
r0_ger['lower'] <- rep(0, dim(r0_ger)[1])
r0_ger['upper'] <- rep(0, dim(r0_ger)[1])
r0_ger$Country <- rep("Germany", dim(r0_ger)[1])

colnames(evalita)[colnames(evalita) == 'r0_ita'] <- 'R_Mean'
evalita['R_Std'] <- rep(0, dim(evalita)[1])
evalita['lower'] <- rep(0, dim(evalita)[1])
evalita['upper'] <- rep(0, dim(evalita)[1])
evalita$Country <- rep("Italy", dim(r0_ger)[1])

non_ci <- rbind(df, evalita)
non_ci <- rbind(non_ci, r0_ger)

simplot <- ggplot(data=non_ci) +  
  geom_line(aes(x=Days, y=R_Mean, color='red')) +
  geom_line(aes(x=Days, y=as.numeric(Rt), col='blue')) +
  geom_ribbon(aes(x=Days, ymin=lower, ymax=upper, fill='red'), alpha=0.3) +
  facet_wrap(~Country, nrow=2, ncol=3) +  ylab('Value of R') +  
  scale_fill_identity(name = 'Confidence Interval', 
    guide = 'legend',labels = c('CI')) +  
  scale_colour_manual(name = 'Rt Type', 
    values =c('red'='red','blue'='blue'), labels = c('Sim Rt', 'Est Rt')) +
  ggsave(paste0(outpath, "All_Countries_withCI_Rt.png"))





###################################
########## Legacy Code ############
###################################

#### Set up distribution on generation time
#GT_obj<-R0::generation.time("gamma", c(6.6, 1.5))
#
## Define time varying effective reproduction number
## important: here you need to specify manually when you would like the true R0 to change
#Ret <- function(date1) {
#  if(date1 <= as.Date("2020-03-15")) {
#    R_val <- 1.5
#  }
#  if (date1 > as.Date("2020-03-15") & date1 <= as.Date("2020-04-15")) {
#    R_val <- 1.25
#  }
#  if (date1 > as.Date("2020-04-15"))  {
#    R_val <- 1.4
#  }
#  return(R_val)
#}
#
## use the R0 to simulate the amount of infected
#routbreak <- function(n=100, Ret, GT_obj, initial_cases = 10) {
#  # Set up time series of incident cases
#  y <- rep(0, n + length(GT_pmf))
#  y[seq_len(length(initial_cases))] <- initial_cases
#  # Outbreak starts on 2020-02-15
#  dates <- as.Date("2020-02-15") + 0:(n-1)
#  # Extract serial interval PMF, ignore support at 0.
#  GT_pmf <- GT_obj$GT[-1]
#
#  # Loop over all time points
#  for (i in 1:n) {
#    date <- dates[i]
#    y[i + 1:length(GT_pmf)] <- y[i] * (Ret(date) * GT_pmf) + y[i + 1:length(GT_pmf)]
#  }
#
#  # Data frame with the result. Assume we start on 15th of Feb
#  res <- data.frame(Date=dates, y=y[1:n])
#
#  #Done
#  return(res)
#}
#
## Generate an outbreak (no stochasticity, just the difference equation)
#out <- routbreak(n=100, Ret=Ret, GT_obj=GT_obj)
#
## Data frame with the true values
#True_val <- NULL
#for (i in 1:length(out$Date)) {
#  date <- out$Date[i]
#  True_val <- c(True_val, Ret(date))
#}
#n1<-100
#names(True_val) <- as.Date("2020-02-15") + 0:(n1-1)
#Data_R_sim <- data.frame(dates=as.Date("2020-02-15") + 0:(n1-1), R_sim=True_val)
