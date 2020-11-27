# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(accelerometry)
library(tidyverse)

###################################
########## Simulating R0 ##########
###################################

### Set up distribution on generation time
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)

GT_obj<-R0::generation.time("gamma", c(6.6, 1.5))

# Define time varying effective reproduction number
# important: here you need to specify manually when you would like the true R0 to change
Ret <- function(date1) {
  if(date1 <= as.Date("2020-03-15")) {
    R_val <- 2
  }
  if (date1 > as.Date("2020-03-15") & date1 <= as.Date("2020-04-15")) {
    R_val <- 1.25
  }
  if (date1 > as.Date("2020-04-15"))  {
    R_val <- 0.9
  }
  return(R_val)
}

# use the R0 to simulate the amount of infected
routbreak <- function(n=100, Ret, GT_obj, initial_cases = 10) {
  # Set up time series of incident cases
  y <- rep(0, n + length(GT_pmf))
  y[seq_len(length(initial_cases))] <- initial_cases
  # Outbreak starts on 2020-02-15
  dates <- as.Date("2020-02-15") + 0:(n-1)
  # Extract serial interval PMF, ignore support at 0.
  GT_pmf <- GT_obj$GT[-1]

  # Loop over all time points
  for (i in 1:n) {
    date <- dates[i]
    y[i + 1:length(GT_pmf)] <- y[i] * (Ret(date) * GT_pmf) + y[i + 1:length(GT_pmf)]
  }

  # Data frame with the result. Assume we start on 15th of Feb
  res <- data.frame(Date=dates, y=y[1:n])

  #Done
  return(res)
}

# Generate an outbreak (no stochasticity, just the difference equation)
out <- routbreak(n=100, Ret=Ret, GT_obj=GT_obj)

# Data frame with the true values
True_val <- NULL
for (i in 1:length(out$Date)) {
  date <- out$Date[i]
  True_val <- c(True_val, Ret(date))
}
n1<-100
names(True_val) <- as.Date("2020-02-15") + 0:(n1-1)
Data_R_sim <- data.frame(dates=as.Date("2020-02-15") + 0:(n1-1), R_sim=True_val)

###################################
########## Estimating R0 ##########
###################################

# Austrian method
est_r0 <- estimate_R(incid = out$y, method = 'parametric_si', config = make_config(list(
                  mean_si = 4.46 , std_si = 2.63)))[[1]]

est_r0 <- est_r0[c('Mean(R)', 'Std(R)')]
est_r0$dates <- as.Date("2020-02-15") + 7:(100-1)
est_r0 <- est_r0 %>%
  rename(R_austr = "Mean(R)")

# put it together
eval <- left_join(Data_R_sim, est_r0[c('R_austr', 'dates')], by = 'dates')

eval_long <- reshape2::melt(eval, id.vars = 'dates')

ggplot(NULL) +
  geom_line(data = eval_long, aes(x = dates, y=value, color = variable))

}

# German method

# R0 with serial time of 4 days
r0_ger <- rep(NA, nrow(out))
for (t in 8:nrow(out)) {
 r0_ger[t] <- sum(out$y[t-0:3]) / sum(out$y[t-4:7])
}
r0_ger <- data.frame(r0_ger)
# bind to evaluation dataframe
r0_ger <- cbind(Data_R_sim, r0_ger)

eval_long <- reshape2::melt(r0_ger, id.vars = 'dates')
ggplot(NULL) +
  geom_line(data = eval_long, aes(x = dates, y=value, color = variable)


params <- list(list(4.46, 2.63), list(5.2, 5.1), list(4, 4))
names(params) <- c('Austria', 'USA', 'France')

full_eval <- list()

for (i in 1:length(params)) {
  mean <- params[[i]][[1]]
  std <- params[[i]][[2]]
  est_r0 <- estimate_R(incid=out$y, method='parametric_si', 
    config=make_config(list(mean_si=mean, std_si=std)))[[1]]
  est_r0 <- est_r0[c('Mean(R)', 'Std(R)')]
  names(est_r0) <- c("R_Mean", 'R_Std')
  est_r0$dates <- as.Date("2020-02-15") + 7:(100-1)
  eval <- left_join(Data_R_sim, est_r0[c('R_Mean', 'R_Std','dates')], 
    by = 'dates')
  eval['lower'] <- eval['R_Mean'] - (1.96 * eval['R_Std'])
  eval['upper'] <- eval['R_Mean'] + (1.96 * eval['R_Std'])
  eval$Country <- rep(names(params)[[i]], dim(eval)[[1]])
  full_eval[[i]] <- eval
}

df <- do.call(rbind, full_eval)

simplot <- ggplot(data=df) +  
  geom_line(aes(x=dates, y=R_Mean, color='red')) +  
  geom_line(aes(x=dates, y=R_sim, col='blue')) +  
  geom_ribbon(aes(x=dates, ymin=lower, ymax=upper, fill='red'), alpha=0.3) +
  facet_wrap(~Country, nrow=1) +  ylab('Value of R') +  
  scale_fill_identity(name = 'Confidence Interval', 
    guide = 'legend',labels = c('CI')) +  
  scale_colour_manual(name = 'R0 Type', 
    values =c('red'='red','blue'='blue'), labels = c('Sim R0', 'Est R0')))

# Italian method
#R0 with rolling window 5 
#Notice that window = 5 is for the MA calculation, 
#while window_interval should 
#represent the generation time, which we take to be 4

out$tot <- function(t) out$y 
r0_itafun <- function(ts, window_interval) {
  data1 <- accelerometry::movingaves(ts, window=5)
  names(data1) <- names(ts)[3:(length(ts)-2)]
  res <- sapply( (1+window_interval):length(data1), function(t) {
    data1[t]/data1[t-window_interval]
  })
  return(res)
}

r0_ita <- r0_itafun(out$y, 4)
r0_ita <- data.frame(r0_ita) #Length of 92 since res only starts at 5th observation and MA window is +-2 

out$Date = as.character(out$Date)

for (t in 7:(length(out$Date)-2)){
r0_ita$dates[t-6] <- out$Date[t]
}

r0_ita$dates = as.Date(r0_ita$dates)
evalita <- left_join(Data_R_sim, r0_ita[c('r0_ita', 'dates')], 
  left.by = 'dates')
names(evalita) <- c("dates", "R_sim", "R_Mean")

names(r0_ger) <- c("dates", "R_sim", "R_Mean")
r0_ger$Country <- rep("Germany", dim(r0_ger)[1])
evalita$Country <- rep("Italy", dim(r0_ger)[1])

ggplot(NULL) +
  geom_line(data = evalita, aes(x = dates, y=r0_ita, color = r0_ita))

non_ci <- rbind(df[c("dates", "R_sim", "R_Mean", "Country")], evalita)
non_ci <- rbind(non_ci, r0_ger)




