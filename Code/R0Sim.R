# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)

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
est_r0_austr <- estimate_R(incid = out$y, method = 'parametric_si', config = make_config(list(
                  mean_si = 4.46 , std_si = 2.63)))[[1]]

est_r0_austr <- est_r0_austr[c('Mean(R)', 'Std(R)')]
est_r0_austr$dates <- as.Date("2020-02-15") + 7:(100-1)
est_r0_austr <- est_r0_austr %>%
  rename(R_austr = "Mean(R)")

# put it together
eval <- left_join(Data_R_sim, est_r0_austr[c('R_austr', 'dates')], by = 'dates')

eval_long <- reshape2::melt(eval, id.vars = 'dates')

ggplot(NULL) +
  geom_line(data = eval_long, aes(x = dates, y=value, color = variable))



# German method

# R0 with serial time of 4 days
r0_ger <- rep(NA, nrow(out))
for (t in 8:nrow(out)) {
 r0_ger[t] <- sum(out$y[t-0:3]) / sum(out$y[t-4:7])
}
r0_ger <- data.frame(r0_ger)


# bind to evaluation dataframe
eval <- eval %>% bind_cols(r0_ger$r0_ger) %>% rename(r0_ger = ...4)

eval_long <- reshape2::melt(eval, id.vars = 'dates')
ggplot(NULL) +
  geom_line(data = eval_long, aes(x = dates, y=value, color = variable))