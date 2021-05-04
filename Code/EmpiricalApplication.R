# Applying our estimators to French data
# Created by: jacobpichelmann
# Created on: 15.04.21

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/SIIFunc.R"))
source(paste0(getwd(), "/Code/Params.R"))

library(urca)

# first application: late summer 2020 (august, sept)
df <- read.csv2(paste0(imppath, 'open_stats_coronavirus.csv')) %>% filter(nom == 'france') %>% select(date, cas) %>%
  mutate(new_cases = cas-lag(cas),
         date = as.Date(date))
colnames(df) <- c('date', 'cum_cases', 'infected_day') # need correct col names for function

# 2020-06-02 must be a data issue as cumulative cases drop by roughly 700

# set global params
params['study_len'] <- 10
params['tau_m'] <- params$study_len
params['sim_mu'] <- 4.8
params['sim_var'] <- 2.8
params['delta'] <- 1

Rt_est_applied <- function(df, params, variant = F) {
  start <- params[['study_len']]
  type <- params[['sim_type']]
  start_variant <- params$start_variant
  tau_m <- params$tau_m

  mean <- params[['sim_mu']]
  var <- params[['sim_var']]

  if (variant) {
    data <- data.frame(matrix(nrow = nrow(df), ncol = 4))
    names(data) <- c('Date', 'Est_Rt', 'Est_Rt1', 'Est_Rt2')
    data$Est_Rt1 <- 1
    data$Est_Rt2 <- 1
    data$Date <- seq(1, nrow(df))
    data[data$Date <= start, ]$Est_Rt1 <- NA
    data[data$Date <= (start_variant + tau_m), ]$Est_Rt2 <- NA

    for (i in start:nrow(data)) {
      if (data[i,]$Date > start) {
        t <- data[i,]$Date

        omega <- rev(gen_distribution(t - 1, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% 1:(t - 1),]$infected_day
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
        I1 <- df[which(df$days == t),]$I1_daily
        I2 <- df[which(df$days == t),]$I2_daily
        I1_vec <- df[df$days %in% 1:(t - 1),]$I1_daily
        I2_vec <- df[df$days %in% 1:(t - 1),]$I2_daily
        data$Est_Rt1[i] <- (I1) / (sum(I1_vec * omega))
        if (t > (start_variant + tau_m)) {
            data$Est_Rt2[i] <- (I2) / (sum(I2_vec * omega))
          }
        }
      }
    }
  else {
    data <- data.frame(matrix(nrow = nrow(df), ncol = 2))
    names(data) <- c('Date', 'Est_Rt')
    data$Date <- seq(1, nrow(df))

    for (i in start:nrow(data)) {
      if (data[i,]$Date > start) {
        t <- data[i,]$Date

        omega <- rev(gen_distribution(t - 1, mean, var, type, 1)$omega)
        I <- df[which(df$days == t),]$infected_day
        I_window <- df[df$days %in% 1:(t - 1),]$infected_day
        data[i,]$Est_Rt <- (I) / (sum(I_window * omega))
      }
    }
  }
  return(data)
}

apply_rt_est <- function(data, start_date, end_date, params, plotname){
  dat <- data %>% filter(date > as.Date(start_date) & date < as.Date(end_date))

  # compute weekly rolling average to deal with weekend bias
  dat <- dat %>% mutate(
    infected_day = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit()
  dat$days <- 1:nrow(dat)

  Rt <- Rt_est_applied(dat, params, variant = F)
  Rt$date <- dat$date
  Rt %>% select(date, Est_Rt) %>% na.omit() %>%
  ggplot() + geom_line(aes(x = date, y = Est_Rt)) +
  theme_minimal() +
  labs(x = '', y = 'Estimated Rt') +
  ggsave(paste0(outpath, plotname, '.png'), width = 10, height = 5)

  return(Rt)
}

forecast_rt <- function(data, Rt_df, start_date, end_date, params, Rt_window, days_ahead, plotname){
  dat <- data %>% filter(date > as.Date(start_date) & date < as.Date(end_date))

  # compute weekly rolling average to deal with weekend bias
  dat <- dat %>% mutate(
    infected_day = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit()
  dat$days <- 1:nrow(dat)

  # take average and forecast next month
  params['R_val'] <- mean(tail(Rt_df, Rt_window)$Est_Rt, na.rm = T)
  params['pop'] <- 67000000 # roughly french population
  params['init_infec'] <- list(tail(dat, params$tau_m)$infected_day)
  params['n_days'] <- days_ahead + params$tau_m

  sim <- si_sim(params)

  actual_cases <- data %>%
    mutate(actual_cases = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit() %>%
    filter(date > tail(Rt_df, 1)$date) %>% head(., days_ahead)

  comp_df <- sim %>% select(days, infected_day) %>% filter(days > params$tau_m) %>%
    cbind(actual_cases %>% select(date, actual_cases))
  colnames(comp_df) <- c('Days', 'Forecast', 'Date', 'Actual')

  comp_df %>% select(Date, Forecast, Actual) %>% reshape2::melt(id.vars = 'Date') %>%
    ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
    theme_minimal() +
    labs(x = '', y = '') +
    theme(legend.title=element_blank()) +
    ggsave(paste0(outpath, plotname, '.png'), width = 10, height = 5)

  return(sim)

}

forecast_rt_var <- function(data, Rt_df, start_date, end_date, params, Rt_window, days_ahead, plotname){
  dat <- data %>% filter(date > as.Date(start_date) & date < as.Date(end_date))

  # take average and forecast next month
  params['R_val'] <- mean(tail(Rt_df, Rt_window)$Est_Rt1, na.rm = T)
  params['R_val_variant'] <- mean(tail(Rt_df, Rt_window)$Est_Rt2, na.rm = T)

  params['pop'] <- 67000000 # roughly french population
  params['init_infec'] <- list(tail(dat, params$tau_m)$I1_daily)
  params['init_infec_var'] <- list(tail(dat, params$tau_m)$I2_daily)
  params['n_days'] <- days_ahead + params$tau_m

  sim <- sii_sim(params)

  actual_cases <- data %>%
      filter(date > tail(Rt_df, 1)$date) %>% head(., days_ahead)
  colnames(actual_cases) <- c('date', 'infected_day_actual', 'I1_actual', 'I2_actual', 'days')

  comp_df <- sim %>% select(days, infected_day, I1_daily, I2_daily) %>% filter(days > params$tau_m) %>%
      cbind(actual_cases %>% select(date, infected_day_actual, I1_actual, I2_actual))
  colnames(comp_df) <- c('Days', 'ForecastTotal', 'ForecastI1', 'ForecastI2', 'Date', 'ActualTotal', 'ActualI1', 'ActualI2')

  comp_df %>% select(Date, ForecastTotal, ForecastI1, ForecastI2, ActualTotal, ActualI1, ActualI2) %>% reshape2::melt(id.vars = 'Date') %>%
      ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
      theme_minimal() +
      labs(x = '', y = '') +
      theme(legend.title=element_blank()) +
      ggsave(paste0(outpath, plotname, '.png'), width = 10, height = 5)

  return(sim)

}


# autumn
Rt_aut <- apply_rt_est(data = df, start_date = '2020-07-15', end_date = '2020-10-01', params = params, plotname = 'Rt_autumn')
Inf_aut <- forecast_rt(data = df, Rt_df = Rt_aut, start_date = '2020-07-15', end_date = '2020-10-01',
                       params = params, Rt_window = 15, days_ahead = 20, plotname = 'Forecast_autumn')

# winter
Rt_win <- apply_rt_est(data = df, start_date = '2020-10-30', end_date = '2020-12-01', params = params, plotname = 'Rt_winter')
Inf_win <- forecast_rt(data = df, Rt_df = Rt_win, start_date = '2020-10-30', end_date = '2020-12-01',
                       params = params, Rt_window = 5, days_ahead = 20, plotname = 'Forecast_winter')



# third application: variant starting 2021
df <- read.csv2(paste0(imppath, 'variants_france.csv')) %>% select(semaine, cl_age90, Nb_tests_PCR_TA_crible, Nb_susp_501Y_V1, Nb_susp_501Y_V2_3, Nb_susp_IND, Nb_susp_ABS)

colnames(df) <- c('week', 'age_grp', 'total_cases', 'brit_cases', 'braz_cases', 'unkw_cases', 'orig_cases')

# helper function to get date from week
substrRight <- function(x){
  substr(x, nchar(x)-10+1, nchar(x))
}

# aggregate age groups and divide by 7 to get daily
df <- df %>% group_by(week) %>%
  summarise(total_cases = sum(total_cases)/7,
          brit_cases = sum(brit_cases)/7,
          braz_cases = sum(braz_cases)/7,
          unkw_cases = sum(unkw_cases)/7,
          orig_cases = sum(orig_cases)/7) %>%
  ungroup() %>%
  mutate(date = as.Date(substrRight(week))) %>%
  select(date, total_cases, brit_cases, braz_cases, unkw_cases, orig_cases)

df_var <- df %>% select(date, total_cases, orig_cases, brit_cases)
colnames(df_var) <- c('date', 'infected_day', 'I1_daily', 'I2_daily')
df_var$days <- 1:nrow(df_var)

params['start_variant'] <- 0
start_date <- '2021-02-17'
end_date <- '2021-03-19'
df_var_sub <- df_var %>% filter(date > as.Date(start_date) & date < as.Date(end_date))


Rt_var <- Rt_est_applied(df_var_sub, params, variant = T)
Rt_var$date <- df_var_sub$date

Rt_var %>% select(Date, Est_Rt, Est_Rt1, Est_Rt2) %>% na.omit() %>%
  reshape2::melt(id.vars = 'Date') %>%
  ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
  ggsave(paste0(outpath, 'Rt_variants.png'))



Inf_var <- forecast_rt_var(data = df_var, Rt_df = Rt_var, start_date = start_date, end_date = end_date,
                           params = params, Rt_window = 15, days_ahead = 15, plotname = 'Forecast_variant')


summary(ur.df(Rt_summer$Est_Rt %>% na.omit(), lags=2, type='drift'))
summary(ur.df(Rt_winter$Est_Rt %>% na.omit(), lags=2, type='drift'))

kpss.test(Rt_summer$Est_Rt %>% na.omit(), null='Level', lshort=TRUE)
kpss.test(Rt_winter$Est_Rt %>% na.omit(), null='Level', lshort=TRUE)