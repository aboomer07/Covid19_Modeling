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


apply_rt_est <- function(data, start_date, end_date, params, plotname){
  dat <- data %>% filter(date > as.Date(start_date) & date < as.Date(end_date))

  # compute weekly rolling average to deal with weekend bias
  dat <- dat %>% mutate(
    infected_day = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit()
  dat$days <- 1:nrow(dat)

  # simulate serial interval study
  serinfect <- samp_pois(params)
  samps <- serinfect$daily

  vals <- serial_ests(samps) # we are not using these values but we need it as an input (sorry ugly code my bad)

  Rt <- Rt_est(dat, vals, params, deterministic = T, correct_bias = F, variant = F)
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
Rt_var <- Rt_est(df_var, vals, params, deterministic = T, correct_bias = F, variant = T, sep_Rt = T)
# still need to fix NA assignment based on simulation data, line 349 in SimFunc.R

Rt_var %>% select(Date, Est_Rt, Est_Rt1, Est_Rt2) %>% na.omit() %>%
  reshape2::melt(id.vars = 'Date') %>%
  ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
  ggsave(paste0(outpath, 'Rt_variants.png'))


summary(ur.df(Rt_summer$Est_Rt %>% na.omit(), lags=2, type='drift'))
summary(ur.df(Rt_winter$Est_Rt %>% na.omit(), lags=2, type='drift'))

kpss.test(Rt_summer$Est_Rt %>% na.omit(), null='Level', lshort=TRUE)
kpss.test(Rt_winter$Est_Rt %>% na.omit(), null='Level', lshort=TRUE)