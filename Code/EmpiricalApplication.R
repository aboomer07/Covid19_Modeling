# Applying our estimators to French data
# Created by: jacobpichelmann
# Created on: 15.04.21

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/SIIFunc.R"))
source(paste0(getwd(), "/Code/Params.R"))



# first application: late summer 2020 (august, sept)
df <- read.csv2(paste0(imppath, 'open_stats_coronavirus.csv')) %>% filter(nom == 'france') %>% select(date, cas) %>%
  mutate(new_cases = cas-lag(cas),
         date = as.Date(date))
colnames(df) <- c('date', 'cum_cases', 'infected_day') # need correct col names for function

# 2020-06-02 must be a data issue as cumulative cases drop by roughly 700

df_sum <- df %>% filter(date > as.Date('2020-07-15') & date < as.Date('2020-10-01'))

# compute weekly rolling average to deal with weekend bias
df_sum <- df_sum %>% mutate(
  infected_day = zoo::rollmean(infected_day, k = 7, fill = NA)
) %>% na.omit()
df_sum$days <- 1:nrow(df_sum)


# simulate serial interval study
serinfect <- samp_pois(params)
samps <- serinfect$daily

vals <- serial_ests(samps) # we are not using these values but we need it as an input (sorry ugly code my bad)

Rt_summer <- Rt_est(df_sum, vals, params, deterministic = T, correct_bias = F, variant = F)
Rt_summer$date = df_sum$date
Rt_summer %>% select(date, Est_Rt) %>% na.omit() %>%
  ggplot() + geom_line(aes(x = date, y = Est_Rt)) +
  ggsave(paste0(outpath, 'Rt_summer.png'))

# take average and forecast next month
params['R_val'] <- mean(Rt_summer$Est_Rt, na.rm = T)
params['pop'] <- 67000000
params['init_infec'] <- list(tail(df_sum, params$tau_m)$infected_day)
params['n_days'] <- 50

sim_summer <- si_sim(params)

pdf(file = paste0(outpath, "SI_plot_summer.pdf"), width=4, height=7)
si_plot_detail(sim_summer)
dev.off()

actual_cases <- df %>%
  mutate(actual_cases = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit() %>%
  filter(date > as.Date('2020-10-01')) %>% head(., 30)

comp_df <- sim_summer %>% select(days, infected_day) %>% filter(days > 20) %>%
  cbind(actual_cases %>% select(date, actual_cases))
colnames(comp_df) <- c('Days', 'Forecast', 'Date', 'Actual')

comp_df %>% select(Date, Forecast, Actual) %>% reshape2::melt(id.vars = 'Date') %>%
  ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
  ggsave(paste0(outpath, 'Forecast_summer.png'))


# second application: winter 2020 (november, december) continuous lockdown (?)

df_win <- df %>% filter(date > as.Date('2020-10-15') & date < as.Date('2021-01-01'))

# compute weekly rolling average to deal with weekend bias
df_win <- df_win %>% mutate(
  infected_day = zoo::rollmean(infected_day, k = 7, fill = NA)
) %>% na.omit()
df_win$days <- 1:nrow(df_win)

Rt_winter <- Rt_est(df_win, vals, params, deterministic = T, correct_bias = F, variant = F)
Rt_winter$date = df_win$date
Rt_winter %>% select(date, Est_Rt) %>% na.omit() %>%
  ggplot() + geom_line(aes(x = date, y = Est_Rt)) +
  ggsave(paste0(outpath, 'Rt_winter.png'))

# take average and forecast next month
params['R_val'] <- mean(Rt_winter$Est_Rt, na.rm = T)
params['pop'] <- 67000000
params['init_infec'] <- list(tail(df_win, params$tau_m)$infected_day)
params['n_days'] <- 50

sim_winter <- si_sim(params)

pdf(file = paste0(outpath, "SI_plot_winter.pdf"), width=4, height=7)
si_plot_detail(sim_winter)
dev.off()

actual_cases <- df %>%
  mutate(actual_cases = zoo::rollmean(infected_day, k = 7, fill = NA)) %>% na.omit() %>%
  filter(date > as.Date('2021-01-01')) %>% head(., 30)

comp_df <- sim_winter %>% select(days, infected_day) %>% filter(days > 20) %>%
  cbind(actual_cases %>% select(date, actual_cases))
colnames(comp_df) <- c('Days', 'Forecast', 'Date', 'Actual')

comp_df %>% select(Date, Forecast, Actual) %>% reshape2::melt(id.vars = 'Date') %>%
  ggplot() + geom_line(aes(x = Date, y = value, color = variable)) +
  ggsave(paste0(outpath, 'Forecast_winter.png'))



# third application: variant starting 2021
df <- read.csv2(paste0(imppath, 'variants_france.csv')) %>% select(semaine, cl_age90, Nb_tests_PCR_TA_crible,
                                                                   Nb_susp_501Y_V1, Nb_susp_501Y_V2_3, Nb_susp_IND,
                                                                   Nb_susp_ABS)

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

sim_winter <- si_sim(params)






