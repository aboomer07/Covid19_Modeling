# Applying our estimators to French data
# Created by: jacobpichelmann
# Created on: 15.04.21

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/Params.R"))



# first application: late summer 2020 (august, sept)
df <- read.csv2(paste0(imppath, 'open_stats_coronavirus.csv')) %>% filter(nom == 'france') %>% select(date, cas) %>%
  mutate(new_cases = cas-lag(cas),
         date = as.Date(date))
colnames(df) <- c('date', 'cum_cases', 'infected_day') # need correct col names for function

# 2020-06-02 must be a data issue as cumulative cases drop by roughly 700

df_sum <- df %>% filter(date > as.Date('2020-07-15') & date < as.Date('2020-10-01'))
df_sum$days <- 1:nrow(df_sum)

# simulate serial interval study
serinfect <- samp_pois(params)
samps <- serinfect$daily

vals <- serial_ests(samps) # we are not using these values but we need it as an input (ugly code my bad)

Rt_summer <- Rt_est(df_sum, vals, params, deterministic = T, correct_bias = F, variant = F)
Rt_summer$date = df_sum$date
Rt_summer %>% select(date, Est_Rt) %>% na.omit() %>%
  ggplot() + geom_line(aes(x = date, y = Est_Rt))

# second application: winter 2020 (november, december) continuous lockdown (?)

df_win <- df %>% filter(date > as.Date('2020-10-15') & date < as.Date('2021-01-01'))
df_win$days <- 1:nrow(df_win)

Rt_winter <- Rt_est(df_win, vals, params, deterministic = T, correct_bias = F, variant = F)
Rt_winter$date = df_win$date
Rt_winter %>% select(date, Est_Rt) %>% na.omit() %>%
  ggplot() + geom_line(aes(x = date, y = Est_Rt))

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

Rt_var <- Rt_est(df_var, vals, params, deterministic = T, correct_bias = F, variant = T, sep_Rt = T)
# still need to fix NA assignment based on simulation data, line 349 in SimFunc.R




