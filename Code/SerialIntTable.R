# Title     : TODO
# Objective : TODO
# Created by: Luca Poll
# Created on: 5/4/2021

serial_int_study <- function(params){

  # setup
  sim_type <- params$sim_type
  R_type <- params$R_type
  study_len <- params$study_len
  df <- data.frame(matrix(nrow = 0, ncol = 12))
  colnames(df) <- c("sim_type", "R_type", "sample_size", "sim_nr",
                   "mu_gamma", "mu_nonpar", "var_gamma", "var_nonpar",
                   "MSE_gamma", "MSE_nonpar", "ME_gamma", "ME_nonpar")


  # loop over sample sizes
  sample_size <- c(10, 20, 50, 100, 500)



  for (s in seq_along(sample_size)){

    print(paste("Sample Size:", sample_size[s]))
    pb <- txtProgressBar(min = 1, max = 1000, style = 3)

    params$num_people <- sample_size[s]

    data <- data.frame(matrix(nrow = 0, ncol = 12))
    colnames(data) <- c("sim_type", "R_type", "sample_size", "sim_nr",
                   "mu_gamma", "mu_nonpar", "var_gamma", "var_nonpar",
                   "MSE_gamma", "MSE_nonpar", "ME_gamma", "ME_nonpar")

    # repeat simulation 1000 times
    for (i in 1:1000){
      setTxtProgressBar(pb, i)

      data[i,]$sim_nr <- i

      # 1) Serial Interval simulation
      samps <- samp_pois(params)$daily

      # 2) Fit gamma distr to serial interval
      vals <- serial_ests(samps)

      # 3) Fit nonparametric distr to serial interval
      vals_nonpar <- serial_ests_nonpara(samps, range = c(1, study_len), 'nsr')

      # 4) Simulate TSI pandemic
      infect <- nour_sim_data(params)

      # 5) Estimate R(t) and discard burn-in (40 days)
      est_gamma <- Rt_est(infect, vals, params)
      est_nonpar <- Rt_est_nonpara(infect, samps,"nsr", params)

      # 6) Compute comparison parameters
      # 6.1 estimated mean of omega hat
      data$mu_gamma[i] <- vals$meanhat
      data$mu_nonpar[i] <- sum(vals_nonpar$x * vals_nonpar$y)

      # 6.2 estimated variance of omega hat
      data$var_gamma[i] <- vals$varhat
      data$var_nonpar[i] <- sum(vals_nonpar$x^2 * vals_nonpar$y) - sum(vals_nonpar$x * vals_nonpar$y)^2

      # 6.3 MSE
      true_Rt <- Rt_gen(params$R_type, params$n_days)[(params$tau_m*2+1):params$n_days]
      gamma_Rt <- est_gamma$Est_Rt[(params$tau_m*2+1):params$n_days]
      nonpar_Rt <- est_nonpar$Est_Rt[(params$tau_m*2+1):params$n_days]

      data$MSE_gamma[i] <- mean((true_Rt - gamma_Rt)^2)
      data$MSE_nonpar[i] <- mean((true_Rt - nonpar_Rt)^2)

      # 6.4 ME
      data$ME_gamma[i] <- mean(true_Rt - gamma_Rt)
      data$ME_nonpar[i] <- mean(true_Rt - nonpar_Rt)
    }



    data$sample_size <- sample_size[s]
    df <- rbind(df, data)
  }


  df$sim_type <- sim_type
  df$R_type <- R_type

  return(df)

}



create_table <- function(params){

  # loop over true distribution types and different shapes of Rt

  # setup
  df <- data.frame(matrix(nrow = 0, ncol = 12))
  colnames(df) <- c("sim_type", "R_type", "sample_size", "sim_nr",
                   "mu_gamma", "mu_nonpar", "var_gamma", "var_nonpar",
                   "MSE_gamma", "MSE_nonpar", "ME_gamma", "ME_nonpar")

  serialint <- c("gamma", "weibull", "norm")
  rtypes <- c("constant", "increasing", "decreasing", "cave", "panic")

  for (s in serialint){

    params$sim_type <- s
    print(paste("True serial interval:", s))

    for (r in rtypes){

      params$R_type <- r
      print(paste("R type:", r))

      data <- serial_int_study(params)

      df <- rbind(df, data)
    }
  }

  return(df)
}


serial_table <- create_table(params)

save(serial_table, file = paste0(imppath, "SerialIntervalTable.R"))


tabl <- serial_table %>% group_by(sim_type, R_type, sample_size) %>%
  summarise(Mean_mu_gamma = mean(mu_gamma),
            "(SD)_mu_gamma" = sd(mu_gamma),
            Mean_mu_nonpar = mean(mu_nonpar),
            "(SD)_mu_nonpar" = sd(mu_nonpar),
            Mean_var_gamma = mean(var_gamma),
            "(SD)_var_gamma" = sd(var_gamma),
            Mean_var_nonpar = mean(var_nonpar),
            "(SD)_var_nonpar" = sd(var_nonpar),
            Mean_MSE_gamma = mean(MSE_gamma),
            "(SD)_MSE_gamma" = sd(MSE_gamma),
            Mean_MSE_nonpar = mean(MSE_nonpar),
            "(SD)_MSE_nonpar" = sd(MSE_nonpar),
            Mean_ME_gamma = mean(ME_gamma),
            "(SD)_ME_gamma" = sd(ME_gamma),
            Mean_ME_nonpar = mean(ME_nonpar),
            "(SD)_ME_nonpar" = sd(ME_nonpar))

tabl_long <- tabl %>% pivot_longer(!c(sim_type, R_type, sample_size), names_to="names",  values_to = "values") %>%
  mutate(measure = str_extract(names, "[^_]+"),
         type = gsub("_", "", str_extract(names, "_(.*?)_")),
         estimator = gsub("^.+_", "", names)) %>%
  select(-names) %>% mutate(R_type = ifelse(R_type == "constant", "Constant",
                                      ifelse(R_type == "increasing", "Increasing",
                                       ifelse(R_type == "decreasing", "Decreasing",
                                        ifelse(R_type == "cave", "U-inverted",
                                         ifelse(R_type == "panic", "Polynomial", "ERROR")))))) %>%
  mutate(R_type = factor(R_type, levels = c("Constant", "Increasing", "Decreasing", "Polynomial", "U-inverted"))) %>%
  filter(type %in% c("mu", "var", "MSE"))

tabl_wide <- tabl_long %>% pivot_wider(names_from = c(sim_type, type, estimator), values_from = values) %>%
  relocate(c(measure, R_type, sample_size), .before = gamma_mu_gamma) %>% arrange(R_type, sample_size)

is.num <- sapply(tabl_wide, is.numeric)
tabl_wide[is.num] <- lapply(tabl_wide[is.num], round, 4)


library(xtable)




print(xtable(tabl_wide, auto = T), include.rownames = F)

print(xtable(tabl_wide, digits = c(rep(0, 4), rep(c(6,6,2,2,3,3), 3))), include.rownames = F)











