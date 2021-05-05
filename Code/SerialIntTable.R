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
  sample_size <- c(10, 20) #, 50, 100, 500)

  for (s in seq_along(sample_size)){

    params$num_people <- sample_size[s]

    data <- data.frame(matrix(nrow = 0, ncol = 12))
    colnames(data) <- c("sim_type", "R_type", "sample_size", "sim_nr",
                   "mu_gamma", "mu_nonpar", "var_gamma", "var_nonpar",
                   "MSE_gamma", "MSE_nonpar", "ME_gamma", "ME_nonpar")

    # repeat simulation 1000 times
    for (i in 1:100){
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



test <- serial_int_study(params)