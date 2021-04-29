params <- list()
params[['sim_mu']] <- 7
params[['sim_var']] <- 2
params[['study_len']] <- 20
params[['num_people']] <- 30
params[['simulations']] <- 1000
params[['sim_type']] <- 'norm'
params[['sim_type']] <- 'gamma'
params[['n_days']] <- 300
params[['delta']] <- 1
params[['n']] <- params[['tau_m']] * params[['delta']]
params[['R_val']] <- c(1.7, 0.9, 1.3)
params[['pop']] <- 6000000
params[['R_val_variant']] <- 1.9
params[['sim_mu_variant']] <- 7
params[['sim_var_variant']] <- 2
params[['sim_type_variant']] <- 'gamma'
params[['tau_m']] <- params[['study_len']]
params[['start_variant']] <- 115
params[['sii_cross']] <- 0
params[['infec_lim']] <- c(0, 1)
params[['Rt_lim']] <- c(0, 3)
params[['init_infec']] <- rep(1, params[['tau_m']] * params[['delta']])
