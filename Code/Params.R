params <- list()
params[['sim_mu']] <- 7
params[['sim_var']] <- 2
params[['study_len']] <- 20
params[['num_people']] <- 30
params[['simulations']] <- 1000
params[['sim_type']] <- 'gamma'
params[['n_days']] <- 300
params[['delta']] <- 24
params[['n']] <- params[['tau_m']] * params[['delta']]
params[['R_val']] <- c(1.4, 1.1, 1.2)
params[['pop']] <- 6000000
params[['R_val_variant']] <- 1.9
params[['sim_mu_variant']] <- 7
params[['sim_var_variant']] <- 2
params[['sim_type_variant']] <- 'gamma'
params[['tau_m']] <- params[['study_len']]
params[['start_variant']] <- 50
params[['sii_cross']] <- 0.1