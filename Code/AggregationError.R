

aggregation_error <- function(params){

	true_mean <- params[["sim_mu"]]
	true_variance <- params[["sim_var"]]
	Rt_true <- Rt_gen(params[["R_type"]], params[["n_days"]]) 

	deltas <- c(1)
	#deltas <- c(1, 2, 4, 12, 24, 48, 96, 288, 1440)
	omegas <- c("weibull", "gamma", "norm")

	mean_hat <- data.frame(matrix(ncol = 3, nrow = length(deltas)))
	names(mean_hat) <- c("weibull", "gamma", "norm")
	var_hat <- data.frame(matrix(ncol = 3, nrow = length(deltas)))
	names(var_hat) <- c("weibull", "gamma", "norm")
	MSE_Rt <- data.frame(matrix(ncol = 3, nrow = length(deltas)))
	names(MSE_Rt) <- c("weibull", "gamma", "norm")

	for (o in seq(seq_along(omegas))){

		params[["sim_type"]] <- omegas[o]

		for (d in seq_along(deltas)){

			params[["delta"]] <- deltas[d]
			params[['init_infec']] <- rep(1, params[['tau_m']] * params[['delta']])

			estimated_params <- params_distribution(params)		
			mean_hat[d,o] <- estimated_params$avg_params$avg_mean_hat
			var_hat[d,o] <- estimated_params$avg_params$avg_var_hat

			si_estimates <- serial_ests(samp_pois(params)$daily)
			incidence <- nour_sim_data(params)
			Rt_est(incidence, si_estimates$vals, params)
			print("ok")
			MSE_Rt[d,o] <- mean((incidence$R_val-Rt_true)^2)


		}
	}

	output <- list(mean_hat = mean_hat, var_hat = var_hat, MSE_R_t = MSE_Rt, Rt = Rt_true)
	return(output)
}


aggregation_error(params)

incidence <- nour_sim_data(params)
si_estimates <- serial_ests(samp_pois(params)$daily)
Rt_est(incidence, si_estimates, params)