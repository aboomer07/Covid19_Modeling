setwd("..")
imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

source(paste0(getwd(), "/Code/Params.R"))
source(paste0(getwd(), "/Code/SimFunc.R"))


aggregation_error <- function(params){

	true_mean <- params[["sim_mu"]]
	true_variance <- params[["sim_var"]]
	study_len <- params[['study_len']]
	Rt_true <- Rt_gen(params[["R_type"]], params[["n_days"]]) 

	deltas <- seq(1:100)
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

			vals <- list(shape = estimated_params$avg_params$avg_shape_hat,
						 rate = estimated_params$avg_params$avg_rate_hat,
						 meanhat = estimated_params$avg_params$avg_mean_hat,
						 varhat = estimated_params$avg_params$avg_var_hat)

			incidence <- nour_sim_data(params)
			Rt_hat <- Rt_est(incidence, vals, params)
			Rt_hat <- Rt_hat$Est_Rt[((2*study_len)+1):length(Rt_hat$Est_Rt)]
			Rt_true_trunc <- Rt_true[((2*study_len)+1):length(Rt_true)]
			MSE_Rt[d,o] <- mean((Rt_hat-Rt_true_trunc)^2)

		}
	}

	output <- list(deltas = deltas, mean_hat = mean_hat, var_hat = var_hat, MSE_Rt = MSE_Rt, Rt = Rt_true)
	return(output)
}


aggregation_error_plots <- function(params){

		errors <- aggregation_error(params)
		par(mfrow = c(1,2))
		plot(errors$mean_hat$weibull~errors$deltas, type = "l", col = "red",
			xlab = expression(paste("1/", Delta)), ylab = expression(hat(mu)),
			ylim = c(params[['sim_mu']]-0.1, max(errors$mean_hat$weibull, errors$mean_hat$gamma, errors$mean_hat$norm)))
		#axis(1, at=errors$deltas, labels=errors$deltas)
		lines(errors$mean_hat$gamma~errors$deltas, col = "blue")
		lines(errors$mean_hat$norm~errors$deltas, col = "darkgoldenrod2")
		abline(h = params[['sim_mu']], col = "red", lty = 2)
		legend("bottomright",
			legend = c("weibull", "gamma", "normal", "true value"), 
			col = c("red", "blue", "darkgoldenrod2", "red"),
			lty = c(1, 1, 1, 2), 
			cex = 0.7,
			xpd=TRUE,
			horiz = TRUE,
			bty = "n",
			inset=c(0,1))

		plot(errors$MSE_Rt$weibull~errors$deltas, type = "l", col = "red",
			xlab = expression(paste("1/", Delta)), ylab = "MSE",
			ylim = c(min(errors$MSE_Rt$weibull, errors$MSE_Rt$gamma, errors$MSE_Rt$norm),
			 max(errors$MSE_Rt$weibull, errors$MSE_Rt$gamma, errors$MSE_Rt$norm)))
		#axis(1, at=errors$deltas, labels=errors$deltas)
		lines(errors$MSE_Rt$gamma~errors$deltas, col = "blue")
		lines(errors$MSE_Rt$norm~errors$deltas, col = "darkgoldenrod2")
		legend("bottomright",
			legend = c("weibull", "gamma", "normal"), 
			col = c("red", "blue", "darkgoldenrod2"),
			lty = c(1, 1, 1), 
			cex = 0.7,
			xpd=TRUE,
			horiz = TRUE,
			bty = "n",
			inset=c(0,1))

		return(errors)
}

#Specify type and run: constant, increasing, decreasing, panic, cave

types <- c("constant", "increasing", "decreasing", "panic", "cave")

for (type in seq_along(types)){
	params[['R_type']] <- types[type]
	pdf(file = paste0(outpath, 'SerialIntervalBias_Rt', types[type], '.pdf'),
		width=10, height=4)
	aggregation_error_plots(params)	
	dev.off()
}




