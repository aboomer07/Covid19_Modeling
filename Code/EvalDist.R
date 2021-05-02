# Objective : Evaluate fitted distributions on serial interval data
# Created by: jacobpichelmann
# Created on: 21.03.21

source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/Params.R"))

######################################################################
#############       Serial Interval Plots        #####################
######################################################################

# Continuous simulation serial interval histogram
serial_hist_cont <- function (sampscont, dist){
	hist(sampscont, breaks = length(unique(sampscont)), col = "lightblue",
		 border = "lightgrey", freq = F, xlab = "Time", main=NULL)
	lines(dist, lwd=1.5)
	text(quantile(sampscont, probs=0.98), max(dist)*1.3,
	 paste("Distribution =", params$sim_type, "\nMean =", params$sim_mu,
	 		 "; Var =", params$sim_var, "\nDelta = 1/", params$delta,
	 		 "\nSample size =", params$num_people))
}

# discretized serial interval histogram
serial_hist_disc <- function (samps){
	omegadaily <- gen_distribution(params$study_len, params$sim_mu, 
		params$sim_var, params$sim_type, delta = 1)
	hist(samps, breaks=seq(min(samps)-0.5, max(samps)+0.5, by=1), freq = F, 
		col="lightblue", xlab = "Day", main = NULL)
	lines(omegadaily$omega, lwd = 1.5)
	text(quantile(samps, probs=0.98), max(omegadaily$omega)*0.9,
	 paste("Distribution =", params$sim_type, "\nMean =", params$sim_mu,
	 	 "; Var =", params$sim_var, "\nDelta =", 1,
	 	 	"\nSample size =", params$num_people))
}



#####################################################################################
############## Serial Interval Simulation and Estimation Plot  ######################
#####################################################################################

SI_plot_distribution <- function(data){

	R_val <- params[['R_val']]; study_len <- params[['study_len']]
  	num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  	sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  	delta <- params[['delta']]; simulations <- params[['simulations']]

	estimates <- data$distribution
	avg_shape_hat <- data$avg_params[1]
	avg_shape_sd <- data$avg_params[2]
	avg_rate_hat <- data$avg_params[3]
	avg_rate_sd <- data$avg_params[4]
	avg_mean_hat <- data$avg_params[5]
	avg_mean_sd <- data$avg_params[6]
	avg_var_hat <- data$avg_params[7]
	avg_var_sd <- data$avg_params[8]

	print(avg_mean_hat)

	#Show distribution of estimates 
	#pdf(file = paste0(outpath, "SerialEst_", sim_type, "_S", simulations, "_Delta", delta, ".pdf"))
	par(mfrow = c(2,2))

	#Alpha hat 
	hist(estimates$shape_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", alpha)))
	abline(v = avg_shape_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	legend("topright", 
	   c(expression(paste("Mean ", hat(alpha))), "Avg 95% CI"),
	   lty = c(1, 2, 1),  
	   col = c("blue", "blue", "red"),
	   cex = 0.75)

	#Beta hat
	hist(estimates$rate_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", beta)))
	abline(v = avg_rate_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(beta))), "Avg 95% CI"), 
		   lty = c(1, 2, 1),  
		   col = c("blue", "blue", "red"),

		   cex = 0.75)
	#Implied mean
	hist(estimates$mean_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", mu)))
	abline(v = avg_mean_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	#abline(v = (avg_mean_hat-1.96*avg_mean_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = (avg_mean_hat+1.96*avg_mean_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = sim_mu, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(mu))), expression(paste("True ", mu))),
		   lty = c(1, 1),  
		   col = c("blue", "red"),
		   cex = 0.75)

	#Implied variance
	hist(estimates$var_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", sigma^2)))
	abline(v = avg_var_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	#abline(v = (avg_var_hat-1.96*avg_var_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = (avg_var_hat+1.96*avg_var_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = sim_var, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(sigma^2))), expression(paste("True ", sigma^2))),
		   lty = c(1, 1),  
		   col = c("blue", "red"),
		   cex = 0.75)

	#Title
	mtext(paste0("Number of simulations = ", simulations,
	 	"\nSample size = ", num_people,
	 	"\nDiscretization = 1/", delta,
	 	"\nUnderlying distribution = ", sim_type), side = 3, line = -4, outer = TRUE, cex = 0.7)
	mtext(paste0("Number of simulations = ", simulations, 
		"\n Sample size = ", num_people,
		"\nDiscretization = 1/", delta,
		"\nUnderlying distribution = ", sim_type), side = 3, line = -24, outer = TRUE, cex = 0.7)
	#dev.off()
}


SI.bias <- function(data){

	R_val <- params[['R_val']]; study_len <- params[['study_len']]
  	num_people <- params[['num_people']]; sim_mu <- params[['sim_mu']]
  	sim_var <- params[['sim_var']]; sim_type <- params[['sim_type']]
  	delta <- params[['delta']]; simulations <- params[['simulations']]

	estimates <- data$distribution
	avg_shape_hat <- data$avg_params[1]
	avg_shape_sd <- data$avg_params[2]
	avg_rate_hat <- data$avg_params[3]
	avg_rate_sd <- data$avg_params[4]
	avg_mean_hat <- data$avg_params[5]
	avg_mean_sd <- data$avg_params[6]
	avg_var_hat <- data$avg_params[7]
	avg_var_sd <- data$avg_params[8]
)


serial_est_plot <- function(study_len, sim_mean, sim_var, sim_type, vals, nonpara = F){

	sim_mean <- params[['sim_mu']]; sim_var <- params[['sim_var']]; 
	sim_type <- params[['sim_type']]; study_len <- params[['study_len']]

	true <- gen_distribution(params$study_len, sim_mean, sim_var, sim_type, 1)

	if (!nonpara){
		est <- dgamma(1:params$study_len, vals$shape, vals$rate)
		df_e <- as.data.frame(est)
	}

	if (nonpara){
		est <- vals # here vals has to b already fitted data obtained from serial_est_nonpara()
		df_e <- as.data.frame(est$y)
	}
	
	plot(true$omega, type = "l", main = NULL, xlab = "Day", ylab = "Density", lwd=2)
	lines(est, col = "red", lwd=2)

	df_t<-as.data.frame(true$omega)
	df<-bind_cols(df_t, df_e)
	MSE<-MSE_est2(df)

	if (!nonpara){
	legend("topright",
		   legend = c(paste0("True Serial Interval (", sim_type,")"),
					  paste0("Estimated Serial Interval (gamma)")),
		   col = c("black", "red"), lty=1, cex=0.9)
	}

	if (nonpara){
	legend("topright",
		   legend = c(paste0("True Serial Interval (", sim_type,")"),
					  paste0("Estimated Serial Interval (Kernel Estimation)")),
		   col = c("black", "red"), lty=1, cex=0.9)
	}
	
	legend("right", legend=c(paste0("MSE: ", MSE)), pt.cex = 0,pch = c(17,19), cex=0.8)
	
}

# second function that does all at once so there is no confusion between samps creation and simulation parameters for
# true dist

serial_est_plot_full <- function(params, nonpara = F, bw = NULL){

	sim_mean <- params[['sim_mu']]; sim_var <- params[['sim_var']]; 
	sim_type <- params[['sim_type']]; study_len <- params[['study_len']]
	R_val <- params[['R_val']]

	true <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1)

	samps <- samp_pois(params)$daily

	if (!nonpara){
		vals <- serial_ests(samps)
		est <- dgamma(1:study_len, vals$shape, vals$rate)
		df_e <- as.data.frame(est)

	}

	if (nonpara){
		est <- serial_ests_nonpara(samps, range = c(0, study_len), bandwidth = bw) # here vals has to b already fitted data obtained from serial_est_nonpara()
		df_e <- as.data.frame(est$y)
	}

	plot(true$omega, type = "l", main = NULL, xlab = "Day", ylab = "Density", lwd=2)
	lines(est, col = "red", lwd=2)

	df_t<-as.data.frame(true$omega)
	df<-bind_cols(df_t, df_e)
	MSE<-MSE_est2(df)

	if (!nonpara){
	legend("topright",
		   legend = c(paste0("True Serial Interval (", sim_type,")"),
					  paste0("Estimated Serial Interval (gamma)")),
		   col = c("black", "red"), lty=1, cex=0.9)
	}

	if (nonpara){
	legend("topright",
		   legend = c(paste0("True Serial Interval (", sim_type,")"),
					  paste0("Estimated Serial Interval (Bandwidth: ", bw, ")")),
		   col = c("black", "red"), lty=1, cex=0.9)
	}

	legend("right", legend=c(paste0("MSE: ", MSE)), pt.cex = 0,pch = c(17,19), cex=0.8)

}

######################################################################
##############         Simulate Incidence       ######################
######################################################################

infections_plot <- function(incid, params, variant=F, ssii=F){

	if (variant) {
		title <- paste0("New infections per day (", params$sim_type, ")")
		if (ssii) {title <- paste0(title, " CrossOver = ", params$sii_cross)}

		layout(matrix(1:2, ncol = 1, nrow =2))
		plot(incid$I1_pct, type = "l", 
		main = title, xlab = "", ylab = "New infections",
		ylim=params$infec_lim)
		lines(incid$I2_pct, col='red')
		legend('topleft',
			legend=c("Original Strain", 'Variant'),
			col=c("black", 'red'), lty=1, cex=0.9)
		plot(incid$R1_val, type = "l", main = "Simulated R(t)", 
		xlab = "Day", ylab = "R(t)", ylim = params$Rt_lim)
		lines(incid$R2_val, main = "Simulated R(t) Variant", 
		xlab = "Day", ylab = "R(t)", col='red')
		legend('bottomleft',
			legend=c("Original Strain", 'Variant'),
			col=c("black", 'red'), lty=1, cex=0.9)
	}
	else {
	layout(matrix(1:2, ncol = 1, nrow =2))
	plot(incid$I_pct, type = "l", 
		main = paste0("New infections per day (", params$sim_type, ")"),
		xlab = "", ylab = "New infections", ylim=params$infec_lim)
	plot(incid$R_val, type = "l", main = "Simulated R(t)", 
		xlab = "Day", ylab = "R(t)", ylim = params$Rt_lim)
	}
}

### Draw daily incidence of cases and corresponding Rt, given type of omega 
### Compare estimate of Rt to its true value 

compare_rt <- function(Rt, params = NULL, variant = F){
	start_variant <- params[['start_variant']] + params[['tau_m']]
	n_days <- params[['n_days']]

	if(variant){
		Rt2 <- c(rep(NA, start_variant), Rt$Rt2[start_variant + 1:n_days])
		plot(Rt$Rt1, type = "l", 
			#main = paste0("R(t) estimation, true distribution: ", params[['sim_type']]),
			 xlab = "Day", ylab = "R(t)",
			col = "red", ylim = c(0.5, 2.5), lwd = 2)
		lines(Rt2, col = "orange", lwd=2)
		lines(Rt$Est_Rt, col = "black", lwd=2)
		abline(v = start_variant, lwd = 2, lty = "dotted")
		legend("topleft", 
			legend = c("True R(t)", "True R(t) Variant", "Estimated R(t)"),
			col = c("red", "orange", "black"), lty=1, cex=0.9, bty="n")
	}

	else{
		plot(Rt$Rt, type = "l", 
			#main = paste0("R(t) estimation, true distribution: ", params[['sim_type']]),
			xlab = "Day", ylab = "R(t)",
			col = "red", ylim = c(0.5, 2.5), lwd = 2)
		lines(Rt$Est_Rt, col = "black", lwd=2)
		legend("topright", legend = c("True R(t)", "Estimated R(t)"),
		       col = c("red", "black"), lty=1, cex=0.9,bty="n")
	}
}

