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
		params$sim_var, params$sim_type, 1)
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


#Always run this and manually change the parameters in the graph... 
true_gamma <- gen_distribution(params$study_len, params$sim_mean, params$sim_var, "gamma", 1)

SI_plot_distribution <- function(data){

	estimates <- data$distribution
	avg_shape_hat <- data$avg_params[1]
	avg_shape_sd <- data$avg_params[2]
	avg_rate_hat <- data$avg_params[3]
	avg_rate_sd <- data$avg_params[4]
	avg_mean_hat <- data$avg_params[5]
	avg_mean_sd <- data$avg_params[6]
	avg_var_hat <- data$avg_params[7]
	avg_var_sd <- data$avg_params[8]

	#Show distribution of estimates 
	pdf(file = paste0(outpath, "SerialEst_", sim_type, "_S", simulations, ".pdf"))
	par(mfrow = c(2,2))
	#Alpha hat 
	hist(estimates$shape_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", alpha)))
	abline(v = avg_shape_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	#abline(v = (avg_shape_hat-1.96*avg_shape_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = (avg_shape_hat+1.96*avg_shape_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = 6.25, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
	   c(expression(paste("Mean ", hat(alpha)))),
	   lty = c(1, 2, 1),  
	   col = c("blue", "blue", "red"),
	   cex = 0.75)
	#Beta hat
	hist(estimates$rate_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", beta)))
	abline(v = avg_rate_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	#abline(v = (avg_rate_hat-1.96*avg_rate_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = (avg_rate_hat+1.96*avg_rate_sd), lwd = 2, lty = "twodash", col = "blue")
	#abline(v = 1.25, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(beta)))), 
		   lty = c(1, 2, 1),  
		   col = c("blue", "blue", "red"),
		   cex = 0.75)
	#Implied mean
	hist(estimates$mean_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", mu)))
	abline(v = avg_mean_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	abline(v = (avg_mean_hat-1.96*avg_mean_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = (avg_mean_hat+1.96*avg_mean_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = sim_mu, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(mu))), "95% CI", expression(paste("True ", mu))),
		   lty = c(1, 6, 1),  
		   col = c("blue", "blue", "red"),
		   cex = 0.75)
	#Implied variance
	hist(estimates$var_hat, nclass = 20, xlab = "", col = "aliceblue",
		main = expression(paste("Estimated ", sigma^2)))
	abline(v = avg_var_hat, lwd = 2, lty = "solid", col = alpha("blue", 0.7))
	abline(v = (avg_var_hat-1.96*avg_var_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = (avg_var_hat+1.96*avg_var_sd), lwd = 2, lty = "twodash", col = "blue")
	abline(v = sim_sig, lwd = 2, lty = "solid", col = alpha("red", 0.7))
	legend("topright", 
		   c(expression(paste("Mean ", hat(sigma^2))), "95% CI", expression(paste("True ", sigma^2))),
		   lty = c(1, 6, 1),  
		   col = c("blue", "blue", "red"),
		   cex = 0.75)
	#Title
	mtext(paste0("Number of simulations = ", simulations,
	 	"\nSample size = ", num_people,
	 	"\nDiscretization = 1/", delta), side = 3, line = -3, outer = TRUE, cex = 0.7)
	mtext(paste0("Number of simulations = ", simulations, 
		"\n Sample size = ", num_people,
		"\nDiscretization = 1/", delta), side = 3, line = -24, outer = TRUE, cex = 0.7)
	dev.off()
}


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

infections_plot <- function(incid){
	infected <- incid$infected_day
	Rt_sim <- incid$R_val

	layout(matrix(1:2, ncol = 1, nrow =2))
	plot(infected, type = "l", main = paste0("New infections per day (", sim_type, ")"),
		 xlab = "", ylab = "New infections")
	plot(Rt_sim, type = "l", main = "Simulated R(t)", xlab = "Day", ylab = "R(t)",
		 ylim = c(0, 2))
}

### Draw daily incidence of cases and corresponding Rt, given type of omega #

# infections_plot(type = 'gamma')

### Compare estimate of Rt to its true value 

compare_rt <- function(Rt){
	plot(Rt$Rt, type = "l", main = paste0("R(t) estimation, true distribution: ", params[['sim_type']]),
		 xlab = "Day", ylab = "R(t)", col = "red", ylim = c(min(c(min(Rt$Est_Rt, na.rm=T), min(Rt$Rt, na.rm = T))),
															max(Rt$Est_Rt, na.rm=T)), lwd = 2)
	lines(Rt$Est_Rt, col = "black", lwd=2)
	legend("topright", legend = c("True R(t)", "Estimated R(t)"),
	       col = c("red", "black"), lty=1, cex=0.9)
}


# evaluate performance of non parametric estimator

## sims is the number of simulations
#R_val <- 1.3; study_len <- 15; num_people <- 40; sim_mu <- 6.6; sim_sig <- 1.1;
#sim_type <- 'gamma'; delta <- 24; bw <- 'nsr'; sims <- 1000
#
#true_dist <- gen_distribution(study_len, sim_mu, sim_sig, sim_type, 1)
#
#est_dist <- nonpara_eval(params, bw)
#plot_nonpara_eval(true_dist, est_dist)

