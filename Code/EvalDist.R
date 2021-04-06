# Objective : Evaluate fitted distributions on serial interval data
# Created by: jacobpichelmann
# Created on: 21.03.21

source(paste0(getwd(), "/Code/SimFunc.R"))


######################################################################
#############       Serial Interval Plots        #####################
######################################################################

# Continuous simulation serial interval histogram
serial_hist_cont <- function (sampscont, dist){
	hist(sampscont, breaks = length(unique(sampscont)), col = "lightblue",
		 border = "lightgrey", freq = F, xlab = "Time", main=NULL)
	lines(dist, lwd=1.5)
	text(quantile(sampscont, probs=0.98), max(dist)*1.3,
	 paste("Distribution =", sim_type, "\nMean =", sim_mean, "; Var =", sim_var,
		   "\nDelta =", delta, "\nStudy length =", study_len,
		   "\nSample size =", num_people))
}

# discretized serial interval histogram
serial_hist_disc <- function (samps){
	omegadaily <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1)
	hist(samps, breaks=seq(min(samps)-0.5, max(samps)+0.5, by=1), freq = F,
	 	 col="lightblue", xlab = "Day", main = NULL)
	lines(omegadaily, lwd = 1.5)
	text(quantile(samps, probs=0.98), max(omegadaily)*0.9,
	 paste("Distribution =", sim_type, "\nMean =", sim_mean, "; Var =", sim_var,
		   "\nStudy length =", study_len, "\nSample size =", num_people))
}



#####################################################################################
############## Serial Interval Simulation and Estimation Plot  ######################
#####################################################################################

#Always run this and manually change the parameters in the graph... 
true_gamma <- gen_distribution(study_len, sim_mean, sim_var, "gamma", 1)

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
	true <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1)

	if (!nonpara){
		est <- dgamma(1:study_len, vals$shape, vals$rate)
		df_e <- as.data.frame(est)

	}

	if (nonpara){
		est <- vals # here vals has to b already fitted data obtained from serial_est_nonpara()
		df_e <- as.data.frame(est$y)
	}
	
	plot(true, type = "l", main = NULL, xlab = "Day", ylab = "Density", lwd=2)
	lines(est, col = "red", lwd=2)

	df_t<-as.data.frame(true)
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

serial_est_plot_full <- function(study_len, sim_mean, sim_var, sim_type, R_val, nonpara = F, bw = NULL){
	true <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1)

	samps <- samp_pois(R_val = 1.6, study_len = study_len,
					   num_people = num_people, sim_mu = sim_mean,
					   sim_sig = sim_var, sim_type = sim_type, delta = delta)$daily

	if (!nonpara){
		vals <- serial_ests(samps)
		est <- dgamma(1:study_len, vals$shape, vals$rate)
		df_e <- as.data.frame(est)

	}

	if (nonpara){
		est <- serial_ests_nonpara(samps, range = c(0, study_len), bandwidth = bw) # here vals has to b already fitted data obtained from serial_est_nonpara()
		df_e <- as.data.frame(est$y)
	}

	plot(true, type = "l", main = NULL, xlab = "Day", ylab = "Density", lwd=2)
	lines(est, col = "red", lwd=2)

	df_t<-as.data.frame(true)
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







### Draw daily incidence of cases and corresponding Rt, given type of omega
#


# infections_plot(type = 'gamma')


### Compare estimate of Rt to its true value 

compare_rt <- function(Rt){
	plot(Rt$Rt, type = "l", main = paste0("R(t) estimation, true distribution: ", sim_type),
		 xlab = "Day", ylab = "R(t)", col = "red", ylim = c(0, max(Rt$Est_Rt, na.rm=T)), lwd = 2)
	lines(Rt$Est_Rt, col = "black", lwd=2)
	legend("topright", legend = c("True R(t)", "Estimated R(t)"),
	       col = c("red", "black"), lty=1, cex=0.9)
}





######################################################################
#############               LEGACY               #####################
######################################################################

# # # # Serial Interval

#
# samp_pois_plot <- function(samps, study_len){
#   data.frame(table(samps)) %>%
#     mutate(samps = as.numeric(as.character(samps))) %>%
#     ggplot(aes(x=samps, y=Freq)) + geom_bar(stat = "identity") + xlab("days") + ylab("count") +
#       scale_x_continuous(limits = c(0,study_len)) + theme_minimal()
# }
#
# samp_pois_plot(samps = samps, study_len = study_len*delta)

### Fit gamma and compare it to true distribution -daily




# # # # - - - - - - - - - - - -

#og_gamma <- gen_distribution(15, sim_mu, sim_var, "gamma")
#new_gamma <- gen_distribution(15, mean(params[params$Distribution == 'gamma', 'True_a']), mean(params[params$Distribution == 'gamma', 'True_b']),"gamma")
#og_weibull <- gen_distribution(15, a_weibull, b_weibull, "weibull")
#new_weibull<- gen_distribution(15, mean(params[params$Distribution == 'weibull', 'True_a']), mean(params[params$Distribution == 'weibull', 'True_b']),"gamma")
#og_norm <- gen_distribution(15, a_norm, b_norm, "norm")
#new_norm <- gen_distribution(15, mean(params[params$Distribution == 'norm', 'True_a']), mean(params[params$Distribution == 'norm', 'True_b']),"gamma")
#og_lnorm <- gen_distribution(15, a_lnorm, b_lnorm, "lnorm")
#new_lnorm <- gen_distribution(15, mean(params[params$Distribution == 'lnorm', 'True_a']), mean(params[params$Distribution == 'lnorm', 'True_b']),"gamma")
#
#jpeg("DistCompare.jpg")
#layout(matrix(1:4, nrow = 2, ncol=2))
#plot(og_gamma, type="l")
#lines(new_gamma, type="l", col="green")
#title("Gamma")
#
#plot(og_weibull, type="l")
#lines(new_weibull, type="l", col="red")
#title("Weibull")
#
#plot(og_norm, type="l")
#lines(new_norm, type="l", col="blue")
#title("Normal")
#
#plot(og_lnorm, type="l")
#lines(new_lnorm, type="l", col="orange")
#title("Log-normal")
#dev.off()