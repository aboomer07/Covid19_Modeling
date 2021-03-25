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



######################################################################
############## Serial Interval Estimation Plot  ######################
######################################################################

serial_est_plot <- function(study_len, sim_mean, sim_var, sim_type, vals){
	omegadaily <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1)
	plot(omegadaily, type = "l", main = NULL, xlab = "Day", ylab = "Density", lwd=2)
	lines(dgamma(1:study_len, vals$shape, vals$rate), col = "red", lwd=2)
	legend("topright",
		   legend = c(paste0("True Serial Interval (", sim_type,")"),
					  paste0("Estimated Serial Interval (gamma)")),
		   col = c("black", "red"), lty=1, cex=0.9)
}


######################################################################
##############         Simulate Incidence       ######################
######################################################################

infections_plot <- function(incid){
	infected <- incid$infective_day
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