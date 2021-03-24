# Objective : Evaluate fitted distributions on serial interval data
# Created by: jacobpichelmann
# Created on: 21.03.21

source(paste0(getwd(), "/Code/SimFunc.R"))

### Simulate serial interval 
samp_pois_plot <- function(samps, study_len){ 
  data.frame(table(samps)) %>%
    mutate(samps = as.numeric(as.character(samps))) %>%
    ggplot(aes(x=samps, y=Freq)) + geom_bar(stat = "identity") + xlab("days") + ylab("count") +
      scale_x_continuous(limits = c(0,study_len)) + theme_minimal()
}

samp_pois_plot(samps = samps, study_len = study_len*delta)

### Fit gamma and compare it to true distribution -daily

dist_estimate <- function(study_len, sim_mean, sim_var, sim_type, vals){
	omegadaily <- gen_distribution(study_len, sim_mean, sim_var, sim_type, 1) 
	plot(omegadaily, type = "l", main = "Distribution of true generation time vs its estimate", xlab = "Time", ylab = "Density")
	lines(dgamma(1:study_len, vals$shape, vals$rate), col = "red")
}

dist_estimate(study_len, sim_mean, sim_var, "weibull", vals)


### Draw daily incidence of cases and corresponding Rt, given type of omega

infections_plot <- function(type){
	samps <- samp_pois(R_val = 1.4, study_len = study_len, num_people = num_people,
				   sim_mu = sim_mean, sim_sig = sim_var, sim_type = type, delta = delta)
	vals <- serial_ests(samps)
	incid <- nour_sim_data(sim_mean, sim_var, type, delta)

	infected = incid$infective_day
	Rt_sim = incid$R_val

	jpeg("IncidenceRt.jpg")
	layout(matrix(1:2, ncol = 1, nrow =2))
	plot(infected, type = "l", main = paste0("New infections per day (", type, ")"), xlab = "Day", ylab = "New infections")
	plot(Rt_sim, type = "l", main = "Simulated R(t)", xlab = "Day", ylab = "Rt", ylim = c(0, 2))
	dev.off()
}

infections_plot(type = 'weibull')


### Compare estimate of Rt to its true value 

compare_rt <- function(simulate_dist, estimate_dist){

	samps <- samp_pois(R_val = 1.4, study_len = study_len, num_people = num_people,
				   sim_mu = sim_mean, sim_sig = sim_var, sim_type = simulate_dist, delta = delta)
	vals <- serial_ests(samps)
	incid <- nour_sim_data(sim_mean, sim_var, simulate_dist, delta)
	Rt <- Rt_est(incid, vals, estimate_dist)
	
	jpeg("RtCompare.jpg")	
	plot(Rt$Est_Rt, type = "l", main = paste0("R(t) estimation (", simulate_dist, ")"), xlab = "Day", ylab = "R(t)")
	lines(Rt$Rt, col = "red")
	dev.off()
	}

compare_rt(simulate_dist = 'weibull', estimate_dist = 'gamma')

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