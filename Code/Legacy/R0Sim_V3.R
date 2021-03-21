# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(accelerometry)
library(tidyverse)

imppath <- paste0(getwd(), '/Data/')

# import data and focus on the US
actual_data <- read.csv2(paste0(imppath, 'EuropeCovidData.csv'))
actual_data <- subset(actual_data, Country.Region == 'US')

# drop first 80 observations (2M) to account for testing catching up
# actual_data <- tail(actual_data, -80)

mean <- 6.6
std <- 1.1
R_steps <- c(1.5, 1.1, 1.3)
R_const <- 1.3

window <- 11

# define serial interval distribution taken from Cereda et al (2020)
gamma_x <- seq(0, (window-1), 1)
gamma_y <- dgamma(gamma_x, mean, std)
gamma <- rev(gamma_y)
plot(gamma_y, type = "l")

# time varying R_t

for (t in 1:length(actual_data$date)){
	if (actual_data$date[t] <= as.Date("2020-05-11")){
		actual_data$R_val_step[t] <- R_steps[1]
}
	else if (actual_data$date[t] > as.Date("2020-06-11") & actual_data$date[t] < as.Date("2020-08-12")) {
		actual_data$R_val_step[t] <- R_steps[2]
	}
		else if (actual_data$date[t] > as.Date("2020-08-11")) {
		actual_data$R_val_step[t] <- R_steps[3]
	}
}

# noisy time varying R_t

for (t in 1:length(actual_data$date)){
	if (actual_data$date[t] <= as.Date("2020-06-11")){
		actual_data$R_val_noisystep[t] <- R_steps[1] + rnorm(1, 0.05, 0.025)
}
	else if (actual_data$date[t] > as.Date("2020-06-11") & actual_data$date[t] < as.Date("2020-08-12")) {
		actual_data$R_val_noisystep[t] <- R_steps[2] + rnorm(1, 0.05, 0.025)
	}
		else if (actual_data$date[t] > as.Date("2020-08-11")) {
		actual_data$R_val_noisystep[t] <- R_steps[3] + rnorm(1, 0.05, 0.025)
	}
}


# generate outbreak according to poisson dist with mean R_t sum(I_t-s w_s)


### CONSTANT Rt
actual_data$infective_new <- 0
# actual_data[1:window, ]$infective_new <- round(seq(1, (10 - (9/window)), (9/window)))
actual_data$infective_new[1:window] <- 10

actual_data$R_val <- R_const

for (t in (window + 1):length(actual_data$date)) {
	x <- actual_data$R_val[t] * sum(actual_data[(t - window + 1):t, ]$infective_new * gamma)

	actual_data$infective_new[t] <- rpois(1, x)
}

#Estimate R via the function in the EpiEstim package
tsi_est <- estimate_R(actual_data$infective_new, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))

jpeg("../../TSISim_constant.jpeg")
par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
plot(actual_data$infective_new, type='l', xlab='Number of Days', ylab = 'Daily Infections')
plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
lines(actual_data$R_val, col = "red")
legend(100, 1, 
	legend = c("Estimated Rt", "Simulated Rt"),
	fill = c("black", "red"))
dev.off()

### STEP Rt
actual_data$infective_new <- 0
# actual_data[1:window, ]$infective_new <- round(seq(1, (10 - (9/window)), (9/window)))
actual_data$infective_new[1:window] <- 10

for (t in (window + 1):length(actual_data$date)) {
	x <- actual_data$R_val_step[t] * sum(actual_data[(t - window + 1):t, ]$infective_new * gamma)

	actual_data$infective_new[t] <- rpois(1, x)
}

#Estimate R via the function in the EpiEstim package
tsi_est <- estimate_R(actual_data$infective_new, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))

jpeg("../../TSISim_step.jpeg")
par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
plot(actual_data$infective_new, type='l', xlab='Number of Days', ylab = 'Daily Infections')
plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
lines(tail(actual_data, -8)$R_val_step, col = "red")
legend(100, 1,
	legend = c("Estimated Rt", "Simulated Rt"),
	fill = c("black", "red"))
dev.off()


### NOISY Step Rt
actual_data$infective_new <- 0
# actual_data[1:window, ]$infective_new <- round(seq(1, (10 - (9/window)), (9/window)))
actual_data$infective_new[1:window] <- 10

for (t in (window + 1):length(actual_data$date)) {
	x <- actual_data$R_val_noisystep[t] * sum(actual_data[(t - window + 1):t, ]$infective_new * gamma)

	actual_data$infective_new[t] <- rpois(1, x)
}

#Estimate R via the function in the EpiEstim package
tsi_est <- estimate_R(actual_data$infective_new, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))

jpeg("../../TSISim_noisystep.jpeg")
par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
plot(actual_data$infective_new, type='l', xlab='Number of Days', ylab = 'Daily Infections')
plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
lines(tail(actual_data, -8)$R_val_noisystep, col = "red")
legend(100, 1,
	legend = c("Estimated Rt", "Simulated Rt"),
	fill = c("black", "red"))
dev.off()