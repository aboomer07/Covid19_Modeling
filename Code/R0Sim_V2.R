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

sir_sim <- read.csv2(paste0(imppath, 'SimulatedSIR.csv'))[-1, ]
sir_sim$Delta <- as.numeric(sir_sim$Delta)
sir_sim$S <- as.numeric(sir_sim$S)
sir_sim$Rt <- as.numeric(sir_sim$Rt)
sir_sim$I <- as.numeric(sir_sim$I)
sir_sim$Reff <- sir_sim$Rt * (sir_sim$S / 100000)

# drop first 80 observations (2M) to account for testing catching up
# actual_data <- tail(actual_data, -80)

mean <- 6.6
std <- 1.1
R0 <- 1.6

window <- 11

# define serial interval distribution taken from Cereda et al (2020)
gamma_x <- seq(1, window, 1)
gamma_y <- dgamma(gamma_x, mean, std)
gamma <- rev(gamma_y)
plot(gamma_y, type = "l")


# Define time varying effective reproduction number
# important: here you need to specify manually when you would like the true R0 to change

for (t in 1:length(actual_data$date)){
	if (actual_data$date[t] <= as.Date("2020-06-11")){
		actual_data$R_val[t] <- R0
}
	else if (actual_data$date[t] > as.Date("2020-06-11")) {
		actual_data$R_val[t] <- R0
	}
}

# cannot do it in one loop cause i suck
for (t in 1:length(actual_data$date)){
	if (actual_data$date[t] > as.Date("2020-09-11")) {
		actual_data$R_val[t] <- R0
	}
}

# time varying R_t

# set.seed(1)

# actual_data$R_val <- actual_data$R_val + cumsum(sample(c(-0.1, 0.1), length(actual_data$R_val), TRUE))

#Generate outbreak accoring to E[I_t]=R+t*sum^{t}{s=1}I{t-s}w_s
#we use a window of 5 days

actual_data$infective <- 0
actual_data[1:window, ]$infective <- round(seq(1, (10 - (9/window)), (9/window)))
actual_data$simulation <- 0

for (t in (window + 1):length(actual_data$date)) {

	x <- sum(actual_data[(t - window + 1):t, ]$infective * gamma)
	actual_data$infective[t] <- x * actual_data$R_val[t]
}

# generate outbreak according to poisson dist with mean R_t sum(I_t-s w_s)

actual_data$infective_new <- 0
# actual_data[1:window, ]$infective_new <- round(seq(1, (10 - (9/window)), (9/window)))
actual_data$infective_new[1:window] <- 10

actual_data$R_val <- 1.6

for (t in (window + 1):length(actual_data$date)) {
	x <- actual_data$R_val[t] * sum(actual_data[(t - window + 1):t, ]$infective_new * gamma)

	actual_data$infective_new[t] <- rpois(1, x)
}

#Estimate R via the function in the EpiEstim package
tsi_est <- estimate_R(actual_data$infective_new, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))

jpeg("TSISim.jpeg")
par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
plot(actual_data$infective_new, type='l', xlab='Number of Days', ylab = 'Daily Infections')
plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
lines(actual_data$R_val, col = "red")
legend(100, 1, 
	legend = c("Estimated Rt", "Simulated Rt"),
	fill = c("black", "red"))
dev.off()

#Estimate R via the function in the EpiEstim package
sir_est <- estimate_R(sir_sim$Delta, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))

jpeg("SIRSim.jpeg")
par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
plot(sir_sim$Delta, type='l', xlab='Number of Days', ylab = 'Daily Infections')
plot(sir_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
lines(sir_sim$Rt, col = "red")
lines(sir_sim$Reff, col = "blue")
legend(100, 3, 
	legend = c("Estimated Rt", "Simulated Rt", "Simulated Reff (Rt * (S/N))"),
	fill = c("black", "red", 'blue'))
dev.off()