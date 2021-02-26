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
outpath <- paste0(getwd(), '/Output/')

# import data and focus on the US
actual_data <- read.csv2(paste0(imppath, 'EuropeCovidData.csv'))
actual_data <- subset(actual_data, Country.Region == 'US')

# drop first 80 observations (2M) to account for testing catching up
# actual_data <- tail(actual_data, -80)

mean <- 6.6
std <- 1.1
R <- c(1.8, 1.6, 1.4, 1.2)

window <- c(8, 11, 14)

experimentR <- function(mean, std, R, window, data = actual_data){
	data <- data

	for (w in window){
		gamma_y <- dgamma(gamma_x, mean, std)
		gamma <- rev(gamma_y)

			for (R in R){
				data$R_val <- R
				data$infective_new <- 0
				data$infective_new[1:window] <- 10
				for (t in (window + 1):length(data$date)) {
					x <- data$R_val[t] * sum(data[(t - window):t, ]$infective_new * gamma)
					data$infective_new[t] <- rpois(1, x)
				}
				tsi_est <- estimate_R(data$infective_new, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))
				jpeg(paste0(outpath, "TSISim_", 'Window' , w, '_R', R,".jpeg"))
				par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
				plot(data$infective_new, type='l', xlab='Number of Days', ylab = 'Daily Infections')
				plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
				lines(data$R_val, col = "red")
				legend(100, 1,
					   legend = c("Estimated Rt", "Simulated Rt"),
					   fill = c("black", "red"))
				dev.off()



}
}}

experimentR(mean, std, R, window)

# new Nour simulation
# as we are dealing in deltas we need a data frame that works in deltas (say hours for now)


tau_m <-  15
R <- c(1.6, 0.9, 1.3)
days <- length(R)*40
delta <- 1/24
N <- tau_m / delta

nour_sim <- function(days = days, tau_m = tau_m, R = R, delta, noisy = F){
	# set up data
	data <- data.frame(matrix(nrow = days*24, ncol = 4))
	colnames(data) <- c('t', 'days', 'R_t', 'infective')
	data$t <- 1:dim(data)[1]
	data$days <- rep(1:days, each = 24)
	data$infective[1:(tau_m*24)] <- round(seq(10, 100, length.out = tau_m*24))
	data$R_t <- rep(R, each = 40)

	if (noisy){
		data$R_t <- data$R_t + rnorm(1, 0.05, 0.025)
	}

	# set up gamma
	alpha <- 4.7
	beta <- 20
	mean <- alpha * beta
	std <- sqrt(alpha*beta^2)

	range <- seq(0, tau_m*24-1, 1)
	gamma_y <- dgamma(range, alpha, rate = 1/beta)
	gamma <- rev(gamma_y)

	# simulate outbreak
	start <- ((tau_m) * 24) + 1

	for (t in start:dim(data)[1]) {
		I_vec <- data[data$t %in% (t-N):(t-1),]$infective
		R_mean <- mean(data[which(data$t==t),]$R_t)
		total_infec <- sum(I_vec * gamma)

		infec <- R_mean * total_infec
		data[which(data$t == t),]$infective <- infec
	}

	# aggregate to daily (avg = /delta)
	daily_infec <- data %>%
	group_by(days) %>%
	summarise(infective_day = round(mean(infective)), R_val = mean(R_t))

	# estimate and export plot
	mean <- mean/24
	std <- std/24
	tsi_est <- estimate_R(daily_infec$infective_day, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))
	jpeg(paste0(outpath, "TSISim_Nour.jpeg"))
	par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
	plot(daily_infec$infective_day, type='l', xlab='Number of Days', ylab = 'Daily Infections')
	plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
	lines(daily_infec$R_val, col = "red")
	legend(50, 1,
	   legend = c("Estimated Rt", "Simulated Rt"),
	   fill = c("black", "red"))
	dev.off()
}






