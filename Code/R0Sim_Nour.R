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

window <-  11
R_val <- c(1.6, 0.9, 1.3)
rep <- 40
n_days <- length(R_val)*rep
delta <- 1/24
n <- window / delta



nour_sim <- function(days = n_days, tau_m = window, R = R_val, N = n, plotname, noisy = F){
	# set up data
	# set up data
	data <- data.frame(matrix(nrow = days*24, ncol = 4))
	colnames(data) <- c('t', 'days', 'R_t', 'infective')
	data$t <- 1:dim(data)[1]
	data$days <- rep(1:days, each = 24)
	data$infective[1:(tau_m*24)] <- round(seq(10, 100, length.out = tau_m*24))
	data$R_t <- rep(R, each = rep*24)

	if (noisy){
		noise <- rnorm(dim(data)[1], 0.15, 0.1)
		data$R_t <- data$R_t + noise
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
	jpeg(paste0(outpath, 'TSISim_Nour_', plotname, ".jpeg"))
	par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
	plot(daily_infec$infective_day, type='l', xlab='Number of Days', ylab = 'Daily Infections')
	plot(tsi_est$R$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
	lines(tail(daily_infec, -window)$R_val, col = "red")
	legend(50, 1,
	   legend = c("Estimated Rt", "Simulated Rt"),
	   fill = c("black", "red"))
	dev.off()
}

nour_sim(R = 1.3, plotname = 'ConstantR')

nour_sim(R = 1.3, noisy = T, plotname = 'NoisyR')

nour_sim(plotname = 'StepR')

nour_sim(noisy = T, plotname = 'StepR_noisy')





