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

days <- 100
tau_m <-  11
R <- 1.6
delta <- 1/24
mean <- 6.6
std <- 1.1
N <- tau_m / delta

# construct data
data <- data.frame(matrix(nrow = days*24, ncol = 4))
colnames(data) <- c('t', 'days', 'R_t', 'infective')
data$t <- 1:dim(data)[1]
data$days <- rep(1:days, each = 24)
data$infective[1:(tau_m*24)] <- 100
data$R_t <- R

# set up gamma
# use delta as input for gamma to have more points
gamma_x <- seq(0, tau_m-1, 1)
gamma_y <- rep(dgamma(gamma_x, mean, std), each=1/delta)
gamma <- rev(gamma_y)

# without any knowledge about R_t usually alpha = 1, beta = 5
# https://www.pluralsight.com/guides/beta-and-gamma-function-implementation-in-r
alpha <- 3
beta <- 20
mean <- alpha * beta
std <- sqrt(alpha*beta^2)

range <- seq(0, tau_m*24-1, 1)
gamma_y <- dgamma(range, alpha, rate = 1/beta)
plot(range, gamma_y, type ="l")
gamma <- rev(gamma_y)

# reconstruct from mean = 4.46, std = 2.63
beta <- (2.63^2)/4.46
alpha <- 4.46/beta

alpha <- alpha
beta <- beta

# just to compare with Austrian estimation, looks very similar
# https://www.ages.at/download/0/0/068cb5fb9f2256d267e1a3dc8d464623760fcc30/fileadmin/AGES2015/Wissen-Aktuell/COVID19/Sch%C3%A4tzung_des_seriellen_Intervalles_von_COVID19_2020-04-08.pdf
range <- seq(0, tau_m+4, 1)
gamma_og <- dgamma(range, alpha, rate = 1/beta)
plot(range, gamma_og, type ="l")


start <- ((tau_m) * 24) + 1

for (t in start:dim(data)[1]) {
	I_vec <- data[data$t %in% (t-N):(t-1),]$infective
	R_mean <- mean(data[which(data$t==t),]$R_t)
	total_infec <- sum(I_vec * gamma)

	infec <- R_mean * total_infec
	data[which(data$t == t),]$infective <- infec
}

I_vec <- data[data$t %in% (265-N):(265-1),]$infective
R_mean <- mean(data[which(data$t==265),]$R_t)
total_infec <- sum(I_vec * gamma)