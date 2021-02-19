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
		gamma_x <- seq(0, (w-1), 1)
		gamma_y <- dgamma(gamma_x, mean, std)
		gamma <- rev(gamma_y)

			for (R in R){
				data$R_val <- R
				data$infective_new <- 0
				data$infective_new[1:window] <- 10
				for (t in (window + 1):length(data$date)) {
					x <- data$R_val[t] * sum(data[(t - window + 1):t, ]$infective_new * gamma)
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

# 100 days --> 2400 rows


days <- 100
w <-  11
R <- 1.4

# construct data
data <- data.frame(matrix(nrow = days*24, ncol = 3))
colnames(data) <- c('t', 'R_t', 'infective')
data$t <- rep(1:days, each = 24)
data$infective[1:(w*24)] <- 10
data$R_t <- R

# set up gamma
# use delta as input for gamma to have more points
gamma_x <- seq(0, (w-1), 1) # this should be 10*24 as our guy is now 240 hours infective (increment: 1/24)
gamma_y <- dgamma(gamma_x, mean, std)
gamma <- rev(gamma_y)

for (t in (w + 1):days) {
	infec <- mean(data[which(data$t==t),]$R_t) * 1/24 * sum(data[which(data$t %in% (t-w:t-1)),]$infective * gamma)
	print(infec)
	data[which(data$t == t),]$infective <- infec
}