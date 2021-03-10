# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(fitdistrplus)

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

window <-  11 # simulation window
R_val <- c(1.6, 0.9, 1.3)
rep <- 60
n_days <- length(R_val)*rep
delta <- 1/24
n <- window / delta
tau_est <- 10 # estimation window


nour_sim <- function(days = n_days, tau_m = window, R = R_val, N = n, plotname, noisy = F, conf = T){
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
	# andy:
	# mean <- 94
	# std <- 43.35897
	# alpha <- (mean^2) / (std^2)
	# beta <- mean / (std^2)

	# jp random:
	#alpha <- 4.7
	#beta <- 20
	#mean <- alpha * beta
	#std <- sqrt(alpha*beta^2)

	# jp new:
	mean <- 6.6*24
	std <- 1.1*24
	beta <- std^2/mean
	alpha <- mean/beta

	print(alpha)
	print(beta)

	range <- seq(0, ((tau_m*24)-1), 1)
	gamma_y <- dgamma(range, shape=alpha, rate=1/beta)
	gamma <- rev(gamma_y)

	print(sum(gamma))

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
	mean <- mean*delta
	std <- std*delta
	tsi_est <- estimate_R(daily_infec$infective_day, method = 'parametric_si', config = make_config(list(mean_si = mean, std_si = std)))$R
	tsi_est['lower'] <- tsi_est$`Mean(R)` - (1.96 * tsi_est$`Std(R)`)
	tsi_est['upper'] <- tsi_est$`Mean(R)` + (1.96 * tsi_est$`Std(R)`)
	newx <- 1:(days-7) # substract weekly window

	if (conf){
		jpeg(paste0(outpath, 'TSISim_Nour_', plotname, ".jpeg"))
		par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
		plot(daily_infec$infective_day, type='l', xlab='Number of Days', ylab = 'Daily Infections')
		plot(tsi_est$`Mean(R)`, type = 'n', ylim = c(0,3), xlab='Number of Days', ylab='Rt')
		polygon(x = c(rev(newx), newx), y = c(rev(tsi_est[, 12]), tsi_est[, 13]), col = 'grey80', border = NA)
		adjustcolor("grey80", alpha.f=0.05)
		lines(tsi_est$`Mean(R)`)
		lines(newx, tsi_est[, 12], lty = 'dashed', col = 'blue')
		lines(newx, tsi_est[, 13], lty = 'dashed', col = 'blue')
		lines(tail(daily_infec, -window)$R_val, col = "red")
		legend(50, 1,
			   legend = c("Estimated Rt", "Simulated Rt"),
			   fill = c("black", "red"))
		dev.off()
	}
	else {
		jpeg(paste0(outpath, 'TSISim_Nour_', plotname, ".jpeg"))
		par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
		plot(daily_infec$infective_day, type='l', xlab='Number of Days', ylab = 'Daily Infections')
		plot(tsi_est$`Mean(R)`, type ="l", ylim = c(0,3), xlab='Number of Days', ylab='Rt')
		lines(tail(daily_infec, -window)$R_val, col = "red")
		legend(50, 1,
			   legend = c("Estimated Rt", "Simulated Rt"),
			   fill = c("black", "red"))
		dev.off()
	}
	return(daily_infec)

}

nour_sim(R = 1.3, plotname = 'ConstantR', conf = F)

nour_sim(R = 1.3, noisy = T, plotname = 'NoisyR', conf = F)

nour_sim(plotname = 'StepR', conf = F)

nour_sim(noisy = T, plotname = 'StepR_noisy', conf = F)


nour_sim(R = 1.3, plotname = 'ConstantR_CI')

nour_sim(R = 1.3, noisy = T, plotname = 'NoisyR_CI')

nour_sim(plotname = 'StepR_CI')

nour_sim(noisy = T, plotname = 'StepR_noisy_CI')



# new estimate R function

discr_si <- function(k, mu, sigma) {

	k <- 1:(k)


  if (sigma < 0) {
    stop("sigma must be >=0.")
  }
  if (mu <= 1) {
    stop("mu must be >1")
  }
  if (any(k < 0)) {
    stop("all values in k must be >=0.")
  }

  a <- ((mu - 1) / sigma)^2
  b <- sigma^2 / (mu - 1)

  cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

  res <- k * cdf_gamma(k, a, b) + 
    (k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
  res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) - 
                          cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
  res <- sapply(res, function(e) max(0, e))

  return(res)
}


nour_sim_manual <- function(days = n_days, tau_m = window, tau_e = tau_est, R = R_val, N = n, plotname, noisy = F) {

	data <- data.frame(matrix(nrow = days*24, ncol = 4))
	colnames(data) <- c('t', 'days', 'R_t', 'infective')
	data$t <- 1:dim(data)[1]
	data$days <- rep(1:days, each = 24)
	data$infective[1:(tau_m*24)] <- round(seq(10, 100, length.out = tau_m*24))
	data$R_t <- rep(R, each = rep*24)

	if (noisy) {
		noise <- rnorm(dim(data)[1], 0.15, 0.1)
		data$R_t <- data$R_t + noise
	}

	# jp new:
	mean <- 6.6*24
	std <- 1.1*24
	beta <- std^2/mean
	alpha <- mean/beta

	print(alpha)
	print(beta)

	range <- seq(0, ((tau_m*24)-1), 1)
	gamma_y <- dgamma(range, shape=alpha, rate=1/beta)
	gamma <- rev(gamma_y)

	print(sum(gamma))

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


	# specify daily gamma again (overwrite okay?
	mean <- 6.6
	std <- 1.1
	beta <- std^2/mean
	alpha <- mean/beta

	# range <- seq(0, ((tau_m)-1), 1)
	# gamma_y <- dgamma(range, alpha, rate=1/beta)
	# gamma <- rev(gamma_y)

	start <- 2

	R_t <- list()
	for (t in start:dim(daily_infec)[1]) {
		gamma <- rev(discr_si(t, mean, std))
		I <- daily_infec[which(daily_infec$days==t),]$infective_day
		I_window <- daily_infec[daily_infec$days %in% 1:(t-1),]$infective_day
		R <- (I)/(sum(I_window * gamma))
		# daily_infec[which(daily_infec$days == t),]$R_est <- R_t
		R_t[t] <- R
	}

	R_t <- unlist(R_t, use.names=FALSE)


	jpeg(paste0(outpath, 'TSISim_Nour_Manual_', plotname, ".jpeg"))
	par(mfrow=c(2,1), tcl=0.5, family="serif", mai=c(1,1,0.3,0.3))
	plot(daily_infec$infective_day, type='l', xlab='Number of Days', ylab = 'Daily Infections')
	plot(R_t, type ="l", ylim = c(0,(R+1)), xlab='Number of Days', ylab='Rt')
	lines(tail(daily_infec, -tau_e)$R_val, col = "red")
	legend(50, 1,
		   legend = c("Estimated Rt", "Simulated Rt"),
		   fill = c("black", "red"))
	dev.off()
	return(daily_infec)

}


mean <- 6.6
std <- 1.1
beta <- std^2/mean
alpha <- mean/beta

T <- 25
data <- data.frame(rep(1.4, T))
names(data) <- "Rt"
data$GamMean <- rep(6.6, T)
data$GamStd <- rep(1.1, T)
data$Lambda <- rep(0, T)

sims <- 500
samples <- c()

#Take random sample from infector-infectee

samp_pois <- function(data) {

	range <- 1:nrow(data)

	#Generate Mean Incidence based on gamma prior
	for (t in range) {

		gamma <- discr_si(t, data[t, ]$GamMean, data[t, ]$GamStd)
		data[t, ]$Lambda <- data[t, ]$Rt * gamma[t]

		#Add up all the random poisson infections
		for (sim in 1:sims) {
			samples <- c(samples, rep(t, rpois(1, data[t, ]$Lambda)))
		}
	}
	return(samples)
}

samps <- samp_pois(data)

# pdf <- rowSums(mat) + 0.01
fit.gamma <- fitdist(samps, distr = "gamma")
fit.weibull <- fitdist(samps, distr = "weibull")
fit.norm <- fitdist(samps, distr = 'norm')
fit.lnorm <- fitdist(samps, distr = 'lnorm')

plot(dweibull(0:25, fit.weibull$estimate[1], fit.weibull$estimate[2]), type='l', col='red', ylim=c(0, 0.5))
lines(dgamma(0:25, fit.gamma$estimate[1], fit.gamma$estimate[2]), col='blue')
lines(dgamma(0:25, alpha, rate=1/beta), col='green')

lower_mean <- fit.gamma$estimate[1] - (1.96 * fit.gamma$sd[1])
upper_mean <- fit.gamma$estimate[1] + (1.96 * fit.gamma$sd[1])

lower_sd <- fit.gamma$estimate[2] - (1.96 * fit.gamma$sd[2])
upper_sd <- fit.gamma$estimate[2] + (1.96 * fit.gamma$sd[2])

gamma_vals <- expand.grid(seq(lower_mean, upper_mean, ((upper_mean - lower_mean)/5)), seq(lower_sd, upper_sd, ((upper_sd - lower_sd)/5)))

lower_mean <- fit.weibull$estimate[1] - (1.96 * fit.weibull$sd[1])
upper_mean <- fit.weibull$estimate[1] + (1.96 * fit.weibull$sd[1])

lower_sd <- fit.weibull$estimate[2] - (1.96 * fit.weibull$sd[2])
upper_sd <- fit.weibull$estimate[2] + (1.96 * fit.weibull$sd[2])

weibull_vals <- expand.grid(seq(lower_mean, upper_mean, ((upper_mean - lower_mean)/5)), seq(lower_sd, upper_sd, ((upper_sd - lower_sd)/5)))

################################################################################
#Simulating Serial Interval
#Simulating Poisson Process at Individual Level, with prior
################################################################################

nour_sim_manual(R = 4, plotname = 'ConstantR')

nour_sim_manual(plotname = 'StepR')

nour_sim_manual(tau_e = 8, plotname = 'StepR_noisy_Tau8')







