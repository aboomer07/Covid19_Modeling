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
library(RColorBrewer)

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

window <-  11 # simulation window
R_val <- c(1.6, 0.9, 1.3)
# rep <- 60
n_days <- 180
delta <- 1/24
n <- window / delta
# tau_est <- 10 # estimation window
sim_b <- 1.1^2/6.6
sim_a <- 6.6/sim_b

gen_distribution <- function(k, a, b, type) {
	k <- 1:k
	if (type == 'gamma') {
		cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

		res <- k * cdf_gamma(k, a, b) + 
		(k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
		res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) - 
		                      cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
		res <- sapply(res, function(e) max(0, e))

		return(res)
	}

	if (type == 'weibull') {
		res <- dweibull(k, a, b)
		return(res)
	}

	if (type == 'norm') {
		res <- dnorm(k, a, b)
		return(res)
	}

	if (type == 'lnorm') {
		res <- dlnorm(k, a, b)
		return(res)
	}

}

# Do a serial interval "study"
samp_pois <- function(R_val, study_len, num_people, sim_a, sim_b, sim_type) {

	data <- data.frame(rep(R_val, (study_len/length(R_val))))
	names(data) <- "Rt"
	data$sim_a <- rep(sim_a, study_len)
	data$sim_b <- rep(sim_b, study_len)
	data$Lambda <- rep(0, study_len)

	sims <- num_people
	samples <- c()

	range <- 1:nrow(data)

	#Generate Mean Incidence based on gamma prior
	for (t in range) {

		dist <- gen_distribution(t, data[t, ]$sim_a, data[t, ]$sim_b, sim_type)
		data[t, ]$Lambda <- data[t, ]$Rt * dist[t]

		#Add up all the random poisson infections
		for (sim in 1:sims) {
			samples <- c(samples, rep(t, rpois(1, data[t, ]$Lambda)))
		}
	}
	return(samples)
}

serial_ests <- function(samps, num_vars, dists) {
	num_vars <- 5
	params <- list()
	dists <- c('gamma', 'weibull', 'norm', 'lnorm')

	for (i in 1:length(dists)) {
		dist <- dists[i]
		fit <- fitdist(samps, distr = dist)

		a_est <- fit$estimate[1]
		b_est <- fit$estimate[2]

		a_sd <- fit$sd[1]
		b_sd <- fit$sd[2]

		lower_a <- a_est - (1.96 * a_sd)
		upper_a <- a_est + (1.96 * a_sd)

		lower_b <- b_est - (1.96 * b_sd)
		upper_b <- b_est + (1.96 * b_sd)

		vals <- expand.grid(seq(lower_a, upper_a, ((upper_a - lower_a)/num_vars)),
			seq(lower_b, upper_b, ((upper_b - lower_b)/num_vars)))

		names(vals) <- c("Param_a", 'Param_b')
		vals$Distribution <- dist
		vals$True_a <- a_est
		vals$True_b <- b_est

		params[[i]] <- vals

	}

	params <- do.call('rbind', params)

	return(params)
}

nour_sim_data <- function(sim_a, sim_b, sim_type, days = n_days, tau_m = window, R = R_val, N = n, noisy = F) {

	data <- data.frame(matrix(nrow = days*24, ncol = 4))
	colnames(data) <- c('t', 'days', 'R_t', 'infective')
	data$t <- 1:dim(data)[1]
	data$days <- rep(1:days, each = 24)
	data$infective[1:(tau_m*24)] <- round(seq(10, 100, length.out = tau_m*24))
	data$R_t <- rep(R, each = (24 * days/length(R)))

	if (noisy) {
		noise <- rnorm(dim(data)[1], 0.15, 0.1)
		data$R_t <- data$R_t + noise
	}

	# jp new:
	sim_b <- sim_b / delta

	print(sim_a)
	print(sim_b)

	dist <- rev(gen_distribution(((tau_m*24)), sim_a, sim_b, sim_type))
	print(sum(dist))

	# simulate outbreak
	start <- ((tau_m) * 24) + 1

	for (t in start:dim(data)[1]) {
		I_vec <- data[data$t %in% (t-N):(t-1),]$infective
		R_mean <- mean(data[which(data$t==t),]$R_t)
		total_infec <- sum(I_vec * dist)

		infec <- R_mean * total_infec
		data[which(data$t == t),]$infective <- infec
	}

	# aggregate to daily (avg = /delta)
	daily_infec <- data %>%
	group_by(days) %>%
	summarise(infective_day = round(mean(infective)), R_val = mean(R_t))

	return(daily_infec)
}

Rt_est <- function(df, vals){
	start <- 2
	data <- data.frame(matrix(nrow = n_days * nrow(vals), ncol = 6))
	names(data) <- c('Date', 'est_a', 'est_b', 'est_type', 'Rt', 'Est_Rt')

	data$Date <- rep(1:n_days, nrow(vals))
	data$est_a <- rep(vals[, 1], each=n_days)
	data$est_b <- rep(vals[, 2], each=n_days)
	data$est_type <- rep(vals[, 3], each=n_days)
	data$Rt <- rep(df['R_val'][[1]], nrow(vals))

	for (i in 1:nrow(data)){ 
		a <- data[i, ]$est_a
		b <- data[i, ]$est_b
		type <- data[i, ]$est_type
		t <- data[i, ]$Date

		dist <- rev(gen_distribution(t-1, a, 1/b, type))
		I <- df[which(df$days==t),]$infective_day
		I_window <- df[df$days %in% 1:(t-1),]$infective_day
		# R <- (I)/(sum(I_window * dist))
		data[i, ]$Est_Rt <- (I)/(sum(I_window * dist))
		}

	return(data) 
}

MSE_est <- function(df){

	df$SSE <- sum((df['Rt'] - df['Est_Rt'])^2)

	MSE_mat <- df %>% group_by(est_a, est_b, est_type) %>%
		summarize(MSE = mean(SSE))

	# start <- 2
	# mat <- matrix(nrow = length(vals[,1]), ncol = 3)
	# type <- vals[1, ]$Distribution

	# for (i in seq_along(vals[,1])){
	# 	for (j in seq_along(vals[,2])){
	# 		mean <- vals[i, 1]/vals[j, 2]
	# 		std <- sqrt(mean*(1/vals[j, 2]))
	# 		R_t <- list()
	# 		for (t in start:dim(df)[1]) {
	# 			gamma <- rev(discr_si(t, mean, std))
	# 			I <- df[which(df$days==t),]$infective_day
	# 			I_window <- df[df$days %in% 1:(t-1),]$infective_day
	# 			R <- (I)/(sum(I_window * gamma))
	# 			# daily_infec[which(daily_infec$days == t),]$R_est <- R_t
	# 			R_t[t] <- R
	# }
	# 		R_t <- unlist(R_t, use.names=FALSE)

	# 		MSE <- mean(((tail(rep(R_val, each = rep),-1)) - R_t)^2)
	# 		mat[i, 1] <- vals[i, 1]
	# 		mat[j, 2] <- vals[j, 2]
	# 		mat[j, 3] <- MSE
	# 	}
	# }
	colnames(MSE_mat) <- c('est_a', 'est_b', 'est_type', 'MSE')
	return(MSE_mat)
}

samps <- samp_pois(c(1.4), 25, 500, sim_a, sim_b, 'gamma')

est_dists <- c("gamma", 'weibull', 'norm', 'lnorm')

params <- serial_ests(samps, 5, est_dists)

dat <- nour_sim_data(sim_a, sim_b, 'gamma')

dist_vals <- params[params['Distribution'] == 'gamma',]

mat <- Rt_est(dat, dist_vals)

MSE_est <- MSE_est(mat)

################################################################################
#Plot the different estimated distributions vs. original gamma
################################################################################

og_gamma <- gen_distribution(15, sim_a, sim_b, "gamma")
new_gamma <- gen_distribution(15, mean(params[params$Distribution == 'gamma', 'True_a']), 1/mean(params[params$Distribution == 'gamma', 'True_b']),"gamma")
new_weibull <- gen_distribution(15, mean(params[params$Distribution == 'weibull', 'True_a']), mean(params[params$Distribution == 'weibull', 'True_b']),"weibull")
new_norm <- gen_distribution(15, mean(params[params$Distribution == 'norm', 'True_a']), mean(params[params$Distribution == 'norm', 'True_b']),"norm")
new_lnorm <- gen_distribution(15, mean(params[params$Distribution == 'lnorm', 'True_a']), mean(params[params$Distribution == 'lnorm', 'True_b']),"lnorm")

layout(matrix(1:4, nrow = 2, ncol=2))
plot(og_gamma, type="l")
lines(new_gamma, type="l", col="green")
title("Gamma")

plot(og_gamma, type="l")
lines(new_weibull, type="l", col="red")
title("Weibull")

plot(og_gamma, type="l")
lines(new_norm, type="l", col="blue")
title("Normal")

plot(og_gamma, type="l")
lines(new_lnorm, type="l", col="orange")
title("Log-normal")

################################################################################
#Plot the estimated R's vs. the Simulated R
################################################################################

test <- expand.grid(seq(0, 10, 2), seq(0, 10, 2))
names(test) <- c("Mean", 'Std')
test$MSE <- seq(1, nrow(test), 1)
test <- rbind(test, test, test, test)
test$Distribution <- rep(c("Gamma", 'Weibull', 'Norm', 'Lnorm'), 
	each=nrow(test)/4)
test <- do.call("rbind", replicate(100, test, simplify = F))
test$Period <- rep(1:100, each=144)
test$R <- 1.5
test$R_est <- ifelse(test$Distribution == "Gamma", 1.6 + rnorm(3600, 0, 0.5),
	ifelse(test$Distribution == "Weibull", 1.4 + rnorm(3600, 0, 0.1),
	ifelse(test$Distribution == "Norm", 1.2 + rnorm(3600, 0, 0.4),
		1.3 + rnorm(3600, 0, 0.1))))

test <- test %>%
	group_by(Mean, Std, Distribution) %>%
	mutate(id = cur_group_id()) %>%
	group_by(Distribution, Period) %>%
	mutate(R_max = max(R_est),
	       R_min = min(R_est))

ggplot() +
	geom_line(data = test, aes(x = Period, y=R_est, group=id, color=Distribution), alpha=0.3) +
	geom_line(data = test, aes(x=Period, y=R), color = "red", size=1) +
	scale_color_brewer(palette="Dark2") +
	theme_minimal()


ggplot() +
	geom_ribbon(data = test, aes(x=Period, ymin = R_min, ymax = R_max, fill=Distribution), alpha=0.6) +
	geom_line(data = test, aes(x=Period, y=R), color = "red", size=1) +
	scale_fill_brewer(palette="Dark2") +
	theme_minimal()

################################################################################
#Plot the MSE vs. the Mean and Standard Dev
################################################################################

test$Std<-as.factor(test$Std)
test$Mean<-as.factor(test$Mean)

testLine <- test %>%
  group_by(Distribution) %>%
  summarize(mean = mean(True_Mean), Std=mean(True_Std))

Std<-test %>%
  ggplot( aes(x=Mean, y=MSE, group=Std, color=Std)) +
  geom_line()+ 
  facet_grid(~ Distribution)+ 
  geom_vline(data=testLine, aes(xintercept = mean), linetype="dotted" )+
  scale_color_brewer(palette="Dark2")
Std

Mean <- test %>%
  ggplot( aes(x=Std, y=MSE, group=Mean, color=Mean)) +
  geom_line()+ 
  facet_grid(~ Distribution)+ 
  geom_vline(data=testLine, aes(xintercept=Std))+
  scale_color_brewer(palette="Dark2")
Mean


