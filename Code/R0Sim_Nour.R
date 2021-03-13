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
library(data.table)
library(np)
library(KernSmooth)

imppath <- paste0(getwd(), '/Data/')
outpath <- paste0(getwd(), '/Output/')

window <-  11 # simulation window
R_val <- c(1.6, 0.9, 1.3)
# rep <- 60
n_days <- 180
delta <- 1/24
n <- window / delta
# tau_est <- 10 # estimation window
sim_b <- 1.1^2/(6.6)
sim_a <- (6.6)/sim_b

#Function to generate model parameters correctly

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

	for (i in 1:length(dists)) {
		dist <- dists[i]
		fit <- fitdist(samps, distr = dist)

		a_est <- fit$estimate[1]
		b_est <- fit$estimate[2]

		a_sd <- fit$sd[1]
		b_sd <- fit$sd[2]

		if (dist == 'gamma') {
			mean <- a_est / b_est
			var <- a_est * (b_est^2)
			mean <- mean - 1
			b_est <- mean / sqrt(var)
			a_est <- mean / b_est
		}

		lower_a <- a_est - (1.96 * a_sd)
		upper_a <- a_est + (1.96 * a_sd)

		lower_b <- b_est - (1.96 * b_sd)
		upper_b <- b_est + (1.96 * b_sd)

		vals <- expand.grid(seq(lower_a, upper_a, ((upper_a - lower_a)/num_vars)),
			seq(lower_b, upper_b, ((upper_b - lower_b)/num_vars)))

		names(vals) <- c("Param_a", 'Param_b')
		de <- data.frame(a_est, b_est)
		names(de) <- c("Param_a", 'Param_b')
		vals <- rbind(vals, de)
		vals$Distribution <- dist
		vals$True_a <- a_est
		vals$True_b <- b_est

		params[[i]] <- vals

	}

	params <- do.call('rbind', params)

	return(params)
}

# estimate serial interval nonparametically
serial_ests_nonpara <- function(samps, correction, range){
	h_nsr <- 1.059*sd(samps)*length(samps)^(-1/5)
	kernel_est <- bkde(samps, bandwidth = h_nsr + correction, range.x = range, gridsize = max(range))
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

Rt_est <- function(df, vals, distribution){
	start <- 10
	data <- data.frame(matrix(nrow = n_days * nrow(vals), ncol = 6))
	names(data) <- c('Date', 'est_a', 'est_b', 'est_type', 'Rt', 'Est_Rt')

	data$Date <- rep(1:n_days, nrow(vals))
	data$est_a <- rep(vals[, 1], each=n_days)
	data$est_b <- rep(vals[, 2], each=n_days)
	data$est_type <- rep(vals[, 3], each=n_days)
	data$Rt <- rep(df['R_val'][[1]], nrow(vals))
	data$True_est_a <- rep(vals[, 4], each=n_days)
	data$True_est_b <- rep(vals[, 5], each=n_days)

	for (i in start:nrow(data)){
		if (data[i, ]$Date <= start) {
			data[i, ]$Est_Rt <- NA
		}
		else {
			a <- data[i, ]$est_a
			b <- data[i, ]$est_b
			type <- data[i, ]$est_type
			t <- data[i, ]$Date

			dist <- rev(gen_distribution(t-1, a, b, type))
			I <- df[which(df$days==t),]$infective_day
			I_window <- df[df$days %in% 1:(t-1),]$infective_day
			data[i, ]$Est_Rt <- (I)/(sum(I_window * dist))
		}
	}

	return(data) 
}

Rt_est_nonpara <- function (df, samps, corrections){
	start <- 10
	data <- data.frame(matrix(nrow = n_days * length(corrections), ncol = 4))
	names(data) <- c('Date', 'Cor_Par', 'Rt', 'Est_Rt')

	data$Date <- rep(1:n_days, length(corrections))
	data$Cor_Par <- rep(corrections, each=n_days) # instead of distribution parameters we can correct/opt bandwidth
	data$Rt <- rep(df['R_val'][[1]], length(corrections))

	for (i in start:nrow(data)){
		if (data[i, ]$Date <= start) {
			data[i, ]$Est_Rt <- NA
		}
		else {
			correction <- data[i, ]$Cor_Par
			t <- data[i, ]$Date

			dist <- rev(serial_ests_nonpara(samps, correction, range = c(1,(t-1)))$y)
			I <- df[which(df$days==t),]$infective_day
			I_window <- df[df$days %in% 1:(t-1),]$infective_day
			data[i, ]$Est_Rt <- (I)/(sum(I_window * dist))
		}
	}
	return(data)
}

MSE_est <- function(df){

	df$SSE <- (df$Rt - df$Est_Rt)^2
	df <- df %>% group_by(est_a, est_b, est_type) %>% summarize(MSE=mean(SSE, na.rm=TRUE), True_est_a=mean(True_est_a), True_est_b=mean(True_est_b))

	return(df)
}

samps <- samp_pois(c(1.4), 25, 2000, sim_a, sim_b, 'gamma')

est_dists <- c("gamma", 'weibull', 'norm', 'lnorm')

dat <- nour_sim_data(sim_a, sim_b, 'gamma')

params <- serial_ests(samps, 5, est_dists)

Rt_mat <- list()
MSE_mat <- list()

for (i in 1:length(est_dists)) {
	dist <- est_dists[i]
	dist_vals <- params[params['Distribution'] == dist, ]

	mat <- Rt_est(dat, dist_vals, dist)
	mat_MSE <- MSE_est(mat)

	Rt_mat[[i]] <- mat
	MSE_mat[[i]] <- mat_MSE
}

Rt_mat <- do.call('rbind', Rt_mat)
MSE_mat <- do.call('rbind', MSE_mat)

################################################################################
#Plot the different estimated distributions vs. original gamma
################################################################################

og_gamma <- gen_distribution(15, sim_a, sim_b, "gamma")
new_gamma <- gen_distribution(15, mean(params[params$Distribution == 'gamma', 'True_a']), mean(params[params$Distribution == 'gamma', 'True_b']),"gamma")
new_weibull <- gen_distribution(15, mean(params[params$Distribution == 'weibull', 'True_a']), mean(params[params$Distribution == 'weibull', 'True_b']),"weibull")
new_norm <- gen_distribution(15, mean(params[params$Distribution == 'norm', 'True_a']), mean(params[params$Distribution == 'norm', 'True_b']),"norm")
new_lnorm <- gen_distribution(15, mean(params[params$Distribution == 'lnorm', 'True_a']), mean(params[params$Distribution == 'lnorm', 'True_b']),"lnorm")

jpeg("DistCompare.jpg")
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
dev.off()

################################################################################
#Plot the estimated R's vs. the Simulated R
################################################################################

Rt_plot <- Rt_mat %>%
	group_by(est_a, est_b, est_type) %>%
	mutate(id = cur_group_id()) %>%
	group_by(est_type, Date) %>%
	mutate(R_max = max(Est_Rt),
	       R_min = min(Est_Rt))

ggplot() +
	geom_line(data = Rt_plot, aes(x = Date, y=Est_Rt, group=id, color=est_type), alpha=0.3) +
	geom_line(data = Rt_plot, aes(x=Date, y=Rt), color = "red", size=1) +
	scale_color_brewer(palette="Dark2") +
	theme_minimal()


ggplot() +
	geom_ribbon(data = Rt_plot, aes(x=Date, ymin = R_min, ymax = R_max, fill=est_type), alpha=0.6) +
	geom_line(data = Rt_plot, aes(x=Date, y=Rt), color = "red", size=1) +
	scale_fill_brewer(palette="Dark2") +
	theme_minimal() +
	ggsave("Rt_Est.png")

################################################################################
#Plot the MSE vs. the Mean and Standard Dev
################################################################################

test <- MSE_mat

test$dev_a <- (test$est_a - test$True_est_a) / test$True_est_a
test$dev_b <- (test$est_b - test$True_est_b) / test$True_est_b

test$b_group <-as.factor(test$dev_b)
test$a_group <-as.factor(test$dev_a)

testLine <- test %>%
  group_by(est_type) %>%
  summarize(True_est_a = mean(True_est_a), True_est_b = mean(True_est_b))


a_gam <- ggplot(data=test[test$est_type == 'gamma',], aes(x=dev_a, y=MSE, group=b_group, color=b_group)) +
	geom_line() + 
	scale_color_brewer(palette="Dark2") +
	ggsave("Gamma_MSE.png")

a_wei <- ggplot(data=test[test$est_type == 'weibull',], aes(x=dev_a, y=MSE, group=b_group, color=b_group)) +
	geom_line() + 
	scale_color_brewer(palette="Dark2") +
	ggsave("Weibull_MSE.png")

a_norm <- ggplot(data=test[test$est_type == 'norm',], aes(x=dev_a, y=MSE, group=b_group, color=b_group)) +
	geom_line() + 
	scale_color_brewer(palette="Dark2") +
	ggsave("Norm_MSE.png")

a_lnorm <- ggplot(data=test[test$est_type == 'lnorm',], aes(x=dev_a, y=MSE, group=b_group, color=b_group)) +
	geom_line() + 
	scale_color_brewer(palette="Dark2") +
	ggsave("LNorm_MSE.png")

est_b <- 
  # ggplot(data=test, aes(x=est_a, y=MSE)) +
  ggplot(data=test[test$est_type == 'gamma',], aes(x=dev_a, y=MSE, group=b_group, color=b_group)) +
  geom_line()+ 
  # facet_grid(~ est_type, scales='free')+ 
  # geom_vline(data=testLine[testLine$est_type == 'gamma', ], aes(xintercept = True_est_a), linetype="dotted") +
  scale_color_brewer(palette="Dark2")
est_b

est_a <- test %>%
  ggplot( aes(x=est_b, y=MSE, group=est_a, color=est_a)) +
  geom_line()+ 
  facet_grid(~ est_type)+ 
  geom_vline(data=testLine, aes(xintercept=est_b))+
  scale_color_brewer(palette="Dark2")
est_a


