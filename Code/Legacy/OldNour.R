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

# nour_sim_manual(R = 4, plotname = 'ConstantR')

# nour_sim_manual(plotname = 'StepR')

# nour_sim_manual(tau_e = 8, plotname = 'StepR_noisy_Tau8')