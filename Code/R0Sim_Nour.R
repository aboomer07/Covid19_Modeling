# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

source("DistSens.R")

######################################################################
######## Set right parameters based on distribution and delta ########
######################################################################

window <-  11 # simulation window
R_val <- c(1.6, 0.9, 1.3)
# rep <- 60
n_days <- 180
delta <- 1/24
n <- window / delta
# tau_est <- 10 # estimation window

# parameters of the different distributions
sim_var <- 1.1
sim_mu <- 6.6

# weib<-weibullparinv(1.96, 8.47, loc = 0) #values from June Young Chun, Gyuseung Baek
a_weibull <- weib$mu
b_weibull <- weib$sigma

a_norm <- 6.6
b_norm <- 1.1

a_lnorm <- 1.5
b_lnorm <- 0.3

param_a <- c(sim_mu, a_weibull, a_norm, a_lnorm)
param_b <- c(sim_var, b_weibull, b_norm, b_lnorm)
dist <- c("gamma", 'weibull', 'norm', 'lnorm' )

mat_param<- cbind(dist, param_a, param_b)
mat_param <- as.data.frame(mat_param)
mat_param$param_a <-as.numeric(mat_param$param_a)
mat_param$param_b <-as.numeric(mat_param$param_b)

est_dists <- c("gamma", 'weibull', 'norm', 'lnorm')

dat <- nour_sim_data(sim_mu, sim_var, 'gamma', 24)

params <- serial_ests(est_mu, est_var, 10, est_dists)

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
#Plot the different estimated distributions vs. original distributions
################################################################################

og_gamma <- gen_distribution(15, sim_mu, sim_var, "gamma")
new_gamma <- gen_distribution(15, mean(params[params$Distribution == 'gamma', 'True_a']), mean(params[params$Distribution == 'gamma', 'True_b']),"gamma")
og_weibull <- gen_distribution(15, a_weibull, b_weibull, "weibull")
new_weibull<- gen_distribution(15, mean(params[params$Distribution == 'weibull', 'True_a']), mean(params[params$Distribution == 'weibull', 'True_b']),"gamma")
og_norm <- gen_distribution(15, a_norm, b_norm, "norm")
new_norm <- gen_distribution(15, mean(params[params$Distribution == 'norm', 'True_a']), mean(params[params$Distribution == 'norm', 'True_b']),"gamma")
og_lnorm <- gen_distribution(15, a_lnorm, b_lnorm, "lnorm")
new_lnorm <- gen_distribution(15, mean(params[params$Distribution == 'lnorm', 'True_a']), mean(params[params$Distribution == 'lnorm', 'True_b']),"gamma")

jpeg("DistCompare.jpg")
layout(matrix(1:4, nrow = 2, ncol=2))
plot(og_gamma, type="l")
lines(new_gamma, type="l", col="green")
title("Gamma")

plot(og_weibull, type="l")
lines(new_weibull, type="l", col="red")
title("Weibull")

plot(og_norm, type="l")
lines(new_norm, type="l", col="blue")
title("Normal")

plot(og_lnorm, type="l")
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


