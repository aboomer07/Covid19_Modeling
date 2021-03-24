# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

source(paste0(getwd(), "/Code/SimFunc.R"))

######################################################################
####################### Set parameters ###############################
######################################################################

window <-  11 # simulation window
R_val <- c(1.6, 0.9, 1.3) # incidence R
n_days <- 180
delta <- 24
n <- window*delta

# parameters of the simulated distribution
sim_mean <- 7
sim_var <- 2
study_len <- 20
num_people <- 500
set.seed(14152118) #It spells Nour in numbers hihihi

######################################################################
############## Simulate Serial Interval Data #########################
######################################################################

samps <- samp_pois(R_val = 1.4, study_len = study_len, num_people = num_people,
				   sim_mu = sim_mean, sim_sig = sim_var, sim_type ='weibull', delta = delta) # check how to work with delta here

######################################################################
############## Estimate Serial Interval ##############################
######################################################################

vals <- serial_ests(samps) # here we obtain the params for Rt_est
vals

#estimates are sensitive to num_people, as should be! 

######################################################################
############## Simulate Incidence ####################################
######################################################################

incid <- nour_sim_data(sim_mean, sim_var, 'weibull', delta) # important! must be same dist as in samp_pois

######################################################################
##################### Estimate Rt ####################################
######################################################################

Rt <- Rt_est(incid, vals, 'gamma')
MSE <- MSE_est(Rt)

################################################################################
#Plot the estimated R's vs. the Simulated R
################################################################################

Rt_plot <- Rt %>%
	group_by(est_a, est_b) %>%
	mutate(id = cur_group_id()) %>%
	group_by(Date) %>%
	mutate(R_max = max(Est_Rt),
	       R_min = min(Est_Rt))

ggplot() +
	geom_line(data = Rt_plot, aes(x = Date, y=Est_Rt, group=id), alpha=0.3) +
	geom_line(data = Rt_plot, aes(x=Date, y=Rt), color = "red", size=1) +
	scale_color_brewer(palette="Dark2") +
	theme_minimal()


ggplot() +
	geom_ribbon(data = Rt_plot, 
		aes(x=Date, ymin = R_min, ymax = R_max), alpha=0.6) +
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


