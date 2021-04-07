# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

imppath <- paste0(getwd(), '/Code/')
outpath <- paste0(getwd(), '/Output/')

source(paste0(getwd(), "/Code/Params.R"))
source(paste0(getwd(), "/Code/SimFunc.R"))
source(paste0(getwd(), "/Code/EvalDist.R"))
source(paste0(getwd(), "/Code/SIIFunc.R"))


set.seed(14152118) #It spells Nour in numbers hihihi

######################################################################
############## Simulate Serial Interval Data #########################
######################################################################

# simulate serial interval study
serinfect <- samp_pois(params)
# get discretized and "continuous" secondary cases
samps <- serinfect$daily
sampscont <- serinfect$samplescont
# dist <- serinfect$dist

# plot serial interval simulation in continuous time
pdf(file = paste0(outpath, "SerialHistCont_", params$sim_type, ".pdf"))
serial_hist_cont(sampscont, dist)
dev.off()

# plot discretized serial interval
pdf(file = paste0(outpath, "SerialHistDisc_", params$sim_type, ".pdf"))
serial_hist_disc(samps)
dev.off()

######################################################################
############# Simulate Distribution of Estimators ####################
######################################################################


SI.simulation <- params_distribution(params)
SI.plot <- SI_plot_distribution(data = SI.simulation)


######################################################################
############## Estimate Serial Interval ##############################
######################################################################

vals <- serial_ests(samps) # here we obtain the params for Rt_est

#estimates are sensitive to num_people, as should be!

# plot estimate
pdf(file = paste0(outpath, "SerialEst_", params$sim_type, ".pdf"))
serial_est_plot(params$study_len, params$sim_mean, params$sim_var, 
	params$sim_type, vals)
dev.off()



######################################################################
############## Simulate Incidence ####################################
######################################################################
# simulate outbreak
incid <- nour_sim_data(params) # important! must be same dist as in samp_pois

incid <- sii_sim(params)

# incid <- sii_sim(params)

# plot outbreak
pdf(file = paste0(outpath, "Infections_", sim_type, ".pdf"))
infections_plot(incid)
dev.off()


######################################################################
##################### Estimate Rt ####################################
######################################################################

Rt <- Rt_est(incid, vals, 'gamma', params, deterministic = T, correct_bias = T, variant = T)

Rt_nonpara <- Rt_est_nonpara(incid, samps, 'nsr', params)

# plot
png(file = paste0(outpath, "CompareRt_Variant", params[['sim_type']], ".png"))
compare_rt(Rt, params, variant = T)
dev.off()

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


