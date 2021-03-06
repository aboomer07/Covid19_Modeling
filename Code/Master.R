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
dist <- serinfect$dist


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

SI.bias <- params_bias(params)

pdf(file = paste0(outpath, "SerialEstBias_", sim_type, "_S", simulations, ".pdf"))
par(mfrow = c(2,2))
plot(SI.bias$"Rt 0.8", type = "l", main = "Rt 0.8",
	xlab = "1/Delta", ylab = expression(hat(mu) - mu))
plot(SI.bias$"Rt 1.25", type = "l", main = "Rt 1.25",
	xlab = "1/Delta", ylab = expression(hat(mu) - mu))
plot(SI.bias$"Rt 1.5", type = "l", main = "Rt 1.5",
	xlab = "1/Delta", ylab = expression(hat(mu) - mu))
plot(SI.bias$"Rt 1.75", type = "l", main = "Rt 1.75",
	xlab = "1/Delta", ylab = expression(hat(mu) - mu))
dev.off()

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
############## Evaluate Serial Interval ##############################
######################################################################

# non parametric
nonpar_sim <- nonpara_eval(params, 'nsr')
plot_nonpara_distplot(nonpar_sim, 'ridge', params)
plot_nonpara_distplot(nonpar_sim, 'violin', params)

true_dist <- gen_distribution(params[['study_len']], params[['sim_mu']], params[['sim_var']], params[['sim_type']], 1)$omega

png(paste0(outpath, 'SerialEst_nonpara.png'))
plot_nonpara_eval(true_dist, nonpar_sim)
dev.off()

######################################################################
############## Simulate Incidence ####################################
######################################################################
# SI simulations
# first: gamma true process, constant Rt
# second: gamma true process, varying Rt
# third: weibull true process, constant Rt
# fourth: wibull true process, varying Rt

# simulate outbreaks
tot_len <- params$n_days * params$delta

params[["R_val"]] <- 1.4
params[['R_val_variant']] <- 1.9
incid_si_gamma_const <- si_sim(params)
incid_sii_gamma_const <- sii_sim(params)

Rt_si_determ_const <- Rt_est(incid_si_gamma_const, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_si_sto_const <- Rt_est(incid_si_gamma_const, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt_sii_determ_const <- Rt_est(incid_sii_gamma_const, vals, params, deterministic = T, correct_bias = T, variant = T)
Rt_sii_sto_const <- Rt_est(incid_sii_gamma_const, vals, params, deterministic = F, correct_bias = T, variant = T)

pdf(file = paste0(outpath, "Outbreak_SI_const_gamma.pdf"), width=4, height=7)
si_plot_detail(incid_si_gamma_const)
dev.off()

pdf(file = paste0(outpath, "CompareRt_SI_determ_const.pdf"), height = 4, width=8)
compare_rt(Rt_si_determ_const, params)
dev.off()
pdf(file = paste0(outpath, "CompareRt_SI_sto_const.pdf"), height = 4, width=8)
compare_rt(Rt_si_sto_const, params)
dev.off()

pdf(file = paste0(outpath, "Outbreak_SII_const_gamma_determ.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_const, Rt_sii_determ_const)
dev.off()
pdf(file = paste0(outpath, "Outbreak_SII_const_gamma_sto.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_const, Rt_sii_sto_const)
dev.off()

params[['R_val']] <- seq(1, 2.2, by=1.2/tot_len)
params[['R_val']] <- params[['R_val']][2:length(params[['R_val']])]
params[['R_val_variant']] <- 1.9
incid_si_gamma_inc <- si_sim(params)
incid_sii_gamma_inc <- sii_sim(params)

Rt_si_determ_inc <- Rt_est(incid_si_gamma_inc, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_si_sto_inc <- Rt_est(incid_si_gamma_inc, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt_sii_determ_inc <- Rt_est(incid_sii_gamma_inc, vals, params, deterministic = T, correct_bias = T, variant = T)
Rt_sii_sto_inc <- Rt_est(incid_sii_gamma_inc, vals, params, deterministic = F, correct_bias = T, variant = T)

pdf(file = paste0(outpath, "Outbreak_SI_inc_gamma.pdf"), width=4, height=7)
si_plot_detail(incid_si_gamma_inc)
dev.off()

pdf(file = paste0(outpath, "CompareRt_SI_determ_inc.pdf"), height = 4, width=8)
compare_rt(Rt_si_determ_inc, params)
dev.off()
pdf(file = paste0(outpath, "CompareRt_SI_sto_inc.pdf"), height = 4, width=8)
compare_rt(Rt_si_sto_inc, params)
dev.off()

pdf(file = paste0(outpath, "Outbreak_SII_inc_gamma_determ.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_inc, Rt_sii_determ_inc)
dev.off()
pdf(file = paste0(outpath, "Outbreak_SII_inc_gamma_sto.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_inc, Rt_sii_sto_inc)
dev.off()

params[['R_val']] <- rev(seq(1, 2.2, by=1.2/tot_len))
params[['R_val']] <- params[['R_val']][2:length(params[['R_val']])]
params$R_val_variant <- 1.9
incid_si_gamma_dec <- si_sim(params)
incid_sii_gamma_dec <- sii_sim(params)

Rt_si_determ_dec <- Rt_est(incid_si_gamma_dec, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_si_sto_dec <- Rt_est(incid_si_gamma_dec, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt_sii_determ_dec <- Rt_est(incid_sii_gamma_dec, vals, params, deterministic = T, correct_bias = T, variant = T)
Rt_sii_sto_dec <- Rt_est(incid_sii_gamma_dec, vals, params, deterministic = F, correct_bias = T, variant = T)

pdf(file = paste0(outpath, "Outbreak_SI_dec_gamma.pdf"), width=4, height=7)
si_plot_detail(incid_si_gamma_dec)
dev.off()

pdf(file = paste0(outpath, "CompareRt_SI_determ_dec.pdf"), height = 4, width=8)
compare_rt(Rt_si_determ_dec, params)
dev.off()
pdf(file = paste0(outpath, "CompareRt_SI_sto_dec.pdf"), height = 4, width=8)
compare_rt(Rt_si_sto_dec, params)
dev.off()

pdf(file = paste0(outpath, "Outbreak_SII_dec_gamma_determ.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_dec, Rt_sii_determ_dec)
dev.off()
pdf(file = paste0(outpath, "Outbreak_SII_dec_gamma_sto.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_dec, Rt_sii_sto_dec)
dev.off()

params[['R_val']] <- (0.5 * sin((1:tot_len)/1000)) + 1.2
# params[['R_val']] <- params[['R_val']][2:length(params[['R_val']])]
params$R_val_variant <- 1.9
incid_si_gamma_sin <- si_sim(params)
incid_sii_gamma_sin <- sii_sim(params)

Rt_si_determ_sin <- Rt_est(incid_si_gamma_sin, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_si_sto_sin <- Rt_est(incid_si_gamma_sin, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt_sii_determ_sin <- Rt_est(incid_sii_gamma_sin, vals, params, deterministic = T, correct_bias = T, variant = T)
Rt_sii_sto_sin <- Rt_est(incid_sii_gamma_sin, vals, params, deterministic = F, correct_bias = T, variant = T)

pdf(file = paste0(outpath, "Outbreak_SI_sin_gamma.pdf"), width=4, height=7)
si_plot_detail(incid_si_gamma_sin)
dev.off()

pdf(file = paste0(outpath, "CompareRt_SI_determ_sin.pdf"), height = 4, width=8)
compare_rt(Rt_si_determ_sin, params)
dev.off()
pdf(file = paste0(outpath, "CompareRt_SI_sto_sin.pdf"), height = 4, width=8)
compare_rt(Rt_si_sto_sin, params)
dev.off()

pdf(file = paste0(outpath, "Outbreak_SII_sin_gamma_determ.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_sin, Rt_sii_determ_sin)
dev.off()
pdf(file = paste0(outpath, "Outbreak_SII_sin_gamma_sto.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_sin, Rt_sii_sto_sin)
dev.off()

params[['R_val']] <- -1 * (2/tot_len * (seq(1:tot_len) - (tot_len/2)))**2 + 2
params$R_val_variant <- 1.8
incid_si_gamma_x2 <- si_sim(params)
incid_sii_gamma_x2 <- sii_sim(params)

Rt_si_determ_x2 <- Rt_est(incid_si_gamma_x2, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_si_sto_x2 <- Rt_est(incid_si_gamma_x2, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt_sii_determ_x2 <- Rt_est(incid_sii_gamma_x2, vals, params, deterministic = T, correct_bias = T, variant = T)
Rt_sii_sto_x2 <- Rt_est(incid_sii_gamma_x2, vals, params, deterministic = F, correct_bias = T, variant = T)

pdf(file = paste0(outpath, "Outbreak_SI_x2_gamma.pdf"), width=4, height=7)
si_plot_detail(incid_si_gamma_x2)
dev.off()

pdf(file = paste0(outpath, "CompareRt_SI_determ_x2.pdf"), height = 4, width=8)
compare_rt(Rt_si_determ_x2, params)
dev.off()
pdf(file = paste0(outpath, "CompareRt_SI_sto_x2.pdf"), height = 4, width=8)
compare_rt(Rt_si_sto_x2, params)
dev.off()

pdf(file = paste0(outpath, "Outbreak_SII_x2_gamma_determ.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_x2, Rt_sii_determ_x2)
dev.off()
pdf(file = paste0(outpath, "Outbreak_SII_x2_gamma_sto.pdf"), width=4, height=7)
sii_plot(incid_sii_gamma_x2, Rt_sii_sto_x2)
dev.off()

# params[["R_val"]] <- c(1.7, 0.9, 2.5)
# incid_si_gamma_var <- si_sim(params)

# params[["R_val"]] <- 1.7
# params[["sim_type"]] <- "weibull"
# incid_si_weibull_const <- si_sim(params)

# params[["R_val"]] <- c(1.7, 0.9, 1.3)
# params[["sim_type"]] <- "weibull"
# incid_si_weibull_var <- si_sim(params)

# plot outbreaks si
# pdf(file = paste0(outpath, "Outbreak_SI_const_gamma.pdf"), width=4, height=7)
# si_plot_detail(incid_si_gamma_const)
# dev.off()

# pdf(file = paste0(outpath, "Outbreak_SI_var_gamma.pdf"), width=4, height=7)
# si_plot_detail(incid_si_gamma_var)
# dev.off()

# pdf(file = paste0(outpath, "Outbreak_SI_const_weibull.pdf"), width=4, height=7)
# si_plot_detail(incid_si_weibull_const)
# dev.off()

# pdf(file = paste0(outpath, "Outbreak_SI_var_weibull.pdf"), width=7, height=7)
# si_plot_detail(incid_si_weibull_var)
# dev.off()



######################################################################
##################### Estimate Rt ####################################
######################################################################

# estimate si for deterministic & stochastic gamma and nonparametrically for
# the four different models
# First model
Rt1_si_determ <- Rt_est(incid_si_gamma_const, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt1_si_sto <- Rt_est(incid_si_gamma_const, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt1_nonpara_si <- Rt_est_nonpara(incid_si_gamma_const, samps, 'nsr', params, correct_bias = T)

# Second model
Rt2_si_determ <- Rt_est(incid_si_gamma_var, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt2_si_sto <- Rt_est(incid_si_gamma_var, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt2_nonpara_si <- Rt_est_nonpara(incid_si_gamma_var, samps, 'nsr', params, correct_bias = T)

# # Third model
# # need to rerun the serial interval estimation with a weibull
vals <- serial_ests(samp_pois(params)$daily)
# Rt3_si_sto <- Rt_est(incid_si_weibull_const, vals, params, deterministic = F, correct_bias = T, variant = F)
# Rt3_nonpara_si <- Rt_est_nonpara(incid_si_weibull_const, samps, 'nsr', params, correct_bias = T)

# Fourth model
Rt4_si_sto <- Rt_est(incid_si_weibull_var, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt4_nonpara_si <- Rt_est_nonpara(incid_si_weibull_var, samps, 'nsr', params, correct_bias = T)



# plot different sims and estimations
# Model 1
pdf(file = paste0(outpath, "CompareRt1_SI_determ.pdf"), height = 4, width=8)
compare_rt(Rt1_si_determ, params)
dev.off()

pdf(file = paste0(outpath, "CompareRt1_SI_sto.pdf"), height = 4, width=8)
compare_rt(Rt1_si_sto, params)
dev.off()

pdf(file = paste0(outpath, "CompareRt1_SI_nonpara.pdf"), height = 4, width=8)
compare_rt(Rt1_nonpara_si, params)
dev.off()

# Model 2
pdf(file = paste0(outpath, "CompareRt2_SI_determ.pdf"), height = 4, width=8)
compare_rt(Rt2_si_determ, params)
dev.off()

pdf(file = paste0(outpath, "CompareRt2_SI_sto.pdf"), height = 4, width=8)
compare_rt(Rt2_si_sto, params)
dev.off()

pdf(file = paste0(outpath, "CompareRt2_SI_nonpara.pdf"), height = 4, width=8)
compare_rt(Rt2_nonpara_si, params)
dev.off()

# Model 3
# png(file = paste0(outpath, "CompareRt3_SI_sto.png"))
# compare_rt(Rt3_si_sto, params)
# dev.off()
#
# png(file = paste0(outpath, "CompareRt3_SI_nonpara.png"))
# compare_rt(Rt3_nonpara_si, params)
# dev.off()

# Model 4
pdf(file = paste0(outpath, "CompareRt4_SI_sto.pdf"), height = 4, width=8)
compare_rt(Rt4_si_sto, params)
dev.off()

pdf(file = paste0(outpath, "CompareRt4_SI_nonpara.pdf"), height = 4, width=8)
compare_rt(Rt4_nonpara_si, params)
dev.off()


######################################################################
##################### SI Framework ###################################
######################################################################

# params$R_val <- c(1.3, 1.05)
params$R_val_variant <- 1.8

params[["R_val"]] <- 1.4
incid_sii_gamma_const <- sii_sim(params)

params[['R_val']] <- seq(0.5, 2, by=1.5/tot_len)
incid_sii_gamma_inc <- sii_sim(params)

params[['R_val']] <- rev(seq(0.5, 2, by=1.5/tot_len))
incid_sii_gamma_dec <- sii_sim(params)

params[['R_val']] <- (0.5 * sin((1:tot_len)/1000)) + 1.2
incid_sii_gamma_sin <- sii_sim(params)

params[['R_val']] <- -1 * (2/tot_len * (seq(1:tot_len) - (tot_len/2)))**2 + 2
incid_sii_gamma_x2 <- sii_sim(params)

# Rt_si <- Rt_est(si_model, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_sii <- Rt_est(sii_model, vals, params, deterministic = T, correct_bias = T, variant = T, sep_Rt = T)
Rt_nsr_sii <- Rt_est_nonpara(incid_sii, samps, 'nsr', params, correct_bias = T, variant = T, sep_Rt = T)
# Rt_ssii <- Rt_est(ssii_model, vals, params, deterministic = T, correct_bias = T, variant = T, sep_Rt = T, sep_S = T)

# si_plot(si_model, Rt_si)
png(paste0(outpath, "/SII_Plot_Rt_Constant.png"))
sii_plot(sii_model, Rt_sii)
dev.off()




