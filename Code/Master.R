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

plot_nonpara_eval(true_dist, nonpar_sim)


######################################################################
############## Simulate Incidence ####################################
######################################################################
# SI simulations
# first: gamma true process, constant Rt
# second: gamma true process, varying Rt
# third: weibull true process, constant Rt
# fourth: wibull true process, varying Rt

# simulate outbreaks
params[["R_val"]] <- 1.7
incid_si_gamma_const <- si_sim(params)

params[["R_val"]] <- c(1.7, 0.9, 1.3)
incid_si_gamma_var <- si_sim(params)

params[["R_val"]] <- 1.7
params[["sim_type"]] <- "weibull"
incid_si_weibull_const <- si_sim(params)

params[["R_val"]] <- c(1.7, 0.9, 1.3)
params[["sim_type"]] <- "weibull"
incid_si_weibull_var <- si_sim(params)

# plot outbreaks si
png(file = paste0(outpath, "Outbreak_SI_const_gamma.png"))
si_plot_detail(incid_si_gamma_const)
dev.off()

png(file = paste0(outpath, "Outbreak_SI_var_gamma.png"))
si_plot_detail(incid_si_gamma_var)
dev.off()

png(file = paste0(outpath, "Outbreak_SI_const_weibull.png"))
si_plot_detail(incid_si_weibull_const)
dev.off()

png(file = paste0(outpath, "Outbreak_SI_var_weibull.png"))
si_plot_detail(incid_si_weibull_var)
dev.off()

incid_sii <- sii_sim(params)

# plot outbreak sii
png(file = paste0(outpath, "Outbreak_SII_", params[['sim_type']], ".png"))
sii_plot(incid_sii)
dev.off()

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

# Third model
# need to rerun the serial interval estimation with a weibull
vals <- serial_ests(samp_pois(params)$daily)
Rt3_si_sto <- Rt_est(incid_si_weibull_const, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt3_nonpara_si <- Rt_est_nonpara(incid_si_weibull_const, samps, 'nsr', params, correct_bias = T)

# Fourth model
Rt4_si_sto <- Rt_est(incid_si_weibull_var, vals, params, deterministic = F, correct_bias = T, variant = F)
Rt4_nonpara_si <- Rt_est_nonpara(incid_si_weibull_var, samps, 'nsr', params, correct_bias = T)

Rt_sii <- Rt_est(incid_sii, vals, params, deterministic = T, correct_bias = T, variant = T)

Rt_nonpara_sii <- Rt_est_nonpara(incid_sii, samps, 'nsr', params, correct_bias = T, variant = T)

# plot different sims and estimations
# Model 1
png(file = paste0(outpath, "CompareRt1_SI_determ.png"))
compare_rt(Rt1_si_determ, params)
dev.off()

png(file = paste0(outpath, "CompareRt1_SI_sto.png"))
compare_rt(Rt1_si_sto, params)
dev.off()

png(file = paste0(outpath, "CompareRt1_SI_nonpara.png"))
compare_rt(Rt1_nonpara_si, params)
dev.off()

# Model 2
png(file = paste0(outpath, "CompareRt2_SI_determ.png"))
compare_rt(Rt2_si_determ, params)
dev.off()

png(file = paste0(outpath, "CompareRt2_SI_sto.png"))
compare_rt(Rt2_si_sto, params)
dev.off()

png(file = paste0(outpath, "CompareRt2_SI_nonpara.png"))
compare_rt(Rt2_nonpara_si, params)
dev.off()

# Model 3
png(file = paste0(outpath, "CompareRt3_SI_sto.png"))
compare_rt(Rt3_si_sto, params)
dev.off()

png(file = paste0(outpath, "CompareRt3_SI_nonpara.png"))
compare_rt(Rt3_nonpara_si, params)
dev.off()

# Model 4
png(file = paste0(outpath, "CompareRt4_SI_sto.png"))
compare_rt(Rt4_si_sto, params)
dev.off()

png(file = paste0(outpath, "CompareRt4_SI_nonpara.png"))
compare_rt(Rt4_nonpara_si, params)
dev.off()

png(file = paste0(outpath, "CompareRt_SI_", params[['sim_type']], ".png"))
compare_rt(Rt_sii, params, variant = T)
dev.off()

png(file = paste0(outpath, "CompareRt_SII_nonpara", params[['sim_type']], ".png"))
compare_rt(Rt_nonpara_sii, params, variant = T)
dev.off()

######################################################################
##################### SI Framework ###################################
######################################################################

si_model <- si_sim(params)
sii_model <- sii_sim(params)
ssii_model <- ssii_sim(params)

Rt_si <- Rt_est(si_model, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_sii <- Rt_est(sii_model, vals, params, deterministic = T, correct_bias = T, variant = T, sep_Rt = T)
Rt_ssii <- Rt_est(ssii_model, vals, params, deterministic = T, correct_bias = T, variant = T, sep_Rt = T, sep_S = T)

si_plot(si_model, Rt_si)
png(paste0(outpath, "/SII_Plot.png"))
sii_plot(sii_model, Rt_sii)
dev.off()


