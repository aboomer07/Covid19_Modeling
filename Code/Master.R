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
# first with constant Rt
# simulate outbreak
incid_si <- si_sim(params)
incid_sii <- sii_sim(params)

# plot outbreak si
png(file = paste0(outpath, "Outbreak_SI_", params[['sim_type']], ".png"))
si_plot_detail(incid_si)
dev.off()

# plot outbreak sii
png(file = paste0(outpath, "Outbreak_SII_", params[['sim_type']], ".png"))
sii_plot(incid_sii)
dev.off()


######################################################################
##################### Estimate Rt ####################################
######################################################################

Rt_si <- Rt_est(incid_si, vals, params, deterministic = T, correct_bias = T, variant = F)
Rt_sii <- Rt_est(incid_sii, vals, params, deterministic = T, correct_bias = T, variant = T)

Rt_nonpara_si <- Rt_est_nonpara(incid_si, samps, 'nsr', params, correct_bias = T)
Rt_nonpara_sii <- Rt_est_nonpara(incid_sii, samps, 'nsr', params, correct_bias = T, variant = T)


# plot different sims and estimations
png(file = paste0(outpath, "CompareRt_SI_", params[['sim_type']], ".png"))
compare_rt(Rt_si, params)
dev.off()

png(file = paste0(outpath, "CompareRt_SI_", params[['sim_type']], ".png"))
compare_rt(Rt_sii, params, variant = T)
dev.off()

png(file = paste0(outpath, "CompareRt_SI_nonpara", params[['sim_type']], ".png"))
compare_rt(Rt_nonpara_si, params)
dev.off()

png(file = paste0(outpath, "CompareRt_SII_nonpara", params[['sim_type']], ".png"))
compare_rt(Rt_nonpara_sii, params, variant = T)
dev.off()


