# execute basic SIR model

from Code.SIRModel import *
from Code.SIRParams import *
from Code.SIRPlot import *

country_list = df['Country/Region'].unique()

# get parameters
params = []
for country in country_list:
    param = get_params(country, '2020-03-20', df, 350)
    params.append(param)

# static versions
for country_params, country in zip(params, country_list):
    SIR = est_sir(**country_params)
    plot_sir(SIR, country_params, country)

# static version multiple R0
R0_vec = [1.6, 3, 6]

for country_params, country in zip(params, country_list):
    SIR = est_sir_multi_r0(country_params, R0_vec)
    plot_multiple_sir(SIR, country_params, country)

# dynamic versions (changing R0 over time)
# dynamic = False gives equivalent result to est_sir :)
# SIR_aut_dyn = est_sir_dyn(**aut, dynamic=True)
# plot_sir(SIR_aut_dyn, aut, 'Austria_Dyn')
# plot_i(SIR_aut_dyn, aut, 'Austria_Dny_Infected')  # plot only infected


################################################################################
# Code for running US SIR model, needed functions are on SIRParams file
################################################################################
# sir_data = get_sir_data(df_us)
# states = ['NY', 'NJ', 'CT']
# counties = list(
#     sir_data[sir_data['StateName'].isin(states)]['countyFIPS'].unique())
# us = get_params_us(counties, '04-20-2020', sir_data, 160)
# SIR_us = est_sir(**us)
# plot_sir(SIR_us, " + ".join(states))
#
# SIRs_us = est_sir_multi_r0(us, R0_vec)
# plot_multiple_sir(SIRs_us, " + ".join(states))
#
# SIR_us_dyn = est_sir_dyn(**us, dynamic=True)
# plot_sir(SIR_us_dyn, "Dynamic " + " + ".join(states))
