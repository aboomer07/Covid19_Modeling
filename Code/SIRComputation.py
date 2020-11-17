# execute basic SIR model

from Code.SIRModel import *
from Code.SIRParams import *
from Code.SIRPlot import *

# get parameters
aut = get_params('Austria', '2020-05-10', df, 160)
ita = get_params('Italy', '2020-03-10', df, 160)

# static versions
SIR_aut = est_sir(**aut)
plot_sir(SIR_aut, 'Austria')

# static version multiple R0
R0_vec = [1.6, 3, 6]

SIRs_aut = est_sir_multi_r0(aut, R0_vec)
plot_multiple_sir(SIRs_aut, 'Austria')

# dynamic versions (changing R0 over time)
# dynamic = False gives equivalent result to est_sir :)
SIR_aut_dyn = est_sir_dyn(**aut, dynamic=True)
plot_sir(SIR_aut_dyn, 'Austria_Dyn')


################################################################################
# Code for running US SIR model, needed functions are on SIRParams file
################################################################################
sir_data = get_sir_data(df_us)
states = ['NY', 'NJ', 'CT']
counties = list(
    sir_data[sir_data['StateName'].isin(states)]['countyFIPS'].unique())
us = get_params_us(counties, '04-20-2020', sir_data, 160)
SIR_us = est_sir(**us)
plot_sir(SIR_us, " + ".join(states))

SIRs_us = est_sir_multi_r0(us, R0_vec)
plot_multiple_sir(SIRs_us, " + ".join(states))

SIR_us_dyn = est_sir_dyn(**us, dynamic=True)
plot_sir(SIR_us_dyn, "Dynamic " + " + ".join(states))
