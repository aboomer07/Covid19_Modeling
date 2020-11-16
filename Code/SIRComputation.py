# execute basic SIR model

from Code.SIRModel import *
from Code.SIRParams import *
from Code.SIRPlot import *

# get parameters
aut = get_params('Austria', '2020-05-10', df)
ita = get_params('Italy', '2020-03-10', df)

# static versions
SIR_aut = est_sir(**aut)
plot_sir(SIR_aut, 'Austria')


# dynamic versions (changing R0 over time)
SIR_aut_dyn = est_sir_dyn(**aut, dynamic=True)  # dynamic = False gives equivalent result to est_sir :)
plot_sir(SIR_aut_dyn, 'Austria_Dyn')

