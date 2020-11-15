# execute basic SIR model

from Code.SIRModel import *
from Code.SIRParams import *
from Code.SIRPlot import *

# get parameters
aut = get_params('Austria', '2020-06-10', df)
ita = get_params('Italy', '2020-06-10', df)

SIR_aut = est_sir(**aut)
plot_sir(SIR_aut, 'Austria')

