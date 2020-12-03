# import libraries
from scipy.integrate import odeint
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()

# If output directory does not already exist, create one
if not os.path.isdir("Output"):
    os.mkdir("Output")
output_dir = curr_dir + "/Output/"

data_dir = curr_dir + "/Data/"


def linear_step_func(x, x0, x1):
    R0 = np.piecewise(x, [
        x < x0,
        (x >= x0) & (x <= x1),
        x > x1],
                     [1.75,
                      1.25,
                      0.8]
                     )
    return R0

def r0_dyn(t, r0=3, mu=0.03, r_bar=1.1):
    R0 = r0 * np.exp(- mu * t) + (1 - np.exp(- mu * t)) * r_bar
    return R0


def est_sir_dyn(N, i_0, r_0, R0, gamma, t, dynamic=False):
    # initialize parameters
    N = N
    i_0 = i_0
    r_0 = r_0
    s_0 = N - i_0 - r_0

    gamma = gamma

    # initial conditions
    x_0 = s_0, i_0, r_0

    if dynamic:
        def R0(t): return r0_dyn(t)
    else:
        R0 = R0

    # define function
    def F(x, t, R0=R0):
        s, i, r = x

        # New exposure of susceptibles
        beta = R0(t) * gamma if callable(R0) else R0 * gamma
        # R0 can be constant or function of time

        # Time derivatives
        ds = -beta * s * i / N
        di = beta * s * i / N - gamma * i
        dr = gamma * i

        return ds, di, dr

    def solve_path(R0, t_vec, x_init=x_0):
        def G(x, t): return F(x, t, R0)
        s_path, i_path, r_path = odeint(G, x_init, t_vec).transpose()

        return s_path, i_path, r_path

    S, I, R = solve_path(R0, t, x_init=x_0)
    SIR = pd.DataFrame({'S': S, 'I': I, 'R': R, 'Days': np.arange(0, len(t))})
    if dynamic:
        SIR['R0'] = r0_dyn(t)
    else:
        SIR['R0'] = R0
    return SIR

t = np.linspace(0, 300, 300)
gamma = 1/16
SIR_step = est_sir_dyn(100000, 100, 0, 1.6, gamma, t, dynamic=True)
SIR_step['Delta'] = SIR_step['S'].diff()*-1

SIR_step.to_csv(data_dir + 'SimulatedSIR_dyn.csv', sep=';', index=False)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Days')
ax1.set_ylabel('R0', color=color)
ax1.plot(SIR_step['Days'], SIR_step['R0'], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Infected', color=color)  # we already handled the x-label with ax1
ax2.plot(SIR_step['Days'], SIR_step['I'], color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('γ: ' + str(round(gamma, 2)), xy=(1, 0), xycoords='axes fraction', fontsize=12,
            xytext=(-5, 5), textcoords='offset points',
            ha='right', va='bottom')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(output_dir + 'SimulatedR0_dyn.png')
plt.close()

SIR = est_sir_dyn(100000, 100, 0, 1.6, gamma, t, dynamic=False)
SIR['Delta'] = SIR['S'].diff()*-1

SIR.to_csv(data_dir + 'SimulatedSIR.csv', sep=';', index=False)


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Days')
ax1.set_ylabel('R0', color=color)
ax1.plot(SIR['Days'], SIR['R0'], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Infected', color=color)  # we already handled the x-label with ax1
ax2.plot(SIR['Days'], SIR['I'], color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('γ: ' + str(round(gamma, 2)), xy=(1, 0), xycoords='axes fraction', fontsize=12,
            xytext=(-5, 5), textcoords='offset points',
            ha='right', va='bottom')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(output_dir + 'SimulatedR0.png')
plt.close()
