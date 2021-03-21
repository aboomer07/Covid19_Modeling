# import libraries
from scipy.integrate import odeint
import pandas as pd
import numpy as np

# Function arguments:
# N = total population
# i_0 = initial number of infected
# r_0 = initial number of recovered
# beta = contact rate (beta = R0*gamma)
# gamma = mean recovery rate
# t = grid of time points in days


def est_sir(N, i_0, r_0, R0, gamma, t):
    N = N
    i_0 = i_0
    r_0 = r_0

    gamma = gamma
    beta = R0*gamma

    # Everyone else, S0, is susceptible to infection initially.
    s_0 = N - i_0 - r_0

    # The SIR model differential equations.
    def F(y, t, N, beta, gamma):
        S, I, R = y
        dSdt = -beta * S * I / N
        dIdt = beta * S * I / N - gamma * I
        dRdt = gamma * I
        return dSdt, dIdt, dRdt

    # Initial conditions vector
    y0 = s_0, i_0, r_0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(F, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T / N
    SIR = pd.DataFrame({'S': S, 'I': I, 'R': R, 'Days': np.arange(0, len(t))})
    return SIR


def est_sir_multi_r0(params, r0_vec):
    SIRs = []
    new_params = params
    for i in r0_vec:
        new_params['R0'] = i
        SIR = est_sir(**new_params)
        SIR['R0'] = i
        SIRs.append(SIR)
    SIRs = pd.concat(SIRs)
    return SIRs


# Dynamic function
# Function arguments:
# x = state vector (array like)
# t is time (scalar)
# !! R0 is effective transmission rate, defaulting to constant

# R0 changing over time (starting at 1.5 and going down to 1.1 at a rate of 0.06)
def r0_dyn(t, r0=3, mu=0.06, r_bar=1.1):
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
    return SIR
