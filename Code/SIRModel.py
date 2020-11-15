# import libraries
from scipy.integrate import odeint
import pandas as pd
import numpy as np

# Function arguments:
# N = total population
# I0 = initial number of infected
# R0 = initial number of recovered
# beta = contact rate
# gamma = mean recovery rate
# t = grid of time points in days

def est_sir(N, I0, R0, beta, gamma, t):
    N = N
    I0 = I0
    R0 = R0
    beta = beta
    gamma = gamma
    t = t

    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0

    # The SIR model differential equations.
    def deriv(y, t, N, beta, gamma):
        S, I, R = y
        dSdt = -beta * S * I / N
        dIdt = beta * S * I / N - gamma * I
        dRdt = gamma * I
        return dSdt, dIdt, dRdt

    # Initial conditions vector
    y0 = S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T
    SIR = pd.DataFrame({'S': S, 'I': I, 'R': R, 'Days': np.arange(0, len(t))})
    return SIR








