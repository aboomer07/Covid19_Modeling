# old notes, delete when not helpful:

from scipy.integrate import odeint
import numpy as np

## build SIR model
# Total population, N.
N = 8000000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 143, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days). ... beta is probability of infection?
beta, gamma = 0.2, 1./14
# A grid of time points (in days)
t = np.linspace(0, 23, 23)

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

