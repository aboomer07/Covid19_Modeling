# standalone script for some descriptive statistics and plots


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from Code.DataFunc import *
from Code.SIRModel import r0_dyn

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()
data_dir = curr_dir + "/Data/"

# If output directory does not already exist, create one
if not os.path.isdir("Output"):
    os.mkdir("Output")
output_dir = curr_dir + "/Output/"

# import data
df = pd.read_csv(data_dir + 'EuropeCovidData.csv', sep=';')
country_list = df['Country/Region'].unique()

df['date'] = pd.to_datetime(df['date'])

# 1. Plot: confirmed cases by country

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=df, x='date', y='confirmed', hue='Country/Region', ax=ax)
ax.set_title('Confirmed Covid-19 Cases by Country')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'CasesByCountry.png')
plt.close()

# 2. Plot: empirical R0, RKI method
R0s = []
for country in country_list:
    R0 = get_R0(country, df)
    # R0 = R0.rename(country)
    R0s.append(R0)

R0s = pd.concat(R0s)
R0s = R0s.replace([np.inf, -np.inf], np.nan)
R0s = R0s.dropna()


fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=R0s, x='date', y='R0', hue='Country/Region', ax=ax)
ax.set_title('R0 by RKI method (rolling window)')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'R0ByCountry.png')
plt.close()

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=R0s[R0s['date'] > '2020-03-10'], x='date', y='R0', hue='Country/Region', ax=ax)
ax.set_title('R0 by RKI method (rolling window)')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'R0ByCountry_Subset_03_10.png')
plt.close()

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=R0s[R0s['date'] > '2020-03-30'], x='date', y='R0', hue='Country/Region', ax=ax)
ax.set_title('R0 by RKI method (rolling window)')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'R0ByCountry_Subset_03_30.png')
plt.close()

# 3. Plot: Compare to degrading R0 for Austria
R0_aut = R0s[(R0s['date'] > '2020-03-30') & (R0s['Country/Region'] == 'Austria')]
t = len(R0s[(R0s['date'] > '2020-03-30') & (R0s['Country/Region'] == 'Austria')].confirmed)
timespan = np.linspace(0, t, t)

R0_aut['R0_sim'] = r0_dyn(timespan, r0=1.5, r_bar=1.1, mu=0.06)
R0_aut = R0_aut.melt(id_vars=['date'], value_vars=['R0', 'R0_sim'], var_name='Type')

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=R0_aut, x='date', y='value', style='Type', ax=ax)
ax.set_title('R0 by RKI method vs time decaying R0')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'R0Sim_Aut.png')
plt.close()


# 3. Plot: Cases in Austria with events

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=df[df['Country/Region'] == 'Austria'], x='date', y='confirmed', ax=ax)
ax.set_title('Confirmed Cases Austria')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'Cases_Aut.png')
plt.close()
