# standalone script for some descriptive statistics and plots


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px
import plotly.graph_objs as go


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
events = pd.read_csv(data_dir + 'CovidEventsAUT.csv', sep=';')


# Data preparation for plots
df_event_aut = df[df['Country/Region'] == 'Austria'].merge(events, on='date', how='left')

country_list = df['Country/Region'].unique()

R0s = []
for country in country_list:
    R0 = get_R0(country, df)
    R0s.append(R0)

R0s = pd.concat(R0s)
R0s = R0s.replace([np.inf, -np.inf], np.nan)
R0s = R0s[R0s['R0'].notna()]

R0s_event_aut = R0s[R0s['Country/Region'] == 'Austria'].merge(events, on='date', how='left')

df['date'] = pd.to_datetime(df['date'])
R0s['date'] = pd.to_datetime(R0s['date'])
df_event_aut['date'] = pd.to_datetime(df_event_aut['date'])
R0s_event_aut['date'] = pd.to_datetime(R0s_event_aut['date'])

# new infections daily
df_event_aut['delta'] = df_event_aut['confirmed'].diff()


# 1. Plot: confirmed cases by country

fig, ax = plt.subplots(ncols=1, nrows=1)
sns.lineplot(data=df, x='date', y='confirmed', hue='Country/Region', ax=ax)
ax.set_title('Confirmed Covid-19 Cases by Country')
fig.set_size_inches(18.5, 10.5)
plt.savefig(output_dir + 'CasesByCountry.png')
plt.close()

# 2. Plot: empirical R0, RKI method

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

# confirmed cases interactive with events
df_event_aut_sub = df_event_aut[df_event_aut['event'].notnull()]
fig = px.line(df_event_aut, x='date', y="confirmed")
fig.add_trace(go.Scatter(mode="markers", x=df_event_aut_sub['date'], y=df_event_aut_sub['confirmed'],
                         hovertext=df_event_aut_sub['event'], name='Event'))
fig.write_html(output_dir + "ConfirmedCasesEvents.html")

# delta cases interactive with events
fig = px.line(df_event_aut, x='date', y="delta")
fig.add_trace(go.Scatter(mode="markers", x=df_event_aut_sub['date'], y=df_event_aut_sub['delta'],
                         hovertext=df_event_aut_sub['event'], name='Event'))
fig.write_html(output_dir + "DeltaCasesEvents.html")

# R0 interactive with events
# restrict to smaller time period to account for unrealistic R0 in first few days
R0s_event_aut = R0s_event_aut[R0s_event_aut['date'] > '2020-03-30']

R0s_event_aut_sub = R0s_event_aut[R0s_event_aut['event'].notnull()]
fig = px.line(R0s_event_aut, x='date', y="R0")
fig.add_trace(go.Scatter(mode="markers", x=R0s_event_aut_sub['date'], y=R0s_event_aut_sub['R0'],
                         hovertext=R0s_event_aut_sub['event'], name='Event'))
fig.write_html(output_dir + "R0Events.html")


