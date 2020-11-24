# import libraries
# import requests
import pandas as pd
import os

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()

# If data directory does not already exist, create one
if not os.path.isdir("Data"):
    os.mkdir("Data")
data_dir = curr_dir + "/Data/"

# 1. Option: full historic data from github
# A. All countries and US total
url = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'
df = pd.read_csv(url, sep=',')
df = df.drop(['Lat', 'Long'], axis=1)

# wide to long
df = df.melt(id_vars=['Province/State', 'Country/Region'],
             var_name='date', value_name='confirmed')
df['date'] = pd.to_datetime(df['date'], format="%m/%d/%y")

# restrict to countries of interest (here US is grouped?)
country_list = ['Austria', 'Germany', 'France', 'Italy', 'US']
df = df[df['Country/Region'].isin(country_list)]

# french regions are not of interest
df = df.groupby(['Country/Region', 'date'])['confirmed'].sum().reset_index()


# save csv
df.to_csv(data_dir + 'EuropeCovidData.csv', sep=';', index=False)

# B. US separate to get county level
url = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
df_us = pd.read_csv(url, sep=',')
df_us = df_us.drop(['Lat', 'Long_', 'UID', 'iso2', 'iso3', 'code3',
                    'Country_Region', 'Combined_Key', 'Admin2'], axis=1).reset_index(drop=True)
df_us = df_us.rename({'Province_State': "State"}, axis=1)
df_us['FIPS'] = df_us['FIPS'].apply(
    lambda x: str(int(x)) if pd.notnull(x) else "NaN")
df_us['FIPS'] = df_us['FIPS'].apply(
    lambda x: "0" + x if len(x) == 4 else x)
# FIPS is the county code, keeping instead of Admin2 which is county name

# wide to long
df_us = df_us.melt(id_vars=['FIPS', 'State'],
                   var_name='date', value_name='confirmed')
df_us['date'] = pd.to_datetime(df_us['date'], format="%m/%d/%y")

# restrict to counties of interest, US specific
# county_list = ['Oregon']
# df_us = df_us[df_us['State'].isin(county_list)]

# save data
df_us.to_csv(data_dir + 'USCovidData.csv', sep=';', index=False)
