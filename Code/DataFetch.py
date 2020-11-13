# import libraries
import requests
import pandas as pd
import os

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()

# If an output directory does not already exist, create one
if not os.path.isdir("Data"):
    os.mkdir("Data")
data_dir = curr_dir + "/Data"

# 1. Option: full historic data from github
# A. All countries and US total
url = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'
df = pd.read_csv(url, sep=',')
df = df.drop(['Lat', 'Long'], axis=1)

# wide to long
df = df.melt(id_vars=['Province/State', 'Country/Region'], var_name='date', value_name='confirmed')
df['date'] = pd.to_datetime(df['date'], format="%m/%d/%y")

# restrict to countries of interest (here US is grouped?)
country_list = ['Austria', 'Germany', 'France', 'Italy', 'US']
df = df[df['Country/Region'].isin(country_list)]

# save csv
df.to_csv(data_dir + 'EuropeCovidData.csv', sep=';', index=False)

# A. US separate to get county level
url = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
df_us = pd.read_csv(url, sep=',')
df_us = df_us.drop(['Lat', 'Long_', 'UID', 'iso2', 'iso3', 'code3', 'FIPS', 'Country_Region', 'Combined_Key'], axis=1).reset_index(drop=True)
# check if UUID or FIPS are useful

# wide to long
df_us = df_us.melt(id_vars=['Admin2', 'Province_State'], var_name='date', value_name='confirmed')
df_us['date'] = pd.to_datetime(df_us['date'], format="%m/%d/%y")

# restrict to counties of interest, US specific
county_list = ['Oregon']
df_us = df_us[df_us['Province_State'].isin(county_list)]

# save data
df_us.to_csv(data_dir + 'USCovidData.csv', sep=';', index=False)

# 2. Option: specific historic data from hopkins

api_response = requests.get('https://covid19api.herokuapp.com/confirmed')

# get latest data for location ID 16: Austria
total_cases_aut = api_response.json()['locations'][16]['history']
