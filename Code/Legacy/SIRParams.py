# specify parameters to be used in model runs

# import libraries
import pandas as pd
import numpy as np
import os
import world_bank_data as wb
import country_converter as coco

# import data
# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("../..")
curr_dir = os.getcwd()
data_dir = curr_dir + "/Data/"

europe = input("Do you want to input European data? True or False: ")
us = input("Do you want to input US County Level Data? True or False: ")

if europe:
    df = pd.read_csv(data_dir + 'EuropeCovidData.csv', sep=';')

# if us:
#     df_us = pd.read_csv(data_dir + "county_data.csv", sep=",")
#     df_us['countyFIPS'] = df_us['countyFIPS'].astype(str).str.zfill(5)


def get_params(country, startdate, data, timespan):
    params = {
        'N': int(wb.get_series('SP.POP.TOTL', mrv=1, country=coco.convert(country, to='ISO3'))),
        'i_0': int(data[(data['Country/Region'] == country) &
                        (data['date'] == startdate)]['confirmed']),
        'r_0': 0,  # need more data to specify
        'R0': 1.6,  # derive this from data
        'gamma': 1/18,  # based on Atkenson's note for now
        't': np.linspace(0, timespan, timespan),
    }
    return params


def get_sir_data(data):
    data = data.copy(deep=True)
    death_cols = [col for col in data.columns if "#Deaths_" in col]
    case_cols = [col for col in data.columns if "#Cases_" in col]
    data = data[['countyFIPS', 'StateName', 'PopulationEstimate2018'] +
                death_cols + case_cols]
    data = data.set_index(['countyFIPS', 'StateName', 'PopulationEstimate2018'])
    data.columns = pd.MultiIndex.from_tuples(
        [(i.split("_")[0], i.split("_")[1]) for i in data.columns])
    data = data.stack(level=1).reset_index().\
        rename({'level_3': "Date", '#Cases': 'confirmed',
                '#Deaths': 'deaths', 'PopulationEstimate2018': 'Pop'}, axis=1)

    return(data)


def get_params_us(counties, startdate, data, timespan):
    data = data.copy(deep=True)
    params = {
        'N': data[(data['countyFIPS'].isin(counties)) & (data['Date'] == startdate)]['Pop'].sum(),
        'i_0': int(data[(data['countyFIPS'].isin(counties)) & (data['Date'] == startdate)]['confirmed'].sum()),
        'r_0': 0,
        'R0': 1.6,
        'gamma': 1/18,  # based on Atkenson's note for now
        't': np.linspace(0, timespan, timespan)
    }
    return params

