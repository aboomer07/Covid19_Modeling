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
    os.chdir("..")
curr_dir = os.getcwd()
data_dir = curr_dir + "/Data/"

df = pd.read_csv(data_dir + 'EuropeCovidData.csv', sep=';')


def get_params(country, startdate, data, timespan):
    params = {
        'N': int(wb.get_series('SP.POP.TOTL', mrv=1, country=coco.convert(country, to='ISO3'))),
        'i_0': int(data[(data['Country/Region'] == country) &
                       (data['date'] == startdate)]['confirmed']),
        'r_0': 0,  # need more data to specify
        'R0': 1.6,  # derive this from data
        'gamma': 0.1,  # derive this from data
        't': np.linspace(0, timespan, timespan)
    }
    return params
