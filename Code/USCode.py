################################################################################
# Import Libraries and Other Scripts for Functions
################################################################################
import os
from DataPrep import *
from Plots import *
################################################################################
# Finish Import Libraries
################################################################################

################################################################################
# Sort of out location of all the folders
################################################################################
if os.getcwd().split("/")[-1] == 'Code':
    data_dir = os.path.abspath("..") + "/Data"

    if not os.path.isdir(os.path.abspath("..") + '/Output'):
        os.mkdir(os.path.abspath("..") + '/Output')
    output_dir = os.path.abspath("..") + '/Output'

    report_dir = os.path.abspath("..") + '/Report'

################################################################################
# Define the set of datasets to pull into dictionary from DataPrep script

#       'df' = full county level dataset with all columns
#       'fips' = mapping of all fips codes to their states
#       'info' = county specific time invariant data
#       'cd' = cumulative cases and deaths per county per day
#       'mobility' = google moblity data per day per county, diff categories
#       'measures' = the policy interventions per county
#       'geo_map' = the mapping from fips codes to polygons for geopandas
#       'full_fips' = the full set of fips codes and dates in case data is null
################################################################################

datasets = ['fips', 'cd', 'geo_map', 'measures', 'mobility']
data_dict = get_data(data_dir, datasets)

################################################################################
# Define the set of plotting params for making choropleth

#       'size' = tuple of figure size
#       'col' = the value column for the choropleth
#       'cmap' = the color scale to use
#       'legend' = boolean for whether to include a legend
#       'title' = the text string of the plot title
#       'title_font' = a font dictionary inclding fontsize and fontweight
#       'date_bar' = a boolean for whether to convert bar labels to string dates
#       'file' = the file path to write graph to
################################################################################

plot_params = {
    'size': (50, 20),
    'col': 'Date',
    'cmap': 'viridis',
    'legend': True,
    'title': 'Week of 100 Cumulative Confirmed Cases Per 100,000 Pop',
    'title_font': {'fontsize': 50, 'fontweight': 15},
    'bar_font': {'fontsize': 30},
    'date_bar': True,
    'file': output_dir + "/mincases.png"
}

################################################################################
# Plot the first date that each county hit 100 cases per 100,000
################################################################################

week = data_dict['cd'].groupby([data_dict['cd'].index, pd.Grouper(
    key='Date', freq='W-MON')])[['CasesCapita', 'DeathsCapita']].\
    agg(np.nanmean).reset_index(level=1)
week['MinDist'] = week['CasesCapita'].apply(
    lambda x: x - 100 if x - 100 >= 0 else np.nan)

first = week.reset_index().set_index('Date').\
    groupby('countyFIPS')['MinDist'].idxmin().\
    reset_index().rename({'MinDist': 'Date'}, axis=1).set_index('countyFIPS')
first['Date'] = first['Date'].apply(
    lambda x: x.timestamp() if pd.notnull(x) else np.nan)
first['State'] = first.index.map(data_dict['fips']['State'])
first = first[~first['State'].isin(
    ['AK', 'HI', 'VI', 'PR'])].set_index('countyFIPS')

first = get_gp_data(first, data_dict['geo_map'])
make_choropleth(first, plot_params)

################################################################################
# Plot the first date that each county hit 1 death per 100,000
################################################################################

plot_params['title'] = 'Week of 1 Cumulative Confirmed Death Per 100,000 Pop'
plot_params['file'] = output_dir + "/mindeaths.png"
week['MinDist'] = week['DeathsCapita'].apply(
    lambda x: x - 1 if x - 1 >= 0 else np.nan)

first = week.reset_index().set_index('Date').\
    groupby('countyFIPS')['MinDist'].idxmin().\
    reset_index().rename({'MinDist': 'Date'}, axis=1).set_index('countyFIPS')
first['Date'] = first['Date'].apply(
    lambda x: x.timestamp() if pd.notnull(x) else np.nan)
first['State'] = first.index.map(data_dict['fips']['State'])
first = first[~first['State'].isin(
    ['AK', 'HI', 'VI', 'PR'])].set_index('countyFIPS')

first = get_gp_data(first, data_dict['geo_map'])
make_choropleth(first, plot_params)

################################################################################
# Plot the date of stay at home measures for each county
################################################################################

plot_params['col'] = 'stay at home'
plot_params['title'] = 'Stay at Home Issuance Date by US County'
plot_params['file'] = output_dir + '/stayhome.png'

stayhome = data_dict['measures'][['stay at home', 'stay at home rollback']]
stayhome = stayhome.applymap(lambda x: (x - datetime.date(1970, 1, 1)).
                             total_seconds() if pd.notnull(x) else np.nan)
stayhome['stayhome_duration'] = (stayhome['stay at home rollback'] -
                                 stayhome['stay at home']) / (3600 * 24)
stayhome['State'] = stayhome.index.map(data_dict['fips']['State'])
stayhome = stayhome[~stayhome['State'].isin(['AK', 'HI', 'VI', 'PR'])]

stayhome = get_gp_data(stayhome, data_dict['geo_map'])
make_choropleth(stayhome, plot_params)

################################################################################
# Plot the date of stay at home rollback measures for each county
################################################################################

plot_params['col'] = 'stay at home rollback'
plot_params['title'] = 'Stay at Home Rollback Date by US County'
plot_params['file'] = output_dir + '/stayhome_end.png'
make_choropleth(stayhome, plot_params)

################################################################################
# Plot the duration of the stay at home orders per county
################################################################################

plot_params['col'] = 'stayhome_duration'
plot_params['title'] = 'Stay at Home Order Duration by US County'
plot_params['file'] = output_dir + '/stayhome_length.png'
plot_params['date_bar'] = False
make_choropleth(stayhome, plot_params)

################################################################################
# Plot the current cases per 100,000
################################################################################

current = week.groupby(week.index).agg(lambda x: x.iloc[x.Date.argmax()])
current['State'] = current.index.map(data_dict['fips']['State'])
current = current[~current['State'].isin(['AK', 'HI', 'VI', 'PR'])]
current = get_gp_data(current, data_dict['geo_map'])

plot_params['col'] = 'CasesCapita'
plot_params['title'] = 'Current Cases per 100,000 Pop'
plot_params['file'] = output_dir + '/curr_cases.png'
make_choropleth(current, plot_params)

################################################################################
# Plot the current deaths per 100,000
################################################################################

plot_params['col'] = 'DeathsCapita'
plot_params['title'] = 'Current Deaths per 100,000 Pop'
plot_params['file'] = output_dir + '/curr_deaths.png'
make_choropleth(current, plot_params)
