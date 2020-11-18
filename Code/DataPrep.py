import pandas as pd
import json
from itertools import product
from shapely.geometry import Polygon, MultiPolygon
import geopandas as gp
import datetime
import numpy as np


def get_data(data_dir, datasets):

    data_dict = {key: None for key in datasets}

    df = pd.read_csv(data_dir + '/county_data.csv', low_memory=False)
    df['countyFIPS'] = df['countyFIPS'].astype(str).str.zfill(5)

    # melt data into long format, wide format pisses me off
    death_cols = [col for col in df.columns if "#Deaths_" in col]
    case_cols = [col for col in df.columns if "#Cases_" in col]
    mobility_tags = ['_Grocery-Pharmacy', '_Parks', '_Residential',
                     '_Retail-Recreation', '_Transit', '_Workplace']
    mobility_cols = [col for col in df.columns if any(
        [i in col for i in mobility_tags])]

    intervention_cols = ['stay at home', '>50 gatherings', '>500 gatherings',
                         'public schools', 'restaurant dine-in', 'entertainment/gym',
                         'federal guidelines', 'foreign travel ban', 'stay at home rollback',
                         '>50 gatherings rollback', '>500 gatherings rollback',
                         'restaurant dine-in rollback', 'entertainment/gym rollback']

    map_cols = ['countyFIPS', 'PopulationEstimate2018']

    county_info = ['CountyName', 'StateName', 'State', 'lat', 'lon',
                   'STATEFP', 'COUNTYFP', 'POP_LATITUDE', 'POP_LONGITUDE',
                   'CensusRegionCode', 'CensusRegionName', 'CensusDivisionCode',
                   'CensusDivisionName', 'FederalRegionCode', 'SSABeneficiaryCode',
                   'CoreBasedStatAreaCode(CBSA)Metropolitan/Micropolitan2018',
                   'CoreBasedStatAreaName(CBSA)Metropolitan/Micropolitan2018',
                   'CBSAIndicatorCode0=Not,1=Metro,2=Micro2018',
                   'CBSACountyStatusCentralorOutlying2018',
                   'MetropolitanDivisionCode2018', 'MetropolitanDivisionName2018', 'CombinedStatisticalAreaCode2018',
                   'CombinedStatisticalAreaName2018', 'Rural-UrbanContinuumCode2013',
                   'UrbanInfluenceCode2013', 'Economic-DependntTypologyCode2015',
                   'Farming-DependentTypologyCode2015', 'Mining-DependentTypologyCode2015',
                   'Manufacturing-DepTypologyCode2015',
                   'Fed/StGovt-DepdntTypolgyCodeFederal/StateGovernment2015',
                   'RecreationTypolpgyCode2015', 'Nonspecializd-DepTypologyCode2015',
                   'LowEducationTypologyCode2015', 'LowEmploymentTypologyCode2015',
                   'HighPovertyTypologyCode2014', 'PersistentPovrtyTypologyCode2014',
                   'PersistentChildPovTypolCodeRelatedChildren2015',
                   'PopulationLossTypologyCode2015', 'RetirementDestnatnTyplgyCode2015',
                   'BEAEconomicAreaCode2004', 'BEAComponentEconomcAreaCode2004',
                   'BEAEconomicAreaName2004', 'BEAComponentEconomcAreaName2004']

    if 'fips' in datasets:
        fips = pd.read_csv(data_dir + "/FullFips.csv", header=0)
        fips['countyFIPS'] = fips['countyFIPS'].astype(str).str.zfill(5)
        fips = fips.set_index('countyFIPS')
        data_dict['fips'] = fips

    if 'info' in datasets:
        info_df = df[['countyFIPS'] + county_info]
        info_df = info_df.set_index("countyFIPS")
        data_dict['info'] = info_df

    if 'cd' in datasets:
        cd_df = df[map_cols + death_cols + case_cols]
        cd_df = cd_df.set_index(map_cols)
        cd_df.columns = pd.MultiIndex.from_tuples(
            [(i.split("_")[0], i.split("_")[1]) for i in cd_df.columns])
        cd_df = cd_df.stack(level=1).reset_index(
            level=[1, 2]).rename({'level_2': "Date", 'PopulationEstimate2018': 'Pop', '#Deaths': 'Deaths', '#Cases': 'Cases'}, axis=1)
        cd_df['Date'] = pd.to_datetime(cd_df['Date'], format='%m-%d-%Y')
        cd_df['DeathsCapita'] = 100000 * (cd_df['Deaths'] / cd_df['Pop'])
        cd_df['CasesCapita'] = 100000 * (cd_df['Cases'] / cd_df['Pop'])
        data_dict['cd'] = cd_df

    if 'mobility' in datasets:
        mobility = df[['countyFIPS', 'StateName'] + mobility_cols]
        mobility = mobility.set_index(['countyFIPS', 'StateName'])
        mobility.columns = pd.MultiIndex.from_tuples(
            [(i.split("_")[0], i.split("_")[1]) for i in mobility.columns])
        mobility = mobility.stack(level=0).reset_index(
            level=[1, 2]).rename({'level_2': "Date"}, axis=1)
        data_dict['mobility'] = mobility

    if 'measures' in datasets:
        measures = df[['countyFIPS'] + intervention_cols]
        measures = measures.set_index('countyFIPS')
        measures = measures.astype('Int32').applymap(
            lambda x: datetime.date.fromordinal(x) if pd.notnull(x) else np.nan)
        data_dict['measures'] = measures

    if 'geo_map' in datasets:
        with open(data_dir + "/geo_map.json", encoding='ISO-8859-1') as f:
            geo_dict = json.load(f)

        for i in range(len(geo_dict['features'])):
            geo_dict['features'][i]['id'] = \
                geo_dict['features'][i]['properties']['STATE'] + \
                geo_dict['features'][i]['properties']['COUNTY']

        geo_map = {}
        for i in range(len(geo_dict['features'])):
            coords = geo_dict['features'][i]['geometry']['coordinates']
            poly_type = geo_dict['features'][i]['geometry']['type']
            geo_id = geo_dict['features'][i]['id']

            if poly_type == 'Polygon':
                geo_map[geo_id] = Polygon(coords[0])

            elif poly_type == 'MultiPolygon':
                geo_map[geo_id] = MultiPolygon([Polygon(i[0]) for i in coords])

            else:
                geo_map[geo_id] = np.nan

        data_dict['geo_map'] = geo_map

    # if 'full_fips' in datasets:
    #     full_fips = pd.DataFrame(
    #         [i for i in product(fips['countyFIPS'].values,
    #                             df['Date'].unique())],
    #         columns=['countyFIPS', 'Date'])
    #     full_fips['StateCode'] = full_fips.\
    #         apply(lambda x: x['countyFIPS'][:2], axis=1)
    #     full_fips['State'] = full_fips['StateCode'].map(codes)
    #     data_dict['full_fips'] = full_fips

    return data_dict


def get_gp_data(data, geo_map):
    data = data.copy(deep=True)

    if 'countyFIPS' not in data.columns:
        data['geometry'] = data.index.map(geo_map)
    else:
        data['geometry'] = data['countyFIPS'].map(geo_map)

    geo_data = gp.GeoDataFrame(data, geometry=data['geometry'], crs="EPSG:4326")

    return geo_data
