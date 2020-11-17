# Additional data functions

import pandas as pd


# compute empirical R0 based on RKI method
# https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121343/r0_technical_note_v3.4.pdf

def get_R0(country, df):
    df.loc[df['Country/Region'] == country, 'R0'] = df.confirmed.rolling(window=4).sum() / df.confirmed.rolling(
        window=4).sum().shift(4)
    df = df[df['Country/Region'] == country]
    return df
