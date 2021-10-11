# Covid-19 Modeling

This project was done in collaboration with Jacob Pichelmann, Camille Calandre, Nikita Marini, and Luca Poll during my masters program at the Toulouse School of Economics.

## Overview
This project was the year long empirical project during my second year at the Toulouse School of Economics. We used a combination of two different types of epidemiological models, SIR and Time Since Infection, to study the efficacy of a popular public health metric of the novel coronavirus, SARS-COV2 (Covid- 19) under different methodological and disease conditions. To this end, we simulated several outbreaks, taking different assumptions for the underlying infectiousness of the disease, and estimating the rate of transmission in these scenarios. We then forecasted the progression of future cases given our estimators. We found that using a non-parametric estimator of the rate of transmission provides robustness when the shape of the underlying distribution of infectiousness is unknown. We also noted that the assumption of stable conditions makes disease outbreak forecasting in a generalized sense more difficult.

## Data Sources
Data for France was gathered from the Sante Publique France [here](https://www.data.gouv.fr/en/organizations/sante-publique-france/). Data for the US was downloaded from the JHU CSSE COVID-19 Dataset on Github [here](https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data).

## Tools
The paper was built using Latex. The coding was done in R using mainly the EpiEstim, fitdistrplus, KernSmooth, zoo, tseries, and forecast packages. 
