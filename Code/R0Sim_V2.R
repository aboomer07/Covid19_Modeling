# Title     : Nikita's R0 Estimation
# Created by: jacobpichelmann
# Created on: 06.02.21

# import libraries
library(EpiEstim)
library(R0)
library(dplyr)
library(ggplot2)
library(accelerometry)
library(tidyverse)

imppath <- paste0(getwd(), '/Data/')

# import data and focus on the US
actual_data <- read.csv2(paste0(imppath, 'EuropeCovidData.csv'))
actual_data <- subset(actual_data, Country.Region == 'US')

# drop first 80 observations (2M) to account for testing catching up
actual_data <- tail(actual_data, -80)


# define serial interval distribution taken from Cereda et al (2020)
gamma_x <- seq(1, 100, 1)
gamma_y <- dgamma(gamma_x, 6.6, 1.1)
plot(gamma_y, type = "l")


# Define time varying effective reproduction number
# important: here you need to specify manually when you would like the true R0 to change

for (t in 1:length(actual_data$date)){
if (actual_data$date[t] <= as.Date("2020-06-11")){
actual_data$R_val[t] <- 1.5
}
else if (actual_data$date[t] > as.Date("2020-06-11")) {
actual_data$R_val[t] <- 1.25
}
}

# cannot do it in one loop cause i suck
for (t in 1:length(actual_data$date)){
if (actual_data$date[t] > as.Date("2020-09-11")) {
actual_data$R_val[t] <- 1.4
}
}

#Generate outbreak accoring to E[I_t]=R+t*sum^{t}{s=1}I{t-s}w_s
#we use a window of 5 days

actual_data$infective <- 0
actual_data$simulation <- 0

for (t in 6:length(actual_data$date)) {
actual_data$infective[t] <- actual_data$confirmed[t-5]*gamma_y[5] +
actual_data$confirmed[t-4]*gamma_y[4] + actual_data$confirmed[t-3]*gamma_y[3] +
actual_data$confirmed[t-2]*gamma_y[2] + actual_data$confirmed[t-1]*gamma_y[1]

actual_data$simulation[t] <- actual_data$R_val[t]*actual_data$infective[t]
}


#Estimate R via the function in the EpiEstim package
est_r0 <- estimate_R(actual_data$simulation, method = 'parametric_si', config = make_config(list(mean_si = 6.6, std_si = 1.1)))
plot(est_r0$R$`Mean(R)`, type ="l", ylim = c(0,4))
lines(actual_data$R_val, col = "red")