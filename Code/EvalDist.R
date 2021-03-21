# Objective : Evaluate fitted distributions on serial interval data
# Created by: jacobpichelmann
# Created on: 21.03.21

source(paste0(getwd(), "/Code/SimFunc.R"))


### Specify true distributions

### Simulate serial interval
samp_pois_plot <- function(samps, study_len){
  data.frame(table(samps)) %>%
    mutate(samps = as.numeric(as.character(samps))) %>%
    ggplot(aes(x=samps, y=Freq)) + geom_bar(stat = "identity") + xlab("days") + ylab("count") +
      scale_x_continuous(limits = c(0,study_len)) + theme_minimal()
}

samp_pois_plot(samps = samps, study_len = study_len)

### Fit gamma

### Evaluate

og_gamma <- gen_distribution(15, sim_mu, sim_var, "gamma")
new_gamma <- gen_distribution(15, mean(params[params$Distribution == 'gamma', 'True_a']), mean(params[params$Distribution == 'gamma', 'True_b']),"gamma")
og_weibull <- gen_distribution(15, a_weibull, b_weibull, "weibull")
new_weibull<- gen_distribution(15, mean(params[params$Distribution == 'weibull', 'True_a']), mean(params[params$Distribution == 'weibull', 'True_b']),"gamma")
og_norm <- gen_distribution(15, a_norm, b_norm, "norm")
new_norm <- gen_distribution(15, mean(params[params$Distribution == 'norm', 'True_a']), mean(params[params$Distribution == 'norm', 'True_b']),"gamma")
og_lnorm <- gen_distribution(15, a_lnorm, b_lnorm, "lnorm")
new_lnorm <- gen_distribution(15, mean(params[params$Distribution == 'lnorm', 'True_a']), mean(params[params$Distribution == 'lnorm', 'True_b']),"gamma")

jpeg("DistCompare.jpg")
layout(matrix(1:4, nrow = 2, ncol=2))
plot(og_gamma, type="l")
lines(new_gamma, type="l", col="green")
title("Gamma")

plot(og_weibull, type="l")
lines(new_weibull, type="l", col="red")
title("Weibull")

plot(og_norm, type="l")
lines(new_norm, type="l", col="blue")
title("Normal")

plot(og_lnorm, type="l")
lines(new_lnorm, type="l", col="orange")
title("Log-normal")
dev.off()