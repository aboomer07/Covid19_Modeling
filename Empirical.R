
#import libraries
library(tidyverse)
library(plotly)
library(zoo)
library(viridis)
library(hrbrthemes)
library(gridExtra)
# moving averages

#import datasets
dt <- read.csv("Data/data_covid_FR.csv")
COVID_events_FRANCE <- read.csv("Data/COVID_events_FRANCE_rest.csv", sep=";")
variants <- read.csv("Data/variants_france.csv", sep=";")

#adjust date in the dataframes
dt$date <- (as.POSIXct(paste(dt$date), format = "%Y-%m-%d",tz="UTC"))
COVID_events_FRANCE$date <-  (as.POSIXct(paste(COVID_events_FRANCE$DATE), format = "%d/%m/%Y",tz="UTC"))
variants$semaine <- substr(variants$semaine, 0, 10)
Sys.setlocale("LC_TIME", "en_US"); weekdays(Sys.Date()+0:6)
variants$semaine <- as.Date(variants$semaine)


#################################### Main events france data: daily data#############################

#get rolling values for daily cases 
dt<-dt %>%
  mutate(roll_7=rollmean(conf_j1, k=7, fill=NA))


COVID_events_FRANCE<-inner_join(dt, COVID_events_FRANCE, by="date")
COVID_events_FRANCE<-COVID_events_FRANCE%>%
  select(date, roll_7, event )%>%
  mutate(yend=roll_7, ystart=yend*0.5)

#GGPLOT
ggp <- ggplot(NULL) +                                 
  geom_line(data = dt, aes(x = date, y = roll_7)) +
  geom_segment(data = COVID_events_FRANCE, aes(x = date, y = ystart, xend = date, yend = yend), color = "red") +
  geom_text(data = COVID_events_FRANCE, aes(x = date, y = ystart * 0.9, label = event), color = "red", size=3)
ggp 

#PLOTLY
fig <- plot_ly(dt, x = ~date, y = ~roll_7, type = 'scatter', mode = 'lines')
fig <- fig %>% add_annotations(x = COVID_events_FRANCE$date,
                               y = COVID_events_FRANCE$roll_7,
                               text = COVID_events_FRANCE$event,
                               xref = "x",
                               yref = "y",
                               showarrow = TRUE,
                               size=0.1,
                               arrowhead = 4,
                               arrowsize = .5,
                               ax = -20,
                               ay = -40)%>% 
  layout(xaxis = list(title="Date", titlefont=list(size=12)), yaxis = list(title="Daily cases (7 day moving average)", titlefont=list(size=12)))


fig

#################################### Main events france data: cumulative data#############################

COVID_events_FRANCE<-inner_join(dt, COVID_events_FRANCE, by="date")
COVID_events_FRANCE<-COVID_events_FRANCE%>%
  select(date, conf, event )%>%
  mutate(yend=conf, ystart=yend*0.5)


fig <- plot_ly(dt, x = ~date, y = ~conf, type = 'scatter', mode = 'lines')
fig <- fig %>% add_annotations(x = COVID_events_FRANCE$date,
                               y = COVID_events_FRANCE$conf,
                               text = COVID_events_FRANCE$event,
                               xref = "x",
                               yref = "y",
                               showarrow = TRUE,
                               size=0.1,
                               arrowhead = 4,
                               arrowsize = .5,
                               ax = -20,
                               ay = -40) %>% 
  layout(xaxis = list(title="Date", titlefont=f), yaxis = list(title="Cumulative cases", titlefont=f))

fig


#################################### Variants #############################

#arrange data 
variants<-variants%>%
  group_by(semaine)%>%
  summarise(UK=sum(Nb_susp_501Y_V1)/7, BR=sum(Nb_susp_501Y_V2_3)/7, IND=sum(Nb_susp_IND)/7, ABS=sum(Nb_susp_ABS)/7, total=sum(Nb_tests_PCR_TA_crible)/7)

UK_var<-variants%>%
  select(semaine, UK)%>%
  mutate(Variant="UK")%>%
  rename(case=UK)

BR_var<-variants%>%
  select(semaine, BR)%>%
  mutate(Variant="BR")%>%
  rename(case=BR)

IND_var<-variants%>%
  select(semaine, IND)%>%
  mutate(Variant="IND")%>%
  rename(case=IND)

ABS_var<-variants%>%
  select(semaine, ABS)%>%
  mutate(Variant="ABS")%>%
  rename(case=ABS)

variants_long<-bind_rows(UK_var,BR_var,IND_var,ABS_var)
variants_long$case<-as.numeric(variants_long$case)

"Bresilian and South-Africa variant (P2 and B.1.351 )"
"British variant (5.1.1.7)"
variants_long$Variant <- factor(variants_long$Variant , levels=c("UK", "BR", "IND","ABS") )

#daily cases by variant

daily<-ggplot(variants_long, aes(x=semaine, y=case, fill=Variant)) + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  theme_ipsum() + 
  labs(title="Daily cases of different lineages", x="date", y="daily cases")+
  guides(color=guide_legend("Variant"))  # add guide properties by aesthetic

  
#percentage of each variant 

variants_long <- variants_long  %>%
  group_by(semaine, Variant) %>%
  summarise(n = sum(case)) %>%
  mutate(percentage = (n / sum(n))*100)

# Plot
perc<-ggplot(variants_long, aes(x=semaine, y=percentage, fill=Variant))  + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  theme_ipsum() + 
  labs(title="Percentage of different lineages", x="date", y="%") +
  guides(color=guide_legend("Variant"))  # add guide properties by aesthetic

pdf("variants.pdf")
grid.arrange(daily, perc, ncol=2)
dev.off()
