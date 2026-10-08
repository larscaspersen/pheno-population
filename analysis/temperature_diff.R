#temperature graph
library(tidyverse)


cka <- read.csv('experimental/cka_clean.csv') |> 
  filter(Year >= 1998 & Year <= 2011,
         Month %in% c(1:5, 10:12)) |> 
  mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by(Month, Year) |> 
  summarise(mean = mean(Tmean),
            median = median(Tmean),
            q_05 = quantile(Tmean, 0.05),
            q_95 = quantile(Tmean, 0.95)) |> 
  ungroup() |> 
  mutate(loc = 'Rhineland (Klein-Altendorf)')

kob <- read.csv('experimental/Ravensburg_hourly_temp.csv') |> 
  filter(Year >= 2003 & Year <= 2022,
         Month %in% c(1:5, 10:12)) |> 
  #mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by(Month, Year) |> 
  summarise(mean = mean(Temp),
            median = median(Temp),
            q_05 = quantile(Temp, 0.05),
            q_95 = quantile(Temp, 0.95)) |> 
  ungroup() |> 
  mutate(loc = 'Lake Constance (Bavendorf)')

kob |> 
  rbind(cka) |> 
  mutate(mon_abb = month.abb[Month],
         mon_abb = factor(mon_abb, levels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May'))) |> 
  ggplot(aes(x=mon_abb)) +
  geom_boxplot(aes(y = mean, fill = loc))


read.csv('experimental/cka_clean.csv') |> 
  filter(Year >= 1998 & Year <= 2011,
         Month %in% c(1:5, 10:12)) |> 
  mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by() |> 
  summarise(mean = mean(Tmean),
            median = median(Tmean),
            q_05 = quantile(Tmean, 0.05),
            q_95 = quantile(Tmean, 0.95)) |> 
  ungroup() |> 
  mutate(loc = 'Rhineland (Klein-Altendorf)')

read.csv('experimental/Ravensburg_hourly_temp.csv') |> 
  filter(Year >= 2003 & Year <= 2022,
         Month %in% c(1:5, 10:12)) |> 
  #mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by() |> 
  summarise(mean = mean(Temp),
            median = median(Temp),
            q_05 = quantile(Temp, 0.05),
            q_95 = quantile(Temp, 0.95)) |> 
  ungroup() |> 
  mutate(loc = 'Lake Constance (Bavendorf)')

## by month

cka_mon <- read.csv('experimental/cka_clean.csv') |> 
  filter(Year >= 1998 & Year <= 2011,
         Month %in% c(1:5, 10:12)) |> 
  mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by(Month) |> 
  summarise(cka = mean(Tmean))

kob_mon <- read.csv('experimental/Ravensburg_hourly_temp.csv') |> 
  filter(Year >= 2003 & Year <= 2022,
         Month %in% c(1:5, 10:12)) |> 
  #mutate(Tmean = (Tmin + Tmax) / 2) |> 
  group_by(Month) |> 
  summarise(kob = mean(Temp)) 

merge(cka_mon, kob_mon, by = 'Month') |> 
  mutate(diff = cka - kob)
