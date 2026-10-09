#load data
library(tidyverse)
library(chillR)

#read temperature and phenology data
source('R/data_io.R')
kob_season <- load_kob_season()
kob_bloom  <- load_kob_bloom()
obs_topaz  <- load_topaz_flowering_cka()
cka_season <- load_cka_season(years = unique(obs_topaz$year))





