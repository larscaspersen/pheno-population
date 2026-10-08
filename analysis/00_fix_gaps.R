#fix kob season
raven <- read.csv('Ravensburg_hourly_temp.csv') 
source('R/temperature_gaps.R')

raven_compl <- raven |> fill_hourly_temperature_gaps( months = c(8:12, 1:5))

raven_filled <- chillR::interpolate_gaps_hourly(raven_compl$data, latitude = 47.782)
write.csv(raven_filled$weather, file = 'Ravensburg_hourly_temp_fixed.csv', row.names = FALSE)