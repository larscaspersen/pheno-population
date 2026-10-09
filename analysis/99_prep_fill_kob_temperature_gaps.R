# Run from the project root: source('analysis/fill_kob_temperature_gaps.R')
# Creates kob_checked (NA-filled data and gap report) and kob_filled (imputed
# temperatures). The input CSV is not overwritten.
source('R/temperature_gaps.R')

latitude_kob <- 47.8 # Approximate Ravensburg latitude; replace with station latitude.
kob_checked <- fill_hourly_temperature_gaps(
  read.csv('Ravensburg_hourly_temp.csv'), months = c(8:12, 1:5))
kob_filled <- kob_checked$data
kob_filled$Temp_imputed <- is.na(kob_filled$Temp)

# Work within seasons so chillR does not reconstruct the omitted summer months.
season <- kob_filled$Year + (kob_filled$Month >= 8)
hour_key <- function(x) paste(x$Year, x$Month, x$Day, x$Hour, sep = '-')
for (yr in unique(season[kob_filled$Temp_imputed])) {
  rows <- which(season == yr)
  weather <- kob_filled[rows, ]
  interpolated <- chillR::interpolate_gaps_hourly(
    hourtemps = weather, latitude = latitude_kob)$weather
  missing <- which(is.na(weather$Temp))
  matched <- match(hour_key(weather[missing, ]), hour_key(interpolated))
  if (anyNA(matched)) stop('Interpolated timestamps could not be matched.')
  # Copy only the missing temperatures, retaining measured values exactly.
  kob_filled$Temp[rows[missing]] <- interpolated$Temp[matched]
}
stopifnot(!anyNA(kob_filled$Temp),
          identical(kob_filled$Temp[!kob_filled$Temp_imputed],
                    kob_checked$data$Temp[!kob_filled$Temp_imputed]))
print(kob_checked$gaps, row.names = FALSE)
message(sum(kob_filled$Temp_imputed), ' temperatures interpolated; ',
        sum(is.na(kob_filled$Temp)), ' NAs remaining.')
