# Identify missing hourly records and insert rows with Temp = NA.
# Run on raw hourly data BEFORE genSeasonList(), which drops the Hour column.
#
# source('R/temperature_gaps.R')
# kob <- read.csv('Ravensburg_hourly_temp.csv')
# checked <- fill_hourly_temperature_gaps(kob, months = c(8:12, 1:5))
# checked$gaps                  # start, end and number of missing hours
# kob_complete <- checked$data # includes NA rows; no interpolation
# checked$inserted_rows         # row indices of newly inserted records
#
# Input: data frame with Date (ISO date string or Date), Hour (0:23), Temp.
# Existing NA temperatures are retained and are not counted as missing records.
# The grid covers whole days from the first through the last observed date.
# months defaults to all months; restrict it explicitly for seasonal datasets.
# Calendar fields Year, Month, Day and JDay are filled when present in the input.
# Other columns are preserved, with NA for newly inserted records.
# UTC is used only for calendar arithmetic: supplied dates/hours are not shifted
# or converted, and every calendar day is treated as having 24 hourly records.
# Output: list(data, gaps, inserted_rows), sorted by date and hour.
fill_hourly_temperature_gaps <- function(data, months = 1:12) {
  if (!is.data.frame(data) || !all(c('Date', 'Hour', 'Temp') %in% names(data)) ||
      anyDuplicated(names(data)) || nrow(data) == 0L) {
    stop('data must be a non-empty data frame with Date, Hour and Temp columns.',
         call. = FALSE)
  }
  if (!(inherits(data$Date, 'Date') || is.character(data$Date))) {
    stop('Date must contain Date values or ISO date strings (YYYY-MM-DD).',
         call. = FALSE)
  }
  dates <- tryCatch(as.Date(data$Date), error = function(e) NULL)
  if (is.null(dates) || anyNA(dates) || any(!is.finite(as.numeric(dates)))) {
    stop('Date contains missing or invalid dates.', call. = FALSE)
  }
  if (!is.numeric(data$Hour) || any(!is.finite(data$Hour)) ||
      any(data$Hour != floor(data$Hour) | data$Hour < 0 | data$Hour > 23)) {
    stop('Hour must contain integer hours from 0 to 23.', call. = FALSE)
  }
  if (!is.numeric(data$Temp)) {
    stop('Temp must be numeric (NA values are allowed).', call. = FALSE)
  }
  if (!is.numeric(months) || !length(months) || any(!is.finite(months)) ||
      any(months != floor(months) | months < 1 | months > 12)) {
    stop('months must contain month numbers from 1 to 12.', call. = FALSE)
  }
  if (any(!as.integer(format(dates, '%m')) %in% months)) {
    stop('months must include every month present in data.', call. = FALSE)
  }
  observed <- as.numeric(as.POSIXct(dates, tz = 'UTC')) + data$Hour * 3600
  if (anyDuplicated(observed)) {
    stop('Duplicate Date/Hour records found; resolve duplicates before filling gaps.',
         call. = FALSE)
  }

  # Start/end at day boundaries to also complete partial first and last days.
  grid <- seq(as.numeric(as.POSIXct(min(dates), tz = 'UTC')),
              as.numeric(as.POSIXct(max(dates), tz = 'UTC')) + 23 * 3600,
              by = 3600)
  timestamps <- as.POSIXct(grid, origin = '1970-01-01', tz = 'UTC')
  keep <- as.integer(format(timestamps, '%m', tz = 'UTC')) %in% months
  grid <- grid[keep]
  timestamps <- timestamps[keep]
  original_rows <- match(grid, observed)
  inserted <- which(is.na(original_rows))
  filled <- data[original_rows, , drop = FALSE]
  rownames(filled) <- NULL

  if (length(inserted)) {
    missing_dates <- as.Date(timestamps[inserted], tz = 'UTC')
    filled$Date[inserted] <- if (inherits(data$Date, 'Date')) missing_dates else
      as.character(missing_dates)
    filled$Hour[inserted] <- as.integer(format(timestamps[inserted], '%H', tz = 'UTC'))
    fields <- c(Year = '%Y', Month = '%m', Day = '%d', JDay = '%j')
    for (field in intersect(names(fields), names(filled))) {
      filled[[field]][inserted] <- as.integer(format(missing_dates, fields[[field]]))
    }
    filled$Temp[inserted] <- NA_real_
    missing <- grid[inserted]
    # Separate gaps across excluded months, even when consecutive in the grid.
    groups <- split(seq_along(missing),
                    cumsum(c(TRUE, diff(missing) != 3600)))
    gaps <- data.frame(
      start = as.POSIXct(vapply(groups, function(i) missing[i[1]], numeric(1)),
                        origin = '1970-01-01', tz = 'UTC'),
      end = as.POSIXct(vapply(groups, function(i) missing[i[length(i)]], numeric(1)),
                      origin = '1970-01-01', tz = 'UTC'),
      n_hours = unname(vapply(groups, length, integer(1))),
      row.names = NULL)
  } else {
    gaps <- data.frame(start = as.POSIXct(character(), tz = 'UTC'),
                       end = as.POSIXct(character(), tz = 'UTC'),
                       n_hours = integer())
  }
  list(data = filled, gaps = gaps, inserted_rows = inserted)
}
