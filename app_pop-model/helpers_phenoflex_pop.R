# App bootstrap: use the shared manuscript helpers, never compile the local kernel.
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(chillR)

app_project_root <- normalizePath(
  if (file.exists('R/run_population_model_evalpheno.R')) '.' else '..',
  winslash = '/', mustWork = TRUE)
for (file in c('data_io.R', 'run_population_model_evalpheno.R',
               'plots_forcing.R', 'plots_flowering.R', 'plots_misc.R')) {
  source(file.path(app_project_root, 'R', file), local = TRUE)
}

app_parameter_files <- list.files(file.path(app_project_root, 'parameters'),
                                  pattern = '\\.(yaml|yml)$', ignore.case = TRUE)
app_default_parameter_file <- 'topaz_population_normal.yaml'

# Populate the controls from either YAML parameterization without a round trip
# through the other one. Inactive chill controls are left empty.
app_parameter_inputs <- function(config) {
  par <- config$par
  names(par)[names(par) == 'pie_c'] <- 'pi_c'
  if ('theta_star' %in% names(par)) {
    par[c('theta_star', 'theta_c')] <- par[c('theta_star', 'theta_c')] - 273.15
  }
  values <- as.list(stats::setNames(rep(NA_real_, 8L),
    c('E0', 'E1', 'A0', 'A1', 'theta_star', 'theta_c', 'tau', 'pi_c')))
  values <- utils::modifyList(values, as.list(par))
  distribution <- config$distribution
  skew <- distribution$add_par
  if (is.null(skew)) skew <- c(0, 0)
  adjustment <- config$scale_yc_budbreak
  if (is.null(adjustment)) adjustment <- 1
  c(values, list(yc_sd = distribution$yc_sd, zc_sd = distribution$zc_sd,
    adjust_zc = adjustment, dist_chill = distribution$dist_chill,
    dist_heat = distribution$dist_heat, skew_chill = skew[1], skew_heat = skew[2],
    n_pop = distribution$n_pop,
    seed = if (is.null(distribution$seed)) NA_real_ else distribution$seed))
}

# The shared forcing loader reads paths relative to the project. Shiny starts
# with the app directory as its working directory; restore it after each read.
app_prepare_obs_data <- function(sheet, start_yday, end_yday) {
  previous_dir <- setwd(app_project_root)
  on.exit(setwd(previous_dir), add = TRUE)
  helper_prepare_obs_data(sheet, start_yday, end_yday)
}

kob_season <- load_kob_season(
  path = file.path(app_project_root, 'data/Ravensburg_hourly_temp_fixed.csv'))
kob_bloom <- load_kob_bloom(
  path = file.path(app_project_root, 'data/Ravensburg_bloom_dates.csv'))
kob_bloom$year <- as.character(kob_bloom$year)

app_select_years <- function(selection, available = names(kob_season)) {
  validation <- c('2007', '2009', '2010', '2016')
  years <- as.character(selection)
  if ('all' %in% selection) return(as.character(available))
  if ('calibration' %in% selection) years <- c(years, setdiff(available, validation))
  if ('validation' %in% selection) years <- c(years, validation)
  as.character(available)[as.character(available) %in% years]
}

app_predict_orchard <- function(seasons, params, budbreak = FALSE) {
  if (budbreak) {
    # Match the manuscript's orchard convention: change the mean heat
    # requirement, retaining its specified population standard deviation.
    params$par['zc'] <- params$par['zc'] * params$scale_yc_budbreak
  }
  bind_rows(lapply(names(seasons), function(year) {
    season <- seasons[[year]]
    result <- run_population_model(season, params, basic_output = TRUE)
    data.frame(year = year, value = vapply(result$bloomindex,
      helper_bloomint_to_jday, numeric(1), x = season))
  }))
}
