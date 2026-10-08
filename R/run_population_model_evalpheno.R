# Package-based alternative to R/run_bud-population-model.R.
# Source this file from the project root; no local C++ compilation is needed.
#
# source('R/run_population_model_evalpheno.R')
# params <- load_population_params('parameters/topaz_population_normal.yaml')
# result <- run_population_model(forcing_obs$temp_df, params,
#                                jday_cut = forcing_obs$jday_cut,
#                                stop_at_zc = FALSE)
# save_population_params(params, 'parameters/topaz_population_updated.yaml')
#
# Uses evalpheno::phenoflex_population() to retain the manuscript output fields
# (bloomindex, x, y, z, exp, z_effec, exp_force_inc, yc_pop, zc_pop).
# The modular predict_population_phenology() API uses different output fields
# and forcing-time conventions, so it is not substituted here.

.population_parameter_vector <- function(par, convert = TRUE) {
  kinetic <- c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1',
               'Tf', 'Tc', 'Tb', 'slope')
  characteristic <- c('yc', 'zc', 's1', 'Tu', 'theta_star', 'theta_c',
                      'tau', 'pie_c', 'Tf', 'Tc', 'Tb', 'slope')
  if (is.list(par)) par <- unlist(par)
  if (!is.numeric(par) || length(par) != 12L || any(!is.finite(par))) {
    stop('par must contain twelve finite numeric PhenoFlex parameters.', call. = FALSE)
  }
  if (is.null(names(par))) {
    # Unnamed vectors retain the original helper's kinetic parameter order.
    names(par) <- kinetic
  }
  # The manuscript's characteristic YAML uses pi_c; evalpheno uses pie_c.
  names(par)[names(par) == 'pi_c'] <- 'pie_c'
  if (anyNA(names(par)) || anyDuplicated(names(par))) {
    stop('Parameter names must be unique and non-missing.', call. = FALSE)
  }
  if (setequal(names(par), kinetic)) return(par[kinetic])
  if (!setequal(names(par), characteristic)) {
    stop('Use the twelve kinetic or characteristic PhenoFlex parameter names.',
         call. = FALSE)
  }
  par <- par[characteristic]
  if (!convert) return(par)
  par <- evalpheno::convert_parameters(par, failure_return = 'NA')
  if (!is.numeric(par) || length(par) != 12L || any(!is.finite(par))) {
    stop('Characteristic-to-kinetic parameter conversion failed.', call. = FALSE)
  }
  stats::setNames(par, kinetic)
}

.normalize_population_params <- function(params) {
  if (!is.list(params) || is.null(params$par)) {
    stop('params must be a configuration list containing par.', call. = FALSE)
  }
  # Keep the original parameterization for saving; convert only when running.
  params$par <- .population_parameter_vector(params$par, convert = FALSE)
  defaults <- list(yc_sd = 0, zc_sd = 0, dist_chill = 'normal',
                   dist_heat = 'normal', add_par = NULL, seed = 12345, n_pop = 100)
  if (is.null(params$distribution)) params$distribution <- list()
  if (!is.list(params$distribution)) {
    stop('distribution must be a list of population settings.', call. = FALSE)
  }
  params$distribution <- utils::modifyList(defaults, params$distribution,
                                            keep.null = TRUE)
  params
}

# Read the existing manuscript YAML format, preserving description, stage
# thresholds, forcing adjustment and any other configuration metadata.
load_population_params <- function(path) {
  if (!requireNamespace('yaml', quietly = TRUE)) {
    stop("Install 'yaml' to read population parameter files.", call. = FALSE)
  }
  .normalize_population_params(yaml::read_yaml(path))
}

# Save named parameters as a YAML mapping, rather than an unnamed sequence.
# Characteristic parameters remain characteristic, allowing a lossless reload.
save_population_params <- function(params, path) {
  if (!requireNamespace('yaml', quietly = TRUE)) {
    stop("Install 'yaml' to save population parameter files.", call. = FALSE)
  }
  params <- .normalize_population_params(params)
  params$par <- as.list(params$par)
  yaml::write_yaml(params, path)
  invisible(path)
}

# Accept either a loaded configuration or a YAML path. Additional named
# arguments override package settings (e.g. n, seed, adjust_zc, population,
# max_days_forcing or forcing_temperature).
run_population_model <- function(temp_df, params, jday_cut = NULL,
                                  basic_output = FALSE, stop_at_zc = TRUE, ...) {
  if (is.character(params) && length(params) == 1L) {
    params <- load_population_params(params)
  }
  params <- .normalize_population_params(params)
  distribution <- params$distribution
  adjustment <- params$scale_yc_budbreak
  if (is.null(adjustment)) adjustment <- 1
  args <- list(par = params$par, temp_df = temp_df, jday_cut = jday_cut,
               yc_sd = distribution$yc_sd, zc_sd = distribution$zc_sd,
               n = distribution$n_pop, dist_chill = distribution$dist_chill,
               dist_heat = distribution$dist_heat,
               add_distpar = distribution$add_par, seed = distribution$seed,
               adjust_zc = adjustment, basic_output = basic_output,
               stop_at_zc = stop_at_zc)
  overrides <- list(...)
  if (length(overrides) &&
      (is.null(names(overrides)) || any(names(overrides) == '') ||
       anyDuplicated(names(overrides)))) {
    stop('Additional model arguments must have unique names.', call. = FALSE)
  }
  do.call(helper_run_pop_model,
          utils::modifyList(args, overrides, keep.null = TRUE))
}

# Compatibility entry point for existing manuscript helper_run_pop_model()
# calls. Source this alternative after other helpers to select the package
# implementation. The package validates weather, cuts and sampled requirements.
helper_run_pop_model <- function(par, yc_sd = 0, zc_sd = 0, jday_cut = NULL,
                                 temp_df, n = 100, adjust_zc = 1,
                                 basic_output = FALSE, stop_at_zc = TRUE,
                                 dist_chill = 'normal', dist_heat = 'normal',
                                 add_distpar = NULL, seed = 12345,
                                 max_days_forcing = 50, forcing_temperature = 23,
                                 population = NULL) {
  if (!requireNamespace('evalpheno', quietly = TRUE)) {
    stop("Install 'evalpheno' to run the population model.", call. = FALSE)
  }
  evalpheno::phenoflex_population(
    temp_df = temp_df, par = .population_parameter_vector(par),
    yc_sd = yc_sd, zc_sd = zc_sd, jday_cut = jday_cut, n = n,
    adjust_zc = adjust_zc, basic_output = basic_output, stop_at_zc = stop_at_zc,
    dist_chill = dist_chill, dist_heat = dist_heat, add_distpar = add_distpar,
    seed = seed, max_days_forcing = max_days_forcing,
    forcing_temperature = forcing_temperature, population = population)
}
