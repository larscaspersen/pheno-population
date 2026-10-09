library(shiny)
library(bslib)
library(patchwork)
options(shiny.sanitize.errors = FALSE)

# Support both source('app_pop-model/app.R') and shiny::runApp('app_pop-model').
helper_path <- if (file.exists('helpers_phenoflex_pop.R')) {
  'helpers_phenoflex_pop.R'
} else {
  'app_pop-model/helpers_phenoflex_pop.R'
}
source(helper_path, local = TRUE)


default_config <- load_population_params(file.path(
  app_project_root, 'parameters', app_default_parameter_file))
default_params <- app_parameter_inputs(default_config)

# Parameter controls and model plots.
ui <- page_sidebar(
  # App title ----
  title = "PhenoFlex Population Model - Forcing Experiment",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    
    checkboxInput(
      inputId = "plot_forcing",
      label = "Show Forcing Plot",
      value = TRUE
    ),
    checkboxInput(
      inputId = "plot_flowering",
      label = "Show Flowering Plot",
      value = FALSE
    ),
    checkboxInput(
      inputId = "plot_budbreak_orchard",
      label = "Show predicted budbreak in orchard",
      value = FALSE
    ),
    checkboxInput(
      inputId = "plot_tempresponse",
      label = "Show Temperature Response Plot",
      value = FALSE
    ),
    numericInput(
      inputId = "yc",
      label = "yc\n(mean chill requirement)",
      min = 10,
      max = 80,
      value = default_params$yc
    ),
    numericInput(
      inputId = "yc_sd",
      label = "yc sd\n(standard deviation of chill requirement)",
      min = 0,
      max = 20,
      value = default_params$yc_sd
    ),
    numericInput(
      inputId = "zc",
      label = "zc\n(mean of heat requirement)",
      min = 100,
      max = 500,
      value = default_params$zc
    ),
    numericInput(
      inputId = "zc_sd",
      label = "zc sd\n(standard deviation of heat requirement)",
      min = 0,
      max = 20,
      value = default_params$zc_sd
    ),
    numericInput(
      inputId = "s1",
      label = "s1\n(slope of transition function between chill and heat accumulation)",
      min = 0.01,
      max = 1.5,
      value = default_params$s1
    ),
    numericInput(
      inputId = "adjust_zc",
      label = "conversion heat flowering to budbreak",
      min = 0.1,
      max = 1,
      value = default_params$adjust_zc
    ),
    selectizeInput(
      inputId = "dist_chill",
      label = "Chill requirement distribution",
      choices = c('normal' = 'normal',
                  'normal skewed' = 'normal_skewed'),
      selected = default_params$dist_chill,
      multiple = FALSE
    ),
    conditionalPanel(
      condition = "input.dist_chill == 'normal_skewed'",
      numericInput('skew_chill', 'Chill distribution skewness', value = default_params$skew_chill)
    ),
    selectizeInput(
      'dist_heat', 'Heat requirement distribution',
      choices = c('normal' = 'normal', 'normal skewed' = 'normal_skewed'),
      selected = default_params$dist_heat, multiple = FALSE
    ),
    conditionalPanel(
      condition = "input.dist_heat == 'normal_skewed'",
      numericInput('skew_heat', 'Heat distribution skewness', value = default_params$skew_heat)
    ),
    numericInput('n_pop', 'Number of buds', value = default_params$n_pop,
                 min = 1, step = 1),
    numericInput('seed', 'Random seed (leave empty for random draws)',
                 value = default_params$seed, min = 0, step = 1),
    conditionalPanel(
      condition = "output.chill_parameterization == 'kinetic'",
      numericInput('E0', 'E0 (activation energy of precursor formation)',
                   value = default_params$E0, min = 0),
      numericInput('E1', 'E1 (activation energy of precursor destruction)',
                   value = default_params$E1, min = 0),
      numericInput('A0', 'A0 (rate coefficient of precursor formation)',
                   value = default_params$A0, min = 0),
      numericInput('A1', 'A1 (rate coefficient of precursor destruction)',
                   value = default_params$A1, min = 0)
    ),
    conditionalPanel(
      condition = "output.chill_parameterization == 'characteristic'",
      numericInput(
        inputId = "theta_star",
        label = "theta_star\n(optimal temperature in °C for chill accumulation)",
        min = 5,
        max = 8,
        value = default_params$theta_star
      ),
      numericInput(
        inputId = "theta_c",
        label = "theta_c\n(critical temperature in °C for chill accumulation)",
        min = 12,
        max = 15,
        value = default_params$theta_c
      ),
      numericInput(
        inputId = "tau",
        label = "tau\n(time interval for chill accumulation under optimal conditions)",
        min = 16,
        max = 48,
        value = default_params$tau
      ),
      numericInput(
        inputId = "pi_c",
        label = "pi_c\n(time interval leading to chill negation)",
        min = 24,
        max = 40,
        value = default_params$pi_c
      )
    ),
    numericInput(
      inputId = "Tf",
      label = "Tf\n(temperature for transition from PDBF to DBF)",
      min = 0,
      max = 10,
      value = default_params$Tf
    ),
    numericInput(
      inputId = "slope",
      label = "slope\n(slope of transition function from PDBF to DBF)",
      min = 0.1,
      max = 15,
      value = default_params$slope
    ),
    numericInput(
      inputId = "Tb",
      label = "Tb\n(base temperature heat accumulation)",
      min = 0,
      max = 10,
      value = default_params$Tb
    ),
    numericInput(
      inputId = "Tu",
      label = "Tu\n(optimal temperature heat accumulation)",
      min = 20,
      max = 35,
      value = default_params$Tu
    ),
    numericInput(
      inputId = "Tc",
      label = "Tc\n(critical temperature heat accumulation)",
      min = 30,
      max = 40,
      value = default_params$Tc
    ),
    selectInput(
      inputId = 'par_select',
      label = 'Parameter YAML file',
      choices = app_parameter_files,
      selected = app_default_parameter_file,
      multiple = FALSE,
      selectize = TRUE,
      width = NULL,
      size = NULL
    ),
    actionButton(
      inputId = "reload_parameters",
      label = "Reload parameters from YAML",
      class = "btn-secondary"
    ),
    selectizeInput(
      inputId = "forcing_exp",
      label = "Choose which forcing experiment to plot",
      choices = c('Kanzi, 2020' = 'K_2020_term+spur',
                  'Kanzi, 2022' = 'K_2022_term+spur',
                   'Topaz, 2020' = 'T_2020_term+spur',
                   'Topaz, 2022' = 'T_2022_term+spur'),
      selected = 'T_2022_term+spur',
      multiple = FALSE
    ),
    numericInput(
      inputId = "exp_start_yday",
      label = "Julian Day of earliest cutting experiment",
      min = 1,
      max = 365,
      value = 300
    ),
    numericInput(
      inputId = "exp_end_yday",
      label = "Julian Day of latest cutting experiment",
      min = 1,
      max = 365,
      value = 55
    ),
    selectizeInput(
      inputId = "years_bloom",
      label = "Choose years for bloom prediction",
      choices = c('all', 'calibration', 'validation', as.character(2004:2022)),
      selected = 'all',
      multiple = TRUE
    )
  ),
  # Combined output for selected plots.
  plotOutput(outputId = "pop_plot")
)

# Model and plots use the shared helpers and the evalpheno population adapter.
server <- function(input, output, session) {
  # Show useful validation messages while users edit parameters.
  model_result <- function(expr) {
    tryCatch(force(expr), error = function(error) {
      validate(need(FALSE, conditionMessage(error)))
    })
  }

  loaded_params <- reactive({
    req(input$par_select)
    input$reload_parameters
    validate(need(input$par_select %in% app_parameter_files,
                  'Select a parameter file from the parameters directory.'))
    model_result(load_population_params(file.path(
      app_project_root, 'parameters', input$par_select)))
  })

  output$chill_parameterization <- reactive({
    if ('E0' %in% names(loaded_params()$par)) 'kinetic' else 'characteristic'
  })
  outputOptions(output, 'chill_parameterization', suspendWhenHidden = FALSE)

  observeEvent(loaded_params(), {
    values <- app_parameter_inputs(loaded_params())
    for (name in names(values)) {
      freezeReactiveValue(input, name)
      if (name %in% c('dist_chill', 'dist_heat')) {
        updateSelectizeInput(session, name, selected = values[[name]])
      } else {
        updateNumericInput(session, name, value = values[[name]])
      }
    }
  })

  par <- reactive({
    req(input$yc, input$zc, input$s1, input$Tu, input$Tf, input$Tc,
        input$Tb, input$slope)
    values <- c(yc = input$yc, zc = input$zc, s1 = input$s1, Tu = input$Tu,
                Tf = input$Tf, Tc = input$Tc, Tb = input$Tb, slope = input$slope)
    if ('E0' %in% names(loaded_params()$par)) {
      req(input$E0, input$E1, input$A0, input$A1)
      values <- c(values, E0 = input$E0, E1 = input$E1,
                  A0 = input$A0, A1 = input$A1)
      validate(need(input$E0 > 0 && input$E1 > input$E0 &&
                      input$A0 > 0 && input$A1 > 0,
                    'Kinetic chill parameters require 0 < E0 < E1 and positive rate coefficients.'))
    } else {
      req(input$theta_star, input$theta_c, input$tau, input$pi_c)
      values <- c(values, theta_star = input$theta_star + 273.15,
                  theta_c = input$theta_c + 273.15,
                  tau = input$tau, pi_c = input$pi_c)
      validate(need(input$theta_star < input$theta_c && input$tau > 0 && input$pi_c > 0,
                    'Chill parameters require theta_star < theta_c and positive time intervals.'))
    }
    validate(
      need(all(is.finite(values)), 'Enter finite parameter values.'),
      need(input$yc > 0 && input$zc > 0 && input$s1 > 0 && input$slope > 0,
           'Requirements and slopes must be positive.'),
      need(input$Tb < input$Tu && input$Tu < input$Tc,
           'Heat temperatures must satisfy Tb < Tu < Tc.'))
    model_result(.population_parameter_vector(values))
  })

  params <- reactive({
    req(input$yc_sd, input$zc_sd, input$adjust_zc, input$dist_chill,
        input$dist_heat, input$skew_chill, input$skew_heat, input$n_pop)
    seed <- input$seed
    if (is.null(seed) || is.na(seed)) seed <- NULL
    validate(
      need(is.finite(input$yc_sd) && input$yc_sd >= 0 &&
             is.finite(input$zc_sd) && input$zc_sd >= 0,
           'Population standard deviations must be non-negative.'),
      need(is.finite(input$adjust_zc) && input$adjust_zc > 0,
           'The budbreak heat multiplier must be positive.'),
      need(is.finite(input$n_pop) && input$n_pop >= 1 && input$n_pop == floor(input$n_pop),
           'The number of buds must be a positive integer.'),
      need(is.null(seed) || (is.finite(seed) && seed >= 0 && seed == floor(seed)),
           'Enter a non-negative integer seed or leave it empty.'))
    config <- loaded_params()
    config$par <- par()
    config$scale_yc_budbreak <- input$adjust_zc
    config$distribution$yc_sd <- input$yc_sd
    config$distribution$zc_sd <- input$zc_sd
    config$distribution$dist_chill <- input$dist_chill
    config$distribution$dist_heat <- input$dist_heat
    skew <- c(input$skew_chill, input$skew_heat)
    config$distribution['add_par'] <- list(
      if (is.null(config$distribution$add_par) && all(skew == 0)) NULL else skew)
    config$distribution['seed'] <- list(seed)
    config$distribution$n_pop <- input$n_pop
    config
  })

  obs_list <- reactive({
    req(input$forcing_exp, input$exp_start_yday, input$exp_end_yday)
    validate(need(all(c(input$exp_start_yday, input$exp_end_yday) %in% 1:365),
                  'Cutting days must be integers between 1 and 365.'))
    observations <- model_result(app_prepare_obs_data(
      sheet = input$forcing_exp, start_yday = input$exp_start_yday,
      end_yday = input$exp_end_yday))
    validate(need(length(observations$jday_cut) > 0,
                  'No cutting experiments fall within the selected days.'))
    observations
  })

  pop_out <- reactive({
    req(input$plot_forcing)
    observations <- obs_list()
    model_result(run_population_model(
      temp_df = observations$temp_df, params = params(),
      jday_cut = observations$jday_cut, stop_at_zc = FALSE))
  })

  selected_seasons <- reactive({
    req(input$years_bloom)
    years <- app_select_years(input$years_bloom)
    validate(need(length(years) > 0, 'Select at least one available orchard year.'))
    kob_season[years]
  })

  pop_bloom <- reactive({
    req(input$plot_flowering)
    model_result(app_predict_orchard(selected_seasons(), params()))
  })

  pop_bb <- reactive({
    req(input$plot_budbreak_orchard)
    model_result(app_predict_orchard(selected_seasons(), params(), budbreak = TRUE))
  })

  output$pop_plot <- renderPlot({
    plots <- list()
    if (isTRUE(input$plot_forcing)) {
      observations <- obs_list()
      plots <- c(plots, list(helper_plot_forcing_exp_test(
        model_res_list = list(pop_out()), obs = observations$exp_obs,
        jday_cut = observations$jday_cut, jday_name = observations$jday_name,
        legend_names = 'Modelled')))
    }
    if (isTRUE(input$plot_flowering)) {
      plots <- c(plots, list(helper_plot_flowering(
        bloom_df = pop_bloom(), obs_df = kob_bloom,
        year_select = names(selected_seasons()))))
    }
    if (isTRUE(input$plot_budbreak_orchard)) {
      plots <- c(plots, list(helper_plot_budbreak_orchard(
        bloom_df = pop_bb(), obs_df = kob_bloom,
        year_select = names(selected_seasons()))))
    }
    if (isTRUE(input$plot_tempresponse)) {
      plots <- c(plots, list(model_result(helper_combined_response_plot(
        par = par(), temp_values = seq(-10, 40, by = 0.1)))))
    }
    validate(need(length(plots) > 0, 'Select at least one plot.'))
    if (length(plots) == 1L) return(plots[[1]])
    if (length(plots) == 3L) {
      return(wrap_plots(plots, design = 'AB\nCC'))
    }
    wrap_plots(plots, ncol = 2)
  })
}

shinyApp(ui = ui, server = server)

