helper_prepare_forcing_data <- function(model_res_list, obs, 
                                        jday_cut, 
                                        jday_name, 
                                        forcing_days = 50, placeholder_fail = 9999,
                                        legend_names=NULL, annotate_performance = TRUE){
  
  # Validate inputs
  if(any(names(model_res_list) %in% c("bloomindex","x","y","z"))) model_res_list <- list(model_res_list)
  
  # Handle legend_names
  n_models <- length(model_res_list)
  if (!is.null(legend_names)) {
    if (length(legend_names) != n_models) {
      stop("'legend_names' must have the same length as 'model_res_list'")
    }
  } else {
    legend_names <- if (n_models == 1) "Modelled" else paste("Model", seq_len(n_models))
  }
  
  # Process each model run and bind together
  mod_ecdf <- purrr::map2(model_res_list, legend_names, function(model_res, leg_name) {
    
    mod_df <- model_res$exp %>%
      as.data.frame() %>%
      mutate(jday = jday_cut) %>%
      pivot_longer(cols = -jday) %>%
      mutate(jday_mod = ifelse(jday > 220, yes = jday - 365, no = jday),
             value = ifelse(value == placeholder_fail, yes = 24 * (forcing_days + 1), no = value),
             value_mod = value / 24,
             jday_fact = factor(jday, levels = jday_cut,
                                labels = jday_name))
    
    purrr::map(levels(mod_df$jday_fact), function(jdf) {
      ecdf_fun <- mod_df %>%
        filter(jday_fact == jdf) %>%
        pull(value_mod) %>%
        ecdf()
      
      data.frame(cumsum_mod = ecdf_fun(seq(0, forcing_days, by = 0.1)),
                 Days        = seq(0, forcing_days, by = 0.1),
                 jday_fact   = jdf,
                 model_run   = leg_name)
    }) %>%
      bind_rows()
    
  }) %>%
    bind_rows()
  
  ssd_df <- NULL
  # SSD per model run x chilling date combination
  if(annotate_performance){
    ssd_df <- purrr::map(legend_names, function(leg_name) {
      obs[, c('Days', 'cumsum', 'jday_fact')] %>%
        merge(mod_ecdf %>% filter(model_run == leg_name),
              by = c('Days', 'jday_fact')) %>%
        group_by(jday_fact) %>%
        summarise(ssd = sum(((cumsum - cumsum_mod) * 100)^2) %>% round(digits = 0),
                  .groups = "drop") %>%
        mutate(model_run = leg_name)
    }) %>%
      bind_rows()
  }
  
  
  # Build curve dataset: one row per Days x jday_fact x source
  curve_df <- mod_ecdf %>%
    rename(cumsum = cumsum_mod) %>%
    mutate(source = model_run) %>%
    select(-model_run) %>%
    rbind(
      cbind(obs[, c("Days", "cumsum", "jday_fact")],
            source = "Observed")
    ) %>%
    mutate(jday_fact = factor(jday_fact, levels = jday_name),
           source    = factor(source, levels = c("Observed", legend_names))) %>%
    filter(Days <= forcing_days)
  
  return(list(curve_df = curve_df, 
              ssd_df = ssd_df))
  
}


summary_experiment <- function(mod_out, 
                               time_point = c(10,20, 42),
                               proportion_budbreak = 0.5,
                               quan = c(0, 0.5, 1),
                               na_label = 9999){
  colnames_df <- c(paste0('budbreak_t', time_point),
                   paste0('days_p', proportion_budbreak))
  res_df <- data.frame()
  
  for(t in time_point){
    x <- apply(mod_out$exp,MARGIN = 1, FUN = function(x) ecdf(x)(t*24))
    res_df <- rbind(res_df, x)
  }
  
  for(p in proportion_budbreak){
    x <- apply(mod_out$exp,MARGIN = 1, FUN = quantile, probs = p)
    res_df <- rbind(res_df, ifelse(x == na_label, yes = na_label, no = round(x  / 24, digits = 1)))
  }
  res_df <- t(res_df)
  row.names(res_df) <- NULL
  colnames(res_df)  <- colnames_df
  return(res_df)
}

prepare_performance_plot_forcing <- function(obs_list,
                                             model_list_calibration, 
                                             model_list_validation,
                                             season_labels = NULL){
  if(is.null(season_labels)){
    season_labels <- c('Calibration', 'Validation')
  }
  
  obs_data <- obs_list[[1]]
  obs_data_val <- obs_list[[2]]
  
  #test to produce a predicted vs observed plot
  prep_list <- helper_prepare_forcing_data(model_res_list = model_list_calibration, 
                                           obs = obs_data$exp_obs %>% 
                                             mutate(jday_fact = factor(jday_fact,
                                                                       levels = levels(jday_fact),
                                                                       labels = obs_data$jday_name_label)), 
                                           jday_cut = obs_data$jday_cut, 
                                           jday_name = obs_data$jday_name_label,
                                           annotate_performance = FALSE,
                                           legend_names = names(model_list_calibration))
  prep_df <- prep_list$curve_df %>% 
    mutate(season = 'Calibration')
  
  #do same for validation season
  prep_list_val <- helper_prepare_forcing_data(model_res_list = model_list_validation, 
                                               obs = obs_data_val$exp_obs %>% 
                                                 mutate(jday_fact = factor(jday_fact,
                                                                           levels = levels(jday_fact),
                                                                           labels = obs_data_val$jday_name_label)), 
                                               jday_cut = obs_data_val$jday_cut, 
                                               jday_name = obs_data_val$jday_name_label,
                                               annotate_performance = FALSE,
                                               legend_names = names(model_list_validation))
  prep_df_val <- prep_list_val$curve_df %>% 
    mutate(season = 'Validation')
  
  prepared_df_combined <- prep_df |> 
    rbind(prep_df_val) |> 
    mutate(season_label = factor(season, levels = c('Calibration', 'Validation'),
                                 labels = season_labels))
  
  forc_perf <- prepared_df_combined |> 
    mutate(source = as.vector(source)) |> 
    filter(source == 'Observed') %>% 
    mutate(observed = cumsum) %>% 
    select(-cumsum, -source) %>% 
    merge(rbind(prep_df[prep_df$source != 'Observed', c('Days', 'jday_fact', 'source', 'cumsum', 'season')],
                prep_df_val[prep_df_val$source != 'Observed', c('Days', 'jday_fact', 'source', 'cumsum', 'season')]),
          by = c('Days', 'jday_fact', 'season')) %>% 
    mutate(month = sub(".*\\([0-9]+ ([A-Za-z]{3}).*", "\\1", jday_fact),
           month = factor(month, levels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb'))) %>% 
    filter(month != 'Oct') %>% 
    group_by(month, jday_fact, season) %>% 
    arrange(Days) %>% 
    mutate(lag_dif = lag(observed, default = 100) - observed) %>% 
    ungroup() %>% 
    filter(lag_dif != 0) %>% 
    mutate(season_label = factor(season,
                                 levels = c('Calibration', 'Validation'),
                                 labels = season_labels),
           month_label = factor(month, levels = levels(month),
                                labels = c(month.name[c(10:12, 1:2)]))) %>% 
    group_by(month_label, source, season) %>% 
    summarise(rmse = chillR::RMSEP(predicted = cumsum*100, observed = observed*100, na.rm = TRUE) %>% round(digits =  1),
              rpiq = chillR::RPIQ(predicted = cumsum, observed = observed, na.rm = TRUE) %>% round(digits =  1),
              mean_bias = mean(cumsum - observed, na.rm = TRUE) * 100  %>% round(digits =  1)) %>% 
    ungroup() %>% 
    pivot_longer(cols = rmse:mean_bias) %>% 
    pivot_wider(names_from = season, values_from = value) %>% 
    #mutate(label = paste0(format(round(Calibration, digits = 1), nsmall = 1), ' (', format(round(Validation,digits=1), nsmall = 1), ')')) %>% 
    mutate(label = paste0(sprintf("%.1f", Calibration), ' (', sprintf("%.1f", Validation), ')')) %>% 
    select(-Validation, -Calibration) %>% 
    pivot_wider(values_from = label) 
  
  return(list(prepared_data = prepared_df_combined |> as_tibble(),
              performance_data = forc_perf))
}
