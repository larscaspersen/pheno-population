#plot the forcing data
helper_plot_forcing_exp_test <- function(model_res_list, obs, 
                                         jday_cut, 
                                         jday_name, 
                                         forcing_days = 50, placeholder_fail = 9999,
                                         legend_names=NULL,annotate_performance = TRUE){
  
  prep_list <- helper_prepare_forcing_data(model_res_list = model_res_list, 
                                           obs = obs,
                                           jday_cut = jday_cut,
                                           jday_name = jday_name, 
                                           forcing_days = forcing_days, 
                                           placeholder_fail = placeholder_fail,
                                           legend_names = legend_names,
                                           annotate_performance = annotate_performance)
  curve_df <- prep_list$curve_df
  ssd_df <- prep_list$ssd_df
  
  
  # Assign linetypes: solid for Observed, dashed variants for model runs
  n_models   <- length(legend_names)
  lty_values <- c("solid", rep(c("dashed", "dotdash", "dotted", "longdash", "twodash"),
                               length.out = n_models))
  names(lty_values) <- c("Observed", legend_names)
  
  p_forcing <- curve_df %>%
    ggplot() +
    geom_line(aes(x = Days, y = cumsum,
                  col = source, linetype = source)) +
    scale_color_discrete(name = "Data Source") +
    scale_linetype_manual(name = "Data Source", values = lty_values) +
    ylab("Share of Buds Reaching Budbreak (%)") +
    xlab("Days After Cutting Under Forcing Conditions") +
    coord_cartesian(xlim = c(0, forcing_days + 1)) +
    scale_y_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.25)) +
    facet_wrap(~jday_fact) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  if (annotate_performance) {
    # Total SSD label: sum across all model runs, or one per run
    total_ssd_label <- ssd_df %>%
      group_by(model_run) %>%
      summarise(total = sum(ssd), .groups = "drop") %>%
      mutate(label = paste0(model_run, ": ", format(total, big.mark = ",", nsmall = 0))) %>%
      pull(label) %>%
      paste(collapse = "  |  ")
    
    p_forcing <- p_forcing +
      ggtitle(paste("Forcing Experiment. Total SSD —", total_ssd_label)) +
      ggrepel::geom_label_repel(
        data = ssd_df,
        aes(x = 1, y = 0.9,
            label = paste0(model_run, "\nSSD: ", format(ssd, nsmall = 0, big.mark = ",")),
            col = model_run),
        inherit.aes = FALSE,
        direction    = "both",
        min.segment.length = Inf,
        box.padding  = 0.8,
        point.padding = 0.8,
        segment.color = "grey40",
        force         = 5,
        show.legend   = FALSE
      )
  }
  
  return(p_forcing)
}

#function to prepare data ready for plotting
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

draw_predicted_observed_forcing <- function(po_prepared,
                                            month_display = c('Nov', 'Dec', 'Jan', 'Feb')
){
  
  mon_int <-  match(month_display, table = month.abb) 
  mon_full <- month.name[mon_int]
  
  
  p_o_forcing <- po_prepared$prepared_data |> 
    filter(source == 'Observed') %>% 
    mutate(observed = cumsum) %>% 
    select(-cumsum, -source) %>% 
    merge(rbind(po_prepared$prepared_data[po_prepared$prepared_data$source != 'Observed', c('Days', 'jday_fact', 'source', 'cumsum', 'season')],
               po_prepared$prepared_data[po_prepared$prepared_data$source != 'Observed', c('Days', 'jday_fact', 'source', 'cumsum', 'season')]),
          by = c('Days', 'jday_fact', 'season')) %>% 
    mutate(month = sub(".*\\([0-9]+ ([A-Za-z]{3}).*", "\\1", jday_fact),
           month = factor(month, levels = month_display)) %>% 
    filter(month %in% month_display) %>% 
    group_by(month, jday_fact, season) %>% 
    arrange(Days) %>% 
    mutate(lag_dif = lag(observed, default = 100) - observed) %>% 
    ungroup() %>% 
    filter(lag_dif != 0) %>% 
    mutate(season_label = factor(season,
                                 levels = c('Calibration', 'Validation'),
                                 labels = c('Calibration (2021-2022)', 'Validation (2019-2020)')),
           month_label = factor(month, levels = levels(month),
                                labels = mon_full)) %>% 
    ggplot(aes(x=observed, y = cumsum)) +
    # geom_jitter(aes(color = source),
    #            alpha = 0.4, width = 0.01, height = 0.01) +
    geom_point(aes(color = source, shape = season_label)) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    scale_y_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.2)) +
    facet_wrap(~month_label) +
    scale_shape_manual(values = c(16, 4), name = "Dataset") +
    scale_color_manual(values = c('steelblue', 'tomato'), name = "Model") +
    theme_bw(base_size = 15) +
    ylab('Predicted Share of Budbreak (%)') +
    xlab('Observed Share of Budbreak (%)') +
    theme(legend.position = "bottom") 
  
  return(p_o_forcing)
  
}
