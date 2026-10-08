source('R/phenoflex_transition-function.R')

draw_plot_processes <- function(obs_data, model_output, 
                                yc_pop, zc_pop, s1){
  
  #plot temperature
  #plot for temperature
  p_temp <- obs_data$temp_df %>% 
    group_by(JDay, Year) %>% 
    summarise(Tmean = mean(Temp)) %>% 
    ungroup() %>% 
    mutate(run_mean = chillR::runn_mean(vec = Tmean, runn_mean = 15),
           yday_plot = ifelse(Year == 2022, yes = JDay, no = JDay - 365)) %>% 
    ggplot(aes(x=yday_plot)) +
    geom_line(aes(y = run_mean)) +
    annotate("label",
             y = Inf,
             x = -Inf, 
             label = "A",
             hjust = -0.2, vjust = 1.2, # slight offset inside plot
             fontface = "bold") +
    xlab('Date') +
    ylab('Daily Mean\nTemperature (°C)') +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122),
                    ylim = c(0, 15)) +
    theme_bw(base_size = 13) +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  #----------------------------------#
  #format output of model
  #----------------------------------#
  
  #format output
  z_out <- model_output$z %>% 
    as.data.frame() %>% 
    mutate(yday = obs_data$temp_df$JDay,
           year = obs_data$temp_df$Year,
           yday_plot = ifelse(year == 2022, yes = yday, no = yday - 365),
           hour = rep(0:23, length(unique(obs_data$temp_df$JDay))),
           h_mod = ((hour + 1)/24)-0.5 ,
           yday_cont = yday_plot + h_mod) %>% 
    select(-yday, -year, -yday_plot, -hour, -h_mod) %>% 
    pivot_longer(cols = -yday_cont) 
  
  #summarize heat
  z_summ <- z_out %>% 
    group_by(yday_cont) %>% 
    summarise(min_z = min(value),
              med_z = median(value),
              max_z = max(value)) %>% 
    ungroup()
  
  
  #format chill
  y_out <-  data.frame(y = model_output$y) %>% 
    mutate(yday = obs_data$temp_df$JDay,
           year = obs_data$temp_df$Year,
           yday_plot = ifelse(year == 2022, yes = yday, no = yday - 365),
           hour = rep(0:23, length(unique(obs_data$temp_df$JDay))),
           h_mod = ((hour + 1)/24)-0.5 ,
           yday_cont = yday_plot + h_mod) %>% 
    select(-yday, -year, -yday_plot, -hour, -h_mod) 
  
  #calculate transition function
  
  s1_pop <- purrr::map(yc_pop, .f = function(yc) get_transition_fun(yc, s1, y_vector = y_out$y)) %>% 
    bind_cols() %>% 
    dplyr::rename_with(~ str_replace_all(.x, "\\.\\.\\.", "V")) %>% 
    mutate(yday_cont = y_out$yday_cont) %>% 
    pivot_longer(cols = -yday_cont)
  
  s1_summ <- s1_pop %>% 
    group_by(yday_cont) %>% 
    summarise(min_py = min(value),
              med_py = median(value),
              max_py = max(value)) %>% 
    ungroup()
  
  #-----------------#
  #calculate spread of chill and heat
  #-----------------#
  yc_spread <- c(min(yc_pop), max(yc_pop))
  zc_spread <- c(min(zc_pop), max(zc_pop))
  index_spread <- c(min(model_output$bloomindex), max(model_output$bloomindex))
  yday_spread <- c(y_out$yday_cont[index_spread[1]], y_out$yday_cont[index_spread[2]])
  
  #spread when yc is met
  index_yc_spread <-  c(which(y_out$y >= yc_spread[1]) %>% min(),
                        which(y_out$y >= yc_spread[2]) %>% min())
  yday_yc_spread <- c(y_out$yday_cont[index_yc_spread[1]], y_out$yday_cont[index_yc_spread[2]])
  
  
  #----------------#
  #draw chill
  #----------------#
  
  p_chill <- y_out %>% 
    ggplot(aes(x = yday_cont)) +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = yc_spread[1], ymax = yc_spread[2],
             alpha = 0.2) +
    annotate("segment",
             x = -97, xend = -97,
             y = yc_spread[1], yend = yc_spread[2],
             arrow = arrow(ends = "both",
                           type = "closed",
                           length = unit(0.2, "cm"))) +
    annotate("text",
             x = -70,
             y = mean(c(yc_spread[1], yend = yc_spread[2])), 
             label = "atop(Chill~requirement,(y[c]))",,parse = TRUE,
             angle = 0) +
    annotate("rect",
             xmin = yday_yc_spread[1], xmax = yday_yc_spread[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.2) +
    annotate("text",
             y = 20,
             x = mean(c(yday_yc_spread[1], yend = yday_yc_spread[2])), 
             label = "atop(Chill~requirement~met, (y == y[c]))",,parse = TRUE) +
    annotate("segment",
             y = 0, yend = 0,
             x = yday_yc_spread[1], xend = yday_yc_spread[2],
             arrow = arrow(ends = "both",
                           type = "closed",
                           length = unit(0.2, "cm"))) +
    geom_line(aes(y = y)) +
    annotate("label",
             y = Inf,
             x = -Inf, 
             label = "B",
             hjust = -0.2, vjust = 1.2, # slight offset inside plot
             fontface = "bold") +
    xlab('Date') +
    ylab('Accumulated Chill\nPortions (y)') +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122)) +
    theme_bw(base_size = 13) +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  #--------------#
  #draw heat
  #--------------#
  p_heat <- z_summ %>% 
    ggplot(aes(x = yday_cont)) +
    #annotate heat requirement
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = zc_spread[1], ymax = zc_spread[2],
             alpha = 0.2) +
    # annotate("segment",
    #          x = -80, xend = -80,
    #          y = zc_spread[1], yend = zc_spread[2],
    #          arrow = arrow(ends = "both",
    #                        type = "closed",
    #                        length = unit(0.2, "cm"))) +
    annotate("text",
             x = -50,
             y = mean(c(zc_spread[1], yend = zc_spread[2])), 
             label = "Heat~requirement~(z[c])",,parse = TRUE) +
    #annotate bloom date
    annotate("rect",
             xmin = yday_spread[1], xmax = yday_spread[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.2) +
    # annotate("segment",
    #          y = 20, yend = 20,
    #          x = yday_spread[1], xend = yday_spread[2],
    #          arrow = arrow(ends = "both",
    #                        type = "closed",
    #                        length = unit(0.2, "cm"))) +
    annotate("text",
             y = 50,
             x = mean(c(yday_spread[1], yday_spread[2])), 
             label = "atop(Bloom, (z == z[c]))",,parse = TRUE) +
    geom_ribbon(aes(ymin = min_z, ymax = max_z),
                fill = '#FF474C', alpha = 0.7) +
    geom_line(aes(y = med_z)) +
    annotate("label",
             y = Inf,
             x = -Inf, 
             label = "D",
             hjust = -0.2, vjust = 1.2, # slight offset inside plot
             fontface = "bold") +
    xlab('Date') +
    ylab('Accumulated\nHeat (z)') +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122), ylim = c(0, 500)) +
    theme_bw(base_size = 13)
  
  #------------------#
  #draw transition function
  #------------------#
  p_py <- s1_summ %>% 
    ggplot(aes(x = yday_cont)) +
    geom_ribbon(aes(ymin = min_py, ymax = max_py),
                fill = 'lightblue', alpha = 0.7) +
    geom_line(aes(y = med_py)) +
    annotate("label",
             y = Inf,
             x = -Inf, 
             label = "C",
             hjust = -0.2, vjust = 1.2, # slight offset inside plot
             fontface = "bold") +
    xlab('Date') +
    ylab(expression(
      "Competence to Grow ("*P[y]*")")) +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122)) +
    theme_bw(base_size = 13) +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  
  #combine plots
  p_temp / p_chill / p_py / p_heat  
  
  
  
}

draw_plot_processes_v2 <- function(
    obs_data,
    model_output,
    yc_pop,
    zc_pop,
    s1,
    
    # Optional stage-specific model
    stage_model_output = NULL,
    yc_stage = NULL,
    zc_stage = NULL,
    s1_stage = NULL,
    
    # Labels
    axis_labels = list(),
    model_labels = c(
      population = "PhenoFlex (population)",
      stage_specific = "PhenoFlex (stage-specific)"
    ),
    legend_title = "Model",
    
    # Colors
    model_colors = c(
      population = "tomato",
      stage_specific = "steelblue"
    ),
    
    # Axis settings
    date_breaks = c(
      275 - 365, 306 - 365, 336 - 365,
      1, 32, 61, 92, 122
    ),
    date_labels = c(
      "Oct", "Nov", "Dec", "Jan",
      "Feb", "Mar", "Apr", "May"
    ),
    x_limits = c(275 - 365, 122),
    temp_limits = c(0, 15),
    heat_limits = c(0, 500)
) {
  
  # ------------------------------------------------------------------
  # Defaults and input checks
  # ------------------------------------------------------------------
  
  default_axis_labels <- list(
    x = "Date",
    temperature = "Daily Mean\nTemperature (°C)",
    chill = "Accumulated Chill\nPortions (y)",
    competence = expression(
      "Competence to Grow (" * P[y] * ")"
    ),
    heat = "Accumulated\nHeat (z)"
  )
  
  axis_labels <- utils::modifyList(
    default_axis_labels,
    axis_labels
  )
  
  required_model_names <- c("population", "stage_specific")
  
  if (!all(required_model_names %in% names(model_labels))) {
    stop(
      "`model_labels` must contain named elements ",
      "`population` and `stage_specific`."
    )
  }
  
  if (!all(required_model_names %in% names(model_colors))) {
    stop(
      "`model_colors` must contain named elements ",
      "`population` and `stage_specific`."
    )
  }
  
  has_stage_model <- !is.null(stage_model_output)
  
  if (
    has_stage_model &&
    any(vapply(
      list(yc_stage, zc_stage, s1_stage),
      is.null,
      logical(1)
    ))
  ) {
    stop(
      "When `stage_model_output` is supplied, you must also provide ",
      "`yc_stage`, `zc_stage`, and `s1_stage`."
    )
  }
  
  pop_name <- unname(model_labels["population"])
  stage_name <- unname(model_labels["stage_specific"])
  
  color_values <- stats::setNames(
    unname(model_colors[required_model_names]),
    unname(model_labels[required_model_names])
  )
  
  # ------------------------------------------------------------------
  # Create continuous time axis
  # ------------------------------------------------------------------
  
  time_index <- obs_data$temp_df %>%
    dplyr::transmute(
      yday = JDay,
      year = Year
    ) %>%
    dplyr::group_by(year, yday) %>%
    dplyr::mutate(
      hour = dplyr::row_number() - 1,
      n_hours = dplyr::n(),
      yday_plot = dplyr::if_else(
        year == 2022,
        as.numeric(yday),
        as.numeric(yday) - 365
      ),
      h_mod = ((hour + 1) / n_hours) - 0.5,
      yday_cont = yday_plot + h_mod
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(yday_cont)
  
  n_time <- nrow(time_index)
  
  if (length(model_output$y) != n_time) {
    stop(
      "The length of `model_output$y` does not match ",
      "the number of temperature records."
    )
  }
  
  # Safely find the first threshold crossing
  first_crossing <- function(x, threshold) {
    index <- which(x >= threshold)
    
    if (length(index) == 0) {
      return(NA_integer_)
    }
    
    index[1]
  }
  
  # ------------------------------------------------------------------
  # Temperature
  # ------------------------------------------------------------------
  
  p_temp <- obs_data$temp_df %>%
    dplyr::group_by(JDay, Year) %>%
    dplyr::summarise(
      Tmean = mean(Temp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Year, JDay) %>%
    dplyr::mutate(
      run_mean = chillR::runn_mean(
        vec = Tmean,
        runn_mean = 15
      ),
      yday_plot = dplyr::if_else(
        Year == 2022,
        as.numeric(JDay),
        as.numeric(JDay) - 365
      )
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = yday_plot, y = run_mean)
    ) +
    ggplot2::geom_line() +
    ggplot2::annotate(
      "label",
      y = Inf,
      x = -Inf,
      label = "A",
      hjust = -0.2,
      vjust = 1.2,
      fontface = "bold"
    ) +
    ggplot2::labs(
      x = axis_labels$x,
      y = axis_labels$temperature
    ) +
    ggplot2::scale_x_continuous(
      breaks = date_breaks,
      labels = date_labels
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = temp_limits
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------------
  # Population-model output
  # ------------------------------------------------------------------
  
  y_out <- dplyr::bind_cols(
    time_index,
    tibble::tibble(y = as.numeric(model_output$y))
  )
  
  z_matrix <- as.matrix(model_output$z)
  
  if (nrow(z_matrix) != n_time) {
    stop(
      "The number of rows in `model_output$z` does not match ",
      "the number of temperature records."
    )
  }
  
  z_summ <- tibble::tibble(
    yday_cont = time_index$yday_cont,
    min_z = apply(z_matrix, 1, min, na.rm = TRUE),
    med_z = apply(z_matrix, 1, median, na.rm = TRUE),
    max_z = apply(z_matrix, 1, max, na.rm = TRUE),
    model = pop_name
  )
  
  py_matrix <- do.call(
    cbind,
    lapply(
      yc_pop,
      function(yc) {
        as.numeric(
          get_transition_fun(
            yc = yc,
            s1 = s1,
            y_vector = y_out$y
          )
        )
      }
    )
  )
  
  s1_summ <- tibble::tibble(
    yday_cont = time_index$yday_cont,
    min_py = apply(py_matrix, 1, min, na.rm = TRUE),
    med_py = apply(py_matrix, 1, median, na.rm = TRUE),
    max_py = apply(py_matrix, 1, max, na.rm = TRUE),
    model = pop_name
  )
  
  # ------------------------------------------------------------------
  # Population requirement and event ranges
  # ------------------------------------------------------------------
  
  yc_spread <- range(yc_pop, na.rm = TRUE)
  zc_spread <- range(zc_pop, na.rm = TRUE)
  
  index_yc_spread <- c(
    first_crossing(y_out$y, yc_spread[1]),
    first_crossing(y_out$y, yc_spread[2])
  )
  
  if (all(is.finite(index_yc_spread))) {
    yday_yc_spread <- time_index$yday_cont[index_yc_spread]
  } else {
    yday_yc_spread <- c(NA_real_, NA_real_)
  }
  
  bloom_indices <- as.integer(model_output$bloomindex)
  
  bloom_indices <- bloom_indices[
    is.finite(bloom_indices) &
      bloom_indices >= 1 &
      bloom_indices <= n_time
  ]
  
  if (length(bloom_indices) > 0) {
    yday_spread <- range(
      time_index$yday_cont[bloom_indices],
      na.rm = TRUE
    )
  } else {
    yday_spread <- c(NA_real_, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # Stage-specific model output
  # ------------------------------------------------------------------
  
  if (has_stage_model) {
    
    if (
      length(stage_model_output$y) != n_time ||
      length(stage_model_output$z) != n_time
    ) {
      stop(
        "The stage-specific model output does not match ",
        "the number of temperature records."
      )
    }
    
    # Rescale stage-specific chill so that its requirement is shown
    # at the center of the population chill-requirement interval.
    yc_target <- mean(yc_spread)
    chill_rescale_factor <- yc_target / yc_stage
    
    stage_out <- tibble::tibble(
      yday_cont = time_index$yday_cont,
      y = as.numeric(stage_model_output$y),
      z = as.numeric(stage_model_output$z),
      model = stage_name
    ) %>%
      dplyr::mutate(
        y_rescaled = y * chill_rescale_factor,
        py = get_transition_fun(
          yc = yc_stage,
          s1 = s1_stage,
          y_vector = y
        )
      )
    
    stage_yc_index <- first_crossing(
      stage_out$y,
      yc_stage
    )
    
    stage_zc_index <- first_crossing(
      stage_out$z,
      zc_stage
    )
    
    stage_yc_segment <- tibble::tibble(
      model = stage_name,
      x = ifelse(
        is.na(stage_yc_index),
        NA_real_,
        stage_out$yday_cont[stage_yc_index]
      ),
      y = yc_target
    )
    
    stage_zc_segment <- tibble::tibble(
      model = stage_name,
      x = ifelse(
        is.na(stage_zc_index),
        NA_real_,
        stage_out$yday_cont[stage_zc_index]
      ),
      y = zc_stage
    )
  }
  
  # ------------------------------------------------------------------
  # Chill panel
  # ------------------------------------------------------------------
  
  p_chill <- y_out %>%
    dplyr::mutate(model = pop_name) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = yday_cont)
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = yc_spread[1],
      ymax = yc_spread[2],
      alpha = 0.2
    ) +
    ggplot2::annotate(
      "text",
      x = -70,
      y = mean(yc_spread),
      label = "atop(Chill~requirement, (y[c]))",
      parse = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = y, colour = model),
      show.legend = has_stage_model
    )
  
  if (all(is.finite(yday_yc_spread))) {
    p_chill <- p_chill +
      ggplot2::annotate(
        "rect",
        xmin = yday_yc_spread[1],
        xmax = yday_yc_spread[2],
        ymin = -Inf,
        ymax = yc_spread[1],
        alpha = 0.2
      ) +
      ggplot2::annotate(
        "rect",
        xmin = yday_yc_spread[1],
        xmax = yday_yc_spread[2],
        ymin = yc_spread[2],
        ymax = Inf,
        alpha = 0.2
      ) +
      ggplot2::annotate(
        "text",
        x = mean(yday_yc_spread),
        y = 20,
        label = "atop(Chill~requirement, met~(y == y[c]))",
        parse = TRUE
      )
  }
  
  if (has_stage_model) {
    p_chill <- p_chill +
      ggplot2::geom_line(
        data = stage_out,
        ggplot2::aes(
          y = y_rescaled,
          colour = model
        )
      )
    
    if (is.finite(stage_yc_segment$x)) {
      p_chill <- p_chill +
        ggplot2::geom_segment(
          data = stage_yc_segment,
          ggplot2::aes(
            x = -Inf,
            xend = x,
            y = y,
            yend = y,
            colour = model
          ),
          inherit.aes = FALSE,
          linetype = "dashed"
        ) +
        ggplot2::geom_segment(
          data = stage_yc_segment,
          ggplot2::aes(
            x = x,
            xend = x,
            y = y,
            yend = -Inf,
            colour = model
          ),
          inherit.aes = FALSE,
          linetype = "dashed"
        )
    }
  }
  
  p_chill <- p_chill +
    ggplot2::annotate(
      "label",
      y = Inf,
      x = -Inf,
      label = "B",
      hjust = -0.2,
      vjust = 1.2,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(
      values = color_values,
      breaks = unname(model_labels),
      name = legend_title
    ) +
    ggplot2::labs(
      x = axis_labels$x,
      y = axis_labels$chill
    ) +
    ggplot2::scale_x_continuous(
      breaks = date_breaks,
      labels = date_labels
    ) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------------
  # Competence-to-grow panel
  # ------------------------------------------------------------------
  
  p_py <- s1_summ %>%
    ggplot2::ggplot(
      ggplot2::aes(x = yday_cont)
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = min_py,
        ymax = max_py
      ),
      fill = unname(model_colors["population"]),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = med_py,
        colour = model
      ),
      show.legend = FALSE
    )
  
  if (has_stage_model) {
    p_py <- p_py +
      ggplot2::geom_line(
        data = stage_out,
        ggplot2::aes(
          y = py,
          colour = model
        ),
        show.legend = FALSE
      )
  }
  
  p_py <- p_py +
    ggplot2::annotate(
      "label",
      y = Inf,
      x = -Inf,
      label = "C",
      hjust = -0.2,
      vjust = 1.2,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(
      values = color_values,
      breaks = unname(model_labels),
      name = legend_title
    ) +
    ggplot2::labs(
      x = axis_labels$x,
      y = axis_labels$competence
    ) +
    ggplot2::scale_x_continuous(
      breaks = date_breaks,
      labels = date_labels
    ) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  # ------------------------------------------------------------------
  # Heat panel
  # ------------------------------------------------------------------
  
  p_heat <- z_summ %>%
    ggplot2::ggplot(
      ggplot2::aes(x = yday_cont)
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = zc_spread[1],
      ymax = zc_spread[2],
      alpha = 0.2
    ) +
    ggplot2::annotate(
      "text",
      x = -50,
      y = mean(zc_spread),
      label = "Heat~requirement~(z[c])",
      parse = TRUE
    )
  
  if (all(is.finite(yday_spread))) {
    p_heat <- p_heat +
      ggplot2::annotate(
        "rect",
        xmin = yday_spread[1],
        xmax = yday_spread[2],
        ymin = -Inf,
        ymax = zc_spread[1],
        alpha = 0.2
      ) +
      ggplot2::annotate(
        "rect",
        xmin = yday_spread[1],
        xmax = yday_spread[2],
        ymin = zc_spread[2],
        ymax = Inf,
        alpha = 0.2
      ) +
      ggplot2::annotate(
        "text",
        x = mean(yday_spread),
        y = 50,
        label = "atop(Bloom, (z == z[c]))",
        parse = TRUE
      )
  }
  
  p_heat <- p_heat +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = min_z,
        ymax = max_z
      ),
      fill = unname(model_colors["population"]),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = med_z,
        colour = model
      ),
      show.legend = FALSE
    )
  
  if (has_stage_model) {
    p_heat <- p_heat +
      ggplot2::geom_line(
        data = stage_out,
        ggplot2::aes(
          y = z,
          colour = model
        ),
        show.legend = FALSE
      )
    
    if (is.finite(stage_zc_segment$x)) {
      p_heat <- p_heat +
        ggplot2::geom_segment(
          data = stage_zc_segment,
          ggplot2::aes(
            x = -Inf,
            xend = x,
            y = y,
            yend = y,
            colour = model
          ),
          inherit.aes = FALSE,
          linetype = "dashed",
          show.legend = FALSE
        ) +
        ggplot2::geom_segment(
          data = stage_zc_segment,
          ggplot2::aes(
            x = x,
            xend = x,
            y = y,
            yend = -Inf,
            colour = model
          ),
          inherit.aes = FALSE,
          linetype = "dashed",
          show.legend = FALSE
        )
    }
  }
  
  p_heat <- p_heat +
    ggplot2::annotate(
      "label",
      y = Inf,
      x = -Inf,
      label = "D",
      hjust = -0.2,
      vjust = 1.2,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(
      values = color_values,
      breaks = unname(model_labels),
      name = legend_title
    ) +
    ggplot2::labs(
      x = axis_labels$x,
      y = axis_labels$heat
    ) +
    ggplot2::scale_x_continuous(
      breaks = date_breaks,
      labels = date_labels
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = heat_limits
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "none")
  
  # ------------------------------------------------------------------
  # Combine panels
  # ------------------------------------------------------------------
  
  patchwork::wrap_plots(
    p_temp,
    p_chill,
    p_py,
    p_heat,
    ncol = 1,
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = if (has_stage_model) "bottom" else "none"
    )
}



summarize_processes <- function(obs_data, model_out, par){
  
  #format chill
  y_out <-  data.frame(y = model_out$y) %>% 
    mutate(yday = obs_data$temp_df$JDay,
           year = obs_data$temp_df$Year,
           yday_plot = ifelse(year == 2022, yes = yday, no = yday - 365),
           hour = rep(0:23, length(unique(obs_data$temp_df$JDay))),
           h_mod = ((hour + 1)/24)-0.5 ,
           yday_cont = yday_plot + h_mod) %>% 
    select(-yday, -year, -yday_plot, -hour, -h_mod) 
  
  #numbers for manuscript
  #mean + range when yc met
  #mean + range when Py is more than 90%
  #mean + range for flowering
  index_spread <- c(min(model_out$bloomindex), quantile(model_out$bloomindex, 0.25),
                    median(model_out$bloomindex), mean(model_out$bloomindex), quantile(model_out$bloomindex, 0.75),
                    max(model_out$bloomindex))
  yday_spread <- c(y_out$yday_cont[index_spread[1]], y_out$yday_cont[index_spread[2]], y_out$yday_cont[index_spread[3]],
                   y_out$yday_cont[index_spread[4]], y_out$yday_cont[index_spread[5]],  y_out$yday_cont[index_spread[6]])
  
  yday_sum <- as.Date(yday_spread, origin = '2021-12-31')
  names(yday_sum) <- c('min', '25%', 'median', 'mean', '75%', 'max')
  yday_span <- yday_spread[3] - yday_spread[1]
  
  #spread when yc is met
  yc_spread <- c(min(model_out$yc_pop), quantile(model_out$yc_pop, 0.25),
                 median(model_out$yc_pop), mean(model_out$yc_pop), quantile(model_out$yc_pop, 0.75),
                 max(model_out$yc_pop))
  
  index_yc_spread <-  c(which(y_out$y >= yc_spread[1]) %>% min(),
                        which(y_out$y >= yc_spread[2]) %>% min(),
                        which(y_out$y >= yc_spread[3]) %>% min(),
                        which(y_out$y >= yc_spread[4]) %>% min(),
                        which(y_out$y >= yc_spread[5]) %>% min(),
                        which(y_out$y >= yc_spread[6]) %>% min())
  yday_yc_spread <- c(y_out$yday_cont[index_yc_spread[1]], y_out$yday_cont[index_yc_spread[2]],
                      y_out$yday_cont[index_yc_spread[3]], y_out$yday_cont[index_yc_spread[4]],
                      y_out$yday_cont[index_yc_spread[5]], y_out$yday_cont[index_yc_spread[6]])
  yday_yc_sum <- as.Date(yday_yc_spread, origin = '2021-12-31')
  names(yday_yc_sum) <- c('min', '25%', 'median', 'mean', '75%', 'max')
  
  
  #get populaiton of transition functions
  s1_pop <- purrr::map(model_out$yc_pop, .f = function(yc) get_transition_fun(yc = yc, s1 = par[3], y_vector = y_out$y)) %>% 
    bind_cols() %>% 
    dplyr::rename_with(~ str_replace_all(.x, "\\.\\.\\.", "V")) %>% 
    mutate(yday_cont = y_out$yday_cont) %>% 
    pivot_longer(cols = -yday_cont)
  
  s1_sum <- s1_pop  %>% 
    group_by(name) %>% 
    summarise(index_09 = which(value >= 0.9) %>% min(),
              yday_09 =  yday_cont[index_09]) %>% 
    pull(yday_09) %>% 
    summary() %>% 
    as.numeric() %>% 
    as.Date(origin = '2021-12-31')
  names(s1_sum) <- c('min', '25%', 'median', 'mean', '75%', 'max')
  
  
  return(list(yc = yday_yc_sum,
              py_90 = s1_sum,
              bloom = yday_sum))
}
