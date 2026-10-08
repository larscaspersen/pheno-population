source('experimental/R/phenoflex_transition-function.R')

get_stage_bud <- function(y, z, yc, zc, s1, 
                          py_endo = 0.1,
                          py_eco = 0.9,
                          p_zc_eco = 0.2){
  
  #transition function
  py <- get_transition_fun(yc = yc,s1 = s1, y_vector = y)
  
  #endormancy
  start_endo <- 1
  end_endo <- which(py > py_endo) %>% min()
  
  #flowering
  start_flower <- which(z >= zc) %>% min()
  end_flower <- length(y)
  
  #transition from endo to ecodormancy
  start_transition <- end_endo + 1
  end_transition <- which(py >= py_endo & py < py_eco) %>% max()
  
  #ecodormancy
  start_eco <- end_transition + 1
  end_eco <- which(z > (p_zc_eco * zc)) %>% min()
  
  #ontogenetic growth
  start_onto <- end_eco +1
  end_onto <- start_flower -1
  
  stage_vec <- rep(NA, length(y))
  stage_vec[start_endo:end_endo] <- 0
  stage_vec[start_transition:end_transition] <- 1
  stage_vec[start_eco:end_eco] <- 2
  stage_vec[start_onto:end_onto] <- 3
  stage_vec[start_flower:end_flower] <- 4
  
  return(stage_vec)
  
}

get_stage_pop <- function(mod_out, yc_pop, zc_pop, s1,
                          py_endo = 0.1,
                          py_eco = 0.9,
                          p_zc_eco = 0.2){
  
  
  stages_list <- purrr::map(1:ncol(mod_out$z), function(i){
    get_stage_bud(y = mod_out$y, 
                  z = mod_out$z[,i], 
                  yc = yc_pop[i], 
                  zc = zc_pop[i], 
                  s1 = s1, 
                  py_endo = py_endo,
                  py_eco = py_eco,
                  p_zc_eco = p_zc_eco)
  })
  
  do.call(cbind, stages_list) %>% 
    return()
  
  
}

add_zero_anchors <- function(df) {
  # df has columns: jday_h_plot, value (label), share, name
  time_step <- min(diff(sort(unique(df$jday_h_plot))))
  
  df %>%
    group_by(name, value) %>%
    arrange(jday_h_plot) %>%
    reframe(
      jday_h_plot = c(min(jday_h_plot) - time_step,   # anchor before first rise
                      jday_h_plot,
                      max(jday_h_plot) + time_step),   # anchor after last drop
      share       = c(0, share, 0),
      .groups = "drop"
    )
}


plot_stages_bars <- function(stages,
                             obs_data,
                             stage_names = c('Endodormancy', 
                                             'Transition:\nEndo- to Ecodormancy',
                                             'Ecodormancy', 
                                             'Growth', 
                                             'Flowering'),
                             fill_scheme = c('#0571b0','grey85', '#d7191c', '#a6d96a', '#7b3294'),
                             next_stage_on_bottom = FALSE){
  
  # Build a lookup: code (0-based) -> display label
  n_codes      <- length(stage_names)
  code_to_label <- setNames(stage_names, as.character(seq_len(n_codes) - 1))
  
  # Unique labels in code order (for factor levels / fill scale)
  unique_labels <- unique(stage_names)
  
  # One color per unique label (take first occurrence color for each)
  unique_colors <- fill_scheme[match(unique_labels, stage_names)]
  names(unique_colors) <- unique_labels
  
  n_run <- ncol(stages) 
  
  df_plot <- stages %>% 
    as.data.frame() %>% 
    mutate(jday = obs_data$temp_df$JDay,
           year = obs_data$temp_df$Year,
           jday_plot = ifelse(year == 2022, yes = jday, no = jday - 365),
           h = rep(0:23, nrow(obs_data$temp_df) / 24),
           jday_h_plot = jday_plot + (h /24) - 0.5) %>% 
    select(-jday, -year, -h) %>% 
    pivot_longer(cols = starts_with("V")) %>% 
    mutate(value = code_to_label[as.character(value)]) %>%
    group_by(jday_h_plot, value) %>% 
    # Dynamically count share per unique label
    summarise(share = n() / n_run,
              .groups = "drop") %>%
    mutate(name = factor(value, levels = rev(unique_labels)))
  
  if(next_stage_on_bottom){
    df_plot <- df_plot %>% 
      mutate(name = factor(value, levels = (unique_labels)))
  }
                                                      
  p_stages <-  df_plot %>%
    # summarise(endo = sum(value == 0) / n(),
    #           trans = sum(value == 1) / n(),
    #           eco = sum(value == 2) / n(),
    #           onto = sum(value == 3) / n(),
    #           flower = sum(value == 4) / n()) %>% 
    # pivot_longer(cols = endo:flower) %>% 
    # mutate(name = factor(name, levels = rev(c('endo', 'trans', 'eco', 'onto', 'flower')),
    #                      labels = rev(c('Endodormancy', 'Transition: Endo- to Ecodormancy',
    #                                     'Ecodormancy',
    #                                     'Ontogenetic Growth',
    #                                     'Flowering')))) %>% 
    ggplot(aes(x = jday_h_plot, y = share, fill = name)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    ylab('Comulative Density of Population') +
    xlab('Date') +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122)) +
    theme_bw(base_size = 12) +
    scale_fill_manual(name = 'Phenological\nStage',
                      values = unique_colors,
                      breaks = unique_labels) +
    scale_y_continuous(expand = c(0,0))+
    guides(fill = guide_legend(nrow = 2)) + 
    theme(legend.position = 'bottom')
  
  return(p_stages)
}




plot_stages_lines <- function(stages,
                              obs_data,
                             stage_names = c('Endodormancy', 
                                             'Transition:\nEndo- to Ecodormancy',
                                             'Ecodormancy', 
                                             'Growth', 
                                             'Flowering'),
                             fill_scheme = c('#0571b0','grey65', '#d7191c', '#a6d96a', '#7b3294' )){
  # Build a lookup: code (0-based) -> display label
  n_codes      <- length(stage_names)
  code_to_label <- setNames(stage_names, as.character(seq_len(n_codes) - 1))
  
  # Unique labels in code order (for factor levels / fill scale)
  unique_labels <- unique(stage_names)
  
  # One color per unique label (take first occurrence color for each)
  unique_colors <- fill_scheme[match(unique_labels, stage_names)]
  names(unique_colors) <- unique_labels
  
  n_run <- ncol(stages)                                                    
  
  p_stages <-  stages %>% 
    as.data.frame() %>% 
    mutate(jday = obs_data$temp_df$JDay,
           year = obs_data$temp_df$Year,
           jday_plot = ifelse(year == 2022, yes = jday, no = jday - 365),
           h = rep(0:23, nrow(obs_data$temp_df) / 24),
           jday_h_plot = jday_plot + (h /24) - 0.5) %>% 
    select(-jday, -year, -h) %>% 
    pivot_longer(cols = starts_with("V")) %>% 
    mutate(value = code_to_label[as.character(value)]) %>%
    group_by(jday_h_plot, value) %>% 
    # Dynamically count share per unique label
    summarise(share = n() / n_run,
              .groups = "drop") %>%
    mutate(name = factor(value, levels = rev(unique_labels))) %>%
    add_zero_anchors() %>% 
    # summarise(endo = sum(value == 0) / n(),
    #           trans = sum(value == 1) / n(),
    #           eco = sum(value == 2) / n(),
    #           onto = sum(value == 3) / n(),
    #           flower = sum(value == 4) / n()) %>% 
    # pivot_longer(cols = endo:flower) %>% 
    # mutate(name = factor(name, levels = rev(c('endo', 'trans', 'eco', 'onto', 'flower')),
    #                      labels = rev(c('Endodormancy', 'Transition: Endo- to Ecodormancy',
    #                                     'Ecodormancy',
    #                                     'Ontogenetic Growth',
    #                                     'Flowering')))) %>% 
    ggplot(aes(x = jday_h_plot, y = share, color = name)) +
    geom_line() +
    ylab('Comulative Density of Population') +
    xlab('Date') +
    scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                       labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
    coord_cartesian(xlim = c(275-365, 122)) +
    theme_bw(base_size = 12) +
    scale_color_manual(name = 'Phenological\nStage',
                      values = unique_colors,
                      breaks = unique_labels) +
    scale_y_continuous(expand = c(0,0.01))+
    guides(color = guide_legend(nrow = 2)) + 
    theme(legend.position = 'bottom')
  
  return(p_stages)
}

