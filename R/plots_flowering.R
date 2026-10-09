
helper_bloomint_to_jday <- function(bloomindex, x, return_fail = 9999){
  if (bloomindex == 0) {
    return(return_fail)
  }
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  if (length(unique(x$Year)) == 2 & x$Year[bloomindex] == min(x$Year)) {
    JDay <- JDay - 365
  }
  n <- length(JDaylist)
  if (n == 1) {
    return(JDay)
  }
  return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
}

helper_plot_flowering <- function(bloom_df, 
                                  obs_df, 
                                  raw_obs = TRUE,
                                  year_select = 'all', 
                                  placeholder_fail = 9999, 
                                  annotate_performance = TRUE){
  
  years <- as.numeric(year_select)
  unique_years <- unique(bloom_df$year)
  val_years <- c(2007, 2009, 2010, 2016)
  
  if(length(year_select) == 1){
    if(year_select == 'all') years <- unique_years
    if(year_select == 'calibration') years <- unique_years[(unique_years %in% val_years) == FALSE]
    if(year_select == 'validation') years <- val_years
  }
  
  strip_fills <- ifelse(years %in% val_years,
                        "#F4CCCC",   # validation (red-ish)
                        "#D9EAD3")   # calibration (green-ish)
  
  
  if(raw_obs){
    obs_df <- obs_df %>% 
      filter(year %in% years) %>% 
      select(year, firstbloom, fullbloom) %>% 
      pivot_longer(cols = -year) %>% 
      mutate(share = ifelse(name == 'firstbloom', yes = 0.1, no = 0.5))
  }
  
  
  bloom_df$value_mod <- ifelse(bloom_df$value == placeholder_fail, yes = NA, no = bloom_df$value)
  
  yday_x_scale <- (floor(min(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE))):(ceiling(max(c(bloom_df$value_mod, obs_df$value), na.rm = TRUE)))
  yday_x_scale_label <- as.Date(yday_x_scale,
                                origin = '2021-12-31') %>% 
    format("%b %d")
  
  yday_x_scale_sub <- yday_x_scale[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  yday_x_scale_label_sub <- yday_x_scale_label[seq.int(from = 1, to = length(yday_x_scale), length.out = 10)]
  
  
  
  p_flower <- bloom_df %>% 
    filter(year %in% years) %>% 
    ggplot() +
    stat_ecdf(aes(x = value,
                  color = 'Modelled'),
              geom = 'step', linewidth = 1.5) +
    ylab('Share of buds flowering (%)') +
    xlab('Date') +
    scale_x_continuous(breaks = yday_x_scale_sub, 
                       labels = yday_x_scale_label_sub,
                       limits = c(min(yday_x_scale), 
                                  max(yday_x_scale))) +
    scale_y_continuous(labels = scales::percent_format(scale = 100),
                       breaks = seq(0, 1, by = 0.25)) +
    geom_point(data = obs_df, aes(x = value, y = share, 
                                  col = 'Observed')) +
    ggh4x::facet_wrap2(~year, strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_fills,
                                                                                               colour = 'black'))) +
    scale_color_manual(values = c('steelblue', 'firebrick', '#E69F00')) +
    scale_shape_manual(values = 2) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(annotate_performance){
    
    #merge predicted and observed
    rmse_df <- bloom_df %>% 
      filter(year %in% years) %>% 
      group_by(year) %>% 
      summarise(firstbloom = quantile(value, 0.1),
                fullbloom = quantile(value, 0.5)) %>% 
      pivot_longer(cols = -year, values_to = 'pred') %>% 
      merge(obs_df, by = c('year', 'name')) %>% 
      group_by(name) %>% 
      summarise(rmse = RMSEP(predicted = pred, observed = value) %>% round(digits = 1)) %>% 
      ungroup() %>% 
      pivot_wider(values_from = rmse)
    
    #calculate difference of predicted and observed
    performance_df <- bloom_df %>% 
      filter(year %in% years) %>% 
      group_by(year) %>% 
      summarise(firstbloom = quantile(value, 0.1),
                fullbloom = quantile(value, 0.5)) %>% 
      pivot_longer(cols = -year, values_to = 'pred') %>% 
      merge(obs_df, by = c('year', 'name')) %>% 
      mutate(diff = pred - value) %>% 
      select(year, name, diff) %>% 
      pivot_wider(names_from = name, values_from = diff)
    
    p_flower <- p_flower  +
      annotate(geom = 'text', x = min(yday_x_scale), y = 1, label = 'Error',
               vjust = 1, hjust = 0) +
      geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.8, 
                                           label = paste('F1:', format(round(firstbloom, digits = 1), nsmall = 1))),
                vjust = 1, hjust = 0) +
      geom_text(data = performance_df, aes(x = min(yday_x_scale), y = 0.6, 
                                           label = paste('F2:', format(round(fullbloom, digits = 1), nsmall = 1))),
                vjust = 1, hjust = 0)+
      ggtitle(paste('Predicted Bloom KOB. RMSE Firstbloom (F1):', format(rmse_df$firstbloom, nsmall = 1),
                    'RMSE Fullbloom (F2):', format(rmse_df$fullbloom, nsmall = 1) ))
    
  }
  return(p_flower)  
}


helper_plot_budbreak_orchard <- function(bloom_df, 
                                         obs_df, 
                                         year_select = 'all', 
                                         placeholder_fail = 9999, 
                                         annotate_performance = TRUE){
  years <- as.numeric(year_select)
  unique_years <- unique(bloom_df$year)
  val_years <- c(2007, 2009, 2010, 2016)
  
  if(length(year_select) == 1){
    if(year_select == 'all') years <- unique_years
    if(year_select == 'calibration') years <- unique_years[(unique_years %in% val_years) == FALSE]
    if(year_select == 'validation') years <- val_years
  }
  
  #subset observations
  obs_df <- obs_df %>% 
    filter(year %in% years) %>% 
    select(year, budbreak) |> 
    mutate(year = as.character(year))
  
  
  bloom_df$value_mod <- ifelse(bloom_df$value == placeholder_fail, yes = NA, no = bloom_df$value)
  
  yday_y_scale <- (floor(min(c(bloom_df$value_mod, obs_df$budbreak), na.rm = TRUE))):(ceiling(max(c(bloom_df$value_mod, obs_df$budbreak), na.rm = TRUE)))
  yday_y_scale_label <- as.Date(yday_y_scale,
                                origin = '2021-12-31') %>% 
    format("%b %d")
  
  yday_y_scale_sub <- yday_y_scale[seq.int(from = 1, to = length(yday_y_scale), length.out = 10)]
  yday_y_scale_label_sub <- yday_y_scale_label[seq.int(from = 1, to = length(yday_y_scale), length.out = 10)]
  
  
  ggplot() +
    geom_point(data = obs_df, aes(x = year, y = budbreak, 
                                  col = 'Observed'))
  
  
  p_budbreak <- ggplot() +
    geom_boxplot(data = bloom_df, aes(x = year, y = value_mod,
                                      col = 'Predicted'),
                 fill = 'grey70') +
    ylab('Date of budbreak (%)') +
    xlab('Year') +
    # scale_y_continuous(breaks = yday_y_scale_sub, 
    #                    labels = yday_y_scale_label_sub,
    #                    limits = c(min(yday_y_scale), 
    #                               max(yday_y_scale))) +
    geom_point(data = obs_df, aes(x = year, y = budbreak, 
                                  col = 'Observed'),
               shape = 25, size = 2, fill = 'firebrick') +
    scale_color_manual(values = c('firebrick', 'black')) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(annotate_performance){
    
    #merge predicted and observed
    rmse_df <- bloom_df %>% 
      filter(year %in% years) %>% 
      group_by(year) |> 
      summarise(budbreak = median(value)) |> 
      pivot_longer(cols = -year, values_to = 'pred') %>% 
      merge(obs_df, by = c('year')) %>% 
      summarise(rmse = RMSEP(predicted = pred, observed = budbreak) %>% round(digits = 1)) %>% 
      ungroup() 
    
    p_budbreak <- p_budbreak  +
      ggtitle(paste('Predicted Budbreak KOB. RMSE (days):', format(rmse_df$rmse, nsmall = 1)))
    
  }
  return(p_budbreak)
}
