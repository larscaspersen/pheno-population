helper_temp_response_df <- function(par, temp_values, hours = 1200L) {
  # Constant-temperature chill over 50 days and normalized hourly heat response.
  # Use the same evalpheno components as the population model.
  chill <- vapply(temp_values, function(temperature) {
    result <- evalpheno::calculate_chill_dynamic(
      temp = rep(temperature, hours + 1L), times = 0:hours,
      E0 = par[['E0']], E1 = par[['E1']], A0 = par[['A0']], A1 = par[['A1']],
      Tf = par[['Tf']], slope = par[['slope']])
    result[nrow(result), 'y']
  }, numeric(1))
  heat <- evalpheno::calculate_heat_gdh_unscaled(
    temp = c(temp_values, tail(temp_values, 1L)),
    times = seq_len(length(temp_values) + 1L),
    Tu = par[['Tu']], Tb = par[['Tb']], Tc = par[['Tc']])
  data.frame(Temperature = temp_values, Chill_response = chill,
             Heat_response = heat)
}

helper_combined_response_plot <- function(par, temp_values){
  
  response_df <- helper_temp_response_df(par, temp_values)
  
  max_chill <- max(response_df$Chill_response)
  if (max_chill <= 0) max_chill <- 1
  
  response_df %>% 
    mutate(Heat_response = Heat_response * max_chill) %>% 
    pivot_longer(cols = c('Chill_response', 'Heat_response')) %>% 
    mutate(name_plot = factor(name, 
                              levels = c('Chill_response', 'Heat_response'),
                              labels = c('Chill Response', 'Heat Response'))) %>% 
    ggplot(aes(x = Temperature, y = value, group = name_plot, col = name_plot)) +
    geom_line(aes(linetype = name_plot), show.legend = FALSE,
              linewidth = 1.5) +
    scale_y_continuous(
      
      # Features of the first axis
      name = "Chill Portions over 50 Days",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis( transform=~./max_chill, name="Heat Response")
    ) +
    theme_bw(base_size = 15) +
    scale_color_manual(values = c('#377eb8',  '#e41a1c')) +
    scale_linetype_manual(values = c('dashed', 'solid')) +
    theme(
      # Primary Y-axis (left)
      axis.text.y.left = element_text(color = '#377eb8'),
      axis.title.y.left = element_text(color = '#377eb8'),
      axis.line.y.left = element_line(color = '#377eb8'),
      
      # Secondary Y-axis (right)
      axis.text.y.right = element_text(color = '#e41a1c'),
      axis.title.y.right = element_text(color = '#e41a1c'),
      axis.line.y.right = element_line(color = '#e41a1c'),
      legend.position = 'none'
    ) 
  
  #rescale heat
}
