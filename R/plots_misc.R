helper_combined_response_plot <- function(par, temp_values){
  
  response_df <- LarsChill::get_temp_response_df(par, temp_values)
  
  max_chill <- max(response_df$Chill_response)
  
  response_df %>% 
    mutate(Heat_response = Heat_response * max_chill) %>% 
    pivot_longer(cols = c('Chill_response', 'Heat_response')) %>% 
    mutate(name_plot = factor(name, 
                              levels = c('Chill_response', 'Heat_response'),
                              labels = c('Chill Response', 'Heat Response'))) %>% 
    ggplot(aes(x = Temperature, y = value, group = name_plot, col = name_plot)) +
    geom_line(aes(linetype = name_plot), show.legend = FALSE,
              size = 1.5) +
    scale_y_continuous(
      
      # Features of the first axis
      name = "Chill Response",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis( transform=~./40, name="Heat Response")
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