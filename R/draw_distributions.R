draw_dist_plots <- function(n, par, yc_sd, zc_sd, dist_chill, dist_heat, add_par=NULL,
                            seed = 12345){
  
  #sample distribudions (from the helpers_phenoflex_pop script)
  dist_out <- helper_sample_dist(n = n, par = par, yc_sd = yc_sd, zc_sd = zc_sd,
                                 dist_chill = dist_chill, dist_heat = dist_heat, add_par = add_par,
                                 seed = seed)
  
  yc_pop <- dist_out$yc_pop 
  zc_pop <- dist_out$zc_pop 
  s1 <- par[3]
  
  
  #labels for plotting
  mu <-  "\u03BC"
  sigma <- "\u03C3"
  val_mu <- format(round(par[1], digits = 0), nsmall = 0)
  val_sigma <- format(round(yc_sd, digits = 0), nsmall = 0)
  
  #histogram of yc
  p_ycpop <- yc_pop %>% 
    data.frame() %>% 
    ggplot(aes(x = .)) +
    geom_histogram() +
    annotate(geom = 'label', x = 40, y = 12, label = paste0('N(', mu,'=', val_mu, ' ,', sigma, '=', val_sigma, ')'),
             hjust = 0, vjust = 1) +
    ylab('Count') +
    xlab(expression(Chill~Requirement~(y[c]))) +
    theme_bw(base_size = 15)
  
  #histogram for zc
  val_mu <- format(round(par[2], digits = 0), nsmall = 0)
  val_sigma <- format(round(zc_sd, digits = 0), nsmall = 0)
  
  p_zcpop <- zc_pop %>% 
    data.frame() %>% 
    ggplot(aes(x = .)) +
    geom_histogram() +
    annotate(geom = 'label', x = 227, y = 12, 
             label =  paste0('N(', mu,'=', val_mu, ' ,', sigma, '=', val_sigma, ')'),
             hjust = 0, vjust = 1) +
    ylab('Count') +
    xlab(expression(Heat~Requirement~(z[c]))) +
    theme_bw(base_size = 15)
  
  #transition function
  
  #span the transition plot between 50% yc and 150% yc
  y_vector <- seq(min(yc_pop) * 0.5, max(yc_pop) * 1.5, 0.5)
  
  #draw the curves for min, max and median yc
  s1_min <- get_transition_fun(yc = min(yc_pop), s1, y_vector = y_vector)
  s1_max <- get_transition_fun(yc = max(yc_pop), s1, y_vector = y_vector)
  s1_med <- get_transition_fun(yc = median(yc_pop), s1, y_vector = y_vector)
  
  p_transition <- data.frame(y_vector = y_vector,
                             min = s1_min,
                             med = s1_med,
                             max = s1_max) %>% 
    ggplot(aes(x = y_vector)) +
    geom_ribbon(aes(ymin = min, ymax = max),
                fill = 'lightblue', alpha = 0.7) +
    annotate(geom = 'label', x = 90, y= 0.3 ,
             label = expression(P[y]==frac(s[y],1+s[y])*","~s[y]==exp(s[1]*y[c]*frac(y-y[c],y))),
             hjust = 0, vjust =1)+
    ylab(expression(
        "Competence to Grow (" * P[y] * ")")) +
    xlab('Accumulated Chill Portions (y)') +
    geom_line(aes(y = med)) +
    theme_bw(base_size = 15)
  
  layout <- "
AB
CC
"
  p_ycpop + p_zcpop + p_transition +
    plot_layout(design = layout)
}
