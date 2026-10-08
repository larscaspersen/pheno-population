Rcpp::sourceCpp("src/phenoflex_pop.cpp")

source('experimental/R/population_sampling.R')

helper_run_pop_model <- function(par, yc_sd, zc_sd, jday_cut, temp_df, n = 100,
                                 adjust_zc = 1, basic_output = FALSE, stop_at_zc = TRUE,
                                 dist_chill = 'normal', dist_heat = 'normal', 
                                 add_distpar = NULL){
  
  #--------------#
  #draw population of yc, zc
  
  yc_pop <- par[1]
  zc_pop <- par[2]
  
  #in case there is a distribution to sample from, sample yc_pop and zc_pop
  if(n > 1){
    dist_out <- helper_sample_dist(n = n, par = par,
                                   yc_sd = yc_sd, zc_sd = zc_sd, dist_chill = dist_chill,
                                   dist_heat = dist_heat, add_par = add_distpar)
    
    yc_pop <- dist_out$yc_pop
    zc_pop <- dist_out$zc_pop
  }
  
  #---------------#
  #identify timepoints of cutting in temperature data
  
  i_cut <-purrr::map_int(jday_cut, function(x){
    floor(median(which(x == temp_df$JDay)))
  })  
  #c++ starts counting at zero, correct for that
  i_cut <- i_cut -1
  
  #----------------#
  #run model
  
  mod_out <- PhenoFlex_pop_slim(temp = temp_df$Temp, 
                                times = seq_along(temp_df$Temp), 
                                yc = yc_pop,
                                zc = zc_pop,
                                i_cut = i_cut,
                                max_days_forcing = 50,
                                forcing_temperature = 23, 
                                s1 = par[3],
                                E0 = par[5],
                                E1 = par[6],
                                A0 = par[7],
                                A1 = par[8],
                                slope = par[12],
                                Tf = par[9],
                                Tb = par[11],
                                Tu = par[4],
                                Tc = par[10],
                                placeholder_fail = 9999,
                                basic_output = basic_output, 
                                adjust_zc_forcing_exp = adjust_zc, 
                                stopatzc = stop_at_zc) 
  
  #attach yc_pop, zc_pop
  mod_out$yc_pop <- yc_pop
  mod_out$zc_pop <- zc_pop
  
  return(mod_out)
  
}