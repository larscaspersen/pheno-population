
source('R/data_io.R') #load data
forcing_obs <- helper_prepare_obs_data()

#read parameters
par_pop <- load_params('parameters/topaz_population_normal.yaml')


#source('R/run_bud-population-model.R')
source('R/run_population_model_evalpheno.R')



#run population model
#it is implemented in the evalpheno package
pop_out <- evalpheno::phenoflex_population(par = par_pop$par,
                                yc_sd = par_pop$distribution$yc_sd,
                                zc_sd = par_pop$distribution$zc_sd, 
                                dist_chill = par_pop$distribution$dist_chill,
                                dist_heat = par_pop$distribution$dist_heat, 
                                add_distpar = par_pop$distribution$add_par,
                                n = par_pop$distribution$n_pop, 
                                adjust_zc = par_pop$scale_yc_budbreak,
                                temp_df = forcing_obs$temp_df, 
                                jday_cut = forcing_obs$jday_cut, 
                                basic_output = FALSE,
                                stop_at_zc = FALSE)

#alternatively, you can also use the convenience function
# run_population_model(temp = forcing_obs$temp_df, 
#                      params = par_pop, 
#                      jday_cut =  forcing_obs$jday_cut, 
#                      basic_output = FALSE,
#                      stop_at_zc = FALSE)


source('R/plots_processes.R') #draw plot

draw_plot_processes(model_output = pop_out, 
                    obs_data = forcing_obs, 
                    yc_pop = pop_out$yc_pop,
                    zc_pop = pop_out$zc_pop, 
                    s1 = par_pop$par[3])

ggsave(filename = 'fig/manuscript/f3_processes.jpeg',
       height = 24,
       width = 20, device = 'jpeg', units = 'cm')
ggsave(filename = 'fig/manuscript/f3_processes.pdf',
       height = 24, width = 20,  units = 'cm', device = grDevices::cairo_pdf)


summarize_processes(obs_data = forcing_obs, 
                    model_out = pop_out, par = par_pop$par)


#---------------------#
#processes of population and single stage model
#---------------------#



par_single <- load_params('parameters/topaz_stage-specific.yaml')
par_single$par <- par_single$par |> LarsChill::convert_parameters()
names(par_single$par) <- LarsChill::phenoflex_parnames_old

#use full bloom for heat
par_single$par[2] <- par_single$zc_stages$full_bloom

pop_out_single <- helper_run_pop_model(par = par_single$par,
                                       yc_sd = par_single$distribution$yc_sd,
                                       zc_sd = par_single$distribution$zc_sd, 
                                       dist_chill = par_single$distribution$dist_chill,
                                       dist_heat = par_single$distribution$dist_heat, 
                                       add_distpar = par_single$distribution$add_par,
                                       n = par_single$distribution$n_pop, 
                                       temp_df = forcing_obs$temp_df, 
                                       jday_cut = forcing_obs$jday_cut, 
                                       basic_output = FALSE,
                                       stop_at_zc = FALSE)

draw_plot_processes_v2(model_output = pop_out, 
                    obs_data = forcing_obs, 
                    yc_pop = pop_out$yc_pop,
                    zc_pop = pop_out$zc_pop, 
                    s1 = par_pop$par[3], 
                    stage_model_output = pop_out_single,
                    yc_stage = par_single$par[1],
                    zc_stage = par_single$par[2],
                    s1_stage = par_single$par[3])
ggsave(filename = 'experimental/fig/manuscript/s2_processes_population-single.png',
       height = 24,
       width = 20, device = 'png', units = 'cm')