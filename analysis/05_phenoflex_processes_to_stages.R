

source('R/data_io.R') #load data
forcing_obs <- helper_prepare_obs_data()

#read parameters
par_pop <- load_params('parameters/topaz_population_normal.yaml')
par_single <- load_params('parameters/topaz_stage-specific.yaml')

#convert from new to old format
par_single$par <- par_single$par |> LarsChill::convert_parameters()
names(par_single$par) <- LarsChill::phenoflex_parnames_old

#use budbreak heat requirement
par_single$par[2] <- par_single$zc_stages$full_bloom

#source('R/run_bud-population-model.R')
source('R/run_population_model_evalpheno.R')

#run population model
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

pop_out_single <-  evalpheno::phenoflex_population(par = par_single$par,
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

#-----------------------#
# development stages ####
#-----------------------#
source('R/process_to_stages.R')

stages <- get_stage_pop(mod_out = pop_out,
                        yc_pop = pop_out$yc_pop, 
                        zc_pop = pop_out$zc_pop, 
                        s1 = par_pop$par[3],
                        p_zc_eco = par_pop$scale_yc_budbreak)


fill_scheme <- c('#0571b0','grey85', '#d7191c', '#a6d96a', '#7b3294' )


plot_stages_pop <- plot_stages_bars(stages,next_stage_on_bottom = TRUE, 
                                    obs_data = forcing_obs,
                                    c('Endodormancy', 
                                      'Transition:\nEndo- to Ecodormancy',
                                      'Ecodormancy', 
                                      'Visible Development', 
                                      'Bloom')) +
  theme_classic(base_size = 13) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     limits = c(0, 1.001),
                     breaks = seq(0, 1, by = 0.25), expand = c(0,0)) +
  theme_classic(base_size = 13) + 
  ylab('Share of Buds in Developmental Stage (%)')+
  theme(legend.position = 'bottom',
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank()) 
plot_stages_pop
ggsave('fig/population_plot-inv-order.jpeg',
       height = 12.5, width =16.5, units = 'cm', device = 'jpeg')


#data for the manuscript text

#need continuous jdays for the summary
yday_cont <- ifelse(forcing_obs$temp_df$Year == max(forcing_obs$temp_df$Year), 
                    yes =  forcing_obs$temp_df$JDay,
                    no = forcing_obs$temp_df$JDay - 365)

stages %>% 
  apply(MARGIN = 2, FUN = function(x){
    #get levels
    x_levels <-unique(x)
    #get index when levels first appear
    index_x <- sapply(x_levels, FUN = function(x_i){
      return(min(which(x == x_i)))
    })
    return(index_x)
  }) %>% 
  as.data.frame() %>% 
  mutate(stage = 1:nrow(.)) %>% 
  pivot_longer(cols = -stage) %>% 
  group_by(stage) %>% 
  summarise(min = yday_cont[min(value)] %>% as.Date(origin = '2021-12-31'),
            med = yday_cont[median(value)]%>% as.Date(origin = '2021-12-31'),
            max = yday_cont[max(value)]%>% as.Date(origin = '2021-12-31'))


#stages for stage-specific model
#make the same plot for single-stage model as well
p_scale <- par_single$zc_stages$budbreak / par_single$zc_stages$full_bloom
par_single$par[2] <- par_single$zc_stages$full_bloom

stages_single <- get_stage_pop(mod_out = pop_out_single,
                               yc_pop = par_single$par[1], 
                               zc_pop = par_single$par[2], 
                               s1 = par_single$par[3],
                               p_zc_eco = p_scale)

plot_stage_single <- plot_stages_bars(stages = stages_single,
                                      obs_data = forcing_obs,
                                      next_stage_on_bottom = TRUE, 
                                      stage_names = c('Endodormancy', 
                                                      'Transition:\nEndo- to Ecodormancy',
                                                      'Ecodormancy', 
                                                      'Visible Development', 
                                                      'Bloom')) +
  theme_classic(base_size = 13) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     limits = c(0, 1.001),
                     breaks = seq(0, 1, by = 0.25), expand = c(0,0)) +
  theme_classic(base_size = 13) + 
  ylab('')+
  theme(legend.position = 'bottom',
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
plot_stage_single


library(patchwork)
library(ggtext)
library(grid)

layout <- "
A
A
A
A
A
B
"
(plot_stages_pop + 
    ggtitle('Bud-Population Model') +
    ylab('Share of Buds in\nDevelopmental Stage (%)')) / (plot_stage_single+
                                                            ggtitle('Stage-Specific Model')) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(design = layout,
              guides = 'collect') &
  theme(legend.position = 'bottom')


ggsave('fig/manuscript/f4_phenostages_both-models.jpeg',
       height = 14, width =18, units = 'cm', device = 'jpeg')
ggsave('fig/manuscript/f4_phenostages_both-models.pdf',
       height = 14, width =18, units = 'cm', device = grDevices::cairo_pdf)
