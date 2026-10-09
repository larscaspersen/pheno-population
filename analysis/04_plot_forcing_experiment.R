source('R/data_io.R') #load data
forcing_obs <- helper_prepare_obs_data()

#read parameters
par_pop <- load_params('parameters/topaz_population_normal.yaml')
par_single <- load_params('parameters/topaz_stage-specific.yaml')

#convert from new to old format
par_single$par <- par_single$par |> evalpheno::characteristic_to_kinetic()
names(par_single$par) <- evalpheno::phenoflex_parnames_kinetic

#use budbreak heat requirement
par_single$par[2] <- par_single$zc_stages$budbreak


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

pop_out_single <- evalpheno::phenoflex_population(par = par_single$par,
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

source('R/plots_forcing.R')

#adjust how the names are shown in the facet boxes
labels <- c(
  "\u2460","\u2461","\u2462","\u2463","\u2464","\u2465","\u2466","\u2467",
  "\u2468","\u2469","\u246A","\u246B","\u246C","\u246D","\u246E","\u246F"
)

cutting_annotation_df <- data.frame(
  yday_plot = ifelse(forcing_obs$jday_cut > 220, yes = forcing_obs$jday_cut -365, no = forcing_obs$jday_cut ),
  label = labels[1:length(forcing_obs$jday_cut)]
)

#overwrite the labels
forcing_obs$jday_name_label <- factor(forcing_obs$jday_name,
                                       levels = forcing_obs$jday_name,
                                       labels = paste0(cutting_annotation_df$label,' (',
                                                       forcing_obs$jday_name, ')'))
#need to correspond to the labels in observations table
forcing_obs$exp_obs <- forcing_obs$exp_obs %>% 
  mutate(jday_fact = factor(jday_fact,
                            levels = levels(jday_fact),
                            labels = forcing_obs$jday_name_label))

#--------------#
#needed to handle the special characters in the plot

library(systemfonts)
library(showtext)

fonts <- system_fonts()

symbol_font <- fonts[
  fonts$family == "Segoe UI Symbol",
]

sysfonts::font_add(
  family = "symbols",
  regular = symbol_font$path[1]
)
#-----------------------------#


plot_forcing <- helper_plot_forcing_exp_test(model_res_list = list(pop_out_single,
                                                   pop_out),
                             obs = forcing_obs$exp_obs, 
                             jday_cut = forcing_obs$jday_cut, 
                             jday_name = forcing_obs$jday_name_label,
                             annotate_performance = FALSE,
                             legend_names = c('PhenoFlex (single-stage)',
                                              'PhenoFlex (population)'))+
  scale_color_manual(values = c("Observed" = "black", "PhenoFlex (single-stage)" = "steelblue", "PhenoFlex (population)" = "tomato"),
                     name = "Data Source") +
  scale_linetype_manual(values = c("Observed" = "solid", "PhenoFlex (single-stage)" = "solid", "PhenoFlex (population)" = "dashed"),
                        name = "Data Source")

p_temp <- forcing_obs$temp_df %>% 
  group_by(JDay, Year) %>% 
  summarise(Tmean = mean(Temp)) %>% 
  ungroup() %>% 
  mutate(run_mean = chillR::runn_mean(vec = Tmean, runn_mean = 15),
         yday_plot = ifelse(Year == 2022, yes = JDay, no = JDay - 365)) %>% 
  ggplot(aes(x=yday_plot)) +
  geom_line(aes(y = run_mean)) +
  xlab('Date') +
  ylab('Daily Mean\nTemperature (°C)') +
  scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                     labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
  coord_cartesian(xlim = c(275-365, 122),
                  ylim = c(-2, 15)) +
  geom_rect(aes(ymax = 0, ymin = -Inf, xmin = -Inf, xmax = 1 ), fill = 'grey80') +
  geom_rect(aes(ymax = 0, ymin = -Inf, xmin = Inf, xmax = 1 ), fill = 'grey60') +
  annotate(geom = 'text', x = -90, y = 0-1.2, label = '2021', hjust = 0) +
  annotate(geom = 'text', x = 122, y = 0-1.2, label = '2022', hjust = 1) +
  geom_segment(data = cutting_annotation_df,
               aes(x = yday_plot, xend = yday_plot,
                   y = -Inf, yend = 9),
               linetype = 'dotted')+
  geom_text(data = cutting_annotation_df,
            aes(label = label,
                x = yday_plot, y = 10),
            size = 5, family = "symbols") +
  annotate(geom = 'text', x = mean(c(min(cutting_annotation_df$yday_plot), 
                                     max(cutting_annotation_df$yday_plot))), 
           y = 13, label = 'Cutting experiments') +
  theme_bw(base_size = 13) 

library(patchwork)

layout <- "
A
B
B
B
B
"

(p_temp / plot_forcing + coord_cartesian(xlim = c(0,42))) +
  plot_layout(design = layout) &
  theme_bw(base_size = 13) &
  theme(legend.position = 'bottom',
        panel.grid.minor = element_blank())

showtext::showtext_opts(dpi = 300)
showtext::showtext_auto()
ggsave('fig/manuscript/f2_forcing_experiment.jpg', device = 'jpeg',
       height = 20, width = 19, units = 'cm')
showtext::showtext_auto(FALSE)
showtext::showtext_opts(dpi = 96)

showtext::showtext_auto()
ggsave('fig/manuscript/f2_forcing_experiment.pdf', 
       height = 20, width = 19, units = 'cm', device = grDevices::cairo_pdf)
showtext::showtext_auto(FALSE)







#----------------#
#validation
#----------------#

forcing_obs_val <- helper_prepare_obs_data(sheet =  'T_2020_term+spur')

#run population model
pop_out_val <- evalpheno::phenoflex_population(par = par_pop$par,
                                yc_sd = par_pop$distribution$yc_sd,
                                zc_sd = par_pop$distribution$zc_sd, 
                                dist_chill = par_pop$distribution$dist_chill,
                                dist_heat = par_pop$distribution$dist_heat, 
                                add_distpar = par_pop$distribution$add_par,
                                n = par_pop$distribution$n_pop, 
                                adjust_zc = par_pop$scale_yc_budbreak,
                                temp_df = forcing_obs_val$temp_df, 
                                jday_cut = forcing_obs_val$jday_cut, 
                                basic_output = FALSE,
                                stop_at_zc = FALSE)

pop_out_single_val <- evalpheno::phenoflex_population(par = par_single$par,
                                       yc_sd = par_single$distribution$yc_sd,
                                       zc_sd = par_single$distribution$zc_sd, 
                                       dist_chill = par_single$distribution$dist_chill,
                                       dist_heat = par_single$distribution$dist_heat, 
                                       add_distpar = par_single$distribution$add_par,
                                       n = par_single$distribution$n_pop, 
                                       temp_df = forcing_obs_val$temp_df, 
                                       jday_cut = forcing_obs_val$jday_cut, 
                                       basic_output = FALSE,
                                       stop_at_zc = FALSE)

#adjust how the names are shown in the facet boxes
labels <- c(
  "\u2460","\u2461","\u2462","\u2463","\u2464","\u2465","\u2466","\u2467",
  "\u2468","\u2469","\u246A","\u246B","\u246C","\u246D","\u246E","\u246F"
)

cutting_annotation_df <- data.frame(
  yday_plot = ifelse(forcing_obs_val$jday_cut > 220, 
                     yes = forcing_obs_val$jday_cut -365,
                     no = forcing_obs_val$jday_cut ),
  label = labels[1:length(forcing_obs_val$jday_cut)]
)

#overwrite the labels
forcing_obs_val$jday_name_label <- factor(forcing_obs_val$jday_name,
                                      levels = forcing_obs_val$jday_name,
                                      labels = paste0(cutting_annotation_df$label,' (',
                                                      forcing_obs_val$jday_name, ')'))
#need to correspond to the labels in observations table
forcing_obs_val$exp_obs <- forcing_obs_val$exp_obs %>% 
  mutate(jday_fact = factor(jday_fact,
                            levels = levels(jday_fact),
                            labels = forcing_obs_val$jday_name_label))

plot_forcing_val <- helper_plot_forcing_exp_test(model_res_list = list(pop_out_single_val,
                                                                   pop_out_val),
                                             obs = forcing_obs_val$exp_obs, 
                                             jday_cut = forcing_obs_val$jday_cut, 
                                             jday_name = forcing_obs_val$jday_name_label,
                                             annotate_performance = FALSE,
                                             legend_names = c('PhenoFlex (single-stage)',
                                                              'PhenoFlex (population)'))+
  scale_color_manual(values = c("Observed" = "black", "PhenoFlex (single-stage)" = "steelblue", "PhenoFlex (population)" = "tomato"),
                     name = "Data Source") +
  scale_linetype_manual(values = c("Observed" = "solid", "PhenoFlex (single-stage)" = "solid", "PhenoFlex (population)" = "dashed"),
                        name = "Data Source")

p_temp <- forcing_obs_val$temp_df %>% 
  group_by(JDay, Year) %>% 
  summarise(Tmean = mean(Temp)) %>% 
  ungroup() %>% 
  mutate(run_mean = chillR::runn_mean(vec = Tmean, runn_mean = 15),
         yday_plot = ifelse(Year == 2020, yes = JDay, no = JDay - 365)) %>% 
  ggplot(aes(x=yday_plot)) +
  geom_line(aes(y = run_mean)) +
  xlab('Date') +
  ylab('Daily Mean\nTemperature (°C)') +
  scale_x_continuous(breaks = c(275-365, 306-365, 336-365, 1, 32, 61, 92, 122), 
                     labels = c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')) +
  coord_cartesian(xlim = c(275-365, 122),
                  ylim = c(-2, 15)) +
  geom_rect(aes(ymax = 0, ymin = -Inf, xmin = -Inf, xmax = 1 ), fill = 'grey80') +
  geom_rect(aes(ymax = 0, ymin = -Inf, xmin = Inf, xmax = 1 ), fill = 'grey60') +
  annotate(geom = 'text', x = -90, y = 0-1.2, label = '2019', hjust = 0) +
  annotate(geom = 'text', x = 122, y = 0-1.2, label = '2020', hjust = 1) +
  geom_segment(data = cutting_annotation_df,
               aes(x = yday_plot, xend = yday_plot,
                   y = -Inf, yend = 9),
               linetype = 'dotted')+
  geom_text(data = cutting_annotation_df,
            aes(label = label,
                x = yday_plot, y = 10),
            size = 5, family = 'symbols') +
  annotate(geom = 'text', x = mean(c(min(cutting_annotation_df$yday_plot), 
                                     max(cutting_annotation_df$yday_plot))), 
           y = 13, label = 'Cutting experiments') +
  theme_bw(base_size = 15) 



layout <- "
A
B
B
B
B
"

(p_temp / plot_forcing_val + coord_cartesian(xlim = c(0,42))) +
  plot_layout(design = layout) &
  theme_bw(base_size = 13) &
  theme(legend.position = 'bottom',
        panel.grid.minor = element_blank())

showtext::showtext_opts(dpi = 300)
showtext::showtext_auto()

ggsave('fig/manuscript/s3_forcing_experiment_validation.jpg', device = 'jpeg',
       height = 20, width = 18, units = 'cm')

showtext::showtext_auto(FALSE)
showtext::showtext_opts(dpi = 96)

showtext::showtext_auto()
ggsave('fig/manuscript/s3_forcing_experiment_validation.pdf', device = grDevices::cairo_pdf,
       height = 20, width = 18, units = 'cm')
showtext::showtext_auto(FALSE)


# #-------------------------#
# # forcing: pred vs obs ####
# #-------------------------#
# 
# source('R/forcing_prep.R')
# 
# po_prepared <- prepare_performance_plot_forcing(obs_list = list(forcing_obs,
#                                                                 forcing_obs_val),
#                                                 model_list_calibration = list('PhenoFlex (stage-specific)' = pop_out_single,
#                                                                               'PhenoFlex (population)' = pop_out), 
#                                                 model_list_validation = list('PhenoFlex (stage-specific)' = pop_out_single_val,
#                                                                              'PhenoFlex (population)' = pop_out_val),
#                                                 season_labels = c('Calibration (2021-2022)', 
#                                                                   'Validation (2019-2020)'))
# 
# pos_x <- 0.3
# pos_x_alternative <- 0.6
# pos_y <- 0.95
# pos_y_alternative <- 0.3
# step_y <- -0.07
# step_x <- -0.1
# txt_size <- 3
# 
# 
# forc_perf <- po_prepared$performance_data |> 
#   mutate(nudge_right = case_when(source == 'PhenoFlex (stage-specific)' ~ 0.13, .default = 0.37),
#          nudge_bottom = case_when(source == 'PhenoFlex (stages-specific)' ~ 0, .default = 0.07),
#          pos_x_plot = case_when(month_label %in% c('January', 'February') ~ pos_x_alternative, .default = pos_x),
#          pos_y_plot = case_when(month_label %in% c('January', 'February') ~ pos_y_alternative, .default = pos_y)) 
# 
# annotation_df <- data.frame(month_label = factor(month.name[c(11:12,1:2)])) |> 
#   mutate( pos_x_plot = case_when(month_label %in% c('January', 'February') ~ pos_x_alternative, .default = pos_x),
#           pos_y_plot = case_when(month_label %in% c('January', 'February') ~ pos_y_alternative, .default = pos_y),
#           rmse_label = 'RMSE (%):',
#           rpiq_label = 'RPIQ (-):',
#           bias_label = 'Mean Bias (%):')
# 
# p_o_forcing <- draw_predicted_observed_forcing(po_prepared = po_prepared, month_display = c('Nov', 'Dec', 'Jan', 'Feb'))
# 
# 
# p_o_forcing +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot +nudge_right, y =  pos_y_plot + (step_y*0),
#                 label = rmse, col = source),
#             show.legend = FALSE, size = txt_size) +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot +nudge_right, y =  pos_y_plot + (step_y*1),
#                 label = rpiq, col = source),
#             show.legend = FALSE, size = txt_size) +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot +nudge_right, y =  pos_y_plot + (step_y*2),
#                 label = mean_bias, col = source), 
#             show.legend = FALSE, size = txt_size) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot, y = pos_y_plot, label = rmse_label),
#             size = txt_size, hjust = 1) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot, y = pos_y_plot + step_y, label = rpiq_label),
#             size = txt_size, hjust = 1) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot, y = pos_y_plot + (step_y*2), label = bias_label),
#             size = txt_size, hjust = 1) +
#   guides(color = guide_legend(nrow = 2, ncol = 1),
#          shape = guide_legend(nrow = 2, ncol = 1)) +
#   theme(panel.grid.minor = element_blank())
# 
# ggsave('experimental/fig/manuscript/s4_pred-obs-forcing.jpeg',
#        height = 18, width = 18, units = 'cm',
#        device = 'jpeg')

# pos_x <- 0.00
# pos_y <- 1.0
# step_y <- -0.07
# step_x <- 0.25
# txt_size <- 3
# 
# forc_perf <- po_prepared$performance_data |> 
#   mutate(nudge_right = case_when(source == 'PhenoFlex (stage-specific)' ~ 0.13, .default = 0.37),
#          nudge_bottom = case_when(source == 'PhenoFlex (stage-specific)' ~ 0, .default = 0.07),
#          pos_x_plot = case_when(month_label %in% c('January', 'February') ~ 0.7, .default = pos_x),
#          pos_y_plot = case_when(month_label %in% c('January', 'February') ~ 0.3, .default = pos_y)) 
# 
# annotation_df <- data.frame(month_label = factor(month.name[c(11:12,1:2)])) |> 
#   mutate( pos_x_plot = case_when(month_label %in% c('January', 'February') ~ 0.7, .default = pos_x),
#           pos_y_plot = case_when(month_label %in% c('January', 'February') ~ 0.3, .default = pos_y),
#           rmse_label = 'RMSE (%):',
#           rpiq_label = 'RPIQ (-):',
#           bias_label = 'Mean Bias (%):')
# 
# 
# #different orientation of text
# 
# p_o_forcing +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot, y =  pos_y_plot+ step_y - nudge_bottom,
#                 label = rmse, col = source),
#             show.legend = FALSE, size = txt_size, hjust = 0) +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot + step_x, y =  pos_y_plot + step_y - nudge_bottom,
#                 label = rpiq, col = source),
#             show.legend = FALSE, size = txt_size, hjust = 0) +
#   geom_text(data = forc_perf,
#             aes(x = pos_x_plot + (step_x*2), y =  pos_y_plot + step_y - nudge_bottom,
#                 label = mean_bias, col = source),
#             show.legend = FALSE, size = txt_size, hjust = 0) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot, y = pos_y_plot, label = rmse_label),
#             size = txt_size, hjust = 0) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot+ step_x, y = pos_y_plot, label = rpiq_label),
#             size = txt_size, hjust = 0) +
#   geom_text(data = annotation_df,
#             aes(x = pos_x_plot+(step_x*2), y = pos_y_plot, label = bias_label),
#             size = txt_size, hjust = 0) +
#   guides(color = guide_legend(nrow = 2, ncol = 1),
#          shape = guide_legend(nrow = 2, ncol = 1)) +
#   theme(panel.grid.minor = element_blank())
# 
# ggsave('experimental/fig/pred-obs-forcing-v2.jpeg',
#        height = 18, width = 18, units = 'cm',
#        device = 'jpeg')

