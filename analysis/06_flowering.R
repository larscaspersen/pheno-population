library(chillR)
library(evalpheno)
library(tidyverse)
library(patchwork)

source('R/data_io.R') #load data
kob_season <- load_kob_season()
kob_obs <- load_kob_bloom()
cka_obs <- load_topaz_flowering_cka()
cka_season <- load_cka_season(years = unique(cka_obs$year))

#read parameters
par_pop <- load_params('parameters/topaz_population_normal.yaml')
par_single <- load_params('parameters/topaz_stage-specific.yaml')

#convert from new to old format
par_single$par <- par_single$par |> evalpheno::convert_parameters()
names(par_single$par) <- evalpheno::phenoflex_parnames_old

#use budbreak heat requirement
par_single$par[2] <- par_single$zc_stages$full_bloom

#source('R/run_bud-population-model.R')
source('R/run_population_model_evalpheno.R')

source('R/plots_flowering.R')

pop_bloom_kob <- purrr::map(kob_season, function(s1){
  bloom <- evalpheno::phenoflex_population(par = par_pop$par, 
                                yc_sd = par_pop$distribution$yc_sd, 
                                zc_sd =  par_pop$distribution$zc_sd, 
                                jday_cut = NULL, 
                                temp_df = s1,
                                dist_chill = par_pop$distribution$dist_chill,
                                dist_heat = par_pop$distribution$dist_heat,
                                add_distpar = par_pop$distribution$add_par,
                                n= par_pop$distribution$n_pop,
                                basic_output = TRUE,
                                stop_at_zc = TRUE) %>% 
    purrr::pluck('bloomindex') %>% 
    purrr::map_dbl(helper_bloomint_to_jday, x = s1) %>% 
    return()
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  setNames(names(kob_season)) %>% 
  pivot_longer(cols = everything(), names_to = 'year')

zc_pop_org <- par_pop$par[2]
par_pop$par[2] <- zc_pop_org * par_pop$scale_yc_budbreak

pop_bb_kob <- purrr::map(kob_season, function(s1){
  bloom <- helper_run_pop_model(par = par_pop$par, 
                                yc_sd = par_pop$distribution$yc_sd, 
                                zc_sd =  par_pop$distribution$zc_sd, 
                                jday_cut = NULL, 
                                temp_df = s1,
                                dist_chill = par_pop$distribution$dist_chill,
                                dist_heat = par_pop$distribution$dist_heat,
                                add_distpar = par_pop$distribution$add_par,
                                n= par_pop$distribution$n_pop,
                                basic_output = TRUE,
                                stop_at_zc = TRUE) %>% 
    purrr::pluck('bloomindex') %>% 
    purrr::map_dbl(helper_bloomint_to_jday, x = s1) %>% 
    return()
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  setNames(names(kob_season)) %>% 
  pivot_longer(cols = everything(), names_to = 'year')


par_pop$par["zc"] <- zc_pop_org
pop_bloom_cka <- purrr::map(cka_season, function(s1){
  bloom <- helper_run_pop_model(par = par_pop$par, 
                                yc_sd = par_pop$distribution$yc_sd, 
                                zc_sd =  par_pop$distribution$zc_sd, 
                                jday_cut = NULL, 
                                temp_df = s1,
                                dist_chill = par_pop$distribution$dist_chill,
                                dist_heat = par_pop$distribution$dist_heat,
                                add_distpar = par_pop$distribution$add_par,
                                n= par_pop$distribution$n_pop,
                                basic_output = TRUE,
                                stop_at_zc = TRUE) %>% 
    purrr::pluck('bloomindex') %>% 
    purrr::map_dbl(helper_bloomint_to_jday, x = s1) %>% 
    return()
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  setNames(names(cka_season)) %>% 
  pivot_longer(cols = everything(), names_to = 'year')


#prediction single model kob, cka
#add prediction of the standard model
par_single$par[2] <- par_single$zc_stages$first_bloom

model <- evalpheno::pheno_model('phenoflex', chill = evalpheno::chill_dynamic('kinetic'))


fistbloom_single_kob <- evalpheno::predict_phenology(model = model,
                                                     weather = kob_season, parameters = par_single$par)

par_single$par[2] <- par_single$zc_stages$full_bloom
fullbloom_single_kob <- evalpheno::predict_phenology(model = model,
                                                     weather = kob_season, parameters = par_single$par)


par_single$par[2] <- par_single$zc_stages$budbreak
budbreak_single_kob <- evalpheno::predict_phenology(model = model,
                                                    weather = kob_season, parameters = par_single$par)


single_kob_df <- data.frame(share = c(rep(0.1, length(fistbloom_single_kob)),
                                    rep(0.5, length(fistbloom_single_kob)),
                                    rep(NA, length(fistbloom_single_kob))),
                          pred = c(fistbloom_single_kob, fullbloom_single_kob, budbreak_single_kob),
                          year = c(names(fistbloom_single_kob), names(fullbloom_single_kob), 
                                   names(budbreak_single_kob)),
                          stage = rep(c('Firstbloom', 'Fullbloom', 'Budbreak'),
                                       each = length(fistbloom_single_kob)))


helper_plot_flowering(bloom_df = pop_bloom_kob, 
                      obs_df = kob_obs, 
                      year_select = 'all',
                      annotate_performance = FALSE) +
  geom_point(data = single_kob_df, 
             aes(x = pred, y = share, col = 'strd'),
             shape = 17) +
  scale_x_continuous(breaks = c(91, 98, 105, 112,119,126),
                     labels = c('Apr-01', 'Apr-08', 'Apr-15','Apr-22', 'Apr-29', 'May-06')) +
  coord_cartesian(xlim = c(95,130)) +
  facet_wrap(~year, ncol = 4) +
  scale_color_manual(name = 'Data Source',
                     breaks = c('Observed', 'strd', 'Modelled'),
                     labels = c('Observed', 'PhenoFlex (single-stage)', 'PhenoFlex (population)'),
                     values = c('black', 'steelblue',  'tomato')) +
  ylab('Share of Buds Flowering (%)') +
  theme(panel.grid.minor = element_blank())
ggsave('fig/manuscript/f5_flowering_kob_population.jpeg',
       height = 16, width = 16, units = 'cm', device = 'jpeg')
ggsave('fig/manuscript/f5_flowering_kob_population.pdf',
       height = 16, width = 16, units = 'cm',
       device = grDevices::cairo_pdf)

#predicted / observed plot
nudge_right_step <- 6

# Keep failed buds in the distribution as right-censored values. A quantile
# beyond the simulated season becomes NA; removing failed buds before taking
# quantiles would change the population fractions used to define each stage.
stage_quantile <- function(value, probability) {
  value[value == 9999] <- Inf
  if (anyNA(value)) return(NA_real_)
  prediction <- quantile(value, probability, names = FALSE)
  if (is.finite(prediction)) prediction else NA_real_
}

kob_obs_long <- kob_obs |>
  transmute(year = as.character(year), budbreak, firstbloom, fullbloom) |>
  pivot_longer(-year, values_to = 'observed') |>
  mutate(name = case_match(name,
                           'budbreak' ~ 'Budbreak',
                           'firstbloom' ~ 'First Bloom (10% open)',
                           'fullbloom' ~ 'Full Bloom (50% open)'))
cka_obs_long <- cka_obs |>
  transmute(year = as.character(year), name, observed = value)
end_of_bloom_share <- unique(cka_obs$share[
  cka_obs$name == 'End of Bloom (most petals fallen)'])
stopifnot(length(end_of_bloom_share) == 1L)

# Raw population objects contain year/value, not observed/predicted/name.
# Summarize them first, then join observations by year AND stage.
pop_bloom_kob_plot <- pop_bloom_kob |>
  group_by(year) |>
  summarise(`First Bloom (10% open)` = stage_quantile(value, 0.1),
            `Full Bloom (50% open)` = stage_quantile(value, 0.5),
            .groups = 'drop') |>
  pivot_longer(-year, values_to = 'predicted') |>
  mutate(year = as.character(year), location = 'Lake Constance') |>
  inner_join(kob_obs_long, by = c('year', 'name'))
pop_bb_sum_plot <- pop_bb_kob |>
  group_by(year) |>
  summarise(predicted = stage_quantile(value, 0.5), .groups = 'drop') |>
  mutate(year = as.character(year), name = 'Budbreak', location = 'Lake Constance') |>
  inner_join(kob_obs_long, by = c('year', 'name'))
pop_bloom_cka_plot <- pop_bloom_cka |>
  group_by(year) |>
  summarise(`First Bloom (10% open)` = stage_quantile(value, 0.1),
            `Full Bloom (50% open)` = stage_quantile(value, 0.5),
            `End of Bloom (most petals fallen)` = stage_quantile(value, end_of_bloom_share),
            .groups = 'drop') |>
  pivot_longer(-year, values_to = 'predicted') |>
  mutate(year = as.character(year), location = 'Rhineland') |>
  inner_join(cka_obs_long, by = c('year', 'name'))

single_bloom_kob_plot <- single_kob_df |>
  transmute(year = as.character(year), predicted = pred,
            name = case_match(stage,
                              'Budbreak' ~ 'Budbreak',
                              'Firstbloom' ~ 'First Bloom (10% open)',
                              'Fullbloom' ~ 'Full Bloom (50% open)'),
            location = 'Lake Constance') |>
  inner_join(kob_obs_long, by = c('year', 'name'))

# The current 05_flowering.R only predicts the single-stage model at KOB.
# Add the missing CKA predictions. No end-of-bloom threshold is defined in the
# single-stage YAML, so that model is compared for first/full bloom only at CKA.
single_model <- evalpheno::pheno_model('phenoflex',
                                       chill = evalpheno::chill_dynamic('kinetic'))
single_cka_par <- par_single$par
single_cka_par['zc'] <- par_single$zc_stages$first_bloom
firstbloom_single_cka <- evalpheno::predict_phenology(single_model,
                                                      weather = cka_season, parameters = single_cka_par)
single_cka_par['zc'] <- par_single$zc_stages$full_bloom
fullbloom_single_cka <- evalpheno::predict_phenology(single_model,
                                                     weather = cka_season, parameters = single_cka_par)
single_bloom_cka_plot <- bind_rows(
  tibble(year = names(firstbloom_single_cka), predicted = unname(firstbloom_single_cka),
         name = 'First Bloom (10% open)'),
  tibble(year = names(fullbloom_single_cka), predicted = unname(fullbloom_single_cka),
         name = 'Full Bloom (50% open)')) |>
  mutate(location = 'Rhineland') |>
  inner_join(cka_obs_long, by = c('year', 'name'))

# Assemble one comparison table for both metrics and figures. Bind rows rather
# than pairing independent bud/year records by their current row position.
comparison_df <- bind_rows(
  bind_rows(pop_bloom_kob_plot, pop_bb_sum_plot, pop_bloom_cka_plot) |>
    mutate(model = 'PhenoFlex (population)'),
  bind_rows(single_bloom_kob_plot, single_bloom_cka_plot) |>
    mutate(model = 'PhenoFlex (single-stage)'))
comparison_complete <- comparison_df |>
  filter(is.finite(observed), is.finite(predicted))
perf_df <- comparison_complete |>
  group_by(model, location, name) |>
  summarise(n_pairs = n(),
            rmse = chillR::RMSEP(predicted = predicted, observed = observed),
            rpiq = chillR::RPIQ(predicted = predicted, observed = observed),
            mean_bias = mean(predicted - observed), .groups = 'drop') |>
  mutate(nudge_right = case_when(
    name == 'Budbreak' ~ 0,
    name == 'First Bloom (10% open)' & location == 'Lake Constance' ~ nudge_right_step,
    name == 'First Bloom (10% open)' ~ 0,
    name == 'Full Bloom (50% open)' & location == 'Lake Constance' ~ nudge_right_step * 2,
    name == 'Full Bloom (50% open)' ~ nudge_right_step,
    name == 'End of Bloom (most petals fallen)' ~ nudge_right_step * 2))
perf_df_kob <- perf_df |> filter(location == 'Lake Constance')
perf_df_cka <- perf_df |> filter(location == 'Rhineland')

mon_start <- c(60, 91, 121)
mon_label <- c('Mar-01', 'Apr-01', 'May-01')
p_o_plot <- ggplot(comparison_complete, aes(x = observed, y = predicted, col = name)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_point() +
  scale_color_manual(name = 'Phenological Stage',
                     values = c('Budbreak' = '#009E73', 'First Bloom (10% open)' = 'tomato',
                                'Full Bloom (50% open)' = 'steelblue',
                                'End of Bloom (most petals fallen)' = '#E69F00'),
                     breaks = c('Budbreak', 'First Bloom (10% open)', 'Full Bloom (50% open)',
                                'End of Bloom (most petals fallen)')) +
  theme_bw(base_size = 15) +
  scale_x_continuous(breaks = mon_start, labels = mon_label) +
  scale_y_continuous(breaks = mon_start, labels = mon_label) +
  guides(color = guide_legend(nrow = 2, ncol = 2)) +
  theme(legend.position = 'bottom', panel.grid.minor = element_blank()) +
  labs(x = 'Observed', y = 'Predicted') +
  facet_grid(model ~ location)

pox_x <- 120
pos_y_start <- 75
step_y <- 5
txt_size <- 3
performance_labels <- perf_df |>
  pivot_longer(c(rmse, rpiq, mean_bias), names_to = 'metric', values_to = 'value') |>
  mutate(label_y = pos_y_start - step_y * (match(metric, c('rmse', 'rpiq', 'mean_bias')) - 1))
p_o_plot_annotated <- p_o_plot +
  annotate('rect', xmin = pox_x - 30, xmax = pox_x + 17,
           ymin = pos_y_start - step_y * 2 - 5, ymax = pos_y_start + 5,
           fill = 'white', color = 'black') +
  geom_text(data = performance_labels,
            aes(x = pox_x + nudge_right, y = label_y, col = name,
                label = format(round(value, 1), nsmall = 1)),
            inherit.aes = FALSE, show.legend = FALSE, size = txt_size) +
  annotate('text', x = pox_x - 4, y = pos_y_start,
           label = 'RMSE (days)', size = txt_size, hjust = 1) +
  annotate('text', x = pox_x - 4, y = pos_y_start - step_y,
           label = 'RPIQ (-)', size = txt_size, hjust = 1) +
  annotate('text', x = pox_x - 4, y = pos_y_start - step_y * 2,
           label = 'Mean Bias (days)', size = txt_size, hjust = 1)

ggsave('fig/pred-obs-cka-kob-topaz_withsinglemode.jpeg', plot = p_o_plot,
       height = 16, width = 18, units = 'cm', device = 'jpeg')
ggsave('fig/manuscript/f6_pred-obs-cka-kob-topaz_withsinglemode_text.jpeg', plot = p_o_plot_annotated,
       height = 16, width = 19, units = 'cm', device = 'jpeg')
ggsave('fig/manuscript/f6_pred-obs-cka-kob-topaz_withsinglemode_text.pdf', 
       plot = p_o_plot_annotated,
       height = 16, width = 19, units = 'cm', device = grDevices::cairo_pdf)


