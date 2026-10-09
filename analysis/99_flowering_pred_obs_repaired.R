# Replacement for the commented predicted/observed block in 05_flowering.R.
# Source after creating pop_bloom_kob, pop_bb_kob, pop_bloom_cka, single_kob_df,
# kob_obs, cka_obs, cka_season and the converted (kinetic) par_single.
# IMPORTANT: restore par_pop$par['zc'] <- zc_pop_org BEFORE computing
# pop_bloom_cka; otherwise that object contains budbreak predictions.
# Requires the tidyverse and chillR used by 00_stetup.R.

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
ggsave('fig/pred-obs-cka-kob-topaz_withsinglemode_text.jpeg', plot = p_o_plot_annotated,
       height = 16, width = 19, units = 'cm', device = 'jpeg')
