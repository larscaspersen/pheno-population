
library(patchwork)

source('R/draw_distributions.R')
source('R/data_io.R')
source('R/population_sampling.R')
source('R/phenoflex_transition-function.R')

par_population <- load_params('parameters/topaz_population_normal.yaml')

draw_dist_plots(n = par_population$distribution$n_pop,
                par = par_population$par,
                yc_sd = par_population$distribution$yc_sd,
                zc_sd = par_population$distribution$zc_sd,
                dist_chill = par_population$distribution$dist_chill,
                dist_heat = par_population$distribution$dist_heat,
                add_par = par_population$distribution$add_par,
                seed = par_population$distribution$seed)
#have to adjust the label if it is not normal
#fine tune boxes, they are all a little bit off now

ggsave('fig/manuscript/s1_distributions.jpeg',
       height = 15, width = 20, units = 'cm', device = 'jpeg')
