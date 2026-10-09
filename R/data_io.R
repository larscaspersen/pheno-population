#process data
helper_prepare_obs_data <- function(sheet = 'T_2022_term+spur',
                                    start_yday = 300,
                                    end_yday = 55){
  
  exp_obs <- readxl::read_excel('Time to budbreak data for Sigma.xlsx', 
                                sheet = sheet) %>% 
    group_by(Data) %>% 
    mutate(cumsum = cumsum(Bubble_size)) %>% 
    mutate(Date = lubridate::parse_date_time(Data, orders = 'dmy'),
           yday = lubridate::yday(Date),
           h_after_cut = Days*24) 
  
  #subsetting
  
  if(start_yday >end_yday){
    exp_obs <- exp_obs[exp_obs$yday >= start_yday | exp_obs$yday <= end_yday, ]
  } else {
    exp_obs <- exp_obs[exp_obs$yday >= start_yday & exp_obs$yday <= end_yday, ]
  }
  
  if(grepl(pattern = '2020', x = sheet)) temp_file <- 'data/hohenheim_aug19-jul20.csv'
  if(grepl(pattern = '2022', x = sheet)) temp_file <- 'data/hohenheim_aug21-jul22.csv'
  #read temperature data
  s <- read.csv(temp_file, sep = ';', dec = ',') %>% 
    mutate(Date = lubridate::dmy(Tag),
           Hour = lubridate::hm(Stunde) %>% hour(),
           JDay = lubridate::yday(Date),
           Year = lubridate::year(Date),
           Temp = AVG_TA200,
    ) %>% 
    select(Temp, JDay, Year)
  
  #-----------------------------------#
  #cutting experiments
  
  #the function expects to get already the position at which the forcing experiment should take place
  jday_cut <- unique(exp_obs$yday)
  jday_name <- lubridate::stamp("03 Nov", orders = '%d %b', locale = 'C')(unique(exp_obs$Date))
  #jday_name <- c('Nov 03', 'Nov 17', 'Dec 01', 'Dec 15', 'Dec 29', 'Jan 12', 'Jan 26', 'Feb 09', 'Feb23')
  i_cut <-purrr::map_int(jday_cut, function(x){
    floor(median(which(x == s$JDay)))
  })  
  #c++ starts counting at zero, correct for that
  i_cut <- i_cut -1
  
  #add missing days to antons table
  miss_df <- exp_obs %>% 
    group_by(Date, yday) %>% 
    summarise(max = max(cumsum),
              miss = 1 - max) %>% 
    ungroup() %>% 
    mutate(h_after_cut = 51 *24,
           Bubble_size = miss,
           Days = 51,
           cumsum = 1) %>% 
    select(Days, Date, yday, h_after_cut, Bubble_size, cumsum)
  
  exp_obs <- exp_obs %>% 
    select(Days, Date, yday, h_after_cut, Bubble_size, cumsum) %>% 
    rbind(miss_df) %>%  
    mutate(jday_fact = factor(yday, levels = jday_cut, 
                              labels = jday_name)) 
  
  return(list(
    exp_obs = exp_obs,
    jday_cut = jday_cut,
    jday_name = jday_name,
    i_cut = i_cut,
    temp_df = s
    
  ) )
}

# R/data_io.R

load_kob_season <- function(path = 'data/Ravensburg_hourly_temp_fixed.csv',
                            years = 2004:2022) {
  read.csv(path) %>%
    genSeasonList(years = years,mrange = c(8,5)) %>%
    setNames(years)
}

load_kob_bloom <- function(path = 'data/Ravensburg_bloom_dates.csv',
                           variety = 'Topaz') {
  read.csv(path) %>%
    filter(variety == !!variety) %>%
    mutate(firstbloom = lubridate::yday(firstbloom),
           fullbloom  = lubridate::yday(fullbloom),
           budbreak   = lubridate::yday(budbreak))
}

load_cka_season <- function(path = 'data/cka_clean.csv',
                            latitude = 50.7,
                            years) {
  read.csv(path) %>%
    stack_hourly_temps(latitude = latitude) %>%
    purrr::pluck('hourtemps') %>%
    chillR::genSeasonList(years = years) %>%
    setNames(years)
}

load_topaz_flowering_cka <- function(path = 'data/topaz_flowering_cka.csv',
                                     end_of_bloom_share = 0.9) {
  read.csv(path) %>%
    mutate(share = case_match(name,
                              'begin_flowering_f5' ~ 0.1,
                              'flowering_f50'      ~ 0.5,
                              'end_flowering'      ~ end_of_bloom_share),
           name = case_match(name,
                             'begin_flowering_f5' ~ 'First Bloom (10% open)',
                             'flowering_f50'      ~ 'Full Bloom (50% open)',
                             'end_flowering'      ~ 'End of Bloom (most petals fallen)'))
}


load_params <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$par <- unlist(cfg$par)   # named numeric vector, e.g. par["yc"]
  cfg
}

# convert a named par list/vector into the raw ordered vector
# your C++/PhenoFlex_pop_slim call expects
par_to_vector <- function(par) {
  c(par["yc"], par["zc"], par["s1"], par["Tu"],
    par["E0"], par["E1"], par["A0"], par["A1"],
    par["Tf"], par["Tc"], par["Tb"], par["slope"])
}
