t_22 <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                              sheet = 'T_2022_term+spur')

helper_formate_forcing_data <- function(exp_obs, cultivar = 'Topaz', ndays_forcing = 42){
  exp_obs <- exp_obs %>% 
    group_by(Data) %>% 
    mutate(cumsum = cumsum(Bubble_size)) %>% 
    ungroup()
  
  exp_obs$Date <- lubridate::parse_date_time(exp_obs$Data, orders = 'dmy')
  exp_obs$yday <- lubridate::yday(exp_obs$Date)
  
  exp_obs_sub <-  exp_obs %>% 
    select(Days, Date, cumsum, yday) %>% 
    filter(Days <= ndays_forcing)
  
  
  exp_obs_sub$cumsum <- round(exp_obs_sub$cumsum, digits = 2)
  
  cultivar <- cultivar
  year <- exp_obs_sub$Date %>% lubridate::year()
  season <- max(year)
  
  exp_obs_sub %>% 
    mutate(cultivar = cultivar,
           year = year,
           season = season) %>% 
    arrange(Date, Days) %>% 
    return()
}


t_22 <- helper_formate_forcing_data(t_22)

t_20 <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                   sheet = 'T_2020_term+spur') %>% 
  helper_formate_forcing_data()

k_20 <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                           sheet = 'K_2020_term+spur') %>% 
  helper_formate_forcing_data(cultivar = 'Kanzi')

k_22 <- readxl::read_excel('experimental/Time to budbreak data for Sigma.xlsx', 
                           sheet = 'K_2022_term+spur') %>% 
  helper_formate_forcing_data(cultivar = 'Kanzi')

rbind(t_20, t_22, k_20, k_22) %>% 
  write.csv('experimental/forcing_data_formatted.csv', row.names = FALSE)




#also format the forcing hourly obs
s <- read.csv('experimental/hohenheim_aug21-jul22.csv', sep = ';', dec = ',') %>% 
  mutate(Date = lubridate::dmy(Tag),
         Hour = lubridate::hm(Stunde) %>% hour(),
         JDay = lubridate::yday(Date),
         Year = lubridate::year(Date),
         Day = lubridate::day(Date),
         Month = lubridate::month(Date),
         Temp = AVG_TA200,
  ) %>% 
  select(Date, Year, Month, Day, Hour, Temp)
s2 <- read.csv('experimental/hohenheim_aug19-jul20.csv', sep = ';', dec = ',') %>% 
  mutate(Date = lubridate::dmy(Tag),
         Hour = lubridate::hm(Stunde) %>% hour(),
         JDay = lubridate::yday(Date),
         Year = lubridate::year(Date),
         Day = lubridate::day(Date),
         Month = lubridate::month(Date),
         Temp = AVG_TA200,
  ) %>% 
  select(Date, Year, Month, Day, Hour, Temp)

rbind(s, s2) %>% 
  write.csv('experimental/temp_forcing_season_complete.csv', row.names = FALSE)
