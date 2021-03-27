# identify bbs stops and routes in the ebird data
# based on several 3 min point counts by the same person, on the same day, 
# within limits of time and distance in the sequence


library(lubridate)
library(sf)
library(tidyverse)
library(lwgeom)

# name of data extraction
data_tag <- "mayjune_201718_bcr27"

# folder to save training and validation datasets
data_proc_folder <- "data_proc/"

# processed data
# they will only be complete counts, so here we use the reduced dataset
eb_zf <- str_glue("{data_proc_folder}/ebd_{data_tag}_zf.csv") %>% 
  read_csv()

# potential bbs checklists: at least 40 complete, 3 min stationary counts 
# from same observer on same day

bbs_candidates <- eb_zf %>%
  filter(all_species_reported, 
         protocol_type == "Stationary", 
         duration_minutes == 3,
         yday(observation_date) > 134) %>% 
  distinct(checklist_id, observer_id, latitude, longitude, 
           observation_date, time_observations_started) %>% 
  group_by(observer_id, observation_date) %>% 
  mutate(n_checklists = n()) %>% 
  ungroup() %>% 
  filter(n_checklists >= 40)


# calculate time and distance between adjacent stops
time_dist_between_stops <- function(df) {
  # time
  df <- df %>% 
    arrange(time_observations_started) %>% 
    mutate(time_bt_stops = c(0, diff(time_observations_started)))
  
  # distance
  df$dist_bt_stops <- df %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    arrange(time_observations_started) %>% 
    mutate(next_pt = geometry[row_number() + 1],
           d = st_distance(geometry, next_pt, by_element = TRUE)) %>% 
    pull(d) %>% 
    as.numeric() %>% 
    head(., length(.) - 1) %>% 
    c(0, .)
    
  
  return(df)
}


bbs_stop_diff <- bbs_candidates %>% 
  arrange(observer_id, observation_date, time_observations_started) %>% 
  group_by(observer_id, observation_date) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map(data, time_dist_between_stops)) %>% 
  unnest(cols = data)

# can now identify routes based on distance or time threshold
bbs_stop_diff$time_bt_stops %>% 
  keep(~ . < 0.5) %>% 
  hist()
bbs_stop_diff$dist_bt_stops %>% 
  keep(~ . < 3000) %>% 
  hist()
# these could be suitable thresholds
t_thresh <- 0.33 # 20 minutes
d_thresh <- 1500 # 1.5 km

bbs_routes <- bbs_stop_diff %>% 
  # apply thresholds
  mutate(new_route = time_bt_stops > t_thresh & 
           dist_bt_stops > d_thresh) %>% 
  group_by(observer_id, observation_date) %>% 
  arrange(observer_id, observation_date, time_observations_started) %>% 
  # assign routes
  mutate(route = cumsum(new_route)) %>% 
  select(-new_route) %>% 
  # stops in route
  group_by(observer_id, observation_date, route) %>% 
  mutate(n_stops = n()) %>% 
  ungroup()

  bbs_routes$route_id <- group_indices(bbs_routes, observer_id, observation_date, route)


# save to file
bbs_routes %>% 
  count(route_id) %>% 
  filter(n > 30) %>% 
  left_join(bbs_routes) %>%
  select(checklist_id, route_id) %>%
  write_csv(path = paste0(data_proc_folder, "bbs_route_id.txt"))

