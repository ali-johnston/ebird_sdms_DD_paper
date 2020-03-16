# identify bbs stops and routes in the ebird data
# based on several 3 min point counts by the same person, on the same day, 
# within limits of time and distance in the sequence


library(lubridate)
library(sf)
library(tidyverse)

# name of data extraction
data_tag <- "mayjune_201718_bcr27"

# folder to save training and validation datasets
data_proc_folder <- "data_proc/"

# processed data
# they will only be complete counts, so here we use the reduced dataset
eb_zf <- str_glue("{data_proc_folder}/ebd_{data_tag}_zf.csv") %>% 
  read_csv()

# potential bbs checklists: at least 40 complete, 3 min stationary counts
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
# note: ali has 3 additional checklists here than are 3 min traveling counts

# identify nearest neighbours within each group
find_nearest_neighbour <- function(df){
  df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
  
  dist_mat <- st_distance(df_sf) %>% round()
  df$min_dist <- apply(dist_mat, 1, function(x){min(x[x > 0])})
  df$n_near1 <- apply(dist_mat, 1, function(x) {sum((x[x > 0]) < 3300)}) - 1
  df$n_near2 <- apply(dist_mat, 1, function(x) {sum((x[x > 0]) < 6700)}) - 1
  df$n_near3 <- apply(dist_mat, 1, function(x) {sum((x[x > 0]) < 10000)}) - 1
  
  return(df)
}

bbs_nn <- bbs_candidates %>% 
  group_by(observer_id, observation_date) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map(data, find_nearest_neighbour)) %>% 
  unnest(cols = data)

# bbs checklists must have at least: 
# 1 checklist within 2 km
# 2 within 2 miles
# 5 within 2 miles
# 9 within 3 miles
dist_thresh <- 2000
threshold1 <- 2
threshold2 <- 5
threshold3 <- 9
bbs_checklists <- bbs_nn %>% 
  filter(min_dist < 2000,
         n_near1 > threshold1,
         n_near2 > threshold2,
         n_near3 > threshold3)

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


# # identify weird cases
# bbs_routes %>% 
#   count(observer_id, observation_date, route) %>% 
#   filter(n < 40) %>% 
#   inner_join(bbs_routes, by = c("observer_id", "observation_date", "route")) %>% 
#   arrange(observer_id, observation_date, time_observations_started)
# # definitely some cases where observer throws in a 3 min bbs-like count 
# # before, after, or between routes, e.g.
# # bbs stop: https://ebird.org/checklist/S36974226
# # bbs-like checklist after route: https://ebird.org/checklist/S36967837


# ggplot(bbs_routes, aes(x = longitude, y = latitude, colour = as.factor(route))) + 
#   geom_point() + scale_colour_hue()

# ggplot(bbs_routes, aes(x = longitude, y = latitude, colour = as.factor(route))) + 
#   geom_point() + scale_colour_hue() + xlim(-88, -85) + ylim(30, 31)

# ggplot(bbs_routes, aes(x = longitude, y = latitude, colour = as.factor(route))) + 
#   geom_point() + scale_colour_hue() + xlim(-84, -82) + ylim(29.7, 31)


# # digging into cases with >50 point counts on a single route
# bbs_routes %>%
#   count(observer_id, observation_date, route) %>% 
#   filter(n > 52) %>% 
#   inner_join(bbs_routes, by = c("observer_id", "observation_date", "route")) %>% 
#   arrange(observer_id, observation_date, time_observations_started)

# bbs_routes %>%
#   filter(observer_id=="obs115524", observation_date == "2017-05-26") %>%
#   ggplot(aes(x = longitude, y = latitude, colour = as.factor(route))) + 
#   geom_point() + scale_colour_hue()

# test <- bbs_routes %>%
#   filter(observer_id=="obs115524", observation_date == "2017-05-26")

# nrow(test)
# test %>% select(longitude, latitude) %>% distinct() %>% nrow()
# # so repeat 3 min point counts at the same places

