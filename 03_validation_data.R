# take the processed ebird data and separate into train, test1, test2
# train = 2018 non-bbs data
# test1 = 2018 bbs data 
# test2 = 2017 non-bbs data


library(tidyverse)
library(lubridate)
library(rgdal)
library(sf)

# custom functions
walk(list.files("R", full.names = TRUE), source)

# folder to save training and validation datasets
data_proc_folder <- "data_proc/"

# folder to save plots of the data
figure_folder <- "figures/raw_data/"
dir.create(figure_folder)

# name of data extraction
data_tag <- "mayjune_201718_bcr27"

# ####################################################################
# READ IN PROCESSED EBIRD DATA

eb_zf_loc <- paste0(data_proc_folder, "ebd_", data_tag, "_zf.csv")
eb_zf <- read.csv(eb_zf_loc) %>%
          select(-species_observed) %>%
          spread(species_code, observation_count) %>%
          filter(yday(observation_date)>134)

# read in badly processed data
eb_all_zf_loc <- paste0(data_proc_folder, "ebd_", data_tag, "_zf_all.csv")
eb_all_zf <- read.csv(eb_all_zf_loc) %>%
          select(-species_observed) %>%
          spread(species_code, observation_count) %>%
          filter(yday(observation_date)>134)



# ####################################################################
# assign the suspected bbs counts

bbs_routes <- read_csv(paste0(data_proc_folder, "bbs_route_id.txt")) %>%
            mutate(bbs = 1)


# --------------------------------------------------------------------
# merge back in with main dataset

eb_zf_type <- eb_zf %>%
            mutate(checklist_id = as.character(checklist_id)) %>%
            left_join(bbs_routes) %>%
            mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
            mutate(type = case_when(
                  bbs==0 & year(observation_date)==2018 ~ "train",
                  bbs==0 & year(observation_date)==2017 ~ "test_2017",
                  bbs==1 & year(observation_date)==2018 ~ "test_bbs",
                  TRUE                                   ~ "other"))

# copy across to the 'bad' dataset
eb_all_zf_type <- eb_all_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(select(eb_zf_type, checklist_id, type)) %>%
        mutate(type = ifelse(is.na(type), 
                            ifelse(year(observation_date)==2018, "train", "other"), 
                            type))


# ####################################################################
# WRITE TO FILES

data_loc <- paste0(data_proc_folder, "data_4_models_", data_tag, ".csv")
write_csv(eb_zf_type, data_loc)

data_all_loc <- paste0(data_proc_folder, "data_all_4_models_", data_tag, ".csv")
write_csv(eb_all_zf_type, data_all_loc)



# ####################################################################
# SUMMARISE NUMBERS

# all data
nrow(eb_all_zf_type)

# split by year
table(year(eb_all_zf_type$observation_date))

# split by year and bbs
eb_all_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(bbs_routes) %>%
        mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
        mutate(yr = year(observation_date)) %>%
        select(bbs, yr) %>%
        table()

# 2017 non-bbs
eb_all_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(bbs_routes) %>%
        mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
        mutate(yr = year(observation_date)) %>%
        filter(yr == 2017, bbs == 0) %>%
        nrow()

# 2017 non-bbs, complete, effort filtered
eb_all_zf_type %>%
        filter(type == "test_2017") %>%
        filter(all_species_reported) %>%
        filter(protocol_type %in% c("Stationary", "Traveling")) %>%
        mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
        filter(effort_distance_km <= 5,
               duration_minutes <= 5 * 60,
               number_observers <= 10) %>%
        nrow()



# ####################################################################
# PLOT DIFFERENT DATASETS

map_proj <- st_crs(102003)
# borders
f_gpkg <- paste0(data_proc_folder, "gis-data.gpkg")
ne_land <- read_sf(f_gpkg, "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf(f_gpkg, "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf(f_gpkg, "ne_state_lines") %>% 
  filter(country_code == "US") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf(f_gpkg, "bcr") %>% 
  filter(bcr_code == 27) %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


# --------------------------------------------------------------------
# project the datasets

train_data <- eb_zf_type %>%
        filter(type=="train")

train_data_proj <- as.matrix(cbind(train_data$longitude, train_data$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)

# spatial subsample the training data
train_data_ss <- eb_zf_type %>%
        filter(all_species_reported) %>%
        filter(type=="train") %>%
        mutate(week = lubridate::week(observation_date)) %>%
        hex_sample(spacing = 5,
                         regime = "together", byvar = "week") %>%
    select(checklist_id, sampling_event_identifier, type) %>%
    mutate(selected = 1) %>%
    left_join(eb_zf_type)

train_data_ss_proj <- as.matrix(cbind(train_data_ss$longitude, train_data_ss$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


# test data bbs
test_data_bbs <- eb_zf_type %>%
        filter(type=="test_bbs")

test_data_bbs_proj <- as.matrix(cbind(test_data_bbs$longitude, test_data_bbs$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)

# test data 2017
test_data_2017 <- eb_zf_type %>%
        filter(type=="test_2017")

test_data_2017_proj <- as.matrix(cbind(test_data_2017$longitude, test_data_2017$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


plot_name <- paste0(figure_folder, "train_test_maps_", Sys.Date(), ".png")
png(plot_name, width = 14, height = 12, units="cm", pointsize=9, res=300)

  par(mfrow=c(2,2), mar = c(0.5, 0.5, 0.5, 0.5))

  # ----- training data
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(train_data_proj, add = TRUE, pch=16, cex=0.2)
  title(main = "A", line = -1)


  # ----- spatial subsampled training data 
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(train_data_ss_proj, add = TRUE, pch=16, cex=0.2)
  title(main = "B", line = -1)



  # ----- 2017 test data 
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(test_data_2017_proj, add = TRUE, pch=16, cex=0.2)
  title(main = "C", line = -1)


  # ----- BBS stops
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(test_data_bbs_proj, add = TRUE, pch=16, cex=0.2)
  title(main = "D", line = -1)

dev.off()







