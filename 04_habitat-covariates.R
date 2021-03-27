# download modis landcover data
# calculate pland within buffer around each checklist
# do good and bad practice checklists at the same time
# create prediction surface

library(raster)
library(exactextractr)
library(fasterize)
library(sf)
library(MODIS)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
# resolve namespace conflicts
select <- dplyr::select


# load data ----


# data folders
data_proc_folder <- "data_proc/"
data_raw_folder <- "data_raw/"

# landcover classes
lc_classes <- read_csv(paste0(data_raw_folder, "/modis/modis_umd_classes.csv"))

# bcr 27 boundary
bcr <- read_sf(file.path(data_proc_folder, "gis-data.gpkg"), "bcr") %>% 
  filter(bcr_code == 27)


# name of data extraction
data_tag <- "mayjune_201718_bcr27"


# ebird data
ebd <- read_csv(paste0(data_proc_folder, "data_4_models_", data_tag, ".csv"))

# bad ebird data
ebd_all <- read_csv(paste0(data_proc_folder, "data_all_4_models_", data_tag, ".csv"))



# modis data ----

# download and mosaic modis mcd12q1 v6 landcover data with umd classes
tif_dir <- "data_raw/modis"
if (length(list.files(tif_dir, "tif$")) < 2) {
  tiles <- getTile(bcr)
  # earliest year of ebird data
  ebd_start_year <- format(min(ebd$observation_date), "%Y.01.01")
  tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                  tileH = tiles@tileH, tileV = tiles@tileV,
                  begin = ebd_start_year, end = "2018.12.31", 
                  job = "modis_umd_bcr27") %>% 
    pluck("MCD12Q1.006") %>% 
    unlist()
  # save tifs in project directory
  if (!dir.exists(tif_dir)) {
    dir.create(tif_dir)
  }
  for (i in seq_along(tifs)) {
    yr <- format(as.Date(names(tifs)[i]), "%Y")
    f <- file.path(tif_dir, "modis_umd_{yr}.tif") %>% 
      str_glue()
    file.copy(tifs[i], f)
  }
}

# load annaul landcover layers
f_tifs <- list.files(tif_dir, "^modis_umd.*tif$", full.names = TRUE)
layer_year <- str_extract(f_tifs, "(?<=modis_umd_)[0-9]{4}") %>% 
  paste0("y", .)
landcover <- stack(f_tifs) %>% 
  setNames(layer_year)


# landscape metrics ----

# calculate pland within neighborhood of every unique checklist location
neighborhood_radius <- 2 * ceiling(max(res(landcover))) # ~ 2 modis cells
agg_factor <- 5 # this will produce a 5x5 cell neighbourhood
# ebird data to sf object
ebd_pts <- bind_rows(ebd, ebd_all) %>% 
  # get unique locations
  distinct(year = as.integer(format(observation_date, "%Y")),
           locality_id, latitude, longitude) %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover))
# buffer every point to create circular neighborhood
ebd_buffs <- ebd_pts %>% 
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(locality_id, geometry))

# function to summarize landcover data for all checklists in a given year
calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[paste0("y", yr)]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}
# extract values within buffer
landcover_buffer <- ebd_buffs %>% 
  mutate(pland = map2(year, data, calculate_pland, lc = landcover)) %>% 
  select(year, pland) %>% 
  unnest(cols = pland) %>% 
  # set values that aren't valid landcover classes to NA
  mutate(landcover = if_else(landcover %in% lc_classes$class,
                             landcover, NA_integer_))

# calculate the percent of each landcover class
pland <- landcover_buffer %>% 
  count(locality_id, year, landcover) %>% 
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))
# tranform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>% 
  mutate(landcover = paste0("pland_", str_pad(landcover, 2, pad = "0"))) %>% 
  arrange(landcover) %>% 
  pivot_wider(names_from = landcover,
              values_from = pland, 
              values_fill = list(pland = 0))

# save
write_csv(pland, paste0(data_proc_folder, "modis_pland_loc-year_", data_tag, ".csv"))

# attach back to ebd data by year and location
ebd_pland <- bind_rows(ebd, ebd_all) %>% 
  distinct(checklist_id, locality_id, observation_date) %>% 
  mutate(year = as.integer(format(observation_date, "%Y"))) %>% 
  select(checklist_id, locality_id, year) %>% 
  inner_join(pland, by = c("locality_id", "year")) %>% 
  select(-locality_id, -year)
write_csv(ebd_pland, paste0(data_proc_folder, "modis_pland_checklists_", data_tag, ".csv"))


# prediction surface ----

# template raster, cell size equal to neighborhood size from buffering
# cells = 1 within BCR
r <- raster(landcover) %>% 
  aggregate(agg_factor) %>% 
  fasterize(st_transform(bcr, crs = projection(.)), .) %>% 
  trim() %>% 
  writeRaster(filename = str_glue(data_proc_folder, "modis_{agg_factor}xagg.tif"), 
              overwrite = TRUE)
# extract cell centers and buffer
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  transmute(id = 1:nrow(.))
r_buffer <- st_buffer(r_centers, dist = neighborhood_radius)
# extract values within buffer, only need 2018
landcover_buffer_curr <- landcover[["y2018"]] %>% 
  exact_extract(r_buffer, progress = FALSE) %>% 
  map(~ count(., landcover = value)) %>% 
  tibble(id = r_buffer$id, data = .) %>% 
  unnest() %>% 
  # set values that aren't valid landcover classes to 0
  mutate(landcover = if_else(landcover %in% lc_classes$class,
                             landcover, NA_integer_))
# calculate the percent of each landcover class
pland <- landcover_buffer_curr %>% 
  count(id, landcover) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>%
  filter(!is.na(landcover))
pland_wide <- pland %>% 
  mutate(landcover = paste0("pland_", str_pad(landcover, 2, pad = "0"))) %>% 
  arrange(landcover) %>% 
  pivot_wider(names_from = landcover, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  mutate(year = 2018L) %>% 
  select(id, year, everything())
# bring in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_wide, by = "id")
# save
write_csv(pland_coords, paste0(data_proc_folder, "modis_pland_prediction-surface.csv"))
