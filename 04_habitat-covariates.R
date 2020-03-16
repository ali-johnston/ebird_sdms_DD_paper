# download modis landcover data
# calculate pland within buffer around each checklist
# do good and bad practice checklists at the same time
# create prediction surface

library(raster)
library(velox)
library(fasterize)
library(sf)
library(MODIS)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
# resolve namespace conflicts
select <- dplyr::select


# load data ----


# data folders
data_proc_folder <- "data_proc/"
data_raw_folder <- "data_raw/"

# landcover classes
lc_classes <- read_csv(paste0(data_raw_folder, "/modis/modis_umd_classes.csv"))

# bcr 27 boundary
bcr <- read_sf(paste0(data_proc_folder, "/gis-data.gpkg"), "bcr") %>% 
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
if (length(list.files(tif_dir, "tif$")) < 9) {
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
f_tifs <- list.files(tif_dir, "^modis_umd", full.names = TRUE)
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
  distinct(locality_id, latitude, longitude) %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover)) %>% 
  # add an id column for joining later
  mutate(id = 1:nrow(.))
# locality id must be unique
stopifnot(n_distinct(ebd_pts$locality_id) == length(ebd_pts$locality_id))
# buffer every point to create circular neighborhood
ebd_buffs <- st_buffer(ebd_pts, dist = neighborhood_radius)
# extract values within buffer, using velox for speed
landcover_vlx <- velox(landcover)
landcover_buffer <- landcover_vlx$extract(ebd_buffs, df = TRUE) %>% 
  setNames(c("id", names(landcover))) %>% 
  gather("year", "landcover", -id) %>% 
  mutate(year = as.integer(str_extract(year, "[0-9]+"))) %>% 
  # set values that aren't valid landcover classes to NA
  mutate(landcover = if_else(landcover %in% lc_classes$class,
                             landcover, NA_integer_))
# covert id back to locality id
landcover_buffer <- ebd_buffs %>% 
  st_set_geometry(NULL) %>% 
  inner_join(landcover_buffer, by = "id") %>% 
  select(-id)

# calculate the percent of each landcover class
pland <- landcover_buffer %>% 
  count(locality_id, year, landcover) %>% 
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n)
# fill in implicit missing values with 0s
pland <- complete(pland, 
                  locality_id, year, 
                  landcover = lc_classes$class, 
                  fill = list(pland = 0)) %>% 
  mutate(landcover = paste0("pland_", str_pad(landcover, 2, pad = "0"))) %>% 
  # transform to wide format
  spread(landcover, pland)
# save
write_csv(pland, paste0(data_proc_folder, "modis_pland_loc-year_", data_tag, ".csv"))

# attach back to ebd data by year and location
ebd_pland <- bind_rows(ebd, ebd_all) %>% 
  distinct(checklist_id, locality_id, observation_date) %>% 
  # modis only for 2001-16, for dates outside use closest year
  mutate(year = as.integer(format(observation_date, "%Y")),
         year = if_else(year > 2016, 2016L, year),
         year = if_else(year < 2001, 2001L, year)) %>% 
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
# extract values within buffer, only need 2016
landcover_vlx_curr <- velox(landcover[["y2016"]])
landcover_buffer_curr <- landcover_vlx_curr$extract(r_buffer, df = TRUE) %>% 
  setNames(c("id", "landcover")) %>% 
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
  # fill in implicit missing values with 0s
  complete(id, landcover = lc_classes$class, fill = list(pland = 0)) %>% 
  mutate(landcover = paste0("pland_", str_pad(landcover, 2, pad = "0"))) %>% 
  # transform to wide format
  spread(landcover, pland) %>% 
  mutate(year = 2016L) %>% 
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
