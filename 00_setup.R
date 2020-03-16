# setup for auk and MODIS packages
# download and prepare gis layers for plotting
# only needs to be run once

library(sf)
library(auk)
library(MODIS)
library(rvest)
library(rnaturalearth)
library(dplyr)
library(readr)

# auk package setup ----

auk_set_ebd_path("******")


# modis setup ----

# package setup
EarthdataLogin(usr = "******", pwd = "******")


# gis data ----

f_ne <- "data_proc/gis-data.gpkg"
# bcrs
tmp_dir <- tempdir()
tmp_bcr <- file.path(tmp_dir, "bcr.zip")
paste0("https://www.birdscanada.org/research/gislab/download/", 
       "bcr_terrestrial_shape.zip") %>% 
  download.file(destfile = tmp_bcr)
unzip(tmp_bcr, exdir = tmp_dir)
bcr <- file.path(tmp_dir, "BCR_Terrestrial_master_International.shp") %>% 
  read_sf() %>% 
  select(bcr_code = BCR, bcr_name = LABEL) %>% 
  mutate(bcr_name = recode(bcr_name, 
                           "Lower Great Lakes/ St. Lawrence Plain" = 
                             "Lower Great Lakes/St. Lawrence Plain"))
list.files(tmp_dir, "bcr", ignore.case = TRUE, full.names = TRUE) %>% 
  unlink()

# political boundaries
# bounding box
ne_bbox <- ne_download(scale = 50, category = "physical",
                       type = "wgs84_bounding_box",
                       returnclass = "sf")
# 15 degree graticules
ne_graticules <- ne_download(scale = 50, category = "physical",
                             type = "graticules_15",
                             returnclass = "sf")
# land border with lakes removed
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  st_set_precision(1e6) %>%
  st_union()
# country lines
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf")
# states, north america
ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(adm0_a3, USA = "US", CAN = "CAN")) %>% 
  select(country = adm0_name, country_code = iso_a2)

# output
unlink(f_ne)
write_sf(ne_bbox, f_ne, "ne_bbox")
write_sf(ne_graticules, f_ne, "ne_graticules")
write_sf(ne_land, f_ne, "ne_land")
write_sf(ne_country_lines, f_ne, "ne_country_lines")
write_sf(ne_state_lines, f_ne, "ne_state_lines")
write_sf(bcr, f_ne, "bcr")