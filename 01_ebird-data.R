# extract and zero-fill ebird data from the raw ebird dataset. 

# three species are included, although only wood thrush is explored in paper
# an incorrectly processed dataset is also prepared here


library(auk)
library(sf)
library(dplyr)
library(purrr)
library(lubridate)
library(readr)
# custom functions
walk(list.files("R", full.names = TRUE), source)

# location of large unzipped eBird data
ebd_data_folder <- "/Volumes/ebird_data/EBD/"
ebd_release <- "May-2019"

# where to save the processed datasets
data_save <- "data_proc/"

# name for this specific data extraction
data_tag <- "mayjune_201718_bcr27"

# zero-filled ebird data ----

# ebd extraction
orig_ebd <- paste0(ebd_data_folder, "ebd_rel", ebd_release, "/ebd_rel", ebd_release, ".txt")
orig_sampling <- paste0(ebd_data_folder, "ebd_sampling_rel", ebd_release, "/ebd_sampling_rel", ebd_release, ".txt")


# intentionally keep incomplete checklists and all protocols

f_ebd_all <- paste0(data_save, "ebd_", data_tag, "_all.txt")
f_sampling_all <- paste0(data_save, "ebd_sampling_", data_tag, "_all.txt")

# process the ebird data keeping only certain species, regions, and dates
if (!file.exists(f_ebd_all)) {
  ebd_filtered <- auk_ebd(orig_ebd, 
                          file_sampling = orig_sampling) %>% 
    auk_species(c("Wood Thrush", "Northern Bobwhite", 
                  "White Ibis")) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-05-01", "*-06-30")) %>% 
    auk_filter(file = f_ebd_all, file_sampling = f_sampling_all)
}

# zero fill
# intentional bad practice: zero-filling incomplete checklists
ebd_zf_all <- auk_zerofill(f_ebd_all, f_sampling_all, 
                           collapse = TRUE, complete = FALSE)

# clean up variables
ebd_zf_all <- ebd_zf_all %>% 
  mutate(
    # use species code
    species_code = ebird_species(scientific_name, "code"),
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(as.character(observation_count)),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started))

# additional filtering - last 10 years of data
# additional effort filters intentionally not applied
ebd_zf_all <- filter(ebd_zf_all, 
    year(observation_date) > 2016,
    year(observation_date) < 2019)

# output
ebd_zf_all <- ebd_zf_all %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         species_code,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebd_zf_all, paste0(data_save, "ebd_", data_tag, "_zf_all.csv"), na = "")




# ####################################################################
# good processing practices
f_ebd <- paste0(data_save, "ebd_", data_tag, ".txt")
f_sampling <- paste0(data_save, "ebd_sampling_", data_tag, ".txt")

# additional filtering
ebd_zf <- ebd_zf_all %>% 
  filter(
    # complete checklists
    all_species_reported == TRUE,

    # only keep stationary or traveling protocols
    protocol_type %in% c("Stationary", "Traveling"),

    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    !is.na(time_observations_started),
    number_observers <= 10)

# output
ebd_zf <- ebd_zf %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         species_code,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebd_zf, paste0(data_save, "ebd_", data_tag, "_zf.csv"), na = "")


