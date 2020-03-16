


# t <- 1

# complete <- bp_runs$complete[t]
# spatial_subsample <- bp_runs$spatial_subsample[t]
# effort_filter <- bp_runs$effort_filter[t]
# effort_covs <- bp_runs$effort_covs[t]
# spacing <- sample_spacing
# data <- ebird_habitat


fit_model_occu <- function(complete, incomplete, 
                         spatial_subsample, 
                         effort_filter, effort_covs,
                         data, spacing, ...) {

  print("------- Running occupancy model -------")
  print(paste("raw data n =", nrow(data)))

  # incomplete checklists only
  if (incomplete) {
    data <- filter(data, !as.logical(all_species_reported))
  }

  # complete checklists only
  if (complete) {
    data <- filter(data, as.logical(all_species_reported))
  }

  # filter on effort covariates
  if (effort_filter){
    data <- data %>%
      mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%      
      filter(protocol_type %in% c("Stationary", "Traveling"),
             effort_distance_km <= 5,
             duration_minutes <= 5 * 60,
             number_observers <= 5) %>% 
      mutate(protocol_type = factor(protocol_type, 
                                    levels = c("Stationary" , 
                                               "Traveling")))
  }
  
  # prepare for unmarked
  # require a minimum of 2 observations per site 
  occ <- filter_repeat_visits(data, min_obs = 2, max_obs = 10,
                              annual_closure = TRUE,
                              date_var = "observation_date",
                              site_vars = c("locality_id")) 
  
  # convert to a wide format for use with unmarked
  # include all detection covariates
  # choose habitat covariates based on prior knowledge of species
  # for wood thrush:
  #   - deciduous (pland_04) and mixed (pland_05) forest; known nesting habitat
  #   - cropland (pland_12) and urban (pland_13); avoided
  # for northern bobwhite:
  #   - savanna (pland_09) and woody savanna (pland_08); preferred
  #   - deciduous forest (pland_04) and urban (pland_13); avoided
  # for white ibis:
  #   - fresh water (pland_00) and wetlands (pland_11); preferred
  #   - savanna (pland_09) and woody savanna (pland_08); avoided
  if (sp_code == "woothr") {
    plands <- paste0("pland_", c("04", "05", "12", "13"))
  } else if (sp_code == "norbob") {
    plands <- paste0("pland_", c("04", "08", "09", "13"))
  } else if (sp_code == "whiibi") {
    plands <- paste0("pland_", c("00", "08", "09", "11"))
  } else {
    stop("species code not valid.")
  }
  o_covs <- c("day_of_year",
              "time_observations_started", 
              "duration_minutes", 
              "effort_distance_km", 
              "number_observers", 
              "protocol_traveling",
              "checklist_calibration_index") %>% 
    intersect(names(occ))
  if(!effort_covs) o_covs <- "day_of_year"
  occ_wide <- occ %>%
      mutate(latitude = round(latitude, 6)) %>%
      mutate(longitude = round(longitude, 6)) %>%
      format_unmarked_occu(site_id = "site", 
                                   response = "species_observed",
                                   site_covs = c("latitude", "longitude",
                                                # habitat covariates
                                                 plands),
                                   obs_covs = o_covs)
  

  # spatial subsampling
  if (spatial_subsample) {
    # generate hexagonal grid
    dggs <- dgconstruct(spacing = spacing)
    occ_wide_cell <- occ_wide %>% 
      mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
    # sample one record (set of repeated observations at a site) per grid cell
    occ_wide <- occ_wide_cell %>% 
      group_by(cell) %>% 
      sample_n(size = 1) %>% 
      ungroup() %>% 
      select(-cell)
  } 
  
  # make unmarked object from the spatially subsampled data
  occ_um <- formatWide(occ_wide, type = "unmarkedFrameOccu")
  print(paste("number of sites = ", nrow(occ_wide)))
  
  # fit an occupancy model
  # four habitat types and seven detection covariates
  # specify detection model
    model_formula <- o_covs %>% 
      paste(collapse = " + ") %>% 
      paste("~", ., "~ .") %>% 
      as.formula()
  
  # add in species specific occupancy formula
  occ_formula <- keep(names(occ_wide), ~ str_detect(., "^pland")) %>% 
    paste(collapse = " + ") %>% 
    paste("~ ", .) %>% 
    as.formula()
  model_formula <- update(model_formula, occ_formula)
  
  # fit model
  occu(model_formula, data = occ_um)
}
