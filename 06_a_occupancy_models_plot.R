# occupancy model

library(auk)
library(sf)
library(raster)
library(dggridR)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(PresenceAbsence)
library(fields)
library(viridis)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(lubridate)

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# read in functions
walk(list.files("R", full.names = TRUE), source)


set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")
# setup spatial sampling regime
sample_spacing <- 5 # approx distance in km used for spatial subsampling cells

date <- Sys.Date()
# date <- "2020-02-20"
run_name <- paste0("occu_", sp_code, "_", date)


# --------------------------------------------------------------------
# set paths
figure_folder <- paste0("figures/occupancy/", run_name, "/")
output_folder <- paste0("output/occupancy/", run_name, "/")

dir.create(figure_folder, recursive = TRUE)
dir.create(output_folder, recursive = TRUE)


# load data ----
data_folder <- "data_proc/"
data_tag <- "mayjune_201718_bcr27"

# ebird data
ebird <- read_csv(paste0(data_folder, "data_all_4_models_", data_tag, ".csv"), na = "") 
species_count <- ebird[,which(colnames(ebird)==sp_code)] %>%
                  as.matrix() %>% as.vector() %>% as.numeric()
species_binary <- ifelse(is.na(species_count), 1, ifelse(species_count==0, 0, 1))
ebird$species_observed <- species_binary
ebird <- ebird %>%
        select(checklist_id, sampling_event_identifier, observer_id, species_observed, latitude, longitude,
                protocol_type, all_species_reported, observation_date, time_observations_started,
                duration_minutes, effort_distance_km, number_observers, type)

# modis covariates
habitat <- read_csv(paste0(data_folder, "modis_pland_checklists_", data_tag, ".csv"), 
                    col_types = cols(
                      .default = col_double(),
                      checklist_id = col_character()))
pred_surface <- read_csv(paste0(data_folder, "modis_pland_prediction-surface.csv"), 
                         col_types = cols(
                           .default = col_double(),
                           id = col_integer(),
                           year = col_integer()))

# combine modis and ebird data
ebird_habitat <- inner_join(ebird, habitat, by = "checklist_id") %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
        mutate(day_of_year = yday(observation_date))
        mutate(locality_id = paste0(round(ebird_habitat$longitude, 6), "_", 
          round(ebird_habitat$latitude, 6)))


# map data ----

map_proj <- st_crs(102003)
# borders
f_gpkg <- paste0(data_folder, "gis-data.gpkg")
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


# define model set ----
mod_set <- tibble(run_name = c("incomplete", "all", "complete", 
                               "sss", "effort", "best_practice"),
                  incomplete =        c(1, 0, 0, 0, 0, 0),
                  complete =          c(0, 0, 1, 1, 1, 1),
                  spatial_subsample = c(0, 0, 0, 1, 1, 1),
                  effort_filter =     c(0, 0, 0, 0, 1, 1),
                  effort_covs =       c(0, 0, 0, 0, 0, 1)) %>% 
  mutate_if(is.numeric, as.logical) %>% 
  mutate(run_id = row_number() + 1) %>% 
  select(run_id, everything())


# define locality within some tolerance ----


# ####################################################################
# RUN THE MODELS

print("RUNNING MODELS")

set.seed(1)
mod_set$models <- pmap(mod_set, fit_model_occu, data = ebird_habitat,
                      spacing = sample_spacing)


print("PREDICT OCCUPANCY")

# prediction ----

predict_raster <- function(model, data, template) {
  pred <- predict(model, newdata = as.data.frame(data), type = "state") %>% 
    select(occ_prob = Predicted, occ_se = SE) %>% 
    bind_cols(pred_surface, .)
  
  # rasterize
  r_pred <- pred %>% 
    select(occ_prob, occ_se, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(template)) %>% 
    rasterize(template)
  r_pred[[-1]]
}

r <- raster("data/modis_5xagg.tif")
mod_set <- mod_set %>% 
  mutate(r_pred = map(models, predict_raster, data = pred_surface, 
                      template = r))

# occupancy probability
f_prob <- str_glue("{output_folder}/07_bad_2_occupancy-model_pred-occ_{sp_code}_", Sys.Date(), ".tif")
r_occ_prob <- map(mod_set$r_pred, ~ .[["occ_prob"]]) %>% 
  stack() %>% 
  writeRaster(f_prob, overwrite = TRUE) %>% 
  setNames(mod_set$run_name)

# occupancy se
f_se <- str_glue("{output_folder}/07_bad_2_occupancy-model_pred-se_{sp_code}_", Sys.Date(), ".tif")
r_occ_se <- map(mod_set$r_pred, ~ .[["occ_se"]]) %>% 
  stack() %>% 
  writeRaster(f_se, overwrite = TRUE) %>% 
  setNames(mod_set$run_name)


# maps ----

r_pred_proj <- projectRaster(r_occ_prob, crs = map_proj$proj4string,
                             method = "ngb")

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
mx <- 1
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("{figure_folder}/occupancy_preds_{sp_code}_{date}.png") %>% 
  png(width = 2400, height = 3000, res = 300)

par(mfrow = c(3, 2), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))
for (i in 2:7) {
  
  r_plot <- r_pred_proj[[i - 1]]
  name <- paste("Model", i)
  
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # probability of detection
  plot(r_plot, col = pal, breaks = brks, maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE) 
  
  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  
  title(main = name, line = -2, cex.main = 1.5, font.main = 1)
  box()
}

# legend
par(new = TRUE, mfrow = c(1, 1), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
lbl_brks <- seq(0, mx, by = 0.1)
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.03, 0.05),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.8),
           legend.args = list(text = NULL,
                              side = 3, col = "black",
                              cex = 1, line = 0))
dev.off()




# --------------------------------------------------------------------
# predict detectability

  # complete checklists only
  pred_data <- filter(ebird_habitat, as.logical(all_species_reported))
  # filter on effort covariates
  pred_data <- pred_data %>%
      mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%      
      filter(protocol_type %in% c("Stationary", "Traveling"),
             effort_distance_km <= 5,
             duration_minutes <= 5 * 60,
             number_observers <= 5) %>% 
      mutate(protocol_type = factor(protocol_type, 
                                    levels = c("Stationary" , 
                                               "Traveling")))

  # prepare for unmarked
  occ <- filter_repeat_visits(pred_data, min_obs = 2, max_obs = 10,
                              annual_closure = TRUE,
                              date_var = "observation_date",
                              site_vars = c("locality_id")) #, "observer_id"))

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
    # generate hexagonal grid
    dggs <- dgconstruct(spacing = sample_spacing)
    occ_wide_cell <- occ_wide %>% 
      mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
    # sample one record (set of repeated observations at a site) per grid cell
    occ_wide <- occ_wide_cell %>% 
      group_by(cell) %>% 
      sample_n(size = 1) %>% 
      ungroup() %>% 
      select(-cell)
  
  # make unmarked object from the spatially subsampled data
  occ_um <- formatWide(occ_wide, type = "unmarkedFrameOccu")
  

det3 <- predict(mod_set$models[[2]], newdata = occ_um, type = "det")
det6 <- predict(mod_set$models[[5]], newdata = occ_um, type = "det")
det7 <- predict(mod_set$models[[6]], newdata = occ_um, type = "det")


str_glue("{figure_folder}/occupancy-model_predictions_detectability_{sp_code}_{date}.png") %>% 
  png(width = 14, height = 15, res = 600, pointsize = 9, units = "cm")

br <- seq(0, 1, by=0.05)
par(mfrow=c(2,1), mar = c(2, 5, 1, 1), oma = c(5, 1, 1, 1))

h_counts <- hist(det6$Predicted, breaks = br, plot = FALSE)$counts
at_y <- pretty(c(0, max(h_counts)), 3)
hist(det6$Predicted, breaks = br, xlab = "", ylab = "Frequency", main = "", xaxt="n", yaxt="n")
axis(side = 1, at = c(0, 0.5, 1), pos = 0)
axis(side = 2, at = at_y, xpd = NA)
text(x = 0, y = max(at_y)*0.95, labels = "a", font = 2, pos = 4, cex = 0.98, xpd = NA)
text(x = 0.03, y = max(at_y)*0.95, labels = "Model 6", pos = 4, cex = 0.98, xpd = NA)

h_counts <- hist(det7$Predicted, breaks = br, plot = FALSE)$counts
at_y <- pretty(c(0, max(h_counts)), 3)
hist(det7$Predicted, breaks = br, xlab = "Estimated detectability", ylab = "Frequency", main = "", xaxt="n", xpd = NA, yaxt="n")
axis(side=1, at=c(0, 0.5, 1), pos = 0)
axis(side = 2, at = at_y, xpd = NA)
text(x = 0, y = max(at_y)*0.95, labels = "b", font = 2, pos = 4, cex = 0.98, xpd = NA)
text(x = 0.03, y = max(at_y)*0.95, labels = "Model 7", pos = 4, cex = 0.98, xpd = NA)

dev.off()


# --------------------------------------------------------------------
# predictability on the 2017 data and aggregate within sites. 

test_2017 <- ebird_habitat %>%
            filter(type == "test_2017") %>%
            as.data.frame()

test_2017$det <- predict(mod_set$models[[6]], newdata = test_2017, type = "det")$Predicted

# define 'sites'
test_2017$site <- paste0(round(test_2017$longitude, 6), "_", round(test_2017$latitude, 6))

# split by sites
test_2017_sp <- split(test_2017, as.factor(test_2017$site))
aggregate_det <- function(x) {1 - prod(1 - x$det)}
site_det <- sapply(test_2017_sp, FUN = aggregate_det)
site_obs <- sapply(test_2017_sp, FUN = function(x){max(x$species_observed)})

site_det_df <- data.frame(site = names(site_det), det = as.numeric(site_det), obs = as.numeric(site_obs))
site_high_det <- site_det_df$site[site_det_df$det>0.95]


# separate out the sites with high detectability (e.g. det>0.9)
test_2017_high_det <- test_2017 %>%
          select(starts_with("pland"), site) %>%
          distinct() %>%
          right_join(site_det_df, by="site") %>%
          filter(det>0.9)

# predict occupancy at each site
test_2017_high_det$occ <- predict(mod_set$models[[6]], newdata = test_2017_high_det, type = "state")$Predicted

par(mfrow=c(2,1))
br <- seq(0, 1, by = 0.05)
hist(test_2017_high_det$occ[test_2017_high_det$obs==0], breaks = br)
hist(test_2017_high_det$occ[test_2017_high_det$obs==1], breaks = br)

# calculate ppms
all_ppms <- map(mod_set$models, predict, newdata = test_2017_high_det, type = "state") %>%
            lapply(FUN=function(x){x$Predicted}) %>%
            map(calculate_ppms, obs = test_2017_high_det$obs) %>%
            bind_rows() %>%
            mutate(model = paste0("Model ", 2:7)) %>% 
            mutate(tss = sensitivity + specificity - 1) 

plot_data <- all_ppms %>%
            select(-threshold) %>%
            gather(metric, value, -model) %>%
            mutate(metric = factor(metric, 
                           levels = c("mse", "auc", "tss",
                                      "kappa", "sensitivity", "specificity"),
                           labels = c("Mean Squared Error (MSE)", "AUC", "TSS",
                                      "Kappa", "Sensitivity", "Specificity")))


# plot ppms
  g_ppm <- ggplot(plot_data) +
    aes(x = model, y = value) +
    geom_point() +
    facet_wrap(~ metric, nrow = 2, scales = "free_y") +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 12, hjust = 0),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "transparent"),
          axis.ticks.y = element_line(),
          panel.grid = element_blank())
  plotname <- str_glue("{figure_folder}/occu-model_2017_high_det_{sp_code}.png") %>%
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)




# number of observed positives

obs <- occ_wide[,2:11]
table(apply(obs, 1, max, na.rm = TRUE))  
 #   0    1 
 # 769 1318 


occ_pred <- predict(mod_set$models[[6]], newdata = occ_um, type = "state")


## find only the zero sites
occ_zero <- occ_wide %>%
            filter(pland_04==0, pland_05==0, pland_12==0, pland_13==0)
# make unmarked object from the spatially subsampled data
occ_zero_um <- formatWide(occ_zero, type = "unmarkedFrameOccu")

mod <- occu(~day_of_year + time_observations_started + duration_minutes + 
    effort_distance_km + number_observers + protocol_traveling ~ 1,
    data = occ_zero_um)
