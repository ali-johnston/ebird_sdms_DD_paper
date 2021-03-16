# Run encounter rate models with different data processing
# Run 1 time with all data for each option
# Predict to entire region
# Plot maps

library(auk)
library(sf)
library(raster)
library(dggridR)
library(ranger)
library(maxnet)
library(scam) 
library(PresenceAbsence)
library(verification)
#devtools::install_github("zmjones/edarf", subdir = "pkg")
library(edarf)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(hexbin)
library(viridis)
library(fields)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(forcats)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# read in functions
walk(list.files("R", full.names = TRUE), source)

#all_species <- c("Swainson's Warbler", "Red-cockaded Woodpecker",
#                  "Henslow's Sparrow",  # crashed out on low prevalence
# "Eastern Whip-poor-will", # crashed out on optimal time of day.
all_species <- c("Wood Thrush", "Chuck-will's-widow")


for(i_species in 1:length(all_species)){

    set.seed(1)
    # set species for analysis
    species <- all_species[i_species]
    sp_code <- ebird_species(species, "code")

    # setup spatial sampling regime
    sample_regime <- "together"     # spatial subsample regime. see R/hex_sample.R for details 
    sample_spacing <- 5             # approx distance in km used for spatial subsampling cells
    calibrate <- TRUE               # whether to calibrate the predictions for random forest

    date <- Sys.Date()
    run_name <- paste0("RF_sample_", sample_regime, "_", sp_code, "_", date)

    # --------------------------------------------------------------------
    # set paths

    # set folder paths for results ----
    figure_folder <- paste0("figures/encounter/", run_name, "/")
    output_folder <- paste0("output/encounter/", run_name, "/")

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
            select(checklist_id, sampling_event_identifier, species_observed, latitude, longitude,
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
            mutate(week = lubridate::week(observation_date)) %>%
            mutate(type_week = paste(type, week, sep="_")) %>%
            mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
            mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
            mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
            mutate(number_observers = as.numeric(as.character(number_observers))) %>%
            mutate(day_of_year = yday(observation_date))


    # map data ----

    map_proj <- st_crs(5070)
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


    # create the model set definitions

    mod_set <- tibble(run_name = c("maxent", "incomplete", 
                                          "all", "complete", 
                                          "sss", "effort", "covs"),
                      maxnet =            c(1, 0, 0, 0, 0, 0, 0),
                      incomplete =        c(0, 1, 0, 0, 0, 0, 0),
                      complete =          c(0, 0, 0, 1, 1, 1, 1),
                      spatial_subsample = c(0, 0, 0, 0, 1, 1, 1),
                      effort_filter =     c(0, 0, 0, 0, 0, 1, 1),
                      effort_covs =       c(0, 0, 0, 0, 0, 0, 1)) %>% 
      mutate_if(is.numeric, as.logical) %>% 
      mutate(run_id = row_number()) %>% 
      select(run_id, everything())



    # spatial subsampling for train and test_2017 ----
    set.seed(1)

    # ####################################################################
    # PREPARE THE DATA

    print("PREPARING DATA")

    # reduce to 75% of all data
    ebird_habitat_sample <- ebird_habitat %>%
        sample_frac(0.75)

      # subsample positives and negatives together for train data
      ebird_subsamp <- ebird_habitat_sample %>%
        hex_sample(spacing = sample_spacing, regime = sample_regime, byvar = "type_week") %>%
        select(checklist_id, sampling_event_identifier, type) %>%
        mutate(selected = 1) %>%
        right_join(ebird_habitat_sample) %>%
        mutate(selected = ifelse(is.na(selected), 0, selected))

      # include all obs and the model function selects only subsampled on request
      model_data <- ebird_subsamp %>%
        filter(type == "train")

    # fit model set ----
    mod_set$models <- pmap(mod_set, fit_model_enc, data = model_data,
                           calibrate = calibrate, calibrate_plot = TRUE)

    # amount of data in each run
    run_counts <- mod_set %>% 
      mutate(n_checklists = unlist(map(models, "n_checklists")),
             n_sightings = unlist(map(models, "n_sightings"))) %>% 
      select(-models)
    str_glue("{output_folder}/RF-model_counts_{sp_code}.csv") %>% 
      write_csv(run_counts, .)


    # ####################################################################
    # PREDICT ACROSS MAP

    predict_raster <- function(model, data, template) {
      # add effort covariates to prediction surface
      data <- data %>% 
        mutate(observation_date = ymd("2018-06-15"),
               day_of_year = yday(observation_date),
               time_observations_started = model$t_max_det[1],
               duration_minutes = 60,
               effort_distance_km = 1,
               number_observers = 1, 
               checklist_calibration_index = 2,
               protocol_type = "Traveling", 
               protocol_traveling = 1)
      
      # predict
      pred <- predict_model_enc(model, data)
      pred_df <- bind_cols(data, prob = pred)
      
      # rasterize
      print("rasterize")
      r_pred <- pred_df %>% 
        select(prob, latitude, longitude) %>% 
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
        st_transform(crs = projection(template)) %>% 
        rasterize(template)
      r_pred[[-1]]
    }

    r <- raster(paste0(data_folder, "/modis_5xagg.tif"))
    r_pred <- map(mod_set$models, predict_raster, 
                  data = pred_surface, template = r) %>% 
      stack()
    r_pred <- str_glue("{output_folder}/enc_RF_pred_{sp_code}_{date}.tif") %>% 
      writeRaster(r_pred, ., overwrite = TRUE) %>% 
      setNames(mod_set$run_name)


    # density plot ----

    # prepare data
    pred_compare <- rasterToPoints(r_pred) %>% 
      as_tibble() %>% 
      drop_na() %>% 
      select(-x, -y) %>% 
      gather("run_name", "mod_type", -covs) %>% 
      inner_join(mod_set %>% select(run_id, run_name), by = "run_name") %>% 
      select(run_id, mod_type, covs) %>%
      arrange(run_id) %>% 
      mutate(run = paste("Model", run_id)) %>%
      filter(run!="Model 7") %>%
      mutate(run = as_factor(run))

    # Function to generate correlation coefficient for the charts
    corr_eqn <- function(x,y, digits = 2) {
      corr_coef <- round(cor(x, y), digits = digits)
      corr_coef <- expression(paste(italic(r)," = ", corr_coef))
      return(corr_coef)
    }

    # plot
    g_density <- ggplot(pred_compare) + 
      aes(x = covs, y = mod_type) +
      geom_hex(aes(fill = stat(count))) + 
      stat_cor(method = "pearson") + 
      scale_x_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) + 
      scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) +
      coord_equal() +
      scale_fill_viridis_c(trans = "log10", labels = scales::comma) + 
      facet_wrap(~ run, nrow = 2) + 
      labs(x = "Model 7 predictions", y = "Model 1-6 predictions") + 
      guides(fill = guide_colorbar(title = "# predictions", 
                                   title.position = "left",
                                   barwidth = 0.5, barheight = 12)) +
      theme_few() +
      theme(legend.position = "right",
            legend.title = element_text(angle = 90, hjust = 0.5))
    str_glue("{figure_folder}/enc_RF-model_scatter_density_{sp_code}.png")  %>% 
      ggsave(g_density, width = 20, height = 12, units = "cm", dpi = 300)


    # maps ----

    r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

    # breaks and palette
    plasma_rev <- rev(plasma(25, end = 0.9))
    gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
    pal <- c(gray_int(4)[2], plasma_rev)
    # mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
    mx <- 1
    brks <- seq(0, mx, length.out = length(pal) + 1)

    str_glue("{figure_folder}/RF-model_predictions_{sp_code}_{date}.png") %>% 
      png(width = 2400, height = 3600, res = 300)

      par(mfrow = c(4, 2), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))

      plot(0, 0, col="white", bty="n", xlab="", ylab="", xaxt="n", yaxt="n")

      for (i in seq.int(nlayers(r_pred_proj))) {
        r_plot <- r_pred_proj[[i]]
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
                 smallplot = c(0.25, 0.75, 0.035, 0.055),
                 horizontal = TRUE,
                 axis.args = list(at = lbl_brks, labels = lbl_brks,
                                  fg = "black", col.axis = "black",
                                  cex.axis = 0.75, lwd.ticks = 0.5,
                                  padj = -1.8),
                 legend.args = list(text = NULL,
                                    side = 3, col = "black",
                                    cex = 1, line = 0))
    dev.off()



    # ####################################################################
    # COMPARE ENCOUNTER IN MODEL 2, MODEL 3, MODEL 7

    enc1 <- r_pred$maxent %>% getValues()
    enc2 <- r_pred$incomplete %>% getValues()
    enc3 <- r_pred$all %>% getValues()
    enc4 <- r_pred$complete %>% getValues()
    enc5 <- r_pred$sss %>% getValues()
    enc6 <- r_pred$effort %>% getValues()
    enc7 <- r_pred$covs %>% getValues()

    plot_hist_pretty <- function(x, br = seq(0, 1, by=0.05), 
      xlabel="", ylabel="", letter = 1, name = "", 
      axis_x = FALSE, ...) {
      h_counts <- hist(x, breaks = br, plot = FALSE)$counts
      at_y <- pretty(c(0, max(h_counts)), 3)
      hist(x, breaks = br, main = "", xaxt="n", yaxt="n", xlab=xlabel, ylab=ylabel, xpd = NA)
      axis(side = 2, at = at_y, xpd = NA)
      axis(side = 1, at = c(0, 0.5, 1), pos = 0, labels = rep("", 3))
      if(axis_x) axis(side=1, at=c(0, 0.5, 1), pos = 0)
      text(x = 0, y = max(at_y)*0.9, labels = LETTERS[letter], font = 2, pos = 4, cex = 0.98, xpd = NA)
      text(x = 0.07, y = max(at_y)*0.9, labels = name, pos = 4, cex = 0.98, xpd = NA)
    }


    str_glue("{figure_folder}/encounter-model_predictions_encounter_{sp_code}_{date}.png") %>% 
      png(width = 13, height = 17, res = 600, pointsize = 9, units = "cm")

    par(mfrow=c(4,2), mar = c(2, 5, 2, 1), oma = c(3, 1, 4, 1))

    plot(0, 0, col="white", bty="n", axes = FALSE, xlab = "", ylab = "")
    plot_hist_pretty(enc1, letter = 1, name = "Model 1")

    plot_hist_pretty(enc2, letter = 2, name = "Model 2")
    plot_hist_pretty(enc3, letter = 3, name = "Model 3")

    plot_hist_pretty(enc4, letter = 4, name = "Model 4", ylabel = "Frequency")
    plot_hist_pretty(enc5, letter = 5, name = "Model 5")

    plot_hist_pretty(enc6, letter = 6, name = "Model 6", axis_x = TRUE, xlabel = "Estimated encounter rate")
    plot_hist_pretty(enc7, letter = 7, name = "Model 7", axis_x = TRUE, xlabel = "Estimated encounter rate")

    dev.off()

} # close i_species

