# Run encounter rate models with different data processing
# Run 25 times with different random data samples
# Calculate validation metrics

library(auk)
library(sf)
library(raster)
library(dggridR)
library(ranger)
library(maxnet)
library(scam) 
library(PresenceAbsence)
library(verification)
library(edarf)
library(ggplot2)
library(ggthemes)
library(hexbin)
library(viridis)
library(fields)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(forcats)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# read in functions
walk(list.files("R", full.names = TRUE), source)

set.seed(1)

# set parameters for analysis
species <- "Wood Thrush" #species
sp_code <- ebird_species(species, "code")
sample_regime <- "together" # spatial subsample regime. see R/hex_sample.R for details 
sample_spacing <- 5 # approx distance in km used for spatial subsampling cells
calibrate <- TRUE # whether to calibrate the predictions for random forest
anchor_model <- 3 # model to use as comparison for validation metrics

date <- Sys.Date()
run_name <- paste0("RF_sample_", sample_regime, "_", sp_code, "_", date)


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
ebird_habitat <- inner_join(ebird, habitat, by = "checklist_id")


ebird_habitat <- ebird_habitat %>%
        mutate(week = lubridate::week(observation_date)) %>%
        mutate(type_week = paste(type, week, sep="_")) %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(day_of_year = yday(observation_date))

# ####################################################################
# fit the set of models ----


# define the params for each model combinations ----
# note that this matrix matches table 1 in the paper 

mod_set_master <- tibble(run_name = c("maxent", "incomplete", 
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


nsim <- 25

for(i in 1:nsim){

  print("")
  print("=======================================")
  print(paste("RUNNING FOR SIM", i))
  set.seed(i)

  # ####################################################################
  # PREPARE THE DATA

  print("PREPARING DATA")

  # reduce to 75% of all data (in all subsets)
  ebird_habitat_sample <- ebird_habitat %>%
      sample_frac(0.75)

  # subsample positives and negatives together 
  # for train data and test 2017
  ebird_subsamp <- ebird_habitat_sample %>%
    hex_sample(spacing = sample_spacing, regime = sample_regime, byvar = "type_week") %>%
    select(checklist_id, sampling_event_identifier, type) %>%
    mutate(selected = 1) %>%
    right_join(ebird_habitat_sample) %>%
    mutate(selected = ifelse(is.na(selected), 0, selected))

  # --------------------------------------------------------------------
  # create separate datasets
  
  ebird_test_bbs <- ebird_subsamp %>%
    filter(type == "test_bbs") %>%
    left_join(read_csv(paste0(data_folder, "bbs_route_id.txt")), by="checklist_id")

  ebird_test_2017 <- ebird_subsamp %>%
    filter(type == "test_2017", selected == 1)

  # include all obs and the model function selects only subsampled on request
  model_data <- ebird_subsamp %>%
    filter(type == "train")

  # --------------------------------------------------------------------
  # create balanced test data with 2017 data 
  # subsample positives and negatives separately

  subsamp_pn <- ebird_habitat_sample %>%
    filter(type == "test_2017") %>%
    hex_sample(spacing = sample_spacing, regime = "separate", byvar = "week")

  npos <- sum(subsamp_pn$species_observed)

  ebird_test_2017_bal <- subsamp_pn %>% 
    select(species_observed, checklist_id) %>%
    group_by(species_observed) %>%
    sample_n(npos) %>%
    ungroup() %>%
    left_join(ebird_subsamp)


  # ####################################################################
  # RUN THE MODELS
  
  print("RUNNING MODELS")
  mod_set <- mod_set_master

  mod_set$models <- pmap(mod_set, fit_model_enc, data = model_data,
                       calibrate = calibrate, calibrate_plot = FALSE, 
                       subsample_seed = i)

  # ####################################################################
  # VALIDATE THE MODELS
  print("VALIDATING MODELS")

  ppm_bbs_stop <- mutate(mod_set, ppms = map(models, validate, data = ebird_test_bbs, bbs_combine = FALSE)) %>%
    select(-models) %>% unnest(cols = c(ppms)) %>%
    mutate(val_type = "bbs")

  ppm_bbs_route <- mutate(mod_set, ppms = map(models, validate, data = ebird_test_bbs, bbs_combine = TRUE)) %>%
    select(-models) %>% unnest(cols = c(ppms)) %>%
    mutate(val_type = "bbs_route")

  ppm_2017 <- mutate(mod_set, ppms = map(models, validate, data = ebird_test_2017)) %>%
    select(-models) %>% unnest(cols = c(ppms)) %>%
    mutate(val_type = "2017")

  ppm_2017_bal <- mutate(mod_set, ppms = map(models, validate, data = ebird_test_2017_bal)) %>%
    select(-models) %>% unnest(cols = c(ppms)) %>%
    mutate(val_type = "2017_bal")

  ppms <- rbind(ppm_bbs_stop, ppm_bbs_route, ppm_2017, ppm_2017_bal) %>%
    mutate(tss = sensitivity + specificity - 1) %>%
    mutate(sim_id = i)

  if(i==1) all_ppms <- ppms
  if(i>1) all_ppms <- rbind(all_ppms, ppms)

}

str_glue("{output_folder}/enc_RF_ppms_{sp_code}_{date}.csv") %>% 
  write_csv(all_ppms, .)

# all_ppms <- str_glue("{output_folder}/enc_RF_ppms_{sp_code}_{date}.csv") %>% read_csv()


# 2017 validation plot

# plot comparing ppms
ppm_plot <- all_ppms %>% 
  select_if(~ !is.logical(.)) %>% 
  gather("metric", "value", -run_id, -run_name, -sim_id, -val_type) %>%
  filter(metric != "threshold", metric != "n_checklists", metric != "n_pos") %>% 
  arrange(val_type, sim_id, run_id) %>% 
  mutate(metric_label = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
  mutate(metric_label = factor(metric_label, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS"))) %>%
  mutate(run = paste("Model", run_id),
         run = as_factor(run),
         start = if_else(metric == "AUC", 0.5, 0)) %>%
  filter(!is.na(metric))

for(i in 1:4){

  val_type_plot <- c("bbs", "bbs_route", "2017", "2017_bal")[i]
  plot_data <- filter(ppm_plot, val_type==val_type_plot)

  g_ppm <- ggplot(plot_data) +
    aes(x = run, y = value) +
    geom_boxplot(coef=5) +
    facet_wrap(~ metric_label, nrow = 2, scales = "free_y") +
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
  plotname <- str_glue("{figure_folder}/enc_RF_ppms_{val_type_plot}_{sp_code}.png") %>%
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)

}

for(i in 1:4){

  val_type_plot <- c("bbs", "bbs_route", "2017", "2017_bal")[i]
  plot_data <- filter(ppm_plot, val_type==val_type_plot)

  g_ppm <- ggplot(plot_data) +
    aes(x = run, y = value) +
    geom_boxplot(coef=5) +
    facet_wrap(~ metric_label, nrow = 2, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 12, hjust = 0),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "transparent"),
          axis.ticks.y = element_line(),
          panel.grid = element_blank())
  plotname <- str_glue("{figure_folder}/enc_RF_ppms_unscaled_{val_type_plot}_{sp_code}.png") %>%
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)

}



# --------------------------------------------------------------------
# plot differences from anchor model

ppms_gather <- all_ppms %>%
  select_if(~ !is.logical(.)) %>% 
  select(-threshold) %>%
  gather("metric", "value", -run_id, -run_name, -sim_id, -val_type)

diff_plot <- ppms_gather %>%
  filter(run_id == anchor_model) %>%
  select(metric, sim_id, value, val_type) %>%
  rename(value_anchor = value) %>%
  right_join(ppms_gather) %>%
  mutate(diff = value - value_anchor) %>%
  arrange(val_type, sim_id, run_id) %>% 
  mutate(run = paste("Model", run_id),
         run = as_factor(run),
         start = if_else(metric == "AUC", 0.5, 0)) %>%
  filter(!is.na(metric))



add_grey <- TRUE

for(i in 1:4){

  val_type_plot <- c("bbs", "bbs_route", "2017", "2017_bal")[i]
  plot_data <- diff_plot %>%
      filter(val_type==val_type_plot) %>%
      filter(! metric %in% c("n_checklists", "n_pos")) %>%
      mutate(metric_short = metric) %>%
      mutate(metric = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
      mutate(metric = factor(metric, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS")))

  maxy <- plot_data %>% 
            select(metric, diff) %>% 
            group_by(metric) %>%
            summarise(max_abs = max(abs(diff))) %>%
            ungroup()

  nmod <- length(table(plot_data$run))
  nmetric <- 6

  str_glue("{figure_folder}/enc_RF_ppms_DIFF_{val_type_plot}_anchor{anchor_model}_{sp_code}_grey{add_grey}.png") %>%
    png(width = 21, height = 17, units="cm", res = 600)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))

        # performance metrics
        for (j in 1:nmetric) {
          m <- levels(plot_data$metric)[j]
          maxyy <- maxy$max_abs[j]
          xnames <- rep("", nmod)
          if(j>3) xnames <- paste("Model", 1:nmod)
          boxplot(as.formula("diff ~ run"), 
                  data = plot_data[plot_data$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  boxcol = "white", 
                  col = "white",
                  xlab = "", ylab = "",
                  ylim=c(-1*maxyy, maxyy), las = 2, names = xnames)
          if(!add_grey) abline(h=0, lwd=2, col="grey70")
          if(add_grey) {
            ymin <- -2
            ymax <- 0
            if(j==1) {ymin <- 0; ymax <- 1 }
            polygon(x = c(-1, 10, 10, -1, -1), y = c(ymin, ymin, ymax, ymax, ymin), col="grey78", border = alpha("white", 0))
          }
          par(new=TRUE)
          boxplot(as.formula("diff ~ run"), 
                  data = plot_data[plot_data$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  xlab = "", ylim = c(-1*maxyy, maxyy), names = rep("", nmod), 
                  col = alpha("white", 0.4), 
                  xaxt="n", yaxt="n", 
                  ylab = bquote(Delta~.(levels(plot_data$metric)[j])))
          text(x = 0.5, y = maxyy*0.95, 
               labels = LETTERS[j], 
               font = 2, pos=4)
        }

    dev.off()

}


