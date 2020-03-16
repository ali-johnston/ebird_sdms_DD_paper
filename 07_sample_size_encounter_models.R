# Run encounter rate models 3 and 7 with different sample sizes
# Run 25 times with different random data samples
# Calculate validation metrics

library(auk)
library(sf)
library(raster)
library(dggridR)
library(ranger)
library(scam) 
library(PresenceAbsence)
library(verification)
library(edarf)
library(fields)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)

walk(list.files("R", full.names = TRUE), source)

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")
# setup spatial sampling regime
sample_regime <- "together" # spatial subsample regime. see R/hex_sample.R for details 
sample_spacing <- 5 # approx distance in km used for spatial subsampling cells
calibrate <- TRUE # whether to calibrate the predictions for random forest
anchor_model <- 3 # model to use as comparison for validation metrics


date <- Sys.Date()
# date <- "2020-02-20"
run_name <- paste0("SS_RF_", sample_regime, "_", sp_code, "_", date)


# --------------------------------------------------------------------
# set paths
figure_folder <- paste0("figures/sample_size/", run_name, "/")
output_folder <- paste0("output/sample_size/", run_name, "/")

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
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(day_of_year = yday(observation_date))


# setup bad practice combinations ----
prop_runs <- seq(0.1, 0.9, by=0.2)

mod_set_master <- cross_df(list(good = c(TRUE, FALSE), 
                         prop_data = prop_runs)) %>%
            mutate(maxnet = FALSE) %>%
            mutate(incomplete = FALSE) %>%
            mutate(complete = ifelse(good, TRUE, FALSE)) %>%
            mutate(spatial_subsample = ifelse(good, TRUE, FALSE)) %>%
            mutate(effort_filter = ifelse(good, TRUE, FALSE)) %>%
            mutate(effort_covs = ifelse(good, TRUE, FALSE)) %>%
            mutate(run_id = row_number()) %>% 
            select(run_id, everything())



nsim <- 25

results_file_name <- str_glue("{output_folder}/enc_RF_sample_size_{nsim}_{sp_code}_{date}.csv")

if(!file.exists(results_file_name)){

  for(i in 1:nsim){

    ss <- i
    set.seed(ss)

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
                        calibrate = calibrate, subsample_seed = ss)

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

  results_file_name %>% 
    write_csv(all_ppms, .)

}

all_ppms <- results_file_name %>% 
  read_csv()

# 2017 validation plot

# plot comparing ppms
ppm_plot <- all_ppms %>% 
  mutate(model_name = ifelse(good, "Model 7", "Model 3")) %>%
  select_if(~ !is.logical(.)) %>% 
  gather("metric", "value", -run_id, -sim_id, -val_type, -prop_data, -model_name, -n_checklists, -n_pos) %>%
  filter(metric != "threshold", metric != "n_checklists", metric != "n_pos") %>% 
  arrange(val_type, sim_id, run_id) %>% 
  mutate(metric_short = metric) %>%
  mutate(metric = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
  mutate(metric = factor(metric, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS"))) %>%
  mutate(run = paste("Model", run_id),
         run = as.factor(run),
         start = if_else(metric == "AUC", 0.5, 0)) %>%
  filter(!is.na(metric)) %>%
  mutate(prop_data_cat = factor(prop_data, levels=sort(as.numeric(names(table(all_ppms$prop_data))))))


for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    ppm_best <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 6")
    ppm_bad <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 2")

    str_glue("{figure_folder}/enc_RF_sample_size_{plot_val_type}_{sp_code}.png") %>% 
      png(width = 2200, height = 1500, res = 300)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))
        cols <- c(bad = "grey65", best = "black")

        ppm_min <- c(0, 0.50, 0, 0, 0, 0)
        ppm_max <- c(0.2, 1, 1, 1, 1, 1)

        # performance metrics
        for (i in 1:6) {
          m <- levels(ppm_best$metric)[i]
          brks <- pretty(c(ppm_min[i], ppm_max[i]), n = 2)
          boxplot(as.formula("value ~ prop_data"), 
                  data = ppm_bad[ppm_bad$metric==m,],
                  range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
                  boxcol = "transparent", 
                  border = cols["bad"], 
                  xlim = c(0.5, 11), ylim = range(brks), 
                  xaxt = "n", yaxt = "n", 
                  xlab = "", ylab = levels(ppm_best$metric)[i])
          par(new = TRUE)
          boxplot(as.formula("value ~ prop_data"), 
                  data = ppm_best[ppm_best$metric==m,],
                  range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
                  boxcol = "transparent", 
                  border = cols["best"], 
                  xlim = c(0.15, 10.65), ylim = range(brks), 
                  xaxt = "n", yaxt = "n",
                  xlab = "", ylab = "")
          par(new = TRUE)
          axis(side = 2, at = brks)
          if (i > 3) {
            axis(side = 1, at = c(0, 5, 10), labels = c(0, 0.5, 1))
          }
          text(x = 0.5, y = max(brks) - 0.02 * diff(range(brks)), 
               labels = LETTERS[i], 
               font = 2)

          # add legend
          if(i==1) legend(x = 2, y = 0.15, 
                 lwd = 2, col = cols, 
                 legend = c("Model 2", "Model 6"), cex = 0.8)

        }

    dev.off()

}


metrics <- levels(ppm_plot$metric)

for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    ppm_best <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 7")
    ppm_bad <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 3")

    str_glue("{figure_folder}/enc_RF_sample_size_{plot_val_type}_{sp_code}.png") %>% 
      png(width = 2200, height = 1500, res = 300)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))
        cols <- c(bad = "grey65", best = "black")

        # ppm_min <- c(0, 0.50, 0, 0, 0, 0)
        # ppm_max <- c(0.2, 1, 1, 1, 1, 1)
        # ylimits <- data.frame(metric_short = metrics, ymin = ppm_min, ymax = ppm_max)

        ylimits <- ppm_best %>% 
              rbind(ppm_bad) %>%
              select(metric, value) %>%
              # filter(value < Inf) %>%
              # filter(value > -Inf) %>%
              group_by(metric) %>%
              summarise(pmin = min(value), pmax = max(value)) %>%
              ungroup() %>%
              mutate(ymin = pmin - (pmax - pmin)*0.05) %>%
              mutate(ymax = pmax + (pmax - pmin)*0.05)


        # performance metrics
        for (j in 1:6) {
          m <- metrics[j]
          m_long <- metrics[j]
          brks_y <- pretty(c(ylimits$ymin[j], ylimits$ymax[j]), n = 2)
          b <- boxplot(as.formula("value ~ prop_data"), 
                  data = ppm_bad[ppm_bad$metric==m,], 
                  plot = FALSE,
                  range = 0)
          xlabel <- ""
          if(j==5) xlabel <- "Proportion of dataset"
          ylabel <- m_long

          plot(0, 0, col="white",
            xlim=c(0, 1), ylim=c(brks_y[1], brks_y[length(brks_y)]), 
            xaxt="n", yaxt="n",
            xlab = xlabel, ylab = ylabel, xpd = NA)
          axis(side = 2, at=brks_y)
          if(j>3) axis(side = 1, at = c(0, 0.5, 1))

          x_real <- as.numeric(b$names)
          x_tr <- x_real
          x_wd <- 0.01
          x_adj <- 0.02
          segments(x0 = x_tr - x_adj - x_wd, x1 = x_tr - x_adj + x_wd, y0 = b$stats[3,], y1=b$stats[3,], col=cols[1], lwd=3)
          segments(x0 = x_tr - x_adj, x1 = x_tr - x_adj, y0 = b$stats[1,], y1=b$stats[2,], col=cols[1])
          segments(x0 = x_tr - x_adj, x1 = x_tr - x_adj, y0 = b$stats[4,], y1=b$stats[5,], col=cols[1])

          b_gd <- boxplot(as.formula("value ~ prop_data"), 
                  data = ppm_best[ppm_best$metric==m,], 
                  plot = FALSE,
                  range = 0)

          segments(x0 = x_tr + x_adj - x_wd, x1 = x_tr + x_adj + x_wd, y0 = b_gd$stats[3,], y1=b_gd$stats[3,], col=cols[2], lwd=3)
          segments(x0 = x_tr + x_adj, x1 = x_tr + x_adj, y0 = b_gd$stats[1,], y1=b_gd$stats[2,], col=cols[2])
          segments(x0 = x_tr + x_adj, x1 = x_tr + x_adj, y0 = b_gd$stats[4,], y1=b_gd$stats[5,], col=cols[2])

          text(x = 0, y = max(brks_y) - 0.02 * diff(range(brks_y)), 
               labels = LETTERS[j], 
               font = 2, pos = 4)

          # add legend
          if(j==1) legend(x = 2, y = -1.5, 
                 lwd = 2, col = cols, 
                 legend = c("Model 3", "Model 7"), cex = 0.8)

        }

    dev.off()

}





for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    ppm_best <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 6")
    ppm_bad <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 2")

    str_glue("{figure_folder}/enc_RF_sample_size_no_checklists_{plot_val_type}_{sp_code}.png") %>% 
      png(width = 2200, height = 1500, res = 300)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))
        cols <- c(bad = "grey65", best = "black")

        ppm_min <- c(0, 0.50, 0, 0, 0, 0)
        ppm_max <- c(0.2, 1, 1, 1, 1, 1)

        xlim <- range(c(ppm_best$n_checklists, ppm_bad$n_checklists))
        xlim <- signif(xlim, 2)
        xlim_log10 <- log10(xlim)

        xlim <- c(max(c(0, xlim[1] - (xlim[2]-xlim[1])*0.05)), xlim[2] + (xlim[2]-xlim[1])*0.05)
        xlim_log10 <- c(xlim_log10[1] - (xlim_log10[2]-xlim_log10[1])*0.05, xlim_log10[2] + (xlim_log10[2]-xlim_log10[1])*0.05)

        brks_x <- pretty(c(xlim[1], xlim[2]), n = 3)


        ylimits <- ppm_best %>% 
              rbind(ppm_bad) %>%
              select(metric, value) %>%
              # filter(value < Inf) %>%
              # filter(value > -Inf) %>%
              group_by(metric) %>%
              summarise(pmin = min(value), pmax = max(value)) %>%
              ungroup() %>%
              mutate(ymin = pmin - (pmax - pmin)*0.05) %>%
              mutate(ymax = pmax + (pmax - pmin)*0.05)


        # performance metrics
        for (j in 1:6) {
          m <- levels(ppm_best$metric)[j]
          m_long <- levels(ppm_best$metric)[j]
          brks_y <- pretty(c(ylimits$ymin[j], ylimits$ymax[j]), n = 2)

          xlabel <- ""
          if(j==5) xlabel <- "Number of checklists"
          ylabel <- m_long

          plot(0, 0, col="white",
            xlim=xlim, ylim=c(brks_y[1], brks_y[length(brks_y)]), 
            xaxt="n", yaxt="n",
            xlab = xlabel, ylab = ylabel, xpd = NA)
          axis(side = 2, at=brks_y)
          if(j>3) axis(side = 1, at = brks_x, labels = brks_x)

          s <- ppm_bad$metric==m
          points(ppm_bad$n_checklists[s], ppm_bad$value[s], pch=16, col=cols[1])

          s <- ppm_best$metric==m
          points(ppm_best$n_checklists[s], ppm_best$value[s], pch=16, col=cols[2])

          text(x = 0.5, y = max(brks_y) - 0.02 * diff(range(brks_y)), 
               labels = LETTERS[j], 
               font = 2)

          # add legend
          if(j==1) legend(x = 2, y = -1.5, 
                 lwd = 2, col = cols, 
                 legend = c("Model 3", "Model 7"), cex = 0.8)

        }

    dev.off()

}

# ####################################################################
# NUMBER OF POSITIVES

log_plot <- FALSE

for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    ppm_best <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 7")
    ppm_bad <- filter(ppm_plot, val_type==plot_val_type, model_name == "Model 3")

    str_glue("{figure_folder}/enc_RF_sample_size_no_pos_{plot_val_type}_{sp_code}.png") %>% 
      png(width = 2200, height = 1500, res = 300)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))
        cols <- c(bad = "grey65", best = "black")

        ppm_min <- c(0, 0.50, 0, 0, 0, 0)
        ppm_max <- c(0.2, 1, 1, 1, 1, 1)

        xlim <- range(c(ppm_best$n_pos, ppm_bad$n_pos))
        xlim_log10 <- c(floor(log10(min(xlim))), ceiling(log10(max(xlim))))

        xlimits <- xlim
        if(log_plot) xlimits <- xlim_log10

        brks_x <- pretty(min(xlimits):max(xlimits), n=2)


        ylimits <- ppm_best %>% 
              rbind(ppm_bad) %>%
              select(metric, value) %>%
              # filter(value < Inf) %>%
              # filter(value > -Inf) %>%
              group_by(metric) %>%
              summarise(pmin = min(value), pmax = max(value)) %>%
              ungroup() %>%
              mutate(ymin = pmin - (pmax - pmin)*0.05) %>%
              mutate(ymax = pmax + (pmax - pmin)*0.05)


        # performance metrics
        for (j in 1:6) {
          m <- levels(ppm_best$metric)[j]
          m_long <- levels(ppm_best$metric)[j]
          brks_y <- pretty(c(ylimits$ymin[j], ylimits$ymax[j]), n = 2)
          xlabel <- ""
          if(j==5) xlabel <- "Number of detections"
          ylabel <- m_long

          plot(0, 0, col="white",
            xlim=c(brks_x[1], brks_x[length(brks_x)]), ylim=c(brks_y[1], brks_y[length(brks_y)]), 
            xaxt="n", yaxt="n",
            xlab = xlabel, ylab = ylabel, xpd = NA)
          axis(side = 2, at=brks_y)

          xtick_labels <- signif(10^brks_x,2)
          if(!log_plot) xtick_labels <- brks_x
          if(j>3) axis(side = 1, at = brks_x, labels = xtick_labels)

          s <- ppm_bad$metric==m
          points(ppm_bad$n_pos[s], ppm_bad$value[s], pch=16, col=cols[1])

          s <- ppm_best$metric==m
          points(ppm_best$n_pos[s], ppm_best$value[s], pch=16, col=cols[2])

          text(x = 20, y = max(brks_y) - 0.02 * diff(range(brks_y)), 
               labels = LETTERS[j], 
               font = 2, pos = 4)

          # add legend
          if(j==1) legend(x = 2, y = -1.5, 
                 lwd = 2, col = cols, 
                 legend = c("Model 3", "Model 7"), cex = 0.8)

        }

    dev.off()

}

# --------------------------------------------------------------------
# plot differences


# plot comparing ppms
diff_plot <- all_ppms %>% 
  mutate(model_name = ifelse(good, "Model 7", "Model 3")) %>%
  select(-n_checklists, -n_pos) %>%
  select_if(~ !is.logical(.)) %>% 
  gather("metric", "value", -run_id, -sim_id, -val_type, -prop_data, -model_name) %>%
  filter(metric != "threshold") %>% 
  arrange(val_type, sim_id, run_id) %>%
  mutate(model_code = ifelse(model_name == "Model 7", "mod7", "mod3")) %>%
  select(-model_name, -run_id) %>%
  spread(model_code, value) %>% 
  mutate(anchor_col = anchor_model) %>%
  mutate(diff = ifelse(anchor_col==3, mod7 - mod3, mod3 - mod7)) %>%
  mutate(metric_short = metric) %>%
  mutate(metric = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
  mutate(metric = factor(metric, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS"))) %>%
  mutate(start = if_else(metric == "AUC", 0.5, 0)) %>%
  filter(!is.na(metric)) %>%
  mutate(prop_data_cat = factor(prop_data, levels=sort(as.numeric(names(table(all_ppms$prop_data))))))


  maxy <- diff_plot %>% 
            select(metric, diff) %>% 
            group_by(metric) %>%
            summarise(max_abs = max(abs(diff))) %>%
            ungroup()


for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    diff_d <- filter(diff_plot, val_type==plot_val_type)

    str_glue("{figure_folder}/enc_RF_sample_size_DIFF_{plot_val_type}_anchor{anchor_model}_{sp_code}.png") %>% 
      png(width = 21, height = 17, units="cm", res = 600)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))

        # performance metrics
        for (j in 1:6) {
          m <- levels(diff_d$metric)[j]
          maxyy <- maxy$max_abs[j]
          xlabel <- ""
          if(j==5) xlabel <- "Proportion of original data"
          boxplot(as.formula("diff ~ prop_data"), 
                  data = diff_d[diff_d$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  boxcol = "transparent", 
                  col = "transparent",
                  xaxt = "n", xlim=c(0.5, 5.5), ylim = c(-1*maxyy, maxyy),
                  xlab = xlabel, ylab = "")
          abline(h=0, lwd=2, col="grey70")
          par(new=TRUE)
          boxplot(as.formula("diff ~ prop_data"), 
                  data = diff_d[diff_d$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
#                  boxcol = "transparent", 
                  xaxt = "n", xlim=c(0.5, 5.5), ylim = c(-1*maxyy, maxyy),
                  xlab = "", 
                  ylab = bquote(Delta~.levels(diff_d$metric)[j]))
          if (j > 3) {
            axis(side = 1, at = c(0.5, 3, 5.5), labels = c(0, 0.5, 1))
          }
          text(x = 0.5, y = maxyy*0.95, 
               labels = LETTERS[j], 
               font = 2, pos=4)

        }

    dev.off()

}


# --------------------------------------------------------------------
# allow y-scales to vary between test datasets


  maxy <- diff_plot %>% 
            select(metric, diff, val_type) %>% 
            group_by(metric, val_type) %>%
            summarise(max_abs = max(abs(diff))) %>%
            ungroup()


add_grey <- TRUE

for(i in 1:4){

    plot_val_type <- c("2017", "bbs", "2017_bal", "bbs_route")[i]

    diff_d <- filter(diff_plot, val_type==plot_val_type)

    str_glue("{figure_folder}/enc_RF_sample_size_DIFF_freey_{plot_val_type}_anchor{anchor_model}_{sp_code}_grey{add_grey}.png") %>% 
      png(width = 21, height = 17, units="cm", res = 600)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1)) # 4

        # performance metrics
        for (j in 1:6) {
          m <- levels(diff_d$metric)[j]
          maxyy <- maxy$max_abs[maxy$val_type==plot_val_type & maxy$metric==m]
          xlabel <- ""
          if(j==5) xlabel <- "Proportion of original data"
          boxplot(as.formula("diff ~ prop_data"), 
                  data = diff_d[diff_d$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  boxcol = "transparent", 
                  col = "transparent",
                  xaxt = "n", xlim=c(0.5, 5.5), 
                  yaxt = "n", ylim = c(-1*maxyy, maxyy),
                  xlab = "", ylab = "")
          if(!add_grey) abline(h=0, lwd=2, col="grey70")

          if(add_grey) {
            ymin <- -1
            ymax <- 0
            if(j==1) {ymin <- 0; ymax <- 1 }
            polygon(x = c(-1, 10, 10, -1, -1), y = c(ymin, ymin, ymax, ymax, ymin), col="grey78", border = alpha("white", 0))
          }
          par(new=TRUE)
          boxplot(as.formula("diff ~ prop_data"), 
                  data = diff_d[diff_d$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                 col = alpha("white", 0.4), 
                  xaxt = "n", xlim=c(0.5, 5.5), ylim = c(-1*maxyy, maxyy),
                  xlab = xlabel, xpd = NA, 
                  ylab = bquote(Delta~.(levels(diff_d$metric)[j])))
          if (j > 3) {
            axis(side = 1, at = c(0.5, 3, 5.5), labels = c(0, 0.5, 1))
          }
          if(j==5) text(x = 3, y = maxyy*-1.45, labels = xlabel, xpd = NA)
          text(x = 0.5, y = maxyy*0.95, 
               labels = LETTERS[j], 
               font = 2, pos=4)

        }

    dev.off()

}



