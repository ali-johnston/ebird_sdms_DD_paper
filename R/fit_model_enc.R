
# t <- 2
# maxnet <- bp_runs$maxnet[t]
# complete <- bp_runs$complete[t]
# spatial_subsample <- bp_runs$spatial_subsample[t]
# effort_filter <- bp_runs$effort_filter[t]
# effort_covs <- bp_runs$effort_covs[t]
# data <- model_data
# prop_data <- 1
# calibrate <- TRUE
# subsample_seed <- 1


fit_model_enc <- function(maxnet, incomplete, 
                         complete, spatial_subsample, 
                         effort_filter, effort_covs,
                         data,
                         prop_data = 1, 
                         calibrate = TRUE, 
                         calibrate_plot = FALSE,
                         calibrate_plot_name = NULL,
                         subsample_seed = NULL,
                         n_background = 50000,
                         ...) {

  print("------- Running encounter rate model -------")
  print(paste("raw data n =", nrow(data)))

  if(prop_data < 1) {
    if(is.null(subsample_seed)) subsample_seed <- round(runif(1, 0, 1000))
    set.seed(round(prop_data*1000)+subsample_seed)
    print(paste("Subsampling to", prop_data, "of data"))
    data <- data %>% sample_frac(prop_data)
    print(paste("subsampled to n =", nrow(data)))
    print(head(data$checklist_id))
  }

  if (maxnet) {
    print("Keeping only presences")
    present <- data %>%
      filter(as.logical(species_observed)) %>%
      select(starts_with("pland"))
    # randomly select background points
    bg_n <- n_background
    print(paste(bg_n, "background points"))
    background <- pred_surface %>%
      sample_n(size = bg_n, replace = FALSE) %>%
      select(starts_with("pland"))
    p <- c(rep(1, nrow(present)), rep(0, nrow(background)))
    combined <- bind_rows(present, background)
    print("Running Maxnet model")
    mod <- maxnet(p, combined, maxnet.formula(p, combined, classes = "lq"))
    return(list(model = mod, calibration = NULL, t_max_det = 7,
                n_checklists = nrow(present), 
                n_sightings = nrow(present)))
  } else {

    # incomplete checklists only
    if (incomplete) {
      print("Keeping only incomplete checklists")
      data <- filter(data, !as.logical(all_species_reported))
    }

    # complete checklists only
    if (complete) {
      print("Keeping only complete checklists")
      data <- filter(data, as.logical(all_species_reported))
    }

    # filter on effort covariates
    if (effort_filter){
      print("Filtering by effort")
      data <- data %>%
        filter(protocol_type %in% c("Stationary", "Traveling")) %>%
        mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
        filter(effort_distance_km <= 5,
               duration_minutes <= 5 * 60,
               number_observers <= 10)
    }

    # spatial subsample
    if (spatial_subsample) {
      print("Spatial subsampling")
      data <- filter(data, selected == 1)
    }

    print(paste("processed data n =", nrow(data)))
    
    # select covariates
    if (effort_covs) {
      print("Keeping effort covariates")
      data <- data %>% 
        select(species_observed,
               day_of_year, 
               time_observations_started, duration_minutes,
               effort_distance_km, number_observers, protocol_traveling,
               contains("checklist_calibration_index"),
               starts_with("pland_"))
     
      # ensure there are no rows with NAs
      stopifnot(all(complete.cases(data)))
    } else {
      data <- select(data, species_observed, starts_with("pland_"))
      
      # ensure there are no rows with NAs
      stopifnot(all(complete.cases(data)))
    }

    # define the positive proportion of the training data
    # for balanced random forest
    pos_prop <- mean(as.numeric(as.character(data$species_observed)))
    print(paste("train prevalence =", round(pos_prop, 3)))

    # covert response to factor
    data$species_observed <- factor(data$species_observed, levels=c("0", "1"))


    print("Running random forest")

    # fit model
     mod <- ranger(formula =  species_observed ~ ., 
                  data = data,
                  num.trees = 1000,
                  importance = "impurity",
                  probability = TRUE,
                  mtry = sqrt(ncol(data) - 1),
                  replace = TRUE)
                  # ,sample.fraction = c(pos_prop, pos_prop))
    
    # calibrate
    if(calibrate){
        print("Running calibration")
        logit <- function(p) {log(p/(1-p))}
        mod_pred_train <- predict(mod, data = data, type = "response")
        mod_pred_train <- data.frame(id = nrow(data),
                                     obs = data$species_observed,
                                     pred = mod_pred_train$predictions[,2]) %>% 
          mutate(obs = as.numeric(as.character(obs))) %>%
          drop_na()
        mod_cal <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                        data = mod_pred_train, 
    #                    family = binomial, 
                        gamma = 1.4)

        if(calibrate_plot){
          if(is.null(calibrate_plot_name)) calibrate_plot_name <- paste0("calibration_plot_", Sys.Date(), ".png")
          png(calibrate_plot_name, width = 9, height = 9, units="cm", pointsize = 9, res = 300)
            plot(mod_pred_train$pred, jitter(mod_pred_train$obs) + 0.2 - mod_pred_train$obs*0.4,
              xlim=c(0, 1), ylim=c(0,1), 
              xaxt="n", yaxt= "n",
              xlab = "predicted", ylab = "observed")
            axis(side = 1, at=c(0, 0.5, 1))
            axis(side = 2, at=c(0, 0.5, 1))
            x <- seq(0, 1, by=0.01)
            pred_y <- predict(mod_cal, newdata = data.frame(pred = x))
            pred_y <- ifelse(pred_y<0, 0, ifelse(pred_y>1, 1, pred_y))
            lines(x, pred_y, col="red", lwd=2)
            segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col="yellow")
            segments(x0 = 0, y0 = 0.002, x1 = 1, y1 = 1.002, col="black")
          dev.off()
        }
    }
    if(!calibrate) mod_cal <- NULL
    
    # maximum time of day for detection
    if (effort_covs) {

      print("Finding max time of day")

      # 10 minute sections with sufficient data
      search_hours <- data %>%
        mutate(min_10 = floor(as.numeric(time_observations_started*6))/6) %>%
        count(min_10) %>%
        mutate(pct = n / sum(n)) %>%
        filter(pct > 0.0005 | n>20)

      # coarse search to find approximate time of day
      sub <- data %>%
        mutate(min_10 = floor(as.numeric(time_observations_started*6))/6) %>%
        filter(min_10 >= min(search_hours$min_10), min_10 <= max(search_hours$min_10))

      pred_t <- partial_dependence(mod, vars = "time_observations_started", 
                                   n = c(20, min(c(1000, nrow(sub)))), data = sub) %>%
        rename(prob = `1`) %>%
        top_n(n=3, wt = prob) %>%
        mutate(min_10 = floor(as.numeric(time_observations_started*6))/6)

      # finer resolution search to find accurate time of day
      sub <- data %>%
        mutate(min_10 = floor(as.numeric(time_observations_started*6))/6) %>%
        filter(min_10 <= max(pred_t$min_10), min_10 >= min(pred_t$min_10))

      t_max_det <- partial_dependence(mod, vars = "time_observations_started", 
                                   n = c(20, nrow(sub)), data = sub) %>%
        rename(prob = `1`) %>%
        top_n(n=1, wt=prob) %>%
        select(time_observations_started) %>%
        as.matrix() %>% as.numeric()

      print(t_max_det)

    } 
    if(!effort_covs) {
      t_max_det <- 7
    }
    
    return(list(model = mod, calibration = mod_cal, t_max_det = t_max_det,
                n_checklists = nrow(data), 
                n_sightings = sum(as.numeric(as.character(data$species_observed)))))
  }
}

