  
validate <- function(model, data, bbs_combine = FALSE) {

  pred_vec <- predict_model_enc(model, data)
  pred <- data.frame(id = 1:nrow(data), 
                     obs = data$species_observed,
                     pred = pred_vec)

  pred <- drop_na(pred)

  if(bbs_combine){
    pred_df <- data %>%
        cbind(pred = pred_vec) %>%
        mutate(obs = species_observed) %>%
        group_by(route_id) %>%
        summarise(mean_pred = mean(pred), 
                  mean_obs = mean(obs),
                  max_pred = 1 - prod(1-pred),
                  max_obs = max(obs)) %>%
        ungroup()

    pred <- pred_df %>%
      rename(pred = max_pred, obs = max_obs) %>%
      mutate(id = 1:nrow(pred_df)) %>%
      select(id, obs, pred)

  }
  
  # validation metrics
  # mean squared error (mse)
  mse <- mean((pred$obs - pred$pred)^2, na.rm = TRUE)
  # pick threshold to maximize kappa
  thresh <- optimal.thresholds(as.matrix(pred), opt.methods = "MaxKappa")[1, 2]
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
  pa_metrics <- presence.absence.accuracy(as.matrix(pred), threshold = thresh, 
                                          na.rm = TRUE, st.dev = FALSE)
  
  # summarise the performance of this model
  tibble(mse = round(mse, 5),
         auc = round(pa_metrics$AUC, 5),
         # threshold required
         threshold = round(thresh, 5),
         kappa = round(pa_metrics$Kappa, 5), 
         sensitivity = round(pa_metrics$sensitivity, 5),
         specificity = round(pa_metrics$specificity, 5),
         n_checklists = model$n_checklists, 
         n_pos = model$n_sightings
  )
}
