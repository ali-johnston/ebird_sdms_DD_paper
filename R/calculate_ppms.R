  
calculate_ppms <- function(pred, obs) { 

  data <- data.frame(id = 1:length(pred), obs = obs, pred = pred) 

  # validation metrics
  print("calculate ppms")
  # mean squared error (mse)
  mse <- mean((obs - pred)^2, na.rm = TRUE)
  # pick threshold to maximize kappa
  thresh <- optimal.thresholds(data, opt.methods = "MaxKappa")[1, 2]
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
  pa_metrics <- presence.absence.accuracy(data, threshold = thresh, 
                                          na.rm = TRUE, st.dev = FALSE)
  
  # summarise the performance of this model
  tibble(mse = round(mse, 5),
         auc = round(pa_metrics$AUC, 5),
         # threshold required
         threshold = round(thresh, 5),
         kappa = round(pa_metrics$Kappa, 5), 
         sensitivity = round(pa_metrics$sensitivity, 5),
         specificity = round(pa_metrics$specificity, 5)
  )
}
