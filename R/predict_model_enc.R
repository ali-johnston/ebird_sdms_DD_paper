  
predict_model_enc <- function(model, data) {
  # predict on test data set
  if (inherits(model$model, "maxnet")) {
    # predict on standardised fixed test dataset
    pred <- predict(model$model, newdata = data, 
                    type = "logistic", clamp = FALSE)
  } else {
    pred_rf <- predict(model$model, data = data, type = "response")
    pred <- pred_rf$predictions[,2]

    if(!is.null(model$calibration)){
      pred_cal <- predict(model$calibration, 
                        newdata = data.frame(pred = pred), 
                        type = "response")
      pred_cal <- ifelse(pred_cal<0, 0, ifelse(pred_cal>1, 1, pred_cal))
      pred <- pred_cal
    }
  } 
  return(pred)
}
