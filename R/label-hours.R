label_hours <- function(x) {
  stopifnot(is.numeric(x))
  
  label_hours_single <- function(y) {
    if (y > 24 * 365) {
      paste(round(y / 24 / 365, 1), "years")
    } else if (y > 24 * 7) {
      paste(round(y / 24 / 7, 1), "weeks")
    } else if (y > 24) {
      paste(round(y / 24, 1), "days")
    } else {
      paste(round(y, 1), "hours")
    }
  }
  vapply(x, label_hours_single, "")
}

