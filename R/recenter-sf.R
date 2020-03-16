recenter_sf <- function(x, crs) {
  stopifnot(inherits(x, c("sf", "sfc")))
  # find central meridian
  center <- as.numeric(stringr::str_extract(crs, "(?<=lon_0=)[-.0-9]+"))
  if (is.na(center) || length(center) != 1 || !is.numeric(center)) {
    stop("CRS has no central meridian term (+lon_0).")
  }
  
  # edge is offset from center by 180 degrees
  edge <- ifelse(center < 0, center + 180, center - 180)
  
  # create an very narrow sliver to clip out
  delta <- 1e-5
  clipper <- sf::st_bbox(c(xmin = edge - delta, xmax = edge + delta,
                           ymin = -90, ymax = 90),
                         crs = 4326)
  clipper <- sf::st_as_sfc(clipper)
  clipper <- suppressWarnings(smoothr::densify(clipper, max_distance = 1e-3))
  clipper <- sf::st_transform(clipper, crs = sf::st_crs(x))
  
  # cut then project
  x_proj <- sf::st_difference(x, clipper)
  sf::st_transform(x_proj, crs = crs)
}