hex_sample <- function(x, spacing = 5, 
                       regime = c("separate", "positive", "negative", "together"),
                       byvar = NULL) {
  stopifnot(is.data.frame(x), 
            c("observation_date", "longitude", "latitude", byvar) %in% names(x))
  stopifnot(is.numeric(spacing), length(spacing) == 1, spacing > 0)
  regime <- match.arg(regime)
  
  x$byvar <- as.character(as.vector(as.matrix(x[,byvar])))

  # generate hexagonal grid
  dggs <- dggridR::dgconstruct(spacing = spacing)

  if(regime!="together"){
    # get hexagonal cell id and week number for each checklist
    x_cell <- x %>% 
      mutate(cell = dggridR::dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
             obs = as.logical(species_observed))

    if (regime == "separate") {
      # sample detections and non-detections separately 
      x_ss <- x_cell %>% 
        dplyr::group_by(obs, byvar, cell) %>% 
        dplyr::sample_n(size = 1) %>% 
        dplyr::ungroup()
    } else if (regime == "positive") {
      # sample only detections  
      x_ss <- x_cell %>% 
        dplyr::filter(obs) %>% 
        dplyr::group_by(byvar, cell) %>% 
        dplyr::sample_n(size = 1) %>% 
        dplyr::ungroup()
      x_ss <- dplyr::bind_rows(x_ss, dplyr::filter(x_cell, !obs))
    } else if (regime == "negative") {
      # sample only non-detections  
      x_ss <- x_cell %>% 
        dplyr::filter(!obs) %>% 
        dplyr::group_by(byvar, cell) %>% 
        dplyr::sample_n(size = 1) %>% 
        dplyr::ungroup()
      x_ss <- dplyr::bind_rows(x_ss, dplyr::filter(x_cell, obs))
    } 
    x_ss_ret <- dplyr::select(x_ss, -obs, -cell, -byvar)

  }

  # subsample ignorant of whether detection or non-detection
  if(regime == "together") {
    # get hexagonal cell id and week number for each checklist
    x_cell <- x %>% 
      mutate(cell = dggridR::dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)

    x_ss <- x_cell %>% 
      dplyr::group_by(byvar, cell) %>% 
      dplyr::sample_n(size = 1) %>% 
      dplyr::ungroup()

    x_ss_ret <- dplyr::select(x_ss, -cell, -byvar)
  }
  return(x_ss_ret)
}
