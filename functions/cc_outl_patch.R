
### Add a patch to cc_outl function so the output is based on Occurrence_ID instead of row.names() in each subset dataframe per taxa

cc_outl_patch <- function (x, lon = "decimalLongitude", lat = "decimalLatitude", 
                           species = "species", method = "quantile", mltpl = 5, tdi = 1000, 
                           value = "clean", sampling_thresh = 0, verbose = TRUE, min_occs = 7, 
                           thinning = FALSE, thinning_res = 0.5) 
{
  match.arg(value, choices = c("clean", "flagged", "ids"))
  match.arg(method, choices = c("distance", "quantile", "mad"))
  
  ## Patch function by adding an Occurrence_ID variable
  x$Occurrence_ID_for_CoordinatesCleaner <- 1:nrow(x)
  
  init_ids <- rownames(x)
  rownames(x) <- NULL
  if (verbose) {
    message("Testing geographic outliers")
  }
  # Split in list of occurrence per taxa
  splist <- split.data.frame(x, f = as.character(x[[species]]))

  splist[[10]]
  
  # Check for duplicates
  test <- lapply(splist, function(k) {
    duplicated(k[, c(species, lon, lat)])
  })
  
  # Count the number of non-duplicated occurrences
  test <- lapply(test, "!")
  test <- as.vector(unlist(lapply(test, "sum")))
  
  # Keep only taxa with the minimum number of occurrences
  splist <- splist[test >= min_occs]
  if (any(test < min_occs)) {
    warning(sprintf("Species with fewer than %o unique records will not be tested.", 
                    min_occs))
  }
  
  if (all(test < min_occs)) # If all taxa have less occurrences than the minimum threshold, provide outputs directly
  {
    # Ouput of the function depend on the "value" argument. Either the cleaned dataset, the vector of binary flags, or the row.names of outliers
    switch(value, clean = return(x), flagged = return(rep(TRUE, nrow(x))), ids = return(init_ids))
  
  } else { # Run the actual test
    
    # Extract number of non-duplicated records per taxa
    record_numbers <- unlist(lapply(splist, nrow))
    
    # Use raster approximation in case of taxa with more than 10000 records
    if (any(record_numbers >= 10000) | thinning) {
      warning("Using raster approximation.")
      ras <- ras_create(x = x, lat = lat, lon = lon, thinning_res = thinning_res)
    }
    
    flags <- lapply(splist, .flagging_patch, thinning = thinning, 
                    lon = lon, lat = lat, method = method, mltpl = mltpl, 
                    ras = ras, sampling_thresh = sampling_thresh, tdi = tdi)
    flags <- as.numeric(as.vector(unlist(flags)))
    flags <- flags[!is.na(flags)]
    
    if (sampling_thresh > 0) {
      pts <- terra::vect(x[flags, c(lon, lat)], geom = c(lon, 
                                                         lat))
      if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
        stop("package 'rnaturalearth' not found. Needed for sampling_cor = TRUE", 
             call. = FALSE)
      }
      if (!requireNamespace("rgbif", quietly = TRUE)) {
        stop("package 'rgbif' not found. Needed for sampling_cor = TRUE", 
             call. = FALSE)
      }
      if (inherits(try(rgbif::occ_count(country = "DEU"), 
                       silent = TRUE), "try-error")) {
        warnings("Could not retrive records number from GBIF, skipping sampling correction")
      }
      else {
        ref <- terra::vect(rnaturalearth::ne_countries(scale = "medium", 
                                                       returnclass = "sf"))
        area <- data.frame(country = ref$iso_a2, area = terra::expanse(ref))
        area <- area[!is.na(area$area), ]
        area <- area[!is.na(area$country), ]
        nrec <- vapply(area$country, FUN = function(k) {
          rgbif::occ_count(country = k)
        }, FUN.VALUE = 1)
        nrec <- data.frame(country = area$country, recs = unlist(nrec), 
                           row.names = NULL)
        nrec_norm <- dplyr::left_join(nrec, area, by = "country")
        nrec_norm$norm <- log(nrec_norm$recs/(nrec_norm$area/1e+06/100))
        ext_over <- terra::extract(ref, pts)
        out <- !is.na(ext_over[!duplicated(ext_over[, 
                                                    1]), 2])
        country <- ext_over[out, "iso_a2"]
        thresh <- stats::quantile(nrec_norm$norm, probs = sampling_thresh)
        s_flagged <- nrec_norm$norm[match(country, nrec_norm$country)]
        s_flagged <- s_flagged > thresh
        s_flagged[is.na(s_flagged)] <- FALSE
        flags <- flags[s_flagged]
      }
    }
    
    out <- rep(TRUE, nrow(x))
    out[flags] <- FALSE
    
    if (verbose) {
      if (value == "ids") {
        message(sprintf("Flagged %s records.", length(flags)))
      }
      else {
        if (value == "clean") {
          message(sprintf("Removed %s records.", sum(!out)))
        }
        else {
          message(sprintf("Flagged %s records.", sum(!out)))
        }
      }
    }
    
    # Ouput of the function depend on the "value" argument. Either the cleaned dataset, the vector of binary flags, or the row.names of outliers
    switch(value, clean = return(x[out, ]), flagged = return(out), ids = return(init_ids[flags]))
  }
}


.flagging_patch <- function(k, thinning, lon, lat, method, mltpl, ras,
                      sampling_thresh, tdi) {
  
  # Set raster flag inc ase thinning= TRUE,
  # or this particualr species has 100000 or more records
  if (nrow(k) >= 10000 | thinning) {
    raster_flag <- TRUE
  } else {
    raster_flag <- FALSE
  }
  
  # calculate the distance matrix between all points for the outlier tests
  ## for small datasets and without thinning, 
  ## simply a distance matrix using geospheric distance
  if (raster_flag) { 
    # raster approximation for large datasets and thinning
    # get a distance matrix of raster midpoints, with the row 
    # and colnames giving the cell IDs
    
    # if the raster distance is used due to large dataset and 
    # not for thinning, account for the number of points per gridcell
    if (thinning) {
      dist_obj <- ras_dist(x = k, 
                           lat = lat, 
                           lon = lon,
                           ras = ras, 
                           weights = FALSE)
      
      # the point IDS
      pts <-  dist_obj$pts
      
      # the distance matrix 
      dist <- dist_obj$dist
      
    } else {
      dist_obj <-  ras_dist(x = k, 
                            lat = lat, 
                            lon = lon,
                            ras = ras, 
                            weights = TRUE)
      
      # the point IDS
      pts <-  dist_obj$pts
      
      # the distance matrix 
      dist <- dist_obj$dist
      
      # a weight matrix to weight each distance by the number of points in it
      wm <- dist_obj$wm
    }
  } else { 
    #distance calculation
    dist <- geosphere::distm(k[, c(lon, lat)], 
                             fun = geosphere::distHaversine) / 1000
    
    # set diagonale to NA, so it does not influence the mean
    dist[dist == 0] <- NA
  }
  
  # calculate the outliers for the different methods
  ##  distance method useing absolute distance
  if (method == "distance") {
    # Drop an error if the sampling correction is activated
    if (sampling_thresh > 0) {
      stop("Sampling correction impossible for method 'distance'" )
    }
    
    # calculate the minimum distance to any next point
    mins <- apply(dist, 1, min, na.rm = TRUE)
    
    # outliers based on absolute distance threshold
    out <- which(mins > tdi)
  } else { # for the other methods the mean must be claculated 
    # depending on if the raster method is used
    # get row means
    if (raster_flag & !thinning) {
      # get mean distance to all other points
      mins <- apply(dist, 1, sum, na.rm = TRUE) / rowSums(wm, na.rm = TRUE)
    } else {
      # get mean distance to all other points
      mins <-  apply(dist, 1, mean, na.rm = TRUE)
    }
  }
  
  ## the quantile method
  if (method == "quantile") {
    # identify the quaniles
    quo <- quantile(mins, c(0.25, 0.75), na.rm = TRUE)
    
    # flag all upper outliers
    out <- which(mins > (quo[2] + stats::IQR(mins) * mltpl))
  }
  
  ## the mad method
  if (method == "mad") {
    # get the median
    quo <- stats::median(mins, na.rm = TRUE)
    
    # get the mad stat
    tester <- stats::mad(mins, na.rm = TRUE)
    
    # Identify outliers
    out <- which(mins > (quo + tester * mltpl))
  }
  
  # Identify the outlier points depending on 
  # if the raster approximation was used
  if (raster_flag) {
    # create output object
    if (length(out) == 0) {
      ret <- NA
    }
    if (length(out) > 0) {
      ret <- rownames(k)[which(pts %in% gsub("X", "", names(out)))]
    }
  }else{
    # create output object
    if (length(out) == 0) {
      ret <- NA
    }
    if (length(out) > 0) {
      # ret <- rownames(k)[out] # If using row.names, they must be different in each taxa-specific data.frame !!!
      
      ## Patch function by extracting the Occurrence_ID_for_CoordinatesCleanerinstead of row.names
      ret <- as.character(k$Occurrence_ID_for_CoordinatesCleaner[out])
    }
  }
  return(ret)
}
