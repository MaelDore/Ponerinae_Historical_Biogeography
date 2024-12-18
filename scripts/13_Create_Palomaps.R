##### Script 13: Create paleomaps of bioregions #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Create paleomaps delineating current bioregions in the past

###

### Inputs

# Coastlines and plate polygons/boundaries at each time step

###

### Outputs

# Paleomaps delineating current bioregions at each time step

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

# remotes::install_github("macroecology/paleoMap")

library(paleoMap)
library(rgdal)
library(sf)
library(sp)
library(tidyverse)
library(openxlsx)
library(rgplates)
library(chronosphere)
library(sfheaders)
library(magick)

# Custom function to retrieve paleomaps from gplates.org
source("./functions/get_Paleomap.R")

### 1.2/ Define starting time frame ####

## Load phylogeny
Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")

# Extract root age
root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_short_names)) ; root_age

## Use rather epoch boundary:
 # Early Creataceous start at 145My

start_time <- 145


### 1.3/ Functions to manipulate sf objects ####

## Function to select only the largest POLYGON out of a MULTIPOLYGON
keep_only_largest_polygon <- function (sf)
{
  geom <- sf$geometry
  geom_class <- class(geom)
  
  if ("sfc_MULTIPOLYGON" %in% geom_class)
  {
    # Cast in POLYGONs
    sf_casted <- st_cast(x = sf, to = "POLYGON")
    # Compute area
    sf_casted$area <- st_area(sf_casted)
    # Identify largest polygon
    ID_max <- which.max(sf_casted$area)
    # Extract only the largest
    sf_selected <- sf_casted[ID_max,] %>% 
      dplyr::select(-area)
    
    # Export output with only the largest polygon
    return(sf_selected)
    
  } else {
    
    if ("sfc_POLYGON" %in% geom_class)
    {
      # If already a single POLYGON, do nothing
      warning("Object is already a single POLYGON\n")
      return(sf)
      
    } else {
      stop("Object is neither a MULTIPOLYGON, or a POLYGON\n")
    }
  }
}


## Function to remove consecutive duplicated coordinates in POLYGONs
remove_consecutive_duplicated_coordinates_from_polygons <- function (sf, round_coordinates = T, digits = 1) 
{
  sp <- as(sf, 'Spatial')
  
  nb_polygons <- length(sp@polygons[[1]]@Polygons)
  
  # Initiate list of polygons ot remove
  polygons_to_remove <- c()
  
  # Loop per polygons
  for (i in 1:nb_polygons)
  {
    if (round_coordinates)
    {
      coords_i <- as.data.frame(round(sp@polygons[[1]]@Polygons[[i]]@coords, digits = digits))
    } else {
      coords_i <- as.data.frame(sp@polygons[[1]]@Polygons[[i]]@coords)
    }
    
    # Merge coordinates
    merged_data <- apply(X = coords_i, MARGIN = 1, FUN = paste, collapse = "_")
    rle(merged_data)$values
    
    # Use rle to identify consecutive duplicates
    RLE_output <- rle(merged_data)
    
    # Identify the rows you want to keep
    row_ID <- cumsum(c(1, RLE_output$lengths[-length(RLE_output$lengths)]))
    
    # Extract only non-consecutive duplicates
    coords_i <- coords_i[row_ID, ]
    # Convert to matrix
    coords_i <- as.matrix(coords_i)
    
    # Check if there are still at least 4 points
    if (nrow(coords_i) < 4)
    {
      polygons_to_remove <- c(polygons_to_remove , i)
    }
    
    # Store back in sp object
    sp@polygons[[1]]@Polygons[[i]]@coords <- coords_i
  }
  
  # Remove polygons that do not have enough points
  if (length(polygons_to_remove > 0))
  {
    sp@polygons[[1]]@Polygons <- sp@polygons[[1]]@Polygons[-polygons_to_remove]
    sp@polygons[[1]]@plotOrder <- sp@polygons[[1]]@plotOrder[-polygons_to_remove]
  }
  
  # Convert back to sf
  output_sf <- st_as_sf(sp)
  
  return(output_sf)
}

## Function to remove consecutive duplicated coordinates in POINTS
remove_consecutive_duplicated_coordinates_from_points <- function (sf, round_coordinates = T, digits = 1) 
{
  sp <- as(sf, 'Spatial')
  
  if (round_coordinates)
  {
    coords_i <- as.data.frame(round(sp@coords, digits = digits))
  } else {
    coords_i <- as.data.frame(sp@coords)
  }
  
  # Merge coordinates
  merged_data <- apply(X = coords_i, MARGIN = 1, FUN = paste, collapse = "_")
  rle(merged_data)$values
  
  # Use rle to identify consecutive duplicates
  RLE_output <- rle(merged_data)
  
  # Identify the rows you want to keep
  row_ID <- cumsum(c(1, RLE_output$lengths[-length(RLE_output$lengths)]))
  
  # Extract only non-consecutive duplicates
  coords_i <- coords_i[row_ID, ]
  # Convert to matrix
  coords_i <- as.matrix(coords_i)
  
  # Store back in sp object
  sp@coords <- coords_i
  
  # Adjust @data
  sp@data <- sp@data[row_ID, ]
  
  # Convert back to sf
  output_sf <- st_as_sf(sp)
  
  return(output_sf)
}


## Function to find the sequence of neighboring lines according to proximity
find_neighboring_sequence <- function (sf, starting_feature = 1, verbose = T)
{
  # Detect number of features to order
  nb_features <- nrow(sf)
  
  # Assign feature ID
  sf <- sf %>%
    mutate(feature_ID = row_number())
  
  # Initiate seq_ID
  sf$seq_ID <- NA
  sf$seq_ID[starting_feature] <- 1
  
  # Loop until all features are ordered
  i <- 1
  current_ID <- starting_feature
  sf_i <- sf
  
  while (i < nb_features)
  {
    # # Compute nearest feature to current feature
    # sf_i <- sf_i %>% 
    #   mutate(nearest_ID = st_nearest_feature(x = sf_i))
    
    # # Find ID of next feature
    # next_ID <- sf_i$feature_ID[sf_i$nearest_ID[current_ID]]
    
    # Extract data with and wihtout current feature 
    feature_i <- sf_i[current_ID, ]
    sf_i <- sf_i[-current_ID, ]
    # row.names(sf_i) <- NULL
    
    # Compute nearest feature to current feature
    feature_i <- feature_i %>% 
      mutate(nearest_ID = st_nearest_feature(x = ., y = sf_i))
    
    # Find ID of next feature
    next_ID <- sf_i$feature_ID[feature_i$nearest_ID[1]]
    
    # Assign order of next feature
    sf$seq_ID[next_ID] <- i+1
    
    # Update current ID for the next iteration
    current_ID <- which(sf_i$feature_ID == next_ID)
    
    # Increment index i
    i <- i + 1
    
    # Print progress
    if (i %% 100 == 0)
    {
      cat(paste0(Sys.time(), " - Order defined for feature n°", i, "/", nb_features,"\n"))
    }
  }
  
  return(sf)
}

## Function to remove outlier points in a sequence
identify_outliers_from_sequence <- function (sf, quantile_neighbors = 0.02, quantile_dist = 0.98, remove = T, plot = T, main = NULL)
{
  # Extract nb of points and neighbors to use
  nb_points <- nrow(sf)
  nb_neighbors <- round(nb_points*0.05, 0)
  
  # Compute pairwise distances
  pairwise_dist <- st_distance(sf)
  diag(pairwise_dist) <- NA
  
  # Loop per points
  mean_dist <- c()
  for (i in 1:nb_points)
  {
    # Identify neighbors based on sequence
    seq_ID_i <- sf$seq_ID[i]
    min_seq <- (seq_ID_i - round(nb_neighbors/2,0) + nb_points) %% nb_points
    max_seq <- (seq_ID_i + round(nb_neighbors/2,0) + nb_points) %% nb_points
    
    if (min_seq < max_seq)
    {
      neigh_i <- min_seq:max_seq
    } else {
      neigh_i <- c(1:max_seq, min_seq:nb_points)
    }
    
    # hist(pairwise_dist[,i])
    # Extract distance to neighbors
    dist_neigh_i <- pairwise_dist[neigh_i,i]
    # Compute mean distance to neighbors
    mean_dist_i <- mean(dist_neigh_i, na.rm = T)
    # Store mean dist
    mean_dist[i] <- mean_dist_i
  }
  
  # Find outliers as the top X% quantile
  Q_threshold <- quantile(x = mean_dist, p = quantile_dist)
  sf$outliers <- mean_dist > Q_threshold
  
  # Plot
  if (plot)
  {
    if (!is.null(main))
    {
      main <- unique(st_drop_geometry(sf[, main]))[1,1, drop = T]
    }
    
    plot(sf[, "outliers"], main = main)
    plot(sf[, "seq_ID"], main = main)
  }
  
  # Remove outliers
  sf <- sf %>% 
    filter(!outliers)
  
  return(sf)
}


## Function to remove outlier points according to distance from initial position
identify_outliers_by_distance_to_initial_position <- function (sf, sf_init, quantile_dist_max = 0.99, quantile_dist_min = 0.10, remove = T, plot = T, main = NULL)
{
  # Filter initial_sf to keep only points that are found in sf
  current_points <- sf$seq_ID
  sf_init_filtered <- sf_init %>% 
    filter(seq_ID %in% current_points)
  
  # Compute pairwise distances
  pairwise_dist <- st_distance(x = sf, y = sf_init_filtered, by_element = T)
  # Find max_outliers as the top X% quantile
  Qmax_threshold <- quantile(x = pairwise_dist, p = quantile_dist_max)
  sf$max_outliers <- pairwise_dist > Qmax_threshold
  # Find min_outliers as the bottom X% quantile
  Qmin_threshold <- quantile(x = pairwise_dist, p = quantile_dist_min)
  sf$min_outliers <- pairwise_dist < Qmin_threshold
  
  # Plot
  if (plot)
  {
    if (!is.null(main))
    {
      main <- unique(st_drop_geometry(sf[, main]))[1,1, drop = T]
    }
    
    plot(sf[, "seq_ID"], main = main)
    plot(sf[, "max_outliers"], main = main)
    plot(sf[, "min_outliers"], main = main)
    
  }
  
  # Remove outliers
  sf <- sf %>% 
    filter(!max_outliers) %>% 
    filter(!min_outliers)
  
  return(sf)
}

## Function to reconstruct initial POLYGONS from ordered POINTS
reconstruct_initial_polygon_with_sfheaders <- function (sf)
{
  sf_crs <- st_crs(sf)
  sp <- as(sf, "Spatial") 
  sp_coords_df <- as.data.frame(sp@coords)
  names(sp_coords_df)[1:2] <- c("Longitude_dec", "Latitude_dec")
  sp_data_df <- as.data.frame(sp@data)
  sp_df <- cbind(sp_coords_df, sp_data_df)
  sp_df <- sp_df %>% 
    arrange(Subregion, seq_ID)
  sf <- sfheaders::sf_polygon(
    obj = sp_df,
    x = "Longitude_dec",
    y = "Latitude_dec",
    polygon_id = "Subregion",
    close = T,
    keep = T
  ) 
  sf::st_crs(sf) <- st_crs(sf_crs)
  return(sf)
}

## Function to split set of points that are crossing the anti-meridian in two
split_by_anti_meridian <- function (sf_points, group_var = "Subregion")
{
  # Extract groups
  group_var <- st_drop_geometry(sf_points[, group_var, drop = T])
  groups_list <- unique(group_var)
  
  # Create East/West hemisphere polygons
  West_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = -180, xmax = 0, ymax = 90, ymin = -90), crs = st_crs(4326))))
  East_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = 0, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))))
  EW_polygons <- rbind(East_hemisphere, West_hemisphere) %>% 
    rename(geometry = x) %>% 
    mutate(Hemisphere = c("East", "West"))
  
  # Initiate final df
  sf_points_fixed <- data.frame()
  
  # Loop per groups
  for (i in seq_along(groups_list))
  {
    # i <- 5
    
    # Extract group and points
    group_i <- groups_list[i]
    sf_points_i <- sf_points[group_var == group_i, ]
    
    # Split in two groups
    sf_use_s2(FALSE)
    sf_points_split_i <- st_intersection(x = sf_points_i, y = EW_polygons)
    sf_use_s2(TRUE)
    
    # Do something only if points in the two groups
    if (length(table(sf_points_split_i$Hemisphere)) == 2)
    {
      # Extract points per groups
      West_points <- sf_points_split_i[sf_points_split_i$Hemisphere == "West", ]
      East_points <- sf_points_split_i[sf_points_split_i$Hemisphere == "East", ]
      
      # plot(West_points[, "seq_ID"])
      # plot(East_points[, "seq_ID"])
      
      ## Fix the West polygon
      
      # Subtract 360° to East points
      East_points_translated <- East_points
      st_geometry(East_points_translated) <- st_geometry(East_points_translated) + c(-360, 0)
      st_crs(East_points_translated) <- st_crs(East_points)
      
      # Create polygon using the sequence
      Fixed_West <- rbind(West_points, East_points_translated) %>% 
        arrange(seq_ID)
      Fixed_West_polygon <- reconstruct_initial_polygon_with_sfheaders(Fixed_West)
      # plot(Fixed_West[, "seq_ID"])
      # plot(Fixed_West_polygon[, "Subregion"])
      
      # Extract the intersect in the West hemisphere
      sf_use_s2(FALSE)
      Fixed_West_polygon <- st_make_valid(st_buffer(Fixed_West_polygon, 0))
      Fixed_West_polygon <- keep_only_largest_polygon(Fixed_West_polygon)
      Fixed_West_polygon <- st_intersection(x = Fixed_West_polygon, y = EW_polygons[2,])
      sf_use_s2(TRUE)
      # plot(Fixed_West_polygon[, "Subregion"])
      
      # Recast to POINTS and update the sequence
      Fixed_West_points <- Fixed_West_polygon %>% 
        st_cast(to = "POINT") %>% 
        mutate(seq_ID = row_number()) %>% 
        select(-Hemisphere, -Hemisphere.1) %>% 
        mutate(Subregion = paste0(group_i, "_West"))
      
      # plot(Fixed_West_points[, "seq_ID"])
      
      ## Fix the East polygon
      
      # Add 360° to West points
      West_points_translated <- West_points
      st_geometry(West_points_translated) <- st_geometry(West_points_translated) + c(360, 0)
      st_crs(West_points_translated) <- st_crs(West_points)
      
      # Create polygon using the sequence
      Fixed_East <- rbind(East_points, West_points_translated) %>% 
        arrange(seq_ID)
      Fixed_East_polygon <- reconstruct_initial_polygon_with_sfheaders(Fixed_East)
      # plot(Fixed_East[, "seq_ID"])
      # plot(Fixed_East_polygon[, "Subregion"])
      
      # Extract the intersect in the East hemisphere
      sf_use_s2(FALSE)
      Fixed_East_polygon <- st_make_valid(st_buffer(Fixed_East_polygon, 0))
      Fixed_East_polygon <- keep_only_largest_polygon(Fixed_East_polygon)
      Fixed_East_polygon <- st_intersection(x = Fixed_East_polygon, y = EW_polygons[1,])
      sf_use_s2(TRUE)
      # plot(Fixed_East_polygon[, "Subregion"])
      
      # Recast to POINTS and update the sequence
      Fixed_East_points <- Fixed_East_polygon %>% 
        st_cast(to = "POINT") %>% 
        mutate(seq_ID = row_number()) %>% 
        select(-Hemisphere, -Hemisphere.1) %>% 
        mutate(Subregion = paste0(group_i, "_East"))
      
      # plot(Fixed_East_points[, "seq_ID"])
      
      ## Merge the fixed East and West points with updated Subregion names
      sf_points_fixed_i <- rbind(Fixed_East_points, Fixed_West_points)
      # plot(sf_points_fixed_i[, c("Subregion", "seq_ID")])
    } else {
      sf_points_fixed_i <- sf_points_split_i %>%
        select(-Hemisphere)
    }
    ## Merge all the fixed East and West points per groups
    sf_points_fixed <- rbind(sf_points_fixed, sf_points_fixed_i)
  }
  return(sf_points_fixed)
}


##### 2/ Download paleomaps from gplates.org ####

### 2.1/ Explore different models ####

# Available models for coastlines, and associated time frame
  # See https://gwsdoc.gplates.org/models
  # MULLER2022: 0 - 1000 My
  # MERDITH2021: 0 - 1000 My
  # MULLER2019: 0 - 250 My
  # SETON2012: 0-200

# Try several models to compare maps
test_sp <- get_Paleomap(ma = c(50, 100, 150), model = "MULLER2019", show.plates = T, do.plot = T)
test_sp <- get_Paleomap(ma = c(100 ,130), model = "SETON2012", show.plates = T, do.plot = T)
test_sp <- get_Paleomap(ma = c(100, 130), model = "MULLER2019", show.plates = T, do.plot = T)
test_sp <- get_Paleomap(ma = c(100, 130), model = "MULLER2022", show.plates = T, do.plot = T)
test_sp <- get_Paleomap(ma = c(100, 130), model = "MERDITH2021", show.plates = T, do.plot = T)
test_sp <- get_Paleomap(ma = c(100, 130), model = "MULLER2016", show.plates = T, do.plot = T)

# Choose the "best" model and check availability of all time frames
# No model agrees that Palearctics were not connected...

## Select the most recent model: MULLER2022
test_sp <- get_Paleomap(ma = c(0, 77), model = "MULLER2022", get.plate.polygons = T, get.plate.boundaries = T, show.plates = T, do.plot = T)
Paleomaps_sp <- get_Paleomap(ma = 0:start_time, model = "MULLER2022", get.plate.polygons = T, get.plate.boundaries = F, show.plates = F, do.plot = F)

str(test_sp, max.level = 2)
str(Paleomaps_sp, max.level = 2)

test_sp$Plate_polygons$Plate_polygons_0My_MULLER2022
Paleomaps_sp$Plate_polygons$Plate_polygons_0My_MULLER2022

## Save sp objects
# saveRDS(object = test_sp, file = "./outputs/Paleomaps/test_sp.rds")
saveRDS(object = Paleomaps_sp, file = "./outputs/Paleomaps/Paleomaps_sp.rds")


### 2.2/ Convert outputs to sf ####

## Load sp objects
Paleomaps_sp <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sp.rds")


Paleomaps_sf <- Paleomaps_sp
for (i in seq_along(Paleomaps_sf))
{
  # i <- 1
  
  # Extract focal type of maps
  maps_type_i <- Paleomaps_sp[[i]]
  
  # Convert to sf
  maps_type_sf_i <- lapply(X = maps_type_i, FUN = st_as_sf)
  
  # Add shape_ID
  maps_type_sf_i <- lapply(X = maps_type_sf_i, FUN = function (x) { y <- x ; y$shape_ID <- as.factor(1:nrow(x)) ; return(y) })
  
  # Replace with converted maps
  Paleomaps_sf[[i]] <- maps_type_sf_i
}

nrow(Paleomaps_sf$Plate_polygons$Plate_polygons_0My_MULLER2022)
nrow(Paleomaps_sf$Plate_polygons$Plate_polygons_77My_MULLER2022)

# Quick plot
g <- ggplot(data = Paleomaps_sf$Coastlines$Coastlines_10My_MULLER2022) + 
  
  geom_sf(data = Paleomaps_sf$Coastlines$Coastlines_10My_MULLER2022,
          fill = "limegreen", col = "black") +
  
  geom_sf(data = Paleomaps_sf$Plate_polygons$Plate_polygons_10My_MULLER2022,
          mapping = aes(fill = shape_ID), alpha = 0.2, col = "red")
  
  # geom_sf(data = Paleomaps_sf$Plate_boundaries$Coastlines_10My_MULLER2022,
  #         col = "purple")
  
print(g)


### 2.3/ Merge/union coastline polygons to delete internal boundaries ####

# Extract coastline maps
coastline_maps <- Paleomaps_sf[["Coastlines"]]
for (i in seq_along(coastline_maps))
{
  # i <- 1
  
  # Extract map for Time i
  coastline_map_i <- coastline_maps[[i]]
  
  plot(coastline_map_i[, "shape_ID"])
  sf::sf_use_s2(FALSE)
  coastline_map_union_i <- st_make_valid(st_buffer(x = st_union(st_make_valid(coastline_map_i)), dist = 0))
  sf::sf_use_s2(TRUE)
  
  coastline_map_union_i <- coastline_map_union_i %>% 
    st_as_sf() %>% 
    # dplyr::summarise() %>% 
    dplyr::mutate(age = coastline_map_i$age[1],
                  model = coastline_map_i$model[1],
                  name = coastline_map_i$name[1],
                  shape_ID = row_number())
                  # area = st_area(.))
  
  plot(coastline_map_union_i[, "shape_ID"])
  
  # Store unioned map
  coastline_maps[[i]] <- coastline_map_union_i  
  
}
Paleomaps_sf[["Coastlines"]] <- coastline_maps


### 2.4/ Save paleomaps in sf format

# saveRDS(object = test_sf, file = "./outputs/Paleomaps/test_sf.rds")
saveRDS(object = Paleomaps_sf, file = "./outputs/Paleomaps/Paleomaps_sf.rds")


##### 3/ Attribute bioregions to plate polygons #####
  
## Need to attribute bioregions to plate polygons

# Use maps to light plates one by one
# Manually provide classification for each plate
# Leave NA for plates encompassing multiple bioregions

## Load paleomaps in sf format
# test_sf <- readRDS(file = "./outputs/Paleomaps/test_sf.rds")
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")


### 3.1/ Extract plate centroids ####

Paleomaps_sf$Plate_centroids <- Paleomaps_sf$Plate_polygons
for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  Plate_polygons_i <- Paleomaps_sf$Plate_polygons[[i]]
  
  sf::sf_use_s2(FALSE)
  Plate_centroids_i <- st_centroid(Plate_polygons_i)
  sf::sf_use_s2(TRUE)
  
  # Extract WGS84 coordinates
  coordinates <- st_as_text(Plate_centroids_i$geometry)
  longitude_dec <- str_remove(string = coordinates, pattern = ".* \\(")
  Plate_centroids_i$longitude_dec <- as.numeric(str_remove(string = longitude_dec, pattern = " .*"))
  latitude_dec <- str_remove(string = coordinates, pattern = ".* \\(")
  latitude_dec <- str_remove(string = latitude_dec, pattern = ".* ")
  Plate_centroids_i$latitude_dec <- as.numeric(str_remove(string = latitude_dec, pattern = "\\)"))
  
  # Extract Mollweide coordinates
  Plate_centroids_Mollweide_i <- st_transform(x = Plate_centroids_i, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  coordinates_Mollweide <- st_as_text(Plate_centroids_Mollweide_i$geometry)
  longitude_dec_Mollweide <- str_remove(string = coordinates_Mollweide, pattern = ".* \\(")
  Plate_centroids_i$longitude_dec_Mollweide <- as.numeric(str_remove(string = longitude_dec_Mollweide, pattern = " .*"))
  latitude_dec_Mollweide <- str_remove(string = coordinates_Mollweide, pattern = ".* \\(")
  latitude_dec_Mollweide <- str_remove(string = latitude_dec_Mollweide, pattern = ".* ")
  Plate_centroids_i$latitude_dec_Mollweide <- as.numeric(str_remove(string = latitude_dec_Mollweide, pattern = "\\)"))
  
  # Store output
  Paleomaps_sf$Plate_centroids[[i]] <- Plate_centroids_i
  
}

# Save updates with centroids
# saveRDS(object = Paleomaps_sf, file = "./outputs/Paleomaps/Paleomaps_sf.rds")


### 3.2/ Plot plate labels ####

pdf(file = "./outputs/Paleomaps/Plates_all_ages_MULLER2022.pdf", width = 16, height = 8)

for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[i]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[i]]
  # Plate_boundaries_i <- test_sf$Plate_boundaries[[i]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[i]]
  
  # Extract age and model
  age_i <- Coastlines_i$age[1]
  model_i <- Coastlines_i$model[1]
  
  # Create plot
  plates_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "limegreen", col = "black") +
    
    geom_sf(data = Plate_polygons_i,
            mapping = aes(fill = shape_ID), alpha = 0.2, col = "red") +
    
    # # WGS84 labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec, y = latitude_dec), col = "black") +
  
    # Mollweide labels
    geom_label(data = Plate_centroids_i,
               mapping = aes(label = shape_ID, x = longitude_dec_Mollweide, y = latitude_dec_Mollweide), col = "black") +
    
    
    # geom_sf(data = Plate_boundaries_i,
    #         col = "purple") +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Tectonic plates in ", model_i, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 12, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(plates_plot_i)
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Tectonic plates plotted for model ",model_i," at age ", age_i," My - n°", i, "/", length(Paleomaps_sf$Plate_polygons),"\n"))
  }
}

dev.off()


### 3.3/ Export dataframes of bioregion membership of plates to fill ####

Plate_bioregion_membership_df_list <- list()

for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  # Extract data
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[i]]
  nb_plates <- nrow(Plate_polygons_i)
  
  # Extract age and model
  age_i <- Plate_polygons_i$age[1]
  model_i <- Plate_polygons_i$model[1]
  
  # Generate df
  Plate_bioregion_membership_df_i <- data.frame(Plate_ID = 1:nb_plates, Bioregion = NA)
  
  # Store df
  old_names <- names(Plate_bioregion_membership_df_list)
  Plate_bioregion_membership_df_list[[i]] <- Plate_bioregion_membership_df_i
  new_names <- c(old_names, paste0(model_i, "_", age_i, "My"))
  names(Plate_bioregion_membership_df_list) <- new_names
}


## Export df as Excel file to fill by looking at tectonic plate plots

openxlsx::write.xlsx(x = Plate_bioregion_membership_df_list, file = "./outputs/Paleomaps/Plate_bioregion_membership_df_list_to_fill.xlsx")


### 3.4/ Import bioregion membership of plates ####

## 3.4.1/ Import filled df
Plate_bioregion_membership_df_list_filled <- Plate_bioregion_membership_df_list
for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  Plate_bioregion_membership_df_list_filled[[i]] <- openxlsx::read.xlsx(xlsxFile = "./outputs/Paleomaps/Plate_bioregion_membership_df_list_filled.xlsx", sheet = i, na.strings = F)
}

## 3.4.2/ Join data in sf objects

for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  # Extract data
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[i]]
  Plate_bioregion_membership_df_i <- Plate_bioregion_membership_df_list_filled[[i]]
  
  # Join data in sf object
  Plate_polygons_i$shape_ID <- as.numeric(Plate_polygons_i$shape_ID)
  # Plate_polygons_i <- Plate_polygons_i %>%
  #   select(-Bioregion)
  Plate_polygons_i <- dplyr::left_join(x = Plate_polygons_i, y = Plate_bioregion_membership_df_i, by = join_by("shape_ID" == "Plate_ID"))
  Plate_polygons_i$shape_ID <- as.factor(Plate_polygons_i$shape_ID)
  
  # Store output
  Paleomaps_sf$Plate_polygons[[i]] <- Plate_polygons_i
}

# Paleomaps_sf$Plate_polygons[[i]] <- Paleomaps_sf$Plate_polygons[[i]] %>%
#   select(-Bioregion)

# Save updates with bioregion membership for plaques
# saveRDS(object = Paleomaps_sf, file = "./outputs/Paleomaps/Paleomaps_sf.rds")


## 3.4.3/ Plot result

# plot(Paleomaps_sf$Plate_polygons[[i]][,"Bioregion"])

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA

pdf(file = "./outputs/Paleomaps/Plates_all_ages_MULLER2022_with_bioregions.pdf", width = 16, height = 8)

for (i in seq_along(Paleomaps_sf$Plate_polygons))
{
  # i <- 1
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[i]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[i]]
  # Plate_boundaries_i <- Paleomaps_sf$Plate_boundaries[[i]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[i]]
  
  # Extract age and model
  age_i <- Coastlines_i$age[1]
  model_i <- Coastlines_i$model[1]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Plate_polygons_i$Bioregion)]
  
  # Order bioregion factors
  Plate_polygons_i$Bioregion <- factor(x = Plate_polygons_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  
  
  # Create plot
  plates_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Plate_polygons_i,
            mapping = aes(fill = Bioregion), alpha = 0.5, col = "grey50") +
    
    # # WGS84 labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec, y = latitude_dec), col = "black") +
    
    # Mollweide labels
    geom_label(data = Plate_centroids_i,
               mapping = aes(label = shape_ID, x = longitude_dec_Mollweide, y = latitude_dec_Mollweide),
               col = "#00000066", alpha = 0.5) +
    
    
    # geom_sf(data = Plate_boundaries_i,
    #         col = "purple") +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Tectonic plates in ", model_i, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 12, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(plates_plot_i)
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Tectonic plates with Bioregions plotted for model ",model_i," at age ", age_i," My - n°", i, "/", length(Paleomaps_sf$Plate_polygons),"\n"))
  }
}

dev.off()



##### 4/ Refine bioregion attribution within plates using reconstructed subplates ####

### 4.1/ Aggregate adm1 in subregions (Bioregions x plaques) ####

# Load worldwide adm0 GeoBoundaries map
geoBoundaries_Countries_sf <- readRDS(file = "./input_data/geoBoundaries/geoBoundaries_Countries_sf.rds")
View(geoBoundaries_Countries_sf)

# Load Bioregion x adm1 sf
Bioregions_sf_adm1_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level.rds")
View(Bioregions_sf_adm1_level)

plot(Bioregions_sf_adm1_level[, "Bioregion"])

# Use countries/adm1 to divide Bioregion in Subregions as Bioregion x plates boundaries

Bioregions_sf_adm1_level_with_subregions <- Bioregions_sf_adm1_level
Bioregions_sf_adm1_level_with_subregions$Subregion <- NA

table(Bioregions_sf_adm1_level_with_subregions$Bioregion)

## Nearctic: Mexico portion vs. Greenland vs. the rest

# Mexican
Mexican_Nearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Nearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB == "MEX")
Bioregions_sf_adm1_level_with_subregions$Subregion[Mexican_Nearctic_indices] <- "Mexican_Nearctic"

# Greenland
Greenland_Nearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Nearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB == "GRL")
Bioregions_sf_adm1_level_with_subregions$Subregion[Greenland_Nearctic_indices] <- "Greenland_Nearctic"

# American (the rest)
American_Nearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Nearctic" & !Bioregions_sf_adm1_level_with_subregions$Subregion %in% c("Mexican_Nearctic", "Greenland_Nearctic"))
Bioregions_sf_adm1_level_with_subregions$Subregion[American_Nearctic_indices] <- "American_Nearctic"

# View(Bioregions_sf_adm1_level_with_subregions[Bioregions_sf_adm1_level_with_subregions$Bioregion == "Nearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB == "USA", ])

# Correct Puerto Rico & United States Virgin Islands for Caribbean_Neotropics
Puerto_Rico_indices <- which(Bioregions_sf_adm1_level_with_subregions$shapeName == "Puerto Rico")
USVI_indices <- which(Bioregions_sf_adm1_level_with_subregions$shapeName == "United States Virgin Islands")
Bioregions_sf_adm1_level_with_subregions$Bioregion[c(Puerto_Rico_indices, USVI_indices)] <- "Neotropics"
Bioregions_sf_adm1_level_with_subregions$Subregion[c(Puerto_Rico_indices, USVI_indices)] <- "Caribbean_Neotropics"

# Correct American Samoa, Guam & Marianne Islands for Oceanian_Australasia
American_Samoa_indices <- which(Bioregions_sf_adm1_level_with_subregions$shapeName == "American Samoa")
Guam_indices <- which(Bioregions_sf_adm1_level_with_subregions$shapeName == "Guam")
Mariana_indices <- which(Bioregions_sf_adm1_level_with_subregions$shapeName == "Commonwealth of the Northern Mariana Islands")
Bioregions_sf_adm1_level_with_subregions$Bioregion[c(American_Samoa_indices, Guam_indices, Mariana_indices)] <- "Australasia"
Bioregions_sf_adm1_level_with_subregions$Subregion[c(American_Samoa_indices, Guam_indices, Mariana_indices)] <- "Oceanian_Australasia"

## Neotropics: Mesoamerica vs. Caribbean vs. South America

# Mesoamerican: Mexico (South), Guatemala, Belize, Honduras, Salvador, Nicaragua, Costa Rica, Panama
Mesoamerican_Neotropics_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Neotropics" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MEX", "GTM", "BLZ", "HND", "SLV", "NIC", "CRI", "PAN"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Mesoamerican_Neotropics_indices] <- "Mesoamerican_Neotropics"

# South_American: Colombia, Venezuela, Guyana, Suriname, (Guyane Française), Ecuador, Brazil, Peru, Bolivia, Paraguay, Chili, Uruguay, Argentina, Falkland Islands, Isla Brazilera
South_American_Neotropics_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Neotropics" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("COL", "VEN", "GUY", "SUR", "ECU", "BRA", "PER", "BOL", "PRY", "CHL", "URY", "ARG", "117", "120"))
Bioregions_sf_adm1_level_with_subregions$Subregion[South_American_Neotropics_indices] <- "South_American_Neotropics"

# Caribbean: the rest
# View(Bioregions_sf_adm1_level_with_subregions[Bioregions_sf_adm1_level_with_subregions$Bioregion == "Neotropics", ])
Caribbean_Neotropics_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Neotropics" & !Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MEX", "GTM", "BLZ", "HND", "SLV", "NIC", "CRI", "PAN", "COL", "VEN", "GUY", "SUR", "ECU", "BRA", "PER", "BOL", "PRY", "CHL", "URY", "ARG", "117", "120"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Caribbean_Neotropics_indices] <- "Caribbean_Neotropics"

## Afrotropics: Magalasy region vs. Continent

# Magalasy = Madagascar, Comoros, Mauritius, (Mayotte), (Reunion), Seychelles
Magalasy_Afrotropics_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Afrotropics" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MDG", "COM", "MUS", "SYC"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Magalasy_Afrotropics_indices] <- "Magalasy_Afrotropics"

# Continental = the rest
Continental_Afrotropics_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Afrotropics" & !Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MDG", "COM", "MUS", "SYC"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Continental_Afrotropics_indices] <- "Continental_Afrotropics"

## Western Palearctic: North Africa vs. Arabic peninsula vs. South Europa (Portugal/Spain, Italy, Coastal Balkans, Turkey) vs. the rest

# North_African = Morocco, Western Sahara, Algeria, Tunisia, Libya, Egypt, Abyei
North_African_Western_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Western Palearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MAR", "ESH", "DZA", "TUN", "LBY", "EGY", "111"))
Bioregions_sf_adm1_level_with_subregions$Subregion[North_African_Western_Palearctic_indices] <- "North_African_Western_Palearctic"

# Arabic_peninsula = Saudi Arabia, Yemen, Oman, UAE, Qatar, Bahrain, Kuwait, Iraq, Jordan, Israel, Lebanon, Syria, West Bank, Gaza Strip, Sanafir & Tiran Is.
Arabic_Peninsula_Western_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Western Palearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("SAU", "YEM", "OMN", "ARE", "QAT", "BHR", "KWT", "IRQ", "JOR", "ISR", "LBN", "SYR", "129", "118", "126"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Arabic_Peninsula_Western_Palearctic_indices] <- "Arabic_Peninsula_Western_Palearctic"

# South_Europa: Portugal, Spain, Andorra, Italy, Vatican, San Marino, Malta, Slovenia, Croatia, Bosnia-Herzegovina, Montenegro, Albania, Greece, Turkey, Cyprus, No Man's Land, Dragonja
South_Europa_Western_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Western Palearctic" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("PRT", "ESP", "AND", "ITA", "SMR", "VAT", "MLT", "SVN", "HRV", "BIH", "MNE", "ALB", "GRC", "TUR", "CYP", "124", "115"))
Bioregions_sf_adm1_level_with_subregions$Subregion[South_Europa_Western_Palearctic_indices] <- "South_Europa_Western_Palearctic"

# NE_Europa: the rest
NE_Europa_Western_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Western Palearctic" & !Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("MAR", "ESH", "DZA", "TUN", "LBY", "EGY", "111", "SAU", "YEM", "OMN", "ARE", "QAT", "BHR", "KWT", "IRQ", "JOR", "ISR", "LBN", "SYR", "129", "118", "126", "PRT", "ESP", "AND", "ITA", "SMR", "VAT", "MLT", "SVN", "HRV", "BIH", "MNE", "ALB", "GRC", "TUR", "CYP", "124", "115"))
Bioregions_sf_adm1_level_with_subregions$Subregion[NE_Europa_Western_Palearctic_indices] <- "NE_Europa_Western_Palearctic"

# plot(Bioregions_sf_adm1_level_with_subregions[Bioregions_sf_adm1_level_with_subregions$Bioregion == "Western Palearctic", "Subregion"])

## Eastern Palearctic: Far_East (China, Mongolia, Siberia, Japan, Koreas, Himalayan Provinces) vs. Central (Central Asia + Middle East + Oural border)

# Far_East = Palearctic China, Mongolia, Siberia Russia, Japan, Koreas, Himalayan Provinces (Palearctic India and Pakistan), Aksai Chin, Demchok, Dramana-Shakatoe, Kalapani, Siachen-Saltoro, Liancourt Rocks
Siberian_adm1 <- c("Krasnoyarsk Krai", "Khakassia", "Tuva", "Irkutsk Oblast", "Buryatia", "Zabaykalsky Krai", "Amur Oblast", "Sakha Republic", "Jewish Autonomous Oblast", "Khabarovsk Krai", "Primorsky Krai", "Sakhalin Oblast", "Magadan Oblast", "Kamchatka Krai", "Chukotka Autonomous Okrug")
# Himalayan_adm1 <- c("Azad Kashmir", "Gilgit-Baltistan", "Him?chal Pradesh", "Himachal Pradesh", "Jammu and Kashm?r", "Jammu and Kashmir", "Ladakh", "Lad?kh")
Far_East_Eastern_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Eastern Palearctic" & (Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("CHN", "MNG", "JPN", "PRK", "KOR", "IND", "PAK", "112", "113", "114", "116", "119", "121", "123") | Bioregions_sf_adm1_level_with_subregions$shapeName %in% Siberian_adm1))
Bioregions_sf_adm1_level_with_subregions$Subregion[Far_East_Eastern_Palearctic_indices] <- "Far_East_Eastern_Palearctic"

# Central = Central Asia + Middle East + Oural border
Central_Eastern_Palearctic_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Eastern Palearctic" & !(Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("CHN", "MNG", "JPN", "PRK", "KOR", "IND", "PAK", "112", "113", "114", "116", "119", "121", "123") | Bioregions_sf_adm1_level_with_subregions$shapeName %in% Siberian_adm1))
Bioregions_sf_adm1_level_with_subregions$Subregion[Central_Eastern_Palearctic_indices] <- "Central_Eastern_Palearctic"

## Indo-Malaya: Indian subcontinent (India/Pakistan/Bhoutan/Nepal/Sri Lanka/Bangladesh/Maldives) vs. SE Archipelago vs. Philippines vs. Continental SE-Asia

# Indian_subcontinent = India, Pakistan, Bhoutan, Nepal, Sri Lanka, Bangladesh, Maldives, Aksai Chin, Demchok, Dramana-Shakatoe, Kalapani, Siachen-Saltoro
Indian_subcontinent_Indomalaya_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Indomalaya" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("IND", "PAK", "BTN", "NPL", "LKA", "BGD", "MDV", "112", "113", "114", "116", "119", "121"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Indian_subcontinent_Indomalaya_indices] <- "Indian_subcontinent_Indomalaya"

# Philippines = Philippines, Spartly Islands
Philippinian_Indomalaya_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Indomalaya" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("PHL", "128"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Philippinian_Indomalaya_indices] <- "Philippinian_Indomalaya"

# SE_Archipelago = IndoMalayan Indonesia, Brunei, Bornean Malaysia
Bornean_Malaysia_adm1 <- c("Sabah", "Sarawak")
SE_Archipelago_Indomalaya_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Indomalaya" & (Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("IDN", "BRN") | Bioregions_sf_adm1_level_with_subregions$shapeName %in% Bornean_Malaysia_adm1))
Bioregions_sf_adm1_level_with_subregions$Subregion[SE_Archipelago_Indomalaya_indices] <- "SE_Archipelago_Indomalaya"

# Continental SE-Asia = the rest (includes Taiwan, Paracels Islands, Okinawa, Senkakus Islands)
Continental_SE_Asia_Indomalaya_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Indomalaya" & !(Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("IND", "PAK", "BTN", "NPL", "LKA", "BGD", "MDV", "112", "113", "114", "116", "119", "121", "PHL", "128", "IDN", "BRN") | Bioregions_sf_adm1_level_with_subregions$shapeName %in% Bornean_Malaysia_adm1))
Bioregions_sf_adm1_level_with_subregions$Subregion[Continental_SE_Asia_Indomalaya_indices] <- "Continental_SE_Asia_Indomalaya"

## Australasia: Australia vs. PNG/Indonesia/Timor vs. N-Z vs. Oceania

# East_Wallacea = Australasian Indonesia, PNG, Timor Leste
East_Wallacea_Australasia_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Australasia" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("IDN", "PNG", "TLS"))
Bioregions_sf_adm1_level_with_subregions$Subregion[East_Wallacea_Australasia_indices] <- "East_Wallacea_Australasia"

# Australia = Australia, (Cocos/Keeling, Christmas)
Australian_Australasia_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Australasia" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("AUS"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Australian_Australasia_indices] <- "Australian_Australasia"

# Zealandia = N-Z, (Norfolk), (New Caledonia)
Zealandia_Australasia_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Australasia" & Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("NZL"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Zealandia_Australasia_indices] <- "Zealandia_Australasia"

# Oceania = the rest
Oceanian_Australasia_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Australasia" & !Bioregions_sf_adm1_level_with_subregions$ISO_A3_for_GB %in% c("IDN", "PNG", "TLS", "AUS", "NZL"))
Bioregions_sf_adm1_level_with_subregions$Subregion[Oceanian_Australasia_indices] <- "Oceanian_Australasia"

# View(Bioregions_sf_adm1_level_with_subregions[Bioregions_sf_adm1_level_with_subregions$Subregion == "Oceanian_Australasia", ])

## Antarctica = ATA

# Antarctica: Antarctica
Antarctica_indices <- which(Bioregions_sf_adm1_level_with_subregions$Bioregion == "Antarctica")
Bioregions_sf_adm1_level_with_subregions$Subregion[Antarctica_indices] <- "Antarctica"


## Check missing adm1 with no Subregions
table(is.na(Bioregions_sf_adm1_level_with_subregions$Subregion))
table(Bioregions_sf_adm1_level_with_subregions$Subregion)
table(Bioregions_sf_adm1_level_with_subregions$Bioregion)

## Plot result

pdf(file = "./outputs/Paleomaps/Subregions_map.pdf", width = 12, height = 6)

ggplot(data = Bioregions_sf_adm1_level_with_subregions) +
  
  geom_sf(mapping = aes(fill = Subregion)) +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  # Adjust labels
  xlab("") +
  ylab("") +
  
  # Adjust title
  ggtitle(label = paste0("Subregions")) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 8, face = "bold"),
        legend.box.margin = margin(l = 5),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

dev.off()


## Save before aggregation
saveRDS(object = Bioregions_sf_adm1_level_with_subregions, file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level_with_subregions.rds")


## Aggregate adm1 to subregion level

# Load adm1 with subregions membership
Bioregions_sf_adm1_level_with_subregions <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level_with_subregions.rds")

# Issue for shapes crossing the anti-meridian?

sf_use_s2(FALSE)
Bioregions_sf_Subregions_level <- Bioregions_sf_adm1_level_with_subregions %>%
  group_by(Subregion) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Bioregions_sf_Subregions_level <- st_buffer(st_make_valid(Bioregions_sf_Subregions_level), 0)
sf_use_s2(TRUE)

plot(Bioregions_sf_Subregions_level[, "Subregion"])

## Save after aggregation
saveRDS(object = Bioregions_sf_Subregions_level, file = "./input_data/geoBoundaries/Bioregions_sf_Subregions_level.rds")


## 4.2/ Extant subregions to (sub-)plaques using Voronoï cells ####

# Voronoi cells around the Subregions to get subplaques covering the whole Earth

## Load Subregions_sf
Bioregions_sf_Subregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Subregions_level.rds")

## 4.2.1/ Simplify and clean the shapes to fasten computation ####

# Remove tiny shapes < 50 km²

Bioregions_sf_Subregions_level_casted <- st_cast(Bioregions_sf_Subregions_level, to = 'POLYGON')
nrow(Bioregions_sf_Subregions_level) # 23 Subregions
nrow(Bioregions_sf_Subregions_level_casted) # 20033 initial polygons

# View(Bioregions_sf_Subregions_level_casted) 

sf_use_s2(FALSE)
Bioregions_sf_Subregions_level_casted <- Bioregions_sf_Subregions_level_casted %>% 
  mutate(area = units::set_units(st_area(.), km^2)) %>%  # in km²
  filter(area >= units::as_units(50, "km^2")) # Filter our small polygons
Bioregions_sf_Subregions_level_casted <- st_buffer(st_make_valid(Bioregions_sf_Subregions_level_casted), 0)
sf_use_s2(TRUE)

nrow(Bioregions_sf_Subregions_level_casted) # 1634 polygons after removal of small polygons

plot(Bioregions_sf_Subregions_level_casted[, "Subregion"])

# Check length of each polygons to detect errors (should be >4)
nb_points_in_polygons <- c()
for (i in 1:nrow(Bioregions_sf_Subregions_level_casted))
{
  nb_points_in_polygons[i] <- length(st_cast(Bioregions_sf_Subregions_level_casted$geometry[i], to = "POINT"))
}

table(nb_points_in_polygons)
sum(nb_points_in_polygons) # 3.2 M points (too many)

# Simplify shapes to reduce number of points
sf_use_s2(FALSE)
Bioregions_sf_Subregions_level_cleaned <- sf::st_simplify(x = Bioregions_sf_Subregions_level_casted, preserveTopology = T, dTolerance = 1)
Bioregions_sf_Subregions_level_cleaned <- st_make_valid(st_buffer(Bioregions_sf_Subregions_level_cleaned, 0))
Bioregions_sf_Subregions_level_cleaned <- st_crop(x = Bioregions_sf_Subregions_level_cleaned, y = st_bbox(Bioregions_sf_Subregions_level))
sf_use_s2(TRUE)

plot(Bioregions_sf_Subregions_level_cleaned[, "Subregion"])

# Check total nb of points
nb_points_in_polygons <- c()
for (i in 1:nrow(Bioregions_sf_Subregions_level_cleaned))
{
  nb_points_in_polygons[i] <- length(st_cast(Bioregions_sf_Subregions_level_cleaned$geometry[i], to = "POINT"))
}
sum(nb_points_in_polygons) # 9,816 points

## 4.2.2/ Cast into POINTS to run tessellation on all points at once #### 

Bioregions_sf_Subregions_level_cleaned_points <- st_cast(Bioregions_sf_Subregions_level_cleaned, to = "POINT")
nrow(Bioregions_sf_Subregions_level_cleaned_points) # 9,816 points

## Error in attribution of Subregions
# Correct by using intersections with initial polygons

sf_use_s2(FALSE)
Bioregions_sf_Subregions_level_cleaned_points <- Bioregions_sf_Subregions_level_cleaned_points %>%
  select(-Subregion) %>% # Remove Subregion
  st_join(y = Bioregions_sf_Subregions_level_cleaned[, "Subregion"],
          join = st_intersects, 
          left = TRUE)
sf_use_s2(TRUE)

nrow(Bioregions_sf_Subregions_level_cleaned_points) # 10,859 points because some points are duplicated at joint borders

table(is.na(Bioregions_sf_Subregions_level_cleaned_points$Subregion))
table(Bioregions_sf_Subregions_level_cleaned_points$Subregion)

## Plot border points per subregions to check for errors

pdf(file = "./outputs/Paleomaps/Subregions_border_points.pdf", width = 12, height = 6)

Subregions_list <- unique(Bioregions_sf_Subregions_level_cleaned_points$Subregion)
for (i in seq_along(Subregions_list))
{
  # i <- 1
  
  Subregion_i <- Subregions_list[[i]]
  cleaned_points_i <- Bioregions_sf_Subregions_level_cleaned_points %>% 
    filter(Subregion == Subregion_i)
  
  ggplot_Subregion_border_points_i <- ggplot(data = Bioregions_sf_Subregions_level) +
    
    geom_sf(fill = NA) +
    
    geom_sf(data = cleaned_points_i,
            col = "red", alpha = 0.5) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Subregion = ", Subregion_i)) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 8, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(ggplot_Subregion_border_points_i)
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion border points plotted for Subregion = ",Subregion_i," - n°", i, "/", length(Subregions_list),"\n"))
  }
  
}

dev.off()

# Save border points used to create the Voronoi cells
saveRDS(Bioregions_sf_Subregions_level_cleaned_points, file = "./outputs/Paleomaps/Bioregions_sf_Subregions_level_cleaned_points.rds")


## 4.2.3/ Run Derichlet tessellation with st_voronoi to obtain Voronoi cells ####

# # Convert to Mollweide
# Bioregions_sf_Subregions_level_cleaned_points_Mollweide <- st_transform(x = Bioregions_sf_Subregions_level_cleaned_points, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

# # Test without subregions crossing the anti-meridian
# test <- Bioregions_sf_Subregions_level_cleaned_points_Mollweide %>%
#   filter(!(Subregion %in% c("Far_East_Eastern_Palearctic", "Oceanian_Australasia", "Antarctica", "Zealandia_Australasia")))

## Union in a single MULTIPOINT object
# MULTIPOINT_for_Voronoi_Mollweide <- sf::st_union(Bioregions_sf_Subregions_level_cleaned_points_Mollweide)
MULTIPOINT_for_Voronoi <- sf::st_union(Bioregions_sf_Subregions_level_cleaned_points)

## Run Derichlet tessellation
sf_use_s2(FALSE)
# Voronoi_polygons_for_subplaques <- st_voronoi(test, dTolerance = 0, bOnlyEdges = FALSE)
Voronoi_polygons_for_subplaques <- st_voronoi(MULTIPOINT_for_Voronoi, dTolerance = 0, bOnlyEdges = FALSE)
# Voronoi_polygons_for_subplaques <- st_voronoi(MULTIPOINT_for_Voronoi_Mollweide, dTolerance = 0, bOnlyEdges = FALSE)
sf_use_s2(TRUE)

## Convert to POLYGON
Voronoi_polygons_for_subplaques <- Voronoi_polygons_for_subplaques %>%
  st_cast() %>% 
  st_as_sf() %>% 
  select(-area)


## Intersect Voronoi cells with initial points to attribute Subregions

sf_use_s2(FALSE)
Voronoi_polygons_for_subplaques <- Voronoi_polygons_for_subplaques %>%
  st_join(y = Bioregions_sf_Subregions_level_cleaned_points,
          join = st_intersects, 
          left = TRUE) %>%
  rename(geometry = x)
Voronoi_polygons_for_subplaques <- st_crop(x = Voronoi_polygons_for_subplaques, y = st_bbox(Bioregions_sf_Subregions_level))
Voronoi_polygons_for_subplaques <- st_make_valid(st_buffer(Voronoi_polygons_for_subplaques, 0))
sf_use_s2(TRUE)

nrow(Voronoi_polygons_for_subplaques)

View(Voronoi_polygons_for_subplaques)
plot(Voronoi_polygons_for_subplaques)


## Aggregate Voronoi cells per Subregions

sf_use_s2(FALSE)
Voronoi_polygons_for_subplaques <- Voronoi_polygons_for_subplaques %>%
  group_by(Subregion) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Voronoi_polygons_for_subplaques <- st_make_valid(st_buffer(Voronoi_polygons_for_subplaques, 0))
Voronoi_polygons_for_subplaques <- st_crop(x = Voronoi_polygons_for_subplaques, y = st_bbox(Bioregions_sf_Subregions_level))
sf_use_s2(TRUE)

nrow(Voronoi_polygons_for_subplaques)

plot(Voronoi_polygons_for_subplaques[, "Subregion"])


## Plot result

pdf(file = "./outputs/Paleomaps/Subregions_map_with_subplaques.pdf", width = 12, height = 6)

ggplot(data = Bioregions_sf_Subregions_level) +
  
  geom_sf(fill = NA) +
  
  geom_sf(data = Voronoi_polygons_for_subplaques,
          mapping = aes(fill = Subregion), alpha = 0.5) +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  # Adjust labels
  xlab("") +
  ylab("") +
  
  # Adjust title
  ggtitle(label = paste0("Subregions")) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 8, face = "bold"),
        legend.box.margin = margin(l = 5),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

dev.off()


## Save Voronoi polygons
#saveRDS(Voronoi_polygons_for_subplaques, file = "./outputs/Paleomaps/Voronoi_polygons_for_subplaques.rds")


## 4.2.4 / Correct internal borders using intersection with Subregions ####

# Rule: Subregion identity > Subplaque identity

## Load Subregions_sf
Bioregions_sf_Subregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Subregions_level.rds")

# Create a negative polygon when falling outside of Subregions
sf_use_s2(FALSE)
Bioregions_sf_Subregion_negative <- Bioregions_sf_Subregions_level %>% 
  nngeo::st_remove_holes() %>% 
  st_union()
Bioregions_sf_Subregion_negative <- st_difference(x = st_as_sfc(st_bbox(Bioregions_sf_Subregions_level)), y = Bioregions_sf_Subregion_negative) %>% 
  st_as_sf() %>% 
  mutate(Subregion = "NA") %>% 
  rename(geometry = x) %>%
  select(Subregion, geometry)
sf_use_s2(TRUE)

plot(Bioregions_sf_Subregion_negative[, "Subregion"])

# Add negative to Subregions_sf
Bioregions_sf_Subregions_level_with_NA <- rbind(Bioregions_sf_Subregions_level, Bioregions_sf_Subregion_negative)

# Intersect Subregions and subplaques, keeping subplaques with no matches
sf_use_s2(FALSE)
Voronoi_polygons_for_subplaques_corrected <- Voronoi_polygons_for_subplaques %>% 
  rename(Subregion_plaque = Subregion) %>% 
  st_intersection(y = Bioregions_sf_Subregions_level_with_NA)
sf_use_s2(TRUE)

# Apply Rule: Subregion identity > Subplaque identity
Voronoi_polygons_for_subplaques_corrected$Subregion_corrected <- Voronoi_polygons_for_subplaques_corrected$Subregion
Voronoi_polygons_for_subplaques_corrected$Subregion_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "NA"] <- Voronoi_polygons_for_subplaques_corrected$Subregion_plaque[Voronoi_polygons_for_subplaques_corrected$Subregion == "NA"]

# Aggregate with corrected Subregions
sf_use_s2(FALSE)
Voronoi_polygons_for_subplaques_corrected <- Voronoi_polygons_for_subplaques_corrected %>%
  group_by(Subregion_corrected) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup() %>%
  rename(Subregion = Subregion_corrected)
Voronoi_polygons_for_subplaques_corrected <- st_make_valid(st_buffer(Voronoi_polygons_for_subplaques_corrected, 0))
Voronoi_polygons_for_subplaques_corrected <- st_crop(x = Voronoi_polygons_for_subplaques_corrected, y = st_bbox(Bioregions_sf_Subregions_level))
sf_use_s2(TRUE)

view(Voronoi_polygons_for_subplaques_corrected)

plot(Voronoi_polygons_for_subplaques_corrected[, "Subregion_corrected"])

# Remove portion of the Mexican_Nearctic overlapping on the MesoAmerican_Neotropics
Mexican_Nearctic_plaque <- Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "Mexican_Nearctic", ]
Mesoamerican_Neotropics_plaque <- Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "Mesoamerican_Neotropics", ]
sf_use_s2(FALSE)
Mexican_Nearctic_plaque <- st_difference(Mexican_Nearctic_plaque, Mesoamerican_Neotropics_plaque$geometry)
sf_use_s2(TRUE)
Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "Mexican_Nearctic", ] <- Mexican_Nearctic_plaque

# Remove portion of the South_Europa_Western_Palearctic overlapping on the NE_Europa_Western_Palearctic
South_Europa_Western_Palearctic_plaque <- Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "South_Europa_Western_Palearctic", ]
NE_Europa_Western_Palearctic_plaque <- Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "NE_Europa_Western_Palearctic", ]
sf_use_s2(FALSE)
South_Europa_Western_Palearctic_plaque <- st_difference(South_Europa_Western_Palearctic_plaque, NE_Europa_Western_Palearctic_plaque$geometry)
sf_use_s2(TRUE)
Voronoi_polygons_for_subplaques_corrected[Voronoi_polygons_for_subplaques_corrected$Subregion == "South_Europa_Western_Palearctic", ] <- South_Europa_Western_Palearctic_plaque


## Plot corrected result

pdf(file = "./outputs/Paleomaps/Subregions_map_with_subplaques_corrected.pdf", width = 12, height = 6)

ggplot(data = Bioregions_sf_Subregions_level) +
  
  geom_sf(fill = NA) +
  
  geom_sf(data = Voronoi_polygons_for_subplaques_corrected,
          mapping = aes(fill = Subregion), alpha = 0.5) +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  # Adjust labels
  xlab("") +
  ylab("") +
  
  # Adjust title
  ggtitle(label = paste0("Subregions")) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 8, face = "bold"),
        legend.box.margin = margin(l = 5),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

dev.off()

## Save Voronoi polygons corrected
saveRDS(Voronoi_polygons_for_subplaques_corrected, file = "./outputs/Paleomaps/Voronoi_polygons_for_subplaques_corrected.rds")


## Alternative framework to get voronoï cells

# 1/ Get minimum convex polygons (MCP)
?st_convex_hull
# 2/ Extract points per subregions MCP
# 3/ Run Dirichlet tesselation to get Voronoï cells on all points
# 4/ Merge cells from the same subregions


### 4.3/ Project subplaques in the past with the API ####

# https://portal.gplates.org/service/d3_demo/?view=feature_collection

# ## Load Voronoi polygons
# Voronoi_polygons_for_subplaques_corrected <- readRDS(file = "./outputs/Paleomaps/Voronoi_polygons_for_subplaques_corrected.rds")

## Try with the Magalasy subplaque & the Australian subplaque

# plot(Voronoi_polygons_for_subplaques_corrected[, "Subregion"])
# 
# Magalasy_subplaque <- Voronoi_polygons_for_subplaques_corrected %>% 
#   filter(Subregion == "Magalasy_Afrotropics")
# Australian_subplaque <- Voronoi_polygons_for_subplaques_corrected %>% 
#   filter(Subregion == "Australian_Australasia")
# 
# # Remove holes to keep a single polygon
# Magalasy_subplaque <- nngeo::st_remove_holes(Magalasy_subplaque)
# Australian_subplaque <- nngeo::st_remove_holes(Australian_subplaque)
# 
# 
# Magalasy_subplaque <- keep_only_largest_polygon(Magalasy_subplaque)
# Australian_subplaque <- keep_only_largest_polygon(Australian_subplaque)
# 
# plot(Magalasy_subplaque[, "Subregion"])
# plot(Australian_subplaque[, "Subregion"])
# 
# # Simplify shape
# sf_use_s2(FALSE)
# Magalasy_subplaque <- st_simplify(Magalasy_subplaque, dTolerance = 2)
# sf_use_s2(TRUE)
# 
# plot(Magalasy_subplaque[, "Subregion"])
# 
# ## Round coordinates using conversion in sp object
# Magalasy_subplaque_rounded <- remove_consecutive_duplicated_coordinates_from_polygons(sf = Magalasy_subplaque, round_coordinates = T, digits = 1)
# 
# plot(Magalasy_subplaque[, "Subregion"])
# plot(Magalasy_subplaque_rounded[, "Subregion"])
# 
# Australian_subplaque_rounded <- remove_consecutive_duplicated_coordinates_from_polygons(sf = Australian_subplaque, round_coordinates = T, digits = 1)
# 
# plot(Australian_subplaque[, "Subregion"])
# plot(Australian_subplaque_rounded[, "Subregion"])
# 
# ## Merge plaques
# Magalasy_Australian_subplaque_rounded <- rbind(Magalasy_subplaque_rounded, Australian_subplaque_rounded)
# plot(Magalasy_Australian_subplaque_rounded[, "Subregion"])
# 
# ## Export GeoJSON
# 
# Magalasy_subplaque_GeoJSON <- geojsonsf::sf_geojson(sf = Magalasy_subplaque)
# Magalasy_subplaque_GeoJSON <- geojsonsf::sf_geojson(sf = Magalasy_subplaque_rounded)
# 
# Australian_subplaque_GeoJSON <- geojsonsf::sf_geojson(sf = Australian_subplaque_rounded)
# 
# Magalasy_Australian_subplaque_GeoJSON <- geojsonsf::sf_geojson(sf = Magalasy_Australian_subplaque_rounded)
# 
# 
# writeLines(text = Magalasy_subplaque_GeoJSON, con = "./outputs/Paleomaps/Magalasy_subplaque_GeoJSON.geojson", sep = "\\")
# writeLines(text = Australian_subplaque_GeoJSON, con = "./outputs/Paleomaps/Australian_subplaque_GeoJSON.geojson", sep = "\\")
# writeLines(text = Magalasy_Australian_subplaque_GeoJSON, con = "./outputs/Paleomaps/Magalasy_Australian_subplaque_GeoJSON.geojson", sep = "\\")
# 
# ## Get reconstruct shape from the API
# 
# age <- 30
# model <- "MULLER2022"
# 
# # URL for API request
# url <- base::paste0('https://gws.gplates.org/reconstruct/reconstruct_feature_collection/?feature_collection=', Magalasy_subplaque_GeoJSON, '&time=', age, '&model=', model)
# url <- base::paste0('https://gws.gplates.org/reconstruct/reconstruct_feature_collection/?feature_collection=', Australian_subplaque_GeoJSON, '&time=', age, '&model=', model)
# url <- base::paste0('https://gws.gplates.org/reconstruct/reconstruct_feature_collection/?feature_collection=', Magalasy_Australian_subplaque_GeoJSON, '&time=', age, '&model=', model)
# 
# #try to get the requested map of the model
# #print error message if map is not available
# #err: boolean if there was an error getting the map
# err <- FALSE
# shape <- tryCatch(
#   {
#     # Old version with rgdal
#     # rgdal::readOGR(url, verbose = FALSE)
#     
#     # New version with sf objects and geojsonsf
#     geojsonsf::geojson_sf(geojson = url)
#     
#   }, error = function(e) {
#     err <- TRUE
#     message(base::paste0("There is no map for ", age, " mya in ", model, " model. Please check the spelling, the age and the model you chose."))
#     stop()
#   }
# )
# 
# plot(shape)
# 
# 
# test_static <- geojsonsf::geojson_sf(geojson = "https://gws.gplates.org/reconstruct/static_polygons/?&time=0&model=MULLER2022")
# plot(test_static)


### 4.4/ Project subplaques in the past with the rgplates package ####

## Load Voronoi polygons
Voronoi_polygons_for_subplaques_corrected <- readRDS(file = "./outputs/Paleomaps/Voronoi_polygons_for_subplaques_corrected.rds")


GPlates_path <- "C:/Program Files/GPlates/GPlates 2.5.0/"

?rgplates::platemodel
?chronosphere::datasets
?chronosphere::fetch

## 4.4.1/ Register model features to be used locally by the GPlates Desktop App ####

# local_model <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/1000_0_rotfile.rot"),
#                                     features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

local_model_optimized <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/optimisation/1000_0_rotfile_MantleOptimised.rot"),
                                           features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

# local_model_no_net <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/optimisation/no_net_rotation_model.rot"),
#                                     features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

## 4.4.2/ Obtain reconstruct from polygons ####

# # Results = MULTILINESTRING without order to reconstruct polygons...
# Subplaques_reconstruct_120My <- rgplates::reconstruct(x = Voronoi_polygons_for_subplaques_corrected,
#                                                      age = 120,
#                                                      # model = local_model,
#                                                      model = local_model_optimized,
#                                                      # model = local_model_no_net,
#                                                      # from = 0,      # Age of the provided object (only for the online recunstruction)
#                                                      listout = T,   # Should the output be a list
#                                                      verbose = T,
#                                                      path.gplates = paste0(GPlates_path,"gplates.exe"),
#                                                      cleanup = T, # Should temp files be removed?
#                                                      dir = "./outputs/Paleomaps/Temp_maps/", # Directory for temp files
#                                                      gmeta = T)  # To keep metadata produced by GPlates in the sf output
# 
# 
# View(Subplaques_reconstruct_120My)
# plot(Voronoi_polygons_for_subplaques_corrected[, "Subregion"])
# plot(Subplaques_reconstruct_120My[, "Subregion"])

## 4.4.3/ Keep only the largest polygon per subplaques ####

## Special case of Subregions crossing the antimeridian
# Need to identity the East and West polygons
# Extract them as independent POLYGON rather than MULTIPOLYGON

Subregions_to_split_by_antimeridian <- c("American_Nearctic", "Far_East_Eastern_Palearctic", "Oceanian_Australasia", "Zealandia_Australasia")

# Create East/West hemisphere polygons
West_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = -180, xmax = 0, ymax = 90, ymin = -90), crs = st_crs(4326))))
East_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = 0, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))))
EW_polygons <- rbind(East_hemisphere, West_hemisphere) %>% 
  rename(geometry = x) %>% 
  mutate(Hemisphere = c("HEast", "HWest"))

Voronoi_polygons_for_subplaques_split <- Voronoi_polygons_for_subplaques_corrected
for (i in seq_along(Subregions_to_split_by_antimeridian))
{
  # i <- 1
  
  # Extract subplaque
  Subregion_i <- Subregions_to_split_by_antimeridian[i]
  Subplaque_i <- Voronoi_polygons_for_subplaques_corrected %>% 
    filter(Subregion == Subregion_i)
  
  # Cast in POLYGONS
  Subplaques_i <- Subplaque_i %>% 
    st_cast(to = "POLYGON")
  
  # Sort by East/West Hemispheres
  sf_use_s2(FALSE)
  Subplaques_split_i <- st_intersection(x = Subplaques_i, y = EW_polygons)
  sf_use_s2(TRUE)
  
  # plot(Subplaques_split_i[, "Hemisphere"])
  
  # Extract the largest polygon in each Hemisphere
  Subplaques_split_per_H_i <- data.frame()
  for (j in c("HEast", "HWest"))
  {
    # j <- "East"
    
    Subplaques_split_j <- Subplaques_split_i %>% 
      filter(Hemisphere == j) %>% 
      st_cast(to = "MULTIPOLYGON")
    
    # plot(Subplaques_split_j[, "Hemisphere"])
    
    sf_use_s2(FALSE)
    Subplaque_split_j <- keep_only_largest_polygon(Subplaques_split_j) 
    sf_use_s2(TRUE)
    
    # plot(Subplaque_split_j[, "Hemisphere"])
    
    # Store result
    Subplaques_split_per_H_i <- rbind(Subplaques_split_per_H_i, Subplaque_split_j)
  }
  
  # Rename output
  Subplaques_split_per_H_i <- Subplaques_split_per_H_i %>% 
    mutate(Subregion = paste0(Subregion, "_",Hemisphere)) %>% 
    select(-Hemisphere)

  # Store result
  Voronoi_polygons_for_subplaques_split <- Voronoi_polygons_for_subplaques_split %>% 
    filter(Subregion != Subregion_i) %>%
    rbind(Subplaques_split_per_H_i)
}

View(Voronoi_polygons_for_subplaques_split)

## Normal case = select only the biggest POLYGON

Voronoi_polygons_for_subplaques_single <- Voronoi_polygons_for_subplaques_split
for (i in 1:nrow(Voronoi_polygons_for_subplaques_split))
{
  # i <- 1
  
  Subplaque_i <- Voronoi_polygons_for_subplaques_split[i, ]
  
  sf_use_s2(FALSE)
  Subplaque_i <- keep_only_largest_polygon(Subplaque_i) 
  sf_use_s2(TRUE)
  
  # Store result
  Voronoi_polygons_for_subplaques_single[i, ] <- Subplaque_i
}

plot(Voronoi_polygons_for_subplaques_single[, "Subregion"])

# Save Subplaques with unique split polygons
saveRDS(object = Voronoi_polygons_for_subplaques_single, file = "./outputs/Paleomaps/Voronoi_polygons_for_subplaques_single.rds")

## 4.4.4/ Cast to POINTS to be able to reconstruct polygons ####

# Cast to POINTS
Subplaques_list <- Voronoi_polygons_for_subplaques_single$Subregion
Voronoi_points_for_subplaques_single <- data.frame()
for(i in seq_along(Subplaques_list))
{
  # i <- 1
  
  # Extract subplaque
  Subplaque_i <- Subplaques_list[i]
  
  # Cast the polygon in points
  Subplaque_points_i <- Voronoi_polygons_for_subplaques_single %>% 
    filter(Subregion == Subplaque_i) %>% 
    st_cast(to = "POINT") %>% 
    group_by(Subregion) %>%
    mutate(seq_ID = row_number()) %>% 
    ungroup()

  # Store output
  Voronoi_points_for_subplaques_single <- rbind(Voronoi_points_for_subplaques_single, Subplaque_points_i)
}

table(Voronoi_points_for_subplaques_single$Subregion)

plot(Voronoi_polygons_for_subplaques_single[Voronoi_polygons_for_subplaques_single$Subregion == "North_African_Western_Palearctic", "Subregion"])
plot(Voronoi_points_for_subplaques_single[Voronoi_points_for_subplaques_single$Subregion == "North_African_Western_Palearctic", "seq_ID"])

##  4.4.5/ Round and cleaned redundant coordinates to save time ####

Subplaques_list <- Voronoi_polygons_for_subplaques_single$Subregion
Voronoi_points_for_subplaques_single_rounded <- data.frame()
for (i in seq_along(Subplaques_list))
{
  # i <- 1
  
  # Extract subplaque
  Subplaque_i <- Subplaques_list[i]
  
  # Extract points
  Subplaque_points_i <- Voronoi_points_for_subplaques_single %>% 
    filter(Subregion == Subplaque_i)
  
  # Round coordinate and remove consecutive duplicates
  Subplaque_points_i <- remove_consecutive_duplicated_coordinates_from_points(sf = Subplaque_points_i, round_coordinates = T, digits = 1)
  
  # Store results
  Voronoi_points_for_subplaques_single_rounded <- rbind(Voronoi_points_for_subplaques_single_rounded, Subplaque_points_i)
}

# Register seq_ID of points within polygons
Voronoi_points_for_subplaques_single_rounded <- Voronoi_points_for_subplaques_single_rounded %>% 
  group_by(Subregion) %>%
  mutate(seq_ID = row_number()) %>% 
  ungroup()

table(Voronoi_points_for_subplaques_single_rounded$Subregion)
table(Voronoi_points_for_subplaques_single_rounded$seq_ID)

plot(Voronoi_polygons_for_subplaques_single[Voronoi_polygons_for_subplaques_single$Subregion == "North_African_Western_Palearctic", "Subregion"])
plot(Voronoi_points_for_subplaques_single_rounded[Voronoi_points_for_subplaques_single_rounded$Subregion == "North_African_Western_Palearctic", "seq_ID"])

# Save Subplaques border points to provide as input to GPLates
saveRDS(Voronoi_points_for_subplaques_single_rounded, file = "./outputs/Paleomaps/Voronoi_points_for_subplaques_single_rounded.rds")


## 4.4.6/ Obtain reconstruct from POINTS to preserve order and reconstruct polygons ####

# POINTS allows to preserve order and reconstruct polygons

## Set ages to recover

ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Get the reconstructed points
Subplaques_points_reconstruct_all_ages <- rgplates::reconstruct(x = Voronoi_points_for_subplaques_single_rounded,
                                                            age = ages_list,
                                                            # model = local_model,
                                                            model = local_model_optimized,
                                                            # model = local_model_no_net,
                                                            # from = 0,      # Age of the provided object (only for the online recunstruction)
                                                            listout = T,   # Should the output be a list
                                                            verbose = T,
                                                            path.gplates = paste0(GPlates_path,"gplates.exe"),
                                                            cleanup = T, # Should temp files be removed?
                                                            dir = "./outputs/Paleomaps/Temp_maps/", # Directory for temp files
                                                            gmeta = T)  # To keep metadata produced by GPlates in the sf output
names(Subplaques_points_reconstruct_all_ages) <- paste0("Subplaques_",model,"_",ages,"My")

View(Subplaques_points_reconstruct_all_ages[[1]])
plot(Voronoi_points_for_subplaques_single_rounded[, "Subregion"])
plot(Subplaques_points_reconstruct_all_ages[[1]][, "Subregion"])
plot(Subplaques_points_reconstruct_all_ages[[2]][, "Subregion"])
plot(Subplaques_points_reconstruct_all_ages[[3]][, "Subregion"])

## Loop per age to reorder points as initial order and remove subregions with less than 4 points remaining

for (i in seq_along(Subplaques_points_reconstruct_all_ages))
{
  # Extract subplaques_sf for age i
  Subplaques_points_reconstruct_age_i <- Subplaques_points_reconstruct_all_ages[[i]] 
  
  # Reorder points as initial order
  Subplaques_points_reconstruct_age_i <- Subplaques_points_reconstruct_age_i %>% 
    arrange(Subregion, seq_ID)
  
  # Remove subregions with less than 4 points remaining
  table(Subplaques_points_reconstruct_age_i$Subregion)
  subplaques_to_remove <- names(table(Subplaques_points_reconstruct_age_i$Subregion))[table(Subplaques_points_reconstruct_age_i$Subregion) < 4]
  Subplaques_points_reconstruct_age_i <- Subplaques_points_reconstruct_age_i %>% 
    filter(!(Subregion %in% subplaques_to_remove))
  
  # Store updated file
  Subplaques_points_reconstruct_all_ages[[i]] <- Subplaques_points_reconstruct_age_i
}

# Save the raw reconstruct subplaques per ages
saveRDS(Subplaques_points_reconstruct_all_ages, file = "./outputs/Paleomaps/Subplaques_points_reconstruct_all_ages.rds")


### 4.5/ Clean reconstructed point borders of subplaques ####

## 4.5.1/ Split subplaques crossing the antimeridian ####

# Do it only for Subregion suspected to cross the antimeridian, and NOT the meridian !
antimeridian_subplaques <- c("American_Nearctic_HEast", "American_Nearctic_HWest", "East_Wallacea_Australasia", "Far_East_Eastern_Palearctic_HEast", "Oceanian_Australasia_HEast", "Oceanian_Australasia_HWest", "Zealandia_Australasia_HEast", "Zealandia_Australasia_HWest")

# Loop per ages
Subplaques_points_reconstruct_all_ages_split <- Subplaques_points_reconstruct_all_ages
for (i in seq_along(Subplaques_points_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subplaques_points_reconstruct_all_ages)[i]
    
  # Extract subplaques_sf for age i
  Subplaques_points_reconstruct_age_i <- Subplaques_points_reconstruct_all_ages[[i]] 
  
  # Extract only subplaques crossing the antimeridian
  Subplaques_points_reconstruct_age_i_antimeridian <- Subplaques_points_reconstruct_age_i %>% 
    filter(Subregion %in% antimeridian_subplaques)
  table(Subplaques_points_reconstruct_age_i_antimeridian$Subregion)
  
  # Split plaques
  Subplaques_points_reconstruct_age_i_antimeridian <- suppressMessages(suppressWarnings(split_by_anti_meridian(sf_points = Subplaques_points_reconstruct_age_i_antimeridian,
                                                                                        group_var = "Subregion")))
  # Store them back
  Subplaques_points_reconstruct_age_i_split <- Subplaques_points_reconstruct_age_i %>% 
    filter(!(Subregion %in% antimeridian_subplaques)) %>% 
    rbind(Subplaques_points_reconstruct_age_i_antimeridian)
  
  plot(Subplaques_points_reconstruct_age_i_split[, "Subregion"], main = main_i)
  
  # Store updated file
  Subplaques_points_reconstruct_all_ages_split[[i]] <- Subplaques_points_reconstruct_age_i_split
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subplaques split for age = ",age_i," My - n°", i, "/", length(Subplaques_points_reconstruct_all_ages),"\n"))
  }
}  
  
# Save the split reconstruct subplaques per ages
saveRDS(Subplaques_points_reconstruct_all_ages_split, file = "./outputs/Paleomaps/Subplaques_points_reconstruct_all_ages_split.rds")


## 4.5.2/ Filter distance outlier points before reconstruction ####

## Filter by anomaly in distance vs. initial position  (only if not split by anti-meridian)

# Select subplaques to clean for distance outliers
Subplaques_list <- unique(Subplaques_points_reconstruct_120My_split$Subregion)
Suplaques_split_list <- Subplaques_list[str_detect(string = Subplaques_list, pattern = "_East$|_West$")]
Subplaques_no_split_list <- setdiff(Subplaques_list, Suplaques_split_list)

Subplaques_needing_distance_outliers_cleaning <- c("Arabic_Peninsula_Western_Palearctic", "Indian_subcontinent_Indomalaya")  

# Loop per ages
Subplaques_points_reconstruct_all_ages_no_dist_outliers <- Subplaques_points_reconstruct_all_ages_split
for (i in seq_along(Subplaques_points_reconstruct_all_ages_split))
{
  # i <- 3
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subplaques_points_reconstruct_all_ages_split)[i]
  
  # Extract subplaques_sf for age i
  Subplaques_points_reconstruct_age_i_split <- Subplaques_points_reconstruct_all_ages_split[[i]] 
  
  # Initiate new sf_object
  Subplaques_points_reconstruct_age_i_no_dist_outliers <- data.frame()
  
  # Loop per selected subplaques
  for (j in seq_along(Subplaques_needing_distance_outliers_cleaning))
  {
    # j <- 1
    
    # Extract subplaque
    Subplaque_j <- Subplaques_needing_distance_outliers_cleaning[j]
    
    # Extract points from reconstruct time
    Subplaque_points_j <- Subplaques_points_reconstruct_age_i_split %>% 
      filter(Subregion == Subplaque_j)
    
    # Extract points from initial time
    sf_init_j <- Voronoi_points_for_subplaques_single_rounded[Voronoi_points_for_subplaques_single_rounded$Subregion == Subplaque_j, ]
    
    if (nrow(Subplaque_points_j) > 4)
    {
      # Remove outliers
      Subplaque_points_j <- identify_outliers_by_distance_to_initial_position(sf = Subplaque_points_j, 
                                                                              sf_init = sf_init_j,
                                                                              quantile_dist_max = 0.99,
                                                                              quantile_dist_min = 0.20,
                                                                              remove = T, plot = T,
                                                                              main = "Subregion")
      # Remove outlier columns
      Subplaque_points_j <- Subplaque_points_j %>% 
        select(-max_outliers, -min_outliers)
    }
    
    # Store results for age i
    Subplaques_points_reconstruct_age_i_no_dist_outliers <- rbind(Subplaques_points_reconstruct_age_i_no_dist_outliers, Subplaque_points_j)
  }
  
  # Store cleaned subplaques back
  Subplaques_points_reconstruct_age_i_cleaned <- Subplaques_points_reconstruct_age_i_split %>% 
    filter(!(Subregion %in% Subplaques_needing_distance_outliers_cleaning)) %>% 
    rbind(Subplaques_points_reconstruct_age_i_no_dist_outliers)
  
  # plot(Subplaques_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
  
  # Store updated file
  Subplaques_points_reconstruct_all_ages_no_dist_outliers[[i]] <- Subplaques_points_reconstruct_age_i_cleaned
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subplaques cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subplaques_points_reconstruct_all_ages_split),"\n"))
  }
}  

# Save the reconstruct subplaques per ages cleaned from distance outliers
saveRDS(Subplaques_points_reconstruct_all_ages_no_dist_outliers, file = "./outputs/Paleomaps/Subplaques_points_reconstruct_all_ages_no_dist_outliers.rds")


## 4.5.3/ Filter sequence outlier points before reconstruction ####

## Filter by sequence anomaly as points having abnormal distance to their neighboring points in the sequence

# Select subplaques to clean for sequence outliers
Subplaques_list <- unique(Subplaques_points_reconstruct_all_ages_no_dist_outliers[[1]]$Subregion)
Subplaques_needing_sequence_outliers_cleaning <- c("American_Nearctic_HWest_West", "Continental_SE_Asia_Indomalaya", "Far_East_Eastern_Palearctic_HEast_East") 

# Loop per ages
Subplaques_points_reconstruct_all_ages_cleaned <- Subplaques_points_reconstruct_all_ages_no_dist_outliers
for (i in seq_along(Subplaques_points_reconstruct_all_ages_no_dist_outliers))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subplaques_points_reconstruct_all_ages_no_dist_outliers)[i]
  
  # Extract subplaques_sf for age i
  Subplaques_points_reconstruct_age_i_no_dist_outliers <- Subplaques_points_reconstruct_all_ages_no_dist_outliers[[i]] 
  
  # Initiate new sf_object
  Subplaques_points_reconstruct_age_i_no_sequence_outliers <- data.frame()
  
  # Loop per subplaques
  for (j in seq_along(Subplaques_needing_sequence_outliers_cleaning))
  {
    # j <- 1
    
    # Extract subplaque
    Subplaque_j <- Subplaques_needing_sequence_outliers_cleaning[j]
    
    # Extract points
    Subplaque_points_j <- Subplaques_points_reconstruct_age_i_no_dist_outliers %>% 
      filter(Subregion == Subplaque_j)
    
    if (nrow(Subplaque_points_j) > 4)
    {
      # Remove outliers
      Subplaque_points_j <- identify_outliers_from_sequence(sf = Subplaque_points_j, 
                                                            quantile_neighbors = 0.01, quantile_dist = 0.985,
                                                            remove = F, plot = T,
                                                            main = "Subregion")
      # Remove outlier columns
      Subplaque_points_j <- Subplaque_points_j %>% 
        select(-outliers)
    }
    
    # Store results
    Subplaques_points_reconstruct_age_i_no_sequence_outliers <- rbind(Subplaques_points_reconstruct_age_i_no_sequence_outliers, Subplaque_points_j)
  }
  
  # Store cleaned subplaques back
  Subplaques_points_reconstruct_age_i_cleaned <- Subplaques_points_reconstruct_age_i_no_dist_outliers %>% 
    filter(!(Subregion %in% Subplaques_needing_sequence_outliers_cleaning)) %>% 
    rbind(Subplaques_points_reconstruct_age_i_no_sequence_outliers)
  
  # plot(Subplaques_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
  
  # Store updated file
  Subplaques_points_reconstruct_all_ages_cleaned[[i]] <- Subplaques_points_reconstruct_age_i_cleaned
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subplaques cleaned from sequence outliers for age = ",age_i," My - n°", i, "/", length(Subplaques_points_reconstruct_all_ages_no_dist_outliers),"\n"))
  }
}

# table(Subplaques_points_reconstruct_all_ages_no_dist_outliers[[1]]$Subregion)
# table(Subplaques_points_reconstruct_all_ages_cleaned[[1]]$Subregion)

# Save the reconstruct subplaques per ages cleaned from sequence outliers
saveRDS(Subplaques_points_reconstruct_all_ages_cleaned, file = "./outputs/Paleomaps/Subplaques_points_reconstruct_all_ages_cleaned.rds")


## 4.5.4/ Reconstruct polygons from border points ####

# Loop per ages
Subplaques_polygons_reconstruct_all_ages <- Subplaques_points_reconstruct_all_ages_cleaned
for (i in seq_along(Subplaques_polygons_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subplaques_points_reconstruct_all_ages_cleaned)[i]
  
  # Extract subplaques_sf for age i
  Subplaques_points_reconstruct_age_i_cleaned <- Subplaques_points_reconstruct_all_ages_cleaned[[i]] 
  
  # Reconstruct polygons using ordered points
  Subplaques_polygons_reconstruct_age_i <- reconstruct_initial_polygon_with_sfheaders(sf = Subplaques_points_reconstruct_age_i_cleaned)
  # Remove holes
  Subplaques_polygons_reconstruct_age_i <- nngeo::st_remove_holes(Subplaques_polygons_reconstruct_age_i)
  
  # Plot results
  # plot(Subplaques_polygons_reconstruct_age_i[, "Subregion"])
  
  # Store updated file
  Subplaques_polygons_reconstruct_all_ages[[i]] <- Subplaques_polygons_reconstruct_age_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subplaques polygons reconstructed for age = ",age_i," My - n°", i, "/", length(Subplaques_points_reconstruct_all_ages_cleaned),"\n"))
  }
}

# Save the reconstructed subplaque polygons per ages
saveRDS(Subplaques_polygons_reconstruct_all_ages, file = "./outputs/Paleomaps/Subplaques_polygons_reconstruct_all_ages.rds")


## Add to Paleomaps_sf !!


# Alternative: Use minimum convex hull on points?


### 4.6/ Plot reconstructed subplaques ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")


pdf(file = paste0("./outputs/Paleomaps/Plates_with_Subregion_borders_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Subplaques_polygons_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[j]]
  # Plate_boundaries_i <- Paleomaps_sf$Plate_boundaries[[j]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[j]]
  Subplaques_polygons_reconstruct_i <- Subplaques_polygons_reconstruct_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Plate_polygons_i$Bioregion)]
  
  # Order bioregion factors
  Plate_polygons_i$Bioregion <- factor(x = Plate_polygons_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  
  # Plot
  Subplates_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Subplaques_polygons_reconstruct_i,
            mapping = aes(fill = Subregion), alpha = 0.9,
            show.legend = T) +
    
    # geom_sf(data = Plate_polygons_i,
    #         mapping = aes(fill = Bioregion), alpha = 0.2, col = "grey50") +
    
    # # WGS84 labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec, y = latitude_dec), col = "black") +
    
    # # Mollweide labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec_Mollweide, y = latitude_dec_Mollweide),
    #            col = "#00000066", alpha = 0.5) +
    
    
    # geom_sf(data = Plate_boundaries_i,
    #         col = "purple") +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # # Adjust fill legend
    # scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Tectonic plates in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Subplates_plot_i)
}
dev.off()


### 4.7/ Merge subplaques per Bioregions ####

Bioregions_list <- c("Antarctica", "Afrotropics", "Australasia", "Nearctic", "Neotropics", "Indomalaya", "Eastern_Palearctic", "Western_Palearctic")

Bioregions_polygons_reconstruct_all_ages <- Subplaques_polygons_reconstruct_all_ages
# Loop per ages
for (i in seq_along(Subplaques_polygons_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subplaques_polygons_reconstruct_all_ages)[i]
  
  # Extract subplaques_sf for age i
  Subplaques_polygons_reconstruct_age_i <- Subplaques_polygons_reconstruct_all_ages[[i]] 
  Subplaques_polygons_reconstruct_age_i$Bioregion <- NA
  
  # Extract subplaques/subregion
  all_subregions_i <- unique(Subplaques_polygons_reconstruct_age_i$Subregion)
  
  # Loop per Bioregions
  for (j in seq_along(Bioregions_list))
  {
    # j <- 3
    
    # Extract Bioregion
    Bioregion_j <- Bioregions_list[j]
    
    # Detect associated subregions
    Subregions_j <- all_subregions_i[str_detect(string = all_subregions_i, pattern = Bioregion_j)]
    Subplaques_polygons_reconstruct_age_i$Bioregion[Subplaques_polygons_reconstruct_age_i$Subregion %in% Subregions_j] <- Bioregion_j
    
  }
  
  # Merge all associated subregions together
  sf_use_s2(FALSE)
  Subplaques_polygons_reconstruct_age_i <- st_make_valid(Subplaques_polygons_reconstruct_age_i)
  Bioregions_polygons_reconstruct_age_i <- Subplaques_polygons_reconstruct_age_i %>% 
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() 
  sf_use_s2(TRUE)
  
  # Plot results
  plot(Bioregions_polygons_reconstruct_age_i[, "Bioregion"], main = main_i)
  
  # Store updated file
  Bioregions_polygons_reconstruct_all_ages[[i]] <- Bioregions_polygons_reconstruct_age_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion polygons aggregated for age = ",age_i," My - n°", i, "/", length(Subplaques_polygons_reconstruct_all_ages),"\n"))
  }
}

# Save the reconstructed Bioregion polygons per ages
saveRDS(Bioregions_polygons_reconstruct_all_ages, file = "./outputs/Paleomaps/Bioregions_polygons_reconstruct_all_ages.rds")


### 4.8/ Interpolate membership for empty areas using the closest polygon ####

# Voronoi cells around the Bioplaques to get reconstructed Bioplaques covering the whole Earth

## Loop per age
Voronoi_polygons_for_Bioplaques_all_ages <- Bioregions_polygons_reconstruct_all_ages
for (i in seq_along(Bioregions_polygons_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Bioregions_polygons_reconstruct_all_ages)[i]
  
  # Extract bioregion reconstructed polygons
  Bioregions_polygons_reconstruct_age_i <- Bioregions_polygons_reconstruct_all_ages[[i]]
  
  ## 4.8.1/ Cast into POINTS to run tessellation on all points at once #### 
  
  Bioregions_polygons_reconstruct_age_i_points <- Bioregions_polygons_reconstruct_age_i %>% 
    st_collection_extract(type = "POLYGON") %>%
    st_cast(to = "POINT")
  
  # table(Bioregions_polygons_reconstruct_age_i_points$Bioregion)
  
  ## 4.8.2/ Run Derichlet tessellation with st_voronoi to obtain Voronoi cells ####
  
  ## Union in a single MULTIPOINT object
  MULTIPOINT_for_Voronoi_age_i <- sf::st_union(Bioregions_polygons_reconstruct_age_i_points)
  
  ## Run Derichlet tessellation
  # sf_use_s2(FALSE)
  MULTIPOINT_for_Voronoi_age_i <- st_make_valid(MULTIPOINT_for_Voronoi_age_i)
  Voronoi_polygons_for_Bioplaques_age_i <- st_voronoi(MULTIPOINT_for_Voronoi_age_i, dTolerance = 0, bOnlyEdges = FALSE)
  # sf_use_s2(TRUE)
  
  ## Convert to POLYGON
  Voronoi_polygons_for_Bioplaques_age_i <- Voronoi_polygons_for_Bioplaques_age_i %>%
    st_cast() %>% 
    st_as_sf()
  
  ## Intersect Voronoi cells with initial points to attribute Bioregions
  
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioplaques_age_i <- Voronoi_polygons_for_Bioplaques_age_i %>%
    st_join(y = Bioregions_polygons_reconstruct_age_i_points,
            join = st_intersects, 
            left = TRUE) %>%
    rename(geometry = x)
  Voronoi_polygons_for_Bioplaques_age_i <- st_crop(x = Voronoi_polygons_for_Bioplaques_age_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  Voronoi_polygons_for_Bioplaques_age_i <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioplaques_age_i, 0))
  sf_use_s2(TRUE)
  
  # nrow(Voronoi_polygons_for_Bioplaques_age_i)
  
  # View(Voronoi_polygons_for_Bioplaques_age_i)
  # plot(Voronoi_polygons_for_Bioplaques_age_i)
  
  
  ## Aggregate Voronoi cells per Bioregions
  
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioplaques_age_i <- Voronoi_polygons_for_Bioplaques_age_i %>%
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()
  Voronoi_polygons_for_Bioplaques_age_i <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioplaques_age_i, 0))
  Voronoi_polygons_for_Bioplaques_age_i <- st_crop(x = Voronoi_polygons_for_Bioplaques_age_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  # nrow(Voronoi_polygons_for_Bioplaques_age_i)
   
  # plot(Voronoi_polygons_for_Bioplaques_age_i[, "Bioregion"])
  
  ## 4.8.3/ Correct internal borders using intersection with Bioregions ####
  
  # Rule: Bioregion identity > Bioplaque identity
  
  # Create a negative polygon when falling outside of Bioregions
  sf_use_s2(FALSE)
  Bioregions_sf_Bioregion_negative_age_i <- Bioregions_polygons_reconstruct_age_i %>% 
    st_collection_extract(type = "POLYGON") %>%
    nngeo::st_remove_holes() %>% 
    st_union()
  Bioregions_sf_Bioregion_negative_age_i <- st_difference(x = st_as_sfc(st_bbox(Bioregions_polygons_reconstruct_age_i)), y = Bioregions_sf_Bioregion_negative_age_i) %>% 
    st_as_sf() %>% 
    mutate(Bioregion = "NA") %>% 
    rename(geometry = x) %>%
    select(Bioregion, geometry)
  sf_use_s2(TRUE)
  
  # plot(Bioregions_sf_Bioregion_negative_age_i[, "Bioregion"])
  
  # Add negative to Bioregions_sf
  Bioregions_polygons_reconstruct_age_i_with_NA <- rbind(Bioregions_polygons_reconstruct_age_i, Bioregions_sf_Bioregion_negative_age_i)
  
  # Intersect Bioregions and Bioplaques, keeping Bioplaques with no matches
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioplaques_age_i_corrected <- Voronoi_polygons_for_Bioplaques_age_i %>% 
    rename(Bioregion_plaque = Bioregion) %>% 
    st_intersection(y = Bioregions_polygons_reconstruct_age_i_with_NA)
  sf_use_s2(TRUE)
  
  # Apply Rule: Bioregion identity > Bioplaque identity
  Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion_corrected <- Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion
  Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion_corrected[Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion == "NA"] <- Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion_plaque[Voronoi_polygons_for_Bioplaques_age_i_corrected$Bioregion == "NA"]
  
  # Aggregate with corrected Bioregions
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioplaques_age_i_corrected <- Voronoi_polygons_for_Bioplaques_age_i_corrected %>%
    group_by(Bioregion_corrected) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() %>%
    rename(Bioregion = Bioregion_corrected)
  Voronoi_polygons_for_Bioplaques_age_i_corrected <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioplaques_age_i_corrected, 0))
  Voronoi_polygons_for_Bioplaques_age_i_corrected <- st_crop(x = Voronoi_polygons_for_Bioplaques_age_i_corrected, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  ## 4.8.4/ Correct other overlap by applying a priority rule ####
  
  # Define priority sequence
  priority_sequence <- c("Neotropics", "Nearctic", "Autralasia", "Afrotropics", "Western_Palearctic", "Indomalaya", "Eastern_Palearctic")
  
  Voronoi_polygons_for_Bioplaques_age_i_masked <- Voronoi_polygons_for_Bioplaques_age_i_corrected
  Voronoi_polygons_for_Bioplaques_age_i_mask <- Voronoi_polygons_for_Bioplaques_age_i_corrected
  for (j in seq_along(priority_sequence))
  {
    # j <- 1
    Bioregion_j <- rev(priority_sequence)[j]
    
    # Extract bioplaque
    Bioregion_full_j <- Voronoi_polygons_for_Bioplaques_age_i_mask %>% 
      filter(Bioregion == Bioregion_j)
    
    # Merge all other bioregions to create mask
    sf_use_s2(FALSE)
    mask_j <- Voronoi_polygons_for_Bioplaques_age_i_mask %>% 
      filter(Bioregion != Bioregion_j) %>% 
      st_union()
    sf_use_s2(TRUE)
    
    # Apply mask on Bioregion polygon
    sf_use_s2(FALSE)
    Bioregion_masked_j <- Bioregion_full_j %>%
      st_difference(y = mask_j) %>%
      st_as_sf() %>%
      mutate(Bioregion = Bioregion_j)
    sf_use_s2(TRUE)
    
    # plot(Bioregion_full_j, col = "red")
    # plot(Bioregion_masked_j, col = "limegreen")
    
    # Remove focal bioregion from the mask for next iteration
    Voronoi_polygons_for_Bioplaques_age_i_mask <- Voronoi_polygons_for_Bioplaques_age_i_mask %>% 
      filter(Bioregion != Bioregion_j)
    
    # Store updated mask version
    Voronoi_polygons_for_Bioplaques_age_i_masked[Voronoi_polygons_for_Bioplaques_age_i_masked$Bioregion == Bioregion_j, ] <- Bioregion_masked_j
  }
  
  # plot(Voronoi_polygons_for_Bioplaques_age_i_corrected[, "Bioregion"], main = main_i)
  # plot(Voronoi_polygons_for_Bioplaques_age_i_masked[, "Bioregion"], main = main_i)
   
  # view(Voronoi_polygons_for_Bioplaques_age_i_corrected)
  
  plot(Voronoi_polygons_for_Bioplaques_age_i_masked[, "Bioregion"], main = main_i)
  
  # Store results
  Voronoi_polygons_for_Bioplaques_all_ages[[i]] <- Voronoi_polygons_for_Bioplaques_age_i_masked
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion polygons extended for age = ",age_i," My - n°", i, "/", length(Bioregions_polygons_reconstruct_all_ages),"\n"))
  }
}
  
## Save Bioregions Voronoi polygons = Bioplaques for all ages
saveRDS(Voronoi_polygons_for_Bioplaques_all_ages, file = "./outputs/Paleomaps/Voronoi_polygons_for_Bioplaques_all_ages.rds")



### 4.9/ Plot bioregion plates ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
# bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
bioregion_names_reconstruct <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern_Palearctic", "Western_Palearctic")
names(colors_list_for_areas) <- bioregion_names_reconstruct
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
# bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names_reconstruct, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA


pdf(file = paste0("./outputs/Paleomaps/Plates_with_reconstructed_extended_Bioregion_borders_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Voronoi_polygons_for_Bioplaques_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[j]]
  # Plate_boundaries_i <- Paleomaps_sf$Plate_boundaries[[j]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[j]]
  Voronoi_polygons_for_Bioplaques_i <- Voronoi_polygons_for_Bioplaques_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Plate_polygons_i$Bioregion)]
  bioregion_names_present_reconstruct <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Bioregions_polygons_reconstruct_i$Bioregion)]
  
  # Order bioregion factors
  Plate_polygons_i$Bioregion <- factor(x = Plate_polygons_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  Voronoi_polygons_for_Bioplaques_i$Bioregion <- factor(x = Voronoi_polygons_for_Bioplaques_i$Bioregion, levels = bioregion_names_present_reconstruct, labels = bioregion_names_present_reconstruct)
  
  # Plot
  Bioregions_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Voronoi_polygons_for_Bioplaques_i,
            mapping = aes(fill = Bioregion), alpha = 0.5,
            show.legend = T) +
    
    # geom_sf(data = Plate_polygons_i,
    #         mapping = aes(fill = Bioregion), alpha = 0.2, col = "grey50") +
    
    # # WGS84 labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec, y = latitude_dec), col = "black") +
    
    # # Mollweide labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec_Mollweide, y = latitude_dec_Mollweide),
    #            col = "#00000066", alpha = 0.5) +
  
  
  # geom_sf(data = Plate_boundaries_i,
  #         col = "purple") +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
    
    # Adjust fill legend
    # scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present_reconstruct), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present_reconstruct]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Tectonic plates in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Bioregions_plot_i)
}
dev.off()


##### 5/ Intersect reconstructed plaques and plate identities ####

# Rules: 
# Plate identities > Reconstructed plaques
# If NA for Plate identities, use Reconstructed plaques


## Load Bioregions Voronoi polygons = Reconstructed bioplaques for all ages
Voronoi_polygons_for_Bioplaques_all_ages <- readRDS(file = "./outputs/Paleomaps/Voronoi_polygons_for_Bioplaques_all_ages.rds")

## Load paleomaps in sf format with Bioregion membership
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")


## 5.1/ Intersect plates with Bioregion membership and reconstructed Bioplaques ####

## Loop per age
Extended_Bioregions_plates_all_ages <- Voronoi_polygons_for_Bioplaques_all_ages
for (i in seq_along(Extended_Bioregions_plates_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Extended_Bioregions_plates_all_ages)[i]
  
  # Extract plates with Bioregion membership
  sf_use_s2(FALSE)
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[i]] %>% 
    st_make_valid() %>%
    rename(Bioregion_membership = Bioregion)
  # Extract bioregion reconstructed polygons
  Voronoi_polygons_for_Bioplaques_age_i <- Voronoi_polygons_for_Bioplaques_all_ages[[i]] %>% 
    st_make_valid() %>%
    rename(Reconstructed_bioregion = Bioregion)
  sf_use_s2(TRUE)
  
  # Adjust names of bioregions to ensure match
  Voronoi_polygons_for_Bioplaques_age_i$Reconstructed_bioregion[Voronoi_polygons_for_Bioplaques_age_i$Reconstructed_bioregion == "Eastern_Palearctic"] <- "Eastern Palearctic"
  Voronoi_polygons_for_Bioplaques_age_i$Reconstructed_bioregion[Voronoi_polygons_for_Bioplaques_age_i$Reconstructed_bioregion == "Western_Palearctic"] <- "Western Palearctic"
  
  # Intersect Bioregions and Bioplaques
  sf_use_s2(FALSE)
  Extended_Bioregions_plates_i <- Plate_polygons_i %>% 
    st_intersection(y = Voronoi_polygons_for_Bioplaques_age_i)
  sf_use_s2(TRUE)
  
  # Apply Rule: Bioregion identity > Reconstructed Bioplaque identity
  Extended_Bioregions_plates_i$Bioregion <- Extended_Bioregions_plates_i$Bioregion_membership
  Extended_Bioregions_plates_i$Bioregion[Extended_Bioregions_plates_i$Bioregion_membership == "NA"] <- Extended_Bioregions_plates_i$Reconstructed_bioregion[Extended_Bioregions_plates_i$Bioregion_membership == "NA"]
  
  # Aggregate with corrected Bioregions
  sf_use_s2(FALSE)
  Extended_Bioregions_plates_i <- Extended_Bioregions_plates_i %>%
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()
  Extended_Bioregions_plates_i <- st_make_valid(st_buffer(Extended_Bioregions_plates_i, 0))
  Extended_Bioregions_plates_i <- st_crop(x = Extended_Bioregions_plates_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  # View(Extended_Bioregions_plates_i)
  plot(Extended_Bioregions_plates_i[, "Bioregion"], main = main_i)
  
  # Store results
  Extended_Bioregions_plates_all_ages[[i]] <- Extended_Bioregions_plates_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion polygons corrected for age = ",age_i," My - n°", i, "/", length(Extended_Bioregions_plates_all_ages),"\n"))
  }
}

## Save extended bioregion plates as intersection of manually curated and automatically reconstructed bioregions
saveRDS(Extended_Bioregions_plates_all_ages, file = "./outputs/Paleomaps/Extended_Bioregions_plates_all_ages.rds")


## 5.2/ Plot resulting extended bioregion plates ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")
## Load extended bioregion plates as intersection of manually curated and automatically reconstructed bioregions
Extended_Bioregions_plates_all_ages <- readRDS(file = "./outputs/Paleomaps/Extended_Bioregions_plates_all_ages.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA


pdf(file = paste0("./outputs/Paleomaps/Plates_with_curated_Bioregion_borders_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Extended_Bioregions_plates_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Extended_Bioregions_plates_i <- Extended_Bioregions_plates_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Extended_Bioregions_plates_i$Bioregion)]
  
  # Order bioregion factors
  Extended_Bioregions_plates_i$Bioregion <- factor(x = Extended_Bioregions_plates_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  
  # Plot
  Bioregions_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Extended_Bioregions_plates_i,
            mapping = aes(fill = Bioregion), alpha = 0.5,
            show.legend = T) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Bioregion plates in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Bioregions_plot_i)
}
dev.off()



##### 6/ Refine bioregion attribution within plates using reconstructed bioregion convex/concave hulls ####

## Alternative: Reconstruct the subregions, not their plaque (because plaque borders tend to be subducted. Less terrestrial boundaries!)
# However, still need a single polygon per subregion, so use convex_hull and remove overlaps between subregions_hulls

## Load Subregions_sf
Bioregions_sf_Subregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Subregions_level.rds")

### 6.1/ Split subregions crossing the antimeridian ####

# Do it only for Subregion suspected to cross the antimeridian, and NOT the meridian !
antimeridian_subregions <- c("American_Nearctic", "Far_East_Eastern_Palearctic", "Oceanian_Australasia", "Zealandia_Australasia")

# Create East/West hemisphere polygons
West_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = -180, xmax = 0, ymax = 90, ymin = -90), crs = st_crs(4326))))
East_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = 0, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))))
EW_polygons <- rbind(East_hemisphere, West_hemisphere) %>% 
  rename(geometry = x) %>% 
  mutate(Hemisphere = c("HEast", "HWest"))

Subregions_polygons_split <- Bioregions_sf_Subregions_level
for (i in seq_along(antimeridian_subregions))
{
  # i <- 1
  
  # Extract subregion
  Subregion_i <- antimeridian_subregions[i]
  Subregion_multipolygon_i <- Bioregions_sf_Subregions_level %>% 
    filter(Subregion == Subregion_i)
  
  # Cast in POLYGONS
  Subregion_polygons_i <- Subregion_multipolygon_i %>% 
    st_cast(to = "POLYGON")
  
  # Sort by East/West Hemispheres
  sf_use_s2(FALSE)
  Subregion_split_i <- st_intersection(x = Subregion_polygons_i, y = EW_polygons)
  sf_use_s2(TRUE)
  
  # plot(Subregion_split_i[, "Hemisphere"])
  
  # Rename output
  Subregion_split_i <- Subregion_split_i %>% 
    mutate(Subregion = paste0(Subregion, "_",Hemisphere)) %>% 
    select(-Hemisphere)
  
  # Cast back to a MULTIPOLYGON
  sf_use_s2(FALSE)
  Subregion_split_i <- Subregion_split_i %>%
    group_by(Subregion) %>%
    summarize(geometry = st_union(geometry),
              Subregion = first(Subregion))
  sf_use_s2(TRUE)
  
  # Replace previous subregions with its split version in the output
  Subregions_polygons_split <- Subregions_polygons_split %>% 
    filter(Subregion != Subregion_i) %>%
    rbind(Subregion_split_i)
}

plot(Subregions_polygons_split[, "Subregion"])


### 6.2/ Create hull shapes ###

## 6.2.1/ Run concave hulls on split subregions ####

Subregion_hulls_sf <- Subregions_polygons_split %>% 
  # st_convex_hull() %>% 
  st_concave_hull(ratio = 0.1, allow_holes = F)

plot(Subregions_polygons_split[, "Subregion"])
plot(Subregion_hulls_sf[, "Subregion"])


## 6.2.2/ Refine internal borders by removing overlaps between subregions ####

# Rule: Subregion identity > Subregion hulls

# Create a negative polygon when falling outside of Subregions
sf_use_s2(FALSE)
Subregions_polygons_split_negative <- Subregions_polygons_split %>% 
  st_collection_extract(type = "POLYGON") %>%
  nngeo::st_remove_holes() %>% 
  st_union()
Subregions_polygons_split_negative <- st_difference(x = st_as_sfc(st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))),
                                                    # x = st_as_sfc(st_bbox(Subregions_polygons_split)),
                                                    y = Subregions_polygons_split_negative) %>% 
  st_as_sf() %>% 
  mutate(Subregion = "NA") %>% 
  rename(geometry = x) %>%
  select(Subregion, geometry)
sf_use_s2(TRUE)

# plot(Subregions_polygons_split_negative[, "Subregion"])

# Add negative to Subregions_sf
Subregions_polygons_split_with_NA <- rbind(Subregions_polygons_split, Subregions_polygons_split_negative)

# Intersect Subregions and Subregion hulls
sf_use_s2(FALSE)
Subregion_hulls_sf_corrected <- Subregion_hulls_sf %>% 
  rename(Subregion_hull = Subregion) %>% 
  st_intersection(y = Subregions_polygons_split_with_NA)
sf_use_s2(TRUE)

View(Subregion_hulls_sf_corrected)

# Apply Rule: Subregion identity > Subregion hulls
Subregion_hulls_sf_corrected$Subregion_corrected <- Subregion_hulls_sf_corrected$Subregion
Subregion_hulls_sf_corrected$Subregion_corrected[Subregion_hulls_sf_corrected$Subregion == "NA"] <- Subregion_hulls_sf_corrected$Subregion_hull[Subregion_hulls_sf_corrected$Subregion == "NA"]

plot(Subregion_hulls_sf_corrected[, "Subregion_corrected"])

# Aggregate with corrected Subregions
sf_use_s2(FALSE)
Subregion_hulls_sf_corrected <- Subregion_hulls_sf_corrected %>%
  group_by(Subregion_corrected) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup() %>%
  rename(Subregion = Subregion_corrected)
Subregion_hulls_sf_corrected <- st_make_valid(st_buffer(Subregion_hulls_sf_corrected, 0))
Subregion_hulls_sf_corrected <- st_crop(x = Subregion_hulls_sf_corrected,
                                        y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
                                        # y = st_as_sfc(st_bbox(Subregions_polygons_split)))
sf_use_s2(TRUE)

plot(Subregion_hulls_sf_corrected[, "Subregion"])

## 6.2.3/ Keep only the largest polygon in each subregion hull ####

Subregion_hulls_sf_single <- Subregion_hulls_sf_corrected
for (i in 1:nrow(Subregion_hulls_sf_single))
{
  # i <- 1
  
  Subregion_hull_i <- Subregion_hulls_sf_single[i, ]
  
  sf_use_s2(FALSE)
  Subregion_hull_i <- keep_only_largest_polygon(Subregion_hull_i) 
  sf_use_s2(TRUE)
  
  # Store result
  Subregion_hulls_sf_single[i, ] <- Subregion_hull_i
}

plot(Subregion_hulls_sf_single[, "Subregion"])

## 6.2.4/ Fix the Antarctica issue ####

# Create an Antartica patch
Antartica_patch <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = -180, xmax = 180, ymax = -88, ymin = -90), crs = st_crs(4326)))) %>% 
  rename(geometry = x) %>%
  mutate(Subregion = "Antarctica")

Antarctica_current <- Subregion_hulls_sf_single %>% 
  filter(Subregion == "Antarctica")

# Merge with current Antarctica
sf_use_s2(FALSE)
Antarctica_patched <- rbind(Antarctica_current, Antartica_patch) %>% 
  st_union() %>% 
  nngeo::st_remove_holes() %>%
  # st_crop(st_bbox(Subregion_hulls_sf_single)) %>%
  st_as_sf() %>% 
  rename(geometry = x) %>%
  mutate(Subregion = "Antarctica")
sf_use_s2(TRUE)

plot(Antarctica_patched, col = "grey")

# Replace Subregion hull
Subregion_hulls_sf_single <- Subregion_hulls_sf_single %>%
  filter(Subregion != "Antarctica")
Subregion_hulls_sf_single <- rbind(Subregion_hulls_sf_single, Antarctica_patched)  

plot(Subregion_hulls_sf_single[, "Subregion"])

## Save Subregions hulls
saveRDS(object = Subregion_hulls_sf_single, file = "./outputs/Paleomaps/Subregion_hulls_sf_single.rds")


### 6.3/ Plot subregion hulls ####

pdf(file = "./outputs/Paleomaps/Subregion_hulls_map.pdf", width = 12, height = 6)

ggplot(data = Bioregions_sf_Subregions_level) +
  
  geom_sf(fill = NA) +
  
  geom_sf(data = Subregion_hulls_sf_single,
          mapping = aes(fill = Subregion), alpha = 0.5) +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  # Adjust labels
  xlab("") +
  ylab("") +
  
  # Adjust title
  ggtitle(label = paste0("Subregion hulls")) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 8, face = "bold"),
        legend.box.margin = margin(l = 5),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

dev.off()


### 6.4/ Project reconstructed subregion hulls using GPlates ####

## 6.4.1/ Register model features ####

GPlates_path <- "C:/Program Files/GPlates/GPlates 2.5.0/"

# local_model <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/1000_0_rotfile.rot"),
#                                     features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

local_model_optimized <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/optimisation/1000_0_rotfile_MantleOptimised.rot"),
                                              features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

# local_model_no_net <- rgplates::platemodel(rotation = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/optimisation/no_net_rotation_model.rot"),
#                                     features = c("static_polygons" = paste0(GPlates_path, "GeoData/FeatureCollections/AltPlateReconstructions/Muller_etal_2022/shapes_static_polygons_Merdith_et_al.gpml")))

## 6.4.2/ Cast subregion hulls to POINTS to be able to reconstruct polygons ####

# Cast to POINTS
Subregions_list <- Subregion_hulls_sf_single$Subregion
Subregion_hulls_points_sf <- data.frame()
for(i in seq_along(Subregions_list))
{
  # i <- 1
  
  # Extract Subregion
  Subregion_i <- Subregions_list[i]
  
  # Cast the polygon in points
  Subregion_hull_points_i <- Subregion_hulls_sf_single %>% 
    filter(Subregion == Subregion_i) %>% 
    st_cast(to = "POINT") %>% 
    group_by(Subregion) %>%
    mutate(seq_ID = row_number()) %>% 
    ungroup()
  
  # Store output
  Subregion_hulls_points_sf <- rbind(Subregion_hulls_points_sf, Subregion_hull_points_i)
}

table(Subregion_hulls_points_sf$Subregion)

plot(Subregion_hulls_points_sf[Subregion_hulls_points_sf$Subregion == "North_African_Western_Palearctic", "Subregion"])
plot(Subregion_hulls_points_sf[Subregion_hulls_points_sf$Subregion == "North_African_Western_Palearctic", "seq_ID"])

## 6.4.3/ Round and clean redundant coordinates ####

Subregions_list <- Subregion_hulls_sf_single$Subregion
Subregion_hulls_points_sf_rounded <- data.frame()
for (i in seq_along(Subregions_list))
{
  # i <- 1
  
  # Extract Subregion
  Subregion_i <- Subregions_list[i]
  
  # Extract points
  Subregion_hull_points_i <- Subregion_hulls_points_sf %>% 
    filter(Subregion == Subregion_i)
  
  # Round coordinate and remove consecutive duplicates
  Subregion_hull_points_i <- remove_consecutive_duplicated_coordinates_from_points(sf = Subregion_hull_points_i,
                                                                                   round_coordinates = T, digits = 1)
  
  # Store results
  Subregion_hulls_points_sf_rounded <- rbind(Subregion_hulls_points_sf_rounded, Subregion_hull_points_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hull points rounded for Subregion = ",Subregion_i," - n°", i, "/", length(Subregions_list),"\n"))
  }
}

# Register seq_ID of points within polygons
Subregion_hulls_points_sf_rounded <- Subregion_hulls_points_sf_rounded %>% 
  group_by(Subregion) %>%
  mutate(seq_ID = row_number()) %>% 
  ungroup()

table(Subregion_hulls_points_sf$Subregion)
table(Subregion_hulls_points_sf_rounded$Subregion)
table(Subregion_hulls_points_sf_rounded$seq_ID)

plot(Subregion_hulls_points_sf[Subregion_hulls_points_sf$Subregion == "North_African_Western_Palearctic", "Subregion"])
plot(Subregion_hulls_points_sf_rounded[Subregion_hulls_points_sf_rounded$Subregion == "North_African_Western_Palearctic", "seq_ID"])

# Save Subregion hulls border points to provide as input to GPLates
saveRDS(Subregion_hulls_points_sf_rounded, file = "./outputs/Paleomaps/Subregion_hulls_points_sf_rounded.rds")

## 6.4.4/ Obtain reconstruct from POINTS ####

# POINTS allows to preserve order and reconstruct polygons

## Set ages to recover

ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Get the reconstructed points
Subregion_hulls_points_reconstruct_all_ages <- rgplates::reconstruct(x = Subregion_hulls_points_sf_rounded,
                                                                     age = ages_list,
                                                                     # model = local_model,
                                                                     model = local_model_optimized,
                                                                     # model = local_model_no_net,
                                                                     # from = 0,      # Age of the provided object (only for the online recunstruction)
                                                                     listout = T,   # Should the output be a list
                                                                     verbose = T,
                                                                     path.gplates = paste0(GPlates_path,"gplates.exe"),
                                                                     cleanup = T, # Should temp files be removed?
                                                                     dir = "./outputs/Paleomaps/Temp_maps/", # Directory for temp files
                                                                     gmeta = T)  # To keep metadata produced by GPlates in the sf output
names(Subregion_hulls_points_reconstruct_all_ages) <- paste0("Subplaques_",model,"_",ages_list,"My")

View(Subregion_hulls_points_reconstruct_all_ages[[1]])
plot(Subregion_hulls_points_sf_rounded[, "Subregion"])
plot(Subregion_hulls_points_reconstruct_all_ages[[1]][, "Subregion"])
plot(Subregion_hulls_points_reconstruct_all_ages[[2]][, "Subregion"])
plot(Subregion_hulls_points_reconstruct_all_ages[[3]][, "Subregion"])

## Loop per age to reorder points as initial order and remove subregions with less than 4 points remaining

for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages))
{
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i <- Subregion_hulls_points_reconstruct_all_ages[[i]] 
  
  # Reorder points as initial order
  Subregion_hulls_points_reconstruct_age_i <- Subregion_hulls_points_reconstruct_age_i %>% 
    arrange(Subregion, seq_ID)
  
  # Remove subregions with less than 4 points remaining
  table(Subregion_hulls_points_reconstruct_age_i$Subregion)
  Subregions_to_remove <- names(table(Subregion_hulls_points_reconstruct_age_i$Subregion))[table(Subregion_hulls_points_reconstruct_age_i$Subregion) < 4]
  Subregion_hulls_points_reconstruct_age_i <- Subregion_hulls_points_reconstruct_age_i %>% 
    filter(!(Subregion %in% Subregions_to_remove))
  
  # Store updated file
  Subregion_hulls_points_reconstruct_all_ages[[i]] <- Subregion_hulls_points_reconstruct_age_i
}

# Save the raw reconstruct subregion hulls per ages
saveRDS(Subregion_hulls_points_reconstruct_all_ages, file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages.rds")


### 6.5/ Clean reconstructed subregion hulls ####

## 6.5.1/ Split subregion hulls crossing the antimeridian ####

# Do it only for Subregion suspected to cross the antimeridian, and NOT the meridian !
antimeridian_subregions <- c("American_Nearctic_HEast", "American_Nearctic_HWest", "East_Wallacea_Australasia", "Far_East_Eastern_Palearctic_HEast", "Far_East_Eastern_Palearctic_HWest", "Oceanian_Australasia_HEast", "Oceanian_Australasia_HWest", "Zealandia_Australasia_HEast", "Zealandia_Australasia_HWest")

# Loop per ages
Subregion_hulls_points_reconstruct_all_ages_split <- Subregion_hulls_points_reconstruct_all_ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i <- Subregion_hulls_points_reconstruct_all_ages[[i]] 
  
  # Extract only Subregion hulls crossing the antimeridian
  Subregion_hulls_points_reconstruct_age_i_antimeridian <- Subregion_hulls_points_reconstruct_age_i %>% 
    filter(Subregion %in% antimeridian_subregions)
  table(Subregion_hulls_points_reconstruct_age_i_antimeridian$Subregion)
  
  # Split plaques
  Subregion_hulls_points_reconstruct_age_i_antimeridian <- suppressMessages(suppressWarnings(split_by_anti_meridian(sf_points = Subregion_hulls_points_reconstruct_age_i_antimeridian,
                                                                                                                    group_var = "Subregion")))
  # Store them back
  Subregion_hulls_points_reconstruct_age_i_split <- Subregion_hulls_points_reconstruct_age_i %>% 
    filter(!(Subregion %in% antimeridian_subregions)) %>% 
    rbind(Subregion_hulls_points_reconstruct_age_i_antimeridian)
  
  table(Subregion_hulls_points_reconstruct_age_i_split$Subregion)
  plot(Subregion_hulls_points_reconstruct_age_i_split[, "Subregion"], main = main_i)
  
  # Store updated file
  Subregion_hulls_points_reconstruct_all_ages_split[[i]] <- Subregion_hulls_points_reconstruct_age_i_split
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls split for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages),"\n"))
  }
}  

# Save the split reconstruct Subregion hulls per ages
saveRDS(Subregion_hulls_points_reconstruct_all_ages_split, file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages_split.rds")

## 6.5.2/ Filter distance outlier points before reconstruction ####

## Filter by anomaly in distance vs. initial position  (only if not split by anti-meridian)

Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers <- Subregion_hulls_points_reconstruct_all_ages_split

## Make tuned version for each subregion
  # Indian_subcontinent_Indomalaya
  # Arabic_Peninsula_Western_Palearctic
  # Continental_SE_Asia_Indomalaya
  # Central_Eastern_Palearctic

## Filter for Indian_subcontinent_Indomalaya

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers))
{
  # i <- 3
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_dist_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Indian_subcontinent_Indomalaya"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    # Extract points from initial time
    sf_init <- Subregion_hulls_points_sf_rounded[Subregion_hulls_points_sf_rounded$Subregion == Subregion_hull_name, ]
    
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_dist_outliers <- identify_outliers_by_distance_to_initial_position(sf = Subregion_hull_points, 
                                                                                                  sf_init = sf_init,
                                                                                                  quantile_dist_max = 1.0,
                                                                                                  quantile_dist_min = 0.465,
                                                                                                  remove = T, plot = T,
                                                                                                  main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points_no_dist_outliers %>% 
        select(-max_outliers, -min_outliers)
    } else {
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_dist_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers),"\n"))
  }
}  


## Filter for Arabic_Peninsula_Western_Palearctic

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers))
{
  # i <- 3
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_dist_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Arabic_Peninsula_Western_Palearctic"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    # Extract points from initial time
    sf_init <- Subregion_hulls_points_sf_rounded[Subregion_hulls_points_sf_rounded$Subregion == Subregion_hull_name, ]
    
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_dist_outliers <- identify_outliers_by_distance_to_initial_position(sf = Subregion_hull_points, 
                                                                                                  sf_init = sf_init,
                                                                                                  quantile_dist_max = 1.00,
                                                                                                  quantile_dist_min = 0.12,
                                                                                                  remove = T, plot = T,
                                                                                                  main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points_no_dist_outliers %>% 
        select(-max_outliers, -min_outliers)
    } else {
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_dist_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers),"\n"))
  }
}  

## Filter for Continental_SE_Asia_Indomalaya

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers))
{
  # i <- 5
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_dist_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Continental_SE_Asia_Indomalaya"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    # Extract points from initial time
    sf_init <- Subregion_hulls_points_sf_rounded[Subregion_hulls_points_sf_rounded$Subregion == Subregion_hull_name, ]
    
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_dist_outliers <- identify_outliers_by_distance_to_initial_position(sf = Subregion_hull_points, 
                                                                                                  sf_init = sf_init,
                                                                                                  quantile_dist_max = 0.58,
                                                                                                  quantile_dist_min = 0.00,
                                                                                                  remove = T, plot = T,
                                                                                                  main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points_no_dist_outliers %>% 
        select(-max_outliers, -min_outliers)
    } else {
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_dist_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers),"\n"))
  }
}  

## Filter for Central_Eastern_Palearctic

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers))
{
  # i <- 135
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_dist_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Central_Eastern_Palearctic"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    # Extract points from initial time
    sf_init <- Subregion_hulls_points_sf_rounded[Subregion_hulls_points_sf_rounded$Subregion == Subregion_hull_name, ]
    
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_dist_outliers <- identify_outliers_by_distance_to_initial_position(sf = Subregion_hull_points, 
                                                                                                  sf_init = sf_init,
                                                                                                  quantile_dist_max = 0.93,
                                                                                                  quantile_dist_min = 0.00,
                                                                                                  remove = T, plot = T,
                                                                                                  main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points_no_dist_outliers %>% 
        select(-max_outliers, -min_outliers)
    } else {
      Subregion_hull_points_no_dist_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_dist_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers),"\n"))
  }
}  


## Save the reconstruct Subregion hulls per ages cleaned from distance outliers
saveRDS(Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers, file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers.rds")


## 6.5.3/ Filter sequence outlier points before reconstruction ####

## Filter by sequence anomaly as points having abnormal distance to their neighboring points in the sequence

Subregion_hulls_points_reconstruct_all_ages_cleaned <- Subregion_hulls_points_reconstruct_all_ages_no_dist_outliers

## Make tuned version for each subregion
  # Far_East_Eastern_Palearctic_HEast_East
  # Far_East_Eastern_Palearctic_HEast_West
  
## Filter for Far_East_Eastern_Palearctic_HEast_East

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_cleaned))
{
  # i <- 3
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_cleaned)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_cleaned[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_sequence_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Far_East_Eastern_Palearctic_HEast_East"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_sequence_outliers <- identify_outliers_from_sequence(sf = Subregion_hull_points, 
                                                                                    quantile_neighbors = 0.01,
                                                                                    quantile_dist = 0.985,
                                                                                    remove = T, plot = T,
                                                                                    main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_sequence_outliers <- Subregion_hull_points_no_sequence_outliers %>% 
        select(-outliers)
    } else {
      Subregion_hull_points_no_sequence_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_sequence_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_cleaned[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_cleaned),"\n"))
  }
}  

## Filter for Far_East_Eastern_Palearctic_HEast_West

# Loop per ages
for (i in seq_along(Subregion_hulls_points_reconstruct_all_ages_cleaned))
{
  # i <- 5
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_cleaned)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_unfiltered <- Subregion_hulls_points_reconstruct_all_ages_cleaned[[i]] 
  
  # Initiate new sf_object
  Subregion_hulls_points_reconstruct_age_i_no_sequence_outliers <- data.frame()
  
  # Extract Subregion hull name
  Subregion_hull_name <- "Far_East_Eastern_Palearctic_HEast_West"
  
  # Extract points from reconstruct time
  Subregion_hull_points <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
    filter(Subregion == Subregion_hull_name)
  
  if (nrow(Subregion_hull_points) != 0)
  {
    if (nrow(Subregion_hull_points) > 4)
    {
      # Remove outliers
      Subregion_hull_points_no_sequence_outliers <- identify_outliers_from_sequence(sf = Subregion_hull_points, 
                                                                                    quantile_neighbors = 0.05,
                                                                                    quantile_dist = 0.70,
                                                                                    remove = T, plot = T,
                                                                                    main = "Subregion")
      # Remove outlier columns
      Subregion_hull_points_no_sequence_outliers <- Subregion_hull_points_no_sequence_outliers %>% 
        select(-outliers)
    } else {
      Subregion_hull_points_no_sequence_outliers <- Subregion_hull_points
    }
    
    # Store cleaned Subregion hulls back
    Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_age_i_unfiltered %>% 
      filter(!(Subregion %in% Subregion_hull_name)) %>% 
      rbind(Subregion_hull_points_no_sequence_outliers)
    
    # plot(Subregion_hulls_points_reconstruct_age_i_cleaned[, "Subregion"], main = main_i)
    
    # Store updated file
    Subregion_hulls_points_reconstruct_all_ages_cleaned[[i]] <- Subregion_hulls_points_reconstruct_age_i_cleaned
  }
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls cleaned from distance outliers for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_cleaned),"\n"))
  }
}  



# Save the reconstruct Subregion hull points per ages cleaned from sequence outliers
saveRDS(Subregion_hulls_points_reconstruct_all_ages_cleaned, file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages_cleaned.rds")


## 6.5.4/ Reconstruct Subregion hull polygons from border points ####

# ## Load the split reconstruct Subregion hulls per ages
# Subregion_hulls_points_reconstruct_all_ages_split <- readRDS(file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages_split.rds")
# # If skipping the cleaning steps
# Subregion_hulls_points_reconstruct_all_ages_cleaned <- Subregion_hulls_points_reconstruct_all_ages_split

## Load the reconstruct Subregion hull points per ages cleaned from sequence outliers
Subregion_hulls_points_reconstruct_all_ages_cleaned <- readRDS(file = "./outputs/Paleomaps/Subregion_hulls_points_reconstruct_all_ages_cleaned.rds")

## Prepare patches to fix issues with Alaska, Siberia, India, and Myanamar

Alaska_patch <- data.frame(Latitude_dec = c(78.0, 78.0, 74.0, 74.0, 55.0, 55.0, 78.0),
                                 Longitude_dec = c(-120.0, -40.0, -40.0, -20.0, -20.0, -120.0, -120.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Alaska_patch_Nearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Alaska_patch) <- st_crs(4326)

Siberia_East_patch <- data.frame(Latitude_dec = c(89.0, 89.0, 65.0, 65.0, 89.0),
                           Longitude_dec = c(70.0, 120.0, 120.0, 70.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Siberia_East_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Siberia_East_patch) <- st_crs(4326)

Siberia_West_patch <- data.frame(Latitude_dec = c(89.0, 89.0, 85.0, 85.0, 89.0),
                                 Longitude_dec = c(-120.0, -20.0, -20.0, -120.0, -120.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Siberia_West_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Siberia_West_patch) <- st_crs(4326)

Indian_early_patch_1 <- data.frame(Latitude_dec = c(27.5, 20.5, 0, 0, 27.5),
                                 Longitude_dec = c(65.0, 90.0, 90.0, 55.0, 65.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_early_patch_1) <- st_crs(4326)

Indian_early_patch_1_2 <- data.frame(Latitude_dec = c(25.5, 18.5, 0, 0, 25.5),
                                   Longitude_dec = c(65.0, 90.0, 90.0, 55.0, 65.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_early_patch_1_2) <- st_crs(4326)

Indian_early_patch_2 <- data.frame(Latitude_dec = c(33.0, 30.0, 0, 0, 20, 33.0),
                                   Longitude_dec = c(65.0, 70.0, 70.0, 55.0, 58.0, 65.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_2_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_early_patch_2) <- st_crs(4326)

Indian_early_patch_2_2 <- data.frame(Latitude_dec = c(30.0, 28.0, 0, 0, 30.0),
                                   Longitude_dec = c(65.0, 70.0, 70.0, 55.0, 65.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_2_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_early_patch_2_2) <- st_crs(4326)

Indian_early_patch_3 <- data.frame(Latitude_dec = c(25.0, 25.0, 0, 0, 25.0),
                                   Longitude_dec = c(55.0, 70.0, 70.0, 55.0, 55.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_3_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_early_patch_3) <- st_crs(4326)

Indian_late_patch <- data.frame(Latitude_dec = c(-35.0, -35.0, -55.0, -55.0, -35.0),
                                Longitude_dec = c(35.0, 70.0, 70.0, 35.0, 35.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Indian_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Indian_late_patch) <- st_crs(4326)

Myanmar_early_patch <- data.frame(Latitude_dec = c(31.0, 29.0, 21.0, 24.0, 31.0),
                                  Longitude_dec = c(85.0, 95.0, 95.0, 85.0, 85.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Myanmar_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Myanmar_early_patch) <- st_crs(4326)

Myanmar_middle_patch <- data.frame(Latitude_dec = c(33.0, 31.0, 23.0, 23.0, 33.0),
                                   Longitude_dec = c(85.0, 95.0, 95.0, 85.0, 85.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Myanmar_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Myanmar_middle_patch) <- st_crs(4326)

Myanmar_late_patch <- data.frame(Latitude_dec = c(34.0, 34.0, 23.0, 23.0, 34.0),
                                 Longitude_dec = c(85.0, 96.0, 97.0, 85.0, 85.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Myanmar_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Myanmar_late_patch) <- st_crs(4326)

Myanmar_late_patch_2 <- data.frame(Latitude_dec = c(36.0, 36.0, 26.0, 26.0, 36.0),
                                   Longitude_dec = c(90.0, 101.0, 101.0, 90.0, 90.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Myanmar_patch_2_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Myanmar_late_patch_2) <- st_crs(4326)

NW_India_very_early_patch <- data.frame(Latitude_dec = c(37.0, 32.0, 10.0, 10.0, 20.0, 37.0, 37.0),
                                        Longitude_dec = c(72.0, 77.0, 77.0, 70.0, 63.0, 71.0, 72.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "NW_India_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(NW_India_very_early_patch) <- st_crs(4326)

NW_India_very_early_patch_2 <- data.frame(Latitude_dec = c(35.0, 30.0, 10.0, 10.0, 20.0, 30.0, 35.0),
                                        Longitude_dec = c(70.0, 77.0, 77.0, 70.0, 63.0, 65.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "NW_India_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(NW_India_very_early_patch_2) <- st_crs(4326)

NW_India_early_patch <- data.frame(Latitude_dec = c(34.0, 29.0, 10.0, 10.0, 20.0, 34.0, 34.0),
                                   Longitude_dec = c(70.0, 77.0, 77.0, 70.0, 60.0, 68.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "NW_India_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(NW_India_early_patch) <- st_crs(4326)

NW_India_early_patch_2 <- data.frame(Latitude_dec = c(32.0, 27.0, 10.0, 10.0, 20.0, 32.0, 32.0),
                                   Longitude_dec = c(70.0, 77.0, 77.0, 70.0, 60.0, 65.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "NW_India_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(NW_India_early_patch_2) <- st_crs(4326)

NW_India_middle_patch <- data.frame(Latitude_dec = c(28.0, 27.0, 20.0, 20.0, 28.0),
                                    Longitude_dec = c(63.0, 70.0, 70.0, 62.0, 63.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "NW_India_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(NW_India_middle_patch) <- st_crs(4326)

Pakistan_middle_patch <- data.frame(Latitude_dec = c(30.0, 30.0, 25.0, 25.0, 30.0),
                                    Longitude_dec = c(57.0, 63.0, 63.0, 57.0, 57.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Pakistan_middle_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Pakistan_middle_patch) <- st_crs(4326)

Himalaya_very_early_patch <- data.frame(Latitude_dec = c(34.0, 32.0, 28.0, 29.0, 34.0),
                                   Longitude_dec = c(70.0, 88.0, 88.0, 77.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_very_early_patch) <- st_crs(4326)

Himalaya_very_early_patch_2 <- data.frame(Latitude_dec = c(32.0, 30.0, 26.0, 27.0, 32.0),
                                        Longitude_dec = c(70.0, 88.0, 88.0, 77.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_very_early_patch_2) <- st_crs(4326)

Himalaya_early_patch <- data.frame(Latitude_dec = c(28.0, 28.0, 25.0, 25.0, 28.0),
                                   Longitude_dec = c(75.0, 90.0, 90.0, 75.0, 75.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_early_patch) <- st_crs(4326)

Himalaya_early_patch_2 <- data.frame(Latitude_dec = c(28.0, 28.0, 21.0, 21.0, 28.0),
                                   Longitude_dec = c(75.0, 88.0, 88.0, 75.0, 75.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_2_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_early_patch_2) <- st_crs(4326)

Himalaya_middle_patch_0 <- data.frame(Latitude_dec = c(25.0, 23.0, 19.0, 18.0, 25.0),
                                   Longitude_dec = c(75.0, 87.0, 87.0, 75.0, 75.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_middle_patch_0) <- st_crs(4326)

Himalaya_middle_patch_0_2 <- data.frame(Latitude_dec = c(23.0, 23.0, 17.0, 17.0, 23.0),
                                      Longitude_dec = c(87.0, 89.0, 89.0, 87.0, 87.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_2_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_middle_patch_0_2) <- st_crs(4326)


Himalaya_middle_patch_1 <- data.frame(Latitude_dec = c(25.0, 23.0, 18.0, 18.0, 25.0),
                                    Longitude_dec = c(75.0, 87.0, 88.0, 75.0, 75.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_middle_patch_1) <- st_crs(4326)

Himalaya_middle_patch_2 <- data.frame(Latitude_dec = c(25.0, 23.0, 18.0, 18.0, 25.0),
                                    Longitude_dec = c(75.0, 87.0, 89.0, 75.0, 75.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_middle_patch_2) <- st_crs(4326)

Himalaya_late_patch <- data.frame(Latitude_dec = c(25.0, 20.0, 15.0, 15.0, 25.0),
                                      Longitude_dec = c(70.0, 88.0, 88.0, 70.0, 70.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "Himalaya_patch_Eastern_Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(Himalaya_late_patch) <- st_crs(4326)

South_Himalaya_early_patch <- data.frame(Latitude_dec = c(29.0, 29.0, 25.0, 25.0, 29.0),
                                        Longitude_dec = c(79.0, 93.0, 93.0, 79.0, 79.0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Subregion = "South_Himalaya_early_patch_Indomalaya",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders()
st_crs(South_Himalaya_early_patch) <- st_crs(4326)


## Loop per ages
Subregion_hulls_polygons_reconstruct_all_ages <- Subregion_hulls_points_reconstruct_all_ages_cleaned
for (i in seq_along(Subregion_hulls_polygons_reconstruct_all_ages))
{
  # i <- 130
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_points_reconstruct_all_ages_cleaned)[i]
  
  # Extract Subregion hulls points for age i
  Subregion_hulls_points_reconstruct_age_i_cleaned <- Subregion_hulls_points_reconstruct_all_ages_cleaned[[i]] 
  
  # Reconstruct polygons using ordered points
  Subregion_hulls_polygons_reconstruct_age_i <- reconstruct_initial_polygon_with_sfheaders(sf = Subregion_hulls_points_reconstruct_age_i_cleaned)
  # Remove holes
  Subregion_hulls_polygons_reconstruct_age_i <- nngeo::st_remove_holes(Subregion_hulls_polygons_reconstruct_age_i)
  
  # Copy metadata from 1st sf line (except Subregion and seq_ID)
  metadata <- st_drop_geometry(Subregion_hulls_polygons_reconstruct_age_i[1, ]) %>% 
    select(-Subregion, -seq_ID)
  
  ## Add patches to solve issues with Alaska and Siberia
  if (age_i %in% as.character(105:145))
  {
    Alaska_patch_to_bind <- cbind(Alaska_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Alaska_patch_to_bind)
    
    Siberia_East_patch_to_bind <- cbind(Siberia_East_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Siberia_East_patch_to_bind)
    
    Siberia_West_patch_to_bind <- cbind(Siberia_West_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Siberia_West_patch_to_bind)
  }
  
  ## Add patches to solve issues with India
  if (age_i %in% as.character(21:29))
  {
    Indian_early_patch_1_to_bind <- cbind(Indian_early_patch_1, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_early_patch_1_to_bind)
  }
  if (age_i %in% as.character(30:43))
  {
    Indian_early_patch_1_2_to_bind <- cbind(Indian_early_patch_1_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_early_patch_1_2_to_bind)
  }
  if (age_i %in% as.character(21:29))
  {
    Indian_early_patch_2_to_bind <- cbind(Indian_early_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_early_patch_2_to_bind)
  }
  if (age_i %in% as.character(30:32))
  {
    Indian_early_patch_2_2_to_bind <- cbind(Indian_early_patch_2_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_early_patch_2_2_to_bind)
  }
  if (age_i %in% as.character(33:43))
  {
    Indian_early_patch_3_to_bind <- cbind(Indian_early_patch_3, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_early_patch_3_to_bind)
  }
  if (age_i %in% as.character(133:145))
  {
    Indian_late_patch_to_bind <- cbind(Indian_late_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Indian_late_patch_to_bind)
  }

  ## Add patches to solve issues with Myanmar
  if (age_i %in% as.character(52:90))
  {
    Myanmar_early_patch_to_bind <- cbind(Myanmar_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Myanmar_early_patch_to_bind)
  }
  if (age_i %in% as.character(91:108))
  {
    Myanmar_middle_patch_to_bind <- cbind(Myanmar_middle_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Myanmar_middle_patch_to_bind)
  }
  if (age_i %in% as.character(109:114))
  {
    Myanmar_late_patch_to_bind <- cbind(Myanmar_late_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Myanmar_late_patch_to_bind)
  }
  if (age_i %in% as.character(115:145))
  {
    Myanmar_late_patch_2_to_bind <- cbind(Myanmar_late_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Myanmar_late_patch_2_to_bind)
  }

  ## Add patches to solve issues with South_Himalaya
  if (age_i %in% as.character(1:5))
  {
    South_Himalaya_early_patch_to_bind <- cbind(South_Himalaya_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, South_Himalaya_early_patch_to_bind)
  }
  
  ## Add patches to solve issues with Himalaya
  if (age_i %in% as.character(9:12))
  {
    Himalaya_very_early_patch_to_bind <- cbind(Himalaya_very_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_very_early_patch_to_bind)
  }
  if (age_i %in% as.character(13:20))
  {
    Himalaya_very_early_patch_2_to_bind <- cbind(Himalaya_very_early_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_very_early_patch_2_to_bind)
  }
  if (age_i %in% as.character(27:29))
  {
    Himalaya_early_patch_to_bind <- cbind(Himalaya_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_early_patch_to_bind)
  }
  if (age_i %in% as.character(30:45))
  {
    Himalaya_early_patch_2_to_bind <- cbind(Himalaya_early_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_early_patch_2_to_bind)
  }
  if (age_i %in% as.character(46:62))
  {
    Himalaya_middle_patch_0_to_bind <- cbind(Himalaya_middle_patch_0, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_middle_patch_0_to_bind)
  }
  if (age_i %in% as.character(55:60))
  {
    Himalaya_middle_patch_0_2_to_bind <- cbind(Himalaya_middle_patch_0_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_middle_patch_0_2_to_bind)
  }
  if (age_i %in% as.character(63:75))
  {
    Himalaya_middle_patch_1_to_bind <- cbind(Himalaya_middle_patch_1, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_middle_patch_1_to_bind)
  }
  if (age_i %in% as.character(76:80))
  {
    Himalaya_middle_patch_2_to_bind <- cbind(Himalaya_middle_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_middle_patch_2_to_bind)
  }
  if (age_i %in% as.character(81:115))
  {
    Himalaya_late_patch_to_bind <- cbind(Himalaya_late_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Himalaya_late_patch_to_bind)
  }
  
  ## Add patches to solve issues with NW_India
  if (age_i %in% as.character(1:4))
  {
    NW_India_very_early_patch_to_bind <- cbind(NW_India_very_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, NW_India_very_early_patch_to_bind)
  }
  if (age_i %in% as.character(5:8))
  {
    NW_India_very_early_patch_2_to_bind <- cbind(NW_India_very_early_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, NW_India_very_early_patch_2_to_bind)
  }
  if (age_i %in% as.character(9:12))
  {
    NW_India_early_patch_to_bind <- cbind(NW_India_early_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, NW_India_early_patch_to_bind)
  }
  if (age_i %in% as.character(13:20))
  {
    NW_India_early_patch_2_to_bind <- cbind(NW_India_early_patch_2, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, NW_India_early_patch_2_to_bind)
  }
  if (age_i %in% as.character(33:38))
  {
    NW_India_middle_patch_to_bind <- cbind(NW_India_middle_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, NW_India_middle_patch_to_bind)
    
    Pakistan_middle_patch_to_bind <- cbind(Pakistan_middle_patch, metadata)
    # Add patch to the list of subregion hulls
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(Subregion_hulls_polygons_reconstruct_age_i, Pakistan_middle_patch_to_bind)
  }
  
  
  # Plot results
  # plot(Subregion_hulls_polygons_reconstruct_age_i[, "Subregion"])
  
  # Store updated file
  Subregion_hulls_polygons_reconstruct_all_ages[[i]] <- Subregion_hulls_polygons_reconstruct_age_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls polygons reconstructed for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_points_reconstruct_all_ages_cleaned),"\n"))
  }
}

# Save the reconstructed Subregion hull polygons per ages
saveRDS(Subregion_hulls_polygons_reconstruct_all_ages, file = "./outputs/Paleomaps/Subregion_hulls_polygons_reconstruct_all_ages.rds")


### 6.6/ Plot reconstructed subregion hulls ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")


pdf(file = paste0("./outputs/Paleomaps/Plates_with_Subregion_hulls_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Subregion_hulls_polygons_reconstruct_all_ages))
# for (i in 1:65)
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[j]]
  # Plate_boundaries_i <- Paleomaps_sf$Plate_boundaries[[j]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[j]]
  Subregion_hulls_polygons_reconstruct_i <- Subregion_hulls_polygons_reconstruct_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Plate_polygons_i$Bioregion)]
  
  # Order bioregion factors
  Plate_polygons_i$Bioregion <- factor(x = Plate_polygons_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  
  # Plot
  Subregion_hulls_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Subregion_hulls_polygons_reconstruct_i,
            mapping = aes(fill = Subregion), alpha = 0.9,
            show.legend = T) +
    
    # geom_sf(data = Plate_polygons_i,
    #         mapping = aes(fill = Bioregion), alpha = 0.2, col = "grey50") +
    
    # # WGS84 labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec, y = latitude_dec), col = "black") +
    
    # # Mollweide labels
    # geom_label(data = Plate_centroids_i,
    #            mapping = aes(label = shape_ID, x = longitude_dec_Mollweide, y = latitude_dec_Mollweide),
  #            col = "#00000066", alpha = 0.5) +
  
  
  # geom_sf(data = Plate_boundaries_i,
  #         col = "purple") +
  
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
    
    # # Adjust fill legend
    # scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Subregion hulls in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Subregion_hulls_plot_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Subregion hulls plotted for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_polygons_reconstruct_all_ages),"\n"))
  }
}
dev.off()


### 6.7/ Merge reconstructed subregion hulls by Bioregions ####

Bioregions_list <- c("Antarctica", "Afrotropics", "Australasia", "Nearctic", "Neotropics", "Indomalaya", "Eastern_Palearctic", "Western_Palearctic")

Bioregion_hulls_reconstruct_all_ages <- Subregion_hulls_polygons_reconstruct_all_ages
names(Bioregion_hulls_reconstruct_all_ages) <- paste0("Bioregion_hulls_",model,"_",ages_list,"My")
# Loop per ages
for (i in seq_along(Subregion_hulls_polygons_reconstruct_all_ages))
{
  # i <- 41
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Subregion_hulls_polygons_reconstruct_all_ages)[i]
  
  # Extract Subregion hull polygons for age i
  Subregion_hulls_polygons_reconstruct_age_i <- Subregion_hulls_polygons_reconstruct_all_ages[[i]] 
  Subregion_hulls_polygons_reconstruct_age_i$Bioregion <- NA
  
  # Extract Subregions
  all_subregions_i <- unique(Subregion_hulls_polygons_reconstruct_age_i$Subregion)
  
  ## 6.7.1/ Need to deal with the patch overlay before aggregating to Bioregions!!
  
  # Detect patches
  all_patches_i <- all_subregions_i[str_detect(string = all_subregions_i, pattern = "patch")]
  
  if (length(all_patches_i) > 0)
  {
    # Extract patches
    all_patches_polygons_i <- Subregion_hulls_polygons_reconstruct_age_i %>% 
      filter(Subregion %in% all_patches_i) 
    all_patches_union_polygons_i <- all_patches_polygons_i %>%
      st_union()
    # plot(all_patches_polygons_i[, "Subregion"])
    # Remove patch
    no_patches_polygons_i <- Subregion_hulls_polygons_reconstruct_age_i %>% 
      filter(!(Subregion %in% all_patches_i))
    # plot(no_patches_polygons_i[, "Subregion"])
    # Substract patches areas
    sf_use_s2(FALSE)
    sub_patches_polygons_i <- no_patches_polygons_i %>% 
      st_make_valid() %>%
      st_difference(x = ., y = all_patches_union_polygons_i) %>% 
      st_make_valid()
    sf_use_s2(TRUE)
    # plot(sub_patches_polygons_i[, "Subregion"])
    # Put back patches
    Subregion_hulls_polygons_reconstruct_age_i <- rbind(sub_patches_polygons_i, all_patches_polygons_i)
    # plot(Subregion_hulls_polygons_reconstruct_age_i[, "Subregion"])
  }
  
  ## 6.7.2/ Aggregate to Bioregions
  
  # Loop per Bioregions
  for (j in seq_along(Bioregions_list))
  {
    # j <- 7
    
    # Extract Bioregion
    Bioregion_j <- Bioregions_list[j]
    
    # Detect associated subregions
    Subregions_j <- all_subregions_i[str_detect(string = all_subregions_i, pattern = Bioregion_j)]
    Subregion_hulls_polygons_reconstruct_age_i$Bioregion[Subregion_hulls_polygons_reconstruct_age_i$Subregion %in% Subregions_j] <- Bioregion_j
    
  }
  
  # Merge all associated subregions together
  sf_use_s2(FALSE)
  Subregion_hulls_polygons_reconstruct_age_i <- st_make_valid(Subregion_hulls_polygons_reconstruct_age_i)
  Bioregion_hulls_reconstruct_age_i <- Subregion_hulls_polygons_reconstruct_age_i %>% 
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() 
  sf_use_s2(TRUE)
  
  # Plot results
  plot(Bioregion_hulls_reconstruct_age_i[, "Bioregion"], main = main_i)
  
  # Store updated file
  Bioregion_hulls_reconstruct_all_ages[[i]] <- Bioregion_hulls_reconstruct_age_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion hulls aggregated for age = ",age_i," My - n°", i, "/", length(Subregion_hulls_polygons_reconstruct_all_ages),"\n"))
  }
}

# Save the reconstructed Bioregion hulls per ages
saveRDS(Bioregion_hulls_reconstruct_all_ages, file = "./outputs/Paleomaps/Bioregion_hulls_reconstruct_all_ages.rds")


### 6.8/ Plot reconstructed Bioregion hulls ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
colors_list_for_areas_reconstruct <- colors_list_for_areas

bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA

bioregion_names_reconstruct <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern_Palearctic", "Western_Palearctic")
names(colors_list_for_areas_reconstruct) <- bioregion_names_reconstruct
colors_list_for_areas_reconstruct_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas_reconstruct, "grey90")
bioregion_names_reconstruct_with_Antarctica_and_NA <- c("Antarctica", bioregion_names_reconstruct, "NA")
names(colors_list_for_areas_reconstruct_with_Antarctica_and_NA) <- bioregion_names_reconstruct_with_Antarctica_and_NA

pdf(file = paste0("./outputs/Paleomaps/Plates_with_Bioregion_hulls_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Bioregion_hulls_reconstruct_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Bioregion_hulls_reconstruct_i <- Bioregion_hulls_reconstruct_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_reconstruct_present <- bioregion_names_reconstruct_with_Antarctica_and_NA[bioregion_names_reconstruct_with_Antarctica_and_NA %in% unique(Bioregion_hulls_reconstruct_i$Bioregion)]
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_reconstruct_with_Antarctica_and_NA %in% unique(Bioregion_hulls_reconstruct_i$Bioregion)]
  
  # Order bioregion factors
  Bioregion_hulls_reconstruct_i$Bioregion <- factor(x = Bioregion_hulls_reconstruct_i$Bioregion, levels = bioregion_names_reconstruct_present, labels = bioregion_names_present)
  
  # Plot
  Bioregion_hulls_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Bioregion_hulls_reconstruct_i,
            mapping = aes(fill = Bioregion), alpha = 0.9,
            show.legend = T) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Bioregion hulls in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Bioregion_hulls_plot_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion hulls plotted for age = ",age_i," My - n°", i, "/", length(Bioregion_hulls_reconstruct_all_ages),"\n"))
  }
}
dev.off()



### 6.9/ Interpolate membership for empty areas using the closest polygon = Extended bioregion hulls ####

# Voronoi cells around the Bioregion hulls to get extended Bioregion plates covering the whole Earth

## Loop per age
Voronoi_polygons_for_Bioregion_hulls_all_ages <- Bioregion_hulls_reconstruct_all_ages
for (i in seq_along(Voronoi_polygons_for_Bioregion_hulls_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Voronoi_polygons_for_Bioregion_hulls_all_ages)[i]
  
  # Extract bioregion reconstructed hulls
  Bioregion_hulls_reconstruct_age_i <- Voronoi_polygons_for_Bioregion_hulls_all_ages[[i]]
  
  ## 6.9.1/ Cast into POINTS to run tessellation on all points at once #### 
  
  Bioregion_hulls_reconstruct_age_i_points <- Bioregion_hulls_reconstruct_age_i %>% 
    st_collection_extract(type = "POLYGON") %>%
    st_cast(to = "POINT")
  
  # table(Bioregion_hulls_reconstruct_age_i_points$Bioregion)
  
  ## 6.9.2/ Run Derichlet tessellation with st_voronoi to obtain Voronoi cells ####
  
  ## Union in a single MULTIPOINT object
  MULTIPOINT_for_Voronoi_age_i <- sf::st_union(Bioregion_hulls_reconstruct_age_i_points)
  
  ## Run Derichlet tessellation
  # sf_use_s2(FALSE)
  MULTIPOINT_for_Voronoi_age_i <- st_make_valid(MULTIPOINT_for_Voronoi_age_i)
  Voronoi_polygons_for_Bioregion_hulls_age_i <- st_voronoi(MULTIPOINT_for_Voronoi_age_i, dTolerance = 0, bOnlyEdges = FALSE) 
  # sf_use_s2(TRUE)
  
  ## Convert to POLYGON
  Voronoi_polygons_for_Bioregion_hulls_age_i <- Voronoi_polygons_for_Bioregion_hulls_age_i %>%
    st_cast() %>% 
    st_as_sf()
  
  ## Intersect Voronoi cells with initial points to attribute Bioregions
  
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioregion_hulls_age_i <- Voronoi_polygons_for_Bioregion_hulls_age_i %>%
    st_join(y = Bioregion_hulls_reconstruct_age_i_points,
            join = st_intersects, 
            left = TRUE) %>%
    rename(geometry = x)
  Voronoi_polygons_for_Bioregion_hulls_age_i <- st_crop(x = Voronoi_polygons_for_Bioregion_hulls_age_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  Voronoi_polygons_for_Bioregion_hulls_age_i <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioregion_hulls_age_i, 0))
  sf_use_s2(TRUE)
  
  # nrow(Voronoi_polygons_for_Bioregion_hulls_age_i)
  
  # View(Voronoi_polygons_for_Bioregion_hulls_age_i)
  # plot(Voronoi_polygons_for_Bioregion_hulls_age_i)
  
  
  ## Aggregate Voronoi cells per Bioregions
  
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioregion_hulls_age_i <- Voronoi_polygons_for_Bioregion_hulls_age_i %>%
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()
  Voronoi_polygons_for_Bioregion_hulls_age_i <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioregion_hulls_age_i, 0))
  Voronoi_polygons_for_Bioregion_hulls_age_i <- st_crop(x = Voronoi_polygons_for_Bioregion_hulls_age_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  # nrow(Voronoi_polygons_for_Bioregion_hulls_age_i)
  
  # plot(Voronoi_polygons_for_Bioregion_hulls_age_i[, "Bioregion"])
  
  ## 6.9.3/ Correct internal borders using intersection with Bioregions ####
  
  # Rule: Bioregion identity > Extended Bioregion identity
  
  # Create a negative polygon when falling outside of Bioregions
  sf_use_s2(FALSE)
  Bioregion_hulls_reconstruct_negative_age_i <- Bioregion_hulls_reconstruct_age_i %>% 
    st_collection_extract(type = "POLYGON") %>%
    nngeo::st_remove_holes() %>% 
    st_union()
  Bioregion_hulls_reconstruct_negative_age_i <- st_difference(x = st_as_sfc(st_bbox(Bioregion_hulls_reconstruct_age_i)),
                                                              y = Bioregion_hulls_reconstruct_negative_age_i) %>% 
    st_as_sf() %>% 
    mutate(Bioregion = "NA") %>% 
    rename(geometry = x) %>%
    select(Bioregion, geometry)
  sf_use_s2(TRUE)
  
  # plot(Bioregion_hulls_reconstruct_negative_age_i[, "Bioregion"])
  
  # Add negative to Bioregions_sf
  Bioregion_hulls_reconstruct_age_i_with_NA <- rbind(Bioregion_hulls_reconstruct_age_i, Bioregion_hulls_reconstruct_negative_age_i)
  
  # Intersect Bioregions and Extended Bioregions
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected <- Voronoi_polygons_for_Bioregion_hulls_age_i %>% 
    rename(Bioregion_extended = Bioregion) %>% 
    st_intersection(y = Bioregion_hulls_reconstruct_age_i_with_NA)
  sf_use_s2(TRUE)
  
  # Apply Rule: Bioregion identity > Extended Bioregion identity
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion_corrected <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion_corrected[Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion == "NA"] <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion_extended[Voronoi_polygons_for_Bioregion_hulls_age_i_corrected$Bioregion == "NA"]
  
  # Aggregate with corrected Bioregions
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected %>%
    group_by(Bioregion_corrected) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() %>%
    rename(Bioregion = Bioregion_corrected)
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected <- st_make_valid(st_buffer(Voronoi_polygons_for_Bioregion_hulls_age_i_corrected, 0))
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected <- st_crop(x = Voronoi_polygons_for_Bioregion_hulls_age_i_corrected, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  # Remove holes
  sf_use_s2(FALSE)
  Voronoi_polygons_for_Bioregion_hulls_age_i_corrected <- nngeo::st_remove_holes(Voronoi_polygons_for_Bioregion_hulls_age_i_corrected)
  sf_use_s2(TRUE)
  
  ## 6.9.4/ Correct other overlap by applying a priority rule ####

  # Define priority sequence
  priority_sequence <- c("Neotropics", "Nearctic", "Autralasia", "Afrotropics", "Western_Palearctic", "Eastern_Palearctic", "Indomalaya")

  Voronoi_polygons_for_Bioregion_hulls_age_i_masked <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected
  Voronoi_polygons_for_Bioregion_hulls_age_i_mask <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected
  for (j in seq_along(priority_sequence))
  {
    # j <- 1
    Bioregion_j <- rev(priority_sequence)[j]

    # Extract Bioregion
    Bioregion_full_j <- Voronoi_polygons_for_Bioregion_hulls_age_i_mask %>%
      filter(Bioregion == Bioregion_j)

    # Merge all other bioregions to create mask
    sf_use_s2(FALSE)
    mask_j <- Voronoi_polygons_for_Bioregion_hulls_age_i_mask %>%
      filter(Bioregion != Bioregion_j) %>%
      st_union()
    sf_use_s2(TRUE)

    # Apply mask on Bioregion polygon
    sf_use_s2(FALSE)
    Bioregion_masked_j <- Bioregion_full_j %>%
      st_difference(y = mask_j) %>%
      st_as_sf() %>%
      mutate(Bioregion = Bioregion_j)
    sf_use_s2(TRUE)

    # plot(Bioregion_full_j, col = "red")
    # plot(Bioregion_masked_j, col = "limegreen")

    # Remove focal bioregion from the mask for next iteration
    Voronoi_polygons_for_Bioregion_hulls_age_i_mask <- Voronoi_polygons_for_Bioregion_hulls_age_i_mask %>%
      filter(Bioregion != Bioregion_j)

    # Store updated mask version
    Voronoi_polygons_for_Bioregion_hulls_age_i_masked[Voronoi_polygons_for_Bioregion_hulls_age_i_masked$Bioregion == Bioregion_j, ] <- Bioregion_masked_j
  }
  
  # plot(Voronoi_polygons_for_Bioregion_hulls_age_i_corrected[, "Bioregion"], main = main_i)
  # plot(Voronoi_polygons_for_Bioregion_hulls_age_i_masked[, "Bioregion"], main = main_i)
  
  # view(Voronoi_polygons_for_Bioregion_hulls_age_i_corrected)
  
  plot(Voronoi_polygons_for_Bioregion_hulls_age_i_masked[, "Bioregion"], main = main_i)
  
  # Store results
  # Voronoi_polygons_for_Bioregion_hulls_all_ages[[i]] <- Voronoi_polygons_for_Bioregion_hulls_age_i_masked
  Voronoi_polygons_for_Bioregion_hulls_all_ages[[i]] <- Voronoi_polygons_for_Bioregion_hulls_age_i_corrected
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion hulls extended for age = ",age_i," My - n°", i, "/", length(Voronoi_polygons_for_Bioregion_hulls_all_ages),"\n"))
  }
}

## Save Bioregions Voronoi polygons = Extended bioregion hulls for all ages
saveRDS(Voronoi_polygons_for_Bioregion_hulls_all_ages, file = "./outputs/Paleomaps/Voronoi_polygons_for_Bioregion_hulls_all_ages.rds")


### 6.10/ Plot extended bioregion hulls as plates ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
# bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
bioregion_names_reconstruct <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern_Palearctic", "Western_Palearctic")
names(colors_list_for_areas) <- bioregion_names_reconstruct
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
# bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names_reconstruct, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA


pdf(file = paste0("./outputs/Paleomaps/Plates_with_reconstructed_extended_Bioregion_hulls_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Voronoi_polygons_for_Bioregion_hulls_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[j]]
  # Plate_boundaries_i <- Paleomaps_sf$Plate_boundaries[[j]]
  Plate_centroids_i <- Paleomaps_sf$Plate_centroids[[j]]
  Voronoi_polygons_for_Bioregion_hulls_i <- Voronoi_polygons_for_Bioregion_hulls_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Plate_polygons_i$Bioregion)]
  bioregion_names_present_reconstruct <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Voronoi_polygons_for_Bioregion_hulls_i$Bioregion)]
  
  # Order bioregion factors
  Plate_polygons_i$Bioregion <- factor(x = Plate_polygons_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  Voronoi_polygons_for_Bioregion_hulls_i$Bioregion <- factor(x = Voronoi_polygons_for_Bioregion_hulls_i$Bioregion, levels = bioregion_names_present_reconstruct, labels = bioregion_names_present_reconstruct)
  
  # Plot
  Bioregions_hulls_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Voronoi_polygons_for_Bioregion_hulls_i,
            mapping = aes(fill = Bioregion), alpha = 0.5,
            show.legend = T) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    # scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present_reconstruct), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present_reconstruct]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Tectonic plates in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Bioregions_hulls_plot_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion hulls plotted for age = ",age_i," My - n°", i, "/", length(Voronoi_polygons_for_Bioregion_hulls_all_ages),"\n"))
  }
}
dev.off()



##### 7/ Intersect extended Bioregion hulls and plate identities ####

# Rules: 
# Plate identities > Extended Bioregion hulls
# If NA for Plate identities, use Extended Bioregion hulls

## Load Bioregions hulls Voronoi polygons = Extended Bioregion hulls for all ages
Voronoi_polygons_for_Bioregion_hulls_all_ages <- readRDS(file = "./outputs/Paleomaps/Voronoi_polygons_for_Bioregion_hulls_all_ages.rds")

## Load paleomaps in sf format with Bioregion membership
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")


## 7.1/ Intersect plates with Bioregion membership and Extended Bioregion hulls ####

## Loop per age
Corrected_Bioregions_plates_all_ages <- Voronoi_polygons_for_Bioregion_hulls_all_ages
names(Corrected_Bioregions_plates_all_ages) <- paste0("Bioregions_plates_",model,"_",ages_list,"My")
for (i in seq_along(Corrected_Bioregions_plates_all_ages))
{
  # i <- 5
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Corrected_Bioregions_plates_all_ages)[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract plates with Bioregion membership
  sf_use_s2(FALSE)
  Plate_polygons_i  <- Paleomaps_sf$Plate_polygons[[j]] %>% 
    st_make_valid() %>%
    rename(Bioregion_membership = Bioregion)
  # Extract extended Bioregion hulls
  Voronoi_polygons_for_Bioregion_hulls_age_i <- Voronoi_polygons_for_Bioregion_hulls_all_ages[[i]] %>% 
    st_make_valid() %>%
    rename(Reconstructed_bioregion = Bioregion)
  sf_use_s2(TRUE)
  
  # Adjust names of bioregions to ensure match
  Voronoi_polygons_for_Bioregion_hulls_age_i$Reconstructed_bioregion[Voronoi_polygons_for_Bioregion_hulls_age_i$Reconstructed_bioregion == "Eastern_Palearctic"] <- "Eastern Palearctic"
  Voronoi_polygons_for_Bioregion_hulls_age_i$Reconstructed_bioregion[Voronoi_polygons_for_Bioregion_hulls_age_i$Reconstructed_bioregion == "Western_Palearctic"] <- "Western Palearctic"
  
  # plot(Plate_polygons_i[, "Bioregion_membership"])
  # plot(Voronoi_polygons_for_Bioregion_hulls_age_i[, "Reconstructed_bioregion"])
  
  # Intersect Bioregions and extended Bioregion hulls
  sf_use_s2(FALSE)
  Corrected_Bioregions_plates_i <- Plate_polygons_i %>% 
    st_intersection(y = Voronoi_polygons_for_Bioregion_hulls_age_i)
  sf_use_s2(TRUE)
  
  # Apply Rule: Bioregion identity > Reconstructed Bioplaque identity
  Corrected_Bioregions_plates_i$Bioregion <- Corrected_Bioregions_plates_i$Bioregion_membership
  Corrected_Bioregions_plates_i$Bioregion[Corrected_Bioregions_plates_i$Bioregion_membership == "NA"] <- Corrected_Bioregions_plates_i$Reconstructed_bioregion[Corrected_Bioregions_plates_i$Bioregion_membership == "NA"]
  
  # Aggregate with corrected Bioregions
  sf_use_s2(FALSE)
  Corrected_Bioregions_plates_i <- Corrected_Bioregions_plates_i %>%
    group_by(Bioregion) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()
  Corrected_Bioregions_plates_i <- st_make_valid(st_buffer(Corrected_Bioregions_plates_i, 0))
  Corrected_Bioregions_plates_i <- st_crop(x = Corrected_Bioregions_plates_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)))
  sf_use_s2(TRUE)
  
  # View(Corrected_Bioregions_plates_i)
  plot(Corrected_Bioregions_plates_i[, "Bioregion"], main = main_i)
  
  # Store results
  Corrected_Bioregions_plates_all_ages[[i]] <- Corrected_Bioregions_plates_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion polygons corrected for age = ",age_i," My - n°", i, "/", length(Corrected_Bioregions_plates_all_ages),"\n"))
  }
}

## Save corrected bioregion plates as intersection of manually curated and automatically reconstructed bioregions from hulls
saveRDS(Corrected_Bioregions_plates_all_ages, file = "./outputs/Paleomaps/Corrected_Bioregions_plates_all_ages.rds")


## 7.2/ Plot corrected Bioregion plates ####

# Plot as overlay on plate with bioregion membership

# Set ages to plot
ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA


pdf(file = paste0("./outputs/Paleomaps/Corrected_Bioregion_plates_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Corrected_Bioregions_plates_all_ages))
{
  # i <- 1
  
  # Extract age
  age_i <- ages_list[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract data
  Coastlines_i <- Paleomaps_sf$Coastlines[[j]]
  Corrected_Bioregions_plates_i <- Corrected_Bioregions_plates_all_ages[[i]]
  
  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Corrected_Bioregions_plates_i$Bioregion)]

  # Order bioregion factors
  Corrected_Bioregions_plates_i$Bioregion <- factor(x = Corrected_Bioregions_plates_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)

  # Plot
  Bioregions_plates_plot_i <- ggplot(data = Coastlines_i) + 
    
    geom_sf(data = Coastlines_i,
            fill = "grey90", col = "black") +
    
    geom_sf(data = Corrected_Bioregions_plates_i,
            mapping = aes(fill = Bioregion), alpha = 0.5,
            show.legend = T) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +

    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Bioregion plates in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Bioregions_plates_plot_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion corrected plates plotted for age = ",age_i," My - n°", i, "/", length(Corrected_Bioregions_plates_all_ages),"\n"))
  }
}
dev.off()



##### 8/ Attribute bioregions to coastlines ####

### 8.1/ Intersect Coastlines and bioregions extended plates ####

## Load paleomaps in sf format
Paleomaps_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_sf.rds")

## Load corrected bioregion plates as intersection of manually curated and automatically reconstructed bioregions from hulls
Corrected_Bioregions_plates_all_ages <- readRDS(file = "./outputs/Paleomaps/Corrected_Bioregions_plates_all_ages.rds")

# Initiate new object
Paleomaps_with_bioregions_sf <- list()

## Loop per age
for (i in seq_along(Corrected_Bioregions_plates_all_ages))
{
  # i <- 5
  
  # Extract age
  age_i <- ages_list[i]
  main_i <- names(Corrected_Bioregions_plates_all_ages)[i]
  
  # Find matching index
  # j <- age_i + 1
  ages_names <- str_remove(string = names(Paleomaps_sf$Coastlines), pattern = "Coastlines_")
  ages_names <- str_remove(string = ages_names, pattern = paste0("My_", model))
  j <- which(ages_names == age_i)
  
  # Extract coastlines
  sf_use_s2(FALSE)
  Coastlines_i  <- Paleomaps_sf$Coastlines[[j]] %>% 
    st_make_valid()
  
  # Extract corrected Bioregion plates
  Corrected_Bioregions_plates_age_i <- Corrected_Bioregions_plates_all_ages[[i]] %>% 
    st_make_valid()
  sf_use_s2(TRUE)
  
  # plot(Coastlines_i[, "shape_ID"])
  # plot(Corrected_Bioregions_plates_age_i[, "Bioregion"])
  
  # Intersect Coastlines and extended Bioregion hulls
  sf_use_s2(FALSE)
  Coastlines_with_Bioregions_i <- Coastlines_i %>% 
    st_intersection(y = Corrected_Bioregions_plates_age_i) %>%
    rename(geometry = x)
  sf_use_s2(TRUE)
  
  # View(Coastlines_with_Bioregions_i)
  plot(Coastlines_with_Bioregions_i[, "Bioregion"], main = main_i)
  
  # Store results
  Paleomaps_with_bioregions_sf[[i]] <- Coastlines_with_Bioregions_i
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Bioregion coastlines detected for age = ",age_i," My - n°", i, "/", length(Corrected_Bioregions_plates_all_ages),"\n"))
  }
}
names(Paleomaps_with_bioregions_sf) <- paste0("Bioregions_",model,"_",ages_list,"My")

## Save Bioregion coastlines
saveRDS(Paleomaps_with_bioregions_sf, file = "./outputs/Paleomaps/Paleomaps_with_bioregions_sf.rds")


### 8.2/ Plot Bioregion coastlines in time ####

# Set ages to plot
ages_list <- c(0:start_time)
# ages_list <- c(0, 20, 50, 100, 145)

## Load Bioregion coastlines in sf format
Paleomaps_with_bioregions_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_with_bioregions_sf.rds")

# Load color scheme for bioregions
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica_and_NA <- c("grey50", colors_list_for_areas, "grey90")
bioregion_names_with_Antarctica_and_NA <- c("Antarctica", bioregion_names, "NA")
names(colors_list_for_areas_with_Antarctica_and_NA) <- bioregion_names_with_Antarctica_and_NA

## Create polygons of extreme values to ensure all the bbox is plotted
E_patch <- data.frame(Latitude_dec = c(0, 0, 0.01, 0.01, 0),
                                   Longitude_dec = c(179.9, 180.0, 180.0, 179.9, 179.9)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Australasia",
         Subregion = "Australasia",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(E_patch) <- st_crs(4326)

W_patch <- data.frame(Latitude_dec = c(0, 0, 0.01, 0.01, 0),
                       Longitude_dec = c(-179.9, -180.0, -180.0, -179.9, -179.9)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Australasia",
         Subregion = "Australasia",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(W_patch) <- st_crs(4326)


## Loop per ages
pdf(file = paste0("./outputs/Paleomaps/Bioregion_coastlines_MULLER2022_all_ages.pdf"),
    width = 16, height = 8)

for (i in seq_along(Paleomaps_with_bioregions_sf))
{
  # i <- 1
  
  # Extract data
  Coastlines_with_bioregions_i <- Paleomaps_with_bioregions_sf[[i]] 
  
  # Extract age
  age_i <- str_remove(string = names(Paleomaps_with_bioregions_sf)[i], pattern = paste0("Bioregions_", model, "_"))
  age_i <- str_remove(string = age_i, pattern = paste0("My"))
  
  
  # Copy metadata from 1st sf line (except Bioregion)
  metadata <- st_drop_geometry(Coastlines_with_bioregions_i[1, ]) %>% 
    select(-Bioregion)
  
  # Add extreme points to ensure all the bbox is plotted
  E_patch_to_bind <- cbind(E_patch, metadata)
  W_patch_to_bind <- cbind(W_patch, metadata)
  # Add patches to the list of bioregions
  Coastlines_with_bioregions_i <- rbind(Coastlines_with_bioregions_i, E_patch_to_bind, W_patch_to_bind)

  # Adjust list of bioregions
  bioregion_names_present <- bioregion_names_with_Antarctica_and_NA[bioregion_names_with_Antarctica_and_NA %in% unique(Coastlines_with_bioregions_i$Bioregion)]
  
  # Order bioregion factors
  Coastlines_with_bioregions_i$Bioregion <- factor(x = Coastlines_with_bioregions_i$Bioregion, levels = bioregion_names_present, labels = bioregion_names_present)
  
  # Plot
  Coastlines_with_bioregions_plot_i <- ggplot(data = Coastlines_with_bioregions_i) + 
    
    geom_sf(data = Coastlines_with_bioregions_i,
            mapping = aes(fill = Bioregion), alpha = 0.5, col = "black",
            show.legend = T) +
    
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Adjust fill legend
    scale_fill_manual("Bioregions", labels = c(bioregion_names_present), values = colors_list_for_areas_with_Antarctica_and_NA[bioregion_names_present]) +
    
    # Adjust labels
    xlab("") +
    ylab("") +
    
    # Adjust title
    ggtitle(label = paste0("Bioregions in ", model, " model\n",
                           age_i, " My")) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
          legend.text = element_text(size = 6, face = "bold"),
          legend.box.margin = margin(l = 5),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(Coastlines_with_bioregions_plot_i)
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Plot bioregion coastlines for age = ",age_i," My - n°", i, "/", length(Paleomaps_with_bioregions_sf),"\n"))
  }
}
dev.off()


### 8.3/ Make a GIF in forward time ####

source("./functions/image_resize_and_write_gif.R")

pdf_pointer_bioregions_coastlines <- magick::image_read_pdf(path = paste0("./outputs/Paleomaps/Bioregion_coastlines_MULLER2022_all_ages.pdf"),
                                                  pages = NULL, density = 150)

# Get info on PDF pages
magick::image_info(pdf_pointer_bioregions_coastlines)

# Reverse order for forward time
pdf_pointer_bioregions_coastlines <- rev(pdf_pointer_bioregions_coastlines)

image_resize_and_write_gif(image = pdf_pointer_bioregions_coastlines,
                           path =  paste0("./outputs/Paleomaps/Bioregion_coastlines_MULLER2022_all_ages_forward.gif"),
                           delay = 1/5, # Time between frames in seconds
                           width = 1200, height = 600,
                           loop = FALSE,
                           progress = TRUE)



