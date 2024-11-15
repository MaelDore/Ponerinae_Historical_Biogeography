##### Script 11: Map Species richness #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Map current species richness
  # Overall
  # Per Genera
# Map occurrences
# Plot SR gradients along long/lat bands

###

### Inputs

# Occurrence data
# Bioregion maps

###

### Outputs

# Maps of current species richness
  # Overall
  # Per Genera
# Maps of occurrences
# Binary table of taxa P/A in long/lat bands
# Plot of Species Richness along long/lat bands

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(phytools)
library(ape)
library(raster)
library(sf)
library(spatstat)
library(alphahull)
library(geosphere)
library(nngeo)
library(qpdf)

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load bioregions map ###

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

# Convert to Mollweide
Bioregions_sf_Bioregions_level_Mollweide <- st_transform(x = Bioregions_sf_Bioregions_level, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
 
# Save Bioregions map projected in Mollweide
saveRDS(object = Bioregions_sf_Bioregions_level_Mollweide, file = "./outputs/Species_richness_maps/Bioregions_sf_Bioregions_level_Mollweide.rds")

### 1.3/ Load overall bioregion richness data ###

# Bioregion_richness_data <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
Bioregion_richness_data <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

## Seems to miss some taxa (maybe the recent one having appeared within the last time step)

### 1.4/ Load Ponerinae occurrence data ####

Biogeographic_database_Ponerinae_curated_no_duplicates <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_curated_no_duplicates.rds")


##### 2/ Compute overall species richness map ####

### 2.1/ Define raster grid ####

# Create grid in WGS84
worldwide_grid_WGS84 <- raster::raster(xmn = -180, ymn = -60, xmx = 180, ymx = 85,
                                       resolution = 0.25, # Resolution to quarter degree cells = 30km
                                       crs = "+proj=longlat +datum=WGS84")
# Fill terrestrial lands with zeros
Bioregions_sp_Bioregions_level <- as(Bioregions_sf_Bioregions_level, 'Spatial')
terrestrial_bg_WGS84 <- raster::rasterize(x = Bioregions_sp_Bioregions_level, y = worldwide_grid_WGS84, field = 0, background = NA)

plot(terrestrial_bg_WGS84)

# Convert to Mollweide
terrestrial_bg_Mollweide <- raster::projectRaster(from = terrestrial_bg_WGS84,
                                                  crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
                                                  method = "ngb")
plot(terrestrial_bg_Mollweide)

## Save raster maps
saveRDS(object = terrestrial_bg_WGS84, file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")
saveRDS(object = terrestrial_bg_Mollweide, file = "./outputs/Species_richness_maps/terrestrial_bg_Mollweide.rds")

### 2.2/ Rasterize (and thin) occurrences ####

# Extract list of taxa to use
taxa_list_full <- unique(Biogeographic_database_Ponerinae_curated_no_duplicates$Current_name)
taxa_list_full <- taxa_list_full[order(taxa_list_full)]
taxa_list_in_analyses <- unique(Biogeographic_database_Ponerinae_curated_no_duplicates$Current_name[Biogeographic_database_Ponerinae_curated_no_duplicates$Current_status == "valid" | Biogeographic_database_Ponerinae_curated_no_duplicates$In_phylogeny])
taxa_list_in_analyses <- taxa_list_in_analyses[order(taxa_list_in_analyses)] 
  
## Loop per taxa in analyses
Ponerinae_occurrences_species_stack_WGS84 <- stack()
for (i in seq_along(taxa_list_in_analyses))
{
  # i <- 1
  
  # Extract taxa name
  taxa_i <- taxa_list_in_analyses[i]
  
  # Extract occurrences
  occ_df_i <- Biogeographic_database_Ponerinae_curated_no_duplicates %>% 
    filter(Current_name == taxa_i)
  
  # Convert to sf/sp
  occ_sf_i <- st_as_sf(x = occ_df_i, coord = c("Longitude_dec", "Latitude_dec"), crs = "+proj=longlat +datum=WGS84")
  occ_sp_i <- as(occ_sf_i, 'Spatial')
  
  # Rasterize
  raster_i <- raster::rasterize(x = occ_sp_i, y = terrestrial_bg_WGS84, field = 1, update = TRUE)
  
  # Remove non-terrestrial
  # raster_i <- raster::mask(x = raster_i, mask = terrestrial_bg_WGS84)
  
  # Check non-null
  nb_cells <- sum(raster_i[], na.rm = T)
  if (nb_cells < 1)
  {
    stop(paste0(taxa_i, " has no records"))
  }
  
  # Store in stack
  Ponerinae_occurrences_species_stack_WGS84 <- addLayer(x = Ponerinae_occurrences_species_stack_WGS84, raster_i) 
  names(Ponerinae_occurrences_species_stack_WGS84) <- taxa_list_in_analyses[1:i]
  
  # Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Occurrence raster created for taxa n°", i, "/", length(taxa_list_in_analyses),"\n"))
  }
}

plot(Ponerinae_occurrences_species_stack_WGS84[[1:9]])

## Convert to Mollweide
Ponerinae_occurrences_species_stack_Mollweide <- raster::projectRaster(from = Ponerinae_occurrences_species_stack_WGS84,
                                                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
                                                                       method = "ngb")
plot(Ponerinae_occurrences_species_stack_Mollweide[[1:9]])

# Load data in RAM
Ponerinae_occurrences_species_stack_WGS84 <- readAll(Ponerinae_occurrences_species_stack_WGS84)
Ponerinae_occurrences_species_stack_Mollweide <- readAll(Ponerinae_occurrences_species_stack_Mollweide)

## Save occurrence rasters
saveRDS(object = Ponerinae_occurrences_species_stack_WGS84, file = "./outputs/Species_richness_maps/Ponerinae_occurrences_species_stack_WGS84.rds")
saveRDS(object = Ponerinae_occurrences_species_stack_Mollweide, file = "./outputs/Species_richness_maps/Ponerinae_occurrences_species_stack_Mollweide.rds")

Ponerinae_occurrences_species_stack_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occurrences_species_stack_WGS84.rds")
Ponerinae_occurrences_species_stack_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occurrences_species_stack_Mollweide.rds")

### 2.3/ Compute raw species richness (no interpolation) ####

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")

Ponerinae_species_richness_raw_WGS84 <- raster::calc(x = Ponerinae_occurrences_species_stack_WGS84, fun = sum, na.rm = T)
Ponerinae_species_richness_raw_Mollweide <- raster::calc(x = Ponerinae_occurrences_species_stack_Mollweide, fun = sum, na.rm = T)

plot(Ponerinae_species_richness_raw_WGS84, col = pal_bl_red_Mannion)
plot(Ponerinae_species_richness_raw_Mollweide, col = pal_bl_red_Mannion)

# Load data in RAM
Ponerinae_species_richness_raw_WGS84 <- readAll(Ponerinae_species_richness_raw_WGS84)
Ponerinae_species_richness_raw_Mollweide <- readAll(Ponerinae_species_richness_raw_Mollweide)

## Save raw species richness rasters
saveRDS(object = Ponerinae_species_richness_raw_WGS84, file = "./outputs/Species_richness_maps/Ponerinae_species_richness_raw_WGS84.rds")
saveRDS(object = Ponerinae_species_richness_raw_Mollweide, file = "./outputs/Species_richness_maps/Ponerinae_species_richness_raw_Mollweide.rds")


# ### 2.4/ Run Inverse Weighted Distance interpolation on overall distribution ####
# 
# ## 2.4.1/ Create point pattern object
# 
# library(spatstat)
# ?idw
# 
# Ponerinae_species_richness_sp_WGS84 <- raster::rasterToPoints(x = Ponerinae_species_richness_raw_WGS84, fun = function(x) {x > 0}, spatial = TRUE)
# plot(Ponerinae_species_richness_sp_WGS84)
# 
# class(Ponerinae_species_richness_sp_WGS84)
# 
# obs_window_WGS84 <- owin(xrange = bbox(Ponerinae_species_richness_raw_WGS84)[1,], 
#                          yrange = bbox(Ponerinae_species_richness_raw_WGS84)[2,])
# 
# Ponerinae_species_richness_ppp_WGS84 <- spatstat::ppp(x = Ponerinae_species_richness_sp_WGS84@coords[,1],
#                                                       y = Ponerinae_species_richness_sp_WGS84@coords[,2],
#                                                       marks = Ponerinae_species_richness_sp_WGS84@data$layer,
#                                                       window = obs_window_WGS84)
# 
# ## 2.4.2/ Run IWD over different power laws
# 
# power_list <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
# 
# Ponerinae_species_richness_IDW_stack_WGS84 <- stack()
# names(Ponerinae_species_richness_IDW_stack_WGS84) <- NULL
# for (i in seq_along(power_list))
# {
#   # i <- 1
#   
#   # Extract power
#   power_i <- power_list[i]
#   
#   # Run IDW
#   Ponerinae_species_richness_IDW_WGS84 <- spatstat.explore::idw(X = Ponerinae_species_richness_ppp_WGS84,
#                                                                 power = power_i, at = "pixels")
#   # Convert to raster
#   Ponerinae_species_richness_IDW_WGS84 <- raster::raster(Ponerinae_species_richness_IDW_WGS84)
#   Ponerinae_species_richness_IDW_WGS84 <- raster::resample(x = Ponerinae_species_richness_IDW_WGS84, y = Ponerinae_species_richness_raw_WGS84)
#   
#   # Remove non-terrestrial pixels
#   Ponerinae_species_richness_IDW_WGS84 <- raster::mask(x = Ponerinae_species_richness_IDW_WGS84, mask = terrestrial_bg_WGS84)
#   
#   # Adjust names
#   if (i == 1)
#   {
#     old_names <- NULL
#   } else {
#     old_names <- names(Ponerinae_species_richness_IDW_stack_WGS84)
#   }
#   
#   # Stack results
#   Ponerinae_species_richness_IDW_stack_WGS84 <- addLayer(x = Ponerinae_species_richness_IDW_stack_WGS84, Ponerinae_species_richness_IDW_WGS84)
#   new_names <- c(old_names, paste0("Power_", power_i))
#   names(Ponerinae_species_richness_IDW_stack_WGS84) <- new_names
#   
#   # Print progress
#   if (i %% 1 == 0)
#   {
#     cat(paste0(Sys.time(), " - IDW created for Power = ", power_i, " - n°", i, "/", length(power_list),"\n"))
#   }
# }
# 
# plot(Ponerinae_species_richness_IDW_stack_WGS84,
#      col = heat.colors(20))
# 
# ## Result is not convincing...


### 2.5/ Get species range using alpha-hull buffers ####

library(alphahull)

# Load functions to convert alpha-hull in spatial polygons
source("./functions/alpha_functions.R")

## Loop per taxa in analyses
Ponerinae_species_alpha_hull_stack_WGS84 <- stack()
Ponerinae_species_alpha_hull_stack_Mollweide <- stack()
Q80_all_taxa <- c()
for (i in seq_along(taxa_list_in_analyses))
# for (i in taxa_to_update)
{
  # i <- 5
  
  # Extract taxa name
  taxa_i <- taxa_list_in_analyses[i]
  
  # Extract occurrences
  occ_df_i <- Biogeographic_database_Ponerinae_curated_no_duplicates %>% 
    filter(Current_name == taxa_i)
  
  # Convert to sf/sp
  occ_sf_i <- st_as_sf(x = occ_df_i, coord = c("Longitude_dec", "Latitude_dec"), crs = "+proj=longlat +datum=WGS84")
  occ_sp_i <- as(occ_sf_i, 'Spatial')
  
  ## 2.5.1/ Compute the buffer distance as the 80% quantile of the distance to the closest occurrence points among occurrence points for each taxa ####
  
  occ_dist_i = geosphere::distm(x = occ_sp_i@coords)/1000 # Geometric distance on WGS84 ellipsoid, in km
  
  if (nrow(occ_dist_i) >= 5) # only if at least 5 points
  {
    # Remove diagonal (distance to themselves = 0 km)
    diag(occ_dist_i) <- NA
    # Compute distance to closest other occurrence point
    occ_dist_min_i <- apply(X = occ_dist_i, MARGIN = 2, FUN = min, na.rm = T)
    
    # Get the 80% quantile
    Q80_i <- round(quantile(x = occ_dist_min_i, probs = 0.80),1)
    # Use a maximum buffer of 500km to avoid poorly informed ranges to encompass large areas
    Q80_i <- min(500, Q80_i)
  
  } else { # If less than 5 points, apply 1° buffer
    Q80_i <- 0
  }

  Q80_all_taxa <- c(Q80_all_taxa, Q80_i)
  
  ## 2.5.2/ Generate alpha-hull buffer ####

  if (nrow(occ_dist_i) >= 3) # Only if at least 3 points
  {
    # Projection of the spatial object of occurrences in Mollweide
    proj_occ_sp_i <- spTransform(x = occ_sp_i,
                                 CRSobj = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

    tryCatch( # Try an alpha-hull first
      {
        range_sp_i <- ahull_to_SPLDF(ahull(proj_occ_sp_i@coords, alpha = 1000)) # alpha diameter = 1000km
      },
      error = function(e) # If fails, do an alpha shape instead
      {
        cat(paste0("Alpha-hull failed for ", taxa_i, ", use alpha-shape instead \n"))
        cat("ERROR :",conditionMessage(e), "\n")   # Display the error message but do not stop the function
        range_sp_i <<- ashape_to_SPLDF(ashape(proj_occ_sp_i@coords, alpha = 1000)) # alpha diameter = 1000km
      })

    # plot(range_sp_i)

    # Save crs of the newly created alpha-hull
    range_sp_i@proj4string <- proj_occ_sp_i@proj4string

    # Convert to sf
    range_sf_i <- st_as_sf(range_sp_i)

  } else { # If less than 3 points, only apply 1° buffer around occurrences
    range_sf_i <- st_transform(x = occ_sf_i, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  }

  # Extend the shape with a buffer based on the quality of presence sampling (Q80_i)
  # As a minimum, use a 1° buffer = 111.32 km
  range_sf_i <- st_buffer(x = range_sf_i, dist = max(111.32, Q80_i))

  # Remove holes
  range_sf_i <- nngeo::st_remove_holes(x = range_sf_i, max_area = 0)
  # range_sp_i <- spatialEco::remove.holes(x = range_sp_i)

  # Retransposition in WGS84
  range_sf_WGS84_i <- st_transform(x = range_sf_i, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  # range_sp_WGS84_i <- spTransform(range_sp_i, CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  # Rasterization
  range_raster_WGS84_i <- raster::rasterize(x = range_sf_WGS84_i, y = terrestrial_bg_WGS84, field = 1, update = TRUE)
  range_raster_Mollweide_i <- raster::rasterize(x = range_sf_i, y = terrestrial_bg_Mollweide, field = 1, update = TRUE)

  # Remove non-terrestrial areas
  range_raster_WGS84_i <- raster::mask(x = range_raster_WGS84_i, mask = terrestrial_bg_WGS84)
  range_raster_Mollweide_i <- raster::mask(x = range_raster_Mollweide_i, mask = terrestrial_bg_Mollweide)

  # plot(range_raster_WGS84_i)
  # plot(range_raster_Mollweide_i)

  # Store range maps
  Ponerinae_species_alpha_hull_stack_WGS84 <- addLayer(Ponerinae_species_alpha_hull_stack_WGS84, range_raster_WGS84_i)
  Ponerinae_species_alpha_hull_stack_Mollweide <- addLayer(Ponerinae_species_alpha_hull_stack_Mollweide, range_raster_Mollweide_i)
  names(Ponerinae_species_alpha_hull_stack_WGS84) <- taxa_list_in_analyses[1:i]
  names(Ponerinae_species_alpha_hull_stack_Mollweide) <- taxa_list_in_analyses[1:i]
  
  # # Replace updated range maps
  # Ponerinae_species_alpha_hull_stack_WGS84[[i]] <- range_raster_WGS84_i
  # Ponerinae_species_alpha_hull_stack_Mollweide[[i]] <- range_raster_Mollweide_i
  
  # Print progress
  if (i %% 100 == 0)
  {
    # Save temporary final stacks
    saveRDS(Ponerinae_species_alpha_hull_stack_WGS84, file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))
    saveRDS(Ponerinae_species_alpha_hull_stack_Mollweide, file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_Mollweide.rds"))
    
    cat(paste0(Sys.time(), " - Buffered range map created for taxa n°", i, "/", length(taxa_list_in_analyses),"\n"))
  }
}

# names(Ponerinae_species_alpha_hull_stack_WGS84) <- taxa_list_in_analyses
# names(Ponerinae_species_alpha_hull_stack_Mollweide) <- taxa_list_in_analyses

# Explore taxa-specific buffer size used
names(Q80_all_taxa) <- taxa_list_in_analyses
# names(Q80_all_taxa) <- taxa_to_update
hist(Q80_all_taxa)
summary(Q80_all_taxa)
Large_buffer <- Q80_all_taxa[Q80_all_taxa > 500]
Large_buffer[order(Large_buffer, decreasing = T)]
taxa_to_update <- which(taxa_list_in_analyses %in% names(Large_buffer))

plot(Ponerinae_species_alpha_hull_stack_WGS84[["Anochetus_jonesi"]])
plot(Ponerinae_species_alpha_hull_stack_Mollweide[["Anochetus_jonesi"]])


plot(Ponerinae_species_alpha_hull_stack_WGS84[[1:9]])
plot(Ponerinae_species_alpha_hull_stack_Mollweide[[1:9]])

# Load data in RAM
Ponerinae_species_alpha_hull_stack_WGS84 <- readAll(Ponerinae_species_alpha_hull_stack_WGS84)
Ponerinae_species_alpha_hull_stack_Mollweide <- readAll(Ponerinae_species_alpha_hull_stack_Mollweide)

# Save final stacks
saveRDS(Ponerinae_species_alpha_hull_stack_WGS84, file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))
saveRDS(Ponerinae_species_alpha_hull_stack_Mollweide, file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_Mollweide.rds"))

### 2.6/ Compute species richness from alpha-hull ranges ####

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")

# Compute species richness as sum of binary ranges
Ponerinae_species_richness_WGS84 <- raster::calc(x = Ponerinae_species_alpha_hull_stack_WGS84, fun = sum, na.rm = T)
Ponerinae_species_richness_Mollweide <- raster::calc(x = Ponerinae_species_alpha_hull_stack_Mollweide, fun = sum, na.rm = T)

# Load data in RAM
Ponerinae_species_richness_WGS84 <- readAll(Ponerinae_species_richness_WGS84)
Ponerinae_species_richness_Mollweide <- readAll(Ponerinae_species_richness_Mollweide)

# Overlay to terrestrial areas
temp <- terrestrial_bg_WGS84
temp[Ponerinae_species_richness_WGS84@data@values > 0] <- Ponerinae_species_richness_WGS84@data@values[Ponerinae_species_richness_WGS84@data@values > 0]
Ponerinae_species_richness_WGS84 <- temp

temp <- terrestrial_bg_Mollweide
temp[Ponerinae_species_richness_Mollweide@data@values > 0] <- Ponerinae_species_richness_Mollweide@data@values[Ponerinae_species_richness_Mollweide@data@values > 0]
Ponerinae_species_richness_Mollweide <- temp

# Quick plot
plot(Ponerinae_species_richness_WGS84, col = pal_bl_red_Mannion)
plot(Ponerinae_species_richness_Mollweide, col = pal_bl_red_Mannion)

## Save raw species richness rasters
saveRDS(object = Ponerinae_species_richness_WGS84, file = "./outputs/Species_richness_maps/Ponerinae_species_richness_WGS84.rds")
saveRDS(object = Ponerinae_species_richness_Mollweide, file = "./outputs/Species_richness_maps/Ponerinae_species_richness_Mollweide.rds")

### 2.7/ Compute species richness per bioregions from Biogeographic table ####

# Include all taxa in the counts

Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

Ponerinae_sp_richness_per_bioregions_df <- Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA %>% 
  filter(Current_name %in% taxa_list_in_analyses)  %>%
  pivot_longer(cols = c(Afrotropics, Australasia, Indomalaya, Nearctic, Neotropics, `Eastern Palearctic`, `Western Palearctic`),
               names_to = "Bioregion") %>%
  group_by(Bioregion) %>%
  summarize(sp_richness_all_taxa = sum(value)) %>% 
  ungroup() 

# Update Bioregion_richness_data with counts from all taxa, including the one missing from the phylogeny

Bioregion_richness_data <- left_join(Bioregion_richness_data, Ponerinae_sp_richness_per_bioregions_df, by = join_by(bioregions == Bioregion)) %>% 
  mutate(node_labels_all_taxa = paste0(node_ID,"\n",sp_richness_all_taxa))

### 2.8/ Map species richness with ggplot ####

## 2.8.1/ Convert occurrence and bioregion data to Mollweide ####

# Occurrence data
Ponerinae_occ_df <- Biogeographic_database_Ponerinae_curated_no_duplicates %>% 
  filter(Current_name %in% taxa_list_in_analyses)
Ponerinae_occ_sf_WGS84 <- st_as_sf(x = Ponerinae_occ_df, coord = c("Longitude_dec", "Latitude_dec"), crs = "+proj=longlat +datum=WGS84")
Ponerinae_occ_sf_Mollweide <- st_transform(x = Ponerinae_occ_sf_WGS84, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

Mollweide_coordinates <- st_as_text(Ponerinae_occ_sf_Mollweide$geometry)
longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
Ponerinae_occ_sf_Mollweide$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
Ponerinae_occ_sf_Mollweide$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))

# Save occurrence sf files
saveRDS(object = Ponerinae_occ_sf_WGS84, file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_WGS84.rds")
saveRDS(object = Ponerinae_occ_sf_Mollweide, file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_Mollweide.rds")

# Bioregion data
Bioregion_richness_sf <- st_as_sf(x = Bioregion_richness_data, coords = c("longitude", "latitude"), remove = F, crs = st_crs(Ponerinae_occ_sf_WGS84))
Bioregion_richness_sf_Mollweide <- st_transform(Bioregion_richness_sf, crs = st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
Mollweide_coordinates <- st_as_text(Bioregion_richness_sf_Mollweide$geometry)
longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
Bioregion_richness_sf_Mollweide$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
Bioregion_richness_sf_Mollweide$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[Bioregion_richness_sf$node_ID]
# colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
# colors_list_for_areas <- colors_list_for_areas[Bioregion_richness_sf$node_ID]
names(colors_list_for_areas) <- Bioregion_richness_sf$bioregions

# Adjust order of bioregions
Bioregion_richness_sf$bioregions <- factor(Bioregion_richness_sf$bioregions, levels = Bioregion_richness_sf$bioregions, labels = Bioregion_richness_sf$bioregions)
Bioregion_richness_sf_Mollweide$bioregions <- factor(Bioregion_richness_sf_Mollweide$bioregions, levels = Bioregion_richness_sf_Mollweide$bioregions, labels = Bioregion_richness_sf_Mollweide$bioregions)

# Save Bioregion data files
saveRDS(object = Bioregion_richness_sf, file = "./outputs/Species_richness_maps/Bioregion_richness_sf.rds")
saveRDS(object = Bioregion_richness_sf_Mollweide, file = "./outputs/Species_richness_maps/Bioregion_richness_sf_Mollweide.rds")


## 2.7.2/ GGplot ####

# Load Bioregion data files
Bioregion_richness_sf_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregion_richness_sf_Mollweide.rds")

# Convert raster to SpatialPixelsDataFrame
Ponerinae_species_richness_WGS84_spdf <- as(Ponerinae_species_richness_WGS84, "SpatialPixelsDataFrame")
Ponerinae_species_richness_WGS84_spdf <- as.data.frame(Ponerinae_species_richness_WGS84_spdf)
colnames(Ponerinae_species_richness_WGS84_spdf) <- c("value", "x", "y")

Ponerinae_species_richness_Mollweide_spdf <- as(Ponerinae_species_richness_Mollweide, "SpatialPixelsDataFrame")
Ponerinae_species_richness_Mollweide_spdf <- as.data.frame(Ponerinae_species_richness_Mollweide_spdf)
colnames(Ponerinae_species_richness_Mollweide_spdf) <- c("value", "x", "y")

## Plot without node labels for species richness
Ponerinae_species_richness_ggplot <- ggplot(data = Bioregion_richness_sf_Mollweide) +
  
  # Plot species richness as raster background
  geom_tile(data = Ponerinae_species_richness_Mollweide_spdf,
            aes(x = x, y = y, fill = value), alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion) +
  
  # Plot bioregion sf maps
  geom_sf(data = Bioregions_sf_Bioregions_level_Mollweide,
          fill = NA,
          colour = "black",
          alpha = 0.0) +

  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +

  # Add title
  ggtitle(label =  paste0("Species richness of Ponerinae ants")) +
   
  # # Adjust legend aesthetics
  # guides(color = "none",
  #        size = "none",
  #        fill = guide_legend(order = 1)) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
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

pdf(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_richness_map.pdf"),
    width = 10, height = 5)

print(Ponerinae_species_richness_ggplot)

dev.off()


## Plot with node labels for species richness
Ponerinae_species_richness_with_labels_ggplot <- ggplot(data = Bioregion_richness_sf_Mollweide) +
  
  # Plot species richness as raster background
  geom_tile(data = Ponerinae_species_richness_Mollweide_spdf,
            aes(x = x, y = y, fill = value), alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion) +
  
  # Plot bioregion sf maps
  geom_sf(data = Bioregions_sf_Bioregions_level_Mollweide,
          fill = NA,
          colour = "black",
          alpha = 0.0) +
  
  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  # Plot nodes
  geom_point(data = Bioregion_richness_sf_Mollweide,
             mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide,
                           size = mean_counts,
                           # size = sp_richness_all_taxa,
                           col = bioregions),
             alpha = 0.5, show.legend = F) +

  # Adjust color legend
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +

  # Adjust size legend
  scale_size_continuous("Species richness",
                        range = c(5, 30)) +

  # Add node labels
  ggnewscale::new_scale(new_aes = "size") +
  geom_text(data = Bioregion_richness_sf_Mollweide,
            mapping = aes(label = round(mean_counts, 0),
                          # label = round(sp_richness_all_taxa, 0),
                          size = mean_counts,
                          # size = sp_richness_all_taxa,
                          x = longitude_Mollweide, y = latitude_Mollweide),
            # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
            color = "black",
            fontface = 2, show.legend = F) +

  # Adjust size for labels
  scale_size_continuous("Species\nrichness",
                        range = c(3, 10)) +

  # Add title
  ggtitle(label =  paste0("Species richness of Ponerinae ants")) +

  # # Adjust legend aesthetics
  # guides(color = "none",
  #        size = "none",
  #        fill = guide_legend(order = 1)) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
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

pdf(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_richness_map_with_labels.pdf"),
    width = 10, height = 5)

print(Ponerinae_species_richness_with_labels_ggplot)

dev.off()


##### 3/ Compute species richness map per genera ####

# Load stack of taxa ranges
Ponerinae_species_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))
Ponerinae_species_alpha_hull_stack_Mollweide <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_Mollweide.rds"))

# Load occurrence sf files
Ponerinae_occ_sf_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_WGS84.rds")
Ponerinae_occ_sf_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_Mollweide.rds")

# Extract all genera
Biogeographic_database_Ponerinae_curated_no_duplicates$Genus <- str_split(string = Biogeographic_database_Ponerinae_curated_no_duplicates$Current_name, pattern = "_", simplify = T)[,1]
taxa_list_in_analyses <- unique(Biogeographic_database_Ponerinae_curated_no_duplicates$Current_name[Biogeographic_database_Ponerinae_curated_no_duplicates$Current_status == "valid" | Biogeographic_database_Ponerinae_curated_no_duplicates$In_phylogeny])
genera_list <- unique(str_split(string = taxa_list_in_analyses, pattern = "_", simplify = T)[,1])
genera_list <- genera_list[order(genera_list)]

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")

### 3.1/ Compute species richness from alpha-hull ranges ####

Ponerinae_genera_sp_richness_WGS84 <- stack()
Ponerinae_genera_sp_richness_Mollweide <- stack()
## Loop per genera
for (i in seq_along(genera_list))
{
  # i <- 8
  
  # Extract genus
  genus_i <- genera_list[i]
  
  # Extract associated taxa
  taxa_i <- taxa_list_in_analyses[str_detect(string = str_split(string = taxa_list_in_analyses, pattern = "_", simplify = T)[,1], pattern = genus_i)]
  
  # Subset stack of ranges
  genus_alpha_hull_stack_WGS84_i <- raster::subset(Ponerinae_species_alpha_hull_stack_WGS84, subset = taxa_i)
  genus_alpha_hull_stack_Mollweide_i <- raster::subset(Ponerinae_species_alpha_hull_stack_Mollweide, subset = taxa_i)
  
  if (length(taxa_i) > 1)
  {
    # Compute species richness as sum of binary ranges
    genus_richness_WGS84_i <- raster::calc(x = genus_alpha_hull_stack_WGS84_i, fun = sum, na.rm = T)
    genus_richness_Mollweide_i <- raster::calc(x = genus_alpha_hull_stack_Mollweide_i, fun = sum, na.rm = T)
    
  } else {
    genus_richness_WGS84_i <- genus_alpha_hull_stack_WGS84_i
    genus_richness_Mollweide_i <- genus_alpha_hull_stack_Mollweide_i
  }
  
  # Load data in RAM
  genus_richness_WGS84_i <- readAll(genus_richness_WGS84_i)
  genus_richness_Mollweide_i <- readAll(genus_richness_Mollweide_i)
  
  # Overlay to terrestrial areas
  temp <- terrestrial_bg_WGS84
  temp[genus_richness_WGS84_i@data@values > 0] <- genus_richness_WGS84_i@data@values[genus_richness_WGS84_i@data@values > 0]
  genus_richness_WGS84_i <- temp
  
  temp <- terrestrial_bg_Mollweide
  temp[genus_richness_Mollweide_i@data@values > 0] <- genus_richness_Mollweide_i@data@values[genus_richness_Mollweide_i@data@values > 0]
  genus_richness_Mollweide_i <- temp
  
  # # Quick plot
  # plot(genus_richness_WGS84_i, col = pal_bl_red_Mannion)
  # plot(genus_richness_Mollweide_i, col = pal_bl_red_Mannion)
  
  # Store genera richness maps
  Ponerinae_genera_sp_richness_WGS84 <- addLayer(Ponerinae_genera_sp_richness_WGS84, genus_richness_WGS84_i)
  Ponerinae_genera_sp_richness_Mollweide <- addLayer(Ponerinae_genera_sp_richness_Mollweide, genus_richness_Mollweide_i)
  names(Ponerinae_genera_sp_richness_WGS84) <- genera_list[1:i]
  names(Ponerinae_genera_sp_richness_Mollweide) <- genera_list[1:i]
  
  # Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Species richness map created for Genus n°", i, "/", length(genera_list),"\n"))
  }
}

# Plot
plot(Ponerinae_genera_sp_richness_WGS84[[genera_list[1:9]]], col = pal_bl_red_Mannion)
plot(Ponerinae_genera_sp_richness_Mollweide[[genera_list[1:9]]], col = pal_bl_red_Mannion)

## Save stacks of species richness of Genera
saveRDS(object = Ponerinae_genera_sp_richness_WGS84, file = "./outputs/Species_richness_maps/Ponerinae_genera_sp_richness_WGS84.rds")
saveRDS(object = Ponerinae_genera_sp_richness_Mollweide, file = "./outputs/Species_richness_maps/Ponerinae_genera_sp_richness_Mollweide.rds")


### 3.2/ Compute bioregion richness data ####

# Load Bioregion data files
Bioregion_richness_sf_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregion_richness_sf_Mollweide.rds")

## 3.2.1/ Counts according to Historical biogeographic table

Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

Genera_sp_richness_per_bioregions_df <- Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA %>% 
  filter(Current_name %in% taxa_list_in_analyses)

Genera_sp_richness_per_bioregions_df$Genus <- str_split(string = Genera_sp_richness_per_bioregions_df$Current_name, pattern = "_", simplify = T)[,1]
Genera_sp_richness_per_bioregions_df <- Genera_sp_richness_per_bioregions_df %>%
  pivot_longer(cols = c(Afrotropics, Australasia, Indomalaya, Nearctic, Neotropics, `Eastern Palearctic`, `Western Palearctic`),
               names_to = "Bioregion") %>%
  group_by(Genus, Bioregion) %>%
  summarize(sp_richness = sum(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Bioregion, values_from = sp_richness)

# Save counts of species per Genus, per Bioregion
saveRDS(object = Genera_sp_richness_per_bioregions_df, file = "./outputs/Species_richness_maps/Genera_sp_richness_per_bioregions_df.rds")

Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")


## Loop per genera
Bioregion_genus_richness_sf_list <- list()
for (i in seq_along(genera_list))
{
  # i <- 3
  
  # Extract genus
  genus_i <- genera_list[i]
  
  # Extract genus richness
  genus_richness_i <- Genera_sp_richness_per_bioregions_df %>% 
    filter(Genus == genus_i) %>% 
    as.vector() %>%
    unlist()
  genus_richness_i <- as.numeric(genus_richness_i[-1])
  names(genus_richness_i) <- names(Genera_sp_richness_per_bioregions_df)[-1]
  
  # Add genus richness to Bioregion_richness_sf_Mollweide
  Bioregion_richness_sf_Mollweide_i <- Bioregion_richness_sf_Mollweide
  Bioregion_richness_sf_Mollweide_i$genus_richness <- genus_richness_i[match(x = Bioregion_richness_sf_Mollweide_i$bioregions, table = names(genus_richness_i))]
  
  # Adjust order of bioregions
  bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
  Bioregion_richness_sf_Mollweide_i$bioregions <- factor(Bioregion_richness_sf_Mollweide_i$bioregions, levels = bioregion_names, labels = bioregion_names)
  
  # Store in final list
  Bioregion_genus_richness_sf_list <- append(x = Bioregion_genus_richness_sf_list, values = list(Bioregion_richness_sf_Mollweide_i))
  names(Bioregion_genus_richness_sf_list) <- genera_list[1:i]
  
  # Print progress
  if (i %% 10 == 0)
  {
    # Save Bioregion data files
    saveRDS(object = Bioregion_genus_richness_sf_list, file = "./outputs/Species_richness_maps/Bioregion_genus_richness_sf_list.rds")
    
    cat(paste0(Sys.time(), " - Bioregion richness computed for Genus ",genus_i," - n°", i, "/", length(genera_list),"\n"))
  }
}



# ## 3.2.2/ Counts according to range maps
# 
# ## Loop per genera
# Bioregion_genus_richness_sf_list <- list()
# for (i in seq_along(genera_list))
# {
#   # i <- 3
#   
#   # Extract genus
#   genus_i <- genera_list[i]
#   
#   # Extract associated taxa
#   taxa_i <- taxa_list_in_analyses[str_detect(string = str_split(string = taxa_list_in_analyses, pattern = "_", simplify = T)[,1], pattern = genus_i)]
#   
#   # Subset stack of ranges
#   genus_alpha_hull_stack_WGS84_i <- raster::subset(Ponerinae_species_alpha_hull_stack_WGS84, subset = taxa_i)
#   
#   # # Extract data for each Bioregion
#   # genus_richness_sf_i <- raster::extract(x = genus_alpha_hull_stack_WGS84_i,
#   #                                        y = Bioregions_sf_Bioregions_level,
#   #                                        fun = max, na.rm = TRUE,
#   #                                        sp = TRUE)
#   
#   # Extract presence/absence per bioregions
#   genus_richness_data_i <- data.frame()
#   for (j in 1:nrow(Bioregions_sf_Bioregions_level))
#   {
#     # j <- 1
#     
#     # Extract bioregion
#     Bioregion_j <- Bioregions_sf_Bioregions_level$Bioregion[j]
#     
#     # Extract bioregion data
#     Bioregions_sf_j <- Bioregions_sf_Bioregions_level %>% 
#       filter(Bioregion == Bioregion_j)
#     
#     # Extract Genus richness data for Bioregion j
#     Bioregions_data_j <- raster::extract(x = genus_alpha_hull_stack_WGS84_i,
#                                          y = Bioregions_sf_j,
#                                          fun = max, na.rm = TRUE)
#     genus_richness_data_i <- rbind(genus_richness_data_i, Bioregions_data_j)
#   }
#   
#   # Compute bioregion richness
#   genus_richness_data_i <- rowSums(genus_richness_data_i)
#   names(genus_richness_data_i) <- Bioregions_sf_Bioregions_level$Bioregion
#   
#   # Add genus richness to Bioregion_richness_sf_Mollweide
#   Bioregion_richness_sf_Mollweide_i <- Bioregion_richness_sf_Mollweide
#   Bioregion_richness_sf_Mollweide_i$genus_richness <- genus_richness_data_i[match(x = Bioregion_richness_sf_Mollweide_i$bioregions, table = names(genus_richness_data_i))]
#   
#   # Adjust order of bioregions
#   bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
#   Bioregion_richness_sf_Mollweide_i$bioregions <- factor(Bioregion_richness_sf_Mollweide_i$bioregions, levels = bioregion_names, labels = bioregion_names)
#   
#   # Store in final list
#   Bioregion_genus_richness_sf_list <- append(x = Bioregion_genus_richness_sf_list, values = list(Bioregion_richness_sf_Mollweide_i))
#   names(Bioregion_genus_richness_sf_list) <- genera_list[1:i]
#   
#   # Print progress
#   if (i %% 1 == 0)
#   {
#     # Save Bioregion data files
#     saveRDS(object = Bioregion_genus_richness_sf_list, file = "./outputs/Species_richness_maps/Bioregion_genus_richness_sf_list.rds")
#     
#     cat(paste0(Sys.time(), " - Bioregion richness computed for Genus ",genus_i," - n°", i, "/", length(genera_list),"\n"))
#   }
# }
# 
# # Save Bioregion data files
# saveRDS(object = Bioregion_genus_richness_sf_list, file = "./outputs/Species_richness_maps/Bioregion_genus_richness_sf_list.rds")


### 3.3/ Map species richness per Genera with ggplot ####

# Load Bioregion data files
Bioregion_genus_richness_sf_list <- readRDS(file = "./outputs/Species_richness_maps/Bioregion_genus_richness_sf_list.rds")

## Adjust node range size

max_sp_richness <- max(apply(X = Genera_sp_richness_per_bioregions_df[,-1], MARGIN = 2, FUN = max))

# Function to adjust range size of nodes based on the range scale used for overall data
rescale_range_size <- function(x, data_min = 0, data_max = max_sp_richness, range_min = 5, range_max = 30) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}


# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[Bioregion_richness_sf$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[Bioregion_richness_sf$node_ID]
names(colors_list_for_areas) <- Bioregion_richness_sf$bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
colors_list_for_areas <- colors_list_for_areas[bioregion_names]

## Loop per genera
for (i in seq_along(genera_list))
{
  # i <- 3
  
  # Extract genus
  genus_i <- genera_list[i]
  
  # Extract Bioregion richness data
  Bioregion_genus_richness_sf_i <- Bioregion_genus_richness_sf_list[[i]]
  
  # Extract min/max node size to adjust range of points
  min_node_size <- min(Bioregion_genus_richness_sf_i$genus_richness)
  min_range_size <- rescale_range_size(x = min_node_size)
  max_node_size <- max(Bioregion_genus_richness_sf_i$genus_richness)
  max_range_size <- rescale_range_size(x = max_node_size)
  
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract Richness raster
  Genus_species_richness_raster_Mollweide_i <- raster::subset(Ponerinae_genera_sp_richness_Mollweide, subset = i)
  
  # Convert raster to SpatialPixelsDataFrame
  Genus_species_richness_Mollweide_spdf_i <- as(Genus_species_richness_raster_Mollweide_i, "SpatialPixelsDataFrame")
  Genus_species_richness_Mollweide_spdf_i <- as.data.frame(Genus_species_richness_Mollweide_spdf_i)
  colnames(Genus_species_richness_Mollweide_spdf_i) <- c("value", "x", "y")
  
  ### Run with and without labels
  
  # GGplot
  Genus_species_richness_ggplot <- ggplot(data = Bioregion_genus_richness_sf_i) +
    
    # Plot species richness as raster background
    geom_tile(data = Genus_species_richness_Mollweide_spdf_i,
              aes(x = x, y = y, fill = value), alpha = 1.0) +
    
    # Adjust color scheme and legend
    scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion) +
    
    # Plot bioregion sf maps
    geom_sf(data = Bioregions_sf_Bioregions_level_Mollweide,
            fill = NA,
            colour = "black",
            alpha = 0.0) +
    
    # Adjust CRS
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Plot nodes
    geom_point(data = Bioregion_genus_richness_sf_i,
               mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide,
                             size = genus_richness, col = bioregions),
               alpha = 0.5, show.legend = F) +

    # Adjust color legend
    scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +

    # Adjust size legend
    scale_size_continuous("Species richness",
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +

    # Add node labels
    ggnewscale::new_scale(new_aes = "size") +
    geom_text(data = Bioregion_genus_richness_sf_i,
              mapping = aes(label = round(genus_richness, 0),
                            x = longitude_Mollweide, y = latitude_Mollweide,
                            size = genus_richness),
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +

    # Adjust size for labels
    scale_size_continuous("Species\nrichness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
  
    # Add title
    ggtitle(label =  paste0("Species richness of ",genus_i," ants")) +
    
    # # Adjust legend aesthetics
    # guides(color = "none",
    #        size = "none",
    #        fill = guide_legend(order = 1)) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
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
  
  pdf(file = paste0("./outputs/Species_richness_maps/Genera_maps/",genus_i,"_species_richness_map.pdf"),
      width = 10, height = 5)
  
  print(Genus_species_richness_ggplot)
  
  dev.off()
  
  # Print progress
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Species richness mapped for Genus ",genus_i," - n°", i, "/", length(genera_list),"\n"))
  }
}


### 3.4/ Merge in a single PDF ####

all_ggplots_path <- list.files(path = "./outputs/Species_richness_maps/Genera_maps/", pattern = "_species_richness_map.pdf", full.names = T)
nb_ggplots <- length(all_ggplots_path) ; nb_ggplots

qpdf::pdf_combine(input = all_ggplots_path, output = paste0("./outputs/Species_richness_maps/all_Genera_species_richness_maps.pdf"))



##### 4/ Map occurrences per bioregions ####

# ## Load Bioregion sf map in Mollweide
# Bioregions_sf_Bioregions_level_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregions_sf_Bioregions_level_Mollweide.rds")

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

## Load occurrence data (in Mollweide)
Ponerinae_occ_sf_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_WGS84.rds")
Ponerinae_occ_sf_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_occ_sf_Mollweide.rds")

# Set color scheme for areas/bioregions (Use the BSM color scheme)
# colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
# areas_list <- c("A", "U", "E", "I", "R", "N", "W")
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[areas_list]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Adjust order of Bioregions
Bioregions_sf_Bioregions_level_Mollweide$Bioregion <- factor(Bioregions_sf_Bioregions_level_Mollweide$Bioregion, levels = bioregion_names, labels = bioregion_names)

plot(Bioregions_sf_Bioregions_level[, "Bioregion"])

table(Ponerinae_occ_sf_Mollweide$Bioregion_7_PaleA)

Ponerinae_occurrences_ggplot <- ggplot(data = Ponerinae_occ_sf_WGS84) +
  
  # Plot bioregion sf maps
  geom_sf(data = Bioregions_sf_Bioregions_level,
          fill = "grey95",
          colour = "black",
          alpha = 1.0) +
  
  # # Plot occurrences as sf object
  # geom_sf(data = Ponerinae_occ_sf_Mollweide,
  #         mapping = aes(fill = Bioregion),
  #         size = 2,
  #         colour = NA,
  #         alpha = 0.5) +
  
  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow/points outside of World map
           expand = FALSE) + # To avoid extended margins around map
  
  # Plot occurrences as data points
  geom_point(data = Ponerinae_occ_sf_Mollweide,
             mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide,
                           col = Bioregion_7_PaleA),
             size = 2, alpha = 0.5, show.legend = T) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Add title
  ggtitle(label =  paste0("Occurrences of Ponerinae ants")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(alpha = 1.0))) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(# panel.background = element_rect(fill = NA),
        panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_line(colour = "grey20", linetype = "dashed", linewidth = 0.2), # Plot graticules
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

pdf(file = paste0("./outputs/Species_richness_maps/Ponerinae_occurrences_map.pdf"),
    width = 10, height = 5)

print(Ponerinae_occurrences_ggplot)

dev.off()



##### 5/ Plot longitudinal gradient of species richness #####

longitude_scale <- seq(from = -180, to = 180, by = 1)

### 5.1/ Extract list of taxa with range encompassing longitude axes ####

# Load alpha-hull range data
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

## Loop per Longitudinal bands
longitude_binary_presence_df <- data.frame()
for (i in seq_along(longitude_scale))
{
  # i <- 1
  # i <- 70
  
  # Extract longitudinal band value
  longitude_i <- longitude_scale[i]
  
  # Convert into extent
  xmin <- max(-180, longitude_i - 0.5)
  xmax <- min(180, longitude_i + 0.5)
  longitude_extent_i <- raster::extent(c(xmin = xmin, xmax = xmax, ymin = -90, ymax = 90))
  
  # Crop alpha_hulls
  longitude_band_i <- raster::crop(x = Ponerinae_alpha_hull_stack_WGS84, y = longitude_extent_i)
  longitude_band_i <- raster::readAll(longitude_band_i)
  
  # Compute occurrence frequency per species
  taxa_freq_i <- apply(X = longitude_band_i@data@values, MARGIN = 2, FUN = sum, na.rm = T)
  
  # Identify taxa present in the longitudinal band
  taxa_binary_i <- taxa_freq_i > 0
  # table(taxa_binary_i)
  
  # Store result
  longitude_binary_presence_df <- rbind(longitude_binary_presence_df, taxa_binary_i)
  names(longitude_binary_presence_df) <- names(taxa_binary_i)
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Taxa presence detected for Longitude = ",longitude_i,"° - n°", i, "/", length(longitude_scale),"\n"))
  }
  
}
row.names(longitude_binary_presence_df) <- longitude_scale

# Remove alpha hull ranges to save RAM space
rm(Ponerinae_alpha_hull_stack_WGS84) ; gc() 

## Save df for taxa presence along longitudinal bands 
saveRDS(object = longitude_binary_presence_df, file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")


### 5.2/ Compute Species richness along longitudinal bands ####

# Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

longitudinal_bands_df <- data.frame(longitude_dec = longitude_scale)

longitudinal_bands_df$species_richness <- apply(X = longitude_binary_presence_df, MARGIN = 1, FUN = sum)

# Save longitudinal bands df
saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df.rds")

### 5.3/ Compute Species richness along longitudinal bands standardized by degrees of terrestrial lands ####

# Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

# Load terrestrial background
terrestrial_bg_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")

# Load longitudinal bands df
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df.rds")

## 5.3.1/ Compute degrees of terrestrial lands per longitudinal bands

## Loop per Longitudinal bands
terrestrial_land_degrees_per_longitudinal_bands <- c()
for (i in seq_along(longitude_scale))
{
  # i <- 1
  # i <- 70
  
  # Extract longitudinal band value
  longitude_i <- longitude_scale[i]
  
  # Convert into extent
  xmin <- max(-180, longitude_i - 0.5)
  xmax <- min(180, longitude_i + 0.5)
  longitude_extent_i <- raster::extent(c(xmin = xmin, xmax = xmax, ymin = -90, ymax = 90))
  
  # Crop terrestrial areas
  longitude_band_i <- raster::crop(x = terrestrial_bg_WGS84, y = longitude_extent_i)
  
  # Compute proportion of terrestrial lands
  terrestrial_prop_i <- sum(!is.na(longitude_band_i@data@values)) / length(longitude_band_i@data@values)
  
  # Convert to degrees
  terrestrial_degrees_i <- terrestrial_prop_i * 90
  
  # Store result
  terrestrial_land_degrees_per_longitudinal_bands <- c(terrestrial_land_degrees_per_longitudinal_bands, terrestrial_degrees_i)
  names(terrestrial_land_degrees_per_longitudinal_bands)[i] <- longitude_i
  
}
terrestrial_land_degrees_per_longitudinal_bands

# Inform longitudinal bands df
longitudinal_bands_df$terrestrial_degrees <- terrestrial_land_degrees_per_longitudinal_bands

## 5.3.2/ Compute standardized species richness

longitudinal_bands_df$species_richness_std <- round(longitudinal_bands_df$species_richness / longitudinal_bands_df$terrestrial_degrees, 1)

hist(longitudinal_bands_df$species_richness)
hist(longitudinal_bands_df$species_richness_std)

# Save longitudinal bands df
saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df.rds")


### 5.4/ Plot Species richness along longitudinal bands ####

# Load longitudinal bands df with SR data and median immigration ages
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df.rds")

# Create df with fake data for legend
SR_per_latitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                 y = c(50, 50),
                                                 data_type = c("Raw", "Std"))

hist(longitudinal_bands_df$species_richness)
hist(longitudinal_bands_df$species_richness_std)

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")

## 5.4.1/ GGplot with standardized richness ####

pdf(file = "./outputs/Species_richness_maps/SR_per_longitudinal_bands_ggplot_with_std.pdf", height = 8, width = 12)

SR_per_longitudinal_bands_ggplot <- ggplot(data = longitudinal_bands_df) +
  
  
  # Plot standardized species richness
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = species_richness_std * 20, x = longitude_dec),
            col = "dodgerblue",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot species richness
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = species_richness, x = longitude_dec),
            col = "orange",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = SR_per_latitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            size = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Raw", "Std"), labels = c("Raw", "Std"),
                     values = c("orange", "dodgerblue")) +
  
  # Set Y-axis
  scale_y_continuous(
    # Features of the first axis
    name = "Species richness",
    # Add a second axis and specify its transformation depending on the first axis
    sec.axis = sec_axis(transform = ~./20, name = "Standardized species richness")
  ) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Species richness across longitudes")) +
  
  # Set axes labels
  xlab("Longitude  [°]") +
  ylab("Species richness") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        legend.position.inside = c(0.16, 0.50),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y.left = element_text(margin = margin(r = 12)),
        axis.title.y.right =  element_text(margin = margin(l = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(SR_per_longitudinal_bands_ggplot)

dev.off()


## 5.4.2/ GGplot without standardized richness ####

SR_per_longitudinal_bands_ggplot <- ggplot(data = longitudinal_bands_df) +
  
  # Plot species richness as bar plot
  geom_col(data = longitudinal_bands_df,
           mapping = aes(y = species_richness, x = longitude_dec, fill = species_richness),
           show.legend = F,
           col = NA,
           width = 1.0,
           alpha = 1.0,
           linewidth = 0.0) +
  
  # Plot species richness as a line
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = species_richness, x = longitude_dec),
            col = "grey20",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion[1:180]) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Longitude", breaks = c(-120, -60, 0, 60, 120), labels = c("120°W", "60°W", "0°", "60°E", "120°E")) +
  
  # Set plot title +
  ggtitle(label = paste0("Species richness across longitudes")) +
  
  # Set axes labels
  xlab("Longitude  [°]") +
  ylab("Species richness") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y.left = element_text(margin = margin(r = 12)),
        axis.title.y.right =  element_text(margin = margin(l = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot

pdf(file = "./outputs/Species_richness_maps/SR_per_longitudinal_bands_ggplot.pdf", height = 8, width = 12)

print(SR_per_longitudinal_bands_ggplot)

dev.off()


SR_per_latitudinal_bands_ggplot <- ggplot(data = latitudinal_bands_df) +
  
  # Plot species richness as bar plot
  geom_col(data = latitudinal_bands_df,
           mapping = aes(x = species_richness, y = latitude_dec, fill = species_richness),
           show.legend = F,
           orientation = "y",
           col = NA,
           width = 1.0,
           alpha = 1.0,
           linewidth = 0.0) +
  
  # Plot species richness as line
  geom_line(data = latitudinal_bands_df,
            mapping = aes(x = species_richness, y = latitude_dec),
            orientation = "y",
            col = "grey20",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion[1:160]) +
  
  # Adjust label on Latitude axis
  scale_y_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  
  # Set plot title +
  ggtitle(label = paste0("Species richness across latitudes")) +
  
  # Set axes labels
  ylab("Latitude  [°]") +
  xlab("Species richness") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x.top = element_text(margin = margin(b = 12)),
        axis.title.x.bottom = element_text(margin = margin(t = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot

pdf(file = "./outputs/Species_richness_maps/SR_per_latitudinal_bands_ggplot.pdf", height = 8, width = 7)

print(SR_per_latitudinal_bands_ggplot)

dev.off()



##### 6/ Plot latitudinal gradient of species richness #####

latitude_scale <- seq(from = -60, to = 60, by = 1)

### 6.1/ Extract list of taxa with range encompassing longitude axes ####

# Load alpha-hull range data
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

## Loop per latitudinal bands
latitude_binary_presence_df <- data.frame()
for (i in seq_along(latitude_scale))
{
  # i <- 1
  # i <- 70
  
  # Extract latitudinal band value
  latitude_i <- latitude_scale[i]
  
  # Convert into extent
  ymin <- latitude_i - 0.5
  ymax <- latitude_i + 0.5
  latitude_extent_i <- raster::extent(c(xmin = -180, xmax = 180, ymin = ymin, ymax = ymax))
  
  # Crop alpha_hulls
  latitude_band_i <- raster::crop(x = Ponerinae_alpha_hull_stack_WGS84, y = latitude_extent_i)
  # latitude_band_i <- raster::readAll(latitude_band_i)
  
  # Compute occurrence frequency per species
  taxa_freq_i <- apply(X = latitude_band_i@data@values, MARGIN = 2, FUN = sum, na.rm = T)
  
  # Identify taxa present in the latitudinal band
  taxa_binary_i <- taxa_freq_i > 0
  # table(taxa_binary_i)
  
  # Store result
  latitude_binary_presence_df <- rbind(latitude_binary_presence_df, taxa_binary_i)
  names(latitude_binary_presence_df) <- names(taxa_binary_i)
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Taxa presence detected for Latitude = ",latitude_i,"° - n°", i, "/", length(latitude_scale),"\n"))
  }
  
}
row.names(latitude_binary_presence_df) <- latitude_scale

# Remove alpha hull ranges to save RAM space
rm(Ponerinae_alpha_hull_stack_WGS84) ; gc() 

## Save df for taxa presence along latitudinal bands 
saveRDS(object = latitude_binary_presence_df, file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")


### 6.2/ Compute Species richness along latitudinal bands ####

# Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

latitudinal_bands_df <- data.frame(latitude_dec = latitude_scale)

latitudinal_bands_df$species_richness <- apply(X = latitude_binary_presence_df, MARGIN = 1, FUN = sum)

# Save latitudinal bands df
saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df.rds")

### 6.3/ Compute Species richness along latitudinal bands standardized by degrees of terrestrial lands ####

# Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

# Load terrestrial background
terrestrial_bg_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")

# Load latitudinal bands df
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df.rds")

## 6.3.1/ Compute degrees of terrestrial lands per latitudinal bands

## Loop per latitudinal bands
terrestrial_land_degrees_per_latitudinal_bands <- c()
for (i in seq_along(latitude_scale))
{
  # i <- 1
  # i <- 70
  
  # Extract latitudinal band value
  latitude_i <- latitude_scale[i]
  
  # Convert into extent
  ymin <- latitude_i - 0.5
  ymax <- latitude_i + 0.5
  latitude_extent_i <- raster::extent(c(xmin = -180, xmax = 180, ymin = ymin, ymax = ymax))
  
  # Crop terrestrial areas
  latitude_band_i <- raster::crop(x = terrestrial_bg_WGS84, y = latitude_extent_i)
  
  # Compute proportion of terrestrial lands
  terrestrial_prop_i <- sum(!is.na(latitude_band_i@data@values)) / length(latitude_band_i@data@values)
  
  # Convert to degrees
  terrestrial_degrees_i <- terrestrial_prop_i * 90
  
  # Store result
  terrestrial_land_degrees_per_latitudinal_bands <- c(terrestrial_land_degrees_per_latitudinal_bands, terrestrial_degrees_i)
  names(terrestrial_land_degrees_per_latitudinal_bands)[i] <- latitude_i
  
}
terrestrial_land_degrees_per_latitudinal_bands

# Inform latitudinal bands df
latitudinal_bands_df$terrestrial_degrees <- terrestrial_land_degrees_per_latitudinal_bands

## 6.3.2/ Compute standardized species richness

latitudinal_bands_df$species_richness_std <- round(latitudinal_bands_df$species_richness / latitudinal_bands_df$terrestrial_degrees, 1)

hist(latitudinal_bands_df$species_richness)
hist(latitudinal_bands_df$species_richness_std)

# Save latitudinal bands df
saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df.rds")


### 6.4/ Plot Species richness along latitudinal bands ####

# Load latitudinal bands df with SR data and median immigration ages
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df.rds")

# Create df with fake data for legend
SR_per_latitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                 y = c(50, 50),
                                                 data_type = c("Raw", "Std"))

hist(latitudinal_bands_df$species_richness)
hist(latitudinal_bands_df$species_richness_std)

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")

## 6.4.1/ GGplot with standardized richness ####

pdf(file = "./outputs/Species_richness_maps/SR_per_latitudinal_bands_ggplot_with_std.pdf", height = 8, width = 7)

SR_per_latitudinal_bands_ggplot <- ggplot(data = latitudinal_bands_df) +
  
  
  # Plot standardized species richness
  geom_line(data = latitudinal_bands_df,
            mapping = aes(x = species_richness_std * 20, y = latitude_dec),
            orientation = "y",
            col = "dodgerblue",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot species richness
  geom_line(data = latitudinal_bands_df,
            mapping = aes(x = species_richness, y = latitude_dec),
            orientation = "y",
            col = "orange",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = SR_per_latitudinal_bands_legend_df,
            mapping = aes(y = x, x = y, col = data_type),
            size = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Raw", "Std"), labels = c("Raw", "Std"),
                     values = c("orange", "dodgerblue")) +
  
  # Set SR axes
  scale_x_continuous(
    # Features of the first axis
    name = "Species richness",
    # Add a second axis and specify its transformation depending on the first axis
    sec.axis = sec_axis(transform = ~./20, name = "Standardized species richness")
  ) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Species richness across latitudes")) +
  
  # Set axes labels
  ylab("Latitude  [°]") +
  xlab("Species richness") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        legend.position.inside = c(0.16, 0.50),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x.top = element_text(margin = margin(b = 12)),
        axis.title.x.bottom = element_text(margin = margin(t = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(SR_per_latitudinal_bands_ggplot)

dev.off()


## 6.4.2/ GGplot without standardized richness ####

SR_per_latitudinal_bands_ggplot <- ggplot(data = latitudinal_bands_df) +
  
  # Plot species richness as bar plot
  geom_col(data = latitudinal_bands_df,
           mapping = aes(x = species_richness, y = latitude_dec, fill = species_richness),
           show.legend = F,
           orientation = "y",
           col = NA,
           width = 1.0,
           alpha = 1.0,
           linewidth = 0.0) +
  
  # Plot species richness as line
  geom_line(data = latitudinal_bands_df,
            mapping = aes(x = species_richness, y = latitude_dec),
            orientation = "y",
            col = "grey20",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion[1:160]) +
  
  # Adjust label on Latitude axis
  scale_y_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  
  # Set plot title +
  ggtitle(label = paste0("Species richness across latitudes")) +
  
  # Set axes labels
  ylab("Latitude  [°]") +
  xlab("Species richness") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x.top = element_text(margin = margin(b = 12)),
        axis.title.x.bottom = element_text(margin = margin(t = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot

pdf(file = "./outputs/Species_richness_maps/SR_per_latitudinal_bands_ggplot.pdf", height = 8, width = 7)

print(SR_per_latitudinal_bands_ggplot)

dev.off()


