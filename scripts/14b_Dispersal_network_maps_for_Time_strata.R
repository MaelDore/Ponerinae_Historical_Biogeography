##### Script 14b: Create dispersal network maps per geological epochs with paleomaps  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Create dispersal network maps per geological epochs with paleomaps (not just current map as background)

###

### Inputs

# Paleomaps of bioregions
# Dispersal network with Bioregion richness

###

### Outputs

# PDFs and GIF of the dispersal network maps per geological epochs

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(rgdal)
library(sf)
library(sp)
library(tidyverse)
library(magick)
library(deeptime)

# Use the Paleomaps of the middle time
# Use the richness of the final time (max_counts)
# Use the sum of events within the window time for edges

### 1.2/ Load paleomaps ####

## Load Bioregion coastlines in time
Paleomaps_with_bioregions_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_with_bioregions_sf.rds")

## Extract time series
Paleomaps_time_series <- unlist(lapply(X = Paleomaps_with_bioregions_sf, FUN = function (x) { x$age[1] }))

### 1.3/ Load Time-strata info ####

# Get full labels of geological epochs
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")
time_strata_names <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cretaceous", "Early Cretaceous")

# Load geological epoch time boundaries
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
time_boundaries <- c(0, DEC_J_fit$inputs$timeperiods)

# Adjust earliest limit to the root age
Ponerinae_MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_short_names))
time_boundaries[length(time_boundaries)] <- root_age  

# Extract smoothed time for geological epochs as the rounded average age
strata_width <- time_boundaries[-1] - time_boundaries[-length(time_boundaries)]
smooth_time_per_strata <- time_boundaries[-length(time_boundaries)] + strata_width/2
names(smooth_time_per_strata) <- time_strata_names

### 1.4/ Load network data with fixed coordinates ####

## Load node metadata per time strata
# nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_time_strata.rds")
nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata.rds")

# Load edge metadata
# all_dispersal_events_all_strata_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))
all_dispersal_events_all_strata_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))

### 1.5/ Function to generate polgyons from data.frame ####

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


##### 2/ Find matches between time-strata and sliding windows #####

find_closest_value <- function(x, y)
{ 
  ID <- which(abs(y - x) == min(abs(y - x)))
  target <- y[first(ID)]
  return(target)
}

find_ID_closest_value <- function(x, y)
{ 
  ID <- which(abs(y - x) == min(abs(y - x)))
  return(first(ID))
}

# Window width of the sliding windows
window_width <- 5

# Load node metadata for sliding windows with updated coords
# nodes_metadata_per_sliding_windows_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))
nodes_metadata_per_sliding_windows_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))

# Extract smoothed time for sliding windows
smooth_time_per_sliding_windows <- unlist(lapply(X = nodes_metadata_per_sliding_windows_updated_coords, FUN = function (x) {unique(x$time)} ))

# Find the closest match for each time stratum
smooth_time_per_strata
smooth_time_per_strata_updated <- sapply(X = smooth_time_per_strata, FUN = find_closest_value, y = smooth_time_per_sliding_windows)
matching_ID <- sapply(X = smooth_time_per_strata, FUN = find_ID_closest_value, y = smooth_time_per_sliding_windows)

# Build df to convert from strata to sliding windows
strata_to_sliding_windows_df <- data.frame(stratum_ID = 1:length(time_strata_names),
                                           stratum_name = time_strata_names,
                                           stratum_full_name = time_strata_full_names,
                                           stratum_time = smooth_time_per_strata,
                                           sliding_window_ID = matching_ID,
                                           sliding_window_time = smooth_time_per_strata_updated)
strata_to_sliding_windows_df

# Save df to convert from time strata to sliding window
saveRDS(object = strata_to_sliding_windows_df, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/strata_to_sliding_windows_df.rds")


##### 3/ Update node metadata with moving coordinates and bioregion source labels ####

# Load node metadata for sliding windows with updated coords
# nodes_metadata_per_sliding_windows_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))
nodes_metadata_per_sliding_windows_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))

nodes_metadata_per_time_strata


# Initiate new object with updated coordinates for moving continents
nodes_metadata_per_time_strata_updated_coords <- nodes_metadata_per_time_strata
# Loop per time strata to update
for (i in 1:nrow(strata_to_sliding_windows_df))
{
  # i <- 1
  
  # Extract metadata of time stratum
  nodes_metadata_df_i <- nodes_metadata_per_time_strata[[i]]
  # Remove previous fixed coordinates
  nodes_metadata_df_i <- nodes_metadata_df_i %>% 
    dplyr::select(-latitude, -longitude)
  
  # Extract ID of matching sliding window
  matching_ID_i <- strata_to_sliding_windows_df$sliding_window_ID[i]
  
  # Extract metadata of associated sliding window
  matching_metadata_i <- nodes_metadata_per_sliding_windows_updated_coords[[matching_ID_i]]
  # Remove sliding window data
  matching_metadata_i <- matching_metadata_i %>% 
    dplyr::select(-mean_counts, -mean_percentages, -node_size, -node_labels)
  
  # Merge metadata
  nodes_metadata_df_i <- left_join(x = nodes_metadata_df_i, y = matching_metadata_i)
  
  # Store updated metadata
  nodes_metadata_per_time_strata_updated_coords[[i]] <- nodes_metadata_df_i
  
}
nodes_metadata_per_time_strata_updated_coords


# Save updated node metadata df
saveRDS(object = nodes_metadata_per_time_strata_updated_coords, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata_updated_coords.rds")


##### 4/ Update edge metadata with moving coordinates and bioregion source labels ####

# Load edge metadata of mean counts per sliding windows with updated coords
# all_dispersal_events_all_sliding_windows_df_list_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My_updated_coords.rds"))
all_dispersal_events_all_sliding_windows_df_list_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My_updated_coords.rds"))

# Initiate new object with updated coordinates for moving continents
all_dispersal_events_all_strata_df_list_updated_coords <- all_dispersal_events_all_strata_df_list
# Loop per time strata to update
for (i in 1:nrow(strata_to_sliding_windows_df))
{
  # i <- 1
  
  # Extract metadata of time stratum
  all_dispersal_events_df_i <- all_dispersal_events_all_strata_df_list[[i]]
  
  # Extract ID of matching sliding window
  matching_ID_i <- strata_to_sliding_windows_df$sliding_window_ID[i]
  
  # Extract metadata of associated sliding window
  matching_metadata_i <- all_dispersal_events_all_sliding_windows_df_list_updated_coords[[matching_ID_i]]
  # Remove sliding window data
  matching_metadata_i <- matching_metadata_i %>% 
    dplyr::select(-mean_counts, -mean_cum_counts, -mean_perc, -edge_width)
  
  # Merge metadata
  all_dispersal_events_df_i <- left_join(x = all_dispersal_events_df_i, y = matching_metadata_i)
  
  # Store updated metadata
  all_dispersal_events_all_strata_df_list_updated_coords[[i]] <- all_dispersal_events_df_i
  
}
all_dispersal_events_all_strata_df_list_updated_coords

# Save updated edge metadata df
saveRDS(object = all_dispersal_events_all_strata_df_list_updated_coords, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list_updated_coords.rds")


##### 5/ Plot dispersal network per time strata ####

### 5.1/ Load data ####

# Load df to convert from time strata to sliding window
strata_to_sliding_windows_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/strata_to_sliding_windows_df.rds")

## Load overall node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

# Load node metadata per time strata with updated coordinates
nodes_metadata_per_time_strata_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata_updated_coords.rds"))

## Load overall edge metadata
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")

# Load edge metadata of mean counts per time strata with updated coords
all_dispersal_events_all_strata_df_list_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list_updated_coords.rds"))


## Create polygons of extreme values to ensure all the bbox is plotted
E_patch <- data.frame(Latitude_dec = c(0, 0, 0.01, 0.01, 0),
                      Longitude_dec = c(179.9, 180.0, 180.0, 179.9, 179.9)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Australasia",
         Subregion = "Australasia",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  dplyr::select(-Subregion, -seq_ID)
st_crs(E_patch) <- st_crs(4326)

W_patch <- data.frame(Latitude_dec = c(0, 0, 0.01, 0.01, 0),
                      Longitude_dec = c(-179.9, -180.0, -180.0, -179.9, -179.9)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Australasia",
         Subregion = "Australasia",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  dplyr::select(-Subregion, -seq_ID)
st_crs(W_patch) <- st_crs(4326)

N_patch <- data.frame(Latitude_dec = c(90, 90, 89.9, 89.9, 90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Eastern Palearctic",
         Subregion = "Eastern Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  dplyr::select(-Subregion, -seq_ID)
st_crs(N_patch) <- st_crs(4326)

S_patch <- data.frame(Latitude_dec = c(-90, -90, -89.9, -89.9, -90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Antarctica",
         Subregion = "Antarctica",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  dplyr::select(-Subregion, -seq_ID)
st_crs(S_patch) <- st_crs(4326)


### 5.2/ Set color scheme ####

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica <- c("grey90", colors_list_for_areas)
bioregion_names_with_Antarctica <- c("Antarctica", bioregion_names)
names(colors_list_for_areas_with_Antarctica) <- bioregion_names_with_Antarctica


### 5.3/ Set size of nodes/edges ####

# Adjust minimal number of events to be displayed for each time strata
# min_counts_threshold <- 1
min_counts_threshold_list <- c(5, 3, 1, 1, 1, 1, 1)

# Extract min/max number of events to adjust scales
# max_node_counts <- max(unlist(lapply(X = nodes_metadata_per_time_strata_updated_coords, FUN = function (x) { x$max_counts } )))
# max_edge_counts <- max(unlist(lapply(X = all_dispersal_events_all_strata_df_list_updated_coords, FUN = function (x) { x$mean_counts } )))
max_node_counts <- max(nodes_metadata$mean_counts)
max_edge_counts <- max(all_dispersal_events_overall_df$mean_counts)

# Set breaks for edge linewidth
breaks_linewidth_list <- list(Stratum_1 = c(6, 20, 45), Stratum_2 = c(5, 15, 30), Stratum_3 = c(1, 2, 3), Stratum_4 = c(1, 2, 5), Stratum_5 = c(1, 2), Stratum_6 = c(1, 2), Stratum_7 = c(1))

# Function to adjust range size of nods based on the range scale used for overall data
rescale_range_size <- function(x, data_min = 0, data_max = max_node_counts, range_min = 5, range_max = 30) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}

# Function to adjust range linewidth of edges based on the range scale used for overall data
rescale_range_linewidth <- function(x, data_min = 5, data_max = max_edge_counts, range_min = 1, range_max = 6) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}


### 5.4/ Create geoscale ####

# Create fake time-serie data
root_to_tips_time_series <- root_age:0
ts <- rep(40, length(root_to_tips_time_series))
fake_ts_data_df <- data.frame(age = root_to_tips_time_series,
                              y = ts)
# Create custom epochs
epochs_custom <- data.frame(
  name = c("PPHo", "Miocene", "Oligo.", "Eocene", "Paleo.", "Late Cretaceous", "Early Cret."),
  max_age = c(5.3, 23.0, 33.9, 55.8, 66.0, 100.5, root_age),
  min_age = c(0, 5.3, 23.0, 33.9, 55.8, 66.0, 100.5),
  color = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"))

# Create ggplot of the geoscale
geol_scale_plot <- ggplot(fake_ts_data_df) +
  
  geom_line(aes(x = age, y = y), col = "white") +
  
  scale_x_reverse("Age (Mya)") +
  
  # geom_polygon(data = time_marker_polygon,
  #              mapping= aes(x = x, y = y),
  #              fill = "red",
  #              col = "black", lwd = 0.5) +
  
  ylab("") +
  
  deeptime::coord_geo(dat = list(epochs_custom), lwd = 0.3,
                      size = "auto", fontface = "bold",
                      height = unit(1.5, "line"),
                      xlim = c(root_age, 0), ylim = c(0, 40),
                      clip = "off") +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        plot.margin = margin(0.4, 3.5, 0.4, 2.2, "cm"), # trbl
        axis.line = element_line(linewidth = 1.0, color = NA),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(7, "pt"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, color = "black", face = "bold", margin = margin(t = 5, b = 5)),
        axis.title = element_text(size = 12, color = "black", face = "bold", margin = margin(t = 10, b = 5)))

### 5.5/ Plot ggplot for each time stratum ####

for (i in 1:nrow(strata_to_sliding_windows_df))
{
  # i <- 1
  
  # Extract name
  strata_full_name <- strata_to_sliding_windows_df$stratum_full_name[i]
  strata_name <- strata_to_sliding_windows_df$stratum_name[i]
  
  # Extract age
  age_i <- strata_to_sliding_windows_df$stratum_time[i]
  # Extract paleomap ID
  j <- which(Paleomaps_time_series == round(age_i, 0))
  
  ## 5.5.1/ Extract data for time i and adjust display scales ####
  
  # Extract the minimum number of events to display an edge
  # min_counts_threshold_i <- min_counts_threshold
  min_counts_threshold_i <- min_counts_threshold_list[i]
  
  # Extract node metadata
  nodes_metadata_time_strata_i <- nodes_metadata_per_time_strata_updated_coords[[i]]
  # Extract max node size to adjust range of points
  # max_node_size <- max(nodes_metadata_time_strata_i$mean_counts)
  max_node_size <- max(nodes_metadata_time_strata_i$max_counts)
  max_range_size <- rescale_range_size(x = max_node_size)
  # Extract min node size to adjust range of points
  # min_node_size <- min(nodes_metadata_time_strata_i$mean_counts)
  min_node_size <- min(nodes_metadata_time_strata_i$max_counts)
  min_range_size <- rescale_range_size(x = min_node_size)
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract edge metadata
  all_dispersal_events_strata_i_df <- all_dispersal_events_all_strata_df_list_updated_coords[[i]]
  # Extract max node size to adjust range of linewidth
  max_edge_linewidth <- max(all_dispersal_events_strata_i_df$mean_counts)
  max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  # Extract min node size to adjust range of linewidth
  min_edge_linewidth <- min(all_dispersal_events_strata_i_df$mean_counts[all_dispersal_events_strata_i_df$mean_counts >= min_counts_threshold_i])
  min_range_linewidth <- rescale_range_linewidth(x = min_edge_linewidth)
  # Extract manual breaks for edge linewidth per stratum
  # breaks_linewidth_i <- breaks_linewidth_list[[i]]
  breaks_linewidth_i <- c(max_edge_linewidth, max_edge_linewidth, max_edge_linewidth) # Provide 3 valid entries that will be replaced in the legend
  
  # If min counts of edges = max counts of edges, need a fix by adding 0.001 to an edge count
  if (min_edge_linewidth == max_edge_linewidth)
  {
    min_egde_indices <- which(all_dispersal_events_strata_i_df$mean_counts == min_edge_linewidth)
    edge_fix_index <- sample(x = min_egde_indices, size = 1)
    all_dispersal_events_strata_i_df$mean_counts[edge_fix_index] <- all_dispersal_events_strata_i_df$mean_counts[edge_fix_index] + 0.001
    max_edge_linewidth <- max(all_dispersal_events_strata_i_df$mean_counts)
    max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  }
  
  # If number of edges to display is 1, create duplicates of that edge, and add 0.001 to the edge count
  min_egde_indices <- which(all_dispersal_events_strata_i_df$mean_counts >= min_counts_threshold_i)
  if (length(min_egde_indices) == 1)
  {
    fix_df <- all_dispersal_events_strata_i_df[min_egde_indices, ]
    fix_df$mean_counts <- fix_df$mean_counts + 0.001
    all_dispersal_events_strata_i_df <- rbind(all_dispersal_events_strata_i_df, fix_df)
    breaks_linewidth_i <- c(fix_df$mean_counts, fix_df$mean_counts, fix_df$mean_counts)
  }
  # If number of egde to display is null, created two loop edges with different edge values and set linewidth range to zero
  if (length(min_egde_indices) == 0)
  {
    fix_df <- all_dispersal_events_strata_i_df[c(1,9), ]
    fix_df$mean_counts <- c(1, 1.001)
    all_dispersal_events_strata_i_df <- rbind(all_dispersal_events_strata_i_df, fix_df)
    breaks_linewidth_i <- c(1, 1, 1)
  }
  
  # Extract Paleomaps for age j
  Paleomaps_with_bioregions_sf_i <- Paleomaps_with_bioregions_sf[[j]]
  
  # Copy metadata from 1st sf line (except Bioregion)
  metadata <- st_drop_geometry(Paleomaps_with_bioregions_sf_i[1, ]) %>% 
    dplyr::select(-Bioregion)
  
  # Add E/W patches to ensure all the bbox is plotted ans stable
  E_patch_to_bind <- cbind(E_patch, metadata)
  W_patch_to_bind <- cbind(W_patch, metadata)
  N_patch_to_bind <- cbind(N_patch, metadata)
  S_patch_to_bind <- cbind(S_patch, metadata)
  # Add patches to the list of bioregions
  Paleomaps_with_bioregions_sf_i <- rbind(Paleomaps_with_bioregions_sf_i, E_patch_to_bind, W_patch_to_bind)
  Paleomaps_with_bioregions_sf_i <- rbind(Paleomaps_with_bioregions_sf_i, N_patch_to_bind, S_patch_to_bind)
  
  # # Crop to final fixed bbox (does not fix the wobbling...)
  # suppressMessages(sf_use_s2(FALSE))
  # Paleomaps_with_bioregions_sf_i <- suppressMessages(st_crop(x = Paleomaps_with_bioregions_sf_i, y = st_bbox(obj = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))))
  # suppressMessages(sf_use_s2(TRUE))
  
  ## 5.5.2/ Adjust color scheme ####
  
  # Adjust bioregions factors in node/edge data
  nodes_metadata_time_strata_i$bioregions <- factor(nodes_metadata_time_strata_i$bioregions, levels = bioregion_names, labels = bioregion_names)
  all_dispersal_events_strata_i_df$source_labels <- factor(all_dispersal_events_strata_i_df$source_labels, levels = bioregion_names, labels = bioregion_names)
  
  # Adjust bioregions factors in Paleomaps
  Paleomaps_with_bioregions_sf_i$Bioregion <- factor(Paleomaps_with_bioregions_sf_i$Bioregion, levels = bioregion_names_with_Antarctica, labels = bioregion_names_with_Antarctica)
  
  ## 5.5.3/ Plot map and network ####
  
  # Plot PDF
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata_moving_continents/all_dispersal_events_mean_counts_time_stratum_",i,"_moving_continents_ggplot.pdf"),
      width = 9, height = 6)
  
  all_dispersal_events_time_strata_i_ggplot <- ggplot(data = nodes_metadata_time_strata_i) +
    
    # Plot bioregion maps
    geom_sf(data = Paleomaps_with_bioregions_sf_i,
            mapping = aes(fill = Bioregion),
            colour = "black",
            alpha = 1.0) +
    
    # Adjust fill color scheme for bioregions
    scale_fill_manual("Bioregions", breaks = bioregion_names_with_Antarctica, labels = bioregion_names_with_Antarctica, values = colors_list_for_areas_with_Antarctica) +
    
    # Adjust CRS
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Plot nodes
    geom_point(data = nodes_metadata_time_strata_i,
               mapping = aes(x = Longitude_Mollweide, y = Latitude_Mollweide,
                             # size = mean_counts),
                             size = max_counts),
               alpha = 0.3, show.legend = F) +
    
    # Adjust size legend
    scale_size_continuous("Species richness", 
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +
    
    # Plot shades of vertices/edges
    geom_curve(data = all_dispersal_events_strata_i_df[all_dispersal_events_strata_i_df$mean_counts >= min_counts_threshold_i, ],
               aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                   xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                   linewidth = mean_counts*1.5),
               color = "black",
               arrow = arrow(length = unit(0.03, "npc"), 
                             type = "closed"), # Describes arrow head (open or closed)
               angle = 90, # Anything other than 90 or 0 can look unusual
               alpha = 0.3) + 
    
    # Plot vertices/edges for dispersal
    geom_curve(data = all_dispersal_events_strata_i_df[all_dispersal_events_strata_i_df$mean_counts >= min_counts_threshold_i, ],
               aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                   xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                   color = source_labels, linewidth = mean_counts),
               arrow = arrow(length = unit(0.03, "npc"), 
                             type = "closed"), # Describes arrow head (open or closed)
               angle = 90, # Anything other than 90 or 0 can look unusual
               alpha = 1.0) + 
    
    # Adjust color scheme for arrows
    scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
    
    # Adjust edge width legend
    scale_linewidth_continuous("Dispersal events",
                               breaks = breaks_linewidth_i,
                               labels = c("1", "15", "40"),
                               # range = c(5, 30),
                               range = c(min_range_linewidth, max_range_linewidth*1.3)) +
    
    # Add node labels
    ggnewscale::new_scale(new_aes = "size") +
    geom_text(data = nodes_metadata_time_strata_i,
              mapping = aes(# label = round(mean_counts, 0),
                            label = round(max_counts, 0),
                            x = Longitude_Mollweide, y = Latitude_Mollweide,
                            # size = mean_counts),
                            size = max_counts),
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +
    
    # Adjust size for labels
    scale_size_continuous("Species richness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
    
    # Add title
    ggtitle(label =  paste0("Dispersal events between bioregions\n", strata_full_name)) +
    
    # Adjust legend aesthetics
    guides(color = "none",
           size = "none",
           fill = guide_legend(order = 1),
           linewidth = guide_legend(order = 2,
                                    override.aes = list(linewidth = c(0.77, 1.55, 2.95), # Values for 1, 15, 40
                                                        arrow(length = unit(0.2, "cm"))))
    ) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.0, 0.5, "cm"), # trbl
          # legend.position = c(1.1, 0.9),
          legend.box.margin = margin(l = 15, b = 15),
          legend.key.size = unit(1.0, 'lines'),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  # print(all_dispersal_events_time_strata_i_ggplot)
  
  ## 5.5.4/ Plot geoscale ####
  
  ## Add time marker
  current_time <- age_i
  
  time_marker_polygon_i <- data.frame(
    x = c(current_time-1, current_time, current_time+1, current_time-1),
    y = c(100, 0, 100, 100))
  
  geol_scale_i_plot <- geol_scale_plot +
    
    geom_polygon(data = time_marker_polygon_i,
                 mapping= aes(x = x, y = y),
                 fill = "red",
                 col = "black", lwd = 0.5)
  
  # print(geol_scale_i_plot)
  
  ## 5.5.5/ Arrange map and scale in a faceted plot ####
  
  gridExtra::grid.arrange(
    grobs = list(all_dispersal_events_time_strata_i_ggplot,
                 NULL,
                 geol_scale_i_plot),
    widths = 1,  # Width of columns
    heights = c(10, 0.1, 2.8),
    nrow = 3,
    ncol = 1)
  
  dev.off()
  
  ### Print progress
  if (i %% 0 == 0)
  {
    cat(paste0(Sys.time(), " - Map created for ",strata_full_name," - Stratum n°", i, "/", nrow(strata_to_sliding_windows_df),"\n"))
  }
}


##### 6/ Aggregate plots ####


## 6.1/ Aggregate all ggplot plots in a single pdf ####

all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata_moving_continents/", pattern = "all_dispersal_events_mean_counts_time_stratum_", full.names = T)

all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = "ggplot.pdf")]
nb_ggplots <- length(all_ggplots_path)

qpdf::pdf_combine(input = all_ggplots_path, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_counts_all_strata_moving_continents_ggplot.pdf"))

## 6.2/ Aggregate all ggplot in forward timeline + overall ####

overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_updated_ggplot.pdf"

qpdf::pdf_combine(input = c(rev(all_ggplots_path), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_counts_all_strata_moving_continents_ggplot_forward.pdf"))


### 6.3/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

# pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"),
pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"),
                                      pages = NULL, density = 100)
magick::image_info(pdf_pointer)

image_resize_and_write_gif(image = pdf_pointer,
                           # path =  paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.gif"),
                           path =  paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.gif"),
                           delay = 2, # Time between frames in seconds
                           width = 900, height = 400,
                           loop = FALSE,
                           progress = TRUE)


## Add arrow to legend (Does not work with guides. Need to add on PPT)
## Add bg color on sphere (does not clip out of the sphere... Looks better without grid and bg)

# Add Color contour in igraphs on areas that are adjacent at a given time: Section 4.3 and 4.4
# Adjust coordinates in each stratum to show continental drift?

# Create bioregion shapefile that can overlay the PALEOMAPS?

### Make an animation by plotting every geological epoch?
# Use gganimate for paths to make them fade in_out?

