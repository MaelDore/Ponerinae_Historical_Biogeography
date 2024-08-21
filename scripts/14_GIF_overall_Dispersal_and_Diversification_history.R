##### Script 14: Create GIF of overall Dispersal and Diversification history  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Create a GIF of the overall Dispersal and Diversification history of Ponerinae

###

### Inputs

# Paleomaps of bioregions
# Dispersal network with Bioregion richness

###

### Outputs

# GIF of the overall Dispersal and Diversification history of Ponerinae

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


### 1.2/ Load paleomaps ####

## Load Bioregion coastlines in time
Paleomaps_with_bioregions_sf <- readRDS(file = "./outputs/Paleomaps/Paleomaps_with_bioregions_sf.rds")

## Extract time series
Paleomaps_time_series <- unlist(lapply(X = Paleomaps_with_bioregions_sf, FUN = function (x) { x$age[1] }))

# ## Load the reconstructed Bioregion hulls per ages
# Bioregion_hulls_reconstruct_all_ages <- readRDS(file = "./outputs/Paleomaps/Bioregion_hulls_reconstruct_all_ages.rds")
# # Need to extract centroids


### 1.3/ Define the sliding windows/time series ####

# Set sliding window characteristics
window_width <- 5 # Width = 5 My (Must be a multiple of window_steps)
window_steps <- 1 # Steps = 1 My (Must be a divider of window_width)

# Extract root age
# root_age <- max(DEC_J_clado_events_tables[[1]]$time_bp)
root_age <- 113

# Generate start and end times of sliding windows
start_time_list <- seq(from = 0, to = root_age - window_width, by = window_steps)
end_time_list <- seq(from = window_width, to = root_age, by = window_steps)
mean_time_list <- (start_time_list + end_time_list) / 2

nodes_metadata_time_series <- floor(mean_time_list)


### 1.4/ Load network data ####

## Load overall node metadata
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/nodes_metadata.rds")
## Load node metadata per sliding windows
nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))

## Load overall edge metadata
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/all_dispersal_events_overall_df.rds")
# Load edge metadata per sliding windows
all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))


### 1.5 Adjust color scheme for a less shiny output ####

colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]

original_BSM_color_scheme <- c("#0000FF", "#00FFFF", "#00CD00", "#FFFF00", "#FFA500", "#FF0000", "#EE82EE")
new_lighter_color_scheme <- c("royalblue3", "darkslategray2", "limegreen", "#ffff42", "orange1", "firebrick2", "plum1")
names(new_lighter_color_scheme) <- c("A", "U", "I", "R", "N", "E", "W")

saveRDS(new_lighter_color_scheme, file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")

### 1.6/ Function to generate polgyons from data.frame ####

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


##### 2/ Compute centroids of bioregions through time #####

# Create East/West hemisphere pôlygons
West_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = -180, xmax = 0, ymax = 90, ymin = -90), crs = st_crs(4326))))
East_hemisphere <- st_as_sf(x = st_as_sfc(x = st_bbox(c(xmin = 0, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))))
EW_polygons <- rbind(East_hemisphere, West_hemisphere) %>% 
  rename(geometry = x) %>% 
  mutate(Hemisphere = c("East", "West"))

Bioregion_centroids_sf <- Paleomaps_with_bioregions_sf
## Loop per ages
for (i in seq_along(Paleomaps_with_bioregions_sf))
{
  # i <- 1
  
  # Extract age
  age_i <- Paleomaps_time_series[i]
  
  # Extract paleomaps of bioregions
  Bioregion_coastlines_i  <- Paleomaps_with_bioregions_sf[[i]]
  # plot(Bioregion_coastlines_i[, "Bioregion"])
  
  ### 2.1/ Compute centroid objects
  
  # Get centroids
  sf_use_s2(FALSE)
  Bioregion_centroids_i <- Bioregion_coastlines_i %>% 
    st_centroid()
  sf_use_s2(TRUE)
  
  ### 2.2/ Deal with the special case of Australasia and Eastern Palearctic (Centroids should account only for the Eastern Hemisphere)
  
  # Split in two Hemispheres
  sf_use_s2(FALSE)
  Bioregion_centroids_split_i <- st_intersection(x = Bioregion_coastlines_i, y = EW_polygons) %>%
    st_centroid() %>%
    filter(Bioregion %in% c("Australasia", "Eastern Palearctic")) %>%
    filter(Hemisphere == "East") %>%
    select(-Hemisphere)
  sf_use_s2(TRUE)
  
  # Replace centroids with East centroids for Australasia and Eastern Palearctic
  Bioregion_centroids_i <- Bioregion_centroids_i %>% 
    filter(!(Bioregion %in% c("Australasia", "Eastern Palearctic")))
  Bioregion_centroids_i <- rbind(Bioregion_centroids_i, Bioregion_centroids_split_i)
  
  ### 2.3/ Extract WGS84 coordinates 
  
  Bioregion_centroids_sp_i <- as(Bioregion_centroids_i, 'Spatial')
  WGS_coords_i <- as.data.frame(round(Bioregion_centroids_sp_i@coords, digits = 1))
  names(WGS_coords_i) <- c("Longitude_WGS84", "Latitude_WGS84")
  
  ### 2.4/ Convert to Mollweide and extract coordinates
  Bioregion_centroids_Mollweide_i <- Bioregion_centroids_i %>% 
    st_transform(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  Bioregion_centroids_Mollweide_sp_i <- as(Bioregion_centroids_Mollweide_i, 'Spatial')
  Mollweide_coords_i <- as.data.frame(round(Bioregion_centroids_Mollweide_sp_i@coords, digits = 0))
  names(Mollweide_coords_i) <- c("Longitude_Mollweide", "Latitude_Mollweide")
  
  # Inform sf object
  Bioregion_centroids_i <- cbind(Bioregion_centroids_i, WGS_coords_i, Mollweide_coords_i)
  
  # Store results
  Bioregion_centroids_sf[[i]] <- Bioregion_centroids_i
}


##### 3/ Inform nodes and edges medatata with centroid data ####

### 3.1/ Inform node metadata ####

## For the sliding window
nodes_metadata_per_sliding_windows_updated_coords <- nodes_metadata_per_sliding_windows
## Loop per ages
for (i in seq_along(nodes_metadata_time_series))
{
  # i <- 1
  
  # Extract age
  age_i <- nodes_metadata_time_series[i]
  j <- which(Paleomaps_time_series == age_i)
  
  # Extract node metadata
  nodes_metadata_per_sliding_windows_i  <- nodes_metadata_per_sliding_windows[[i]]
  # Extract centroid data
  Bioregion_centroids_i <- Bioregion_centroids_sf[[j]] %>% 
    st_drop_geometry() %>%
    rename(bioregions = Bioregion) %>%
    select(bioregions,
           Longitude_WGS84, Latitude_WGS84,
           Longitude_Mollweide, Latitude_Mollweide)
  
  # Remove previous coordinates info
  nodes_metadata_per_sliding_windows_i <- nodes_metadata_per_sliding_windows_i %>% 
    select(-latitude, -longitude) %>% 
    left_join(Bioregion_centroids_i)
  
  # Store updated results
  nodes_metadata_per_sliding_windows_updated_coords[[i]] <- nodes_metadata_per_sliding_windows_i
}

# Save node metadata with updated coords
saveRDS(object = nodes_metadata_per_sliding_windows_updated_coords, file = paste0("./outputs/Network_events_counts/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))


## For T = 0 overall total

# Extract centroid data
Bioregion_centroids_T0 <- Bioregion_centroids_sf[[1]] %>% 
  st_drop_geometry() %>%
  rename(bioregions = Bioregion) %>%
  select(bioregions,
         Longitude_WGS84, Latitude_WGS84,
         Longitude_Mollweide, Latitude_Mollweide)

# Remove previous coordinates info
nodes_metadata_updated_coords <- nodes_metadata %>% 
  select(-latitude, -longitude) %>% 
  left_join(Bioregion_centroids_T0)

# Save node metadata for T=0 with updated coords
saveRDS(object = nodes_metadata_updated_coords, file = paste0("./outputs/Network_events_counts/nodes_metadata_updated_coords.rds"))


### 3.2/ Inform edge metadata ####

# Conversion df from Bioregion labels to CAPITALS
bioregion_conversion_to_labels_df <- data.frame(area = c("A", "U", "I", "R", "N", "E", "W"),
                                                bioregion = c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern_Palearctic", "Western_Palearctic"),
                                                label = c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic"))

## For the sliding window
all_dispersal_events_all_sliding_windows_df_list_updated_coords <- all_dispersal_events_all_sliding_windows_df_list
## Loop per ages
for (i in seq_along(nodes_metadata_time_series))
{
  # i <- 1
  
  # Extract age
  age_i <- nodes_metadata_time_series[i]
  j <- which(Paleomaps_time_series == age_i)
  
  # Extract edge metadata
  all_dispersal_events_all_sliding_windows_df_i  <- all_dispersal_events_all_sliding_windows_df_list[[i]]
  # Extract centroid data
  Bioregion_centroids_i <- Bioregion_centroids_sf[[j]] %>% 
    st_drop_geometry() %>%
    left_join(y = bioregion_conversion_to_labels_df[, c("label", "area")], by = join_by("Bioregion" == "label")) %>%
    select(area,
           Longitude_WGS84, Latitude_WGS84,
           Longitude_Mollweide, Latitude_Mollweide)
  
  # Update coordinates info for the source
  all_dispersal_events_all_sliding_windows_df_i <- all_dispersal_events_all_sliding_windows_df_i %>% 
    select(-source_latitude, -source_longitude, -source_latitude_Mollweide, -source_longitude_Mollweide) %>% 
    left_join(Bioregion_centroids_i, by = join_by("source" == "area")) %>% 
    rename(source_longitude = Longitude_WGS84,
           source_latitude = Latitude_WGS84,
           source_longitude_Mollweide = Longitude_Mollweide,
           source_latitude_Mollweide = Latitude_Mollweide)
  
  # Update coordinates info for the dest(ination)
  all_dispersal_events_all_sliding_windows_df_i <- all_dispersal_events_all_sliding_windows_df_i %>% 
    select(-dest_latitude, -dest_longitude, -dest_latitude_Mollweide, -dest_longitude_Mollweide) %>% 
    left_join(Bioregion_centroids_i, by = join_by("dest" == "area")) %>% 
    rename(dest_longitude = Longitude_WGS84,
           dest_latitude = Latitude_WGS84,
           dest_longitude_Mollweide = Longitude_Mollweide,
           dest_latitude_Mollweide = Latitude_Mollweide)
  
  # Store updated results
  all_dispersal_events_all_sliding_windows_df_list_updated_coords[[i]] <- all_dispersal_events_all_sliding_windows_df_i
}

# Save edge metadata with updated coords
saveRDS(object = all_dispersal_events_all_sliding_windows_df_list_updated_coords, file = paste0("./outputs/Network_events_counts/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My_updated_coords.rds"))

## For T = 0 overall total

# Extract centroid data
Bioregion_centroids_T0 <- Bioregion_centroids_sf[[1]] %>% 
  st_drop_geometry() %>%
  left_join(y = bioregion_conversion_to_labels_df[, c("label", "area")], by = join_by("Bioregion" == "label")) %>%
  select(area,
         Longitude_WGS84, Latitude_WGS84,
         Longitude_Mollweide, Latitude_Mollweide)

# Update coordinates info for the source
all_dispersal_events_overall_df_updated_coords <- all_dispersal_events_overall_df %>% 
  select(-source_latitude, -source_longitude, -source_latitude_Mollweide, -source_longitude_Mollweide) %>% 
  left_join(Bioregion_centroids_T0, by = join_by("source" == "area")) %>% 
  rename(source_longitude = Longitude_WGS84,
         source_latitude = Latitude_WGS84,
         source_longitude_Mollweide = Longitude_Mollweide,
         source_latitude_Mollweide = Latitude_Mollweide)

# Update coordinates info for the dest(ination)
all_dispersal_events_overall_df_updated_coords <- all_dispersal_events_overall_df_updated_coords %>% 
  select(-dest_latitude, -dest_longitude, -dest_latitude_Mollweide, -dest_longitude_Mollweide) %>% 
  left_join(Bioregion_centroids_T0, by = join_by("dest" == "area")) %>% 
  rename(dest_longitude = Longitude_WGS84,
         dest_latitude = Latitude_WGS84,
         dest_longitude_Mollweide = Longitude_Mollweide,
         dest_latitude_Mollweide = Latitude_Mollweide)

# Save edge metadata for T=0 with updated coords
saveRDS(object = all_dispersal_events_overall_df_updated_coords, file = paste0("./outputs/Network_events_counts/all_dispersal_events_overall_df_updated_coords.rds"))


##### 4/ Map overall network for T=0 #####

### 4.1/ Prepare stuff ####

## Load overall node metadata for T=0 with updated coordinates
nodes_metadata_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/nodes_metadata_updated_coords.rds"))

## Load overall edge metadata for T=0 with updated coords
all_dispersal_events_overall_df_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/all_dispersal_events_overall_df_updated_coords.rds"))

# Set minimum nb of dispersal events to display
min_counts_threshold <- 5

# Set color scheme for areas/bioregions (Use the BSM color scheme)
# colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica <- c("grey90", colors_list_for_areas)
bioregion_names_with_Antarctica <- c("Antarctica", bioregion_names)
names(colors_list_for_areas_with_Antarctica) <- bioregion_names_with_Antarctica

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

N_patch <- data.frame(Latitude_dec = c(90, 90, 89.9, 89.9, 90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Eastern Palearctic",
         Subregion = "Eastern Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(N_patch) <- st_crs(4326)

S_patch <- data.frame(Latitude_dec = c(-90, -90, -89.9, -89.9, -90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Antarctica",
         Subregion = "Antarctica",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(S_patch) <- st_crs(4326)


## Extract Paleomaps for T=0

Paleomaps_with_bioregions_sf_T0 <- Paleomaps_with_bioregions_sf[[1]]

# Copy metadata from 1st sf line (except Bioregion)
metadata <- st_drop_geometry(Paleomaps_with_bioregions_sf_T0[1, ]) %>% 
  select(-Bioregion)

# Add E/W patches to ensure all the bbox is plotted ans stable
E_patch_to_bind <- cbind(E_patch, metadata)
W_patch_to_bind <- cbind(W_patch, metadata)
N_patch_to_bind <- cbind(N_patch, metadata)
S_patch_to_bind <- cbind(S_patch, metadata)

# Add patches to the list of bioregions
Paleomaps_with_bioregions_sf_T0 <- rbind(Paleomaps_with_bioregions_sf_T0, E_patch_to_bind, W_patch_to_bind)
Paleomaps_with_bioregions_sf_T0 <- rbind(Paleomaps_with_bioregions_sf_T0, N_patch_to_bind, S_patch_to_bind)

## Create geoscale

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



## GGplot
pdf(file = "./outputs/Network_events_counts/all_dispersal_events_overall_updated_ggplot.pdf",
    width = 9, height = 6)

### 4.2/ Create network map plot ####
all_dispersal_events_overall_ggplot <- ggplot(data = nodes_metadata_updated_coords) +
  
  # Plot bioregion maps
  geom_sf(data = Paleomaps_with_bioregions_sf_T0,
          mapping = aes(fill = Bioregion),
          colour = "black",
          alpha = 1.0) +
  
  # Adjust fill color scheme for bioregions
  scale_fill_manual("Bioregions", breaks = bioregion_names_with_Antarctica, labels = bioregion_names_with_Antarctica, values = colors_list_for_areas_with_Antarctica) +
  
  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of Wrodl map
           expand = FALSE) +
  
  # Plot nodes
  geom_point(mapping = aes(x = Longitude_Mollweide, y = Latitude_Mollweide, size = mean_counts),
             alpha = 0.3, show.legend = F) +
  
  # Adjust size for points
  scale_size_continuous("Species richness", range = c(5.94, 30)) +
  
  # Plot shades of vertices/edges
  geom_curve(data = all_dispersal_events_overall_df_updated_coords[all_dispersal_events_overall_df_updated_coords$mean_counts >= min_counts_threshold, ],
             aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                 xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                 linewidth = mean_counts*1.5),
             color = "black",
             arrow = arrow(length = unit(0.03, "npc"),
                           type = "closed"), # Describes arrow head (open or closed)
             angle = 90, # Anything other than 90 or 0 can look unusual
             alpha = 0.3) +
  
  # Plot vertices/edges
  geom_curve(data = all_dispersal_events_overall_df_updated_coords[all_dispersal_events_overall_df_updated_coords$mean_counts >= min_counts_threshold, ],
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
                             breaks = c(5, 25, 75),
                             range = c(1, 6)) +
  
  # Add node labels
  ggnewscale::new_scale(new_aes = "size") +
  geom_text(mapping = aes(label = round(mean_counts, 0),
                          x = Longitude_Mollweide, y = Latitude_Mollweide,
                          size = mean_counts),
            # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
            color = "black",
            fontface = 2, show.legend = F) +
  
  # Adjust size for labels
  scale_size_continuous("Species richness", range = c(3, 10)) +
  
  # Add title
  ggtitle(label = paste0("Dispersal events between bioregions\n",
                         "Overall time (",root_age,"-0 Mya)")) +
  
  # Adjust legend aesthetics
  guides(color = "none",
         size = "none",
         fill = guide_legend(order = 1),
         linewidth = guide_legend(order = 2,
                                  override.aes = list(arrow = 0.2))
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

# print(all_dispersal_events_overall_ggplot)

### 4.3/ Create geoscale plot ####

## Add time marker
current_time <- 0

time_marker_polygon <- data.frame(
  x = c(current_time-1, current_time, current_time+1, current_time-1),
  y = c(100, 0, 100, 100))

geol_scale_T0_plot <- geol_scale_plot +
  
  geom_polygon(data = time_marker_polygon,
               mapping= aes(x = x, y = y),
               fill = "red",
               col = "black", lwd = 0.5)

# print(geol_scale_T0_plot)

### 4.4/ Arrange map and scale in a faceted plot ####

gridExtra::grid.arrange(
  grobs = list(all_dispersal_events_overall_ggplot,
               NULL,
               geol_scale_T0_plot),
  widths = 1,  # Width of columns
  heights = c(10, 0.1, 2.8),
  nrow = 3,
  ncol = 1)


dev.off()



##### 5/ Map networks over moving bioregion maps #####

### 5.1/ Load data ####

## Load node metadata per sliding windows  with updated coordinates
nodes_metadata_per_sliding_windows_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/nodes_metadata_per_sliding_windows_",window_width,"My_updated_coords.rds"))

# Load edge metadata of mean counts per sliding windows with updated coords
all_dispersal_events_all_sliding_windows_df_list_updated_coords <- readRDS(file = paste0("./outputs/Network_events_counts/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My_updated_coords.rds"))

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

N_patch <- data.frame(Latitude_dec = c(90, 90, 89.9, 89.9, 90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Eastern Palearctic",
         Subregion = "Eastern Palearctic",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(N_patch) <- st_crs(4326)

S_patch <- data.frame(Latitude_dec = c(-90, -90, -89.9, -89.9, -90),
                      Longitude_dec = c(0, 0.01, 0.01, 0, 0)) %>% 
  st_as_sf(coords = c("Longitude_dec", "Latitude_dec")) %>%
  mutate(Bioregion = "Antarctica",
         Subregion = "Antarctica",
         seq_ID = row_number()) %>%
  reconstruct_initial_polygon_with_sfheaders() %>%
  select(-Subregion, -seq_ID)
st_crs(S_patch) <- st_crs(4326)


### 5.2/ Set color scheme ####

# Set color scheme for areas/bioregions (Use the BSM color scheme)
# colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[c("A", "U", "I", "R", "N", "E", "W")]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas_with_Antarctica <- c("grey90", colors_list_for_areas)
bioregion_names_with_Antarctica <- c("Antarctica", bioregion_names)
names(colors_list_for_areas_with_Antarctica) <- bioregion_names_with_Antarctica


### 5.3/ Set size of nodes/edges ####


# Set minimum nb of dispersal events to display
# min_counts_threshold <- 1
# min_counts_threshold_list <- rep(min_counts_threshold, length(all_dispersal_events_all_sliding_windows_df_list))

min_counts_threshold_early <- 1 # From 110 to 20My
min_counts_threshold_middle <- 2 # From 19 to 8My
min_counts_threshold_late <- 3  # From 7 to 2My
time_threshold_1 <- 20 # Transition at 20My
time_threshold_2 <- 8 # Transition at 8My

min_counts_threshold_list <- rep(min_counts_threshold_early, length(nodes_metadata_time_series))
names(min_counts_threshold_list) <- nodes_metadata_time_series

min_counts_threshold_list[as.numeric(names(min_counts_threshold_list)) < time_threshold_1] <- min_counts_threshold_middle
min_counts_threshold_list[as.numeric(names(min_counts_threshold_list)) < time_threshold_2] <- min_counts_threshold_late

# Extract min/max number of events
max_node_counts <- max(nodes_metadata$mean_counts)
max_edge_counts <- max(all_dispersal_events_overall_df$mean_counts)

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


### 5.5/ Plot ggplot for each sliding window ####

for (i in seq_along(nodes_metadata_time_series))
{
  # i <- 1
  # i <- 31
  # i <- 32
  # i <- 33
  # i <- 94
  
  # Extract age
  age_i <- nodes_metadata_time_series[i]
  j <- which(Paleomaps_time_series == age_i)
  
  ## 5.5.1/ Extract data for time i and adjust display scales ####
  
  # Extract the minimum number of events to display an edge
  # min_counts_threshold_i <- min_counts_threshold
  min_counts_threshold_i <- min_counts_threshold_list[i]
  
  # Extract node metadata
  nodes_metadata_sliding_windows_i <- nodes_metadata_per_sliding_windows_updated_coords[[i]]
  # Extract max node size to adjust range of points
  max_node_size <- max(nodes_metadata_sliding_windows_i$mean_counts)
  max_range_size <- rescale_range_size(x = max_node_size)
  # Extract min node size to adjust range of points
  min_node_size <- min(nodes_metadata_sliding_windows_i$mean_counts)
  min_range_size <- rescale_range_size(x = min_node_size)
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract mean time of the sliding window
  mean_time <- mean(nodes_metadata_sliding_windows_i$time)
  round_time <- floor(mean_time)
  # start_time <- mean_time + window_width/2
  # end_time <- mean_time - window_width/2
  
  # Extract edge metadata
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_all_sliding_windows_df_list_updated_coords[[i]]
  # Extract max node size to adjust range of linewidth
  max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_counts)
  max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  # Extract min node size to adjust range of linewidth
  min_edge_linewidth <- min(all_dispersal_events_sliding_windows_i_df$mean_counts[all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i])
  min_range_linewidth <- rescale_range_linewidth(x = min_edge_linewidth)
  # Extract manual breaks for edge linewidth per stratum
  # breaks_linewidth_i <- breaks_linewidth_list[[i]]
  breaks_linewidth_i <- c(max_edge_linewidth, max_edge_linewidth, max_edge_linewidth) # Provide 3 valid entries that will be replaced in the legend
  
  # If min counts of edges = max counts of edges, need a fix by adding 0.001 to an edge count
  if (min_edge_linewidth == max_edge_linewidth)
  {
    min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_counts == min_edge_linewidth)
    edge_fix_index <- sample(x = min_egde_indices, size = 1)
    all_dispersal_events_sliding_windows_i_df$mean_counts[edge_fix_index] <- all_dispersal_events_sliding_windows_i_df$mean_counts[edge_fix_index] + 0.001
    max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_counts)
    max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  }
  
  # If number of edges to display is 1, create duplicates of that edge, and add 0.001 to the edge count
  min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i)
  if (length(min_egde_indices) == 1)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[min_egde_indices, ]
    fix_df$mean_counts <- fix_df$mean_counts + 0.001
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(fix_df$mean_counts, fix_df$mean_counts, fix_df$mean_counts)
  }
  # If number of egde to display is null, created two loop edges with different edge values and set linewidth range to zero
  if (length(min_egde_indices) == 0)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[c(1,9), ]
    fix_df$mean_counts <- c(1, 1.001)
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(1, 1, 1)
  }
  
  # Extract Paleomaps for age j
  Paleomaps_with_bioregions_sf_i <- Paleomaps_with_bioregions_sf[[j]]
  
  # Copy metadata from 1st sf line (except Bioregion)
  metadata <- st_drop_geometry(Paleomaps_with_bioregions_sf_i[1, ]) %>% 
    select(-Bioregion)
  
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
  nodes_metadata_sliding_windows_i$bioregions <- factor(nodes_metadata_sliding_windows_i$bioregions, levels = bioregion_names, labels = bioregion_names)
  all_dispersal_events_sliding_windows_i_df$source_labels <- factor(all_dispersal_events_sliding_windows_i_df$source_labels, levels = bioregion_names, labels = bioregion_names)
  
  # Adjust bioregions factors in Paleomaps
  Paleomaps_with_bioregions_sf_i$Bioregion <- factor(Paleomaps_with_bioregions_sf_i$Bioregion, levels = bioregion_names_with_Antarctica, labels = bioregion_names_with_Antarctica)
  
  ## 5.5.3/ Plot map and network ####
  
  # Plot PDF
  pdf(file = paste0("./outputs/Network_events_counts/plots_for_sliding_windows_moving_continents/all_dispersal_events_mean_counts_sliding_windows_",i,"_",window_width,"My_moving_continents_ggplot.pdf"),
      width = 9, height = 6)
  
  all_dispersal_events_sliding_windows_i_ggplot <- ggplot(data = nodes_metadata_sliding_windows_i) +
    
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
    geom_point(data = nodes_metadata_sliding_windows_i,
               mapping = aes(x = Longitude_Mollweide, y = Latitude_Mollweide, size = mean_counts),
               alpha = 0.3, show.legend = F) +
    
    # Adjust size legend
    scale_size_continuous("Species richness", 
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +
    
    # Plot shades of vertices/edges
    geom_curve(data = all_dispersal_events_sliding_windows_i_df[all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i, ],
               aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                   xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                   linewidth = mean_counts*1.5),
               color = "black",
               arrow = arrow(length = unit(0.03, "npc"), 
                             type = "closed"), # Describes arrow head (open or closed)
               angle = 90, # Anything other than 90 or 0 can look unusual
               alpha = 0.3) + 
    
    # Plot vertices/edges for dispersal
    geom_curve(data = all_dispersal_events_sliding_windows_i_df[all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i, ],
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
    geom_text(data = nodes_metadata_sliding_windows_i,
              mapping = aes(label = round(mean_counts, 0),
                            x = Longitude_Mollweide, y = Latitude_Mollweide,
                            size = mean_counts),
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +
    
    # Adjust size for labels
    scale_size_continuous("Species richness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
    
    # Add title
    ggtitle(label =  paste0("Dispersal events between bioregions\nTime = ",round_time," Mya")) +
    
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
  
  # print(all_dispersal_events_sliding_windows_i_ggplot)
  
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
    grobs = list(all_dispersal_events_sliding_windows_i_ggplot,
                 NULL,
                 geol_scale_i_plot),
    widths = 1,  # Width of columns
    heights = c(10, 0.1, 2.8),
    nrow = 3,
    ncol = 1)
  
  dev.off()
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Map created for ",round_time," My - Window n°", i, "/", length(nodes_metadata_time_series),"\n"))
  }
  
}


##### 6/ Create aggregated PDF and GIF #####

### 6.1/ Aggregate all ggplot plots in a single pdf ####

## 6.1.1/ With overlap (all windows)

all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/plots_for_sliding_windows_moving_continents/", pattern = "all_dispersal_events_mean_counts_sliding_windows_", full.names = T)
all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = paste0(window_width,"My"))]
all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = "ggplot.pdf")]
nb_ggplots <- length(all_ggplots_path)

# Reorder paths in numerical order
indices <- as.character(1:nb_ggplots)
all_ggplots_path_prefix <- str_remove(string = all_ggplots_path, pattern = "_mean_counts_sliding_windows_.*")
all_ggplots_path_suffix <- paste0("_",window_width,"My_moving_continents_ggplot.pdf")
all_ggplots_path_reordered <- paste0(all_ggplots_path_prefix, "_mean_counts_sliding_windows_", indices, all_ggplots_path_suffix)

qpdf::pdf_combine(input = all_ggplots_path_reordered, output = paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot.pdf"))

## 6.1.2/ Without overlap

window_coverage <- window_width / window_steps
indices_no_overlap <- seq(from = 1, to = nb_ggplots, by = (window_coverage + 1))

all_ggplots_path_reordered_no_overlap <- paste0(all_ggplots_path_prefix[1], "_mean_counts_sliding_windows_", indices_no_overlap, all_ggplots_path_suffix[1])

qpdf::pdf_combine(input = all_ggplots_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot_no_overlap.pdf"))


### 6.2/ Aggregate all ggplot in forward timeline + overall ####

## 6.2.1/ With overlap (all windows)

overall_ggplot_path <- "./outputs/Network_events_counts/all_dispersal_events_overall_updated_ggplot.pdf"
qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot_forward.pdf"))

## 6.2.2/ Without overlap

qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered_no_overlap), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot_forward_no_overlap.pdf"))


### 6.3/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

window_width <- 5
fps <- 5
# fps <- 10

pdf_pointer_mean_counts <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot_forward.pdf"),
                                                 pages = NULL, density = 150)
magick::image_info(pdf_pointer_mean_counts)

image_resize_and_write_gif(image = pdf_pointer_mean_counts,
                           path =  paste0("./outputs/Network_events_counts/all_dispersal_events_mean_counts_all_sliding_windows_",window_width,"My_moving_continents_ggplot_forward_",fps,"fps.gif"),
                           delay = 1/fps, # Time between frames in seconds
                           width = 1350, height = 900,
                           loop = FALSE,
                           progress = TRUE)


## Add arrow to legend (Does not work with guides. Need to add on PPT)
  # Try using manual polygons?


## Create bioregion shapefile that can overlay the PALEOMAPS?
  # Use the PALEOMAPS model instead of MULLER2022


## All can be segregated per clades!
# Across clades to compare clades dynamics
# Within clades, to zoom on a particular subclade of interest


