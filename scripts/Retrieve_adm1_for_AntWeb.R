##### Script XX: Retrieve adm1 for AntWeb #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

##### 1/ Load stuff ####

# Intall rgeos from Archive
# install.packages("./packages/rgeos_0.6-4.tar.gz", repos = NULL, type = "source")

# Install lwgeom from Source
# install.packages("lwgeom", type = "source")

### 1.1/ Load packages ####

library(readxl)
library(openxlsx)  # Use Rccp. No need of Java
library(sf)
library(tidyverse)

### 1.2/ Load data ####

## Load occurrences
Missing_Adm1s <- read_excel("input_data/Biogeographic_data/Missing Adm1s.xlsx")
Missing_Adm1s$Occurrence_ID <- 1:nrow(Missing_Adm1s)

## Convert to sf
Missing_Adm1s$Latitude_dec <- Missing_Adm1s$verbatimLatitude
Missing_Adm1s$Longitude_dec <- Missing_Adm1s$verbatimLongitude
Missing_Adm1s <- st_as_sf(Missing_Adm1s, coords = c("Longitude_dec", "Latitude_dec"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Missing_Adm1s) <- sp::CRS('+init=EPSG:4326')
# Remove copies of Long/Lat of still present
# Missing_Adm1s <- Missing_Adm1s %>% select(-c("Longitude_dec", "Latitude_dec"))


## Load shp file for Countries
geoBoundaries_Countries_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM0/", layer = "geoBoundariesCGAZ_ADM0")
# plot(geoBoundaries_Countries_sf["shapeName"])
table(geoBoundaries_Countries_sf$shapeGroup)

# Create replicate to look for intersecting polygons
geoBoundaries_Countries_sf_Intersect <- geoBoundaries_Countries_sf %>% 
  rename(GB_Country_Intersect = shapeName) %>% 
  select(GB_Country_Intersect)

# Create replicate to look for Nearest polygons
geoBoundaries_Countries_sf_Nearest <- geoBoundaries_Countries_sf %>% 
  rename(GB_Country_Nearest = shapeName) %>% 
  select(GB_Country_Nearest)


## Load shp file for adm1
geoBoundaries_adm1_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM1/", layer = "geoBoundariesCGAZ_ADM1")
# plot(geoBoundaries_adm1_sf["shapeName"])
table(geoBoundaries_adm1_sf$shapeGroup)

# Create replicate to look for intersecting polygons
geoBoundaries_adm1_sf_Intersect <- geoBoundaries_adm1_sf %>% 
  rename(GB_adm1_Intersect = shapeName) %>% 
  select(GB_adm1_Intersect)

# Create replicate to look for Nearest polygons
geoBoundaries_adm1_sf_Nearest <- geoBoundaries_adm1_sf %>% 
  rename(GB_adm1_Nearest = shapeName) %>% 
  select(GB_adm1_Nearest)


##### 2/ Find intersecting and Nearest polygons #####

# Step needed to limit the number of run of st_intersection which would be too big to run as a loop per entries

out <- data.frame() # Store output
batch_size <- 100   # Number of entries per batch
nb_batches <- (nrow(Missing_Adm1s) %/% batch_size) + is.numeric((nrow(Missing_Adm1s) %% batch_size) > 0)
seq_i <- seq(from = 1, to = batch_size * nb_batches + 1, by = batch_size)

sf_use_s2(FALSE)
for (i in seq_i)
{
  # i <- 3701
  
  # Get indices of batch entries
  batch_i <- seq(from = i, by = 1, length.out = batch_size)
  
  # Adjust indices for the last batch
  if (max(batch_i) > nrow(Missing_Adm1s)) { batch_i <- min(batch_i):nrow(Missing_Adm1s) }
  
  # Find intersecting polygons for Countries
  out_batch_Countries_Intersect <- suppressMessages(suppressWarnings(sf::st_intersection(Missing_Adm1s[batch_i, ], geoBoundaries_Countries_sf_Intersect[, "GB_Country_Intersect"])))
  out_batch_Countries_Intersect <- out_batch_Countries_Intersect %>% 
    select(Occurrence_ID, GB_Country_Intersect)
  
  # Find intersecting polygons for adm1
  out_batch_adm1_Intersect <- suppressMessages(suppressWarnings(sf::st_intersection(Missing_Adm1s[batch_i, ], geoBoundaries_adm1_sf_Intersect[, "GB_adm1_Intersect"])))
  out_batch_adm1_Intersect <- out_batch_adm1_Intersect %>% 
    select(Occurrence_ID, GB_adm1_Intersect)
    
  # Find nearest polygons for Countries
  out_batch_Countries_Nearest_ID <- suppressMessages(suppressWarnings(sf::st_nearest_feature(Missing_Adm1s[batch_i, ], geoBoundaries_Countries_sf_Nearest[, "GB_Country_Nearest"])))
  out_batch_Countries_Nearest <- geoBoundaries_Countries_sf_Nearest[out_batch_Countries_Nearest_ID, ]
  out_batch_Countries_Nearest$Occurrence_ID <- batch_i
  
  # Find nearest polygons for adm1
  out_batch_adm1_Nearest_ID <- suppressMessages(suppressWarnings(sf::st_nearest_feature(Missing_Adm1s[batch_i, ], geoBoundaries_adm1_sf_Nearest[, "GB_adm1_Nearest"])))
  out_batch_adm1_Nearest <- geoBoundaries_adm1_sf_Nearest[out_batch_adm1_Nearest_ID, ]
  out_batch_adm1_Nearest$Occurrence_ID <- batch_i
  
  # Keep only the first match per occurrence
  out_batch_Countries_Intersect <- out_batch_Countries_Intersect %>% 
    arrange(Occurrence_ID) %>% 
    group_by(Occurrence_ID) %>% 
    mutate(Duplicates_counter = row_number(Occurrence_ID)) %>% 
    filter(Duplicates_counter == 1) %>% 
    select(-Duplicates_counter) %>%
    ungroup()
  out_batch_adm1_Intersect <- out_batch_adm1_Intersect %>% 
    arrange(Occurrence_ID) %>% 
    group_by(Occurrence_ID) %>% 
    mutate(Duplicates_counter = row_number(Occurrence_ID)) %>% 
    filter(Duplicates_counter == 1) %>% 
    select(-Duplicates_counter) %>%
    ungroup()

  # Join batch output for multiple levels, including entry for occurrences with no match (keeping track with ID)
  out_batch <- suppressMessages(left_join(x = st_drop_geometry(Missing_Adm1s[batch_i, ]), y = st_drop_geometry(out_batch_Countries_Intersect)))
  out_batch <- suppressMessages(left_join(x = out_batch, y = st_drop_geometry(out_batch_Countries_Nearest)))
  out_batch <- suppressMessages(left_join(x = out_batch, y = st_drop_geometry(out_batch_adm1_Intersect)))
  out_batch <- suppressMessages(left_join(x = out_batch, y = st_drop_geometry(out_batch_adm1_Nearest)))

  # Merge batch output with global output
  out <- rbind(out, out_batch)
  
  # Print progress
  cat(paste0(Sys.time(), " - Occurrence n°",max(batch_i), " / ", nrow(Missing_Adm1s),"\n"))
  
}
sf_use_s2(TRUE)

# Check that the final numbers of row are equal
nrow(out)
nrow(Missing_Adm1s)

View(out)

# If the output file seems valid, replace initial file
Missing_Adm1s <- out

## Convert back to sf
Missing_Adm1s$Latitude_dec <- Missing_Adm1s$verbatimLatitude
Missing_Adm1s$Longitude_dec <- Missing_Adm1s$verbatimLongitude
Missing_Adm1s <- st_as_sf(Missing_Adm1s, coords = c("Longitude_dec", "Latitude_dec"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Missing_Adm1s) <- sp::CRS('+init=EPSG:4326')
# Remove copies of Long/Lat of still present
# Missing_Adm1s <- Missing_Adm1s %>% select(-c("Longitude_dec", "Latitude_dec"))

## Save Biogeographic database with all curated coordinates
saveRDS(object = Missing_Adm1s, file = "./input_data/Biogeographic_data/Missing_Adm1s.rds")


##### 3/ Flag for inconsistencies #####

### 3.1/ Flag entries for which TW data does not match GB data ####

# May deserve inspection to check if the change is valid

# For Countries
Missing_Adm1s$Mismatch_TW_Country <- replace_na(Missing_Adm1s$`TW country` != Missing_Adm1s$GB_Country_Intersect, replace = F)
table(Missing_Adm1s$Mismatch_TW_Country)
View(Missing_Adm1s[Missing_Adm1s$Mismatch_TW_Country, ])

# For adm1
Missing_Adm1s$Mismatch_TW_adm1 <- replace_na(Missing_Adm1s$`TW stateProvince` != Missing_Adm1s$GB_adm1_Intersect, replace = F)
table(Missing_Adm1s$Mismatch_TW_adm1)
View(Missing_Adm1s[Missing_Adm1s$Mismatch_TW_adm1, ])

### 3.2/ Flag entries for which the Intersecting shape is different from the Nearest shape (no intersecting shape) ####

# Three cases:
  # 1/ The coordinates falls in a waterbody and is close to its true location => Can use nearest polygon. Should correct coordinates.
  # 2/ The coordinates falls in a waterbody and is a mistake => CANNOT use nearest polygon. Need to correct coordinates.
  # 3/ The coordinates falls on land, but this land is not included in geoBoundaries shape files => Should remain a true NA. Coordinates are fine.

# For Countries
Missing_Adm1s$Mismatch_Nearest_Country <- replace_na(Missing_Adm1s$GB_Country_Intersect != Missing_Adm1s$GB_Country_Nearest, replace = T)
table(Missing_Adm1s$Mismatch_Nearest_Country, is.na(Missing_Adm1s$GB_Country_Intersect))
View(Missing_Adm1s[Missing_Adm1s$Mismatch_Nearest_Country, ])

# For adm1
Missing_Adm1s$Mismatch_Nearest_adm1 <- replace_na(Missing_Adm1s$GB_adm1_Intersect != Missing_Adm1s$GB_adm1_Nearest, replace = T)
table(Missing_Adm1s$Mismatch_Nearest_adm1, is.na(Missing_Adm1s$GB_adm1_Intersect))
View(Missing_Adm1s[Missing_Adm1s$Mismatch_Nearest_adm1, ])


##### 4/ Export results #####

saveRDS(Missing_Adm1s, file = "./input_data/Biogeographic_data/Missing_Adm1s.rds")
write.xlsx(x = st_drop_geometry(Missing_Adm1s), file = "./input_data/Biogeographic_data/Missing_Adm1s.xlsx")
