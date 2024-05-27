##### Script 02: Curate biogeographic information #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate the specimen-level database for biogeographic information
  # Detect errors in coordinates and manually curated them
  # Use maps to visually detect outliers and correct them
  # Draw random coordinates when missing
  # Clean duplicates
# Aggregated information at taxa-level and update the taxa-level summary df
  # Nb of occurrences
  # Presence/Absence in each Bioregion per taxa
     # Create multiple Bioregion schemes

###

### Inputs

# Specimen-level database
# Taxa-level summary database
# NaturalEarth administrative shapefile to provide country shp for CoordinatesCleaner
# GeoBoundaries shapefile to provide new polygon membership for both GABI and AntWeb
# Bentity2 shapefile from GABI

###

### Sources

# AntWeb. Version 8.101. California Academy of Science, online at https://www.antweb.org. Accessed 14 November 2023.
# GABI: Data release 1.0 of 18 January 2020.
   # Guénard, B., Weiser, M., Gomez, K., Narula, N., Economo, E.P. (2017) The Global Ant Biodiversity Informatics (GABI) database: synthesizing data on the geographic distributions of ant species. Myrmecological News 24: 83-89.

###

### Outputs

# Taxa-level database for taxonomic, ecologic and biogeographic information
# Curated Specimen-level database
# Maps of dubious records per taxa
# Genus-level summary table for biogeographic occcurrences

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

# Intall rgeos from Archive
# install.packages("./packages/rgeos_0.6-4.tar.gz", repos = NULL, type = "source")

# Install lwgeom from Source
# install.packages("lwgeom", type = "source")

### 1.1/ Load packages ####
library(tidyverse)
library(readxl)
library(openxlsx)  # Use Rccp. No need of Java
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(phytools)
library(ape)
library(rnaturalearth)
library(rnaturalearthdata)
library(CoordinateCleaner)
library(mapview)
library(sf)
library(sp)
library(rgeos)
library(units)
library(RColorBrewer)
library(htmlwidgets)
library(htmltools)
library(lwgeom)


### 1.2/ Load databases to clean ####

## Load the curated databases
GABI_database_Ponerinae <- readRDS(file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database_Ponerinae.rds")
AntWeb_database_curated <- readRDS(file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")
Biogeographic_database_Ponerinae <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")

## Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

### 1.3/ Save initial geographic information (to compare if updated) ####

Biogeographic_database_Ponerinae <- Biogeographic_database_Ponerinae %>% 
  mutate(Latitude_dec_initial = Latitude_dec) %>%
  mutate(Longitude_dec_initial = Longitude_dec) %>%
  rename(Country_initial = Country) %>%
  mutate(adm1_initial = adm1) %>%
  mutate(adm2_initial = adm2) %>%
  mutate(Locality_initial = Locality) %>%
  mutate(bentity2_name_initial = bentity2_name) %>%
  mutate(Elevation_initial = Elevation)
  
##### 2/ Check dubious records with CoordinateCleaner #####

?CoordinateCleaner::cc_aohi()    # Flag for coordinates falling into Artificial Hotspots listed by Park et al., 2022
?CoordinateCleaner::cc_cap()     # Flag capitals
?CoordinateCleaner::cc_cen()     # Flag country centroids
?CoordinateCleaner::cc_coun()    # Flag for mismatch between country ID and GPS coordinates
?CoordinateCleaner::cc_dupl()    # Flag duplicates based on name and coordinates
?CoordinateCleaner::cc_equ()     # Flag for equal long/lat
?CoordinateCleaner::cc_inst()    # Flag for coordinates in the vicinity of biodiversity institutions
?CoordinateCleaner::cc_outl()    # Flag for outliers based on mean or minimal distance to other points per taxa
?CoordinateCleaner::cc_sea()     # Flag for coordinates outside of landmass
?CoordinateCleaner::cc_urb()     # Flag for coordinates within urban areas
?CoordinateCleaner::cc_val()     # Flag for coordinates within erroneous and NA coordinates
?CoordinateCleaner::cc_zero()    # Flag for coordinates with null lat or long

?CoordinateCleaner::clean_coordinates # Run all tests


### 2.1/ Assign country ISO codes used to detect mismatch with CoordinateCleaner ####

?CoordinateCleaner::countryref

## Extract country metadata

# # From CoordinateCleaner
# Countries_metadata <- CoordinateCleaner::countryref

# From Natural Earth (since the polygon layer used for matching is from NE, better use directly metadata from NE)
# Countries_NE_sf <- rnaturalearth::ne_countries('medium', returnclass = "sf")
Countries_NE_sf <- rnaturalearth::ne_load(file_name = "ne_10m_admin_0_countries", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/Natural_Earth_quick_start/10m_cultural/")
plot(Countries_NE_sf["ISO_A3"])
Countries_NE_metadata <- st_drop_geometry(Countries_NE_sf[, c("SUBUNIT", "ISO_A3")])

table(Countries_NE_metadata$ISO_A3)

# Remove NA
# Countries_NE_metadata <- Countries_NE_metadata[complete.cases(Countries_NE_metadata), ]

#Countries_NE_metadata[!complete.cases(Countries_NE_metadata), ]
Countries_NE_metadata[Countries_NE_metadata$ISO_A3 == "-99", ]
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("France"), "ISO_A3"]) # Includes departments such as Guadeloupe, Martinique, French Guyana, Reunion, and Mayotte
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Norway"), "ISO_A3"]) # Includes Svaldbard and other islands
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Kosovo", "Republic of Serbia"), "ISO_A3"]) 
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Akrotiri Sovereign Base Area", "Dhekelia Sovereign Base Area", "Cyprus", "Cyprus No Mans Area", "Northern Cyprus"), "ISO_A3"])  # British military bases and Turkish area in Cyprus
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Ashmore and Cartier Islands", "Australia", "East Timor"), "ISO_A3"]) # North of Australia
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Bajo Nuevo Bank (Petrel Is.)", "Serranilla Bank", "Colombia", "Jamaica"), "ISO_A3"]) # Caribbean island/reef. Part of Columbia but disputed.
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Baykonur Cosmodrome", "Kazakhstan"), "ISO_A3"]) # In Kazakhstan
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Bir Tawil", "Egypt", "Sudan"), "ISO_A3"]) # Unclaimed territory at the border of Egypt and Sudan
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Brazilian Island", "Uruguay", "Brazil", "Argentina"), "ISO_A3"]) # Disputed Island in Rio Uruguay at the border of Uruguay, Brazil, and Argentina
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Clipperton Island", "Mexico"), "ISO_A3"]) 
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Australia"), "ISO_A3"]) # SE of Australia
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Indian Ocean Territories", "Australia"), "ISO_A3"]) # NW of Australia. Cocos Islands + Christmas Islands
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Scarborough Reef", "Spratly Islands", "Philippines"), "ISO_A3"]) # East of the Philippines. Disputed with China and Taiwan.
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("India", "Siachen Glacier", "Pakistan"), "ISO_A3"]) # Disputed glacier at the border of Pakistan and India
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Somaliland", "Somalia"), "ISO_A3"]) # Norht of Somalia
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Southern Patagonian Ice Field", "Chile", "Argentina"), "ISO_A3"]) # Disputed glacier at the border of Chile and Argentina
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("US Naval Base Guantanamo Bay", "Cuba"), "ISO_A3"]) # US Military base in Cuba

# # Replace NA in "medium" scale
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "Kosovo"] <- "XXK"
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "Ashmore and Cartier Is."] <- "ATC" # North of Australia
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "N. Cyprus"] <- "CYN"
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "Indian Ocean Ter."] <- "CCI" # Christmas Island and Cocos Islands
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "Siachen Glacier"] <- "SCG" # India-Pakistan border
# Countries_NE_metadata$iso_a3[Countries_NE_metadata$name == "Somaliland"] <- "SOL"

# Replace -99 in "small" scale
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Akrotiri Sovereign Base Area"] <- "ASB" # British military bases in Cyprus
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Ashmore and Cartier Islands"] <- "ATC" # North of Australia
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Bajo Nuevo Bank (Petrel Is.)"] <- "BNB" # Caribbean island/reef. Part of Columbia but disputed.
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Baykonur Cosmodrome"] <- "BKC" # In Kazakhstan
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Bir Tawil"] <- "BTW" # Unclaimed territory at the border of Egypt and Sudan
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Brazilian Island"] <- "BIS" # Disputed Island in Rio Uruguay at the border of Uruguay, Brazil, and Argentina
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Clipperton Island"] <- "CPT" # Disputed Island in Rio Uruguay at the border of Uruguay, Brazil, and Argentina
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Coral Sea Islands"] <- "CSI" # SE of Australia
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Cyprus No Mans Area"] <- "CNM"
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Dhekelia Sovereign Base Area"] <- "DSB" # British military bases in Cyprus
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "France"] <- "FRA" # Includes departments such as Guadeloupe, Martinique, French Guyana, Reunion, and Mayotte
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Indian Ocean Territories"] <- "CCI" # Christmas Island and Cocos Islands
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Kosovo"] <- "KOS"
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Northern Cyprus"] <- "CPN" # Turkish part of Cyprus
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Norway"] <- "NOR" # Includes Svaldbard and other islands
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Scarborough Reef"] <- "SCB" # East of the Philippines. Disputed with China and Taiwan.
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Serranilla Bank"] <- "SNB" # Caribbean island/reef. Part of Columbia but disputed.
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Siachen Glacier"] <- "SCG" # India-Pakistan border
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Somaliland"] <- "SOL" # North of Somalia
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Southern Patagonian Ice Field"] <- "SPG" # Disputed glacier at the border of Chile and Argentina
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "Spratly Islands"] <- "SPI" # East of the Philippines. Disputed with China and Taiwan.
Countries_NE_metadata$ISO_A3[Countries_NE_metadata$SUBUNIT == "US Naval Base Guantanamo Bay"] <- "GUT" # US Military base in Cuba

# Countries_NE_metadata[duplicated(Countries_NE_metadata$ISO_A3),]

## Save Countries_NE_sf with updated metadata to use as reference for clean_coordinates
Countries_NE_sf$ISO_A3 <- Countries_NE_metadata$ISO_A3
saveRDS(Countries_NE_sf, file = "./input_data/NaturalEarth_maps/Countries_NE_sf.rds")

# Extract only ISO3 entities
# Countries_metadata_IOS3 <- (Countries_metadata[!duplicated(Countries_metadata$iso3),])
View(Countries_NE_metadata)

# Create new field for ISO3 names
Biogeographic_database_Ponerinae$Country_ISO3_name <- Biogeographic_database_Ponerinae$Country_initial

# Check for names that do not match
# unique(Biogeographic_database_Ponerinae$Country_ISO3_name[!(Biogeographic_database_Ponerinae$Country_ISO3_name %in% Countries_NE_metadata$name)])
unique(Biogeographic_database_Ponerinae$Country_ISO3_name[!(Biogeographic_database_Ponerinae$Country_ISO3_name %in% Countries_NE_metadata$SUBUNIT)])

plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("French Southern and Antarctic Lands", "Madagascar"), "ISO_A3"]) # Includes Kerguelen, Europa island and other stuff
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("New Zealand"), "ISO_A3"]) # Includes Tokelau
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Netherlands"), "ISO_A3"]) # Includes Bonaire, Oranjestad and Saba
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Spain"), "ISO_A3"]) # Includes the Canarias
plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("United Kingdom"), "ISO_A3"]) # Includes the Canarias

# Repair mismatches
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "USA"] <- "United States"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Czech Republic"] <- "Czechia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Zaire"] <- "Democratic Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cote d'Ivoire"] <- "Ivory Coast"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "French Guiana"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Serbia"] <- "Republic of Serbia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "England"] <- "United Kingdom"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo (Brazzaville)"] <- "Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent"] <- "Saint Vincent and the Grenadines"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Christmas Island"] <- "Indian Ocean Territories"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Trinidad"] <- "Trinidad and Tobago"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bahamas"] <- "The Bahamas"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guinea-Bissau [Portuguese Guinea]"] <- "Guinea-Bissau"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea"] <- "Papua New Guinea" # May also be Indonesia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Martinique"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cape Verde"] <- "Cabo Verde"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mariana Islands"] <- "Northern Mariana Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mayotte"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Borneo"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia MULT"] <- "Malaysia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia, Federated States of"] <- "Federated States of Micronesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Great Britain"] <- "United Kingdom"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "?Unknown"] <- "Lesser Antilles"  # Could be many different Islands
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guadeloupe Islands"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "U.S. Virgin Islands"] <- "United States Virgin Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guadeloupe"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea ERR"] <- "Papua New Guinea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "British West Indies"] <- "British West Indies" # Cannot define suitable ISO-3 entities
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea [Dutch New Guinea]"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Channel Islands"] <- "Guernsey" # Could also be Jersey
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea [German New Guinea]"] <- "Papua New Guinea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Tokelau"] <- "New Zealand"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cocos (Keeling) Islands"] <- "Indian Ocean Territories"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hispaniola"] <- "Dominican Republic"  # Could also be in Haiti
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia MULT"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts & Nevis"] <- "Saint Kitts and Nevis"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "US Virgin Islands"] <- "United States Virgin Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia ERR"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea ERR"] <- "Papua New Guinea" # May also be Indonesia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Brunei Darussalam"] <- "Brunei"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands ERR"] <- "Solomon Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Comoros Islands"] <- "Comoros"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "N. Mariana Is."] <- "Northern Mariana Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Wallis and Futuna Islands"] <- "Wallis and Futuna"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea IND"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia [Dutch New Guinea]"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Sao Tome and Principe"] <- "São Tomé and Principe"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Netherlands Antilles"] <- "Aruba" # Could also be Curaçao or Bonaire (Netherlands)
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "UK"] <- "United Kingdom"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea IND/PAPUA"] <- "Papua New Guinea" # May also be Indonesia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Seychelles Islands"] <- "Seychelles"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia"] <- "Federated States of Micronesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands MULT"] <- "Solomon Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bonaire, Sint Eustatius and Saba"] <- "Netherlands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hong Kong"] <- "Hong Kong S.A.R."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo"] <- "Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo [French Congo]"] <- "Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia Peninsular"] <- "Malaysia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macedonia"] <- "North Macedonia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia ERR"] <- "Malaysia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Swaziland"] <- "eSwatini"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macau"] <- "Macao S.A.R"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "?Unknown (Somalia? Or Ethiopia?)"] <- "Somaliland" # Could be Ethiopia/Somalia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Democratic Republic of Congo"] <- "Democratic Republic of the Congo"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Port of Entry"] <- NA
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "null"] <- NA
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Unknown"] <- NA
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Rhodesia"] <- "Zimbabwe"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Timor-Leste"] <- "East Timor"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macaronesia"] <- "Spain"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Reunion"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Europa Island"] <- "French Southern and Antarctic Lands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Panam┬╖"] <- "Panama"

# Biogeographic_database_Ponerinae[Biogeographic_database_Ponerinae$Country_initial == "Macaronesia", ]

# Add ISO-3 code to Biogeographic database
# Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_metadata$iso3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_metadata$name)]
# Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_NE_metadata$iso_a3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_NE_metadata$name)]
Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_NE_metadata$ISO_A3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_NE_metadata$SUBUNIT)]

table(Biogeographic_database_Ponerinae$Country_ISO3_code)

View(Biogeographic_database_Ponerinae[(is.na(Biogeographic_database_Ponerinae$Country_ISO3_code)), ])

## Save specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")

### 2.2/ Remove "Port of Entry data" ####

table(Biogeographic_database_Ponerinae$Country_initial == "Port of Entry")

Biogeographic_database_Ponerinae <- Biogeographic_database_Ponerinae[!(Biogeographic_database_Ponerinae$Country_initial == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags[!(Biogeographic_database_Ponerinae_flags$Country_initial == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_NA_coords <- Biogeographic_database_Ponerinae_NA_coords[!(Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_no_NA_coords <- Biogeographic_database_Ponerinae_no_NA_coords[!(Biogeographic_database_Ponerinae_no_NA_coords$Country_initial == "Port of Entry"), ]

### 2.3/ Repair broken coordinates ####

# Detect invalid coordinates
Biogeographic_database_Ponerinae$valid_coordinates <- cc_val(x = Biogeographic_database_Ponerinae,
                                                            lon = "Longitude_dec", lat = "Latitude_dec", 
                                                            value = "flagged", verbose = T)

# Flag entries with no coordinates
all_NA_coords <- apply(X = Biogeographic_database_Ponerinae[, c("Latitude_dec", "Longitude_dec")], MARGIN = 1, FUN = function (x) {all(is.na(x))})

# 25 cases of incomplete coordinates with either Latitude or Longitude missing
View(Biogeographic_database_Ponerinae[!Biogeographic_database_Ponerinae$valid_coordinates & !all_NA_coords, ])

table(Biogeographic_database_Ponerinae$Latitude_dec[!Biogeographic_database_Ponerinae$valid_coordinates])
table(Biogeographic_database_Ponerinae$Longitude_dec[!Biogeographic_database_Ponerinae$valid_coordinates])

## Repair them using localities using known mean and sd of latitude for the same locality

# Gorongosa in Mozambique, repaired using known mean and sd of latitude for the same locality
Gorongosa_n <- sum((Biogeographic_database_Ponerinae$adm2 == "Gorongosa") & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec), na.rm = T)
Gorongosa_mean <- mean(Biogeographic_database_Ponerinae$Latitude_dec[Biogeographic_database_Ponerinae$adm2 == "Gorongosa"], na.rm = T)
Gorongosa_sd <- sd(Biogeographic_database_Ponerinae$Latitude_dec[Biogeographic_database_Ponerinae$adm2 == "Gorongosa"], na.rm = T)
Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae$adm2 == "Gorongosa"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- rnorm(n = Gorongosa_n, mean = Gorongosa_mean, sd = Gorongosa_sd)

# Gangwon Province in South Korea, repaired using known mean and sd of latitude for the same locality
Gangwon_n <- sum((Biogeographic_database_Ponerinae$adm1 == "Gangwon Province") & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec), na.rm = T)
Gangwon_mean <- mean(Biogeographic_database_Ponerinae$Latitude_dec[Biogeographic_database_Ponerinae$adm1 == "Gangwon Province"], na.rm = T)
Gangwon_sd <- sd(Biogeographic_database_Ponerinae$Latitude_dec[Biogeographic_database_Ponerinae$adm1 == "Gangwon Province"], na.rm = T)
Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae$adm1 == "Gangwon Province"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- rnorm(n = Gangwon_n, mean = Gangwon_mean, sd = Gangwon_sd)

# Texas Tech University Center repaired using unique latitude known for the same locality
table(Biogeographic_database_Ponerinae$Latitude_dec[Biogeographic_database_Ponerinae$Locality == "Texas Tech University Center, Junction"])
Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae$Locality == "Texas Tech University Center, Junction"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 30.471308
Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae$Locality == "Texas Tech University Center, Junction, Lubbock"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 30.471308
Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae$Locality == "Texas Tech University Center, Junction:Lubbock"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 30.471308

# Parque Nacional do Jaú, Castanho, Baixo Jaú using unique latitude known for the same locality
table(Biogeographic_database_Ponerinae$Latitude_dec[str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Baixo Jaú")])
Biogeographic_database_Ponerinae$Latitude_dec[(replace_na(data = str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Baixo Jaú"), replace = F)) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- -1.95

# Kruger Nat. Park in South Africa, repaired using known mean and sd of latitude for the same locality
Kruger_n <- sum((str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Kruger")) & !is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec), na.rm = T)
Kruger_values <- Biogeographic_database_Ponerinae$Longitude_dec[str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Kruger")]
Kruger_values <- Kruger_values[!is.na(Kruger_values) & (Kruger_values > 30)]
Kruger_mean <- mean(Kruger_values, na.rm = T)
Kruger_sd <- sd(Kruger_values, na.rm = T)
Biogeographic_database_Ponerinae$Longitude_dec[replace_na(str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Kruger"), replace = F) & !is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- rnorm(n = Kruger_n, mean = Kruger_mean, sd = Kruger_sd)


## If no available entries with GPS coordinates, repair them using GPS coordinates found from Google Maps.

# Buxa Tiger Reserve in India
Biogeographic_database_Ponerinae$Latitude_dec[(Biogeographic_database_Ponerinae$Locality == "Buxa Tiger Reserve") & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 26.616

# Carle Woods (Carl R. Hansen Woods) in Illinois, USA
Biogeographic_database_Ponerinae$Latitude_dec[(Biogeographic_database_Ponerinae$Locality == "Carle Woods") & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 42.053
Biogeographic_database_Ponerinae$Latitude_dec[(Biogeographic_database_Ponerinae$Locality == "Carle Woods") & replace_na(data = (Biogeographic_database_Ponerinae$Latitude_dec == 0), replace = F) & !is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 42.053

# Asuncion Botanical Garden in Paraguay
Biogeographic_database_Ponerinae$Longitude_dec[(str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Asuncion Botanical Garden")) & !is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- -57.576

# # Oronoque River, Guyana
# Biogeographic_database_Ponerinae$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Oronoque River"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- -57.422
# Biogeographic_database_Ponerinae$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Oronoque River"), replace = F) & is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 2.6

# North Transvaal in South Africa
Biogeographic_database_Ponerinae$Longitude_dec[(str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "North Transvaal")) & !is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- 29.789


View(Biogeographic_database_Ponerinae[!Biogeographic_database_Ponerinae$valid_coordinates & !all_NA_coords, ])
# No more missing coordinates

# Reflag coordinates invalid coordinates
Biogeographic_database_Ponerinae$valid_coordinates <- cc_val(x = Biogeographic_database_Ponerinae,
                                                             lon = "Longitude_dec", lat = "Latitude_dec", 
                                                             value = "flagged", verbose = T)

table(Biogeographic_database_Ponerinae$Latitude_dec[!Biogeographic_database_Ponerinae$valid_coordinates])
table(Biogeographic_database_Ponerinae$Longitude_dec[!Biogeographic_database_Ponerinae$valid_coordinates])
# Empty because only NA

## Save specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")

### 2.4/ Curate entries with zeros ####

#View(Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$.zer, ])

## Set to NA all points with (0,0) coordinates

all_zero_coords <- (Biogeographic_database_Ponerinae$Longitude_dec == 0) & (Biogeographic_database_Ponerinae$Latitude_dec == 0)

Biogeographic_database_Ponerinae$Latitude_dec[all_zero_coords] <- NA
Biogeographic_database_Ponerinae$Longitude_dec[all_zero_coords] <- NA

# Remaining entries should be credible despite Latitude == 0
View(Biogeographic_database_Ponerinae[replace_na(data = (Biogeographic_database_Ponerinae$Latitude_dec == 0), replace = F), ])

# Reflag coordinates invalid coordinates
Biogeographic_database_Ponerinae$valid_coordinates <- cc_val(x = Biogeographic_database_Ponerinae,
                                                             lon = "Longitude_dec", lat = "Latitude_dec", 
                                                             value = "flagged", verbose = T)

## Save specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")


### 2.5/ Run global analysis ####

## Load specimen database for Biogeographic data on Ponerinae
Biogeographic_database_Ponerinae <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")


## 2.5.1/ Prepare landmass polygons ####

## Load more detailed reference for landmasses polygons
# Default is scale 110/small which left aside many small islands

map_NE_lands_sf <- rnaturalearth::ne_load(file_name = "ne_10m_land", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_10m_land/")
map_NE_minor_islands_sf <- rnaturalearth::ne_load(file_name = "ne_10m_minor_islands", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_10m_minor_islands/")
map_NE_all_lands_sf <- rbind(x = map_NE_lands_sf, y = map_NE_minor_islands_sf)

plot(map_NE_lands_sf["featurecla"])
plot(map_NE_minor_islands_sf["featurecla"])
plot(map_NE_all_lands_sf["featurecla"])

## Remove small islands (< 5 km²) to save computation time!

map_NE_lands_sf_cast <- st_as_sf(st_cast(map_NE_lands_sf, "POLYGON")) # Disaggregate MUTLIPOLYGON into POLYGONS
map_NE_lands_sf_cast$area <- st_area(map_NE_lands_sf_cast)
map_NE_lands_sf_cast$area <- set_units(map_NE_lands_sf_cast$area, km^2)

plot(map_NE_lands_sf_cast["area"])
plot(map_NE_lands_sf_cast[map_NE_lands_sf_cast$area > as_units(5, "km^2"), "area"])

map_NE_lands_sf_cast <- map_NE_lands_sf_cast[map_NE_lands_sf_cast$area > as_units(5, "km^2"), "area"]

## Apply buffer ourselves because terra version in CoordinateClean::cc_sea leads to error

# Consider lat/long to be Cartesian spatial coordinates.
# Not strictly correct. Lead to distance deformation and issue with buffer crossing -180/180 line...
# Distance in degrees. 0.05 ≈ 5.55 km at the equator

sf_use_s2(FALSE) # Need to desactivate s2 features that cannot deal with antimeridian crossing and coordinates beyond lat/long limits
map_NE_lands_sf_buffer <- st_buffer(x = map_NE_lands_sf_cast, dist = 0.05)
sf_use_s2(TRUE)

# map_NE_lands_sf_buffer <- st_make_valid(st_buffer(x = map_NE_lands_sf, dist = 0.1))
plot(map_NE_lands_sf["featurecla"])
plot(map_NE_lands_sf_buffer["area"])

## Adjust boundaries
# Need to limit to a 180/90 rectangle (extent beyond when adding buffer to Antarctica for instance)
st_bbox(map_NE_lands_sf)
st_bbox(map_NE_lands_sf_buffer)

geo_limits <- st_geometry(st_sfc(st_point(c(-180,-90)), st_point(c(180,90))), crs = st_crs(map_NE_lands_sf_buffer))
st_crs(geo_limits) <- st_crs(map_NE_lands_sf_buffer)
sf_use_s2(FALSE) # Need to desactivate s2 features that cannot deal with antimeridian crossing and coordinates beyond lat/long limits
map_NE_lands_sf_buffer <- st_crop(x = map_NE_lands_sf_buffer, y = geo_limits)
sf_use_s2(TRUE)

plot(map_NE_lands_sf_buffer)

# Check validity
table(st_is_valid(map_NE_lands_sf_buffer))
# plot(map_NE_lands_sf_buffer[!st_is_valid(map_NE_lands_sf_buffer), ])

## Save landmass buffer
saveRDS(map_NE_lands_sf_buffer, file = "./input_data/NaturalEarth_maps/map_NE_lands_sf_buffer.rds")

# ### sp version ###
# 
# ## Apply buffer ourselves because terra version in CoordinateClean::cc_sea leads to error
# 
# # Apply buffer via RGEOS because sf version sucks
# map_NE_lands_sp <- rnaturalearth::ne_load(file_name = "ne_10m_land", returnclass = 'sp', destdir = "./input_data/NaturalEarth_maps/ne_10m_land/")
# 
# plot(map_NE_lands_sp)
# 
# # Consider lat/long to be Cartesian spatial coordinates. 
# # Not strictly correct. Lead to distance deformation and issue with buffer crossing -180/180 line...
# # Better option with st_buffer in sf, which use spherical geometry, but does not seems to work...
# # Convert to a unique mutipolygon (which may fasten computation?)
# # Distance in degrees. 0.05 ≈ 5.55 km at the equator
# map_NE_lands_sp_buffer <- rgeos::gBuffer(spgeom = map_NE_lands_sp, width = 0.05)
# 
# plot(map_NE_lands_sp_buffer)
# 
# ## Adjust boundaries
# # Need to limit to a 180/90 rectangle (extent beyond when adding buffer to Antarctica for instance)
# st_bbox(map_NE_lands_sp)
# st_bbox(map_NE_lands_sp_buffer)
# 
# # geo_limits <- st_geometry(st_sfc(st_point(c(-180,-90)), st_point(c(180,90))), crs = st_crs(map_NE_lands_sf_buffer))
# # st_crs(geo_limits) <- st_crs(map_NE_lands_sf_buffer)
# # map_NE_lands_sf_buffer <- st_intersection(x = map_NE_lands_sf_buffer, y = geo_limits)
# 
# geo_limits <- SpatialPolygons(Srl = list(Polygons(srl = list(Polygon(coords = rbind(c(-180,-90), c(-180,90), c(180,90), c(180, -90), c(-180,-90)), hole = as.logical(NA))), ID = "bbox_180_90")), proj4string = raster::crs(map_NE_lands_sp_buffer))
# map_NE_lands_sp_buffer <- gIntersection(spgeom1 = map_NE_lands_sp_buffer, spgeom2 = geo_limits)
# 
# ## Simplify shapes to accelerate computation (If applied earlier, make gIntersection crash for some reason...)
# map_NE_lands_sp_buffer <- gSimplify(map_NE_lands_sp_buffer, tol = 0.1, topologyPreserve = T)
# 
# plot(map_NE_lands_sp_buffer)
# 
# # map_NE_lands_sf_buffer <- st_buffer(x = map_NE_lands_sf, dist = 0.01)
# # # map_NE_lands_sf_buffer <- st_make_valid(st_buffer(x = map_NE_lands_sf, dist = 0.1))
# # plot(map_NE_lands_sf["featurecla"])
# # plot(map_NE_lands_sf_buffer["featurecla"])
# 
# ## Convert to sf and clean out invalid features
# 
# map_NE_lands_sf_buffer <- st_as_sf(map_NE_lands_sp_buffer)
# plot(map_NE_lands_sf_buffer)
# 
# # Some features are invalid. Need to remove them?
# st_is_valid(map_NE_lands_sf_buffer)
# # map_NE_lands_sf_buffer_cast <- st_cast(map_NE_lands_sf_buffer, "POLYGON") # Disaggregate MUTLIPOLYGON into POLYGONS
# # map_NE_lands_sf_buffer_cast$valid <- st_is_valid(map_NE_lands_sf_buffer_cast) # Check validity of each polygon
# # 
# # plot(map_NE_lands_sf_buffer_cast["valid"])
# 
# ## Merge polygons to avoid overlaps due to buffering
# sf_use_s2(FALSE)
# map_NE_lands_sf_buffer <- st_union(x = st_make_valid(map_NE_lands_sf_buffer))
# sf_use_s2(TRUE)
# 
# st_is_valid(map_NE_lands_sf_buffer)
# st_bbox(map_NE_lands_sf_buffer)
# 
# plot(map_NE_lands_sf["featurecla"])
# plot(map_NE_lands_sf_buffer)
# 
# ## Remove islands smaller than 10 km² (after buffer is applied so way smaller in practice) to save computation time
# 
# map_NE_lands_sf_buffer <- st_as_sf(st_cast(map_NE_lands_sf_buffer, "POLYGON")) # Disaggregate MUTLIPOLYGON into POLYGONS
# st_geometry(map_NE_lands_sf_buffer) <- "geometry"
# 
# 
# names(map_NE_lands_sf_buffer) <- c("geometry")
# test <- map_NE_lands_sf_buffer
# test$area <- st_area(test)
# test$area <- set_units(test$area, km^2)
# 
# plot(test)
# plot(test[test$area > as_units(80, "km^2"), ])
# 
# ### sp version ###

### 2.5.2/ Run global test ####

## Load landmass buffer
map_NE_lands_sf_buffer <- readRDS(file = "./input_data/NaturalEarth_maps/map_NE_lands_sf_buffer.rds")

## Load better quality country polygons
# Countries_NE_sf <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
Countries_NE_sf <- readRDS(file = "./input_data/NaturalEarth_maps/Countries_NE_sf.rds")

## Extract entries with invalid coordinates
Biogeographic_database_Ponerinae_NA_coords <- Biogeographic_database_Ponerinae[!Biogeographic_database_Ponerinae$valid_coordinates, ]
Biogeographic_database_Ponerinae_no_NA_coords <- Biogeographic_database_Ponerinae[Biogeographic_database_Ponerinae$valid_coordinates, ]

## Save entries with NA data to retrieve later
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

## Obtain flags for each test
# Biogeographic_database_Ponerinae_flags_rerun <- clean_coordinates(x = Biogeographic_database_Ponerinae_flags_rerun,
# Biogeographic_database_Ponerinae_flags_rerun <- clean_coordinates(x = Biogeographic_database_Ponerinae_flags,
Biogeographic_database_Ponerinae_flags <- clean_coordinates(x = Biogeographic_database_Ponerinae_no_NA_coords,
                           lon = "Longitude_dec", lat = "Latitude_dec",
                           countries = "Country_ISO3_code",
                           species = "Current_name",
                           value = "spatialvalid",         # To define the type of output. Can be a spatial df ("spatialvalid"), a logical vector recording entries with at least one flag ("flagged"), a df with dubious coordinates removed ("clean").
                           verbose = T,
                           # List of tests to run
                           tests = c(#"aohi",          # Coordinates found in Artificial Hotspots
                                     "capitals",      # Coordinates found in capitals
                                     "centroids",     # Coordinates found at centroids of countries
                                     "countries",     # Test if occurrences are found within a list of countries or a custom SpatialPolygon
                                     "duplicates",    # Check for duplicates (identical coordinates) within taxa
                                     #"equal",         # Check for coordinates with latitude = longitude
                                     # "gbif",         # Check if the coordinates is within a 1km radius around the GBIF headquarter
                                     "institutions",  # Test if occurrences are found arouf known institution locations
                                     "outliers",      # Test for minimum/mean distance to other occurrences per taxa
                                     # "range"         # Test if coordinates fall within species range extracted from IUCN
                                     "seas"           # Test if fall in the ocean
                                     # "urban"          # Test if fall within urban areas
                                     # "zeros"          # Test for coordinates with null latitude or longitude
                                     ), 
                           outliers_method = "distance",  # Use minimum distance to detect outlier
                           outliers_td = 1000,            # Minimum distance to detect outlier (in km)
                           # outliers_mtp = 5.             # Multiplier for the IRQ of mean distance to flag as outlier
                           outliers_size = 3,             # Minimum number of records for a given taxa to test for outliers
                           # seas_scale = 10,               # Adjust resolution of landmass polygons used as reference. 10 = highest resolution.
                           seas_ref = map_NE_lands_sf_buffer,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                           # seas_ref = map_NE_all_lands_sf,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                           inst_rad = 200,                # Radius around Biodiversity institutions (in meters)
                           capitals_rad = 3000,           # Radius around capitals (in meters)
                           centroids_rad = 200,           # Radius around country centroids (in meters)
                           country_ref = Countries_NE_sf,  # sf SpatialPolygons to use as reference for countries
                           country_refcol = "ISO_A3",     # Name of ISO3 column
                           country_buffer = 5000,         # Buffer around country polygons (in meters)
                           # seas_buffer = 5000,            # Buffer around landmass limits (in meters) # Does not work properly with scale 10 map. Make our own prior
                           aohi_rad = 3000)               # Buffer around AOHI coordinates (in meters)

# In case of successful rerun, replace previous run
# Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags_rerun

# Check results for each test in specific columns in the output
View(Biogeographic_database_Ponerinae_flags)

## Save raw output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")
# saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags_backup.rds")

### 2.6/ Curate records for each category ####

## Load flagged output from CoordinateCleaner
Biogeographic_database_Ponerinae_flags <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")

### 2.6.1/ Prepare data ####

# Convert Flagged df to sf 
Biogeographic_database_Ponerinae_flags$Latitude_copy <- Biogeographic_database_Ponerinae_flags$Latitude_dec
Biogeographic_database_Ponerinae_flags$Longitude_copy <- Biogeographic_database_Ponerinae_flags$Longitude_dec
Biogeographic_database_Ponerinae_flags <- st_as_sf(Biogeographic_database_Ponerinae_flags, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_flags) <- sp::CRS('+init=EPSG:4326')

# Design factor level for type of flags
# Hierarchy of priorities: "zeros", "equal", "countries", "seas", "centroids", "capitals", "aohi", "institutions", "outliers", "urban", "duplicates"
Biogeographic_database_Ponerinae_flags$flag_type <- "OK"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.urb] <- "urban"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst] <- "institutions"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.aohi] <- "aohi"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap] <- "capitals"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cen] <- "centroids"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.sea] <- "seas"
Biogeographic_database_Ponerinae_flags$flag_type[!as.logical(Biogeographic_database_Ponerinae_flags$.con)] <- "countries"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.equ] <- "equal"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.zer] <- "zeros"

Biogeographic_database_Ponerinae_flags$flag_type <- as.factor(Biogeographic_database_Ponerinae_flags$flag_type)  

table(Biogeographic_database_Ponerinae_flags$flag_type)

plot(Biogeographic_database_Ponerinae_flags["flag_type"])

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")
# saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags_no_curation.rds")

### 2.6.2/ Curate records falling in wrong country ####

table(Biogeographic_database_Ponerinae_flags$.con)

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "countries", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality))

test <- Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "countries", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality)
table(test$Source) # Mostly GABI records

## 2.6.2.1/ False positive: case of records flagged as wrong, but that are actually right (tiny islands not in the shapefiles) ####

# Australia, Islands in Torres Strait, Queensland
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- "Murray Island"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- TRUE
# Australia, Islands along Queensland coast. All records are within/around tiny islands and are correct
# plot(Biogeographic_database_Ponerinae_flags[(Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "Queensland"), "adm2"])
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = ((Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "Queensland")), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = ((Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "Queensland")), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00425260")] <- "Lizard Island"
# Australia, Reversby Island, South Australia
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = ((Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "South Australia")), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = ((Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "South Australia")), replace = F)] <- TRUE
# Australia, Imperieuse Reef, West of Broome, Western Australia
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00660331")] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-031706")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00660331")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-031706")] <- TRUE
# Australia, Lord Howe Island, East of New South Wales
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00761142")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00761142")] <- TRUE
# Australia, Deception Bay, Groote Eylandt, Northern territory
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00674639")] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-066527")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00674639")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-066527")] <- TRUE

# Brazil, small islands along the coast of Sao Paulo
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$adm2 %in% c("Ilha Queimada", "Ilha de Búzios", "Ilhabela")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$adm2 %in% c("Ilha Queimada", "Ilha de Búzios", "Ilhabela")] <- TRUE

# Ecudaor, Galapagos Islands, Seymour Norte
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0173306")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0173306")] <- TRUE

# Equatorial Guinea, Isla de N. Grande
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("mncn-ent-146093")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("mncn-ent-146093")] <- TRUE

# Federated States of Micronesia: Woleai Atoll, Nama Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Woleai"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nama Is"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Woleai"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nama Is"), replace = F)] <- TRUE

# France, Port-Cros, Ile de Bagaud
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01913833", "GABI_01913918")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01913833", "GABI_01913918")] <- TRUE

# Cocos (Keeling) Islands / Indian Ocean Territories
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00760735", "GABI_00760736", "GABI_00758961", "GABI_00758962", "GABI_00758963", "GABI_00758964")] <- "Cocos (Keeling) Islands"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00760735", "GABI_00760736", "GABI_00758961", "GABI_00758962", "GABI_00758963", "GABI_00758964")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00760735", "GABI_00760736", "GABI_00758961", "GABI_00758962", "GABI_00758963", "GABI_00758964")] <- TRUE

# Indonesia, Lampung, Krakatau Islands
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Krakatau Island"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Krakatau Island"), replace = F)] <- TRUE
# Indonesia, Lampung, Krakatau Islands, Pulau Sertung
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sertung"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sertung"), replace = F)] <- TRUE
# Indonesia, Lampung, (Krakatau Islands), Rakata Besar
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rakata Besar"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rakata Besar"), replace = F)] <- TRUE
# Indonesia, Nusa Tenggara Timur, Pulau Papagarang, Manggarai Barat
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pulau Papagarang"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pulau Papagarang"), replace = F)] <- TRUE
# Indonesia, Aceh, Pulau Banyak Barat, Pulau Panjang 
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Panjang"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Panjang"), replace = F)] <- TRUE

# Italia, Isola Palmarola and Ponza
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00823282", "GABI_00818898")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00823282", "GABI_00818898")] <- TRUE

# Japan, Oshima District, Ukechima Island
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01936354", "GABI_01936357")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01936354", "GABI_01936357")] <- TRUE
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01936354", "GABI_01936357")] <- "Ukechima Island"
# Japan, Kagoshima District, Minamisatsuma, Uji Island
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01755335")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01755335")] <- TRUE
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01755335")] <- "Uji Island"
# Japan, Kagoshima District, Mishima, Takeshima Island
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00349134")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00349134")] <- TRUE
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00349134")] <- "Takeshima Island"
# Japan, Ehime District, Aoshima Island
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0915202")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0915202")] <- TRUE

# Madagascar, Nosy Ratsy, Nosy Manampaho and Nosy Mangabe on the NE coast
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Ratsy"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Manampao"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Ratsy"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Manampao"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00473283")] <- "Masoala National Park, Nosy Mangabe"
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Mangabe"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mangabé Island"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Maroantsetra"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nosy Mangabe"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mangabé Island"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Maroantsetra"), replace = F)] <- TRUE

# Malaysia, Borneo, Sabah, Mamutik Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mamutik Island"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mamutik Island"), replace = F)] <- TRUE

# Mauritius, Round Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Round Island"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Round Island"), replace = F)] <- TRUE

# Palau, Mecherchar Island
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02144550", "GABI_00679176", "GABI_00676605")] <- "Mecherchar Island"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02144550", "GABI_00679176", "GABI_00676605")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02144550", "GABI_00679176", "GABI_00676605")] <- TRUE
# Palau, Sonsorol, Merir Island
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Merir Island"
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

# PNG, Manus, Hermit group, Luf Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Luf Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Luf Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
# PNG, Morobe Province, Tami Islands
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tami") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Tami Islands"
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tami") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tami") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

# Peru, Madre de Dios, P.N. Pampas de Heath
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pampas de Heath") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pampas de Heath") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

# Saint Vincent and the Grenadines, many small islands: Baliceaux Island, Battowia Island, Catholic Island, Petit Canouan, Petit Saint Vincent Island, Savan Island
Biogeographic_database_Ponerinae_flags$.con[replace_na((Biogeographic_database_Ponerinae_flags$Locality %in% c("Baliceaux Island", "Battowia Island", "Catholic Island", "Petit Canouan", "Petit Saint Vincent Island", "Savan Island")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na((Biogeographic_database_Ponerinae_flags$Locality %in% c("Baliceaux Island", "Battowia Island", "Catholic Island", "Petit Canouan", "Petit Saint Vincent Island", "Savan Island")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

# Samoa, Nu'utele Island
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "Nu'utele Island"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "TRUE"

# Seychelles, Aldabra Atoll and Cosmoledo Island
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-015992", "casent0172609")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-015992", "casent0172609")] <- "TRUE"

# Sierra Leone, Banana Islands
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02315552")] <- "Banana Islands"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02315552")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02315552")] <- "TRUE"

# Solomon Islands, Malaupaina Island
Biogeographic_database_Ponerinae_flags$Locality[(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 161.97)] <- "Malaupaina Island"
Biogeographic_database_Ponerinae_flags$.con[(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 161.97)] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 161.97)] <- "TRUE"
# Solomon Islands, Reef Islands, Nibanga Temau
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-046469")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-046469")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- "TRUE"
# Solomon Islands, Santa Cruz Group, Vanikoro Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
# Solomon Islands, Santa Cruz Group, Anuta Island 
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
# Solomon Islands, Reef Islands, Matema Island
Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matema Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matema Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "TRUE"
# Solomon Island, Unnamed islands North of Santa Isabel
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00759029")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00759029")] <- "TRUE"

# South Korea, Maemuldo Islands
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00352044", "GABI_00352045")] <- "Maemuldo Islands"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00352044", "GABI_00352045")] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00352044", "GABI_00352045")] <- "TRUE"

# Spain, Galicia: Isla de Cíes and Isla de Ons
Biogeographic_database_Ponerinae_flags$.con[replace_na((Biogeographic_database_Ponerinae_flags$adm2_initial %in% c("Isla de Cíes", "Isla de Ons")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na((Biogeographic_database_Ponerinae_flags$adm2_initial %in% c("Isla de Cíes", "Isla de Ons")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

## Convert back to logical the "countries" and "sea" flags
Biogeographic_database_Ponerinae_flags$.con <- as.logical(Biogeographic_database_Ponerinae_flags$.con)
Biogeographic_database_Ponerinae_flags$.sea <- as.logical(Biogeographic_database_Ponerinae_flags$.sea)

# geoBoundaries_adm1_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM1/", layer = "geoBoundariesCGAZ_ADM1")
# plot(geoBoundaries_adm1_sf["shapeName"])

# geoBoundaries_adm2_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM2/", layer = "geoBoundariesCGAZ_ADM2")
# plot(geoBoundaries_adm2_sf["shapeName"]) 
# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$shapeGroup == "AUS", "shapeName"])


## 2.6.2.2/ True positive: Erroneous coordinates to correct (close to border/imprecision, error of sign, error on one digit) ####

# Antigua and Barbuda => Wrong longitude. -61.82 instead of -62.82
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02105751", "GABI_02105752", "GABI_02105753")] <- -61.82

# Argentina, Iguazu National Park. Latitude should be negative
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Iguazu National Park", "2.9 km SW of San ignacio")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Iguazu National Park", "2.9 km SW of San ignacio")])
# Argentina, San Andresito => Wrong longitude. -54.5 instead of -55.9
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "San Androcito"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -54.5
# Argentina, Reserva El Loro Hablador 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "El Loro Hablador"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -62.35
# Argentina, Parque Nacional Chaco => Wrong longitude. -59.61 instead of -41.04
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Parque Nacional Chaco"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -59.61

# Cocos/Keeling Islands. Wrong ISO country (currently in Australia) + wrong rounded coordinates
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- 96.90
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -12.116666670
Biogeographic_database_Ponerinae_flags$adm1[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "Cocos (Keeling) Islands"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "Indian Ocean Territories"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "CCI"

# Australia, Queensland => Wrong longitude. Fall in ocean. 145.65 instead of 145.76666
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00651269")] <- 145.65
# Australia, Northern Territory, Howard Springs => Wrong rounded latitude falling North in the Indian Ocean. -12.50 instead of -12.00
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Howard Springs") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -12.50
# Australia, North Western Australia => Wrong latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("TCMBW")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("TCMBW")])
# Australia, North Western Australia => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$bentity2_name_initial %in% c("South Western Australia")] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$bentity2_name_initial %in% c("South Western Australia")])

# Barbados => Wrong longitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")])
Biogeographic_database_Ponerinae_flags$bentity2_name[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")] <- "Lesser Antilles"
# Barbados => wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -59.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == 13.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -59.58
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -59.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == 13.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- 13.20

# Brazil => Wrong longitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Brazil")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Brazil")])
# Brazil, Ceara => Wrong Latitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Ceara")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Ceara")])
# Brazil, Espirito Santo, Linhares => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = ((Biogeographic_database_Ponerinae_flags$Locality_initial == "Linhares") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")), replace = F)] <- -40.04
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = ((Biogeographic_database_Ponerinae_flags$Locality_initial == "Linhares") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")), replace = F)] <- -19.41
# Brazil, Mato Grosso do Sul, Serra da Bodoquena National Park => Wrong coordinates falling in Paraguay
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Da Mata farm")), replace = F)] <- -56.74
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kadiweu reserve")), replace = F)] <- -57.25
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Salobra river")), replace = F)] <- -56.76
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Salobra river")), replace = F)] <- -20.86
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sta Maria farm-Perdido river")), replace = F)] <- -56.92
# Brazil, Minas Gerais, Juiz de Fora => Wrong latitude. Latitude should be negative. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Juiz de Fora")), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Juiz de Fora")), replace = F)])
# Brazil, Rio de Janeiro => Wrong latitude ending in the ocean. -22.919 instead of 23.919
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Rio de Janeiro") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -22.919
# Brazil, Rondônia, Porto Velho => Wrong longitude ending in the ocean. -63.904 instead of -93.904
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Port Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -63.904
# Brazil, Rondônia, Porto Velho => Latitude and Longitude are reversed!
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)]
# Brazil, Santa Catarina => Wrong signs for latitude/longitude. Should be negative + one records need to be shift back in the landmass
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Santa Catarina") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -48.553
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Santa Catarina") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.158
# Brazil, Sao Paulo, Mogi das Cruzes => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mogi das Cruzes") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -46.142
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mogi das Cruzes") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.532
# Brazil, Sao Paulo, Parque Leon Feffer => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Parque Leon Feffer") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -46.223
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Parque Leon Feffer") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.529

# Cameroon, Mbale Mejo to Ekingli => Wrong rounded coordinates
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c002505")] <- 11.717
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c002505")] <- 3.603

# China, Guangxi, Fangchenggang City, Mt. Shiwandashan => Wrong latitude, falling into China Sea => 21.91 instead of 21.18
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02357907")] <- 21.91
# China, Yunnan => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Yunnan")] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Yunnan")])
# China, Yunnan, Bakaxiaozhai, Menglun Town => Wrong coordinates falling in Vietnam
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Menglun Town") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 101.254
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Menglun Town") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 21.931


# Colombia, Antioquia, Salgar => Wrong latitude. 5.963 instead of 68.70200
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Salgar") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 5.963
# Colombia, Atlantico, Piojó => Wrong coordinates. Fall into Caribbean sea North of the coast
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Piojó") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.108
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Piojó") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.749
# Colombia, Bolivar, Los Colorados Venado => Wrong longitude. -75.116 instead of -52.116
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Los Colorados Venado") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.116
# Colombia, Cauca, PNN Gorgona Mancora => Wrong longitude. -78.183 instead of -18.183
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "PNN Gorgona Mancora") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -78.183
# Colombia, Cauca, Cundinamarca => Wrong latitude. 5.000 instead of 50.000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cundinamarca") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 5.000
# Colombia/Venezuela, Cesar, Sierra de Perija, Socorpa Mission. NP is in Venezuela so rounded coordinates are likely wrong. # Need also to change the country.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0178780"] <- -72.867407
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0178780"] <- 9.948540
# Colombia, La Guajira => Wrong latitude. Latitude should be positive. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("La Guajira")] <- abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("La Guajira")])
# Colombia, Magdalena, Santa Marta => Wrong coordinates. Fall into Caribbean sea North of the coast. 11.242 instead of 11.352
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00529958"), replace = F)] <- -74.205
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00529958"), replace = F)] <- 11.242
# Colombia, Narino, Barbacoas, RN Río Ñambí => Wrong latitude. Latitude should be positive. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Narino")] <- abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Narino")])
# Colombia, Putumayo, Puerto Leguízamo => Wrong longitude. -74.983300 instead of -79.983300
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bajo Cusacante") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -74.983300
# Colombia, Risaralda, El Trapiche => Wrong longitude. -75.954167 instead of -6.954167
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "El Trapiche") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.954167
# Colombia, Valle del Cauca, Farallones de Cali National Park => Wrong longitude. -76.656000 instead of -79.656000
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Farallones de Cali National Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -76.656000
# Colombia, Valle del Cauca, Bajo Calima, Buenaventura => Wrong longitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bajo Calima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bajo Calima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])

# Costa Rica, Alajuela => Coordinates reversed AND wrong signs! (Master class...)
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")] <- -1 * Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")]
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")] <- -1 * Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")]
# Costa Rica, Limon, Salsipuedes => Wrong latitude. Fall North to Limon in the Carribean sea. 10.007 instead of 10.07
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Salsipuedes") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.007
# Costa Rica, Puntarenas, Sirena => Wrong coordinates. Fall South in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sirena") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -83.59129
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sirena") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 8.480170
# Costa Rica, Puntarenas, Estacion Biological Las Cruces => Wrong longitude. -82.960 instead of -82.010000
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Estacion Biological Las Cruces") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -82.960
# Costa Rica, La Virgen => Wrong longitude. -83.75 instead of -8.75
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00527157"), replace = F)] <- -83.75

# Croatia, Dalmacia, Sucurac => Wrong coordinates falling in Montenegro
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sucurac") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.431
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sucurac") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 43.553
# Croatia, Hvar Island, Hvar => Wrong Longitude falling in the Adriatic Sea. 16.453 instead of 14.433333
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Hvar") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.453
# Croatia, Dalmacia, Kastela => Wrong coordinates falling in Bosnia & Herzegovina
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kastela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.350
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kastela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 43.555

# Cuba, Cumanayagua, Mina Carlota => Wrong coordinates falling in Mexico !
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cumanayagua, Mina Carlota") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -80.16667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cumanayagua, Mina Carlota") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 22.06667

# DRC, Masosa => Wrong coordinates in Congo. Not sure of the new ones but better than random.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Masosa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 22.323
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Masosa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 2.215
# DRC, Manamama => Wrong coordinates in Congo. Not sure of the new ones but better than random.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Manamama") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 19.232
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Manamama") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -4.076
# Biogeographic_database_Ponerinae_flags$.con[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Manamama") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- T

# DRC/Central African Republic, Bangassou => Wrong rounded latitude (too South) and Wrong country
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bangassou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.739

# Ecuador, Los Rios, Rio Palenque => Wrong latitude. Latitude should be negative.  
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Los Rios")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Los Rios")])
# Ecuador, Sucumbios, Limoncocha => Wrong longitude falling in Colombia. -76.6 instead of -73.6
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00845352"), replace = F)] <- -76.6
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00845352"), replace = F)] <- "Sucumbíos"

# El Salvador, Santa Ana, Montecristo => Wrong coordinates falling in Guatemala. -76.6 instead of -73.6
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_02152843"), replace = F)] <- -89.39258
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_02152843"), replace = F)] <- 14.39603

# Eritrea/Ethiopia, Nefasit => Wrong rounded coordinates + Wrong country: Eritrea instead of Ethiopia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_01942955"), replace = F)] <- 39.065
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_01942955"), replace = F)] <- 15.331
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_01942955"), replace = F)] <- "Nefasit"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01942955")] <- "Eritrea"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01942955")] <- "ERI"

# Micronesia, Tol Island, Mt. Unibot => Wrong coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mt. Unibot") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.628
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mt. Unibot") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.368
# Micronesia, Pohnpei, Ponape Agriculture & Trade School => Wrong rounded coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Ponape Agriculture & Trade School") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 158.308
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Ponape Agriculture & Trade School") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.842
# Micronesia, Woleai Atoll => Wrong coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Woleai Atoll") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 143.913
# Micronesia, Nama Island => Wrong coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nama Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 152.576
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nama Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.992

# French Guyana, Campus Agronomique in Kourou => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("French Guiana")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("French Guiana")])
# French Guyana, Sinnamary, Petit Saut => Wrong coordinates in Brazil
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00530228"), replace = F)] <- -53.050
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00530228"), replace = F)] <- 5.064

# Gabon, Aire d'Exploition Rationnelle de Faune des Monts Doudou => Wrong Latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Monts Doudou")), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Monts Doudou")), replace = F)])

# Greece, Corfú => Wrong latitude falling in the Ionian sea South to Corfu. 39.456667 instead of 39.286667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Corfú")), replace = F)] <- 39.456667
# Greece, Macedonia, Pieria district, Olympus Mts. Litochoro, Faragi Enipeas => Wrong coordinates falling in the Egean Sea 
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469648", "GABI_00469649")), replace = F)] <- 40.105
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469648", "GABI_00469649")), replace = F)] <- 22.490
# Greece, Macedonia, Pieria district, Olympus Mts. Litochoro, Leptokaria-Karia road => Wrong coordinates falling in the Egean Sea 
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469724")), replace = F)] <- 40.03333
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469724")), replace = F)] <- 22.51666
# Greece, Macedonia, Pieria district, Platamonas Castle => Wrong coordinates falling in the Egean Sea 
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469875")), replace = F)] <- 40.005
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00469875")), replace = F)] <- 22.598

# Guyana, Upper Demerara-Berbice, Mabura Hill => Wrong longitude. Longitude should be negative.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("Guyana")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("Guyana")])

# India, Arunachal Pradesh, Lumla => Wrong coordinates falling in Bhoutan
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01662113")), replace = F)] <- 91.722
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01662113")), replace = F)] <- 27.530

# Indonesia, Sulawesi, Central Sulawesi Province, Berdikari => Wrong Latitude falling North in the Celebes Sea. 0.734 instead of 1.134
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Berdikari") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 0.734
# Indonesia, Sulawesi, Central Sulawesi Province, Bulili => Wrong Longitude falling East in the Celebes Sea. 120.0955 instead of 120.955
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bulili") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 120.0955
# Indonesia, East Nusa Tenggara Province, Nangagete => Wrong coordinates, falling SE in the Savu Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_02008543"), replace = F)] <- 122.2136
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_02008543"), replace = F)] <- -8.6771
# Indonesia, Lampung, Sunda Strait, Krakatau Islands, Pulau Panjang => Wrong rounded longitude, falling East in the sea. 105.38 instead of 105.4
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "anic32-066259"), replace = F)] <- 105.38
# Indonesia, Papua, Jayapura => Wrong longitude falling in the PNG side of the border
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "fmnhins0000117848"), replace = F)] <- 140.993
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_02199373"), replace = F)] <- 140.993
# Indonesia, West Papua, Sorong => Wrong coordinates falling North in the Pacific ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0281867"), replace = F)] <- 131.418
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0281867"), replace = F)] <- -1.039
# Indonesia, Aceh, Pulau Panjang => Wrong longitude falling West in the Indian ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-066249", "anic32-066250")), replace = F)] <- 97.416
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-066249", "anic32-066250")), replace = F)] <- "Pulau Panjang, Pulau Banyak Barat, Aceh"
# Indonesia, Palau Laut => Wrong coordinates falling in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pulo Laut") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 107.977
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pulo Laut") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.715

# Ivory Coast, Tai National Park => Wrong longitude falling in Liberia. -7.09642 instead of -7.89642 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Foret de Tai") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -7.09642
# Ivory Coast, Tai Forest => Wrong rounded longitude falling in Liberia. -7.066667 instead of -7.666667
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tai Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -7.066667

# Japan, Iwate Prefecture => Wrong Latitude falling South in the Pacific Ocean. 39.012222 instead of 33.012222
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Iwate Prefecture") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 39.012222
# Japan, Kagoshima Prefecture, Kagoshima City, Jigenji Park => Wrong coordinates falling South in the Pacific ocean 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Jigenji Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 130.503
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Jigenji Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 31.511
# Japan, Kagoshima Prefecture, Sakurajima => Wrong coordinates falling South in the Pacific ocean 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sakurajima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 130.650
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sakurajima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 31.583
# Japan, close to Akusekijima Island => Wrong coordinates falling East of Akusekijima Island
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00577065", "GABI_00577066")), replace = F)] <- 129.604
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00577065", "GABI_00577066")), replace = F)] <- 29.459

# Lesotho, Likhoele => Wrong Latitude falling North in South Africa. -29.81667 instead of -29.183333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Likhoele") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -29.81667
# Lesotho, Mokhotlong => => Wrong coordinates falling East in South Africa.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mokhotlong") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.08333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mokhotlong") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -29.71667
# Lesotho, Qachas Nek => => Wrong latitude falling South in South Africa. -30.1 instead of -30.9
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Qachas Nek") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -30.1

# Madagascar => Wrong Latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Madagascar")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Madagascar")])
# Madagascar, Marojejy National Park => Wrong Longitude falling West in the Mozambique Canal. 49.934722 instead of 46.934722
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00557292"), replace = F)] <- 49.934722

# Malaysia, Kelantan, Melawi Beach => Wrong coordinates falling NW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Melawi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 102.420
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Melawi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.021
# Malaysia, Borneo, Sarawak, Kampong Sega => Wrong coordinates falling in Kalimantan, Indonesia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kampong Segu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 110.238
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kampong Segu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 1.424

# Malaysia, Borneo, Sarawak, Semengoh Forest Reserve => Wrong Latitude falling in Kalimantan, Indonesia. 1.401 instead of 0.8833333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Semengoh Forest Reserve") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 1.401
# Malaysia, Borneo, Sarawak, Gunung Mulu National Park => Wrong Latitude falling in Brunei. 4.037 instead of 4.95
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Gunung Mulu National Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.037
# Malaysia, Borneo, Sarawak, Lambir Hills National Park => Wrong coordinates falling in Kalimantan, Indonesia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lambir Hills") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 114.04
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lambir Hills") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.200
# Borneo, Tobang => Coordinates are centroids of Borneo. Could relate to Tabang in Kalimantan, Indonesia.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 116.017
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 0.575
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "IDN"

# Mali, Baboye => Wrong coordinates falling South in Burkina Faso
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Baboye") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.952
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Baboye") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 14.260

# Mauritius, Round Island => Wrong Longitude falling West in the Indian ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Round Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 57.784

# Mexico, Chiapas, Irlanda => Wrong coordinates falling in Guatemala
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Irlanda") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -93.332
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Irlanda") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.109
# Mexico, Chiapas, Soconusco => Wrong coordinates falling in Guatemala
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01103516", "GABI_01103514"), replace = F)] <- "Soconusco region"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soconusco") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -92.297
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soconusco") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 15.144
# Mexico, Jalisco, Estación Biológica Chamela. Wrong Latitude falling South into the Pacific Ocean 19.50 instead of 19.05
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Estación Biológica Chamela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 19.50
# Mexico, Quintana Roo, Majahual => Wrong coordinates falling NE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Majahual") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -87.71013
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Majahual") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.71277
# Mexico, Quintana Roo, Puerto Morelos => Wrong rounded coordinates falling SE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Puerto Morelos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -86.90278
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Puerto Morelos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 20.84653 
# Mexico, Quintana Roo, Centro de Investigaciones Costeras La Mancha => Wrong coordinates with sign inversion falling in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01752436", "GABI_01752437"), replace = F)] <- -96.38125
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01752436", "GABI_01752437"), replace = F)] <- 19.59494

# Mozambique, Maputo, Delagoa Bay => Wrong rounded longitude falling East in the Bay
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Delagoa Bay") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 32.561

# Namibia, Oidimba Village => Wrong coordinates falling NW into Angola
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Edimba") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.496
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Edimba") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.476
# Namibia, Katima Mulilo => Wrong coordinates falling North into Zambia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Katima Mulilo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 24.26667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Katima Mulilo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.50

# Palau, Ngatpang => Wrong rounded coordinates falling in SW in the Pacific
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Ngatpang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 134.53244
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Ngatpang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.45053
# Palau, Peleliu => Wrong rounded coordinates falling in SW in the Pacific
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 134.24308
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.01506
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Peleliu"

# Panama, STRI Gamboa, Pipeline Road => Wrong Longitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Panama")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Panama")])
# Panama, Nusagandi => Wrong rounded coordinates in Costa Rica!
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nusagandi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -78.97927
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nusagandi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 9.351060

# PNG, Bougainville Island, Buin Village
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buin"), replace = F)] <- 155.688
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buin"), replace = F)] <- -6.745

# PNG, Hermit group, Luf Island => Wrong Longitude falling East in the Pacific Ocean. 145.063 instead of 145.083
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Luf Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 145.063
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00654783"), replace = F)] <- 145.063
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00654783"), replace = F)] <- "Hermit group, Luf Island"
# PNG, Manaus, Lorengau => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lorengau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 147.272
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lorengau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -2.031
# PNG, Central Province, Wanigela => Wrong Latitude falling South in Pacific Ocean
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Wanigela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.062
# PNG, Central Province, Bisianumu, near Sogeri => Wrong Latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bisianumu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bisianumu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])
# PNG, East New Britain Province, Lamas => Wrong coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lamas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.404
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lamas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -5.6098
# PNG, East New Britain Province, Vouvou => Wrong coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Vouvou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.4594
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Vouvou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -5.440
# PNG, Sandaun Province, Aitape => Wrong rounded coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Aitape") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 142.35
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Aitape") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.13333 
# PNG, Morobe Province, Lae => Wrong coordinates falling far away
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lae") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 147.000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lae") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -6.73333
# PNG, Mussau Talumalaus => Wrong rounded coordinates falling South in the Bismark Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mussau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 149.621
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mussau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.436
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 150) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 149.621
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 150) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.436
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 150) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Bismarck Archipelago, Mussau Island"
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 150) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "New Ireland"
# PNG, Western Province, Fly River => Wrong longitude falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Fly River") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 143.629

# Peru, Madre de Dios, Sachavacayoc Center => Wrong Latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sachavacayoc") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sachavacayoc") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])

# Philippines, no locality => Wrong Longitude falling West of the coast of South Luzon. Bring them back on the land assuming minimal error in coordinates. 124.13750 instead of 124.26750
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == 12.877214) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 124.13750

# Puerto Rico, Bosque Estatal Guanica => Wrong Longitude falling South in the Caribbean Sea. 17.97 instead of 17.84
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bosque Estatal Guanica") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 17.97
# Puerto Rico, Carite Lake region => Wrong rounded coordinates falling SE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Carite") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -66.100
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Carite") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.074

# Rwanda, Rubona => Wrong coordinates falling SW in DRC
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Rubona") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.256
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Rubona") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.728

# Saint Vincent and the Grenadines, Savan Island => Wrong Longitude falling East in the Atlantic Ocean. -61.21111 instead of -61.12111
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Savan Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.21111

# Samoa => Wrong Latitude. Should be negative
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Country_ISO3_name, pattern = "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Country_ISO3_name, pattern = "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])
# Samoa => Wrong Latitude falling South in the Pacific ocean. -13.583 instead of -16.583
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00761072", "GABI_00761107"), replace = F)] <- -13.583

# Saudi Arabia, Al Qatif => => Wrong Latitude falling North in the Arabic Sea. 26.533333 instead of 26.933333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Al Qatif") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 26.533333

# Solomon Islands, Isabel Province => Wrong rounded longitude falling West in the Pacific Ocean. 159.115 instead of 159.5
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Isabel Province") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 159.115
# Solomon Islands, Cristobal Island, Wainoni => Wrong Latitude falling North in the Pacific Ocean. -10.566670 instead of -10.066670/-10.100000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Wainoni") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.566670
# Solomon Islands, Santa Cruz Group, Nendo Island => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Locality[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -11) & is.na(Biogeographic_database_Ponerinae_flags$Locality_initial) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Santa Cruz Group, Nendo Island"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -11) & (Biogeographic_database_Ponerinae_flags$Locality == "Santa Cruz Group, Nendo Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 165.93927
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -11) & (Biogeographic_database_Ponerinae_flags$Locality == "Santa Cruz Group, Nendo Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.718873
# Solomon Islands, Reef Islands, Nibanga Temau => Wrong Longitude falling East in the Pacific Ocean. 166.310 instead of 166.167
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-046469")] <- 166.310
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- "Reef Islands, Nibanga Temau"
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- 166.310
# Solomon Islands, Reef Islands, Matema Island => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Matema") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 166.184
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Matema") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.292
# Solomon Islands, Santa Cruz Group, Vanikoro Island => Wrong rounded coordinates falling NW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 166.909
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -11.674
# Solomon Islands, Santa Cruz Group, Anuta Island => Wrong rounded coordinates falling SW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 169.850
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -11.611


# South Africa, East London, Peddie District => Wrong Latitude falling South into the Ocean. -33.03333 instead of -33.125
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "East London") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -33.03333
# South Africa, Rainwoods (Small Holding), Near Witelsbos => Wrong Latitude falling South into the Ocean. -34.000 instead of -34.125000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Witelsbos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -34.000
# South Africa, Nongoma, Zululand => Wrong Latitude falling North in eSwatini. -27.9 instead of -27.1
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nongom") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.9
# South Africa, Ntendeka Forest Reserve => Wrong Latitude falling North in eSwatini. -27.85 instead of -27.15
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Ntendeka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.85
# South Africa, Umtamvuna Nature Reserve => Wrong Latitude falling South in the Indian Ocean. -31.031 instead of -31.228889
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Umtamvuna") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -31.031
# South Africa, Entabeni Forest Reserve => Wrong Latitude falling North in Botswana. -22.93333 instead of -22.066667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Entabeni") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -22.93333
# South Africa, Soutpansberg Mountain => Wrong coordinates falling NW in Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Soutpansberg") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.42801
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Soutpansberg") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.03564
# South Africa, Eastern Transvaal, Nelspruit => Wrong latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Nelspruit"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec_initial[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nelspruit"), replace = F)])
# South Africa, Thursford Farm => Wrong longitude. Should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Thursford"), replace = F)] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec_initial[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Thursford"), replace = F)])

# Zambia, Lusaka => Wrong coordinates falling NW in Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 28.31938
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -15.43585

# South Korea, Site 114 => Wrong Longitude, falling East in the Pacific Ocean. 129.32 instead of 129.42.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "114") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.32
# South Korea, Site 115 => Wrong Longitude, falling East in the Pacific Ocean. 129.319167 instead of 129.419167.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "115") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.319167
# South Korea, Site 116 => Wrong Longitude, falling East in the Pacific Ocean. 129.323889 instead of 129.433889.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "116") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.323889
# South Korea, Site 51 => Wrong Longitude, falling East in the Pacific Ocean. 128.425 instead of 128.535.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "51") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 128.425
# South Korea, Site 58 => Wrong Longitude, falling East in the Pacific Ocean. 128.529167 instead of 128.629167.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "58") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 128.529167
# South Korea, Site 159 => Wrong Longitude, falling East in the Pacific Ocean. 129.345278 instead of 129.445278.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "159") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.345278

# Spain, Illes Balears, Mallorca => Wrong coordinates falling in the Mediterranean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00794898", "GABI_00794899"), replace = F)] <- 2.986
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00794898", "GABI_00794899"), replace = F)] <- 39.615
# Spain, Illes Canarias, El Hierro => Wrong Latitude falling North in the Atlantic Ocean. 27.75 instead of 27.85
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "El Hierro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 27.75
# Spain, Illes Canarias, Gran Canaria Island, Las Palmas => Wrong Longitude falling East in the Atlantic Ocean. -15.428611 instead of -15.318611
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Las Palmas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -15.428611
# Spain, Illes Canarias => Wrong coordinates falling in the middle of the archipelago, in the Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00831742"), replace = F)] <- -15.70393
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00831742"), replace = F)] <- 28.10043
# Spain, Arroyo de la Ermita, Sierra Almijara => Wrong coordinates falling SE in the Mediterranean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Arroyo de la Ermita") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.969
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Arroyo de la Ermita") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 36.894

# São Tomé and Principe, São Tomé Island => Wrong Latitude falling North in the Atlantic Ocean. 0.26492 instead of 0.56492
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "ey20045"), replace = F)] <- 0.26492

# Trinidad & Tobago, "Trinidad_b" => Falling in the Atlantic Ocean between the two islands. Take random coordiantes on Trinidad instead.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_code, pattern = "Trinidad_b") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.236
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_code, pattern = "Trinidad_b") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.416

# USA, California, Cold Canyon => Wrong latitude falling South into the Pacific Ocean. 38.511 instead of 36.511.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cold Canyon") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 38.511
# USA, Florida, Franklin County => Wrong rounded latitude falling South into the Caribbean Sea Ocean. 29.800 instead of 36.511.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Franklin County") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.916666
# USA, Florida, Seminole State Park => Wrong coordinates falling SE into the Carribean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Seminole State Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -81.604
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Seminole State Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 25.976
# USA, Georgia, Cloud Land Canyon => Wrong Latitude falling South into the Carribean Sea. 34.834444 instead of 24.834444.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cloudland Canyon") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 34.834444
# USA, Georgia, Tugaloo => Wrong Latitude falling South into the Carribean Sea. 34.494722 instead of 24.484722.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tugaloo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 34.494722
# USA, Mississippi, Grand Bay Savanna => Wrong Latitude falling South into the Carribean Sea. 30.36 instead of 30.0475.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Grand Bay Savanna") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 30.36
# USA, Texas, Harlingen => Wrong Latitude falling South into the Carribean Sea. 26.165 instead of 23.6.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Harlingen") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 26.165
# USA, Hawaii, Oahu Island => Wrong rounded Latitude falling South into the Carribean Sea. 21.5 instead of 21.0.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Oahu Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 21.5

# Vanuatu, Efate Island, Port Vila => Wrong rounded coordinates falling NW into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Port Fila") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 168.32
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Port Fila") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.736

# Venezuela, Valle del Abismo => Wrong coordinates falling SE in Brazil
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Valle del Abismo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.6
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Valle del Abismo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.40
# Venezeuala, Azevedo District, Cupo => Wrong Longitude falling East in the Caribbean Sea. -66.372 instead of -65.622
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cupo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -66.372
# Venezeuala, Nueva Esparta, Isla Margarita, ca. Fuentidueño => Wrong Longitude falling West in the Caribbean Sea. -63.913056 instead of -64.913056
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Fuentidueño") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -63.913056
# Venezuela, San Critobal, Loma de Pio => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Loma de Pio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.200
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Loma de Pio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.749
# Venezuela, San Antonio del Tachira => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "San Antonio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.442
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "San Antonio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.818
# Venezuela, Tachira, San Cristobal, Santa Teresa => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Santa Teresa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.225
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Santa Teresa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.800

# Viet Nam, Huong Son => Wrong coordinates falling SE in the China Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Huong Son") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 105.426
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Huong Son") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.512
# Viet Nam, Ba Vi National Park => Wrong Longitude falling East in the China Sea. 105.36083 instead of 108.36083
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Ba Vi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 105.36083

# Wallis and Futuna, Wallis Island => Wrong coordinates falling NE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00758995", "GABI_00754808"), replace = F)] <- -176.202
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00758995", "GABI_00754808"), replace = F)] <- -13.285

# Zambia, Mpulungu => Wrong rounded Latitude falling North in Lake Tanganyika. -8.83333 instead of -8.500
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mpulungu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -8.83333

# Zimbabawe, Matopo Hills => Wrong Latitude falling South in South Africa. -20.45 instead of -27.45
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Matopo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -20.45

# eSwatini, Pigg's Peak => Latitude falling North in South Africa. -26.033333 instead of -25.033333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Piggspeak") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -26.033333


Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Salgar")), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Salgar")), replace = F)]

# plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Australia"), "ISO_A3"]) # Includes the Canarias
# plot(Countries_NE_sf[, "ISO_A3"])
# 
# plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Palestine"), "ISO_A3"]) # Includes the Canarias

## 2.6.2.3/ True positive: Need to change country name ####

# American Samoa => Records are in the Western Samoa Islands, not the American Samoa
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "WSM"

# Norfolk Island => Records are in Norfolk Island, not Australia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Norfolk Island")] <- "Norfolk Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Norfolk Island")] <- "NFK"

# Botswana => Records and localities are in the South African side of the border
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00833495", "GABI_01954574")] <- "South Africa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00833495", "GABI_01954574")] <- "ZAF"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c016271", "sam-hym-c013580", "sam-hym-c013581")] <- "South Africa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c016271", "sam-hym-c013580", "sam-hym-c013581")] <- "ZAF"

# Hong Kong instead of China
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Hong Kong")] <- "Hong Kong S.A.R."
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Hong Kong")] <- "HKG"

# Macao instead of China
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Macau")] <- "Macao S.A.R."
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Macau")] <- "MAC"

# Sierra de Perija. NP and coordinates in Venezuela, not Colombia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = paste0("Sierra de Perija", "|","Sierra de Parija")), replace = F)] <- "Venezuela"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = paste0("Sierra de Perija", "|","Sierra de Parija")), replace = F)] <- "VEN"

# Eritrea instead of some Ethiopia records
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "Tesseney"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "Eritrea"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "ERI"

# Clipperton Island instead of France
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Clipperton Island")] <- "Clipperton Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1_initial %in% c("Clipperton Island")] <- "CPT"

# Guyana instead of French Guyana/France for an erroneous record
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00845045")] <- "Guyana"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00845045")] <- "GUY"

# French Southern and Antarctic Lands instead of Europa Island/France for Europa Island
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality_initialy %in% c("Europa Island")] <- "French Southern and Antarctic Lands"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Europa Island")] <- "ATF"

# Jordan instead of Israel. Coordinates falling in the Jordan side of the river.
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Jourdain")] <- "Jordan"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Jourdain")] <- "JOR"

# Sancurar in Ethiopia instead of Kenya
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00477580")] <- "Ethiopia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00477580")] <- "ETH"

# Jericho in Palestine instead of Lebanon
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798847")] <- "Palestine"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798847")] <- "PSE"
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798847")] <- "Jericho"

# Al-Mazra'a ash-Sharqiya in Palestine instead of Lebanon
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798877")] <- "Palestine"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798877")] <- "PSE"
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00798877")] <- "Al-Mazra'a ash-Sharqiya"

# ?Unknown are in Martinique/France
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("?Unknown")] <- "France"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("?Unknown")] <- "FRA"
Biogeographic_database_Ponerinae_flags$adm1[Biogeographic_database_Ponerinae_flags$Country_initial %in% c("?Unknown")] <- "Martinique"

# Namakunde in Angola instead of Namibia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Namakunde")] <- "Angola"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality_initial %in% c("Namakunde")] <- "AGO"
# Erikson's Drift, Kunene River in Angola instead of Namibia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01943756")] <- "Angola"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01943756")] <- "AGO"

# Kongga/Buin Village in PNG instead of Solomon Islands
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buin"), replace = F)] <- "Autonomous Region of Bougainville [North Solomons]"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buin"), replace = F)] <- "Papua New Guinea"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buin"), replace = F)] <- "PNG"

# Sorong in Indonesia, West Papua instead of PNG
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sorong"), replace = F)] <- "West Papua"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sorong"), replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sorong"), replace = F)] <- "IDN"

# Bali, Laba Sari in Indonesia instead of PNG
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "Bali, Laba Sari"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "IDN"

# Castelnuovo in Montenegro instead of Serbia/Croatia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Castelnuovo"), replace = F)] <- "Montenegro"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Castelnuovo"), replace = F)] <- "MNE"

# Sona Mpungu Forest in Democratic Republic of the Congo instead of Republic of the Congo
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sona Mpungu Forest"), replace = F)] <- "Democratic Republic of the Congo"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sona Mpungu Forest"), replace = F)] <- "COD"

# Bangassou in Central African Republic instead of Equatorial Guinea
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bangassou"), replace = F)] <- "Central African Republic"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bangassou"), replace = F)] <- "CAF"

# Fond St Jacques in Martinique/France instead of Saint Lucia
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Fond St Jacques"), replace = F)] <- "Martinique"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Fond St Jacques"), replace = F)] <- "France"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Fond St Jacques"), replace = F)] <- "FRA"

# American Samoa instead of Samoa
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -170.70), replace = F)] <- "Anu'u Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -170.70), replace = F)] <- "American Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -170.70), replace = F)] <- "ASM"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -169.45), replace = F)] <- "Ta'u Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -169.45), replace = F)] <- "American Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == -169.45), replace = F)] <- "ASM"

# Somalia to Somaliland: [Ouarsangueli] Warsangali, Ainabo, Shimbir Beris, Ceelbuh, Ina Dhakool
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00835211"), replace = F)] <- "Ceelbuh"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00833795", "GABI_00833796")), replace = F)] <- "Ina Dhakool"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality %in% c("[Ouarsangueli] Warsangali", "Ainabo", "Shimbir Beris", "Ceelbuh", "Ina Dhakool")), replace = F)] <- "Somaliland"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality %in% c("[Ouarsangueli] Warsangali", "Ainabo", "Shimbir Beris", "Ceelbuh", "Ina Dhakool")), replace = F)] <- "SOL"

# Lusaka in Zambia instead of South Africa
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Zambia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "ZMB"

# Eastern Equatoria, Lotti Forest, SW Slope Achali Mts. in South Sudan instead of Sudan
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Eastern Equatoria") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "South Sudan"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Eastern Equatoria") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "SSD"

# West Caprivi Park, Manywa River in Namibia instead of Zambia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "West Caprivi Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Namibia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "West Caprivi Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "NAM"

# From Lochiel to Mbabane. Point on South African side of the border.
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01952400")] <- "South Africa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01952400")] <- "ZAF"


## 2.6.2.4/ True positive: False coordinates to remove (turn into NA) ####

# Check if locality is credible...

# DRC, Haut-Ubangi => Wrong coordinates in Congo. To draw at random in the adm1.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- F
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Nord Ubangi"

# Fiji => Wrong rounded coordinates falling in the Pacific Ocean # Draw random coordinates on the Fiji islands
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_dec_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- F

# India, West Bengal => Wrong rounded coordinates falling in Bangladesh # Draw random coordinates in the West Bengal region
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0902607"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0902607"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0902607"), replace = F)] <- "West Bengal"
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0902607"), replace = F)] <- F

# Japan, Kagoshima Prefecture => Wrong coordinates falling in the Pacific Ocean # Draw random coordinates in the Kagoshima Prefecture
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0217561"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0217561"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0217561"), replace = F)] <- "Kagoshima Prefecture"
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0217561"), replace = F)] <- F

# Malaysia, Borneo => Coordinates at centroid of the island falling in Kalimantan, Indonesia but registered in Sarawak, Malaysia # Better to draw them randomly and label them as random
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0913723", "casent0903954", "casent0903955")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0913723", "casent0903954", "casent0903955")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0913723", "casent0903954", "casent0903955")), replace = F)] <- F

# New Zealand => Coordinates right in between the two main islands without further indications. # Drawn random coordinates in NZ
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00757929", "GABI_00757931", "GABI_00756746", "GABI_00757928")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00757929", "GABI_00757931", "GABI_00756746", "GABI_00757928")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00757929", "GABI_00757931", "GABI_00756746", "GABI_00757928")), replace = F)] <- F

# Philippines => Wrong rounded coordinates falling in the middle of the archipelago. # Draw random coordinates in the Philippines.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0907422", "casent0915243")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0907422", "casent0915243")), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("casent0907422", "casent0915243")), replace = F)] <- F

## 2.6.2.5/ Update Geometry ####

# Convert Flagged df to sf 
Biogeographic_database_Ponerinae_flags$Latitude_copy <- Biogeographic_database_Ponerinae_flags$Latitude_dec
Biogeographic_database_Ponerinae_flags$Longitude_copy <- Biogeographic_database_Ponerinae_flags$Longitude_dec
Biogeographic_database_Ponerinae_flags <- st_as_sf(Biogeographic_database_Ponerinae_flags, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_flags) <- sp::CRS('+init=EPSG:4326')
# Remove copies of Long/Lat of still present
Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags %>% select(-c("Longitude_copy", "Latitude_copy"))
Biogeographic_database_Ponerinae_flags

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")

table(Biogeographic_database_Ponerinae_flags$flag_type) # Need to rerun evaluation to see if all tagged entries are cleaned

### 2.6.3/ Curate records falling in water areas ####

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "seas", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

## 2.6.3.1/ False positive: case of records flagged as at seas, but that are actually on lands (tiny islands or deltas) ####

# Brazil, Ilha dos Porcos, Amazon Delta
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00753526")] <- TRUE

# Costa Rica, Osa, Sierpe, Isla del Caño
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Osa, Sierpe") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Osa, Sierpe, Isla del Caño"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Isla del Caño") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Isla de Caño") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Italy, Archipelago Pontino, Ventotene Island
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Ventotene") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Japan, Kagoshima Prefecture, Minamisatsuma, Kamino Island
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00619696", "GABI_00617902")] <- TRUE
# Japan, Okinawa Prefecture, Iou-torishima Island
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00349103")] <- TRUE

# Palau, Sonsorol, Sonsorol island
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial < 5.34) & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Sonsorol, Sonsorol Island"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$Latitude_dec_initial > 5.34) & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Sonsorol, Fanna Island"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1_initial, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Seychelles, Aride Island
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Aride") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Spain, Galicia: Isla de Cíes and Isla de Ons
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm2_initial, pattern = "Isla de Ons") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Islas Cíes") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# USA, Keys Islands
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = " Key") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

## Convert back to logical the "countries" and "sea" flags
Biogeographic_database_Ponerinae_flags$.sea <- as.logical(Biogeographic_database_Ponerinae_flags$.sea)


##  2.6.3.2/ True positive: Slightly erroneous coordinates to correct (Fall in internal lakes, close to land, just imprecision) ####

# Brazil, Rio Grande do Sul, Parque Estadual de Itapua, Viamao => Wrong coordinates falling WE in Lake/Lagoa dos Patos
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Itapua") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -51.026
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Itapua") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -30.347

# Colombia, Magdalena, Cienaga Grande de Sta Maria => Wrong coordinates in the middle of the Lake
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cienaga Grande") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -74.324

# Greece, Evrostina => Wrong coordinates falling into the Corinth Gulf
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Evrostina") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- 22.396
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Evrostina") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- 38.073

# Portugal, Leça da Palmeira => Wrong Longitude, falling West in the Atlantic Ocean. -8.701389 instead of -8.761389
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Leça da Palmeira") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -8.701389

# Seychelles, Aride Island => Wrong rounded Latitude, falling North in the Pacific Ocean. -4.212 instead of -4.200
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Aride") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -4.212


# Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cienaga")), replace = F)]
# Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Cienaga")), replace = F)]

# Need to update the geometry to include new coordinates
# Rerun flag check to ensure changes have fixed them
# Modify again the false positives
# Need to check if they have a sea flag, or any other flag, too to change their flag

### 2.6.4/ Curate records in around country centroids ####

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "centroids", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# Check if locality is available to retrieve true coordinates
# Check if uncertainty is too large to be kept (may be okay to keep those for small islands or countries)

## 2.6.4.1/ False positive: Records in localities that are close to country centroids, or centroids of small islands so it does not matter ####

# Cook Islands, Rarotonga
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Rarotonga"), replace = F)] <- TRUE

# Martinique
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country_initial == "Martinique") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Haiti
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country_initial == "Haiti") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Christmas Island
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country_initial == "Christmas Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Santa Lucia, Praslin
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Praslin"), replace = F)] <- TRUE

# Saint Martin, Loterie Farm
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Loterie Farm"), replace = F)] <- TRUE

# Samoa Island
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country_initial == "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Trinidad and Tobago, Lopinot Complex
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Lopinot complex"), replace = F)] <- TRUE


## 2.6.4.2/ True positive: Records that have wrong coordinates according to their locality => To adjust ####

# Cameroun, Mount Cameroun above Buea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Buea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 9.214
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 4.170

# Czech Republic, Radotínské údolí
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Radotínské") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 14.314
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Radotínské") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 49.997

# Equatorial Guinea, Bioko
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bioko") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 8.70
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bioko") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 3.50

# Ethiopia, Gughé highlands, Bonghé Valley
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bonghé") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 37.35
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bonghé") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 6.06
# Ethiopia, Dire Dawa
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Dire-Dawa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 41.857
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Dire-Dawa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 9.608

# Fiji, Sigatoka
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sigatoka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 177.553
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Sigatoka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -18.054
# Fiji, Kings Road
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kings Rd") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 178.437
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kings Rd") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -17.783

# Malawi, Mlanje
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 35.508
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -16.024
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje_2000") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 35.568
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje_2000") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -15.976

# Mexico, Peñon del Marquis
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Peñon del Marquis") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -99.017
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Peñon del Marquis") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 19.374
# Mexico, Pedrigales
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pedrigales") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -99.207
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pedrigales") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 19.322

# Mozambique, Amatongas Forest
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Amatongas") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 33.768
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Amatongas") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -19.182

# Myanmar, Bhamo
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bhamò") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 97.234
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Bhamò") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 24.262

# Taiwan, Akau
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Akau") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.473
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Akau") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 22.713
# Taiwan, Kosempo
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Jiaxian District [Kosempo, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.614
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.121
# Taiwan, Pilam
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Taitung District [Pilam, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.126
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 22.759
# Taiwan, Taihorin
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Dalin [Taihorin, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.455
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.604
# Taiwan, Taihoku/Taipeh
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Taipei [Taihoku, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.538
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 25.038
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taipeh") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.538
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Taipeh") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 25.038
# Taiwan, Suisharyo
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Yuchi [Suisharyo, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.920
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.872

# Tanzania, Tanganyika
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tanganyika") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 29.958
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Tanganyika") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -6.246

# Uganda, Busnia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Busnia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 34.089
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Busnia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 0.467
# uganda, Mabira Forest
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mabira Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 33.009
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Mabira Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 0.459

# Zimbabwe, Matsuma Dam
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Matsuma Dam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 26.282
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Matsuma Dam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -18.733


# Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau")), replace = F)]
# Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau")), replace = F)]

## 2.6.4.3/ True positive: Records with no usable information => To assign to NA for random drawing ####

# Bolivia. -17.00, -65.0
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Bolivia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Bolivia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Bolivia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Brazil. -10.00, -55.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Brazil") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Brazil") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Brazil") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Cameroun. 6.00, 12.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Cameroon") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Cameroun") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Cameroon") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Cameroun") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Cameroon") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Cameroun") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# China. 35.00, 105.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "China") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "China") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "China") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Colombia. 4.00, -72.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Colombia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Colombia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Colombia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Costa Rica. 10.00, -84.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Costa Rica") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Costa Rica") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Costa Rica") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Ecuador. -2.00, -77.50
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ecuador") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ecuador") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ecuador") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Ethiopia. 8.00, 38.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ethiopia") & (Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Abyssinie", "Ethiopia")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ethiopia") & (Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Abyssinie", "Ethiopia")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Ethiopia") & (Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Abyssinie", "Ethiopia")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Gabon. -1.00, 11.75
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Gabon") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Gabon") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Gabon") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Guinea. 11.00, -10.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Guinea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Guinea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Guinea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# India. 20.00, 77.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "India") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "India") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "India") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Iran. 32.00, 53.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Iran") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Iran") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Iran") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Japan. 36.00, 138.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Japan") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Japan") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Japan") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Liberia. 6.50, -9.50
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Liberia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Liberia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Liberia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Mexico. 23.00, -102.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Mexique", "Mexico_Neotropical", "Mexico_Nearctic")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Mexique", "Mexico_Neotropical", "Mexico_Nearctic")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Mexique", "Mexico_Neotropical", "Mexico_Nearctic")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")]  <- FALSE

# Myanmar. 22.00, 98.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Birmania", "Birmanie", "Burmah")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Birmania", "Birmanie", "Burmah")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Locality_code %in% c("Birmania", "Birmanie", "Burmah")) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")]  <- FALSE

# Nigeria. 10.00, 8.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Nigeria") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Nigeria") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Nigeria") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Sierra Leone. 8.50, -11.50
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sierra Leone") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sierra Leone") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sierra Leone") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Solomon Islands. -8.00, 159.00
Biogeographic_database_Ponerinae_flags$adm1[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Solomon Islands") & replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_initial, pattern = "Isabel"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- "Isabel Province"
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Solomon Islands") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Solomon Islands") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Solomon Islands") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Sri Lanka. 7.00, 81.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sri Lanka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sri Lanka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Sri Lanka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Suriname. 4.00, -56.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Suriname") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Suriname") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Suriname") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Taiwan. 23.50, 121.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Taiwan") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Formosa_TW") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Taiwan") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Formosa_TW") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Taiwan") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Formosa_TW") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# United Kingdom. 54.00, -2.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "United Kingdom") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "United Kingdom") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "United Kingdom") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE

# Zimbabwe. -19.00, 29.00
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Zimbabwe") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Rhodesia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Zimbabwe") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Rhodesia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Zimbabwe") & (Biogeographic_database_Ponerinae_flags$Locality_code == "Rhodesia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- FALSE


## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")

### 2.6.4.4/ Update geometry ####

# Convert Flagged df to sf 
Biogeographic_database_Ponerinae_flags$Latitude_copy <- Biogeographic_database_Ponerinae_flags$Latitude_dec
Biogeographic_database_Ponerinae_flags$Longitude_copy <- Biogeographic_database_Ponerinae_flags$Longitude_dec
Biogeographic_database_Ponerinae_flags <- st_as_sf(Biogeographic_database_Ponerinae_flags, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_flags) <- sp::CRS('+init=EPSG:4326')
# Remove copies of Long/Lat of still present
Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags %>% select(-c("Longitude_copy", "Latitude_copy"))
Biogeographic_database_Ponerinae_flags

### 2.6.4.5/ Cut-Paste NA entries in the NA dataset ####

# Extract entries with new NA coords (assign to NA because they were flagged as error)
Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun <- cc_val(x = Biogeographic_database_Ponerinae_flags,
                                                                         lon = "Longitude_dec", lat = "Latitude_dec", 
                                                                         value = "flagged", verbose = T)

table(Biogeographic_database_Ponerinae_flags$valid_coordinates)
table(Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun)
View(x = Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun, ])
View(x = Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun & Biogeographic_database_Ponerinae_flags$valid_coordinates, ])

# Extract valid coordinates to rerun tests
Biogeographic_database_Ponerinae_flags_rerun <- Biogeographic_database_Ponerinae_flags %>% 
  filter(Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun) %>% 
  select(-c(".val", ".cap", ".cen", ".sea", ".con", ".otl", ".inst", ".dpl", ".summary")) %>% 
  st_drop_geometry() # %>% 
# select(-c("Latitude_copy", "Longitude_copy"))
class(Biogeographic_database_Ponerinae_flags_rerun) <- "data.frame"

# Bind new NA entries to the NA dataset
Biogeographic_database_Ponerinae_NA_coords_new <- st_drop_geometry(Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun, ])
Biogeographic_database_Ponerinae_NA_coords <- rbind(Biogeographic_database_Ponerinae_NA_coords, Biogeographic_database_Ponerinae_NA_coords_new[, names(Biogeographic_database_Ponerinae_NA_coords)])

## Save sf dataframe with coordinates with NA
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

### 2.6.4.6/ Rerun evaluation to update flags ####

Biogeographic_database_Ponerinae_flags_rerun <- clean_coordinates(x = Biogeographic_database_Ponerinae_flags_rerun,
                                                            lon = "Longitude_dec", lat = "Latitude_dec",
                                                            countries = "Country_ISO3_code",
                                                            species = "Current_name",
                                                            value = "spatialvalid",         # To define the type of output. Can be a spatial df ("spatialvalid"), a logical vector recording entries with at least one flag ("flagged"), a df with dubious coordinates removed ("clean").
                                                            verbose = T,
                                                            # List of tests to run
                                                            tests = c(#"aohi",          # Coordinates found in Artificial Hotspots
                                                              "capitals",      # Coordinates found in capitals
                                                              "centroids",     # Coordinates found at centroids of countries
                                                              "countries",     # Test if occurrences are found within a list of countries or a custom SpatialPolygon
                                                              "duplicates",    # Check for duplicates (identical coordinates) within taxa
                                                              #"equal",         # Check for coordinates with latitude = longitude
                                                              # "gbif",         # Check if the coordinates is within a 1km radius around the GBIF headquarter
                                                              "institutions",  # Test if occurrences are found arouf known institution locations
                                                              "outliers",      # Test for minimum/mean distance to other occurrences per taxa
                                                              # "range"         # Test if coordinates fall within species range extracted from IUCN
                                                              "seas"           # Test if fall in the ocean
                                                              # "urban"          # Test if fall within urban areas
                                                              # "zeros"          # Test for coordinates with null latitude or longitude
                                                            ), 
                                                            outliers_method = "distance",  # Use minimum distance to detect outlier
                                                            outliers_td = 1000,            # Minimum distance to detect outlier (in km)
                                                            # outliers_mtp = 5.             # Multiplier for the IRQ of mean distance to flag as outlier
                                                            outliers_size = 3,             # Minimum number of records for a given taxa to test for outliers
                                                            # seas_scale = 10,               # Adjust resolution of landmass polygons used as reference. 10 = highest resolution.
                                                            seas_ref = map_NE_lands_sf_buffer,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                                                            # seas_ref = map_NE_all_lands_sf,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                                                            inst_rad = 200,                # Radius around Biodiversity institutions (in meters)
                                                            capitals_rad = 3000,           # Radius around capitals (in meters)
                                                            centroids_rad = 200,           # Radius around country centroids (in meters)
                                                            country_ref = Countries_NE_sf,  # sf SpatialPolygons to use as reference for countries
                                                            country_refcol = "ISO_A3",     # Name of ISO3 column
                                                            country_buffer = 5000,         # Buffer around country polygons (in meters)
                                                            # seas_buffer = 5000,            # Buffer around landmass limits (in meters) # Does not work properly with scale 10 map. Make our own prior
                                                            aohi_rad = 3000)               # Buffer around AOHI coordinates (in meters)

# In case of successful rerun, replace previous run
# Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags_rerun

# Check results for each test in specific columns in the output
View(Biogeographic_database_Ponerinae_flags)


## 2.6.4.7/ Readjust flag types


# Convert Flagged df to sf 
Biogeographic_database_Ponerinae_flags <- st_drop_geometry(Biogeographic_database_Ponerinae_flags)
Biogeographic_database_Ponerinae_flags$Latitude_copy <- Biogeographic_database_Ponerinae_flags$Latitude_dec
Biogeographic_database_Ponerinae_flags$Longitude_copy <- Biogeographic_database_Ponerinae_flags$Longitude_dec
Biogeographic_database_Ponerinae_flags <- st_as_sf(Biogeographic_database_Ponerinae_flags, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_flags) <- sp::CRS('+init=EPSG:4326')


## Rerun all corrections so false positives are fixed


# Design factor level for type of flags
# Hierarchy of priorities: "zeros", "equal", "countries", "seas", "centroids", "capitals", "aohi", "institutions", "outliers", "urban", "duplicates"
Biogeographic_database_Ponerinae_flags$flag_type <- "OK"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.urb] <- "urban"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst] <- "institutions"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.aohi] <- "aohi"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap] <- "capitals"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cen] <- "centroids"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.sea] <- "seas"
Biogeographic_database_Ponerinae_flags$flag_type[!as.logical(Biogeographic_database_Ponerinae_flags$.con)] <- "countries"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.equ] <- "equal"
# Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.zer] <- "zeros"

Biogeographic_database_Ponerinae_flags$flag_type <- as.factor(Biogeographic_database_Ponerinae_flags$flag_type)  

table(Biogeographic_database_Ponerinae_flags$flag_type)

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "countries", ])
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "seas", ])
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "centroids", ])



## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


### 2.6.5/ Curate records in capitals ####

# Look at locality names to check if nothing suspicious that could reflect the location of a specimen in collection rather than from the field.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "capitals", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# All seems fine, so mark them as false positive
# Biogeographic_database_Ponerinae_flags$.cap <- TRUE

# Update flag types
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & !Biogeographic_database_Ponerinae_flags$.inst] <- "institutions"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & Biogeographic_database_Ponerinae_flags$.inst & !Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & Biogeographic_database_Ponerinae_flags$.inst & Biogeographic_database_Ponerinae_flags$.otl] <- "OK"
Biogeographic_database_Ponerinae_flags$flag_type <- as.factor(Biogeographic_database_Ponerinae_flags$flag_type)  

# Mark them as false positive
Biogeographic_database_Ponerinae_flags$.cap <- TRUE

table(Biogeographic_database_Ponerinae_flags$flag_type)

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


### 2.6.6/ Curate records in around institutions ####

# Look at locality names to check if nothing suspicious that could reflect the location of a specimen in collection rather than from the field.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "institutions", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# All seems fine, so mark them as false positive
# Biogeographic_database_Ponerinae_flags$.inst <- TRUE

# Update flag types
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst & !Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst & Biogeographic_database_Ponerinae_flags$.otl] <- "OK"
Biogeographic_database_Ponerinae_flags$flag_type <- as.factor(Biogeographic_database_Ponerinae_flags$flag_type)

# All seems fine, so mark them as false positive
Biogeographic_database_Ponerinae_flags$.inst <- TRUE

table(Biogeographic_database_Ponerinae_flags$flag_type)

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


### 2.6.7/ Curate outliers based on distance ####

## 2.6.7.1/ Plot occurrences of one taxa with distance outliers ####

# Look at locality names to check if outliers are credible or not.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "outliers", ] %>% arrange(Current_name, Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

taxa_with_outliers_dist <- unique(Biogeographic_database_Ponerinae_flags$Current_name[Biogeographic_database_Ponerinae_flags$flag_type == "outliers"])
taxa_with_outliers_dist <- taxa_with_outliers_dist[order(taxa_with_outliers_dist)]
length(taxa_with_outliers_dist) # 205 taxa with distance-based outliers !

# Create Green/Red palette
Green_col <- RColorBrewer::brewer.pal(n = 6, "Greens")[5]
Red_col <- RColorBrewer::brewer.pal(n = 6, "Reds")[5]
pal_GR <- c(Green_col, Red_col)

# Taxa to test
Focal_taxon <- "Anochetus_afr01"
Focal_taxon <- "Anochetus_afr08"

mapviewOptions(fgb = FALSE) # Need this option to enable screenshots! Need to be set before creating the mapview object

# Plot by taxa in case of doubt for the group to evaluate
interactive_map_single_taxon <- mapview::mapview(x = Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Current_name == Focal_taxon, ], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                                 zcol = "flag_type",
                                                 cex = 7,
                                                 burst = T, legend = T,
                                                 # col.regions = pal(nb_flags), # Color to use to plot spatial units in the spatial object
                                                 col.regions = pal_GR, # Color to use to plot spatial units in the spatial object
                                                 layer.name = "Occurrences", # To specify the legend name of this layer
                                                 # map.types = mapviewGetOption("basemaps")
                                                 map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                                 # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                                 # ... # and many more, depending of the class of the spatial object
)

# Plot without title
interactive_map_single_taxon

mapshot(x = interactive_map_single_taxon,
        url = paste0("./maps/Taxa_occurrence_Control_maps/Interactive_map_occurrences_control_",Focal_taxon,".html"), # To save an interactive .html map
        remove_controls = NULL) # c("zoomControl", "layersControl", "homeButton", "scaleBar"))  # To specify which control buttons to remove (or not) from the output

### 2.6.7.2/ Loop to plot maps across taxa with distance-based outliers ####

for (i in seq_along(taxa_with_outliers_dist))
{
  Focal_taxon <- taxa_with_outliers_dist[i]
  
  # Need this option to enable screenshots! Need to be set before creating the mapview object
  mapviewOptions(fgb = FALSE) 
  
  # Plot by taxa in case of doubt for the group to evaluate
  interactive_map_focal_taxon <- mapview::mapview(x = Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Current_name == Focal_taxon, ], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                                   zcol = "flag_type",
                                                   cex = 7,
                                                   burst = T, legend = T,
                                                   # col.regions = pal(nb_flags), # Color to use to plot spatial units in the spatial object
                                                   col.regions = pal_GR, # Color to use to plot spatial units in the spatial object
                                                   layer.name = "Occurrences", # To specify the legend name of this layer
                                                   # map.types = mapviewGetOption("basemaps")
                                                   map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                                   # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                                   # ... # and many more, depending of the class of the spatial object
  )
  
  # Create title HTML object
  tag.map.title <- htmltools::tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 22px;
  }
   "))
  
  title <- tags$div(
    tag.map.title, 
    HTML(str_replace(string = Focal_taxon, pattern = "_", replacement = " ")) # Type title
  )  
  
  # Add title to the leafmap /mapview object
  interactive_map_focal_taxon@map <- interactive_map_focal_taxon@map %>%
    leaflet::addControl(title, position = "bottomleft", className = "map-title")
  
  interactive_map_focal_taxon
  
  # To save a screenshot
  mapshot(x = interactive_map_focal_taxon,
          file = paste0("./maps/Taxa_occurrence_Control_maps/Outliers_distance/occurrence_control_map_",Focal_taxon,".jpeg"), 
          remove_controls = NULL # delay = 0.5
  )
  
  mapviewOptions(fgb = TRUE)
  
  # Print progress
  cat(paste0(Sys.time(), " - Map plotted for ", Focal_taxon, " = Taxon N°",i,"/",length(taxa_with_outliers_dist),"\n"))
}


### 2.7/ Curate outliers found in dubious bioregions ####

## 2.7.1/ Assign bioregion to each entry ####

# Load country_sf with metadata
Countries_NE_sf <- readRDS(file = "./input_data/NaturalEarth_maps/Countries_NE_sf.rds")

# Export metadata to be filled manually for Bioregions
Countries_NE_sf_metadata <- st_drop_geometry(Countries_NE_sf[, c("SUBUNIT", "ISO_A3")])
Countries_NE_sf_metadata$Bioregion <- NA
# write.xlsx(x = Countries_NE_sf_metadata, file = "./input_data/Biogeographic_data/Countries_NE_sf_metadata.xlsx")

# Load metadata with bioregion assignment
Countries_NE_sf_metadata <- openxlsx::read.xlsx(xlsxFile = "./input_data/Biogeographic_data/Countries_NE_sf_metadata.xlsx")

# Assign Bioregion to entries
Biogeographic_database_Ponerinae_flags$Bioregion <- Countries_NE_sf_metadata$Bioregion[match(Biogeographic_database_Ponerinae_flags$Country_ISO3_code, Countries_NE_sf_metadata$ISO_A3)]

# Deal with the particular cases of countries overlapping bioregions

# France: Mayotte, Reunion, Guyane, Martinique, Guadeloupe
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France", ] %>% arrange(Country_initial))
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Palearctic"
Biogeographic_database_Ponerinae_flags$Bioregion[str_detect(string = Biogeographic_database_Ponerinae_flags$Country_initial, pattern =  "Guadeloupe") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial ==  "Martinique", replace = F) | replace_na(data = Biogeographic_database_Ponerinae_flags$adm1 ==  "Martinique", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial ==  "Mayotte", replace = F) | replace_na(data = Biogeographic_database_Ponerinae_flags$adm1 ==  "Mayotte", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial ==  "Reunion", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial ==  "Reunion", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial ==  "French Guiana", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"

# Indonesia: Longitude 125.5°. Not right but good enough for rough first estimates
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Longitude_dec < 125.5) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Indonesia")] <- "Indomalaya"
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Longitude_dec >= 125.5) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Indonesia")] <- "Australasia"

# Mexico: 22° Latitude. Not right but good enough for rough first estimates
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Latitude_dec < 22) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Mexico")] <- "Neotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Latitude_dec >= 22) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Mexico")] <- "Nearctic"

# China: 33° Latitude (not right but enough for our data)
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Latitude_dec < 33) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "China")] <- "Indomalaya"
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Latitude_dec >= 33) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "China")] <- "Palearctic"

# USA, Hawaii: 150° longitude (works because we have no points in Alaska
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Longitude_dec < -150) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "United States")] <- "Australasia"
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Longitude_dec >= -150) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "United States")] <- "Nearctic"

# TAAF, Ile Europa: Afrotropics
Biogeographic_database_Ponerinae_flags$Bioregion[replace_na(data = Biogeographic_database_Ponerinae_flags$Locality == "Europa Island", replace = F)] <- "Afrotropics"

# Netherlands, Saba Island: Neotropics
Biogeographic_database_Ponerinae_flags$Bioregion[replace_na(data = Biogeographic_database_Ponerinae_flags$Country_initial == "Bonaire, Sint Eustatius and Saba", replace = F)] <- "Neotropics"

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


# Plot to check bioregion assignments
nb_bioregions <- length(unique(Biogeographic_database_Ponerinae_flags$Bioregion))
pal_bioregions <- brewer.pal(n = nb_bioregions, name = "Spectral")

interactive_Bioregion_flags <- mapview::mapview(x = Biogeographic_database_Ponerinae_flags[ ,], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                          zcol = "Bioregion",
                                          cex = 5,
                                          burst = T, legend = T,
                                          # col.regions = pal(nb_flags), # Color to use to plot spatial units in the spatial object
                                          col.regions = pal_bioregions, # Color to use to plot spatial units in the spatial object
                                          layer.name = "Bioregions", # To specify the legend name of this layer
                                          # map.types = mapviewGetOption("basemaps")
                                          map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                          # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                          # ... # and many more, depending of the class of the spatial object
)

interactive_Bioregion_flags


## 2.7.2/ Create taxa x bioregions quantitative table ####

# Create table with nb of occurrence per bioregions
Taxa_bioregions_quantitative_table <- Biogeographic_database_Ponerinae_flags %>%
  st_drop_geometry() %>% 
  group_by(Current_name, Bioregion) %>% 
  summarise(nb_occ = n()) %>% 
  pivot_wider(names_from = Bioregion, values_from = nb_occ) %>%
  mutate_all(., ~replace_na(.,0))

# Turn into percentage, then binary test for bioregions hosting less than 20% of occurrences
Taxa_bioregions_quantitative_table$Total_occ <- apply(X = Taxa_bioregions_quantitative_table[, -1], MARGIN = 1, FUN = sum) 
Taxa_bioregions_perc_table <- round(Taxa_bioregions_quantitative_table[, -1]/Taxa_bioregions_quantitative_table$Total_occ * 100, 1)
Taxa_bioregions_binary_test_table <- t(apply(X = Taxa_bioregions_perc_table, MARGIN = 1, FUN = function (x) {(x > 0) & (x <= 20)} ))
Taxa_bioregions_quantitative_table$Outliers_bioregion  <- apply(X = Taxa_bioregions_binary_test_table, MARGIN = 1, FUN = any)

table(Taxa_bioregions_quantitative_table$Outliers_bioregion)

# List taxa with outliers in dubious bioregions
taxa_with_bioregion_outliers <-Taxa_bioregions_quantitative_table$Current_name[Taxa_bioregions_quantitative_table$Outliers_bioregion]
taxa_with_bioregion_outliers <- taxa_with_bioregion_outliers[order(taxa_with_bioregion_outliers)]

### 2.7.3/ Plot distribution maps of suspicious taxa for bioregion-based outliers

# Prepare bioregion palette
nb_bioregions <- length(unique(Biogeographic_database_Ponerinae_flags$Bioregion))
pal_bioregions <- brewer.pal(n = nb_bioregions, name = "Spectral")
Bioregions_ordered <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Palearctic")

## Loop across taxa to plot maps

for (i in seq_along(taxa_with_bioregion_outliers))
{
  Focal_taxon <- taxa_with_bioregion_outliers[i]
  
  # Extract taxon data
  Taxon_data <- Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Current_name == Focal_taxon, ]
  
  # Define bioregion palette according to available bioregions
  Bioregions_available <- unique(Taxon_data$Bioregion)
  Bioregions_indices <- which(Bioregions_ordered %in% Bioregions_available)
  pal_taxon <- pal_bioregions[Bioregions_indices]
  
  # Need this option to enable screenshots! Need to be set before creating the mapview object
  mapviewOptions(fgb = FALSE) 
  
  # Plot by taxa in case of doubt for the group to evaluate
  interactive_map_focal_taxon <- mapview::mapview(x = Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Current_name == Focal_taxon, ], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                                  zcol = "Bioregion",
                                                  cex = 7,
                                                  burst = T, legend = T,
                                                  # col.regions = pal(nb_flags), # Color to use to plot spatial units in the spatial object
                                                  col.regions = pal_taxon, # Color to use to plot spatial units in the spatial object
                                                  layer.name = "Bioregions", # To specify the legend name of this layer
                                                  # map.types = mapviewGetOption("basemaps")
                                                  map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                                  # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                                  # ... # and many more, depending of the class of the spatial object
  )
  
  # Create title HTML object
  tag.map.title <- htmltools::tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 22px;
  }
   "))
  
  title <- tags$div(
    tag.map.title, 
    HTML(str_replace(string = Focal_taxon, pattern = "_", replacement = " ")) # Type title
  )  
  
  # Add title to the leafmap /mapview object
  interactive_map_focal_taxon@map <- interactive_map_focal_taxon@map %>%
    leaflet::addControl(title, position = "bottomleft", className = "map-title")
  
  interactive_map_focal_taxon
  
  # To save a screenshot
  mapshot(x = interactive_map_focal_taxon,
          file = paste0("./maps/Taxa_occurrence_Control_maps/Outliers_bioregions/occurrence_control_map_",Focal_taxon,".jpeg"), 
          remove_controls = NULL # delay = 0.5
  )
  
  mapviewOptions(fgb = TRUE)
  
  # Print progress
  cat(paste0(Sys.time(), " - Map plotted for ", Focal_taxon, " = Taxon N°",i,"/",length(taxa_with_bioregion_outliers),"\n"))
}

### 2.7.3/ Create .xlsx tables to evaluate taxa ####

taxa_with_outliers <- union(taxa_with_outliers_dist, taxa_with_bioregion_outliers)

Evaluation_table_taxa_outliers <- data.frame(Taxa = taxa_with_outliers, Distance_outlier = NA, Bioregion_outlier = NA, Evaluation = NA, Interactive_map_needed = NA, Comments = NA)
Evaluation_table_taxa_outliers$Distance_outlier <- (Evaluation_table_taxa_outliers$Taxa %in% taxa_with_outliers_dist)
Evaluation_table_taxa_outliers$Bioregion_outlier <- (Evaluation_table_taxa_outliers$Taxa %in% taxa_with_bioregion_outliers)

openxlsx::write.xlsx(x = Evaluation_table_taxa_outliers, file = "./input_data/Biogeographic_data/Evaluation_table_taxa_outliers.xlsx")


##### 3/ Deal with NA coordinates #####

## Load occurrences with no GPS coordinates
Biogeographic_database_Ponerinae_NA_coords <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

# Fix coordinates to NA
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec <- NA
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec <- NA

# Indicate nature of occurrences
Biogeographic_database_Ponerinae_NA_coords$Interpolated_coordinates <- TRUE

### 3.1/ Load shape files ####

# Load bentity2 polygon shapes
Bentity2_sf <- st_read(dsn = "./input_data/GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/", layer = "Bentity2_shapefile_fullres")

# Correct mistakes
Bentity2_sf$BENTITY2_N[Bentity2_sf$BENTITY2_N == "Norde de Santander"] <- "Norte de Santander"

# Remove polygons with NA names
table(is.na(Bentity2_sf$BENTITY2_N))


## Load geoBoundaries shapes files for Countries, adm1 and adm2

geoBoundaries_Countries_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM0/", layer = "geoBoundariesCGAZ_ADM0")
# plot(geoBoundaries_Countries_sf["shapeName"])
table(geoBoundaries_Countries_sf$shapeType)

geoBoundaries_adm1_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM1/", layer = "geoBoundariesCGAZ_ADM1")
# plot(geoBoundaries_adm1_sf["shapeName"])
# plot(geoBoundaries_adm1_sf["shapeGroup"])
table(geoBoundaries_adm1_sf$shapeType)

geoBoundaries_adm2_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM2/", layer = "geoBoundariesCGAZ_ADM2")
# plot(geoBoundaries_adm2_sf["shapeName"])
# plot(geoBoundaries_adm2_sf["shapeGroup"])
# plot(geoBoundaries_adm2_sf["shapeType"])
table(geoBoundaries_adm2_sf$shapeType)

# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$shapeGroup == "AUS", "shapeName"])

# Rename data columns
names(geoBoundaries_Countries_sf) <- c("Country_shapeGroup", "Country_shapeType",  "Country_shapeName",  "geometry")
names(geoBoundaries_adm1_sf) <- c("adm1_shapeName", "adm1_shapeID", "adm1_shapeGroup", "adm1_shapeType", "geometry")
names(geoBoundaries_adm2_sf) <- c("adm2_shapeName", "adm2_shapeID", "adm2_shapeGroup", "adm2_shapeType", "geometry")

# Remove polygons with NA names
table(is.na(geoBoundaries_Countries_sf$Country_shapeName))
table(is.na(geoBoundaries_adm1_sf$adm1_shapeName))
table(is.na(geoBoundaries_adm2_sf$adm2_shapeName))

geoBoundaries_adm2_sf <- geoBoundaries_adm2_sf %>% 
  filter(!is.na(geoBoundaries_adm2_sf$adm2_shapeName))

## Load Natural Earth shape files

NE_Countries_sf <- readRDS(file = "./input_data/NaturalEarth_maps/Countries_NE_sf.rds")
NE_subunits_sf <- ne_load(scale = 50, category = "cultural", type = 'map_subunits', returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_50m_admin_0_map_subunits/")
NE_adm1_sf <- ne_load(file_name = "ne_10m_admin_1_states_provinces", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_10m_admin_1_states_provinces/")

# plot(NE_Countries_sf[, "ISO_A3"])
# plot(NE_subunits_sf[, "SUBUNIT"]) # If possible, use subunits because some countries have disjointed shapes such as France with the DOMTOM
# plot(NE_adm1_sf[, "name"])

# Remove polygons with NA names
table(is.na(NE_Countries_sf$ADMIN))
table(is.na(NE_subunits_sf$SUBUNIT))
table(is.na(NE_adm1_sf$name))

NE_adm1_sf <- NE_adm1_sf %>% 
  filter(!is.na(NE_adm1_sf$name))


### Load Geographic Names Server data

GNS_Localities_df <- read.csv(file = "./input_data/GNS_Administrative_Regions/Administrative_Regions.txt", sep = "\t")
GNS_Localities_sf <- st_as_sf(x = GNS_Localities_df, coords = c("long_dd", "lat_dd"))


### 3.2/ Compute areas of shape files to use to estimate uncertainty when drawing random points ####

## Areas of Bentity2_sf
sf_use_s2(FALSE)
Bentity2_sf$area <- st_area(Bentity2_sf)
Bentity2_sf$area <- set_units(Bentity2_sf$area, km^2)
sf_use_s2(TRUE)

# plot(x = Bentity2_sf[ , "BENTITY2_N"])
# plot(x = Bentity2_sf[ , "area"])

# Save Bentity2_sf
saveRDS(object = Bentity2_sf, file = "./input_data/GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/Bentity2_sf.rds")

# Load Bentity2_sf
Bentity2_sf <- readRDS(file = "./input_data/GABI_Data_Release1.0_18012020/Bentity2_shapefile_fullres/Bentity2_sf.rds")


## Areas of geoBoundaries shape files
sf_use_s2(FALSE)
geoBoundaries_Countries_sf$area <- st_area(geoBoundaries_Countries_sf)
geoBoundaries_Countries_sf$area <- set_units(geoBoundaries_Countries_sf$area, km^2)
geoBoundaries_adm1_sf$area <- st_area(geoBoundaries_adm1_sf)
geoBoundaries_adm1_sf$area <- set_units(geoBoundaries_adm1_sf$area, km^2)
geoBoundaries_adm2_sf$area <- st_area(geoBoundaries_adm2_sf)
geoBoundaries_adm2_sf$area <- set_units(geoBoundaries_adm2_sf$area, km^2)
sf_use_s2(TRUE)

# plot(x = geoBoundaries_Countries_sf[ , "area"])
# plot(x = geoBoundaries_adm1_sf[ , "area"])
# plot(x = geoBoundaries_adm2_sf[ , "area"])

# Save GeoBoundaries shape files
saveRDS(object = geoBoundaries_Countries_sf, file = "./input_data/geoBoundaries/geoBoundaries_Countries_sf.rds")
saveRDS(object = geoBoundaries_adm1_sf, file = "./input_data/geoBoundaries/geoBoundaries_adm1_sf.rds")
saveRDS(object = geoBoundaries_adm2_sf, file = "./input_data/geoBoundaries/geoBoundaries_adm2_sf.rds")

# Load GeoBoundaries shape files
geoBoundaries_Countries_sf <- readRDS(file = "./input_data/geoBoundaries/geoBoundaries_Countries_sf.rds")
geoBoundaries_adm1_sf <- readRDS(file = "./input_data/geoBoundaries/geoBoundaries_adm1_sf.rds")
geoBoundaries_adm2_sf <- readRDS(file = "./input_data/geoBoundaries/geoBoundaries_adm2_sf.rds")


## Areas of Natural Earth shape files
sf_use_s2(FALSE)
NE_Countries_sf$area <- st_area(NE_Countries_sf)
NE_Countries_sf$area <- set_units(NE_Countries_sf$area, km^2)
NE_subunits_sf$area <- st_area(NE_subunits_sf)
NE_subunits_sf$area <- set_units(NE_subunits_sf$area, km^2)
NE_adm1_sf$area <- st_area(NE_adm1_sf)
NE_adm1_sf$area <- set_units(NE_adm1_sf$area, km^2)
sf_use_s2(TRUE)

# plot(x = NE_Countries_sf[ , "area"])
# plot(x = NE_subunits_sf[ , "area"])
# plot(x = NE_adm1_sf[ , "area"])

# Save Natural Earth shape files
saveRDS(object = NE_Countries_sf, file = "./input_data/NaturalEarth_maps/NE_Countries_sf.rds")
saveRDS(object = NE_subunits_sf, file = "./input_data/NaturalEarth_maps/NE_subunits_sf.rds")
saveRDS(object = NE_adm1_sf, file = "./input_data/NaturalEarth_maps/NE_adm1_sf.rds")

# Load Natural Earth shape files
NE_Countries_sf <- readRDS(file = "./input_data/NaturalEarth_maps/NE_Countries_sf.rds")
NE_subunits_sf <- readRDS(file = "./input_data/NaturalEarth_maps/NE_subunits_sf.rds")
NE_adm1_sf <- readRDS(file = "./input_data/NaturalEarth_maps/NE_adm1_sf.rds")


### 3.3/ Find matching Localities/adm2/adm1/Bentity2/Country_ISO2_name based on names ####

## 3.3.1/ Correct some entries to find better matches ####

## Create new column for SUBUNIT names

Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name <- Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name

# Hong Kong is a Country, not an adm1 of China
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Hong Kong"] <- "Hong Kong S.A.R." # SUBUNIT in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Hong Kong"] <- "HKG"
# French Guiana is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "French Guiana"] <- "French Guiana" # SUBUNIT in NE_SUBUNITS
# Mayotte is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Mayotte"] <- "Mayotte" # SUBUNIT in NE_SUBUNITS
# Guadeloupe is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Guadeloupe"] <- "Guadeloupe" # SUBUNIT in NE_SUBUNITS
# Martinique is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Martinique"] <- "Martinique" # SUBUNIT in NE_SUBUNITS
# Christmas Island is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Christmas Island"] <- "Christmas Island" # SUBUNIT in NE_SUBUNITS
# Tokelau is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_initial == "Tokelau"] <- "Tokelau" # SUBUNIT in NE_SUBUNITS
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Tokelau"] <- "Tokelau" # SUBUNIT in NE_SUBUNITS
# Pointe-à-Pitre is in Guadeloupe (Subunit) 
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Pointe-à-Pitre", replace = F)] <- "Guadeloupe" # SUBUNIT in NE_SUBUNITS

# Puerto Rico is a Subunit in NE
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Puerto Rico"] <- "Puerto Rico" # SUBUNIT in NE_SUBUNITS

# SUBUNIT New Caledonia has bad polygons => Match on Country Admin Instead
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name == "New Caledonia", replace = F)] <- NA 


## Fix errors with Country names

# Correct error with Northern Mariana Islands
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "N. Mariana Is."] <- "Northern Mariana Islands"
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Northern Mariana Islands", replace = F)] <- "Northern Mariana Islands"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Northern Mariana Islands", replace = F)] <- "MNP"

# Correct errors with Panama
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Panam·"] <- "Panama" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Panama"] <- "PAN" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Panama"] <- "Panama" # SUBUNIT in NE_Subunits

# Occurrences in Philippines (bentity2) are labeled as in Indonesia
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Philippines"] <- "Philippines" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Philippines"] <- "PHL" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Philippines"] <- "Philippines" # SUBUNIT in NE_Subunits

# Occurrences in Torres Strait Islands (bentity2) are labeled as in PNG as they should be in Australia
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Torres Strait Islands"] <- "Australia" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Torres Strait Islands"] <- "AUS" # Country in NE_Countries
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Torres Strait Islands"] <- "Australia" # SUBUNIT in NE_Subunits


## Fix errors with Bentity2

# Remove Bentity names for cases where SUBUNIT and Country admin are a better match

# Lesser Antilles in Bentity2 => match on SUBUNIT instead because it encompasses many countries
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Lesser Antilles", replace = F)] <- NA 
# Azerbaijan/Turkmenistan/Iran in Bentity2 => match on SUBUNIT instead because Bentity2 allows occurrences to be drawn in the Caspian Sea
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Azerbaijan", replace = F)] <- NA 
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Turkmenistan", replace = F)] <- NA 
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Iran", replace = F)] <- NA 
# Wallis and Futuna in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Wallis and Futuna", replace = F)] <- NA 
# Fiji in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Fiji", replace = F)] <- NA 
# Tonga in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Tonga", replace = F)] <- NA 
# Tokelau in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Tokelau", replace = F)] <- NA 
# Cuba in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Cuba", replace = F)] <- NA 
# Samoan Islands in Bentity2 => match on SUBUNIT instead because does not make a difference between American Samoa and Samoa
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Samoan Islands", replace = F)] <- NA 
# Bahamas in Bentity2 => match on SUBUNIT instead because some records from Barbados are attributed to Bahamas...
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Bahamas", replace = F)] <- NA 
# Hispaniola in Bentity2 => match on SUBUNIT instead because does not make a difference between Dominican Republic and Haiti
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Hispaniola", replace = F)] <- NA 
# Serbia in Bentity2 => match on SUBUNIT instead because does not make a difference between Kosovo and Serbia
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Serbia", replace = F)] <- NA 
# Israel and Palestine in Bentity2 => match on SUBUNIT instead because does not make a difference between Israel and Palestine
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Israel and Palestine", replace = F)] <- NA 
# Cameroon line Archipelago in Bentity2 => match on SUBUNIT instead because does not make a difference between Equatorial Guinea and São Tomé and Principe
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Cameroon line Archipelago", replace = F)] <- NA
# Sudan in Bentity2 => match on SUBUNIT instead because does not make a difference between Sudan and South Sudan
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Sudan", replace = F)] <- NA
# Occurrences in Laffarugh should be in Somaliland, not Ethiopia => Remove Bentity2 and use SUBUNIT Somaliland to match
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Laffarugh", replace = F)] <- NA
# Somalia in Bentity2 => match on SUBUNIT instead because does not make a difference between Somalia and Somaliland
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Somalia", replace = F)] <- NA
# Saudi Arabia in Bentity2 => match on SUBUNIT instead because it has bad polygon boundaries with Yemen
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Saudi Arabia", replace = F)] <- NA
# Comoros in Bentity2 => match on SUBUNIT instead because it does not make a difference between Comoros and Mayotte
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Comoros", replace = F)] <- NA
# Mascarene Islands in Bentity2 => match on SUBUNIT instead because it does not make a difference between Mauritius and La Reunion
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Mascarene Islands", replace = F)] <- NA
# Malaysia and Singapore => match on SUBUNIT instead because it does not make a difference between Malaysia and Singapoure
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Malaysia and Singapore", replace = F)] <- NA
# Caroline Islands in Bentity2 => Does not make a difference between Palau and Federated States of Micronesia. Match on SUBUNIT instead
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Caroline Islands", replace = F)] <- NA
# Marshall Islands in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Marshall Islands", replace = F)] <- NA
# Cyprus in Bentity2 => Does not make a difference between Cyprus and Northern Cyprus. Match on SUBUNIT instead
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Cyprus", replace = F)] <- NA
# Eritrea in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Eritrea", replace = F)] <- NA
# Seychelles in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Seychelles", replace = F)] <- NA
# Kiribati in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Kiribati", replace = F)] <- NA
# Panama in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Panama", replace = F)] <- NA
# New Caledonia in Bentity2 => match on SUBUNIT instead because bad polygons
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "New Caledonia", replace = F)] <- NA


# Hong Kong in Bentity2 => Rename SUBUNIT and Country as Hong Kong S.A.R. (HKG) instead of China
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Hong Kong", replace = F)] <- "Hong Kong S.A.R."
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Hong Kong", replace = F)] <- "Hong Kong S.A.R."
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Hong Kong", replace = F)] <- "HKG"

# Borneo in Bentity2 => match on adm1 instead because it does not make a difference between Malaysia and Indonesia
# Use Sarawak (adm1) for Malaysia when missing
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Borneo" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Malaysia" & is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- "Sarawak"
# Use East Kalimantan for Indonesia when missing
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Borneo" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- "East Kalimantan"
# Rename adequately when not missing by deleting the last portion
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Borneo" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- str_remove(string = Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Borneo" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)], pattern = " Province.*")

# Java in Bentity2 => Rename SUBUNIT and Country as Indonesia (IDN) instead of Malaysia
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Java", replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Java", replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Java", replace = F)] <- "IDN"

# Lesser Sunda Islands in Bentity2
# Name adm1 East Nusa Tenggara if missing
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Lesser Sunda Islands" & is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- "East Nusa Tenggara"
# Rename adequately when not missing by deleting the last portion
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Lesser Sunda Islands" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- str_remove(string = Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Lesser Sunda Islands" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)], pattern = " Province.*")

# New Guinea in Bentity2 => Does not make a difference between Papua New Guinea and Indonesia. 
# Match on SUBUNIT for PNG.
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "New Guinea" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Papua New Guinea", replace = F)] <- NA
# Rename SUBUNIT and Country as PNG for adm1 == "East New Britain Province"
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm1 == "East New Britain Province", replace = F)] <- "Papua New Guinea"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm1 == "East New Britain Province", replace = F)] <- "Papua New Guinea"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm1 == "East New Britain Province", replace = F)] <- "PNG"
# Name adm1 as Papua for Indonesia, except cases with adm1 == "West Papua region" (match on anything with "West") to rename as West Papua
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "New Guinea" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & str_detect(string = Biogeographic_database_Ponerinae_NA_coords$adm1, pattern = "West"), replace = F)] <- "West Papua"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "New Guinea" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & !str_detect(string = Biogeographic_database_Ponerinae_NA_coords$adm1, pattern = "West"), replace = F)] <- "Papua"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "New Guinea" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- "Papua"

# Solomon Islands in Bentity2 => Does not make a difference between Solomon Islands and Bougainville Islands (PNG). 
# For Solomon Islands, remove all 'Province' from adm1 to try to match true Solomon Islands occurrence. Otherwise, match on SUBUNIT. 
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Solomon Islands" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- str_remove(string = Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Solomon Islands" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)], pattern = " Province.*")
# For PNG occurrence, name adm1 as Autonomous Region of Bougainville
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Solomon Islands" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Papua New Guinea", replace = F)] <- "Autonomous Region of Bougainville"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Solomon Islands", replace = F)] <- NA

# Mariana Islands in Bentity2 => Does not make a difference between Northern Mariana Islands and Guam. 
# For Northern Mariana Islands, match on SUBUNIT
# For Guam also, but need to rename Locality == "Guam" as Guam (GUM)
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Guam", replace = F)] <- "Guam"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Guam", replace = F)] <- "Guam"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Guam", replace = F)] <- "GUM"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Mariana Islands", replace = F)] <- NA


## Fix errors with adm1

## GB_adm1

# Balearic Islands = Illes Balears
Biogeographic_database_Ponerinae_NA_coords$adm1[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Balearic Islands"] <- "Illes Balears" # adm1 in GB

# Tenerife & La Palma = Canarias
Biogeographic_database_Ponerinae_NA_coords$adm2[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Tenerife"] <- "Tenerife" 
Biogeographic_database_Ponerinae_NA_coords$adm1[Biogeographic_database_Ponerinae_NA_coords$adm1 == "Tenerife"] <- "Canarias" # adm1 in GB
Biogeographic_database_Ponerinae_NA_coords$adm2[Biogeographic_database_Ponerinae_NA_coords$adm1 == "La Palma"] <- "La Palma" 
Biogeographic_database_Ponerinae_NA_coords$adm1[Biogeographic_database_Ponerinae_NA_coords$adm1 == "La Palma"] <- "Canarias" # adm1 in GB

# Colón (BG_adm1) in Panama is wrongly assigned to Honduras => Rename Colón as Colón Province
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Colón" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Panama", replace = F)] <- "Colón Province"
# San Salvador (BG_adm1) in El Salvador is wrongly assigned to Bahamas => Rename San Salvador as Departamento de San Salvador
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "San Salvador" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "El Salvador", replace = F)] <- "Departamento de San Salvador"
# Santiago (BG_adm1) in Costa Rica is wrongly assigned to Dominican Republic => Add Paraiso as adm2 and rename Santiago as Cartago in adm1
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Santiago" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Costa Rica", replace = F)] <- "Paraiso"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Santiago" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Costa Rica", replace = F)] <- "Cartago"
# Saint Thomas (BG_adm1) in United States Virgin Islands is wrongly assigned to Barbados => Use SUBUNIT (not even Bentity2) as United States Virgin Islands do not have adm2 and adm1 but Saint Thomas is the proper name of the island
Biogeographic_database_Ponerinae_NA_coords$Locality[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Saint Thomas", replace = F)] <- "Saint Thomas"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Saint Thomas", replace = F)] <- NA
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Saint Thomas", replace = F)] <- NA
# Santander (BG_adm1) in Cantabria, Spain is wrongly assigned to Colombia => Rename Santander as Cantabria for both adm1 and adm2
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Santander" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Spain", replace = F)] <- "Cantabria"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Santander" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Spain", replace = F)] <- "Cantabria"
# Distrito Federal (BG_adm1) in Venezuela is wrongly assigned to Brazil => Rename Distrito Federal as Distrito Capital
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Distrito Federal" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Venezuela", replace = F)] <- "Distrito Capital"
# Misiones (BG_adm1) in Paraguay is wrongly assigned to Argentina => Rename Misiones as MISIONES
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Misiones" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Paraguay", replace = F)] <- "MISIONES"
# San José (BG_adm1) in Costa Rica is wrongly assigned to Uruguay => Rename San José as Provincia San José
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "San José" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Costa Rica", replace = F)] <- "Provincia San José"
# Luxembourg (BG_adm1) in Belgium is wrongly assigned to Luxembourg => Rename Luxembourg as Luxemburg
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Luxembourg" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Belgium", replace = F)] <- "Luxemburg"
# Sud (BG_adm1) in Cameroon is wrongly assigned to Italia => Rename Sud as South
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Sud" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Cameroon", replace = F)] <- "South"
# Sud-Ouest (BG_adm1) in Cameroon is wrongly assigned to Burkina-Faso => Rename Sud-Ouest as South-West
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Sud-Ouest" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Cameroon", replace = F)] <- "South-West"
# Centre (BG_adm1) in France is wrongly assigned to Cameroon => Rename Centre as Centre-Val de Loire
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Centre" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "France", replace = F)] <- "Centre-Val de Loire"
# Northern (BG_adm1) in PNG is wrongly assigned to Sudan => Rename Northern as Northern (Oro) Province
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Northern" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Papua New Guinea", replace = F)] <- "Northern (Oro) Province"
# Victoria (BG_adm1) in Tanzania is wrongly assigned to Australia => Rename Victoria as Kigoma + adm2 = Uvinza
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Victoria" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Tanzania", replace = F)] <- "Uvinza"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Victoria" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Tanzania", replace = F)] <- "Kigoma"

# Eastern Equatoria (BG_adm1) in South Sudan. Is properly assigned, but country do not match.
# (Re)name Torit District (adm2) as Torit + Change SUBUNIT and Country for South Sudan
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Eastern Equatoria", replace = F)] <- "South Sudan"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Eastern Equatoria", replace = F)] <- "South Sudan"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Eastern Equatoria", replace = F)] <- "SSD"
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Eastern Equatoria" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "South Sudan", replace = F)] <- "Torit"

# Western (BG_adm1) in Kenya is wrongly assigned to Zambia => Rename Western as Kakamega, except when Locality is "Bungoma, Nzoia", then adm1 is Bungoma
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Western" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Kenya" & Biogeographic_database_Ponerinae_NA_coords$Locality == "Bungoma, Nzoia", replace = F)] <- "Bungoma"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Western" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Kenya", replace = F)] <- "Kakamega"
# Central (BG_adm1) in Paraguay is wrongly assigned to Zambia => Rename Central as CENTRAL
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Central" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Paraguay", replace = F)] <- "CENTRAL"
# Central (BG_adm1) in PNG is wrongly assigned to Zambia => Rename Central as Central Province
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Central" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Papua New Guinea", replace = F)] <- "Central Province"
# Eastern (BG_adm1) in Kenya is wrongly assigned to Zambia => Rename Eastern as Isiolo (like the adm2)
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Eastern" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Kenya", replace = F)] <- "Isiolo"

# Bedforeshire (BG_adm1) in UK is wrongly assigned to another part of UK => Rename adm2 as Bedford. Rename adm1 as England
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Bedforeshire" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "United Kingdom", replace = F)] <- "Bedford"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Bedforeshire" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "United Kingdom", replace = F)] <- "England"

# Sikkim (BG_adm1) in India is wrongly assigned to China and Bhoutan because the adm1 shp has a mistake! => Use North District as adm2
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Sikkim" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "India", replace = F)] <- "North District"
# Sikkim (BG_adm1) in Bangladesh is wrongly assigned to India/China => Rename Sikkim as Rangpur, and use Rangpur as adm2 too
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Sikkim" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Bangladesh", replace = F)] <- "Rangpur"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Sikkim" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Bangladesh", replace = F)] <- "Rangpur"

# Okinawa Prefecture (BG_adm1) is wider than Yaeyama Islands (Bentity2) => Match using Bentity2 instead of adm1
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Okinawa Prefecture" & Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Yaeyama Islands", replace = F)] <- NA

# Western Province (BG_adm1) in Kenya is wrongly assigned to PNG => Rename Western Province as Kakamega 
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Western Province" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Kenya", replace = F)] <- "Kakamega"
# Western Province (BG_adm1) in Solomon Islands is wrongly assigned to PNG => Rename Western Province as Western + Name adm2 as Marovo
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Western" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Solomon Islands", replace = F)] <- "Marovo"

# South Korean adm1 => Remove " Province" from adm1 names
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "South Korea" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- str_remove(string = Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "South Korea" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)], pattern = " Province.*")
# Hiroshima Prefecture (BG_adm1) => Rename Hiroshima Prefecture as Hiroshima
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm1 == "Hiroshima Prefecture", replace = F)] <- "Hiroshima"

# Otsu (Locality) is in Kumamoto, Japan => Provide Kumamoto as adm1
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Otsu", replace = F)] <- "Kumamoto"

# Izu Islands are in Tokyo prefecture but better to match them on Bentity2, so remove adm1
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$bentity2_name == "Izu Islands", replace = F)] <- NA

# North Maluku Province [Maluku Utara] (BG_adm1) is North Maluku
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "North Maluku Province [Maluku Utara]", replace = F)] <- "North Maluku"

# Galapagos Islands (GB_adm1) should be renamed Galápagos
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Galapagos Islands", replace = F)] <- "Galápagos"
# Abidjan Department (GB_adm1) should be renamed District Autonome D'Abidjan
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Abidjan Department", replace = F)] <- "District Autonome D'Abidjan"

# Remove "Province" from Indonesia adm1 names
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm1 == "null", replace = F)] <- NA
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)] <- str_remove(string = Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Indonesia" & !is.na(Biogeographic_database_Ponerinae_NA_coords$adm1), replace = F)], pattern = " Province.*")

# Locality with "Kenting National Park" in Taiwan is in adm1 Pingtung County
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_NA_coords$Locality, pattern = "Kenting National Park"), replace = F)] <- "Pingtung County"
# Locality "Djebel Bou-Berak" is in Boumerdès (GB_adm1) 
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$Locality == "Djebel Bou-Berak", replace = F)] <- "Boumerdès"


## NE_adm1

# Veracruz (NE_adm1) in Argentina is wrongly assigned to Mexico => Rename Veracruz (adm1) as Córdoba
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Veracruz" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Argentina", replace = F)] <- "Córdoba"
# Bolivar (NE_adm1) in Colombia is wrongly assigned to Ecuador => Rename Bolivar as Bolívar
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Bolivar" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Colombia", replace = F)] <- "Bolívar"
# Bolivar (NE_adm1) in Venezuela is wrongly assigned to Ecuador => Rename Bolivar as Bolívar
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Bolivar" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Venezuela", replace = F)] <- "Bolívar"
# Macau (NE_adm1) need to be associated with Macao S.A.R (MAC) as a SUBUNIT/Country Admin
Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Macau", replace = F)] <- "Macao S.A.R"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Macau", replace = F)] <- "Macao S.A.R"
Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_code[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm1 == "Macau", replace = F)] <- "MAC"


## Fix errors with adm2

## GB_adm2
# San Bernardino (GB_adm2) in Cordillera (adm1) in Paraguay is wrongly assigned to California, US => Rename San Bernardino as SAN BERNARDINO
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "San Bernardino" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Paraguay", replace = F)] <- "SAN BERNARDINO"
# Furnas (GB_adm2) in Goiás (adm1) in Brazil is wrongly assigned to Nebraska, US => Delete Furnas and use adm1
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Furnas" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Brazil", replace = F)] <- NA
# San Felipe (GB_adm2) in Texas (adm1) in US is wrongly assigned to Mexico => Rename San Felipe as Austin
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "San Felipe" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "United States", replace = F)] <- "Austin"
# San Jeronimo (GB_adm2) in Baja Verapaz (adm1) in Guatemala is wrongly assigned to Honduras => Rename San Jeronimo as San Jerónimo
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "San Jeronimo" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Guatemala", replace = F)] <- "San Jerónimo"
# La Palma (BG_adm2) in Canarias is wrongly assigned to Cuba => Rename La Palma as Palmas, Las
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "La Palma" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Spain", replace = F)] <- "Palmas, Las"
# Tenerife (BG_adm2) in Canarias is wrongly assigned to Colombia => Rename Tenerife as Santa Cruz de Tenerife
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Tenerife" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Spain", replace = F)] <- "Santa Cruz de Tenerife"
# El Cercado (BG_adm2) in Nuevo Leon, Mexico is wrongly assigned to Dominican Republic => Rename El Cercado as Santiago, rename locality as Horsetail Falls, El Cercado
Biogeographic_database_Ponerinae_NA_coords$Locality[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "El Cercado" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Mexico", replace = F)] <- "Horsetail Falls, El Cercado"
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "El Cercado" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Mexico", replace = F)] <- "Santiago"
# Natal (BG_adm2) in KwaZulu-Natal, South Africa is wrongly assigned to Brazil => Remove Natal to use adm1 KwaZulu-Natal instead
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Natal" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "South Africa", replace = F)] <- NA
# Concepcion (BG_adm2) in Bio Bio, Chile is wrongly assigned to Peru => Rename Concepcion as Provincia de Concepción
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Concepcion" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Chile", replace = F)] <- "Provincia de Concepción"
# Florida (BG_adm2) in US is wrongly assigned to Bolivia => Remove Florida (adm2) and use adm1 Florida instead
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Florida" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "United States", replace = F)] <- NA
# VILLETA (BG_adm2) in Colombia is wrongly assigned to Paraguay => Rename VILLETA as Villeta
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "VILLETA" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Colombia", replace = F)] <- "Villeta"
# Santiago (BG_adm2) in Los Leones, Chile is wrongly assigned to Brazil => Rename Santiago as Provincia de Santiago
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Santiago" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Chile", replace = F)] <- "Provincia de Santiago"
# York (BG_adm2) in Ontario, Canada is wrongly assigned to UK => Rename York as Toronto
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "York" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Canada", replace = F)] <- "Toronto"
# Cordoba (BG_adm2) in Veracruz, Mexico is wrongly assigned to Spain => Rename Cordoba (adm2) as Córdoba, and rename Veracruz (adm1) as Veracruz de Ignacio de la Llave
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Cordoba" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Mexico", replace = F)] <- "Veracruz de Ignacio de la Llave"
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Cordoba" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Mexico", replace = F)] <- "Córdoba"
# Baoshan (BG_adm2) in Yunnan, China is wrongly assigned to Taiwan => Rename Baoshan as Baoshanshi
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(Biogeographic_database_Ponerinae_NA_coords$adm2 == "Baoshan" & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "China", replace = F)] <- "Baoshanshi"

# Locality with "Luganville" in Vanuatu is adm2 Luganville and adm1 Sanma
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_NA_coords$Locality, pattern = "Luganville") & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Vanuatu", replace = F)] <- "Luganville"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_NA_coords$Locality, pattern = "Luganville") & Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name == "Vanuatu", replace = F)] <- "Sanma"

# Kepulauan Yapen Regency (GB_adm2) is Kepulauan Yapen
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$adm2 == "Kepulauan Yapen Regency", replace = F)] <- "Kepulauan Yapen"

 
# View(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeGroup == "TZA", ])
# View(geoBoundaries_adm1_sf[geoBoundaries_adm1_sf$adm1_shapeGroup == "TZA", ])


## Save cleaned NA_coords database before drawing new random coordinates
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

# Load NA_coords database
Biogeographic_database_Ponerinae_NA_coords <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")


## 3.3.2/ Find matching localities in Geo Server Names ####

## To messy... Maybe go back to it for things that have no match in adm2 (and are not duplicates)

i <- 1
GNS_Localities_sf$full_name[str_detect(string = GNS_Localities_sf$full_name, pattern = paste0("(?i)",Biogeographic_database_Ponerinae_NA_coords$adm2[i]))]

GNS_Localities_sf$full_name
GNS_Localities_sf$full_nm_nd

## 3.3.3/ Find matching adm2 in GeoBoundaries ####

geoBoundaries_adm2_sf$adm2_shapeName[str_detect(string = geoBoundaries_adm2_sf$adm2_shapeName, pattern = paste0("(?i)",Biogeographic_database_Ponerinae_NA_coords$adm2[i]))]

table(Biogeographic_database_Ponerinae_NA_coords$adm2 %in% geoBoundaries_adm2_sf$adm2_shapeName)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_shp <- NA
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2 <- (Biogeographic_database_Ponerinae_NA_coords$adm2 %in% geoBoundaries_adm2_sf$adm2_shapeName)
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "GB_adm2"

## 3.3.4/ Find matching adm1 in GeoBoundaries ####

table(Biogeographic_database_Ponerinae_NA_coords$adm1 %in% geoBoundaries_adm1_sf$adm1_shapeName)

# Only for entries with no match for adm2
table(Biogeographic_database_Ponerinae_NA_coords$adm1[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% geoBoundaries_adm1_sf$adm1_shapeName)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 <- (Biogeographic_database_Ponerinae_NA_coords$adm1 %in% geoBoundaries_adm1_sf$adm1_shapeName)
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "GB_adm1"

table(Biogeographic_database_Ponerinae_NA_coords$matching_shp)

## 3.3.6/ Find matching adm1 in Natural Earth ####

table(Biogeographic_database_Ponerinae_NA_coords$adm1 %in% NE_adm1_sf$name)

# Only for entries with no match for adm2 and GeoBoundaries_adm1
table(Biogeographic_database_Ponerinae_NA_coords$adm1[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% NE_adm1_sf$name)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1 <- (Biogeographic_database_Ponerinae_NA_coords$adm1 %in% NE_adm1_sf$name)
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "NE_adm1"

table(Biogeographic_database_Ponerinae_NA_coords$matching_shp)

## 3.3.5/ Find matching Bentity2 ####

table(Biogeographic_database_Ponerinae_NA_coords$bentity2_name %in% Bentity2_sf$BENTITY2_N)

# Only for entries with no match for adm2 and adm1
table(Biogeographic_database_Ponerinae_NA_coords$bentity2_name[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% Bentity2_sf$BENTITY2_N)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2 <- (Biogeographic_database_Ponerinae_NA_coords$bentity2_name %in% Bentity2_sf$BENTITY2_N)
table(Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2, Biogeographic_database_Ponerinae_NA_coords$Source)
# All GABI entries have a match
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2 & !Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "Bentity2"

table(Biogeographic_database_Ponerinae_NA_coords$matching_shp)
table(is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp))

## 3.3.6/ Find matching SUBUNIT from Natural Earth ####

table(Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name %in% NE_subunits_sf$SUBUNIT)

# plot(x = NE_subunits_sf[, "SUBUNIT"])

# Only for entries with no match for adm2, adm1, and Bentity2
table(Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% NE_subunits_sf$SUBUNIT)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_SUBUNIT <- (Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name %in% NE_subunits_sf$SUBUNIT)
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_SUBUNIT & !Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2 & !Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "NE_SUBUNIT"

table(Biogeographic_database_Ponerinae_NA_coords$matching_shp)
table(is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp))

## 3.3.7/ Find matching Country ADMIN from Natural Earth ####

table(Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name %in% NE_Countries_sf$ADMIN)

# Only for entries with no match for adm2, adm1, Bentity2, and Country
table(Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% NE_Countries_sf$ADMIN)
table(Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp)] %in% NE_Countries_sf$ADMIN)

# Record level of matching
Biogeographic_database_Ponerinae_NA_coords$matching_Country_Admin <- (Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name %in% NE_Countries_sf$ADMIN)
Biogeographic_database_Ponerinae_NA_coords$matching_shp[Biogeographic_database_Ponerinae_NA_coords$matching_Country_Admin & !Biogeographic_database_Ponerinae_NA_coords$matching_SUBUNIT & !Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2 & !Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1 & !Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2] <- "NE_ADMIN"

table(Biogeographic_database_Ponerinae_NA_coords$matching_shp)
table(is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp))

### 3.4/ Remove entries with no matches ####

# Stuff from the Antilles with no more information
View(Biogeographic_database_Ponerinae_NA_coords[is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp), ])

Biogeographic_database_Ponerinae_NA_coords <- Biogeographic_database_Ponerinae_NA_coords %>% 
  filter(!is.na(Biogeographic_database_Ponerinae_NA_coords$matching_shp))

### 3.5/ Intersect all shapes to be able to select across homonyms based on name of higher level shape

sf_use_s2(FALSE)

# Project CRS to help intersection computation Use Robinson projection (ESRI:54030)? Winkel Tripel projection (ESRI:53042)? 
geoBoundaries_adm2_sf_Robinson <- st_transform(x = geoBoundaries_adm2_sf[, c("adm2_shapeName")], crs = st_crs("ESRI:54030"))
geoBoundaries_adm1_sf_Robinson <- st_transform(x = geoBoundaries_adm1_sf[, c("adm1_shapeName")], crs = st_crs("ESRI:54030"))
NE_adm1_sf_Robinson <- st_transform(x = NE_adm1_sf[, c("name")], crs = st_crs("ESRI:54030"))
Bentity2_sf_Robinson <- st_transform(x = Bentity2_sf[, c("BENTITY2_N")], crs = st_crs("ESRI:54030"))
NE_subunits_sf_Robinson <- st_transform(x = NE_subunits_sf[, c("SUBUNIT")], crs = st_crs("ESRI:54030"))
NE_Countries_sf_Robinson <- st_transform(x = NE_Countries_sf[, c("ADMIN")], crs = st_crs("ESRI:54030"))

# plot(Intersect_all_shp_GB_adm2_levels["adm2_shapeName"])

# Snap to grid of 500 m to avoid issue with intersecting nodes, but only work if data is projected... 
geoBoundaries_adm2_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(geoBoundaries_adm2_sf_Robinson, 500))
geoBoundaries_adm1_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(geoBoundaries_adm1_sf_Robinson, 500))
NE_adm1_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_adm1_sf_Robinson, 500))
Bentity2_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(Bentity2_sf_Robinson, 500))
NE_subunits_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_subunits_sf_Robinson, 500))
NE_Countries_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_Countries_sf_Robinson, 500))

# Add zero buffer may solve issues...
geoBoundaries_adm2_sf_Robinson_grid500 <- st_buffer(geoBoundaries_adm2_sf_Robinson_grid500, 0)
geoBoundaries_adm1_sf_Robinson_grid500 <- st_buffer(geoBoundaries_adm1_sf_Robinson_grid500, 0)
NE_adm1_sf_Robinson_grid500 <- st_buffer(NE_adm1_sf_Robinson_grid500, 0)
Bentity2_sf_Robinson_grid500 <- st_buffer(Bentity2_sf_Robinson_grid500, 0)
NE_subunits_sf_Robinson_grid500 <- st_buffer(NE_subunits_sf_Robinson_grid500, 0)
NE_Countries_sf_Robinson_grid500 <- st_buffer(NE_Countries_sf_Robinson_grid500, 0)


## Intersect all shapes from adm2 to country ADMIN
Intersect_all_shp_GB_adm2_levels <- st_make_valid(st_intersection(geoBoundaries_adm2_sf_Robinson_grid500, geoBoundaries_adm1_sf_Robinson_grid500))
# Add 0° buffer to avoid issue with intersecting nodes
Intersect_all_shp_GB_adm2_levels <- st_buffer(Intersect_all_shp_GB_adm2_levels, 0)
table(sf::st_is_valid(Intersect_all_shp_GB_adm2_levels))
# Intersect with higher shape levels
Intersect_all_shp_GB_adm2_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm2_levels, NE_adm1_sf_Robinson_grid500)), 0)
Intersect_all_shp_GB_adm2_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm2_levels, Bentity2_sf_Robinson_grid500)), 0)
Intersect_all_shp_GB_adm2_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm2_levels, NE_subunits_sf_Robinson_grid500)), 0)
Intersect_all_shp_GB_adm2_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm2_levels, NE_Countries_sf_Robinson_grid500)), 0)

## Intersect all shapes from GB_adm1 to country ADMIN
Intersect_all_shp_GB_adm1_levels <- st_buffer(st_make_valid(st_intersection(geoBoundaries_adm1_sf_Robinson_grid500, NE_adm1_sf_Robinson_grid500)), 0)
# Intersect with higher shape levels
Intersect_all_shp_GB_adm1_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm1_levels, Bentity2_sf_Robinson_grid500)), 0)
Intersect_all_shp_GB_adm1_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm1_levels, NE_subunits_sf_Robinson_grid500)), 0)
Intersect_all_shp_GB_adm1_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_GB_adm1_levels, NE_Countries_sf_Robinson_grid500)), 0)

## Intersect all shapes from NE_adm1 to country ADMIN
Intersect_all_shp_NE_adm1_levels <- st_buffer(st_make_valid(st_intersection(NE_adm1_sf_Robinson_grid500, Bentity2_sf_Robinson_grid500)), 0)
# Intersect with higher shape levels
Intersect_all_shp_NE_adm1_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_NE_adm1_levels, NE_subunits_sf_Robinson_grid500)), 0)
Intersect_all_shp_NE_adm1_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_NE_adm1_levels, NE_Countries_sf_Robinson_grid500)), 0)

## Intersect all shapes from Bentity2 to country ADMIN
Intersect_all_shp_Bentity2_levels <- st_buffer(st_make_valid(st_intersection(Bentity2_sf_Robinson_grid500, NE_subunits_sf_Robinson_grid500)), 0)
# Intersect with higher shape levels
Intersect_all_shp_Bentity2_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_Bentity2_levels, NE_Countries_sf_Robinson_grid500)), 0)

## Intersect all shapes from NE_subunits to country ADMIN
Intersect_all_shp_NE_subunits_levels <- st_buffer(st_make_valid(st_intersection(NE_subunits_sf_Robinson_grid500, NE_Countries_sf_Robinson_grid500)), 0)

## Compute areas
Intersect_all_shp_GB_adm2_levels$area <- st_area(Intersect_all_shp_GB_adm2_levels)
Intersect_all_shp_GB_adm2_levels$area <- set_units(Intersect_all_shp_GB_adm2_levels$area, km^2)
Intersect_all_shp_GB_adm1_levels$area <- st_area(Intersect_all_shp_GB_adm1_levels)
Intersect_all_shp_GB_adm1_levels$area <- set_units(Intersect_all_shp_GB_adm1_levels$area, km^2)
Intersect_all_shp_NE_adm1_levels$area <- st_area(Intersect_all_shp_NE_adm1_levels)
Intersect_all_shp_NE_adm1_levels$area <- set_units(Intersect_all_shp_NE_adm1_levels$area, km^2)
Intersect_all_shp_Bentity2_levels$area <- st_area(Intersect_all_shp_Bentity2_levels)
Intersect_all_shp_Bentity2_levels$area <- set_units(Intersect_all_shp_Bentity2_levels$area, km^2)
Intersect_all_shp_NE_subunits_levels$area <- st_area(Intersect_all_shp_NE_subunits_levels)
Intersect_all_shp_NE_subunits_levels$area <- set_units(Intersect_all_shp_NE_subunits_levels$area, km^2)

## Transform back to WGS84
Intersect_all_shp_GB_adm2_levels <- st_transform(x = Intersect_all_shp_GB_adm2_levels, crs = st_crs(geoBoundaries_adm2_sf))
plot(Intersect_all_shp_GB_adm2_levels["adm2_shapeName"])
Intersect_all_shp_GB_adm1_levels <- st_transform(x = Intersect_all_shp_GB_adm1_levels, crs = st_crs(geoBoundaries_adm2_sf))
plot(Intersect_all_shp_GB_adm1_levels["adm1_shapeName"])
Intersect_all_shp_NE_adm1_levels <- st_transform(x = Intersect_all_shp_NE_adm1_levels, crs = st_crs(geoBoundaries_adm2_sf))
plot(Intersect_all_shp_NE_adm1_levels["name"])
Intersect_all_shp_Bentity2_levels <- st_transform(x = Intersect_all_shp_Bentity2_levels, crs = st_crs(geoBoundaries_adm2_sf))
plot(Intersect_all_shp_Bentity2_levels["BENTITY2_N"])
Intersect_all_shp_NE_subunits_levels <- st_transform(x = Intersect_all_shp_NE_subunits_levels, crs = st_crs(geoBoundaries_adm2_sf))
plot(Intersect_all_shp_NE_subunits_levels["SUBUNIT"])

## Remove empty polygons
Intersect_all_shp_GB_adm2_levels <- Intersect_all_shp_GB_adm2_levels[!st_is_empty(Intersect_all_shp_GB_adm2_levels), , drop = FALSE]
Intersect_all_shp_GB_adm1_levels <- Intersect_all_shp_GB_adm1_levels[!st_is_empty(Intersect_all_shp_GB_adm1_levels), , drop = FALSE]
Intersect_all_shp_NE_adm1_levels <- Intersect_all_shp_NE_adm1_levels[!st_is_empty(Intersect_all_shp_NE_adm1_levels), , drop = FALSE]
Intersect_all_shp_Bentity2_levels <- Intersect_all_shp_Bentity2_levels[!st_is_empty(Intersect_all_shp_Bentity2_levels), , drop = FALSE]
Intersect_all_shp_NE_subunits_levels <- Intersect_all_shp_NE_subunits_levels[!st_is_empty(Intersect_all_shp_NE_subunits_levels), , drop = FALSE]

sf_use_s2(TRUE)

## Save intersecting shp files
saveRDS(object = Intersect_all_shp_GB_adm2_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_GB_adm2_levels.rds")
saveRDS(object = Intersect_all_shp_GB_adm1_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_GB_adm1_levels.rds")
saveRDS(object = Intersect_all_shp_NE_adm1_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_NE_adm1_levels.rds")
saveRDS(object = Intersect_all_shp_Bentity2_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_Bentity2_levels.rds")
saveRDS(object = Intersect_all_shp_NE_subunits_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_NE_subunits_levels.rds")

# Load intersecting shp files
Intersect_all_shp_GB_adm2_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_GB_adm2_levels.rds")
Intersect_all_shp_GB_adm1_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_GB_adm1_levels.rds")
Intersect_all_shp_NE_adm1_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_NE_adm1_levels.rds")
Intersect_all_shp_Bentity2_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_Bentity2_levels.rds")
Intersect_all_shp_NE_subunits_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_NE_subunits_levels.rds")

### 3.5/ Draw random coordinates in new polygon ####

Biogeographic_database_Ponerinae_NA_coords$Interpolated_occurrence <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area <- NA

# List all sf files
all_sf_files <- list(geoBoundaries_adm2_sf[, c("adm2_shapeName", "area")], 
                     geoBoundaries_adm1_sf[, c("adm1_shapeName", "area")], 
                     NE_adm1_sf[, c("name", "area")], 
                     Bentity2_sf[, c("BENTITY2_N", "area")], 
                     NE_subunits_sf[, c("SUBUNIT", "area")], 
                     NE_Countries_sf[, c("ADMIN", "area")])
# List all intersecting sf files (without non-intersecting polygons)
all_intersecting_sf_files <- list(Intersect_all_shp_GB_adm2_levels, 
                                  Intersect_all_shp_GB_adm1_levels, 
                                  Intersect_all_shp_NE_adm1_levels, 
                                  Intersect_all_shp_Bentity2_levels, 
                                  Intersect_all_shp_NE_subunits_levels, 
                                  NE_Countries_sf[, c("ADMIN", "area")])
# List all matching levels
matching_level_list <- c("GB_adm2", "GB_adm1", "NE_adm1", "Bentity2", "NE_SUBUNIT", "NE_ADMIN")

### 3.5.1/ Loop per entry #### 

set.seed(seed = 644668)

Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID <- NA
Biogeographic_database_Ponerinae_NA_coords$multiple_matches <- FALSE # Record cases of entries with multiple name matching, so the shape must be extracted from the file with all intersected levels of shapes

# Biogeographic_database_Ponerinae_NA_coords <- st_drop_geometry(Biogeographic_database_Ponerinae_NA_coords)

sf_use_s2(FALSE)
for (i in 1:nrow(Biogeographic_database_Ponerinae_NA_coords))
# for (i in 34400:nrow(Biogeographic_database_Ponerinae_NA_coords))
# for (i in occ_to_update)
{ 
  # i <- 1
  # i <- 466
  
  # which(Biogeographic_database_Ponerinae_NA_coords$adm1 == "San Salvador")
  

  # fix <- F
  
  # Extract matching level
  matching_lvl <- Biogeographic_database_Ponerinae_NA_coords$matching_shp[i]
  
  # Extract name of matching shape
  focal_occ_shp_names <- Biogeographic_database_Ponerinae_NA_coords[i, c("adm2", "adm1", "adm1", "bentity2_name", "SUBUNIT_name", "Country_ISO3_name")]
  matching_level_binary <- Biogeographic_database_Ponerinae_NA_coords[i, c("matching_GB_adm2", "matching_GB_adm1", "matching_NE_adm1", "matching_Bentity2", "matching_SUBUNIT", "matching_Country_Admin")]
  matching_shp_index <- which(matching_lvl == matching_level_list)
  focal_occ_shp_name <- focal_occ_shp_names[matching_shp_index][1,1]
  
  # Retrieve associated shape
  focal_occ_sf_file <- all_sf_files[[matching_shp_index]]
  focal_occ_shp_index <- which(st_drop_geometry(focal_occ_sf_file[, 1]) == focal_occ_shp_name)
  
  # Need to decide what to do for cases with multiple matches
  if(length(focal_occ_shp_index) != 1)
  {  
    Biogeographic_database_Ponerinae_NA_coords$multiple_matches[i] <- TRUE
    
    # Extract shp files with intersecting all other higher levels
    Intersect_higher_shp_levels <- all_intersecting_sf_files[[matching_shp_index]]
    
    # Check for other matches in the shape files intersecting all higher levels
    focal_occ_all_shp_indices <- which(st_drop_geometry(Intersect_higher_shp_levels[, 1]) == focal_occ_shp_name)
    focal_occ_all_shp_matches <- Intersect_higher_shp_levels[focal_occ_all_shp_indices, ]
    # plot(focal_occ_all_shp_matches[, 1])
    
    # Find the shape(s) with the most matches
    nb_matches <- apply(X = st_drop_geometry(focal_occ_all_shp_matches), MARGIN = 1, FUN = function (x) { sum(apply(X = sapply(X = x, `==`, focal_occ_shp_names), MARGIN = 2, FUN = any, na.rm = TRUE))  })
    focal_occ_all_shp_best_matches <- focal_occ_all_shp_matches [nb_matches == max(nb_matches), ]   
    # In case of equality, use the one with the highest area, because some shapes are likely artefact of poorly matching boundaries across data source
    focal_occ_all_shp_best_match <- arrange(focal_occ_all_shp_best_matches, desc(area))[1, ]
    # plot(focal_occ_all_shp_best_match[, 1])
    
    # Extract the ID of the shp file in Intersect_higher_shp_levels
    focal_occ_shp_index <- which(row.names(focal_occ_all_shp_best_match) == row.names(Intersect_higher_shp_levels))
   
    # Record shp ID and extract shape
    Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[i] <- focal_occ_shp_index
    focal_occ_shp <- focal_occ_all_shp_best_match
    
    # # Florida in the USA, not in Uruguay
    # if ((focal_occ_shp_names["adm1"] == "Florida") & (focal_occ_shp_names["Country_ISO3_name"] == "United States")) { focal_occ_shp_index <- 108 ; fix <- T }
    # if (!fix) { stop(paste0("i = ",i))}

  } else {
    
    # Record shp ID and extract shape
    Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[i] <- focal_occ_shp_index
    focal_occ_shp <- focal_occ_sf_file[focal_occ_shp_index, ]
    
  }
  
  # Sample randomly a point
  focal_occ_point_sf <- suppressMessages(suppressWarnings(sf::st_sample(x = focal_occ_shp, size = 1, exact = TRUE)))
  # Extract coordinates
  focal_occ_longlat <- as.numeric(focal_occ_point_sf[[1]])
  # Record coordinates
  Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[i] <- focal_occ_longlat[1]
  Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[i] <- focal_occ_longlat[2]
  
  # Inform of uncertainty based on area of the polygon
  Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[i] <- focal_occ_shp$area
  
  # Print progress
  if ( i %% 100 == 0) { cat(paste0(Sys.time(), " - Occurrence n°",i, " / ", nrow(Biogeographic_database_Ponerinae_NA_coords),"\n")) }
}
sf_use_s2(TRUE)

# Transform into sf object
Biogeographic_database_Ponerinae_NA_coords$Latitude_copy <- Biogeographic_database_Ponerinae_NA_coords$Latitude_dec
Biogeographic_database_Ponerinae_NA_coords$Longitude_copy <- Biogeographic_database_Ponerinae_NA_coords$Longitude_dec
Biogeographic_database_Ponerinae_NA_coords <- st_as_sf(Biogeographic_database_Ponerinae_NA_coords, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_NA_coords) <- st_crs(geoBoundaries_adm2_sf)
# Remove copies of Long/Lat if still present
# Biogeographic_database_Ponerinae_NA_coords <- Biogeographic_database_Ponerinae_NA_coords %>% select(-c("Longitude_copy", "Latitude_copy"))

# Check results
plot(NE_Countries_sf[, c("ADMIN")])
plot(Biogeographic_database_Ponerinae_NA_coords[, "Country_ISO3_name"], add = T)
plot(Biogeographic_database_Ponerinae_NA_coords[, "Country_ISO3_name"])

test_map <- ggplot() + 
  geom_sf(data = NE_Countries_sf, fill = NA, show.legend = FALSE) + 
  geom_sf(data = Biogeographic_database_Ponerinae_NA_coords, mapping = aes(col = Country_ISO3_name), show.legend = FALSE) +
  theme_minimal()
print(test_map)

# Save database of NA coordinates
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")


### 3.5.2/ Double-check for potential errors by checking country compatibility using CoordinatesCleaner ####

# Load landmass buffer
map_NE_lands_sf_buffer <- readRDS(file = "./input_data/NaturalEarth_maps/map_NE_lands_sf_buffer.rds")

# Scan new random coordinates
Biogeographic_database_Ponerinae_NA_coords_flags <- clean_coordinates(x = st_drop_geometry(Biogeographic_database_Ponerinae_NA_coords),
                                                            lon = "Longitude_dec", lat = "Latitude_dec",
                                                            countries = "Country_ISO3_code",
                                                            species = "Current_name",
                                                            value = "spatialvalid",         # To define the type of output. Can be a spatial df ("spatialvalid"), a logical vector recording entries with at least one flag ("flagged"), a df with dubious coordinates removed ("clean").
                                                            verbose = T,
                                                            # List of tests to run
                                                            tests = c("countries",     # Test if occurrences are found within a list of countries or a custom SpatialPolygon
                                                                      "seas"           # Test if fall in the ocean
                                                            ), 
                                                            # seas_scale = 10,               # Adjust resolution of landmass polygons used as reference. 10 = highest resolution.
                                                            seas_ref = map_NE_lands_sf_buffer,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                                                            # seas_ref = map_NE_all_lands_sf,   # sf SpatialPolygons to use as reference for landmasses. Default is too coarse. 
                                                            country_ref = Countries_NE_sf,  # sf SpatialPolygons to use as reference for countries
                                                            country_refcol = "ISO_A3",     # Name of ISO3 column
                                                            country_buffer = 5000         # Buffer around country polygons (in meters)
                                                            # seas_buffer = 5000,            # Buffer around landmass limits (in meters) # Does not work properly with scale 10 map. Make our own prior
                                                            )

table(Biogeographic_database_Ponerinae_NA_coords_flags$.con)
table(Biogeographic_database_Ponerinae_NA_coords_flags$.sea)
table(Biogeographic_database_Ponerinae_NA_coords_flags$.con, Biogeographic_database_Ponerinae_NA_coords_flags$.sea)

View(Biogeographic_database_Ponerinae_NA_coords_flags[!Biogeographic_database_Ponerinae_NA_coords_flags$.con, ])

Biogeographic_database_Ponerinae_NA_coords_flags$flag_type <- "OK"
Biogeographic_database_Ponerinae_NA_coords_flags$flag_type[!Biogeographic_database_Ponerinae_NA_coords_flags$.con] <- "country"
Biogeographic_database_Ponerinae_NA_coords_flags$flag_type[!Biogeographic_database_Ponerinae_NA_coords_flags$.sea] <- "sea"

table(Biogeographic_database_Ponerinae_NA_coords_flags$flag_type)

## Make an interactive plot to spot issues

# Transform into sf object
Biogeographic_database_Ponerinae_NA_coords_flags$Latitude_copy <- Biogeographic_database_Ponerinae_NA_coords_flags$Latitude_dec
Biogeographic_database_Ponerinae_NA_coords_flags$Longitude_copy <- Biogeographic_database_Ponerinae_NA_coords_flags$Longitude_dec
Biogeographic_database_Ponerinae_NA_coords_flags <- st_as_sf(Biogeographic_database_Ponerinae_NA_coords_flags, coords = c("Longitude_copy", "Latitude_copy"), crs = sp::CRS('+init=EPSG:4326'))
# Set CRS
st_crs(Biogeographic_database_Ponerinae_NA_coords_flags) <- st_crs(geoBoundaries_adm2_sf)

# Save database of NA coords with flags
saveRDS(Biogeographic_database_Ponerinae_NA_coords_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords_flags.rds")


nb_flags <- length(unique(Biogeographic_database_Ponerinae_NA_coords_flags$flag_type[]))
nb_matching_shp <- length(unique(Biogeographic_database_Ponerinae_NA_coords_flags$matching_shp[]))

# Define palette function
pal <- mapviewPalette("mapviewSpectralColors")

# RColorBrewer::display.brewer.all()

# pal <- brewer.pal(n = nb_taxa, name = "Spectral")
pal_flags <- brewer.pal(n = nb_flags, name = "Spectral")
pal_matching_shp <- brewer.pal(n = nb_matching_shp, name = "Spectral")
# pal_4 <- pal[c(1, 4, 5, 6)]
# pal_2 <- pal[c(5, 6)]

plot(Biogeographic_database_Ponerinae_NA_coords_flags[, ]["flag_type"])

### With color per flag types (2 flag types)
interactive_map_flags <- mapview::mapview(x = Biogeographic_database_Ponerinae_NA_coords_flags[Biogeographic_database_Ponerinae_NA_coords_flags$flag_type != "OK", ], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                          zcol = "flag_type",
                                          cex = 5,
                                          burst = T, legend = T,
                                          col.regions = pal_flags, # Color to use to plot spatial units in the spatial object
                                          # col.regions = pal_2, # Color to use to plot spatial units in the spatial object
                                          layer.name = "Occurrences", # To specify the legend name of this layer
                                          # map.types = mapviewGetOption("basemaps")
                                          map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                          # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                          # ... # and many more, depending of the class of the spatial object
)

interactive_map_flags


### With color per matching levels (6 matching levels)
interactive_map_matching_levels <- mapview::mapview(x = Biogeographic_database_Ponerinae_NA_coords_flags[Biogeographic_database_Ponerinae_NA_coords_flags$flag_type != "OK", ], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                          zcol = "matching_shp",
                                          cex = 5,
                                          burst = T, legend = T,
                                          col.regions = pal_matching_shp, # Color to use to plot spatial units in the spatial object
                                          # col.regions = pal_2, # Color to use to plot spatial units in the spatial object
                                          layer.name = "Occurrences", # To specify the legend name of this layer
                                          # map.types = mapviewGetOption("basemaps")
                                          map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                          # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                          # ... # and many more, depending of the class of the spatial object
)

interactive_map_matching_levels


### 3.5.3/ Manually adjust occurrences falling in seas/oceans ####

# Tokelau, adjusting Longitude towards E
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00761219", replace = F)] <- -171.184

# Clipperton Island, adjusting Longitude towards W
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00419641", replace = F)] <- -109.232

# Christmas Island, adjusting towards SW
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00419755", replace = F)] <- 105.686
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00419755", replace = F)] <- -10.450
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00756976", replace = F)] <- 105.626
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00756976", replace = F)] <- -10.482

# Yaeyama Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00618780", replace = F)] <- 24.055
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00878895", replace = F)] <- 24.233
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00370786", replace = F)] <- 24.243
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00370786", replace = F)] <- 124.0212

# Toshima Island, Izu Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00876563", replace = F)] <- 34.52001245
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00451303", replace = F)] <- 34.52511473
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00876846", replace = F)] <- 34.51973486
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00349130", replace = F)] <- 34.52636874
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988676", replace = F)] <- 34.524
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988676", replace = F)] <- 139.275
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988698", replace = F)] <- 34.520
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988698", replace = F)] <- 139.281
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988700", replace = F)] <- 34.518
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988700", replace = F)] <- 139.287

# Kōzushima, Izu Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00876803", replace = F)] <- 34.194
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988688", replace = F)] <- 34.227
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00876857", replace = F)] <- 34.212
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988685", replace = F)] <- 34.2107

# Hachijoko Island, Izu Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00349250", replace = F)] <- 33.123

# Aogashima, Izu Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00365707", replace = F)] <- 32.457
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00365707", replace = F)] <- 139.768

# Torishima Island, Izu Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988694", replace = F)] <- 30.483
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988686", replace = F)] <- 30.482
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988686", replace = F)] <- 140.310

# Nishinoshima Island, Ogasawara Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348195", replace = F)] <- 27.248
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348195", replace = F)] <- 140.875

# Muko Island, Ogasawara Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348158", replace = F)] <- 27.685

# Nakodojima Island, Ogasawara Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00497313", replace = F)] <- 27.628
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00497313", replace = F)] <- 142.178
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00878467", replace = F)] <- 27.626
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00877245", replace = F)] <- 27.6302
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348107", replace = F)] <- 27.629
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348107", replace = F)] <- 142.182
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00347016", replace = F)] <- 27.6226
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00347016", replace = F)] <- 142.185

# Mejima Island, Ogasawara Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348145", replace = F)] <- 26.567

# Minamitori-shima, Ogasawara Islands, Japan
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348119", replace = F)] <- 24.2905
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348119", replace = F)] <- 153.9829
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988973", replace = F)] <- 24.2872
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00988973", replace = F)] <- 153.9811
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348128", replace = F)] <- 24.2831
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$GABI_accession_ID == "GABI_00348128", replace = F)] <- 153.9844

## Rerun the flagging to check if the correction worked

# Save database of NA coordinates
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

# Load database of NA coordinates
# Biogeographic_database_Ponerinae_NA_coords <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")


# View(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeGroup == "DZA", ])
# View(geoBoundaries_adm1_sf[geoBoundaries_adm1_sf$adm1_shapeGroup == "DZA", ])
# 
# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName %in% c("East District", "North District", "West District", "South District"), "adm2_shapeName"])
# 
# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Marovo", "adm2_shapeName"])
# plot(geoBoundaries_adm1_sf[geoBoundaries_adm1_sf$adm1_shapeName == "Sikkim", "adm1_shapeName"])


### 3.6/ Record new polygons ####

### 3.6.1/ Build shp files with all intersecting shp levels, including non-intersecting shp ####

# Use this instead of previous intersecting shp files because it preserve polygons that are not intersecting on all levels
# Will provide an NA if the level does not intersect

sf_use_s2(FALSE)

# Project CRS to help intersection computation Use Robinson projection (ESRI:54030)? Winkel Tripel projection (ESRI:53042)? 
geoBoundaries_adm2_sf_Robinson <- st_transform(x = geoBoundaries_adm2_sf[, c("adm2_shapeName")], crs = st_crs("ESRI:54030"))
geoBoundaries_adm1_sf_Robinson <- st_transform(x = geoBoundaries_adm1_sf[, c("adm1_shapeName")], crs = st_crs("ESRI:54030"))
NE_adm1_sf_Robinson <- st_transform(x = NE_adm1_sf[, c("name")], crs = st_crs("ESRI:54030"))
Bentity2_sf_Robinson <- st_transform(x = Bentity2_sf[, c("BENTITY2_N")], crs = st_crs("ESRI:54030"))
NE_subunits_sf_Robinson <- st_transform(x = NE_subunits_sf[, c("SUBUNIT")], crs = st_crs("ESRI:54030"))
NE_Countries_sf_Robinson <- st_transform(x = NE_Countries_sf[, c("ADMIN")], crs = st_crs("ESRI:54030"))

# plot(Intersect_all_shp_GB_adm2_levels["adm2_shapeName"])

# Snap to grid of 500 m to avoid issue with intersecting nodes, but only work if data is projected... 
geoBoundaries_adm2_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(geoBoundaries_adm2_sf_Robinson, 500))
geoBoundaries_adm1_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(geoBoundaries_adm1_sf_Robinson, 500))
NE_adm1_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_adm1_sf_Robinson, 500))
Bentity2_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(Bentity2_sf_Robinson, 500))
NE_subunits_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_subunits_sf_Robinson, 500))
NE_Countries_sf_Robinson_grid500 <- st_make_valid(lwgeom::st_snap_to_grid(NE_Countries_sf_Robinson, 500))

# Add zero buffer may solve issues...
geoBoundaries_adm2_sf_Robinson_grid500 <- st_buffer(geoBoundaries_adm2_sf_Robinson_grid500, 0)
geoBoundaries_adm1_sf_Robinson_grid500 <- st_buffer(geoBoundaries_adm1_sf_Robinson_grid500, 0)
NE_adm1_sf_Robinson_grid500 <- st_buffer(NE_adm1_sf_Robinson_grid500, 0)
Bentity2_sf_Robinson_grid500 <- st_buffer(Bentity2_sf_Robinson_grid500, 0)
NE_subunits_sf_Robinson_grid500 <- st_buffer(NE_subunits_sf_Robinson_grid500, 0)
NE_Countries_sf_Robinson_grid500 <- st_buffer(NE_Countries_sf_Robinson_grid500, 0)

# Get intersection between GB_adm2 and GB_adm1
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(st_intersection(geoBoundaries_adm2_sf_Robinson_grid500, geoBoundaries_adm1_sf_Robinson_grid500)), 0)
# Get differences to record non-intersecting polygons in GB_adm2
Intersect_all_shp_all_levels_GB_adm2_diff <- st_buffer(st_make_valid(st_difference(geoBoundaries_adm2_sf_Robinson_grid500, st_union(st_geometry(Intersect_all_shp_all_levels)))), 0)
# Get differences to record non-intersecting polygons in GB_adm1
Intersect_all_shp_all_levels_GB_adm1_diff <- st_buffer(st_make_valid(st_difference(geoBoundaries_adm1_sf_Robinson_grid500, st_union(st_geometry(Intersect_all_shp_all_levels)))), 0)
# Bind all polygons (intersecting and non-intersecting)
Intersect_all_shp_all_levels <- dplyr::bind_rows(Intersect_all_shp_all_levels, Intersect_all_shp_all_levels_GB_adm2_diff, Intersect_all_shp_all_levels_GB_adm1_diff)

# plot(Intersect_all_shp_all_levels_GB_adm2_diff[, "adm2_shapeName"])
# plot(Intersect_all_shp_all_levels_GB_adm1_diff[, "adm1_shapeName"])
# plot(Intersect_all_shp_all_levels[, "area"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")
saveRDS(object = Intersect_all_shp_all_levels_GB_adm2_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_GB_adm2_diff.rds")
saveRDS(object = Intersect_all_shp_all_levels_GB_adm1_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_GB_adm1_diff.rds")


# Get intersection between GB_amd2/GB_adm1 and NE_adm1
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_all_levels, NE_adm1_sf_Robinson_grid500)), 0)
# Get differences to record non-intersecting polygons in NE_adm1
Intersect_all_shp_all_levels_NE_adm1_diff <- st_buffer(st_make_valid(st_difference(NE_adm1_sf_Robinson_grid500, st_union(st_geometry(Intersect_all_shp_all_levels)))), 0)
# Bind all polygons (intersecting and non-intersecting)
Intersect_all_shp_all_levels <- dplyr::bind_rows(Intersect_all_shp_all_levels, Intersect_all_shp_all_levels_NE_adm1_diff)

# plot(Intersect_all_shp_all_levels_NE_adm1_diff[, "name"])
# plot(Intersect_all_shp_all_levels[, "area"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")
saveRDS(object = Intersect_all_shp_all_levels_NE_adm1_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_NE_adm1_diff.rds")


# Get intersection between GB_amd2/GB_adm1/NE_adm1 and Bentity2
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_all_levels, Bentity2_sf_Robinson_grid500)), 0)
# Get differences to record non-intersecting polygons in Bentity2
Intersect_all_shp_all_levels_Bentity2_diff <- st_buffer(st_make_valid(st_difference(Bentity2_sf_Robinson_grid500, st_union(st_geometry(Intersect_all_shp_all_levels)))), 0)
# Bind all polygons (intersecting and non-intersecting)
Intersect_all_shp_all_levels <- dplyr::bind_rows(Intersect_all_shp_all_levels, Intersect_all_shp_all_levels_Bentity2_diff)

# plot(Intersect_all_shp_all_levels_Bentity2_diff[, "BENTITY2_N"])
# plot(Intersect_all_shp_all_levels[, "area"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")
saveRDS(object = Intersect_all_shp_all_levels_Bentity2_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_Bentity2_diff.rds")


# Get intersection between GB_amd2/GB_adm1/NE_adm1/Bentity2 and NE_Subunits
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_all_levels, NE_subunits_sf_Robinson_grid500)), 0)
# Get differences to record non-intersecting polygons in NE_subunits
Union_before_Subunits <- st_buffer(st_make_valid(st_union(st_geometry(Intersect_all_shp_all_levels))), 0)
Intersect_all_shp_all_levels_NE_subunits_diff <- st_buffer(st_make_valid(st_difference(NE_subunits_sf_Robinson_grid500, Union_before_Subunits)), 0)
# Bind all polygons (intersecting and non-intersecting)
Intersect_all_shp_all_levels <- dplyr::bind_rows(Intersect_all_shp_all_levels, Intersect_all_shp_all_levels_NE_subunits_diff)

# plot(Union_before_Subunits[, "adm2_shapeName"])
# plot(Intersect_all_shp_all_levels_NE_subunits_diff[, "SUBUNIT"])
# plot(Intersect_all_shp_all_levels[, "area"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")
saveRDS(object = Intersect_all_shp_all_levels_NE_subunits_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_NE_subunits_diff.rds")


# Get intersection between GB_amd2/GB_adm1/NE_adm1/Bentity2/NE_Subunits and NE_Countries
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(st_intersection(Intersect_all_shp_all_levels, NE_Countries_sf_Robinson_grid500)), 0)
# Get differences to record non-intersecting polygons in NE_Countries
Intersect_all_shp_all_levels_NE_Countries_diff <- st_buffer(st_make_valid(st_difference(NE_Countries_sf_Robinson_grid500, st_union(st_geometry(Intersect_all_shp_all_levels)))), 0)
# Bind all polygons (intersecting and non-intersecting)
Intersect_all_shp_all_levels <- dplyr::bind_rows(Intersect_all_shp_all_levels, Intersect_all_shp_all_levels_NE_Countries_diff)

# plot(Intersect_all_shp_all_levels_NE_Countries_diff[, "ADMIN"])
# plot(Intersect_all_shp_all_levels[, "adm2_shapeName"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")
saveRDS(object = Intersect_all_shp_all_levels_NE_Countries_diff, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels_NE_Countries_diff.rds")

## Compute areas
Intersect_all_shp_all_levels$area <- st_area(Intersect_all_shp_all_levels)
Intersect_all_shp_all_levels$area <- set_units(Intersect_all_shp_all_levels$area, km^2)

hist(Intersect_all_shp_all_levels$area)

## Transform back to WGS84
Intersect_all_shp_all_levels <- st_transform(x = Intersect_all_shp_all_levels, crs = st_crs(geoBoundaries_adm2_sf))

## Repair manually invalid shapes that cross the antimeridian
Intersect_all_shp_all_levels_validity_test <- st_is_valid(Intersect_all_shp_all_levels)
table(Intersect_all_shp_all_levels_validity_test) # 3827 shapes to repair
table(is.na(Intersect_all_shp_all_levels_validity_test)) # 4 shapes that cannot be repaired and need to be manually curated
which(is.na(Intersect_all_shp_all_levels_validity_test))

plot(st_geometry(Intersect_all_shp_all_levels)[[82960]])
# Shape 82960 has a sub-polygon with only 1 point. Needs to be removed manually
sapply(st_geometry(Intersect_all_shp_all_levels)[[82960]], function(x) nrow(x[[1]]))
error_ID <- which(sapply(st_geometry(Intersect_all_shp_all_levels)[[82960]], function(x) nrow(x[[1]])) < 2)
geom_temp <- st_geometry(Intersect_all_shp_all_levels)
geom_temp[[82960]][error_ID] <- NULL
st_geometry(Intersect_all_shp_all_levels) = geom_temp
plot(st_geometry(Intersect_all_shp_all_levels)[[82960]])

plot(st_geometry(Intersect_all_shp_all_levels)[[83171]])
# Shape 83171 has a sub-polygon LinearRing with less than 4 points. Needs to be removed manually
sapply(st_geometry(Intersect_all_shp_all_levels)[[83171]], function(x) nrow(x[[1]]))
error_ID <- which(sapply(st_geometry(Intersect_all_shp_all_levels)[[83171]], function(x) nrow(x[[1]])) < 4)
geom_temp <- st_geometry(Intersect_all_shp_all_levels)
geom_temp[[83171]][error_ID] <- NULL
st_geometry(Intersect_all_shp_all_levels) = geom_temp
plot(st_geometry(Intersect_all_shp_all_levels)[[83171]])

plot(st_geometry(Intersect_all_shp_all_levels)[[251632]])
# Shape 251632 has a sub-polygon LinearRing with less than 4 points. Needs to be removed manually
sapply(st_geometry(Intersect_all_shp_all_levels)[[251632]], function(x) nrow(x[[1]]))
error_ID <- which(sapply(st_geometry(Intersect_all_shp_all_levels)[[251632]], function(x) nrow(x[[1]])) < 4)
geom_temp <- st_geometry(Intersect_all_shp_all_levels)
geom_temp[[251632]][error_ID] <- NULL
st_geometry(Intersect_all_shp_all_levels) = geom_temp
plot(st_geometry(Intersect_all_shp_all_levels)[[251632]])

plot(st_geometry(Intersect_all_shp_all_levels)[[262960]])
# Shape 262960 has a sub-polygon with only 1 point. Needs to be removed manually
sapply(st_geometry(Intersect_all_shp_all_levels)[[262960]], function(x) nrow(x[[1]]))[error_ID]
error_ID <- which(sapply(st_geometry(Intersect_all_shp_all_levels)[[262960]], function(x) nrow(x[[1]])) == 1)
geom_temp <- st_geometry(Intersect_all_shp_all_levels)
geom_temp[[262960]][error_ID] <- NULL
st_geometry(Intersect_all_shp_all_levels) = geom_temp
plot(st_geometry(Intersect_all_shp_all_levels)[[262960]])

## Repair automatically other shapes
Intersect_all_shp_all_levels <- st_buffer(st_make_valid(Intersect_all_shp_all_levels), 0)
Intersect_all_shp_all_levels <- st_make_valid(Intersect_all_shp_all_levels) # To run until no more invalid shapes remain
table(st_is_valid(Intersect_all_shp_all_levels))

## Remove empty polygons
Intersect_all_shp_all_levels <- Intersect_all_shp_all_levels[!st_is_empty(Intersect_all_shp_all_levels), , drop = FALSE]

## Add shape ID
Intersect_all_shp_all_levels$Intersect_shp_ID <- 1:nrow(Intersect_all_shp_all_levels)

plot(Intersect_all_shp_all_levels["Intersect_shp_ID"])

# Save Intersecting files for all levels, including non-intersecting polygons
saveRDS(object = Intersect_all_shp_all_levels, file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")

sf_use_s2(TRUE)


### 3.6.2/ Loop per entry ####

# Load Intersecting files for all levels, including non-intersecting polygons
Intersect_all_shp_all_levels <- readRDS(file = "./input_data/Biogeographic_data/Intersect_all_shp_all_levels.rds")

Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2 <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1 <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_Bentity2 <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2 <- 100
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1 <- 100

entries_with_no_matching_shp <- c()

sf_use_s2(FALSE)
# for (i in 1:nrow(Biogeographic_database_Ponerinae_NA_coords))
# for (i in 50095:nrow(Biogeographic_database_Ponerinae_NA_coords))
for (i in entries_with_no_matching_shp)
{ 
  # i <- 2505
  # i <- 4668
  
  ## Extract matching level and shp ID
  matching_lvl <- Biogeographic_database_Ponerinae_NA_coords$matching_shp[i]
  matching_shp_ID <- Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[i]
  matching_shp_index <- which(matching_lvl == matching_level_list)
  
  ## Extract shape with all levels where coordinates are falling
  new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], Intersect_all_shp_all_levels)))
  
  # If extraction through intersection failed, look for closest shape
  if (nrow(new_shp_all_levels) == 0)
  {
    new_shp_all_levels_ID <- suppressMessages(suppressWarnings(sf::st_nearest_feature(Biogeographic_database_Ponerinae_NA_coords[i, ], Intersect_all_shp_all_levels)))
    new_shp_all_levels <- Intersect_all_shp_all_levels[new_shp_all_levels_ID, ]
  }

  # print(new_shp_all_levels)
  
  ## Test if shape retrieved has the same name as the one used to draw the coordinates
  initial_shp_names <- Biogeographic_database_Ponerinae_NA_coords[i, c("adm2", "adm1", "adm1", "bentity2_name", "SUBUNIT_name", "Country_ISO3_name")]
  new_shp_names <- new_shp_all_levels[1, c("adm2_shapeName", "adm1_shapeName",  "name", "BENTITY2_N", "SUBUNIT",  "ADMIN")]
  initial_shp_name <- st_drop_geometry(initial_shp_names)[1, matching_shp_index]
  new_shp_name <- st_drop_geometry(new_shp_names)[1, matching_shp_index]
  match_test <- replace_na(data = initial_shp_name == new_shp_name, replace = FALSE)
    
  ## Skip an entry of no matching shp if found OR if matching level do not match...
  if ((nrow(new_shp_all_levels) == 0) | !match_test)
  {
    # Record ID of entry with an error
    entries_with_no_matching_shp <- c(entries_with_no_matching_shp, i)
    
    if ((nrow(new_shp_all_levels) == 0)) 
    {
      cat(paste0("WARNING: Entry n°", i, " has no match\n"))
    } else {
      cat(paste0("WARNING: Entry n°", i, " has an error in matching shp names. Initial ",matching_lvl," = ",initial_shp_name,". New ",matching_lvl," = ",new_shp_name,"\n"))
    }
  } else { # If no error, proceed to extract shp names
    
    ## Extract shape used to draw coordinates
    if (!Biogeographic_database_Ponerinae_NA_coords$multiple_matches[i]) # Case with a single match from initial shp file
    {
      focal_occ_sf_file <- all_sf_files[[matching_shp_index]]
    } else { # Case with multiple matches so shape from intersecting shp file
      focal_occ_sf_file <- all_intersecting_sf_files[[matching_shp_index]]
    }
    focal_occ_shp <- focal_occ_sf_file[matching_shp_ID, ]
    
    # plot(focal_occ_shp[, "area"])
    
    ## Record names for all levels others than the one used for matching
    # Also record adm1/adm2 uncertainty based on area overlap with the polygon used to draw random coordinates
    if (matching_lvl %in% c("GB_adm2"))
    { 
      # Record adm1
      Biogeographic_database_Ponerinae_NA_coords$adm1[i] <- st_drop_geometry(new_shp_all_levels[, "adm1_shapeName"])[1,1]
      # Record Bentity2
      Biogeographic_database_Ponerinae_NA_coords$bentity2_name[i] <- st_drop_geometry(new_shp_all_levels[, "BENTITY2_N"])[1,1]
      # Record SUBUNIT
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels[, "SUBUNIT"])[1,1]
      # Record Country Admin
      Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[i] <- st_drop_geometry(new_shp_all_levels[, "ADMIN"])[1,1]
      
      # Uncertainty is null because adm2 was retrieved by name
      Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 0
      # Uncertainty is null because adm1 was retrieved from adm2
      Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 0
    }
    if (matching_lvl %in% c("GB_adm1"))
    {
      # Record adm2 if drawn from adm1
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$adm2[i] <- st_drop_geometry(new_shp_all_levels[, "adm2_shapeName"])[1,1]
      
      # Record Bentity2
      Biogeographic_database_Ponerinae_NA_coords$bentity2_name[i] <- st_drop_geometry(new_shp_all_levels[, "BENTITY2_N"])[1,1]
      # Record SUBUNIT
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels[, "SUBUNIT"])[1,1]
      # Record Country Admin
      Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[i] <- st_drop_geometry(new_shp_all_levels[, "ADMIN"])[1,1]
      
      # adm2 Uncertainty is 100% - the percentage of area covered by GB_adm2 within the NE_SUBUNIT area used to draw the random coordinates
      all_sf_files[[1]]$GB_adm2_shp_ID <- 1:nrow(all_sf_files[[1]]) # Add shp ID
      GB_adm2_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[1]])))[, "GB_adm2_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm2_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm2_shp <- all_sf_files[[1]][GB_adm2_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm2_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm2_shp, focal_occ_shp))) # Build the intersection shp between adm2 shp and adm1 shape used to draw the occurrence
        GB_adm2_intersection_shp$area <- st_area(GB_adm2_intersection_shp) # Compute area of intersecting shape
        GB_adm2_intersection_shp$area <- set_units(GB_adm2_intersection_shp$area, km^2) # Convert to km2
        GB_adm2_overlap_perc <- as.numeric(GB_adm2_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- round(100 - GB_adm2_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
      
      # adm1 Uncertainty is null because adm1 was retrieved by name
      Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 0
    }
    if (matching_lvl %in% c("NE_adm1"))
    {
      # Record adm2 and GB adm1 if drawn from NE adm1
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$adm2[i] <- st_drop_geometry(new_shp_all_levels[, "adm2_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$adm1[i] <- st_drop_geometry(new_shp_all_levels[, "adm1_shapeName"])[1,1]
      
      # Record Bentity2
      Biogeographic_database_Ponerinae_NA_coords$bentity2_name[i] <- st_drop_geometry(new_shp_all_levels[, "BENTITY2_N"])[1,1]
      # Record SUBUNIT
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels[, "SUBUNIT"])[1,1]
      # Record Country Admin
      Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[i] <- st_drop_geometry(new_shp_all_levels[, "ADMIN"])[1,1]
      
      # adm2 Uncertainty is 100% - the percentage of area covered by GB_adm2 within the NE_SUBUNIT area used to draw the random coordinates
      all_sf_files[[1]]$GB_adm2_shp_ID <- 1:nrow(all_sf_files[[1]]) # Add shp ID
      GB_adm2_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[1]])))[, "GB_adm2_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm2_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm2_shp <- all_sf_files[[1]][GB_adm2_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm2_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm2_shp, focal_occ_shp))) # Build the intersection shp between adm2 shp and adm1 shape used to draw the occurrence
        GB_adm2_intersection_shp$area <- st_area(GB_adm2_intersection_shp) # Compute area of intersecting shape
        GB_adm2_intersection_shp$area <- set_units(GB_adm2_intersection_shp$area, km^2) # Convert to km2
        GB_adm2_overlap_perc <- as.numeric(GB_adm2_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- round(100 - GB_adm2_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
      
      # adm1 Uncertainty is 100% - the percentage of area covered by GB_adm1 within the Bentity2 area used to draw the random coordinates
      all_sf_files[[2]]$GB_adm1_shp_ID <- 1:nrow(all_sf_files[[2]]) # Add shp ID
      GB_adm1_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[2]])))[, "GB_adm1_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm1_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm1_shp <- all_sf_files[[2]][GB_adm1_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm1_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm1_shp, focal_occ_shp))) # Build the intersection shp between GB_adm1 shp and NE_adm1 shape used to draw the occurrence
        GB_adm1_intersection_shp$area <- st_area(GB_adm1_intersection_shp) # Compute area of intersecting shape
        GB_adm1_intersection_shp$area <- set_units(GB_adm1_intersection_shp$area, km^2) # Convert to km2
        GB_adm1_overlap_perc <- as.numeric(GB_adm1_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- round(100 - GB_adm1_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
    }
    if (matching_lvl %in% c("Bentity2"))
    {
      # Record adm2 and adm1 if drawn from Bentity2
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$adm2[i] <- st_drop_geometry(new_shp_all_levels[, "adm2_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$adm1[i] <- st_drop_geometry(new_shp_all_levels[, "adm1_shapeName"])[1,1]
      
      # Record SUBUNIT
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels["SUBUNIT"])[1,1]
      # Record Country Admin
      Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[i] <- st_drop_geometry(new_shp_all_levels["ADMIN"])[1,1]
      
      # adm2 Uncertainty is 100% - the percentage of area covered by GB_adm2 within the NE_SUBUNIT area used to draw the random coordinates
      all_sf_files[[1]]$GB_adm2_shp_ID <- 1:nrow(all_sf_files[[1]]) # Add shp ID
      GB_adm2_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[1]])))[, "GB_adm2_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm2_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm2_shp <- all_sf_files[[1]][GB_adm2_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm2_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm2_shp, focal_occ_shp))) # Build the intersection shp between adm2 shp and adm1 shape used to draw the occurrence
        GB_adm2_intersection_shp$area <- st_area(GB_adm2_intersection_shp) # Compute area of intersecting shape
        GB_adm2_intersection_shp$area <- set_units(GB_adm2_intersection_shp$area, km^2) # Convert to km2
        GB_adm2_overlap_perc <- as.numeric(GB_adm2_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- round(100 - GB_adm2_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
      
      # adm1 Uncertainty is 100% - the percentage of area covered by GB_adm1 within the Bentity2 area used to draw the random coordinates
      all_sf_files[[2]]$GB_adm1_shp_ID <- 1:nrow(all_sf_files[[2]]) # Add shp ID
      GB_adm1_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[2]])))[, "GB_adm1_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm1_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm1_shp <- all_sf_files[[2]][GB_adm1_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm1_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm1_shp, focal_occ_shp))) # Build the intersection shp between GB_adm1 shp and NE_adm1 shape used to draw the occurrence
        GB_adm1_intersection_shp$area <- st_area(GB_adm1_intersection_shp) # Compute area of intersecting shape
        GB_adm1_intersection_shp$area <- set_units(GB_adm1_intersection_shp$area, km^2) # Convert to km2
        GB_adm1_overlap_perc <- as.numeric(GB_adm1_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- round(100 - GB_adm1_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
    }
    if (matching_lvl %in% c("NE_SUBUNIT"))
    {
      # Record adm2, adm1 and Bentity2 if drawn from SUBUNIT
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_Bentity2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$adm2[i] <- st_drop_geometry(new_shp_all_levels[, "adm2_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$adm1[i] <- st_drop_geometry(new_shp_all_levels[, "adm1_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels["BENTITY2_N"])[1,1]
      
      # Record Country Admin
      Biogeographic_database_Ponerinae_NA_coords$Country_ISO3_name[i] <- st_drop_geometry(new_shp_all_levels[, "ADMIN"])[1,1]
      
      # adm2 Uncertainty is 100% - the percentage of area covered by GB_adm2 within the NE_SUBUNIT area used to draw the random coordinates
      all_sf_files[[1]]$GB_adm2_shp_ID <- 1:nrow(all_sf_files[[1]]) # Add shp ID
      GB_adm2_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[1]])))[, "GB_adm2_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm2_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm2_shp <- all_sf_files[[1]][GB_adm2_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm2_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm2_shp, focal_occ_shp))) # Build the intersection shp between adm2 shp and adm1 shape used to draw the occurrence
        GB_adm2_intersection_shp$area <- st_area(GB_adm2_intersection_shp) # Compute area of intersecting shape
        GB_adm2_intersection_shp$area <- set_units(GB_adm2_intersection_shp$area, km^2) # Convert to km2
        GB_adm2_overlap_perc <- as.numeric(GB_adm2_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- round(100 - GB_adm2_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
      
      # adm1 Uncertainty is 100% - the percentage of area covered by GB_adm1 within the Bentity2 area used to draw the random coordinates
      all_sf_files[[2]]$GB_adm1_shp_ID <- 1:nrow(all_sf_files[[2]]) # Add shp ID
      GB_adm1_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[2]])))[, "GB_adm1_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm1_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm1_shp <- all_sf_files[[2]][GB_adm1_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm1_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm1_shp, focal_occ_shp))) # Build the intersection shp between GB_adm1 shp and NE_adm1 shape used to draw the occurrence
        GB_adm1_intersection_shp$area <- st_area(GB_adm1_intersection_shp) # Compute area of intersecting shape
        GB_adm1_intersection_shp$area <- set_units(GB_adm1_intersection_shp$area, km^2) # Convert to km2
        GB_adm1_overlap_perc <- as.numeric(GB_adm1_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- round(100 - GB_adm1_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
    }
    if (matching_lvl %in% c("NE_ADMIN"))
    {
      # Record adm2, adm1, Bentity2, and SUBUNIT if drawn from Country ADMIN
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$Interpolated_Bentity2[i] <- TRUE
      Biogeographic_database_Ponerinae_NA_coords$adm2[i] <- st_drop_geometry(new_shp_all_levels[, "adm2_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$adm1[i] <- st_drop_geometry(new_shp_all_levels[, "adm1_shapeName"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$bentity2_name[i] <- st_drop_geometry(new_shp_all_levels[, "BENTITY2_N"])[1,1]
      Biogeographic_database_Ponerinae_NA_coords$SUBUNIT_name[i] <- st_drop_geometry(new_shp_all_levels[, "SUBUNIT"])[1,1]
      
      # adm2 Uncertainty is 100% - the percentage of area covered by GB_adm2 within the NE_SUBUNIT area used to draw the random coordinates
      all_sf_files[[1]]$GB_adm2_shp_ID <- 1:nrow(all_sf_files[[1]]) # Add shp ID
      GB_adm2_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[1]])))[, "GB_adm2_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm2_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm2_shp <- all_sf_files[[1]][GB_adm2_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm2_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm2_shp, focal_occ_shp))) # Build the intersection shp between adm2 shp and adm1 shape used to draw the occurrence
        GB_adm2_intersection_shp$area <- st_area(GB_adm2_intersection_shp) # Compute area of intersecting shape
        GB_adm2_intersection_shp$area <- set_units(GB_adm2_intersection_shp$area, km^2) # Convert to km2
        GB_adm2_overlap_perc <- as.numeric(GB_adm2_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[i] <- round(100 - GB_adm2_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
      
      # adm1 Uncertainty is 100% - the percentage of area covered by GB_adm1 within the Bentity2 area used to draw the random coordinates
      all_sf_files[[2]]$GB_adm1_shp_ID <- 1:nrow(all_sf_files[[2]]) # Add shp ID
      GB_adm1_shp_ID <- st_drop_geometry(suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[i, ], all_sf_files[[2]])))[, "GB_adm1_shp_ID"])[1,1] # Extract ID of the shape where the occurrence falls
      # If intersection failed, record an uncertainty of 100%
      if (is.na(GB_adm1_shp_ID)) 
      {
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- 100
      } else { # Else compute uncertainty as the ratio of areas
        GB_adm1_shp <- all_sf_files[[2]][GB_adm1_shp_ID, ] # Extract the shape where the occurrence falls
        GB_adm1_intersection_shp <- suppressMessages(suppressWarnings(sf::st_intersection(GB_adm1_shp, focal_occ_shp))) # Build the intersection shp between GB_adm1 shp and NE_adm1 shape used to draw the occurrence
        GB_adm1_intersection_shp$area <- st_area(GB_adm1_intersection_shp) # Compute area of intersecting shape
        GB_adm1_intersection_shp$area <- set_units(GB_adm1_intersection_shp$area, km^2) # Convert to km2
        GB_adm1_overlap_perc <- as.numeric(GB_adm1_intersection_shp$area / focal_occ_shp$area * 100) # Compute % of area overlap
        Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[i] <- round(100 - GB_adm1_overlap_perc, 1) # Turn into % of area non-overlap = % uncertainty
      }
    }
  }

  # Print progress
  if ( i %% 100 == 0) { cat(paste0(Sys.time(), " - Occurrence n°",i, " / ", nrow(Biogeographic_database_Ponerinae_NA_coords),"\n")) }
}
sf_use_s2(TRUE)

## Save database of NA coordinates
# saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")


### 3.6.3/ Manual correction of entries with wrong matching due to border effects ####

# WARNING: Entry n°4668 has an error in matching shp names. Initial NE_adm1 = Guangzhou Province. New NE_adm1 = NA
# WARNING: Entry n°4669 has an error in matching shp names. Initial NE_adm1 = Guangzhou Province. New NE_adm1 = Guangdong
# WARNING: Entry n°5801 has an error in matching shp names. Initial Bentity2 = San Andres Providencia and Santa Catalina. New Bentity2 = NA
# WARNING: Entry n°5913 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°5922 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°5994 has an error in matching shp names. Initial NE_adm1 = Lord Howe Island. New NE_adm1 = NA
# WARNING: Entry n°6037 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°6054 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°6067 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°6523 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°6530 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°7166 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°7170 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°7173 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°7178 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°7347 has an error in matching shp names. Initial NE_adm1 = Macau. New NE_adm1 = NA
# WARNING: Entry n°7518 has an error in matching shp names. Initial Bentity2 = Shikoku. New Bentity2 = NA
# WARNING: Entry n°8127 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°8222 has an error in matching shp names. Initial NE_adm1 = Lord Howe Island. New NE_adm1 = NA
# WARNING: Entry n°8315 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°8574 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA
# WARNING: Entry n°8597 has an error in matching shp names. Initial GB_adm1 = Incheon. New GB_adm1 = NA


sf_use_s2(FALSE)
## WARNING: Entry n°1502 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[1502, ]
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Locality = Zengcheng
# Correct adm1 = Guangzhou Province
# Correct adm2 = Zengchengshi

# A priori correction
# Locality "Zengcheng" in China is adm2 Zengchengshi and adm1 Guangzhou Province
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Guangzhou Province"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[1502] <- 23.300
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[1502] <- 113.747
Biogeographic_database_Ponerinae_NA_coords$adm1[1502] <- "Guangzhou Province"
Biogeographic_database_Ponerinae_NA_coords$adm2[1502] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[1502] <- "Guangdong"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[1502] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[1502] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[1502] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[1502] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[1502] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[1502] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[1502] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[1502] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[1502] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[1502] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[1502] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[1502] <- 0

## WARNING: Entry n°2505 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[2505, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, only adjust the uncertainty
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[2505] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[2505] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[2505] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[2505] <- 0

## WARNING: Entry n°2566 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[2566, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# A priori correction
# Locality "Rosis" in France is adm2 Hérault
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Rosis", replace = F)] <- "Hérault"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[2566] <- 43.624
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[2566] <- 3.0048
Biogeographic_database_Ponerinae_NA_coords$adm2[2566] <- "Hérault"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[2566] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[2566] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[2566] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Hérault", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[2566] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Hérault")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[2566] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[2566] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[2566] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[2566] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[2566] <- 0

## WARNING: Entry n°2607 has no match => Same location than 2566
error_entry <- Biogeographic_database_Ponerinae_NA_coords[2607, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# A priori correction
# Locality "Rosis" in France is adm2 Hérault
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Rosis", replace = F)] <- "Hérault"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[2607] <- 43.624
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[2607] <- 3.0048
Biogeographic_database_Ponerinae_NA_coords$adm2[2607] <- "Hérault"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[2607] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[2607] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[2607] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Hérault", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[2607] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Hérault")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[2607] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[2607] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[2607] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[2607] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[2607] <- 0

## WARNING: Entry n°3277 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[3277, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[3277] <- 28.533
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[3277] <- -80.803
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[3277] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[3277] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[3277] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[3277] <- 0

## WARNING: Entry n°3284 has no match => Same locality than 3277
error_entry <- Biogeographic_database_Ponerinae_NA_coords[3284, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[3284] <- 28.533
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[3284] <- -80.803
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[3284] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[3284] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[3284] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[3284] <- 0

## WARNING: Entry n°3285 has no match => Same locality than 3277 and 3284
error_entry <- Biogeographic_database_Ponerinae_NA_coords[3285, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[3285] <- 28.533
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[3285] <- -80.803
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[3285] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[3285] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[3285] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[3285] <- 0

## WARNING: Entry n°3350 has no match => Same GB_adm2 than 3277, 3284 and 3285 
error_entry <- Biogeographic_database_Ponerinae_NA_coords[3350, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[3350] <- 28.613
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[3350] <- -80.818
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[3350] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[3350] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[3350] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[3350] <- 0

## WARNING: Entry n°3965 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[3965, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[3965] <- 27.5148
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[3965] <- -97.856
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[3965] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[3965] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[3965] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[3965] <- 0

## WARNING: Entry n°4367 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4367, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into a lake. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4367] <- 28.696
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4367] <- -95.9656
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4367] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4367] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4367] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4367] <- 0

## WARNING: Entry n°4405 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4405, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# A priori correction
# Locality "Talant" in France is adm2 Côte-d'Or and adm1 Bourgogne-Franche-Comté
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Talant", replace = F)] <- "Côte-d'Or"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Talant", replace = F)] <- "Bourgogne-Franche-Comté"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4405] <- 47.338
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4405] <- 4.999
Biogeographic_database_Ponerinae_NA_coords$adm2[4405] <- "Côte-d'Or"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[4405] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4405] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4405] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Hérault", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4405] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Hérault")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4405] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4405] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4405] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4405] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4405] <- 0


## WARNING: Entry n°4478 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4478, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# A priori correction
# Locality "AEC Savannah River Plant" in USA is adm2 Aiken
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "AEC Savannah River Plant", replace = F)] <- "Aiken"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4478] <- 33.559
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4478] <- -81.724
Biogeographic_database_Ponerinae_NA_coords$adm2[4478] <- "Aiken"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[4478] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4478] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4478] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Aiken", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4478] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Aiken")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4478] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4478] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4478] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4478] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4478] <- 0

## WARNING: Entry n°4668 has no match => Same as 1502
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4668, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Locality = Zengcheng
# Correct adm1 = Guangzhou Province
# Correct adm2 = Zengchengshi

# A priori correction
# Locality "Zengcheng" in China is adm2 Zengchengshi and adm1 Guangzhou Province
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Guangzhou Province"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4668] <- 23.300
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4668] <- 113.747
Biogeographic_database_Ponerinae_NA_coords$adm1[4668] <- "Guangzhou Province"
Biogeographic_database_Ponerinae_NA_coords$adm2[4668] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[4668] <- "Guangdong"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[4668] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4668] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[4668] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[4668] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[4668] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4668] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4668] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4668] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4668] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4668] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4668] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4668] <- 0

## WARNING: Entry n°4669 has no match => Same as 1502 & 4668
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4669, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Locality = Zengcheng
# Correct adm1 = Guangzhou Province
# Correct adm2 = Zengchengshi

# A priori correction
# Locality "Zengcheng" in China is adm2 Zengchengshi and adm1 Guangzhou Province
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Zengcheng", replace = F)] <- "Guangzhou Province"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4669] <- 23.300
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4669] <- 113.747
Biogeographic_database_Ponerinae_NA_coords$adm1[4669] <- "Guangzhou Province"
Biogeographic_database_Ponerinae_NA_coords$adm2[4669] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[4669] <- "Guangdong"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[4669] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4669] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[4669] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[4669] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[4669] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4669] <- set_units(st_area(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi", "adm2_shapeName"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4669] <- which(geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4669] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4669] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4669] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4669] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4669] <- 0


## WARNING: Entry n°4879 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4879, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Locality = Hwy 70, 6 miles W Hillsboro
# Correct adm2 = Orange
# Multiple matches for adm2

# A priori correction
# Locality "Hwy 70, 6 miles W Hillsboro" in USA is adm2 Orange
Biogeographic_database_Ponerinae_NA_coords$adm2[replace_na(data = Biogeographic_database_Ponerinae_NA_coords$Locality == "Hwy 70, 6 miles W Hillsboro", replace = F)] <- "Orange"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4879] <- 36.0815
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4879] <- -79.1874
Biogeographic_database_Ponerinae_NA_coords$adm2[4879] <- "Zengchengshi"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[4879] <- "Orange"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[4879] <- "GB_adm2"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4879] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4879] <- set_units(st_area(geoBoundaries_adm2_sf[959,]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4879] <- 138595
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4879] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4879] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4879] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4879] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4879] <- 0

## WARNING: Entry n°4907 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4907, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into the swamp. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4907] <- 29.589
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4907] <- -83.377
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4907] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4907] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4907] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4907] <- 0

## WARNING: Entry n°4929 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4929, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# adm2 and adm1 are correct, but occurrence fall into the river. Adjust coordinates and uncertainty
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4929] <- 45.6375
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4929] <- -84.5007
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4929] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4929] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4929] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4929] <- 0

## WARNING: Entry n°4979 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4979, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# SUBUNIT is correct, but occurrence fall into the Ocean. Adjust coordinates and uncertainty. Provide random adm1 and adm2
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4979] <- 13.3174
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4979] <- -61.1307
Biogeographic_database_Ponerinae_NA_coords$adm1[4979] <- "Charlotte"
Biogeographic_database_Ponerinae_NA_coords$adm2[4979] <- "Charlotte"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[4979] <- "Lesser Antilles"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4979] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[4979] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[4979] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[4979] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4979] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "Saint Vincent and the Grenadines", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4979] <- which(NE_subunits_sf$SUBUNIT == "Saint Vincent and the Grenadines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4979] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4979] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4979] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4979] <- 50
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4979] <- 50


## WARNING: Entry n°4980 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[4980, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# SUBUNIT is correct, but occurrence fall into the Ocean. Adjust coordinates and uncertainty. Provide random adm1 and adm2
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[4980] <- 13.1359
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[4980] <- -59.4635
Biogeographic_database_Ponerinae_NA_coords$adm1[4980] <- "Saint Philip"
Biogeographic_database_Ponerinae_NA_coords$adm2[4980] <- "Saint Philip"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[4980] <- "Lesser Antilles"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[4980] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[4980] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[4980] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[4980] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[4980] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "Barbados", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[4980] <- which(NE_subunits_sf$SUBUNIT == "Barbados")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[4980] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[4980] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[4980] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[4980] <- 90
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[4980] <- 90

## WARNING: Entry n°5035 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5035, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct, but occurrence fall into the Ocean. Adjust coordinates and uncertainty. Provide random adm1 and adm2
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[5035] <- 9.8882
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[5035] <- 126.0363
Biogeographic_database_Ponerinae_NA_coords$adm1[5035] <- "Caraga"
Biogeographic_database_Ponerinae_NA_coords$adm2[5035] <- "Surigao del Norte"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5035] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5035] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5035] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5035] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5035] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5035] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5035] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5035] <- 99
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5035] <- 90

## WARNING: Entry n°5051 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5051, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# American Samoa (SUBUNIT) is correct. No adm1 and adm2 available. Provide Bentity2
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[5051] <- "Samoan Islands"
Biogeographic_database_Ponerinae_NA_coords$matching_bentity2[5051] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5051] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "American Samoa", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5051] <- which(NE_subunits_sf$SUBUNIT == "American Samoa")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5051] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5051] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5051] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5051] <- NA
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5051] <- NA

## WARNING: Entry n°5053 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5053, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# American Samoa (SUBUNIT) is correct. No adm1 and adm2 available. Provide Bentity2
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[5053] <- "Samoan Islands"
Biogeographic_database_Ponerinae_NA_coords$matching_bentity2[5053] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5053] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "American Samoa", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5053] <- which(NE_subunits_sf$SUBUNIT == "American Samoa")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5053] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5053] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5053] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5053] <- NA
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5053] <- NA

## WARNING: Entry n°5056 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5056, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# American Samoa (SUBUNIT) is correct. No adm1 and adm2 available. Provide Bentity2
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[5056] <- "Samoan Islands"
Biogeographic_database_Ponerinae_NA_coords$matching_bentity2[5056] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5056] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "American Samoa", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5056] <- which(NE_subunits_sf$SUBUNIT == "American Samoa")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5056] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5056] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5056] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5056] <- NA
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5056] <- NA

## WARNING: Entry n°5072 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5072, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct. Provide random adm1 and adm2 and estimate uncertainty.
Biogeographic_database_Ponerinae_NA_coords$adm1[5072] <- "Central Visayas"
Biogeographic_database_Ponerinae_NA_coords$adm2[5072] <- "Bohol"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5072] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5072] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5072] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5072] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5072] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5072] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5072] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5072] <- 98
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5072] <- 90

## WARNING: Entry n°5078 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5078, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct, but Locality is the true adm2. Adjust coordinates and uncertainty. Provide adm1 and adm2.
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[5078] <- 8.9633
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[5078] <- 126.120
Biogeographic_database_Ponerinae_NA_coords$adm1[5078] <- "Caraga"
Biogeographic_database_Ponerinae_NA_coords$adm2[5078] <- "Surigao del Sur"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5078] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5078] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5078] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5078] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5078] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5078] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5078] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5078] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5078] <- 0

## WARNING: Entry n°5092 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5092, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# SUBUNIT is correct, but occurrence fall into the Ocean. Adjust coordinates and uncertainty. Provide random adm1 and adm2
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[5092] <- 13.246
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[5092] <- -59.5698
Biogeographic_database_Ponerinae_NA_coords$adm1[5092] <- "Saint Andrew"
Biogeographic_database_Ponerinae_NA_coords$adm2[5092] <- "Saint Andrew"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[5092] <- "Lesser Antilles"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5092] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5092] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_NE_adm1[5092] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[5092] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5092] <- set_units(st_area(NE_subunits_sf[NE_subunits_sf$SUBUNIT == "Barbados", "SUBUNIT"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5092] <- which(NE_subunits_sf$SUBUNIT == "Barbados")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5092] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5092] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5092] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5092] <- 90
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5092] <- 90

## WARNING: Entry n°5095 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5095, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Christmas Island (Bentity2) is correct. Provide adm1 (Other Territories) and adm2 (Christmas Island).
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[5092] <- -10.4722
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[5092] <- 105.6566
Biogeographic_database_Ponerinae_NA_coords$adm1[5095] <- "Other Territories"
Biogeographic_database_Ponerinae_NA_coords$adm2[5095] <- "Christmas Island"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5095] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5095] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5095] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Christmas Island", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5095] <- which(Bentity2_sf$BENTITY2_N == "Christmas Island")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5095] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5095] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5095] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5095] <- 0
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5095] <- 0

## WARNING: Entry n°5129 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5129, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct. Provide random adm1 and adm2 and estimate uncertainty.
Biogeographic_database_Ponerinae_NA_coords$adm1[5129] <- "Caraga"
Biogeographic_database_Ponerinae_NA_coords$adm2[5129] <- "Agusan del Norte"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5129] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5129] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5129] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5129] <- 95
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5129] <- 90

## WARNING: Entry n°5152 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5152, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct. Provide random adm1 and adm2 and estimate uncertainty.
Biogeographic_database_Ponerinae_NA_coords$adm1[5129] <- "Ilocos Region"
Biogeographic_database_Ponerinae_NA_coords$adm2[5129] <- "Pangasinan"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5129] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5129] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5129] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5129] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5129] <- 98
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5129] <- 90

## WARNING: Entry n°5153 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5153, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))

# Philippines (Bentity2) is correct. Provide random adm1 and adm2 and estimate uncertainty.
Biogeographic_database_Ponerinae_NA_coords$adm1[5153] <- "Mimaropa"
Biogeographic_database_Ponerinae_NA_coords$adm2[5153] <- "Palawan"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[5153] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm1[5153] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[5153] <- set_units(st_area(Bentity2_sf[Bentity2_sf$BENTITY2_N == "Philippines", "BENTITY2_N"]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[5153] <- which(Bentity2_sf$BENTITY2_N == "Philippines")
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[5153] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[5153] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[5153] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[5153] <- 98
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[5153] <- 90

## WARNING: Entry n°5216 has no match
error_entry <- Biogeographic_database_Ponerinae_NA_coords[5216, ]
print(error_entry)
new_shp_all_levels <- suppressMessages(suppressWarnings(sf::st_intersection(error_entry, Intersect_all_shp_all_levels)))


## WARNING: Entry n°5221 has no match
## WARNING: Entry n°5222 has no match
## WARNING: Entry n°5224 has no match
## WARNING: Entry n°5241 has no match
## WARNING: Entry n°5242 has no match
## WARNING: Entry n°5261 has no match

# WARNING: Entry n°5303 has no match
# WARNING: Entry n°5400 has no match

# WARNING: Entry n°5405 has no match
# WARNING: Entry n°5415 has no match
# WARNING: Entry n°5422 has no match
# WARNING: Entry n°5424 has no match
# WARNING: Entry n°5456 has no match
# WARNING: Entry n°5481 has no match
# WARNING: Entry n°5495 has no match
# WARNING: Entry n°5496 has no match

# WARNING: Entry n°5537 has no match
# WARNING: Entry n°5559 has no match
# WARNING: Entry n°5569 has no match
# WARNING: Entry n°5573 has no match
# WARNING: Entry n°5574 has no match

# WARNING: Entry n°5601 has no match
# WARNING: Entry n°5605 has no match
# WARNING: Entry n°5642 has no match
# WARNING: Entry n°5654 has no match
# WARNING: Entry n°5665 has no match
# WARNING: Entry n°5682 has no match
# WARNING: Entry n°5683 has no match

# WARNING: Entry n°5706 has no match
# WARNING: Entry n°5708 has no match
# WARNING: Entry n°5709 has no match
# WARNING: Entry n°5711 has no match
# WARNING: Entry n°5714 has no match
# WARNING: Entry n°5715 has no match
# WARNING: Entry n°5716 has no match
# WARNING: Entry n°5717 has no match
# WARNING: Entry n°5720 has no match
# WARNING: Entry n°5721 has no match
# WARNING: Entry n°5725 has no match
# WARNING: Entry n°5726 has no match
# WARNING: Entry n°5742 has no match
# WARNING: Entry n°5764 has no match

# WARNING: Entry n°5801 has no match
# WARNING: Entry n°5802 has no match
# WARNING: Entry n°5803 has no match
# WARNING: Entry n°5810 has no match
# WARNING: Entry n°5814 has no match
# WARNING: Entry n°5833 has no match
# WARNING: Entry n°5845 has no match
# WARNING: Entry n°5846 has no match
# WARNING: Entry n°5853 has no match
# WARNING: Entry n°5863 has no match
# WARNING: Entry n°5872 has no match
# WARNING: Entry n°5878 has no match

# WARNING: Entry n°5913 has no match
# WARNING: Entry n°5922 has no match
# WARNING: Entry n°5935 has no match
# WARNING: Entry n°5946 has no match
# WARNING: Entry n°5966 has no match
# WARNING: Entry n°5967 has no match
# WARNING: Entry n°5988 has no match
# WARNING: Entry n°5994 has no match

# WARNING: Entry n°6026 has no match
# WARNING: Entry n°6037 has no match
# WARNING: Entry n°6053 has no match
# WARNING: Entry n°6054 has no match
# WARNING: Entry n°6059 has no match
# WARNING: Entry n°6063 has no match
# WARNING: Entry n°6067 has no match
# WARNING: Entry n°6078 has no match
# WARNING: Entry n°6100 has no match

# WARNING: Entry n°6110 has no match
# WARNING: Entry n°6117 has no match
# WARNING: Entry n°6118 has no match
# WARNING: Entry n°6129 has no match
# WARNING: Entry n°6187 has no match

# WARNING: Entry n°6242 has no match
# WARNING: Entry n°6249 has no match
# WARNING: Entry n°6289 has no match

# WARNING: Entry n°6309 has no match

# WARNING: Entry n°6443 has no match

# WARNING: Entry n°6503 has no match
# WARNING: Entry n°6511 has no match
# WARNING: Entry n°6520 has no match
# WARNING: Entry n°6523 has no match
# WARNING: Entry n°6527 has no match
# WARNING: Entry n°6530 has no match
# WARNING: Entry n°6543 has no match
# WARNING: Entry n°6571 has no match
# WARNING: Entry n°6588 has no match
# WARNING: Entry n°6593 has no match

# WARNING: Entry n°6661 has no match
# WARNING: Entry n°6670 has no match
# WARNING: Entry n°6681 has no match

# WARNING: Entry n°6709 has no match
# WARNING: Entry n°6711 has no match
# WARNING: Entry n°6714 has no match
# WARNING: Entry n°6729 has no match
# WARNING: Entry n°6737 has no match
# WARNING: Entry n°6740 has no match
# WARNING: Entry n°6752 has no match

# WARNING: Entry n°6804 has no match
# WARNING: Entry n°6815 has no match
# WARNING: Entry n°6843 has no match
# WARNING: Entry n°6863 has no match
# WARNING: Entry n°6887 has no match

# WARNING: Entry n°6935 has no match

# WARNING: Entry n°7032 has no match
# WARNING: Entry n°7060 has no match
# WARNING: Entry n°7076 has no match

# WARNING: Entry n°7121 has no match
# WARNING: Entry n°7122 has no match
# WARNING: Entry n°7134 has no match
# WARNING: Entry n°7146 has no match
# WARNING: Entry n°7151 has no match
# WARNING: Entry n°7153 has no match
# WARNING: Entry n°7164 has no match
# WARNING: Entry n°7166 has no match
# WARNING: Entry n°7169 has no match
# WARNING: Entry n°7170 has no match
# WARNING: Entry n°7173 has no match
# WARNING: Entry n°7178 has no match

# WARNING: Entry n°7239 has no match
# WARNING: Entry n°7250 has no match
# WARNING: Entry n°7262 has no match
# WARNING: Entry n°7281 has no match

# WARNING: Entry n°7323 has no match
# WARNING: Entry n°7326 has no match
# WARNING: Entry n°7341 has no match
# WARNING: Entry n°7346 has no match
# WARNING: Entry n°7347 has no match
# WARNING: Entry n°7351 has no match
# WARNING: Entry n°7361 has no match
# WARNING: Entry n°7363 has no match
# WARNING: Entry n°7382 has no match
# WARNING: Entry n°7398 has no match

# WARNING: Entry n°7439 has no match
# WARNING: Entry n°7447 has no match
# WARNING: Entry n°7451 has no match
# WARNING: Entry n°7469 has no match
# WARNING: Entry n°7485 has no match

# WARNING: Entry n°7504 has no match
# WARNING: Entry n°7507 has no match
# WARNING: Entry n°7518 has no match
# WARNING: Entry n°7542 has no match

# WARNING: Entry n°7649 has no match
# WARNING: Entry n°7650 has no match
# WARNING: Entry n°7693 has no match

# WARNING: Entry n°7706 has no match
# WARNING: Entry n°7788 has no match

# WARNING: Entry n°7841 has no match
# WARNING: Entry n°7855 has no match
# WARNING: Entry n°7863 has no match
# WARNING: Entry n°7883 has no match
# WARNING: Entry n°7895 has no match

# WARNING: Entry n°7934 has no match
# WARNING: Entry n°7939 has no match
# WARNING: Entry n°7940 has no match
# WARNING: Entry n°7957 has no match

# WARNING: Entry n°8091 has no match

# WARNING: Entry n°8127 has no match
# WARNING: Entry n°8128 has no match
# WARNING: Entry n°8138 has no match
# WARNING: Entry n°8140 has no match
# WARNING: Entry n°8157 has no match

# WARNING: Entry n°8222 has no match
# WARNING: Entry n°8226 has no match
# WARNING: Entry n°8259 has no match
# WARNING: Entry n°8272 has no match

# WARNING: Entry n°8308 has no match
# WARNING: Entry n°8315 has no match
# WARNING: Entry n°8322 has no match
# WARNING: Entry n°8324 has no match
# WARNING: Entry n°8328 has no match
# WARNING: Entry n°8334 has no match
# WARNING: Entry n°8336 has no match
# WARNING: Entry n°8347 has no match
# WARNING: Entry n°8348 has no match

# WARNING: Entry n°8405 has no match
# WARNING: Entry n°8407 has no match
# WARNING: Entry n°8473 has no match

# WARNING: Entry n°8512 has no match
# WARNING: Entry n°8521 has no match
# WARNING: Entry n°8530 has no match
# WARNING: Entry n°8574 has no match
# WARNING: Entry n°8584 has no match
# WARNING: Entry n°8585 has no match
# WARNING: Entry n°8597 has no match
# WARNING: Entry n°8600 has no match

# WARNING: Entry n°8614 has no match
# WARNING: Entry n°8626 has no match
# WARNING: Entry n°8629 has no match
# WARNING: Entry n°8632 has no match
# WARNING: Entry n°8643 has no match
# WARNING: Entry n°8648 has no match
# WARNING: Entry n°8655 has no match
# WARNING: Entry n°8660 has no match
# WARNING: Entry n°8671 has no match
# WARNING: Entry n°8673 has no match
# WARNING: Entry n°8678 has no match
# WARNING: Entry n°8682 has no match
# WARNING: Entry n°8685 has no match
# WARNING: Entry n°8689 has no match
# WARNING: Entry n°8697 has no match
# WARNING: Entry n°8698 has no match
# WARNING: Entry n°8699 has no match
# WARNING: Entry n°8700 has no match

# WARNING: Entry n°8701 has no match
# WARNING: Entry n°8702 has no match
# WARNING: Entry n°8703 has no match
# WARNING: Entry n°8704 has no match
# WARNING: Entry n°8714 has no match
# WARNING: Entry n°8716 has no match
# WARNING: Entry n°8728 has no match
# WARNING: Entry n°8743 has no match
# WARNING: Entry n°8764 has no match
# WARNING: Entry n°8765 has no match
# WARNING: Entry n°8767 has no match
# WARNING: Entry n°8769 has no match
# WARNING: Entry n°8774 has no match
# WARNING: Entry n°8776 has no match
# WARNING: Entry n°8777 has no match
# WARNING: Entry n°8778 has no match
# WARNING: Entry n°8783 has no match
# 
# WARNING: Entry n°8814 has no match
# WARNING: Entry n°8817 has no match
# WARNING: Entry n°8818 has no match
# WARNING: Entry n°8820 has no match
# WARNING: Entry n°8821 has no match
# WARNING: Entry n°8825 has no match
# WARNING: Entry n°8827 has no match
# WARNING: Entry n°8829 has no match
# WARNING: Entry n°8831 has no match
# WARNING: Entry n°8833 has no match
# WARNING: Entry n°8834 has no match
# WARNING: Entry n°8835 has no match
# WARNING: Entry n°8843 has no match
# WARNING: Entry n°8846 has no match

# WARNING: Entry n°8849 has an error

# WARNING: Entry n°50094 has an error
Kaliningrad

# Locality = locality unknown in Samland Peninsula
# Correct adm1 = Kaliningrad (Not Kaliningradskaya Oblast as initially)

# A priori correction
# Locality "locality unknown in Samland Peninsula" in USA is adm1 = Kaliningrad
Biogeographic_database_Ponerinae_NA_coords$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_NA_coords$Locality, pattern = "Samland"), replace = F)] <- "Kaliningrad"

# A posteriori correction
Biogeographic_database_Ponerinae_NA_coords$Latitude_dec[50094] <- 54.657
Biogeographic_database_Ponerinae_NA_coords$Longitude_dec[50094] <- 21.171
Biogeographic_database_Ponerinae_NA_coords$adm2[50094] <- "Gvardeysky District"
Biogeographic_database_Ponerinae_NA_coords$bentity2_name[50094] <- "Kaliningradskaya oblast"
Biogeographic_database_Ponerinae_NA_coords$matching_shp[50094] <- "GB_adm1"
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[50094] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$matching_GB_adm2[50094] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$matching_Bentity2[50094] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_area[50094] <- set_units(st_area(geoBoundaries_adm1_sf[2307,]), km^2)
Biogeographic_database_Ponerinae_NA_coords$matching_shp_ID[50094] <- 2307
Biogeographic_database_Ponerinae_NA_coords$multiple_matches[50094] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm2[50094] <- TRUE
Biogeographic_database_Ponerinae_NA_coords$Interpolated_adm1[50094] <- FALSE
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm2[50094] <- 90
Biogeographic_database_Ponerinae_NA_coords$Uncertainty_adm1[50094] <- 0

sf_use_s2(TRUE)



# View(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeGroup == "RUS", ])
# View(geoBoundaries_adm1_sf[geoBoundaries_adm1_sf$adm1_shapeGroup == "RUS", ])
# 
# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName %in% c("Zengchengshi", "Guangzhoushi"), "adm2_shapeName"])
# 
# plot(geoBoundaries_adm2_sf[geoBoundaries_adm2_sf$adm2_shapeName == "Zengchengshi", "adm2_shapeName"])
# plot(geoBoundaries_adm1_sf[geoBoundaries_adm1_sf$adm1_shapeName == "Other Territories", "adm1_shapeName"])




## Check that all occurrences have an adm1/adm2/Bentity2


## Save database of NA coordinates
# saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")


##### 4/ Merge randomly drawn NA coordinates with curated coordinates ####


### 4.1/ Find intersecting polygons for all entries based on coordinates

# Record when info is updated based on coordinates vs. when match what was provided ????


# ## Find intersecting polygons for all entries
# sf_use_s2(FALSE)
# out <- sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords, geoBoundaries_Countries_sf)
# sf_use_s2(TRUE)

## Find intersecting polygons by batches of 100 entries (to keep track of progress)
out <- data.frame() # Store output
batch_size <- 100   # Number of entries per batch
nb_batches <- (nrow(Biogeographic_database_Ponerinae_NA_coords) %/% batch_size) + is.numeric((nrow(Biogeographic_database_Ponerinae_NA_coords) %% batch_size) > 0)
seq_i <- seq(from = 1, to = batch_size * nb_batches + 1, by = batch_size)

for (i in seq_i)
{
  # i <- 1
  
  # Get indices of batch entries
  batch_i <- seq(from = i, by = 1, length.out = batch_size)
  
  # Adjust indices for the last batch
  if (max(batch_i) > nrow(Biogeographic_database_Ponerinae_NA_coords)) { batch_i <- min(batch_i):nrow(Biogeographic_database_Ponerinae_NA_coords) }
  
  # Find intersecting polygons for Countries/adm1/adm2
  out_batch <- suppressMessages(suppressWarnings(sf::st_intersection(Biogeographic_database_Ponerinae_NA_coords[batch_i, ], geoBoundaries_Countries_sf)))
  out_batch <- suppressMessages(suppressWarnings(sf::st_intersection(out_batch, geoBoundaries_adm1_sf)))
  out_batch <- suppressMessages(suppressWarnings(sf::st_intersection(out_batch, geoBoundaries_adm2_sf)))
  
  # Merge batch output with global output
  out <- rbind(out, out_batch)
  
  # Print progress
  cat(paste0(Sys.time(), " - Occurrence n°",max(batch_i), " / ", nrow(Biogeographic_database_Ponerinae_NA_coords),"\n"))
  
}

nrow(out)
nrow(Biogeographic_database_Ponerinae_NA_coords)


##### 5/ Flag and remove duplicates #####

## Create a database without duplicates. Keep AntWeb over GABI when possible. Run taxa-summary stat on without duplicates data


##### 6/ Detect outliers ####

### Rerun distance based and Bioregion-based outliers detection using also the randomly drawn occurrences

### Find categories of taxa (can have multiple categories, but in this case, use the higher one)



### Plot maps per taxa and display them in proper category

# Use shapes for type of records: true vs. interpolated
# use size and transparency for uncertainty of record => Bigger and more transparent ~ Uncertainty_area
# Use colors for outliers: green = ok ; orange = distance ; red = bioregions


##### 7/ Manually curate biogeographic errors based on visual checking of maps ####

# Anochetus_siphneus => Anochetus_katonae_nr2. Also need to change all specimens records from Zambia!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "katonae_nr2"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "Anochetus_katonae_nr2"
# Anochetus_siphneus => Anochetus_katonae_nr3. Rename specimen from DRC for a new unique name (Anochetus_katonae_nr3) not to use in our database!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "katonae_nr3"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "Anochetus_katonae_nr3"

# Odontomachus_brunneus => Odontomachus_brunneus_nr. Also need to change all specimens records from South America/Outside South East USA!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "brunneus_nr"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "Odontomachus_brunneus_nr"



##### 8/ Update records from initial databases too (GABI and AntWeb), not only the merged one


## Also flag changes and provide the list of updates to Brian when the changes are not regarding stuff that was interpolated


##### 9/ Map cleaned dataset vs. original one ####

world.inp <- map_data("world")

ggplot() + geom_map(data = world.inp, map = world.inp, 
                    aes(x = long, y = lat, map_id = region), fill = "grey80") + 
  xlim(min(dat$decimalLongitude, na.rm = T), max(dat$decimalLongitude, na.rm = T)) + 
  ylim(min(dat$decimalLatitude, na.rm = T), max(dat$decimalLatitude, na.rm = T)) + 
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 1) + 
  geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkgreen", size = 1) + 
  coord_fixed() + theme_bw() + theme(axis.title = element_blank())

map_plot(dat_cl)


### Flag for differences between initial coordinates in GABI and AntWeb, and current coordinates
### Flag for randomly drawn coordinates within country/adm1/2 = "New"

# Also, plot before/after with former coordinates (and showing missing as grey shadows)



##### 10/ Extract specimen database for taxa-level information ####

# Add nb of biogeographic records
# Presence/Absence in each Bioregion per taxa
# Create multiple Bioregion schemes

##### 11/ Create the Genus-level summary table #####

# Create a summary table at Genus-level of nb of valid species, nb of valid species included in our phylogeny, nb of morphospecies in our phylogeny, other (not curated) morphospecies from AntWeb
# Add nb of biogeographic records




##### BONUS/ Plot interactive maps with mapview ####

library(mapview)

?mapview
?mapviewGetOption

## Load output from CoordinateCleaner
Biogeographic_database_Ponerinae_flags <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")
Biogeographic_database_Ponerinae_flags_no_curation <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags_no_curation.rds")

Taxa_in_phylo <- unique(Biogeographic_database_Ponerinae_flags$Current_name[Biogeographic_database_Ponerinae_flags$In_phylogeny])
Taxa_with_valid_name <- unique(Biogeographic_database_Ponerinae_flags$Current_name[Biogeographic_database_Ponerinae_flags$Status_AntCat == "valid"])

Taxa_for_analyses <- union(Taxa_in_phylo, Taxa_with_valid_name)
Taxa_valid_not_in_phylo <- setdiff(Taxa_for_analyses, Taxa_in_phylo)

## B.1/ Prepare color palette ####

# Extract number of categories to plot
nb_taxa <- length(unique(Biogeographic_database_Ponerinae_flags$Current_name[]))
nb_flags <- length(unique(Biogeographic_database_Ponerinae_flags$flag_type[]))
# nb_flags <- 7

# Define palette function
pal <- mapviewPalette("mapviewSpectralColors")

RColorBrewer::display.brewer.all()

# pal <- brewer.pal(n = nb_taxa, name = "Spectral")
pal <- brewer.pal(n = nb_flags, name = "Spectral")
# pal_4 <- pal[c(1, 4, 5, 6)]
# pal_2 <- pal[c(5, 6)]
  
plot(Biogeographic_database_Ponerinae_flags[, ]["flag_type"])

## B.2/ Plot per species ####

selected_taxa <- c("Hypoponera_punctatissima", "Hypoponera_ergatandria", "Hypoponera_ug02", "Hypoponera_ug15", "Hypoponera_afrc-za06", "Hypoponera_casc-mz11", "Hypoponera_mg116")
nb_taxa <- length(selected_taxa)
pal <- brewer.pal(n = nb_taxa, name = "Spectral")

### With color per species (2530 species !)
interactive_map_taxa <- mapview::mapview(
                                    # x = Biogeographic_database_Ponerinae_flags[,], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                    x = Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Current_name %in% selected_taxa ,], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                    zcol = "Current_name",
                                    cex = 7,
                                    burst = T, legend = T,
                                    # col.regions = pal(nb_taxa), # Color to use to plot spatial units in the spatial object
                                    col.regions = pal, # Color to use to plot spatial units in the spatial object
                                    layer.name = "Occurrences", # To specify the legend name of this layer
                                    # map.types = mapviewGetOption("basemaps")
                                    map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "Stadia.StamenWatercolor")
                                    # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stadia.StamenWatercolor") # To specify basemaps
                                    # ... # and many more, depending of the class of the spatial object
)

interactive_map_taxa

# To save the interactive map from mapview or leaflet into an interactive .hmtl or a static .png/.jpg/.pdf map
?mapshot

mapshot(x = interactive_map_taxa,
        url = "./maps/Taxa_occurrence_Control_maps/Interactive_map/interactive_map_occurrences_control_Hypoponera_punctatissima_vs_ergatandria.html", # To save an interactive .html map
        remove_controls = NULL) # c("zoomControl", "layersControl", "homeButton", "scaleBar"))  # To specify which control buttons to remove (or not) from the output


## B.3/ Plot per flag types ####

### With color per flag types (11 flag types)
interactive_map_flags <- mapview::mapview(x = Biogeographic_database_Ponerinae_flags[ ,], # sf/tmap/sp/raster/dataframe(with coordinates)... object to plot interactively
                                    zcol = "flag_type",
                                    cex = 5,
                                    burst = T, legend = T,
                                    # col.regions = pal(nb_flags), # Color to use to plot spatial units in the spatial object
                                    col.regions = pal, # Color to use to plot spatial units in the spatial object
                                    # col.regions = pal_2, # Color to use to plot spatial units in the spatial object
                                    layer.name = "Occurrences", # To specify the legend name of this layer
                                    # map.types = mapviewGetOption("basemaps")
                                    map.types = c("OpenStreetMap", "Esri.WorldImagery", "Esri.WorldStreetMap")
                                    # map.types = c("OpenTopoMap", "Esri.WorldImagery", "Esri.WorldStreetMap", "OpenStreetMap", "Stamen.Watercolor") # To specify basemaps
                                    # ... # and many more, depending of the class of the spatial object
)

interactive_map_flags

interactive_map_flags_no_curation



## Save interactive maps in .rds
saveRDS(interactive_map_flags, file = "./maps/interactive_map_occurrences_flags.rds")
saveRDS(interactive_map_flags_no_curation, file = "./maps/interactive_map_occurrences_flags_no_curation.rds")


# To save the interactive map from mapview or leaflet into an interactive .hmtl or a static .png/.jpg/.pdf map
?mapview::mapshot

mapshot(x = interactive_map_flags,
        url = "./maps/interactive_map_occurrences_flags.html", # To save an interactive .html map
        remove_controls = NULL) # c("zoomControl", "layersControl", "homeButton", "scaleBar"))  # To specify which control buttons to remove (or not) from the output


### Plot interactive maps using tmaps and colors for type of erroneous records
# Once I curated the obvious mistakes. Create .rds and .html for each taxa
# Select the one with mistakes needing arbitration to show to AoW people





# Check correspondence of administrative regions across databases

table(AntWeb_database$Adm1)
table(GABI_database$adm1)

length(unique(AntWeb_database$Adm1))
length(unique(GABI_database$adm1))

table(AntWeb_database$Adm2)
table(GABI_database$adm2)

length(unique(AntWeb_database$Adm1))
length(unique(GABI_database$adm1))

table(is.na(AntWeb_database$LocLatitude))
table(is.na(GABI_database$dec_lat))


##### Attribute Natural Earth adm1 entities  #####

### Load Geographic Names Server data

GNS_adm1_df <- read.csv(file = "./input_data/GNS_Administrative_Regions/Administrative_Regions.txt", sep = "\t")

GNS_adm1_sf <- st_as_sf(x = GNS_adm1_df, coords = c("long_dd", "lat_dd"))

plot(GNS_adm1_sf[GNS_adm1_df$desig_cd == "ADM1", "adm1"], pch = 16)

Biogeographic_database_Ponerinae_NA_coords <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")

table(Biogeographic_database_Ponerinae_NA_coords$Source)


Taxa_in_phylo_NA_coords <- unique(Biogeographic_database_Ponerinae_NA_coords$Current_name[Biogeographic_database_Ponerinae_NA_coords$In_phylogeny])
Taxa_with_valid_name_NA_coords <- unique(Biogeographic_database_Ponerinae_NA_coords$Current_name[Biogeographic_database_Ponerinae_NA_coords$Status_AntCat == "valid"])

Taxa_for_analyses_NA_coords <- union(Taxa_in_phylo_NA_coords, Taxa_with_valid_name_NA_coords)
Taxa_valid_not_in_phylo_NA_coords <- setdiff(Taxa_for_analyses, Taxa_in_phylo_NA_coords)



setdiff(Taxa_in_phylo_NA_coords, Taxa_in_phylo)
setdiff(Taxa_for_analyses_NA_coords, Taxa_for_analyses)
setdiff(Taxa_valid_not_in_phylo_NA_coords, Taxa_valid_not_in_phylo)


setdiff()

# Use Natural Earth adm1 entities as the reference = adm1

# For GABI
# Draw random coordinates from bentity2 region + record uncertainty
# Find new polygon region based on new random coordinates

# For AntWeb
# Create table of AW_adm1 (AntWeb/Geographic Names Server) matches in NE_adm1
# May improve results using match with bentity2 polygons (if names are same as adm1 in AntWeb)
# If multiple match, random region will be picked and random coordinates will be drawn within this new NE region

# Clean duplicates
# Assign correspondence in both database
# Merge both databases

# Assign random coordinates =  draw random point within the area of the smallest geounit available (Bentity2 if no adm1 shapefile) 
# Compare area of geographic entity if possible to select the smallest one
# Use a field to inform on the nature of the coordinates: true, random

# map_subunits_sf <- ne_download(scale = 50, category = "cultural", type = 'map_subunits', returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_50m_admin_0_map_subunits/")

map_countries_sf <- ne_load(scale = 50, category = "cultural", type = 'countries', returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_50m_admin_0_countries/")
map_subunits_sf <- ne_load(scale = 50, category = "cultural", type = 'map_subunits', returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_50m_admin_0_map_subunits/")
map_NE_adm1_sf <- ne_load(file_name = "ne_10m_admin_1_states_provinces", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/ne_10m_admin_1_states_provinces/")

plot(map_NE_adm1_sf[map_NE_adm1_sf$admin == "Papua New Guinea", "adm1_code"])   # Sovereign administrative countries
View(map_NE_adm1_sf[map_NE_adm1_sf$admin == "Papua New Guinea", ])

# Target = GPS coordinates and smallest possible geounit for each specimen data

# The geounit will be used to draw a random GPS coordinates to use for gridded analyses
# For all entries with GPS coordinates, find the matching adm1 from Natural Earth
# Check correspondence of name with AW_adm1 and GABI_adm1
# Create a summary table of AW_adm1 and GABI_adm1 matches in NE_adm1
# For all entries without GPS coordinates
# Remove duplicates based on geounit location before to simulate random coordinates
# For GABI data = Drawn random coordinates within the GABI region and assign new NE_region
# Register uncertainty as the mean distance to centroid in the region   
# For AntWeb data = Match adm1 names
# If multiple match, select one randomly => flag this information
# Once the NE_adm1 is attributed, draw random GPS coordinates (or centroid if it falls within the region?) within it => flag this information
# Register uncertainty as the mean distance to centroid in the region
# Extract Bentity2 correspondance for every AntWeb coordinates

# May want to filter out data with too much uncertainty

# Check new randomized data for outliers the same way I did it for data with coordinates #


?ne_download


