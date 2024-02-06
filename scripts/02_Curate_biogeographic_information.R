##### Script 02: Curate biogeographic information #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate the specimen-level database for biographic information
# Aggregated information at taxa-level and update the taxa-level summary df

###

### Inputs

# Specimen-level database
# Taxa-level summary database
# NaturalEarth adiministrative shapefile
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


### 1.2/ Load databases to clean ####

## Load the curated databases
GABI_database_Ponerinae <- readRDS(file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database_Ponerinae.rds")
AntWeb_database_curated <- readRDS(file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")
Biogeographic_database_Ponerinae <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")

## Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")


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


### 2.1/ Assign country ISO codes ####

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


# Create new field for ISO3 names
Biogeographic_database_Ponerinae$Country_ISO3_name <- Biogeographic_database_Ponerinae$Country

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
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Czech Republic"] <- "Czech Rep."
# Biogeographic_database_Ponerinae$Country_ISO3[Biogeographic_database_Ponerinae$Country_ISO3_name == "Myanmar"] <- "Burma"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Laos"] <- "Lao PDR"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Zaire"] <- "Dem. Rep. Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Zaire"] <- "Democratic Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cote d'Ivoire"] <- "Ivory Coast"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cote d'Ivoire"] <- "Côte d'Ivoire"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "French Guiana"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Serbia"] <- "Republic of Serbia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "South Korea"] <- "Korea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "England"] <- "United Kingdom"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "North Korea"] <- "Dem. Rep. Korea"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo (Brazzaville)"] <- "Congo - Brazzaville"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo (Brazzaville)"] <- "Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo (Brazzaville)"] <- "Republic of the Congo"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent"] <- "St. Vincent & Grenadines"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent"] <- "St. Vin. and Gren."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent"] <- "Saint Vincent and the Grenadines"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent and the Grenadines"] <- "St. Vincent & Grenadines"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Vincent and the Grenadines"] <- "St. Vin. and Gren."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Christmas Island"] <- "Indian Ocean Ter."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Christmas Island"] <- "Indian Ocean Territories"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Trinidad"] <- "Trinidad & Tobago"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Trinidad and Tobago"] <- "Trinidad & Tobago"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Trinidad"] <- "Trinidad and Tobago"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bahamas"] <- "The Bahamas"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guinea-Bissau [Portuguese Guinea]"] <- "Guinea-Bissau"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea"] <- "Papua New Guinea" # May also be Indonesia
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Dominican Republic"] <- "Dominican Rep."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Martinique"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cape Verde"] <- "Cabo Verde"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mariana Islands"] <- "Northern Mariana Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mariana Islands"] <- "N. Mariana Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mariana Islands"] <- "Northern Mariana Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Mayotte"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Borneo"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia MULT"] <- "Malaysia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia, Federated States of"] <- "Micronesia (Federated States of)"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia, Federated States of"] <- "Micronesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia, Federated States of"] <- "Federated States of Micronesia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Federated States of Micronesia"] <- "Micronesia (Federated States of)"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Federated States of Micronesia"] <- "Micronesia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia"] <- "Micronesia (Federated States of)"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands"] <- "Solomon Is."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts and Nevis"] <- "St. Kitts & Nevis"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts and Nevis"] <- "St. Kitts and Nevis"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Great Britain"] <- "United Kingdom"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "British Virgin Islands"] <- "British Virgin Is."  # Could be many different Islands
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "?Unknown"] <- "Lesser Antilles"  # Could be many different Islands
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guadeloupe Islands"] <- "Guadeloupe"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guadeloupe Islands"] <- "France"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Guadeloupe"] <- "France"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "U.S. Virgin Islands"] <- "U.S. Virgin Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "U.S. Virgin Islands"] <- "United States Virgin Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Antigua and Barbuda"] <- "Antigua & Barbuda"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Antigua and Barbuda"] <- "Antigua and Barb."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea ERR"] <- "Papua New Guinea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "British West Indies"] <- "British West Indies" # Cannot define suitable ISO-3 entities
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea [Dutch New Guinea]"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Channel Islands"] <- "Guernsey" # Could also be Jersey
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Marshall Islands"] <- "Marshall Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Papua New Guinea [German New Guinea]"] <- "Papua New Guinea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Tokelau"] <- "New Zealand"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cocos (Keeling) Islands"] <- "Indian Ocean Ter."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cocos (Keeling) Islands"] <- "Indian Ocean Territories"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hispaniola"] <- "Dominican Republic"  # Could also be in Haiti
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hispaniola"] <- "Dominican Rep."  # Could also be in Haiti
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Lucia"] <- "St. Lucia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia MULT"] <- "Indonesia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts & Nevis"] <- "St. Kitts & Nevis"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts & Nevis"] <- "St. Kitts and Nevis"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Kitts & Nevis"] <- "Saint Kitts and Nevis"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "US Virgin Islands"] <- "U.S. Virgin Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "US Virgin Islands"] <- "U.S. Virgin Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "US Virgin Islands"] <- "United States Virgin Islands"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia ERR"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea ERR"] <- "Papua New Guinea" # May also be Indonesia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Brunei Darussalam"] <- "Brunei"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands ERR"] <- "Solomon Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands ERR"] <- "Solomon Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Comoros Islands"] <- "Comoros"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Northern Mariana Islands"] <- "N. Mariana Is."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Wallis and Futuna Islands"] <- "Wallis & Futuna"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Wallis and Futuna Islands"] <- "Wallis and Futuna Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Wallis and Futuna Islands"] <- "Wallis and Futuna"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea IND"] <- "Indonesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Indonesia [Dutch New Guinea]"] <- "Indonesia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Sao Tome and Principe"] <- "Sao Tome & Principe"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Sao Tome and Principe"] <- "São Tomé and Principe"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Martin"] <- "Saint Martin (French part)"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Saint Martin"] <- "St-Martin"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Netherlands Antilles"] <- "Aruba" # Could also be Curaçao or Bonaire (Netherlands)
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "UK"] <- "United Kingdom"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Equatorial Guinea"] <- "Eq. Guinea"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "New Guinea IND/PAPUA"] <- "Papua New Guinea" # May also be Indonesia
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Seychelles Islands"] <- "Seychelles"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Micronesia"] <- "Federated States of Micronesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands MULT"] <- "Solomon Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Solomon Islands MULT"] <- "Solomon Is."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bonaire, Sint Eustatius and Saba"] <- "Bonaire"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bonaire, Sint Eustatius and Saba"] <- "Netherlands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hong Kong"] <- "Hong Kong SAR China"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Hong Kong"] <- "Hong Kong S.A.R."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo"] <- "Congo - Brazzaville"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo"] <- "Republic of the Congo"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo [French Congo]"] <- "Congo - Brazzaville"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo [French Congo]"] <- "Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Congo [French Congo]"] <- "Republic of the Congo"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "South Sudan"] <- "S. Sudan"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bosnia and Herzegovina"] <- "Bosnia & Herzegovina"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Bosnia and Herzegovina"] <- "Bosnia and Herz."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Palestine"] <- "Palestinian Territories"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia Peninsular"] <- "Malaysia"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Ivory Coast"] <- "Côte d'Ivoire"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macedonia"] <- "North Macedonia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Malaysia ERR"] <- "Malaysia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Swaziland"] <- "eSwatini"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Central African Republic"] <- "Central African Rep."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macau"] <- "Macao"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macau"] <- "Macao S.A.R"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "?Unknown (Somalia? Or Ethiopia?)"] <- "Somaliland" # Could be Ethiopia/Somalia
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Democratic Republic of Congo"] <- "Congo - Kinshasa"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Democratic Republic of Congo"] <- "Dem. Rep. Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Democratic Republic of Congo"] <- "Democratic Republic of the Congo"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Port of Entry"] <- NA
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "null"] <- NA
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Unknown"] <- NA
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Cook Islands"] <- "Cook Is."
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "United States Virgin Islands"] <- "U.S. Virgin Islands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "United States Virgin Islands"] <- "U.S. Virgin Is."
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Rhodesia"] <- "Zimbabwe"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Timor-Leste"] <- "East Timor"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Macaronesia"] <- "Spain"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Panam·"] <- "Panama"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Europa Island"] <- "French Southern Territories"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Europa Island"] <- "French Southern and Antarctic Lands"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "French Polynesia"] <- "Fr. Polynesia"
Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "Reunion"] <- "France"
# Biogeographic_database_Ponerinae$Country_ISO3_name[Biogeographic_database_Ponerinae$Country_ISO3_name == "French Southern Territories"] <- "Fr. S. Antarctic Lands"

# Biogeographic_database_Ponerinae[Biogeographic_database_Ponerinae$Country == "Macaronesia", ]

# Add ISO-3 code to Biogeographic database
# Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_metadata$iso3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_metadata$name)]
# Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_NE_metadata$iso_a3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_NE_metadata$name)]
Biogeographic_database_Ponerinae$Country_ISO3_code <- Countries_NE_metadata$ISO_A3[match(x = Biogeographic_database_Ponerinae$Country_ISO3_name, table = Countries_NE_metadata$SUBUNIT)]

table(Biogeographic_database_Ponerinae$Country_ISO3_code)

View(Biogeographic_database_Ponerinae[(is.na(Biogeographic_database_Ponerinae$Country_ISO3_code)), ])

## Save specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")

### 2.2/ Remove "Port of Entry data" ####

table(Biogeographic_database_Ponerinae$Country == "Port of Entry")

Biogeographic_database_Ponerinae <- Biogeographic_database_Ponerinae[!(Biogeographic_database_Ponerinae$Country == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_flags <- Biogeographic_database_Ponerinae_flags[!(Biogeographic_database_Ponerinae_flags$Country == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_NA_coords <- Biogeographic_database_Ponerinae_NA_coords[!(Biogeographic_database_Ponerinae_NA_coords$Country == "Port of Entry"), ]
# Biogeographic_database_Ponerinae_no_NA_coords <- Biogeographic_database_Ponerinae_no_NA_coords[!(Biogeographic_database_Ponerinae_no_NA_coords$Country == "Port of Entry"), ]

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

# Oronoque River, Guyana
Biogeographic_database_Ponerinae$Longitude_dec[(str_detect(string = Biogeographic_database_Ponerinae$Locality, pattern = "Oronoque River")) & !is.na(Biogeographic_database_Ponerinae$Latitude_dec) & is.na(Biogeographic_database_Ponerinae$Longitude_dec)] <- -57.422

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

# Load flagged output from CoordinateCleaner
Biogeographic_database_Ponerinae_flags <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")

## 2.6.1/ Prepare data ####

## Load output from CoordinateCleaner
Biogeographic_database_Ponerinae_flags <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")

### Create new variables to record former coordinates and new ones
# Biogeographic_database_Ponerinae_flags$Latitude_initial <- Biogeographic_database_Ponerinae_flags$Latitude_dec
# Biogeographic_database_Ponerinae_flags$Longitude_initial <- Biogeographic_database_Ponerinae_flags$Longitude_dec

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

table(Biogeographic_database_Ponerinae_flags$.con)

## 2.6.2/ Curate records falling in wrong country ####

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "countries", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality))

test <- Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "countries", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality)
table(is.na(test$GABI_accession_ID))

Biogeographic_database_Ponerinae[Biogeographic_database_Ponerinae$Specimen_code == "casent0752451", ] 

## 2.6.2.1/ False positive: case of records flagged as wrong, but that are actually right (tiny islands not in the shapefiles) ####

# Australia, Islands in Torres Strait, Queensland
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- "Murray Island"
Biogeographic_database_Ponerinae_flags$.con[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00005107", "GABI_00005108")] <- TRUE
# Australia, Islands along Queensland coast. All records are within/around tiny islands and are correct
plot(Biogeographic_database_Ponerinae_flags[(Biogeographic_database_Ponerinae_flags$flag_type == "countries") & (Biogeographic_database_Ponerinae_flags$adm1 == "Queensland"), "adm2"])
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
Biogeographic_database_Ponerinae_flags$Locality[(Biogeographic_database_Ponerinae_flags$Latitude_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == 161.97)] <- "Malaupaina Island"
Biogeographic_database_Ponerinae_flags$.con[(Biogeographic_database_Ponerinae_flags$Latitude_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == 161.97)] <- "TRUE"
Biogeographic_database_Ponerinae_flags$.sea[(Biogeographic_database_Ponerinae_flags$Latitude_initial == -10.248) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == 161.97)] <- "TRUE"
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
Biogeographic_database_Ponerinae_flags$.con[replace_na((Biogeographic_database_Ponerinae_flags$adm2 %in% c("Isla de Cíes", "Isla de Ons")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na((Biogeographic_database_Ponerinae_flags$adm2 %in% c("Isla de Cíes", "Isla de Ons")) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- TRUE

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
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality %in% c("Iguazu National Park", "2.9 km SW of San ignacio")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality %in% c("Iguazu National Park", "2.9 km SW of San ignacio")])
# Argentina, San Andresito => Wrong longitude. -54.5 instead of -55.9
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "San Androcito"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -54.5
# Argentina, Reserva El Loro Hablador 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "El Loro Hablador"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -62.35
# Argentina, Parque Nacional Chaco => Wrong longitude. -59.61 instead of -41.04
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Parque Nacional Chaco"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -59.61

# Cocos/Keeling Islands. Wrong ISO country (currently in Australia) + wrong rounded coordinates
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- 96.90
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -12.116666670
Biogeographic_database_Ponerinae_flags$Country_adm1[(Biogeographic_database_Ponerinae_flags$Longitude_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "Cocos (Keeling) Islands"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[(Biogeographic_database_Ponerinae_flags$Longitude_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "Indian Ocean Territories"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[(Biogeographic_database_Ponerinae_flags$Longitude_initial == 96.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -12.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- "CCI"

# Australia, Queensland => Wrong longitude. Fall in ocean. 145.65 instead of 145.76666
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00651269")] <- 145.65
# Australia, Northern Territory, Howard Springs => Wrong rounded latitude falling North in the Indian Ocean. -12.50 instead of -12.00
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Howard Springs") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -12.50
# Australia, North Western Australia => Wrong latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality %in% c("TCMBW")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Locality %in% c("TCMBW")])
# Australia, North Western Australia => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$bentity2_name %in% c("South Western Australia")] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$bentity2_name %in% c("South Western Australia")])

# Barbados => Wrong longitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")])
Biogeographic_database_Ponerinae_flags$bentity2_name[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Barbados")] <- "Lesser Antilles"
# Barbados => wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_initial == -59.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == 13.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- -59.58
Biogeographic_database_Ponerinae_flags$Latitude_dec[(Biogeographic_database_Ponerinae_flags$Longitude_initial == -59.00000) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == 13.000000) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")] <- 13.20

# Brazil => Wrong longitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Brazil")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Brazil")])
# Brazil, Ceara => Wrong Latitude. Longitude should be negative. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Ceara")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Ceara")])
# Brazil, Espirito Santo, Linhares => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = ((Biogeographic_database_Ponerinae_flags$Locality == "Linhares") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")), replace = F)] <- -40.04
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = ((Biogeographic_database_Ponerinae_flags$Locality == "Linhares") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries")), replace = F)] <- -19.41
# Brazil, Mato Grosso do Sul, Serra da Bodoquena National Park => Wrong coordinates falling in Paraguay
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Da Mata farm")), replace = F)] <- -56.74
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kadiweu reserve")), replace = F)] <- -57.25
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Salobra river")), replace = F)] <- -56.76
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Salobra river")), replace = F)] <- -20.86
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sta Maria farm-Perdido river")), replace = F)] <- -56.92
# Brazil, Minas Gerais, Juiz de Fora => Wrong latitude. Latitude should be negative. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Juiz de Fora")), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Juiz de Fora")), replace = F)])
# Brazil, Rio de Janeiro => Wrong latitude ending in the ocean. -22.919 instead of 23.919
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rio de Janeiro") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -22.919
# Brazil, Rondônia, Porto Velho => Wrong longitude ending in the ocean. -63.904 instead of -93.904
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Port Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -63.904
# Brazil, Rondônia, Porto Velho => Latitude and Longitude are reversed!
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- Biogeographic_database_Ponerinae_flags$Latitude_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- Biogeographic_database_Ponerinae_flags$Longitude_initial[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "Porto Velho") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)]
# Brazil, Santa Catarina => Wrong signs for latitude/longitude. Should be negative + one records need to be shift back in the landmass
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Santa Catarina") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -48.553
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Santa Catarina") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.158
# Brazil, Sao Paulo, Mogi das Cruzes => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mogi das Cruzes") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -46.142
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mogi das Cruzes") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.532
# Brazil, Sao Paulo, Parque Leon Feffer => Wrong rounded coordinates ending in the ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Parque Leon Feffer") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -46.223
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Parque Leon Feffer") & Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.529

# Cameroon, Mbale Mejo to Ekingli => Wrong rounded coordinates
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c002505")] <- 11.717
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("sam-hym-c002505")] <- 3.603

# China, Guangxi, Fangchenggang City, Mt. Shiwandashan => Wrong latitude, falling into China Sea => 21.91 instead of 21.18
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02357907")] <- 21.91
# China, Yunnan => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Yunnan")] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Yunnan")])
# China, Yunnan, Bakaxiaozhai, Menglun Town => Wrong coordinates falling in Vietnam
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Menglun Town") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 101.254
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Menglun Town") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 21.931

# Colombia, Atlantico, Piojó => Wrong coordinates. Fall into Caribbean sea North of the coast
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Piojó") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.108
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Piojó") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.749
# Colombia, Bolivar, Los Colorados Venado => Wrong longitude. -75.116 instead of -52.116
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Los Colorados Venado") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.116
# Colombia, Cauca, PNN Gorgona Mancora => Wrong longitude. -78.183 instead of -18.183
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "PNN Gorgona Mancora") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -78.183
# Colombia/Venezuela, Cesar, Sierra de Perija, Socorpa Mission. NP is in Venezuela so rounded coordinates are likely wrong. # Need also to change the country.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0178780"] <- -72.867407
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code == "casent0178780"] <- 9.948540
# Colombia, La Guajira => Wrong latitude. Latitude should be positive. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("La Guajira")] <- abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("La Guajira")])
# Colombia, Magdalena, Santa Marta => Wrong coordinates. Fall into Caribbean sea North of the coast. 11.242 instead of 11.352
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00529958"), replace = F)] <- -74.205
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00529958"), replace = F)] <- 11.242
# Colombia, Narino, Barbacoas, RN Río Ñambí => Wrong latitude. Latitude should be positive. 
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Narino")] <- abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Narino")])
# Colombia, Putumayo, Puerto Leguízamo => Wrong longitude. -74.983300 instead of -79.983300
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bajo Cusacante") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -74.983300
# Colombia, Risaralda, El Trapiche => Wrong longitude. -75.954167 instead of -6.954167
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "El Trapiche") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -75.954167
# Colombia, Valle del Cauca, Farallones de Cali National Park => Wrong longitude. -76.656000 instead of -79.656000
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Farallones de Cali National Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -76.656000
# Colombia, Valle del Cauca, Bajo Calima, Buenaventura => Wrong longitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bajo Calima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bajo Calima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])


# Costa Rica, Alajuela => Coordinates reversed AND wrong signs! (Master class...)
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")] <- -1 * Biogeographic_database_Ponerinae_flags$Latitude_init[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")]
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")] <- -1 * Biogeographic_database_Ponerinae_flags$Longitude_init[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01692539", "GABI_01692601")]
# Costa Rica, Limon, Salsipuedes => Wrong latitude. Fall North to Limon in the Carribean sea. 10.007 instead of 10.07
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Salsipuedes") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.007
# Costa Rica, Puntarenas, Sirena => Wrong coordinates. Fall South in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sirena") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -83.59129
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sirena") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 8.480170
# Costa Rica, Puntarenas, Estacion Biological Las Cruces => Wrong longitude. -82.960 instead of -82.010000
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Estacion Biological Las Cruces") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -82.960
# Costa Rica, La Virgen => Wrong longitude. -83.75 instead of -8.75
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00527157"), replace = F)] <- -83.75

# Croatia, Dalmacia, Sucurac => Wrong coordinates falling in Montenegro
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sucurac") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.431
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sucurac") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 43.553
# Croatia, Hvar Island, Hvar => Wrong Longitude falling in the Adriatic Sea. 16.453 instead of 14.433333
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Hvar") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.453
# Croatia, Dalmacia, Kastela => Wrong coordinates falling in Bosnia & Herzegovina
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kastela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.350
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kastela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 43.555

# Cuba, Cumanayagua, Mina Carlota => Wrong coordinates falling in Mexico !
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cumanayagua, Mina Carlota") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -80.16667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cumanayagua, Mina Carlota") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 22.06667

# DRC, Masosa => Wrong coordinates in Congo. Not sure of the new ones but better than random.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Masosa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 22.323
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Masosa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 2.215
# DRC, Manamama => Wrong coordinates in Congo. Not sure of the new ones but better than random.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Manamana") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 19.232
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Manamana") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -4.076

# DRC/Central African Republic, Bangassou => Wrong rounded latitude (too South) and Wrong country
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bangassou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.739

# Ecuador, Los Rios, Rio Palenque => Wrong latitude. Latitude should be negative.  
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Los Rios")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Los Rios")])
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
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mt. Unibot") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.628
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mt. Unibot") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.368
# Micronesia, Pohnpei, Ponape Agriculture & Trade School => Wrong rounded coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Ponape Agriculture & Trade School") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 158.308
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Ponape Agriculture & Trade School") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.842
# Micronesia, Woleai Atoll => Wrong coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Woleai Atoll") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 143.913
# Micronesia, Nama Island => Wrong coordinates falling in the lagoon
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nama Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 152.576
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nama Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.992

# French Guyana, Campus Agronomique in Kourou => Wrong longitude. Longitude should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country %in% c("French Guiana")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country %in% c("French Guiana")])
# French Guyana, Sinnamary, Petit Saut => Wrong coordinates in Brazil
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00530228"), replace = F)] <- -53.050
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00530228"), replace = F)] <- 5.064

# Gabon, Aire d'Exploition Rationnelle de Faune des Monts Doudou => Wrong Latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Monts Doudou")), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Monts Doudou")), replace = F)])

# Greece, Corfú => Wrong latitude falling in the Ionian sea South to Corfu. 39.456667 instead of 39.286667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Corfú")), replace = F)] <- 39.456667
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
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country %in% c("Guyana")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country %in% c("Guyana")])

# India, Arunachal Pradesh, Lumla => Wrong coordinates falling in Bhoutan
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01662113")), replace = F)] <- 91.722
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01662113")), replace = F)] <- 27.530

# Indonesia, Sulawesi, Central Sulawesi Province, Berdikari => Wrong Latitude falling North in the Celebes Sea. 0.734 instead of 1.134
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Berdikari") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 0.734
# Indonesia, Sulawesi, Central Sulawesi Province, Bulili => Wrong Longitude falling East in the Celebes Sea. 120.0955 instead of 120.955
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bulili") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 120.0955
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
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pulo Laut") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 107.977
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pulo Laut") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.715

# Ivory Coast, Tai National Park => Wrong longitude falling in Liberia. -7.09642 instead of -7.89642 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Foret de Tai") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -7.09642
# Ivory Coast, Tai Forest => Wrong rounded longitude falling in Liberia. -7.066667 instead of -7.666667
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tai Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -7.066667

# Japan, Iwate Prefecture => Wrong Latitude falling South in the Pacific Ocean. 39.012222 instead of 33.012222
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Iwate Prefecture") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 39.012222
# Japan, Kagoshima Prefecture, Kagoshima City, Jigenji Park => Wrong coordinates falling South in the Pacific ocean 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Jigenji Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 130.503
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Jigenji Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 31.511
# Japan, Kagoshima Prefecture, Sakurajima => Wrong coordinates falling South in the Pacific ocean 
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sakurajima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 130.650
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sakurajima") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 31.583
# Japan, close to Akusekijima Island => Wrong coordinates falling East of Akusekijima Island
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00577065", "GABI_00577066")), replace = F)] <- 129.604
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00577065", "GABI_00577066")), replace = F)] <- 29.459

# Lesotho, Likhoele => Wrong Latitude falling North in South Africa. -29.81667 instead of -29.183333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Likhoele") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -29.81667
# Lesotho, Mokhotlong => => Wrong coordinates falling East in South Africa.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mokhotlong") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.08333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mokhotlong") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -29.71667
# Lesotho, Qachas Nek => => Wrong latitude falling South in South Africa. -30.1 instead of -30.9
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Qachas Nek") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -30.1

# Madagascar => Wrong Latitude. Latitude should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Madagascar")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Madagascar")])
# Madagascar, Marojejy National Park => Wrong Longitude falling West in the Mozambique Canal. 49.934722 instead of 46.934722
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00557292"), replace = F)] <- 49.934722

# Malaysia, Kelantan, Melawi Beach => Wrong coordinates falling NW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Melawi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 102.420
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Melawi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 6.021
# Malaysia, Borneo, Sarawak, Kampong Sega => Wrong coordinates falling in Kalimantan, Indonesia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kampong Segu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 110.238
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kampong Segu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 1.424

# Malaysia, Borneo, Sarawak, Semengoh Forest Reserve => Wrong Latitude falling in Kalimantan, Indonesia. 1.401 instead of 0.8833333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Semengoh Forest Reserve") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 1.401
# Malaysia, Borneo, Sarawak, Gunung Mulu National Park => Wrong Latitude falling in Brunei. 4.037 instead of 4.95
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Gunung Mulu National Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.037
# Malaysia, Borneo, Sarawak, Lambir Hills National Park => Wrong coordinates falling in Kalimantan, Indonesia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lambir Hills") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 114.04
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lambir Hills") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.200
# Borneo, Tobang => Coordinates are centroids of Borneo. Could relate to Tabang in Kalimantan, Indonesia.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 116.017
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 0.575
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tobang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "IDN"

# Mali, Baboye => Wrong coordinates falling South in Burkina Faso
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Baboye") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.952
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Baboye") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 14.260

# Mauritius, Round Island => Wrong Longitude falling West in the Indian ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Round Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 57.784

# Mexico, Chiapas, Irlanda => Wrong coordinates falling in Guatemala
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Irlanda") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -93.332
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Irlanda") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.109
# Mexico, Chiapas, Soconusco => Wrong coordinates falling in Guatemala
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01103516", "GABI_01103514"), replace = F)] <- "Soconusco region"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soconusco") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -92.297
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soconusco") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 15.144
# Mexico, Jalisco, Estación Biológica Chamela. Wrong Latitude falling South into the Pacific Ocean 19.50 instead of 19.05
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Estación Biológica Chamela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 19.50
# Mexico, Quintana Roo, Majahual => Wrong coordinates falling NE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Majahual") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -87.71013
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Majahual") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.71277
# Mexico, Quintana Roo, Puerto Morelos => Wrong rounded coordinates falling SE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Puerto Morelos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -86.90278
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Puerto Morelos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 20.84653 
# Mexico, Quintana Roo, Centro de Investigaciones Costeras La Mancha => Wrong coordinates with sign inversion falling in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01752436", "GABI_01752437"), replace = F)] <- -96.38125
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01752436", "GABI_01752437"), replace = F)] <- 19.59494

# Mozambique, Maputo, Delagoa Bay => Wrong rounded longitude falling East in the Bay
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Delagoa Bay") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 32.561

# Namibia, Oidimba Village => Wrong coordinates falling NW into Angola
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Edimba") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 16.496
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Edimba") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.476
# Namibia, Katima Mulilo => Wrong coordinates falling North into Zambia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Katima Mulilo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 24.26667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Katima Mulilo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.50

# Palau, Ngatpang => Wrong rounded coordinates falling in SW in the Pacific
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Ngatpang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 134.53244
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Ngatpang") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.45053
# Palau, Peleliu => Wrong rounded coordinates falling in SW in the Pacific
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 134.24308
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.01506
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Peleliu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Peleliu"

# Panama, STRI Gamboa, Pipeline Road => Wrong Longitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Panama")] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Country_ISO3_name %in% c("Panama")])
# Panama, Nusagandi => Wrong rounded coordinates in Costa Rica!
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nusagandi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -78.97927
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nusagandi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 9.351060

# PNG, Bougainville Island, Buin Village
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buin") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 155.688
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buin") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -6.745

# PNG, Hermit group, Luf Island => Wrong Longitude falling East in the Pacific Ocean. 145.063 instead of 145.083
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Luf Is") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 145.063
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00654783"), replace = F)] <- 145.063
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00654783"), replace = F)] <- "Hermit group, Luf Island"
# PNG, Manus, Lorengau => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lorengau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 147.272
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lorengau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -2.031
# PNG, Central Province, Wanigela => Wrong Latitude falling South in Pacific Ocean
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Wanigela") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.062
# PNG, Central Province, Bisianumu, near Sogeri => Wrong Latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bisianumu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bisianumu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])
# PNG, East New Britain Province, Lamas => Wrong coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lamas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.404
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lamas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -5.6098
# PNG, East New Britain Province, Vouvou => Wrong coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vouvou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 151.4594
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vouvou") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -5.440
# PNG, Sandaun Province, Aitape => Wrong rounded coordinates falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Aitape") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 142.35
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Aitape") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.13333 
# PNG, Morobe Province, Lae => Wrong coordinates falling far away
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lae") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 147.000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lae") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -6.73333
# PNG, Mussau Talumalaus => Wrong rounded coordinates falling South in the Bismark Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mussau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 149.621
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mussau") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.436
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec == 150) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 149.621
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec == -5) & (Biogeographic_database_Ponerinae_flags$Longitude_dec == 149.621) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.436
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec == -1.436) & (Biogeographic_database_Ponerinae_flags$Longitude_dec == 149.621) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Bismarck Archipelago, Mussau Island"
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_dec == -1.436) & (Biogeographic_database_Ponerinae_flags$Longitude_dec == 149.621) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "New Ireland"
# PNG, Western Province, Fly River => Wrong longitude falling East into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Fly River") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 143.629

# Peru, Madre de Dios, Sachavacayoc Center => Wrong Latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sachavacayoc") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sachavacayoc") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])

# Philippines, no locality => Wrong Longitude falling West of the coast of South Luzon. Bring them back on the land assuming minimal error in coordinates. 124.13750 instead of 124.26750
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_initial == 12.877214) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 124.13750

# Puerto Rico, Bosque Estatal Guanica => Wrong Longitude falling South in the Caribbean Sea. 17.97 instead of 17.84
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bosque Estatal Guanica") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 17.97
# Puerto Rico, Carite Lake region => Wrong rounded coordinates falling SE in the Caribbean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Carite") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -66.100
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Carite") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.074

# Rwanda, Rubona => Wrong coordinates falling SW in DRC
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rubona") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.256
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rubona") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1.728

# Saint Vincent and the Grenadines, Savan Island => Wrong Longitude falling East in the Atlantic Ocean. -61.21111 instead of -61.12111
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Savan Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.21111

# Samoa => Wrong Latitude. Should be negative
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Country, pattern = "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Country, pattern = "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)])
# Samoa => Wrong Latitude falling South in the Pacific ocean. -13.583 instead of -16.583
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00761072", "GABI_00761107"), replace = F)] <- -13.583

# Saudi Arabia, Al Qatif => => Wrong Latitude falling North in the Arabic Sea. 26.533333 instead of 26.933333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Al Qatif") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 26.533333

# Solomon Islands, Isabel Province => Wrong rounded longitude falling West in the Pacific Ocean. 159.115 instead of 159.5
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Isabel Province") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 159.115
# Solomon Islands, Cristobal Island, Wainoni => Wrong Latitude falling North in the Pacific Ocean. -10.566670 instead of -10.066670/-10.100000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Wainoni") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.566670
# Solomon Islands, Santa Cruz Group, Nendo Island => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Locality[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_initial == -11) & is.na(Biogeographic_database_Ponerinae_flags$Locality) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Santa Cruz Group, Nendo Island"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_initial == -11) & (Biogeographic_database_Ponerinae_flags$Locality == "Santa Cruz Group, Nendo Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 165.93927
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na((Biogeographic_database_Ponerinae_flags$Latitude_initial == -11) & (Biogeographic_database_Ponerinae_flags$Locality == "Santa Cruz Group, Nendo Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.718873
# Solomon Islands, Reef Islands, Nibanga Temau => Wrong Longitude falling East in the Pacific Ocean. 166.310 instead of 166.167
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$Specimen_code %in% c("anic32-046469")] <- 166.310
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- "Reef Islands, Nibanga Temau"
Biogeographic_database_Ponerinae_flags$Longitude_dec[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00669091")] <- 166.310
# Solomon Islands, Reef Islands, Matema Island => Wrong rounded coordinates falling SE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matema") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 166.184
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matema") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -10.292
# Solomon Islands, Santa Cruz Group, Vanikoro Island => Wrong rounded coordinates falling NW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 166.909
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Vanikoro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -11.674
# Solomon Islands, Santa Cruz Group, Anuta Island => Wrong rounded coordinates falling SW in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 169.850
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Anuda Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -11.611


# South Africa, East London, Peddie District => Wrong Latitude falling South into the Ocean. -33.03333 instead of -33.125
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "East London") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -33.03333
# South Africa, Rainwoods (Small Holding), Near Witelsbos => Wrong Latitude falling South into the Ocean. -34.000 instead of -34.125000
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Witelsbos") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -34.000
# South Africa, Nongoma, Zululand => Wrong Latitude falling North in eSwatini. -27.9 instead of -27.1
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nongom") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.9
# South Africa, Ntendeka Forest Reserve => Wrong Latitude falling North in eSwatini. -27.85 instead of -27.15
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Ntendeka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -27.85
# South Africa, Umtamvuna Nature Reserve => Wrong Latitude falling South in the Indian Ocean. -31.031 instead of -31.228889
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Umtamvuna") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -31.031
# South Africa, Entabeni Forest Reserve => Wrong Latitude falling North in Botswana. -22.93333 instead of -22.066667
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Entabeni") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -22.93333
# South Africa, Soutpansberg Mountain => Wrong coordinates falling NW in Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soutpansberg") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.42801
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Soutpansberg") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -23.03564
# South Africa, Eastern Transvaal, Nelspruit => Wrong latitude. Should be negative.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nelspruit"), replace = F)] <- -1 * abs(Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Nelspruit"), replace = F)])
# South Africa, Thursford Farm => Wrong longitude. Should be positive.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Thursford"), replace = F)] <- abs(Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Thursford"), replace = F)])

# Zambia, Lusaka => Wrong coordinates falling NW in Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 28.31938
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -15.43585

# South Korea, Site 114 => Wrong Longitude, falling East in the Pacific Ocean. 129.32 instead of 129.42.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "114") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.32
# South Korea, Site 115 => Wrong Longitude, falling East in the Pacific Ocean. 129.319167 instead of 129.419167.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "115") & (Biogeographic_database_Ponerinae_flags$Country == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.319167
# South Korea, Site 116 => Wrong Longitude, falling East in the Pacific Ocean. 129.323889 instead of 129.433889.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "116") & (Biogeographic_database_Ponerinae_flags$Country == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.323889
# South Korea, Site 51 => Wrong Longitude, falling East in the Pacific Ocean. 128.425 instead of 128.535.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "51") & (Biogeographic_database_Ponerinae_flags$Country == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 128.425
# South Korea, Site 58 => Wrong Longitude, falling East in the Pacific Ocean. 128.529167 instead of 128.629167.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "58") & (Biogeographic_database_Ponerinae_flags$Country == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 128.529167
# South Korea, Site 159 => Wrong Longitude, falling East in the Pacific Ocean. 129.345278 instead of 129.445278.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "159") & (Biogeographic_database_Ponerinae_flags$Country == "South Korea") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 129.345278

# Spain, Illes Balears, Mallorca => Wrong coordinates falling in the Mediterranean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00794898", "GABI_00794899"), replace = F)] <- 2.986
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00794898", "GABI_00794899"), replace = F)] <- 39.615
# Spain, Illes Canarias, El Hierro => Wrong Latitude falling North in the Atlantic Ocean. 27.75 instead of 27.85
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "El Hierro") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 27.75
# Spain, Illes Canarias, Gran Canaria Island, Las Palmas => Wrong Longitude falling East in the Atlantic Ocean. -15.428611 instead of -15.318611
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Las Palmas") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -15.428611
# Spain, Illes Canarias => Wrong coordinates falling in the middle of the archipelago, in the Atlantic Ocean.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00831742"), replace = F)] <- -15.70393
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00831742"), replace = F)] <- 28.10043
# Spain, Arroyo de la Ermita, Sierra Almijara => Wrong coordinates falling SE in the Mediterranean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Arroyo de la Ermita") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -3.969
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Arroyo de la Ermita") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 36.894

# São Tomé and Principe, São Tomé Island => Wrong Latitude falling North in the Atlantic Ocean. 0.26492 instead of 0.56492
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Specimen_code == "ey20045"), replace = F)] <- 0.26492

# Trinidad & Tobago, "Trinidad_b" => Falling in the Atlantic Ocean between the two islands. Take random coordiantes on Trinidad instead.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_code, pattern = "Trinidad_b") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.236
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality_code, pattern = "Trinidad_b") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 10.416

# USA, California, Cold Canyon => Wrong latitude falling South into the Pacific Ocean. 38.511 instead of 36.511.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cold Canyon") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 38.511
# USA, Florida, Franklin County => Wrong rounded latitude falling South into the Caribbean Sea Ocean. 29.800 instead of 36.511.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Franklin County") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 29.916666
# USA, Florida, Seminole State Park => Wrong coordinates falling SE into the Carribean Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Seminole State Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -81.604
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Seminole State Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 25.976
# USA, Georgia, Cloud Land Canyon => Wrong Latitude falling South into the Carribean Sea. 34.834444 instead of 24.834444.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cloudland Canyon") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 34.834444
# USA, Georgia, Tugaloo => Wrong Latitude falling South into the Carribean Sea. 34.494722 instead of 24.484722.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tugaloo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 34.494722
# USA, Mississippi, Grand Bay Savanna => Wrong Latitude falling South into the Carribean Sea. 30.36 instead of 30.0475.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Grand Bay Savanna") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 30.36
# USA, Texas, Harlingen => Wrong Latitude falling South into the Carribean Sea. 26.165 instead of 23.6.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Harlingen") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 26.165
# USA, Hawaii, Oahu Island => Wrong rounded Latitude falling South into the Carribean Sea. 21.5 instead of 21.0.
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Oahu Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 21.5

# Vanuatu, Efate Island, Port Vila => Wrong rounded coordinates falling NW into the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Port Fila") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 168.32
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Port Fila") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -17.736

# Venezuela, Valle del Abismo => Wrong coordinates falling SE in Brazil
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Valle del Abismo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -61.6
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Valle del Abismo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 4.40
# Venezeuala, Azevedo District, Cupo => Wrong Longitude falling East in the Caribbean Sea. -66.372 instead of -65.622
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cupo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -66.372
# Venezeuala, Nueva Esparta, Isla Margarita, ca. Fuentidueño => Wrong Longitude falling West in the Caribbean Sea. -63.913056 instead of -64.913056
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Fuentidueño") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -63.913056
# Venezuela, San Critobal, Loma de Pio => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Loma de Pio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.200
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Loma de Pio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.749
# Venezuela, San Antonio del Tachira => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "San Antonio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.442
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "San Antonio") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.818
# Venezuela, Tachira, San Cristobal, Santa Teresa => Wrong coordinates falling NW in Colombia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Santa Teresa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -72.225
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Santa Teresa") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 7.800

# Viet Nam, Huong Son => Wrong coordinates falling SE in the China Sea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Huong Son") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 105.426
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Huong Son") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 18.512
# Viet Nam, Ba Vi National Park => Wrong Longitude falling East in the China Sea. 105.36083 instead of 108.36083
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Ba Vi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- 105.36083

# Wallis and Futuna, Wallis Island => Wrong coordinates falling NE in the Pacific Ocean
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00758995", "GABI_00754808"), replace = F)] <- -176.202
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00758995", "GABI_00754808"), replace = F)] <- -13.285

# Zambia, Mpulungu => Wrong rounded Latitude falling North in Lake Tanganyika. -8.83333 instead of -8.500
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mpulungu") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -8.83333

# Zimbabawe, Matopo Hills => Wrong Latitude falling South in South Africa. -20.45 instead of -27.45
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matopo") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -20.45

# eSwatini, Pigg's Peak => Latitude falling North in South Africa. -26.033333 instead of -25.033333
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Piggspeak") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- -26.033333


Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Piggspeak")), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Piggspeak")), replace = F)]

plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Australia"), "ISO_A3"]) # Includes the Canarias
plot(Countries_NE_sf[, "ISO_A3"])

plot(Countries_NE_sf[Countries_NE_sf$SUBUNIT %in% c("Palestine"), "ISO_A3"]) # Includes the Canarias

## 2.5.2.3/ True positive: Need to change country name ####

# American Samoa => Records are in the Western Samoa Islands, not the American Samoa
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00941224", "GABI_00941234")] <- "WSM"

# Norfolk Island => Records are in Norfolk Island, not Australia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Norfolk Island")] <- "Norfolk Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Norfolk Island")] <- "NFK"

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
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(Biogeographic_database_Ponerinae_flags$Locality, pattern = paste0("Sierra de Perija", "|","Sierra de Parija")), replace = F)] <- "Venezuela"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(Biogeographic_database_Ponerinae_flags$Locality, pattern = paste0("Sierra de Perija", "|","Sierra de Parija")), replace = F)] <- "VEN"

# Eritrea instead of some Ethiopia records
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "Tesseney"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "Eritrea"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00824708")] <- "ERI"

# Clipperton Island instead of France
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Clipperton Island")] <- "Clipperton Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$adm1 %in% c("Clipperton Island")] <- "CPT"

# Guyana instead of French Guyana/France for an erroneous record
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00845045")] <- "Guyana"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00845045")] <- "GUY"

# French Southern and Antarctic Lands instead of Europa Island/France for Europa Island
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality %in% c("Europa Island")] <- "French Southern and Antarctic Lands"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality %in% c("Europa Island")] <- "ATF"

# Jordan instead of Israel. Coordinates falling in the Jordan side of the river.
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality %in% c("Jourdain")] <- "Jordan"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality %in% c("Jourdain")] <- "JOR"

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
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Country %in% c("?Unknown")] <- "France"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Country %in% c("?Unknown")] <- "FRA"
Biogeographic_database_Ponerinae_flags$Country[Biogeographic_database_Ponerinae_flags$Country %in% c("?Unknown")] <- "Martinique"

# Namakunde in Angola instead of Namibia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$Locality %in% c("Namakunde")] <- "Angola"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$Locality %in% c("Namakunde")] <- "AGO"
# Erikson's Drift, Kunene River in Angola instead of Namibia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01943756")] <- "Angola"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01943756")] <- "AGO"

# Kongga/Buin Village in PNG instead of Solomon Islands
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buin"), replace = F)] <- "Autonomous Region of Bougainville [North Solomons]"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buin"), replace = F)] <- "Papua New Guinea"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buin"), replace = F)] <- "PNG"

# Sorong in Indonesia, West Papua instead of PNG
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sorong"), replace = F)] <- "West Papua"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sorong"), replace = F)] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sorong"), replace = F)] <- "IDN"

# Bali, Laba Sari in Indonesia instead of PNG
Biogeographic_database_Ponerinae_flags$Locality[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "Bali, Laba Sari"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "Indonesia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_02095758", "GABI_02008560")] <- "IDN"

# Castelnuovo in Montenegro instead of Serbia/Croatia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Castelnuovo"), replace = F)] <- "Montenegro"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Castelnuovo"), replace = F)] <- "MNE"

# Sona Mpungu Forest in Democratic Republic of the Congo instead of Republic of the Congo
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sona Mpungu Forest"), replace = F)] <- "Democratic Republic of the Congo"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sona Mpungu Forest"), replace = F)] <- "COD"

# Bangassou in Central African Republic instead of Equatorial Guinea
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bangassou"), replace = F)] <- "Central African Republic"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bangassou"), replace = F)] <- "CAF"

# Fond St Jacques in Martinique/France instead of Saint Lucia
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Fond St Jacques"), replace = F)] <- "Martinique"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Fond St Jacques"), replace = F)] <- "France"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Fond St Jacques"), replace = F)] <- "FRA"

# American Samoa instead of Samoa
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -170.70), replace = F)] <- "Anu'u Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -170.70), replace = F)] <- "American Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.295) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -170.70), replace = F)] <- "ASM"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -169.45), replace = F)] <- "Ta'u Island"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -169.45), replace = F)] <- "American Samoa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Latitude_initial == -14.23) & (Biogeographic_database_Ponerinae_flags$Longitude_initial == -169.45), replace = F)] <- "ASM"

# Somalia to Somaliland: [Ouarsangueli] Warsangali, Ainabo, Shimbir Beris, Ceelbuh, Ina Dhakool
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID == "GABI_00835211"), replace = F)] <- "Ceelbuh"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = (Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00833795", "GABI_00833796")), replace = F)] <- "Ina Dhakool"
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality %in% c("[Ouarsangueli] Warsangali", "Ainabo", "Shimbir Beris", "Ceelbuh", "Ina Dhakool")), replace = F)] <- "Somaliland"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality %in% c("[Ouarsangueli] Warsangali", "Ainabo", "Shimbir Beris", "Ceelbuh", "Ina Dhakool")), replace = F)] <- "SOL"

# Lusaka in Zambia instead of South Africa
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Zambia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lusaka") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "ZMB"

# Eastern Equatoria, Lotti Forest, SW Slope Achali Mts. in South Sudan instead of Sudan
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Eastern Equatoria") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "South Sudan"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Eastern Equatoria") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "SSD"

# West Caprivi Park, Manywa River in Namibia instead of Zambia
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "West Caprivi Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Namibia"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "West Caprivi Park") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "NAM"

# From Lochiel to Mbabane. Point on South African side of the border.
Biogeographic_database_Ponerinae_flags$Country_ISO3_name[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01952400")] <- "South Africa"
Biogeographic_database_Ponerinae_flags$Country_ISO3_code[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_01952400")] <- "ZAF"


## 2.6.2.4/ True positive: False coordinates to remove (turn into NA) ####

# Check if locality is credible...

# DRC, Haut-Ubangi => Wrong coordinates in Congo. To draw at random in the adm1.
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- F
Biogeographic_database_Ponerinae_flags$adm1[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Haut Ubangi") & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- "Nord Ubangi"

# Fiji => Wrong rounded coordinates falling in the Pacific Ocean # Draw random coordinates on the Fiji islands
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- NA
Biogeographic_database_Ponerinae_flags$valid_coordinates[replace_na(data = (Biogeographic_database_Ponerinae_flags$Longitude_initial == 179) & (Biogeographic_database_Ponerinae_flags$Latitude_initial == -18) & (Biogeographic_database_Ponerinae_flags$flag_type == "countries"), replace = F)] <- F

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



## 2.6.3/ Curate records falling in water areas ####

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "seas", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

## 2.6.3.1/ False positive: case of records flagged as at seas, but that are actually on lands (tiny islands or deltas) ####

# Brazil, Ilha dos Porcos, Amazon Delta
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00753526")] <- TRUE

# Costa Rica, Osa, Sierpe, Isla del Caño
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Osa, Sierpe") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Osa, Sierpe, Isla del Caño"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Isla del Caño") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Isla de Caño") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Italy, Archipelago Pontino, Ventotene Island
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Ventotene") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Japan, Kagoshima Prefecture, Minamisatsuma, Kamino Island
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00619696", "GABI_00617902")] <- TRUE
# Japan, Okinawa Prefecture, Iou-torishima Island
Biogeographic_database_Ponerinae_flags$.sea[Biogeographic_database_Ponerinae_flags$GABI_accession_ID %in% c("GABI_00349103")] <- TRUE

# Palau, Sonsorol, Sonsorol island
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$Latitude_dec < 5.34) & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Sonsorol, Sonsorol Island"
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$Latitude_dec > 5.34) & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- "Sonsorol, Fanna Island"
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm1, pattern = "Sonsorol") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Seychelles, Aride Island
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Aride") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# Spain, Galicia: Isla de Cíes and Isla de Ons
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$adm2, pattern = "Isla de Ons") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Islas Cíes") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

# USA, Keys Islands
Biogeographic_database_Ponerinae_flags$.sea[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = " Key") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- TRUE

## Convert back to logical the "countries" and "sea" flags
Biogeographic_database_Ponerinae_flags$.sea <- as.logical(Biogeographic_database_Ponerinae_flags$.sea)


##  2.6.3.2/ True positive: Slightly erroneous coordinates to correct (Fall in internal lakes, close to land, just imprecision) ####

# Brazil, Rio Grande do Sul, Parque Estadual de Itapua, Viamao => Wrong coordinates falling WE in Lake/Lagoa dos Patos
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Itapua") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -51.026
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Itapua") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -30.347

# Colombia, Magdalena, Cienaga Grande de Sta Maria => Wrong coordinates in the middle of the Lake
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cienaga Grande") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -74.324

# Greece, Evrostina => Wrong coordinates falling into the Corinth Gulf
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Evrostina") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- 22.396
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Evrostina") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- 38.073

# Portugal, Leça da Palmeira => Wrong Longitude, falling West in the Atlantic Ocean. -8.701389 instead of -8.761389
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Leça da Palmeira") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -8.701389

# Seychelles, Aride Island => Wrong rounded Latitude, falling North in the Pacific Ocean. -4.212 instead of -4.200
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Aride") & (Biogeographic_database_Ponerinae_flags$flag_type == "seas"), replace = F)] <- -4.212


Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cienaga")), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Cienaga")), replace = F)]

# Need to update the geometry to include new coordinates
# Rerun flag check to ensure changes have fixed them
# Modify again the false positives
# Need to check if they have a sea flag, or any other flag, too to change their flag

## 2.6.4/ Curate records in around country centroids ####

View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "centroids", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# Check if locality is available to retrieve true coordinates
# Check if uncertainty is too large to be kept (may be okay to keep those for small islands or countries)

## 2.6.4.1/ False positive: Records in localities that are close to country centroids, or centroids of small islands so it does not matter ####

# Cook Islands, Rarotonga
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Rarotonga"), replace = F)] <- TRUE

# Martinique
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country == "Martinique") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Haiti
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country == "Haiti") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Christmas Island
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country == "Christmas Island") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Santa Lucia, Praslin
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Praslin"), replace = F)] <- TRUE

# Saint Martin, Loterie Farm
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Loterie Farm"), replace = F)] <- TRUE

# Samoa Island
Biogeographic_database_Ponerinae_flags$.cen[(Biogeographic_database_Ponerinae_flags$Country == "Samoa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- TRUE

# Trinidad and Tobago, Lopinot Complex
Biogeographic_database_Ponerinae_flags$.cen[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Lopinot complex"), replace = F)] <- TRUE


## 2.6.4.2/ True positive: Records that have wrong coordinates according to their locality => To adjust ####

# Cameroun, Mount Cameroun above Buea
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 9.214
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Buea") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 4.170

# Czech Republic, Radotínské údolí
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Radotínské") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 14.314
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Radotínské") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 49.997

# Equatorial Guinea, Bioko
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bioko") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 8.70
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bioko") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 3.50

# Ethiopia, Gughé highlands, Bonghé Valley
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bonghé") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 37.35
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bonghé") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 6.06
# Ethiopia, Dire Dawa
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Dire-Dawa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 41.857
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Dire-Dawa") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 9.608

# Fiji, Sigatoka
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sigatoka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 177.553
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Sigatoka") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -18.054
# Fiji, Kings Road
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kings Rd") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 178.437
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kings Rd") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -17.783

# Malawi, Mlanje
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 35.508
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -16.024
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje_2000") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 35.568
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (Biogeographic_database_Ponerinae_flags$Locality_code == "Mlanje_2000") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -15.976

# Mexico, Peñon del Marquis
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Peñon del Marquis") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -99.017
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Peñon del Marquis") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 19.374
# Mexico, Pedrigales
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pedrigales") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -99.207
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pedrigales") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 19.322

# Mozambique, Amatongas Forest
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Amatongas") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 33.768
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Amatongas") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -19.182

# Myanmar, Bhamo
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bhamò") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 97.234
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Bhamò") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 24.262

# Taiwan, Akau
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.473
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 22.713
# Taiwan, Kosempo
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Jiaxian District [Kosempo, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.614
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Kosempo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.121
# Taiwan, Pilam
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Taitung District [Pilam, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.126
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Pilam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 22.759
# Taiwan, Taihorin
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Dalin [Taihorin, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.455
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihorin") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.604
# Taiwan, Taihoku/Taipeh
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Taipei [Taihoku, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.538
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taihoku") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 25.038
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taipeh") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 121.538
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Taipeh") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 25.038
# Taiwan, Suisharyo
Biogeographic_database_Ponerinae_flags$Locality[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- "Yuchi [Suisharyo, Formosa]"
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 120.920
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Suisharyo") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 23.872

# Tanzania, Tanganyika
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tanganyika") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 29.958
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Tanganyika") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -6.246

# Uganda, Busnia
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Busnia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 34.089
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Busnia") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 0.467
# uganda, Mabira Forest
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mabira Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 33.009
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Mabira Forest") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 0.459

# Zimbabwe, Matsuma Dam
Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matsuma Dam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- 26.282
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Matsuma Dam") & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids"), replace = F)] <- -18.733


Biogeographic_database_Ponerinae_flags$Longitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau")), replace = F)]
Biogeographic_database_Ponerinae_flags$Latitude_dec[replace_na(data = (str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Akau")), replace = F)]

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
Biogeographic_database_Ponerinae_flags$adm1[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "Solomon Islands") & replace_na(data = str_detect(string = Biogeographic_database_Ponerinae_flags$Locality, pattern = "Isabel"), replace = F) & (Biogeographic_database_Ponerinae_flags$flag_type == "centroids")] <- "Isabel Province"
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


## 2.6.5/ Curate records in capitals ####

# Look at locality names to check if nothing suspicious that could reflect the location of a specimen in collection rather than from the field.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "capitals", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# All seems fine, so mark them as false positive
# Biogeographic_database_Ponerinae_flags$.cap <- TRUE

# Update flag types
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & !Biogeographic_database_Ponerinae_flags$.inst] <- "institutions"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & Biogeographic_database_Ponerinae_flags$.inst & !Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.cap & Biogeographic_database_Ponerinae_flags$.inst & Biogeographic_database_Ponerinae_flags$.otl] <- "OK"

table(Biogeographic_database_Ponerinae_flags$flag_type)

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


## 2.6.6/ Curate records in around institutions ####

# Look at locality names to check if nothing suspicious that could reflect the location of a specimen in collection rather than from the field.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "institutions", ] %>% arrange(Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

# All seems fine, so mark them as false positive
# Biogeographic_database_Ponerinae_flags$.inst <- TRUE

# Update flag types
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst & !Biogeographic_database_Ponerinae_flags$.otl] <- "outliers"
Biogeographic_database_Ponerinae_flags$flag_type[!Biogeographic_database_Ponerinae_flags$.inst & Biogeographic_database_Ponerinae_flags$.otl] <- "OK"

table(Biogeographic_database_Ponerinae_flags$flag_type)

## Save sf output from CoordinateCleaner
saveRDS(Biogeographic_database_Ponerinae_flags, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_flags.rds")


## 2.6.7/ Curate outliers based on distance ####

## 2.6.7.1/ True positive: False coordinates to remove (turn into NA) ####

# Look at locality names to check if outliers are credible or not.
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$flag_type == "outliers", ] %>% arrange(Current_name, Country_ISO3_name, adm1, adm2, Locality, GABI_accession_ID, Specimen_code))

taxa_with_outliers_dist <- unique(Biogeographic_database_Ponerinae_flags$Current_name[Biogeographic_database_Ponerinae_flags$flag_type == "outliers"])
taxa_with_outliers_dist <- taxa_with_outliers_dist[order(taxa_with_outliers_dist)]
length(taxa_with_outliers_dist) # 207 taxa with outliers !

# Create Green/Red palette
Green_col <- RColorBrewer::brewer.pal(n = 6, "Greens")[5]
Red_col <- RColorBrewer::brewer.pal(n = 6, "Reds")[5]
pal_GR <- c(Green_col, Red_col)

# Taxa to test
Focal_taxon <- "Anochetus_afr01"
Focal_taxon <- "Anochetus_afr08"
Focal_taxon <- "Anochetus_africanus"
Focal_taxon <- "Anochetus_bequaerti"
Focal_taxon <- "Anochetus_bispinosus"
Focal_taxon <- "Anochetus_emarginatus"
Focal_taxon <- "Anochetus_fuliginosus"
Focal_taxon <- "Anochetus_ghilianii"
Focal_taxon <- "Anochetus_graeffei"
Focal_taxon <- "Anochetus_isolatus"
Focal_taxon <- "Anochetus_maynei"
Focal_taxon <- "Anochetus_obscuratus"
Focal_taxon <- "Anochetus_orchidicola"
Focal_taxon <- "Anochetus_pubescens"
Focal_taxon <- "Anochetus_punctaticeps"
Focal_taxon <- "Anochetus_risii"
Focal_taxon <- "Cryptopone_ochracea"

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


## List of taxa that looks fine
# Anochetus_bispinosus

## List of dubious taxa to plot
# Anochetus_afr01 => Mostly in Cameroon, Gabon, Congo. Outliers in Angola and Mozambique
# Anochetus_afr08 => Mostly in Cameroon, Gabon, Congo. Outliers in Angola and Mozambique
# Anochetus_africanus => Mostly Equatorial. Some records in Mozambique. Outliers in South Africa.
# Anochetus_altisquamis => Mostly South of Brazil. Outlier is alone in North of Brazil
# Anochetus_bequaerti => Equatorial Africa + SE Africa. Outlier in Ghana
# Anochetus_emarginatus => North of South America. One outlier in Caera that looks fine. One outlier in Cancun, Mexico that looks dubious.
# Anochetus_fuliginosus => Equatorial Africa + in South Africa (weird) + Outlier in Liberia
# Anochetus_isolatus => In PNG, Salomon Islands, but also in the Philippines. Outlier in Palau Seram in Indonesia
# Anochetus_maynei => West Africa. Outlier in Kenya
# Anochetus_obscuratus => Equatorial Africa + East Africa + Outlier in Ivory Coast
# Anochetus_orchidicola => Central America + Outlier in Mexico
# Anochetus_pubescens => In Comoros. 3 outliers on Africa mainland in Kenya, Zimbabwe and South Africa
# Anochetus_punctaticeps => In South Africa. 3 outliers in RDC, Cameroon and Ivory Coast

## List of issues across continents that could affect Biogeographic inferences
# Anochetus_ghilianii => In Morocco and Spain. Outlier in Brazil
# Anochetus_graeffei => Mostly in Indo-Malay and Australasia. Weird records in Seychelles, and Brazil. Outlier in eSwatini. 
# Anochetus_risii => Indo-Malay. Outlier in 

## Loop across taxa to plot maps

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


# Start with screenshots and run interactive version only for problematic cases?


### 2.7/ Curate outliers found in dubious bioregions ####

## 2.7.1/ Assign bioregion to each entry ####

# Load country_sf with metadata
Countries_NE_sf <- readRDS(file = "./input_data/NaturalEarth_maps/Countries_NE_sf.rds")

# Export metadata to be filled manually for Bioregions
Countries_NE_sf_metadata <- st_drop_geometry(Countries_NE_sf[, c("SUBUNIT", "ISO_A3")])
Countries_NE_sf_metadata$Bioregion <- NA
write.xlsx(x = Countries_NE_sf_metadata, file = "./input_data/Biogeographic_data/Countries_NE_sf_metadata.xlsx")

# Load metadata with bioregion assignment
Countries_NE_sf_metadata <- openxlsx::read.xlsx(xlsxFile = "./input_data/Biogeographic_data/Countries_NE_sf_metadata.xlsx")

# Assign Bioregion to entries
Biogeographic_database_Ponerinae_flags$Bioregion <- Countries_NE_sf_metadata$Bioregion[match(Biogeographic_database_Ponerinae_flags$Country_ISO3_code, Countries_NE_sf_metadata$ISO_A3)]

# Deal with the particular cases of countries overlapping bioregions

# France: Mayotte, Reunion, Guyane, Martinique, Guadeloupe
View(Biogeographic_database_Ponerinae_flags[Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France", ] %>% arrange(Country))
Biogeographic_database_Ponerinae_flags$Bioregion[(Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Palearctic"
Biogeographic_database_Ponerinae_flags$Bioregion[str_detect(string = Biogeographic_database_Ponerinae_flags$Country, pattern =  "Guadeloupe") & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country ==  "Martinique", replace = F) | replace_na(data = Biogeographic_database_Ponerinae_flags$adm1 ==  "Martinique", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country ==  "Mayotte", replace = F) | replace_na(data = Biogeographic_database_Ponerinae_flags$adm1 ==  "Mayotte", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country ==  "Reunion", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country ==  "Reunion", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Afrotropics"
Biogeographic_database_Ponerinae_flags$Bioregion[(replace_na(data = Biogeographic_database_Ponerinae_flags$Country ==  "French Guiana", replace = F)) & (Biogeographic_database_Ponerinae_flags$Country_ISO3_name == "France")] <- "Neotropics"

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
Biogeographic_database_Ponerinae_flags$Bioregion[replace_na(data = Biogeographic_database_Ponerinae_flags$Country == "Bonaire, Sint Eustatius and Saba", replace = F)] <- "Neotropics"

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

### 2.7.3/ Plot distribution maps of suspicious taxa

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




### Reintegrate NA entries dataset in the full curated dataset

# Extract entries with new NA coords (assign to NA becase they were flagged as error)
Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun <- cc_val(x = Biogeographic_database_Ponerinae_flags,
                                                                         lon = "Longitude_dec", lat = "Latitude_dec", 
                                                                         value = "flagged", verbose = T)

table(Biogeographic_database_Ponerinae_flags$valid_coordinates)
table(Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun)
View(x = Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun, ])

# Extract valid coordinates to rerun tests
Biogeographic_database_Ponerinae_flags_rerun <- Biogeographic_database_Ponerinae_flags %>% 
  filter(Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun) %>% 
  select(-c(".val", ".cap", ".cen", ".sea", ".con", ".otl", ".inst", ".dpl", ".summary")) %>% 
  st_drop_geometry() # %>% 
  # select(-c("Latitude_copy", "Longitude_copy"))
class(Biogeographic_database_Ponerinae_flags_rerun) <- "data.frame"


Biogeographic_database_Ponerinae_NA_coords_new <- st_drop_geometry(Biogeographic_database_Ponerinae_flags[!Biogeographic_database_Ponerinae_flags$valid_coordinates_rerun, ])
Biogeographic_database_Ponerinae_NA_coords <- rbind(Biogeographic_database_Ponerinae_NA_coords, Biogeographic_database_Ponerinae_NA_coords_new[, names(Biogeographic_database_Ponerinae_NA_coords)])

## Save sf dataframe with coordinates with NA
saveRDS(Biogeographic_database_Ponerinae_NA_coords, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_NA_coords.rds")



### 2.8/ Plot interactive maps with mapview ####

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

## 2.8.1/ Prepare color palette ####

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

## 2.8.2/ Plot per species ####

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


## 2.8.3/ Plot per flag types ####

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



# Check dubious records with CoordinateCleaner. Remove erroneous. Flag dubious (be careful if it is the last data for a taxa)
# Plot maps of dubious stuff. See my script from Ithomiini

## Do an interactive map per taxa using a loop?

# Share all dubious with AoW team


# Flag especially records across continents that represent a low count/percentage => could cause mistakes in biogeographic analyses


# Number of totally cleaned occurrences with no flag
sum(flags$.summary)


# Plot flagged or not occurrences
?plot.spatialvalid

plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
# Need to go into the function to edit the fact they use ugly shape
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude", details = T, pts_size = 2)

# Remove flagged occurrences
dat_cl <- dat_cl[flags$.summary, ]

### 5/ Map cleaned dataset vs. original one #

world.inp <- map_data("world")

ggplot() + geom_map(data = world.inp, map = world.inp, 
                    aes(x = long, y = lat, map_id = region), fill = "grey80") + 
  xlim(min(dat$decimalLongitude, na.rm = T), max(dat$decimalLongitude, na.rm = T)) + 
  ylim(min(dat$decimalLatitude, na.rm = T), max(dat$decimalLatitude, na.rm = T)) + 
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 1) + 
  geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkgreen", size = 1) + 
  coord_fixed() + theme_bw() + theme(axis.title = element_blank())

map_plot(dat_cl)

### 6/ Save the result

dat_cl_with_fossils <- dat_cl
dat_cl_without_fossils <- dat_cl

save(dat_cl_with_fossils, file = "./data/Xenarthra/Xenarthra_America_GBIF_full_cleaned_with_fossils.RData")
write_csv(dat_cl_with_fossils, "./data/Xenarthra/Xenarthra_America_GBIF_full_cleaned_with_fossils.csv")

save(dat_cl_without_fossils, file = "./data/Xenarthra/Xenarthra_America_GBIF_full_cleaned_without_fossils.RData")
write_csv(dat_cl_without_fossils, "./data/Xenarthra/Xenarthra_America_GBIF_full_cleaned_without_fossils.csv")

### Flag for differences between initial coordinates in GABI and AntWeb, and current coordinates
### Flag for randomly drawn coordinates within country/adm1/2

### Remove dubious records from initial databases too (GABI and AntWeb), not only the merged one

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

# Create table of GABI_adm1 (Bentity2) matches in NE_adm1
# Create matching table based on spatial overlay
# Manually curate GABI_adm1 with multiple matches (may be due to border effects)
# Random region will be picked based on random coordinates drawn within the GABI region
# Create table of AW_adm1 (AntWeb/Geographic Names Server) matches in NE_adm1
# Create matching table based on name similarity
# If multiple match, random region will be picked and random coordinates will be drawn within this new NE region

# Simulate only one per duplicate taxa/adm and no need if already present with coordinates

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

plot(map_NE_adm1_sf[, "adm1_code"])   # Sovereign administrative countries

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


##### Flag duplicates #####

## Create a database without duplicates. Keep AntWeb over GABI when possible. Run taxa-summary stat on without duplicates data

##### Extract specimen database for taxa-level information ####

# Add nb of biogeographic records


##### Create the Genus-level summary table #####

# Create a summary table at Genus-level of nb of valid species, nb of valid species included in our phylogeny, nb of morphospecies in our phylogeny, other (not curated) morphospecies from AntWeb
# Add nb of biogeographic records
