##### Script 19: Generate Voucher table for NCBI #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Generate Voucher table needed for raw sequence upload in NCBI SRA - BioSample (only newly sequenced data)
# Generate voucher metadata table for Supplementary Data 1 (includes all taxa in the phylogeny)
  # Includes NCBI SRA - BioSample for raw data, and NCBI GenBank accession ID for UCE contigs and COI

###

### Inputs

# AntWeb specimen database accessed on 30 May 2024
# Phylogeny with metadata including voucher codes
# Source table tracking source of extracts
# Curated occurrence database

###

### Sources

# AntWeb. Version 8.101. California Academy of Science, online at https://www.antweb.org. Accessed on 30 May 2024.
# NCBI - BioSample: Template for submission: https://submit.ncbi.nlm.nih.gov/biosample/template/

###

### Outputs

# Voucher table following template provided by NCBI - BioSample (missing proper biomaterial information)

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(xlsx)      # Need the Java Development Kit (JDK) installed
library(openxlsx)  # Use Rccp. No need of Java
library(phytools)
library(ape)
library(treeio)

### 1.2/ Load databases ####

# Load (curated) AntWeb databases
AntWeb_database_curated <- readRDS(file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")
AntWeb_database <- read_excel("./input_data/AntWeb_data/AntWeb_database_Ponerinae_2024_05_30.xlsx")

# Load updated df for Metadata of samples in the phylogeny
Phylogeny_sample_data_789t <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

# Load curated occurrence database
Biogeographic_database_Ponerinae_curated <- readRDS(file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae_curated.rds")


##### 2/ Generate template and populate with already available metadata ####

NCBI_BioSample_Voucher_table <- Phylogeny_sample_data_789t %>% 
  rename(organism = Current_name,
         specimen_voucher = Specimen_Code) %>%
  mutate(sample_name = paste0(organism, "_", Extraction_Code, "_", specimen_voucher),
         bioproject_accession = NA,
         isolate = "not applicable",
         breed = "not applicable",
         host = "not applicable",
         isolation_source = "not applicable",
         collection_date = NA, 
         tissue = "entire specimen (destructively or non-destructively)",
         sex = NA,
         `infra specific rank` = NA,
         `infra specific name` = NA,
         description = NA)  %>% 
  dplyr::select(Current_status, Subspecies, # To keep to inform infraspecific ranks
                sample_name, bioproject_accession, organism, isolate,	breed,	host,	isolation_source,	collection_date, tissue, sex,	specimen_voucher,	`infra specific rank`, `infra specific name`, description)

# View(NCBI_BioSample_Voucher_table)

# Deal with morphospecies
NCBI_BioSample_Voucher_table$`infra specific rank` <- NCBI_BioSample_Voucher_table$Current_status
table(NCBI_BioSample_Voucher_table$`infra specific rank`)
NCBI_BioSample_Voucher_table$`infra specific rank`[NCBI_BioSample_Voucher_table$`infra specific rank` != "morphotaxon"] <- NA
table(NCBI_BioSample_Voucher_table$`infra specific rank`)

# Deal with subspecies
Subspecies_ID <- which(NCBI_BioSample_Voucher_table$Subspecies)
NCBI_BioSample_Voucher_table$`infra specific name`[Subspecies_ID] <- str_remove(string = NCBI_BioSample_Voucher_table$organism[Subspecies_ID], pattern = ".*_")
NCBI_BioSample_Voucher_table$organism[Subspecies_ID] <- str_remove(string = NCBI_BioSample_Voucher_table$organism[Subspecies_ID], pattern = paste(NCBI_BioSample_Voucher_table$`infra specific name`[Subspecies_ID], collapse = "|"))
NCBI_BioSample_Voucher_table$organism[Subspecies_ID] <- str_remove(string = NCBI_BioSample_Voucher_table$organism[Subspecies_ID], pattern = "_$")

# Remove fields
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  dplyr::select(-Current_status, -Subspecies)

# View(NCBI_BioSample_Voucher_table)


##### 3/ Complement metadata using AntWeb #####

AntWeb_database_to_merge <- AntWeb_database_curated %>% 
  mutate(specimen_voucher = str_to_upper(SpecimenCode)) %>%
  rename(biomaterial_provider = ownedby,
         dev_stage = life_stage,
         collected_by = collectedby,
         identified_by = determinedby,
         altitude = elevation,
         collection_date_temp = other) %>%
  dplyr::select(specimen_voucher, biomaterial_provider, dev_stage, collected_by, identified_by, altitude, collection_date_temp) %>%
  filter(specimen_voucher %in% NCBI_BioSample_Voucher_table$specimen_voucher)

# Remove "null"
AntWeb_database_to_merge$identified_by[AntWeb_database_to_merge$identified_by == "null"] <- NA
AntWeb_database_to_merge$altitude[AntWeb_database_to_merge$altitude == "null"] <- NA
AntWeb_database_to_merge$biomaterial_provider[AntWeb_database_to_merge$biomaterial_provider == "null"] <- NA

# View(AntWeb_database_to_merge)

# Merge with voucher data
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>%
  left_join(y = AntWeb_database_to_merge, by = join_by(specimen_voucher)) %>% 
  dplyr::select(sample_name, bioproject_accession, organism, isolate,	breed,	host,	isolation_source,	collection_date, collection_date_temp, tissue, altitude, biomaterial_provider, collected_by, dev_stage, identified_by, sex,	specimen_voucher,	`infra specific rank`, `infra specific name`, description)

# View(NCBI_BioSample_Voucher_table)

## Deal with collection data

# In other: <collrecorddate>17 Jun 2014</collrecorddate>
NCBI_BioSample_Voucher_table$collection_date_temp <- str_remove(string = NCBI_BioSample_Voucher_table$collection_date_temp, pattern = ".*\\<collrecorddate\\>")
NCBI_BioSample_Voucher_table$collection_date_temp <- str_remove(string = NCBI_BioSample_Voucher_table$collection_date_temp, pattern = "\\<\\/collrecorddate\\>.*")

collection_date_list <- str_split(string = NCBI_BioSample_Voucher_table$collection_date_temp, pattern = " ")
unlist(lapply(X = collection_date_list, FUN = length))
collection_date_list[unlist(lapply(X = collection_date_list, FUN = length)) != 3] <- NA
collection_date_list <- lapply(X = collection_date_list, FUN = paste0, collapse = "-")
collection_date_values <- unlist(collection_date_list)
collection_date_values[collection_date_values == "NA"] <- NA 
collection_date_values[nchar(collection_date_values) > 11] <- NA

NCBI_BioSample_Voucher_table$collection_date <- collection_date_values

# Remove temp fields
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>%
  dplyr::select(-collection_date_temp)

View(NCBI_BioSample_Voucher_table)


##### 4/ Complement metadata using Occurrence database #####

# geo_loc_name = Country: Localityname (without ":") (Use data from curated occurrences)
# lat_lon = "d[d.dddd] N|S d[dd.dddd] W|E", eg, 38.98 N 77.11 W. (Use data from curated occurrences)

Biogeographic_database_to_merge <- Biogeographic_database_Ponerinae_curated %>% 
  sf::st_drop_geometry() %>% 
  mutate(specimen_voucher = str_to_upper(Specimen_code)) %>%
  dplyr::select(specimen_voucher, Latitude_dec, Longitude_dec, Country_ISO3_name, adm1, adm2, Locality) %>%
  filter(specimen_voucher %in% NCBI_BioSample_Voucher_table$specimen_voucher)

# View(Biogeographic_database_to_merge)

## Get coordinates in proper format

# No need to convert in degrees! Only to cut down to two decimals

Biogeographic_database_to_merge$lat_lon <- NA

decimal_to_Biosample_format <- function (lat, long, nb_decimals = 2)
{
  factor <- 10^nb_decimals 
    
  abs_lat <- abs(lat)
  lat_rounded <- round(abs_lat * factor, 0) / factor
  if (abs_lat == lat) { lat_card <- "N" } else { lat_card <- "S" } # Get cardinal

  abs_long <- abs(long)
  long_rounded <- round(abs_long * factor, 0) / factor
  if (abs_long == long) { long_card <- "E" } else { long_card <- "W" } # Get cardinal

  # Aggregate result
  lat_long_decimals <- paste0(lat_rounded, " ", lat_card, " ", long_rounded, " ", long_card)
  return(lat_long_decimals)
}

# decimal_to_Biosample_format(lat = 4.08851, long = 52.67922, nb_decimals = 4)


# Loop per entry
Biogeographic_database_to_merge$lat_lon <- NA
for (i in 1:nrow(Biogeographic_database_to_merge))
{
  Biogeographic_database_to_merge$lat_lon[i] <- decimal_to_Biosample_format(lat = Biogeographic_database_to_merge$Latitude_dec[i],
                                                                           long = Biogeographic_database_to_merge$Longitude_dec[i])
}

# # Convert degrees to long/lat
# decimal_to_degrees_minutes <- function (lat, long)
# {
#   abs_lat <- abs(lat)
#   lat_degrees <- floor(abs_lat) # Get degrees
#   lat_minutes <- round((abs_lat - lat_degrees) * 60, 0) # Get rounded minutes
#   if (abs_lat == lat) { lat_card <- "N" } else { lat_card <- "S" } # Get cardinal
#   
#   abs_long <- abs(long)
#   long_degrees <- floor(abs_long) # Get degrees
#   long_minutes <- round((abs_long - long_degrees) * 60, 0) # Get rounded minutes
#   if (abs_long == long) { long_card <- "E" } else { long_card <- "W" } # Get cardinal
#   
#   # Aggregate result
#   lat_long_deg <- paste0(lat_degrees, ".", lat_minutes, " ", lat_card, " ", long_degrees, ".", long_minutes, " ", long_card)
#   return(lat_long_deg)
# }
# 
# # Loop per entry
# Biogeographic_database_to_merge$lat_lon <- NA
# for (i in 1:nrow(Biogeographic_database_to_merge))
# {
#   Biogeographic_database_to_merge$lat_lon[i] <- decimal_to_degrees_minutes(lat = Biogeographic_database_to_merge$Latitude_dec[i],
#                                                                            long = Biogeographic_database_to_merge$Longitude_dec[i])
# }

## Get locality info in proper format

# Replace ":" to avoid issues
Biogeographic_database_to_merge$Country_ISO3_name <- str_replace_all(string = Biogeographic_database_to_merge$Country_ISO3_name, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$adm1 <- str_replace_all(string = Biogeographic_database_to_merge$adm1, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$adm2 <- str_replace_all(string = Biogeographic_database_to_merge$adm2, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$Locality <- str_replace_all(string = Biogeographic_database_to_merge$Locality, pattern = ":", replacement = " ")

# Merge in one field
Biogeographic_database_to_merge$geo_loc_name <- paste0(Biogeographic_database_to_merge$Country_ISO3_name, ": ", Biogeographic_database_to_merge$adm1, ", ", Biogeographic_database_to_merge$adm2, ", ", Biogeographic_database_to_merge$Locality)
  
# View(Biogeographic_database_to_merge)

# Merge with voucher data
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>%
  left_join(y = Biogeographic_database_to_merge[, c("specimen_voucher", "geo_loc_name", "lat_lon")], by = join_by(specimen_voucher)) %>% 
  dplyr::select(sample_name, bioproject_accession, organism, isolate,	breed,	host,	isolation_source,	collection_date, geo_loc_name, tissue, altitude, biomaterial_provider, collected_by, dev_stage, identified_by, lat_lon,	sex,	specimen_voucher,	`infra specific rank`, `infra specific name`, description) %>%
  arrange(sample_name)

View(NCBI_BioSample_Voucher_table)

## Change Honk Kong S.A.R. for Honk Kong (BioSample rule)
NCBI_BioSample_Voucher_table$geo_loc_name <- str_replace(string = NCBI_BioSample_Voucher_table$geo_loc_name, pattern = "Hong Kong S.A.R.", replacement = "Hong Kong")

## Update location of CASENT0777939 specimen in Mozambique
NCBI_BioSample_Voucher_table$lat_lon[NCBI_BioSample_Voucher_table$specimen_voucher == "CASENT0777939"] <- "12.82 S 39.70 E"
NCBI_BioSample_Voucher_table$geo_loc_name[NCBI_BioSample_Voucher_table$specimen_voucher == "CASENT0777939"] <- "Mozambique: Cabo Delgado, Parque Nacional Quirimbas, Taratibu"


##### 5/ Deal with specimens with no match in databases ####

# Not in my curated Antweb database (probably introduced)
Not_in_curated_AntWeb <- NCBI_BioSample_Voucher_table$specimen_voucher[!(NCBI_BioSample_Voucher_table$specimen_voucher %in% AntWeb_database_to_merge$specimen_voucher)]
Not_in_curated_AntWeb_ID <- which(!(NCBI_BioSample_Voucher_table$specimen_voucher %in% AntWeb_database_to_merge$specimen_voucher))
Not_in_curated_AntWeb

# Not in my curated occurrence database (probably introduced + location errors)
Not_in_Occ_database <- NCBI_BioSample_Voucher_table$specimen_voucher[!(NCBI_BioSample_Voucher_table$specimen_voucher %in% Biogeographic_database_to_merge$specimen_voucher)]
Not_in_Occ_database_ID <- which(!(NCBI_BioSample_Voucher_table$specimen_voucher %in% Biogeographic_database_to_merge$specimen_voucher))
Not_in_Occ_database

# All included in the specimen missing from occurrence data
Not_in_curated_AntWeb %in% Not_in_Occ_database

# Check if found in the non-curated AntWeb database
Not_in_curated_AntWeb %in% str_to_upper(AntWeb_database$SpecimenCode)
Not_in_Occ_database %in% str_to_upper(AntWeb_database$SpecimenCode)

# Get metdata from the non-curated AntWeb database
AntWeb_database_complement_to_merge <- AntWeb_database %>% 
  mutate(specimen_voucher = str_to_upper(SpecimenCode)) %>%
  rename(biomaterial_provider = ownedby,
         dev_stage = life_stage,
         collected_by = collectedby,
         identified_by = determinedby,
         altitude = elevation,
         collection_date_temp = other) %>%
  dplyr::select(specimen_voucher, biomaterial_provider, dev_stage, collected_by, identified_by, altitude, collection_date_temp) %>%
  filter(specimen_voucher %in% Not_in_curated_AntWeb)
# Check order
Not_in_curated_AntWeb == AntWeb_database_complement_to_merge$specimen_voucher

# Remove "null"
AntWeb_database_complement_to_merge$identified_by[AntWeb_database_complement_to_merge$identified_by == "null"] <- NA
AntWeb_database_complement_to_merge$altitude[AntWeb_database_complement_to_merge$altitude == "null"] <- NA
AntWeb_database_complement_to_merge$biomaterial_provider[AntWeb_database_complement_to_merge$biomaterial_provider == "null"] <- NA

# View(AntWeb_database_complement_to_merge)

# Merge with voucher data
NCBI_BioSample_Voucher_table$biomaterial_provider[Not_in_curated_AntWeb_ID] <- AntWeb_database_complement_to_merge$biomaterial_provider
NCBI_BioSample_Voucher_table$dev_stage[Not_in_curated_AntWeb_ID] <- AntWeb_database_complement_to_merge$dev_stage
NCBI_BioSample_Voucher_table$collected_by[Not_in_curated_AntWeb_ID] <- AntWeb_database_complement_to_merge$collected_by
NCBI_BioSample_Voucher_table$identified_by[Not_in_curated_AntWeb_ID] <- AntWeb_database_complement_to_merge$identified_by
NCBI_BioSample_Voucher_table$altitude[Not_in_curated_AntWeb_ID] <- AntWeb_database_complement_to_merge$altitude

# View(NCBI_BioSample_Voucher_table)

## Deal with collection data

# In other: <collrecorddate>17 Jun 2014</collrecorddate>
collection_date_temp <- str_remove(string = AntWeb_database_complement_to_merge$collection_date_temp, pattern = ".*\\<collrecorddate\\>")
collection_date_temp <- str_remove(string = collection_date_temp, pattern = "\\<\\/collrecorddate\\>.*")

collection_date_list <- str_split(string = collection_date_temp, pattern = " ")
unlist(lapply(X = collection_date_list, FUN = length))
collection_date_list[unlist(lapply(X = collection_date_list, FUN = length)) != 3] <- NA
collection_date_list <- lapply(X = collection_date_list, FUN = paste0, collapse = "-")
collection_date_values <- unlist(collection_date_list)
collection_date_values[collection_date_values == "NA"] <- NA 
collection_date_values[nchar(collection_date_values) > 11] <- NA

NCBI_BioSample_Voucher_table$collection_date[Not_in_curated_AntWeb_ID] <- collection_date_values

## Add biogeographic info

Biogeographic_database_complement_to_merge <- AntWeb_database %>% 
  mutate(specimen_voucher = str_to_upper(SpecimenCode)) %>%
  dplyr::select(specimen_voucher, decimal_latitude, decimal_longitude, country, adm1, adm2, localityname) %>%
  filter(specimen_voucher %in% Not_in_Occ_database)

# Remove "null"
Biogeographic_database_complement_to_merge$adm1[Biogeographic_database_complement_to_merge$adm1 == "null"] <- NA
Biogeographic_database_complement_to_merge$adm2[Biogeographic_database_complement_to_merge$adm2 == "null"] <- NA

# View(Biogeographic_database_complement_to_merge)

## Get coordinates in proper format

# Loop per entry
Biogeographic_database_complement_to_merge$lat_lon <- NA
for (i in 1:nrow(Biogeographic_database_complement_to_merge))
{
  Biogeographic_database_complement_to_merge$lat_lon[i] <- decimal_to_degrees_minutes(lat = Biogeographic_database_complement_to_merge$decimal_latitude[i],
                                                                              long = Biogeographic_database_complement_to_merge$decimal_longitude[i])
}

## Get locality info in proper format

# Replace ":" to avoid issues
Biogeographic_database_complement_to_merge$country <- str_replace_all(string = Biogeographic_database_complement_to_merge$country , pattern = ":", replacement = " ")
Biogeographic_database_complement_to_merge$adm1 <- str_replace_all(string = Biogeographic_database_complement_to_merge$adm1, pattern = ":", replacement = " ")
Biogeographic_database_complement_to_merge$adm2 <- str_replace_all(string = Biogeographic_database_complement_to_merge$adm2, pattern = ":", replacement = " ")
Biogeographic_database_complement_to_merge$localityname <- str_replace_all(string = Biogeographic_database_complement_to_merge$localityname, pattern = ":", replacement = " ")

# Merge in one field
Biogeographic_database_complement_to_merge$geo_loc_name <- paste0(Biogeographic_database_complement_to_merge$country, ": ", Biogeographic_database_complement_to_merge$adm1, ", ", Biogeographic_database_complement_to_merge$adm2, ", ", Biogeographic_database_complement_to_merge$localityname)

# View(Biogeographic_database_complement_to_merge)

# Merge with voucher data
NCBI_BioSample_Voucher_table$lat_lon[Not_in_Occ_database_ID] <- Biogeographic_database_complement_to_merge$lat_lon
NCBI_BioSample_Voucher_table$geo_loc_name[Not_in_Occ_database_ID] <- Biogeographic_database_complement_to_merge$geo_loc_name

View(NCBI_BioSample_Voucher_table)

##### 6/ Deal with the outgroups #####

# Load the 792t phylogeny with metadata
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata <- read.beast.newick(file = "./input_data/Phylogenies/ponerinae-792t-spruce-75p-iqtree-swscmerge-mfp_v2_v2.tre")

Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label[!(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label %in% NCBI_BioSample_Voucher_table$sample_name)]

## Add them manually!
write.xlsx(x = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

# Reload modified table
NCBI_BioSample_Voucher_table <- read.xlsx(xlsxFile = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

# Update error in manual dates
NCBI_BioSample_Voucher_table$collection_date[NCBI_BioSample_Voucher_table$sample_name == "Amblyopone_australis_D0872_CASENT0106229"] <- "13-Jan-1999"
NCBI_BioSample_Voucher_table$collection_date[NCBI_BioSample_Voucher_table$sample_name == "Paraponera_clavata_EX1573_CASENT0633292"] <- "8-Jun-2011"
NCBI_BioSample_Voucher_table$collection_date[NCBI_BioSample_Voucher_table$sample_name == "Proceratium_google_MAMI0434_CASENT0035028"] <- "12-Mar-2003"

##### 7/ Curate dev_stage and sex #####

# NCBI_BioSample_Voucher_table <- readRDS(file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

NCBI_BioSample_Voucher_table$caste_info <- NCBI_BioSample_Voucher_table$dev_stage
NCBI_BioSample_Voucher_table$dev_stage <- NA

# adult worker; female
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("Adult worker", "1 adult worker", "1 W", "1W", "1 worker", "1 worker (head detached)", "1 worker (missing abdomen)", "1w", "2 worker", "2 workers", "3 workers", "4 workers", "adult worker", "worker")] <- "adult worker"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("Adult worker", "1 adult worker", "1 W", "1W", "1 worker", "1 worker (head detached)", "1 worker (missing abdomen)", "1w", "2 worker", "2 workers", "3 workers", "4 workers", "adult worker", "worker")] <- "female"
  
# alate queen; female
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("1 alate queen", "1aq", "alate queen")] <- "alate queen"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("1 alate queen", "1aq", "alate queen")] <- "female"

# queen; female
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("queen", "1eq")] <- "queen"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("queen", "1eq")] <- "female"

# dealate queen; female
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("1 dealate queen", "1dq", "dealate queen")] <- "dealate queen"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("1 dealate queen", "1dq", "dealate queen")] <- "female"

# adult male; male
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("1 male", "male")] <- "adult male"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("1 male", "male")] <- "male"

# missing; missing
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("1 dealate queen, 1 worker")] <- "missing"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("1 dealate queen, 1 worker")] <- "missing"

# missing; missing
NCBI_BioSample_Voucher_table$dev_stage[NCBI_BioSample_Voucher_table$caste_info %in% c("null")] <- "missing"
NCBI_BioSample_Voucher_table$sex[NCBI_BioSample_Voucher_table$caste_info %in% c("null")] <- "missing"

NCBI_BioSample_Voucher_table$caste_info[is.na(NCBI_BioSample_Voucher_table$dev_stage)]
table(is.na(NCBI_BioSample_Voucher_table$dev_stage))
table(is.na(NCBI_BioSample_Voucher_table$sex))

# Remove temporary field
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  dplyr::select(-caste_info)

##### 8/ Add description field ####

NCBI_BioSample_Voucher_table$description <- paste0("The DNA voucher specimen ", NCBI_BioSample_Voucher_table$specimen_voucher," is deposited at ", NCBI_BioSample_Voucher_table$biomaterial_provider," with data available on www.antweb.org. The voucher is the same specimen as the one used for extraction, or is by default from the same nest/series/locality as the sequenced specimen.")

### Save/Export the Voucher table
saveRDS(object = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
write.xlsx(x = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

## Reorder by inversing for "Neoponera_metanotalis" vs. "Neoponera_metanotalis_2"

##### 9/ Filter to keep only newly sequenced data in the Biosample table ####

# Load Biosample voucher table
# NCBI_BioSample_Voucher_table <- readRDS(file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
NCBI_BioSample_Voucher_table <- openxlsx::read.xlsx(xlsxFile = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

# Save the non-filtered version to build the SD1 table with all taxa
SD1_Voucher_specimens_metadata <- NCBI_BioSample_Voucher_table
saveRDS(object = SD1_Voucher_specimens_metadata, file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

# NCBI_BioSample_Voucher_table <- readRDS(file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

# test <- SD1_Voucher_specimens_metadata %>% 
#   filter(Voucher_specimen_code %in% NCBI_BioSample_Voucher_table$specimen_voucher)
# 
# # Loop per entry
# test$lat_lon <- NA
# for (i in 1:nrow(test))
# {
#   test$lat_lon[i] <- decimal_to_Biosample_format(lat = test$Latitude_dec[i],
#                                                  long = test$Longitude_dec[i])
# }
# 
# 
# write.xlsx(x = test, file = "./input_data/Molecular_data/test.xlsx")

# Load Source table
AoW_Ponerinae_sequence_sources <- openxlsx::read.xlsx(xlsxFile = "./input_data/Molecular_data/AoW ponerine sequence source.xlsx")

table(AoW_Ponerinae_sequence_sources$status)

# Mismatch for outgroups (3) and subspecies (2)
AoW_Ponerinae_sequence_sources$Current_name[which(!(AoW_Ponerinae_sequence_sources$Current_name %in% NCBI_BioSample_Voucher_table$organism))]

AoW_Ponerinae_sequence_sources <- AoW_Ponerinae_sequence_sources %>% 
  arrange(Current_name)

# Fixed subspecies mismatch by using full_name
NCBI_BioSample_Voucher_table$Full_name <- paste0(NCBI_BioSample_Voucher_table$organism, "_", NCBI_BioSample_Voucher_table$infra.specific.name)
NCBI_BioSample_Voucher_table$Full_name <- str_remove(string = NCBI_BioSample_Voucher_table$Full_name, pattern = "_NA$")

table(AoW_Ponerinae_sequence_sources$Current_name %in% NCBI_BioSample_Voucher_table$Full_name)

# Mismatch order for "Neoponera_metanotalis" vs. "Neoponera_metanotalis_2"
AoW_Ponerinae_sequence_sources$Current_name[which(!(AoW_Ponerinae_sequence_sources$Current_name == NCBI_BioSample_Voucher_table$Full_name))]

NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  left_join(y = AoW_Ponerinae_sequence_sources[, c("Current_name", "status")], by = join_by("Full_name" == "Current_name"))

# Remove non-newly sequenced data
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  filter(status == "sequenced") %>%
  dplyr::select(-status, - Full_name)

# Remove "_" in organism names
NCBI_BioSample_Voucher_table$organism <- str_replace(string = NCBI_BioSample_Voucher_table$organism, pattern = "_", replacement = " ")

# Add missing in missing collection dates
NCBI_BioSample_Voucher_table$collection_date[is.na(NCBI_BioSample_Voucher_table$collection_date)] <- "missing"

# Detect non-ASCII characters
detect_non_ASCII_characters_from_df <- function (df)
{
  all_chars <- paste(unlist(df), collapse = "")
  all_non_ASCII_chars <- unique(unlist(str_extract_all(string = all_chars, pattern = "[^\t\n\r\x20-\x7E]")))
  cat(paste0(length(all_non_ASCII_chars), " unique non-ASCII characters detected: ", paste(all_non_ASCII_chars, collapse = " ")))
}
detect_non_ASCII_characters_from_df(df = NCBI_BioSample_Voucher_table)

# Replace non-ASCII characters
# Run str_replace on all columns

### Save/Export the Voucher table
saveRDS(object = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
write.xlsx(x = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

# Need to be pasted in the formatted version "_AOW_Ponerinae"


## Info
# `infra specific rank` = to tag morphospecies
# `infra specific name` = to name subspecies

test <- SD1_Voucher_specimens_metadata

# Loop per entry
test$lat_lon <- NA
for (i in 1:nrow(test))
{
  test$lat_lon[i] <- decimal_to_Biosample_format(lat = test$Latitude_dec[i],
                                                 long = test$Longitude_dec[i])
}


##### 10/ Create Supplementary data for all taxa ####

# Load Source table
# SD1_Voucher_specimens_metadata <- openxlsx::read.xlsx(xlsxFile = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.xlsx")
SD1_Voucher_specimens_metadata <- readRDS(file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

# Deal with dates
SD1_Voucher_specimens_metadata$collection_date
date_vector <- lubridate::as_date(x = SD1_Voucher_specimens_metadata$collection_date, origin = "1899-12-30")
formatted_date_vector <- format(date_vector, "%Y-%m-%d")
SD1_Voucher_specimens_metadata$collection_date <- formatted_date_vector

### 10.1/ Merge with Source table ####

# Load Source table
AoW_Ponerinae_sequence_sources <- openxlsx::read.xlsx(xlsxFile = "./input_data/Molecular_data/AoW ponerine sequence source.xlsx")

# Mismatch for subspecies (2)
AoW_Ponerinae_sequence_sources$Current_name[which(!(AoW_Ponerinae_sequence_sources$Current_name %in% SD1_Voucher_specimens_metadata$organism))]

AoW_Ponerinae_sequence_sources <- AoW_Ponerinae_sequence_sources %>% 
  arrange(Current_name)

# Fixed subspecies mismatch by using full_name
SD1_Voucher_specimens_metadata$Full_name <- paste0(SD1_Voucher_specimens_metadata$organism, "_", SD1_Voucher_specimens_metadata$infra.specific.name)
SD1_Voucher_specimens_metadata$Full_name <- str_remove(string = SD1_Voucher_specimens_metadata$Full_name, pattern = "_NA$")

table(AoW_Ponerinae_sequence_sources$Current_name %in% SD1_Voucher_specimens_metadata$Full_name)

# Mismatch order for "Neoponera_metanotalis" vs. "Neoponera_metanotalis_2"
AoW_Ponerinae_sequence_sources$Current_name[which(!(AoW_Ponerinae_sequence_sources$Current_name == SD1_Voucher_specimens_metadata$Full_name))]

# Merge and rename fields
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  left_join(y = AoW_Ponerinae_sequence_sources, by = join_by("Full_name" == "Current_name")) %>% 
  rename(Specimen_code_Sources = Specimen_code) %>% 
  rename(Extraction_code_Sources = Extraction_code) %>% 
  rename(Outgroup = outgroup) %>% 
  rename(Status_data = status) %>% 
  rename(Reference = ref) %>%
  rename(Taxa_name = Full_name) %>%
  rename(Voucher_specimen_code = specimen_voucher) %>%
  rename(Collection_date = collection_date) %>%
  rename(Locality = geo_loc_name) %>%
  rename(Caste = dev_stage) %>%
  rename(Sex = sex) %>%
  rename(Host_institution = biomaterial_provider) %>%
  rename(Collected_by = collected_by) %>%
  rename(Identified_by = identified_by) %>%
  rename(Status_taxa = infra.specific.rank) %>%
  rename(Description = description) %>%
  rename(BioSample_ID = Biosample) %>%
  rename(BioProject_ID = BioProject) %>%
  dplyr::select(Taxa_name, Status_taxa, Outgroup, Status_data, Voucher_specimen_code, Specimen_code_Sources, Extraction_code_Sources, note, BioProject_ID, BioSample_ID, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, lat_lon, Collected_by, Identified_by, Description)


### 10.2/ Merge with Macroevolution_taxa_database  ####

Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

# Extraction code to double check with Sources's codes
# Subspecies

SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  left_join(y = Ponerinae_Macroevolution_taxa_database[, c("Current_name", "Specimen_phylogeny_Extraction_code", "Subspecies")], by = join_by("Taxa_name" == "Current_name")) %>% 
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Specimen_code_Sources, Extraction_code_Sources, Specimen_phylogeny_Extraction_code, note, BioProject_ID, BioSample_ID, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, lat_lon, Collected_by, Identified_by, Description)


### 10.3/ Curate fields ####

## Fill Subspecies for outgroups
SD1_Voucher_specimens_metadata$Subspecies[is.na(SD1_Voucher_specimens_metadata$Subspecies)] <- F
table(SD1_Voucher_specimens_metadata$Subspecies)

## Fill Status_taxa for valid names
SD1_Voucher_specimens_metadata$Status_taxa[is.na(SD1_Voucher_specimens_metadata$Status_taxa)] <- "valid"
table(SD1_Voucher_specimens_metadata$Status_taxa)

## Check Specimen codes
table(SD1_Voucher_specimens_metadata$Voucher_specimen_code == SD1_Voucher_specimens_metadata$Specimen_code_Sources)
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(-Specimen_code_Sources)

## Check Extraction codes
# Fill missing codes for outgroups
SD1_Voucher_specimens_metadata$Specimen_phylogeny_Extraction_code[is.na(SD1_Voucher_specimens_metadata$Specimen_phylogeny_Extraction_code)] <- SD1_Voucher_specimens_metadata$Extraction_code_Sources[is.na(SD1_Voucher_specimens_metadata$Specimen_phylogeny_Extraction_code)]
table(SD1_Voucher_specimens_metadata$Extraction_code_Sources == SD1_Voucher_specimens_metadata$Specimen_phylogeny_Extraction_code)
View(SD1_Voucher_specimens_metadata[!(SD1_Voucher_specimens_metadata$Extraction_code_Sources == SD1_Voucher_specimens_metadata$Specimen_phylogeny_Extraction_code), c("Extraction_code_Sources", "Specimen_phylogeny_Extraction_code", "note")])
# Mismatch for extract with failing attempt ("a" added to the label) # Use the codes from the phylogeny
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  rename(Extraction_code = Specimen_phylogeny_Extraction_code) %>%
  dplyr::select(-Extraction_code_Sources, -note)

## Convert latitude/longitude back to decimals

# Convert long/lat to degrees
degrees_minutes_to_decimal <- function (lat_long)
{
  lat_long_vector <- unlist(str_split(string = lat_long, pattern = " |\\."))
  
  if (lat_long_vector[3] == "S") { sign_lat <- "-" } else { sign_lat <- "" } # Get sign
  lat_degrees <- lat_long_vector[1] # Get degrees
  lat_minutes <- round(as.numeric(lat_long_vector[2]) * 100 / 60, 0) # Get rounded minutes
  
  if (lat_long_vector[6] == "W") { sign_long <- "-" } else { sign_long <- "" } # Get sign
  long_degrees <- lat_long_vector[4] # Get degrees
  long_minutes <- round(as.numeric(lat_long_vector[5]) * 100 / 60, 0) # Get rounded minutes
  
  # Aggregate result
  lat_dec <- as.numeric(paste0(sign_lat, lat_degrees, ".", lat_minutes))
  long_dec <- as.numeric(paste0(sign_long, long_degrees, ".", long_minutes))
  return(list(lat_dec = lat_dec, long_dec = long_dec))
}


SD1_Voucher_specimens_metadata$Latitude_decimal <- SD1_Voucher_specimens_metadata$Longitude_decimal <- NA
for (i in 1:nrow(SD1_Voucher_specimens_metadata))
{
  # i <- 1
  
  lat_long_i <- degrees_minutes_to_decimal(SD1_Voucher_specimens_metadata$lat_lon[i])
  
  SD1_Voucher_specimens_metadata$Latitude_decimal[i] <- lat_long_i$lat_dec
  SD1_Voucher_specimens_metadata$Longitude_decimal[i] <- lat_long_i$long_dec
}
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Extraction_code, BioProject_ID, BioSample_ID, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, Latitude_decimal, Longitude_decimal, Collected_by, Identified_by, Description)


## Add additional fields for sequence references
SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID <- NA
SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID <- NA
SD1_Voucher_specimens_metadata$COX1_GenBank_ID <- NA
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Extraction_code, BioProject_ID, BioSample_ID, Read_sequences_SRA_SRR_ID, UCE_TLS_GenBank_ID, COX1_GenBank_ID, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, Latitude_decimal, Longitude_decimal, Collected_by, Identified_by, Description)

# Save updated table
saveRDS(object = SD1_Voucher_specimens_metadata, file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")
write.xlsx(x = SD1_Voucher_specimens_metadata, file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.xlsx")


### 10.4/ Add GenBank metadata from previous studies ####

# SD1_Voucher_specimens_metadata <- read.xlsx(xlsxFile = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.xlsx")
SD1_Voucher_specimens_metadata <- readRDS(file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

# Load Metadata from Longino & Branstetter 2020 (DOI: 10.1093/isd/ixaa004)
ixaa004_data <- read.xlsx(xlsxFile = "./input_data/Molecular_data/ixaa004_suppl_supplementary_tables.xlsx", sheet = "S4-NCBI-Accession#s")

# Match on BioSample ID
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  left_join(y = ixaa004_data[, c("BioSample#", "SRA.SRR#", "TLS.GenBank#.[UCEs]", "GenBank#.[COI]")], by = join_by("BioSample_ID" == "BioSample#")) %>%
  rename(Read_sequences_SRA_SRR_ID_2 = `SRA.SRR#`) %>%
  rename(UCE_TLS_GenBank_ID_2 = `TLS.GenBank#.[UCEs]`) %>%
  rename(COX1_GenBank_ID_2 = `GenBank#.[COI]`) %>%
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Extraction_code, BioProject_ID, BioSample_ID, Read_sequences_SRA_SRR_ID, Read_sequences_SRA_SRR_ID_2, UCE_TLS_GenBank_ID, UCE_TLS_GenBank_ID_2, COX1_GenBank_ID, COX1_GenBank_ID_2, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, Latitude_decimal, Longitude_decimal, Collected_by, Identified_by, Description)

# Check that ID matches
View(SD1_Voucher_specimens_metadata)

# Remove previous data
SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID <- SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID_2
SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID <- SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID_2
SD1_Voucher_specimens_metadata$COX1_GenBank_ID <- SD1_Voucher_specimens_metadata$COX1_GenBank_ID_2
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(-Read_sequences_SRA_SRR_ID_2, -UCE_TLS_GenBank_ID_2, -COX1_GenBank_ID_2)

# Load Metadata from Branstetter & Longino (2022) (DOI: 10.1093/isd/ixab031)
ixab031_data <- read.xlsx(xlsxFile = "./input_data/Molecular_data/ixab031_suppl_supplementary_tables.xlsx", sheet = "S4-NCBI-Accessions")
ixab031_data <- ixab031_data %>% 
  filter(!is.na(ixab031_data$`BioSample#`))

# Match on BioSample ID
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  left_join(y = ixab031_data[, c("BioSample#", "SRA.SRR#", "TLS.GenBank#.[UCEs]", "GenBank#.[COI]")], by = join_by("BioSample_ID" == "BioSample#")) %>%
  rename(Read_sequences_SRA_SRR_ID_2 = `SRA.SRR#`) %>%
  rename(UCE_TLS_GenBank_ID_2 = `TLS.GenBank#.[UCEs]`) %>%
  rename(COX1_GenBank_ID_2 = `GenBank#.[COI]`) %>%
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Extraction_code, BioProject_ID, BioSample_ID, Read_sequences_SRA_SRR_ID, Read_sequences_SRA_SRR_ID_2, UCE_TLS_GenBank_ID, UCE_TLS_GenBank_ID_2, COX1_GenBank_ID, COX1_GenBank_ID_2, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, Latitude_decimal, Longitude_decimal, Collected_by, Identified_by, Description)

# Check that ID matches
View(SD1_Voucher_specimens_metadata)

# Inform new metadata
SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID[!is.na(SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID_2)] <- SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID_2[!is.na(SD1_Voucher_specimens_metadata$Read_sequences_SRA_SRR_ID_2)]
SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID[!is.na(SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID_2)] <- SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID_2[!is.na(SD1_Voucher_specimens_metadata$UCE_TLS_GenBank_ID_2)]
SD1_Voucher_specimens_metadata$COX1_GenBank_ID[!is.na(SD1_Voucher_specimens_metadata$COX1_GenBank_ID_2)] <- SD1_Voucher_specimens_metadata$COX1_GenBank_ID_2[!is.na(SD1_Voucher_specimens_metadata$COX1_GenBank_ID_2)]
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(-Read_sequences_SRA_SRR_ID_2, -UCE_TLS_GenBank_ID_2, -COX1_GenBank_ID_2)

# SD1_Voucher_specimens_metadata <- read.xlsx(xlsxFile = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata_formatted.xlsx")

# Load Metadata from Camacho et al., 2025
SuppInfoS2_data <- read.xlsx(xlsxFile = "./input_data/Molecular_data/SuppInfoS2_Camacho et al.2024_Zooregions in Madagascar_final_14.i.25.xlsx", sheet = "Table S2.1 - Sequencing info")

SuppInfoS2_data <- SuppInfoS2_data %>% 
  filter(SuppInfoS2_data$Voucher.Code %in% SD1_Voucher_specimens_metadata$Voucher_specimen_code) %>% 
  filter(SuppInfoS2_data$Data.Reference == "This study")

View(SD1_Voucher_specimens_metadata[SD1_Voucher_specimens_metadata$Voucher_specimen_code %in% SuppInfoS2_data$Voucher.Code, ])

# Match on BioSample ID
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  left_join(y = SuppInfoS2_data[, c("Voucher.Code", "BioSample#")], by = join_by("Voucher_specimen_code" == "Voucher.Code")) %>%
  rename(BioSample_ID_2 = `BioSample#`) %>%
  dplyr::select(Taxa_name, Status_taxa, Subspecies, Outgroup, Status_data, Voucher_specimen_code, Extraction_code, BioProject_ID, BioSample_ID, BioSample_ID_2, Read_sequences_SRA_SRR_ID, UCE_TLS_GenBank_ID, COX1_GenBank_ID, Reference, Reference_DOI, Host_institution, Caste, Sex, Collection_date, Locality, Latitude_decimal, Longitude_decimal, Collected_by, Identified_by, Description)

# Check that ID matches
View(SD1_Voucher_specimens_metadata)

# Inform new metadata
SD1_Voucher_specimens_metadata$BioSample_ID[!is.na(SD1_Voucher_specimens_metadata$BioSample_ID_2)] <- SD1_Voucher_specimens_metadata$BioSample_ID_2[!is.na(SD1_Voucher_specimens_metadata$BioSample_ID_2)]
SD1_Voucher_specimens_metadata <- SD1_Voucher_specimens_metadata %>% 
  dplyr::select(-BioSample_ID_2)


# Save updated table
saveRDS(object = SD1_Voucher_specimens_metadata, file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")
write.xlsx(x = SD1_Voucher_specimens_metadata, file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.xlsx")
  


# Need to be pasted in the Supplementary Data file as SD1