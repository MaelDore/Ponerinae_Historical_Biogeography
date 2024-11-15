##### Script 19: Generate Voucher table for NCBI #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Generate Voucher table needed for sequence upload in NCBI - BioSample

###

### Inputs

# AntWeb specimen database accessed on 30 May 2024
# Phylogeny with metadata including voucher codes
# Curated occurrence database

###

### Sources

# AntWeb. Version 8.101. California Academy of Science, online at https://www.antweb.org. Accessed on 30 May 2024.
# NCBI - BioSample: Template for submission: https://submit.ncbi.nlm.nih.gov/biosample/template/

###

### Outputs

# Voucher table following template provided by NCBI - BioSample

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
  mutate(specimen_voucher = str_to_upper(Specimen_code)) %>%
  dplyr::select(specimen_voucher, Latitude_dec, Longitude_dec, Country_ISO3_name, adm1, adm2, Locality) %>%
  filter(specimen_voucher %in% NCBI_BioSample_Voucher_table$specimen_voucher)

View(Biogeographic_database_to_merge)

## Get coordinates in proper format

# Convert long/lat to degrees
decimal_to_degrees_minutes <- function (lat, long)
{
  abs_lat <- abs(lat)
  lat_degrees <- floor(abs_lat) # Get degrees
  lat_minutes <- round((abs_lat - lat_degrees) * 60, 0) # Get rounded minutes
  if (abs_lat == lat) { lat_card <- "N" } else { lat_card <- "S" } # Get cardinal
  
  abs_long <- abs(long)
  long_degrees <- floor(abs_long) # Get degrees
  long_minutes <- round((abs_long - long_degrees) * 60, 0) # Get rounded minutes
  if (abs_long == long) { long_card <- "E" } else { long_card <- "W" } # Get cardinal
  
  # Aggregate result
  lat_long_deg <- paste0(lat_degrees, ".", lat_minutes, " ", lat_card, " ", long_degrees, ".", long_minutes, " ", long_card)
  return(lat_long_deg)
}

# Loop per entry
Biogeographic_database_to_merge$lat_lon <- NA
for (i in 1:nrow(Biogeographic_database_to_merge))
{
  Biogeographic_database_to_merge$lat_lon[i] <- decimal_to_degrees_minutes(lat = Biogeographic_database_to_merge$Latitude_dec[i],
                                                                           long = Biogeographic_database_to_merge$Longitude_dec[i])
}

## Get locality info in proper format

# Replace ":" to avoid issues
Biogeographic_database_to_merge$Country_ISO3_name <- str_replace_all(string = Biogeographic_database_to_merge$Country_ISO3_name, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$adm1 <- str_replace_all(string = Biogeographic_database_to_merge$adm1, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$adm2 <- str_replace_all(string = Biogeographic_database_to_merge$adm2, pattern = ":", replacement = " ")
Biogeographic_database_to_merge$Locality <- str_replace_all(string = Biogeographic_database_to_merge$Locality, pattern = ":", replacement = " ")

# Merge in one field
Biogeographic_database_to_merge$geo_loc_name <- paste0(Biogeographic_database_to_merge$Country_ISO3_name, ": ", Biogeographic_database_to_merge$adm1, ", ", Biogeographic_database_to_merge$adm2, ", ", Biogeographic_database_to_merge$Locality)
  
View(Biogeographic_database_to_merge)

# Merge with voucher data
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>%
  left_join(y = Biogeographic_database_to_merge[, c("specimen_voucher", "geo_loc_name", "lat_lon")], by = join_by(specimen_voucher)) %>% 
  dplyr::select(sample_name, bioproject_accession, organism, isolate,	breed,	host,	isolation_source,	collection_date, geo_loc_name, tissue, altitude, biomaterial_provider, collected_by, dev_stage, identified_by, lat_lon,	sex,	specimen_voucher,	`infra specific rank`, `infra specific name`, description) %>%
  arrange(sample_name)

View(NCBI_BioSample_Voucher_table)


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

View(Biogeographic_database_complement_to_merge)

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

View(Biogeographic_database_complement_to_merge)

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
test <- read.xlsx(xlsxFile = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

##### 7/ Add description field ####

NCBI_BioSample_Voucher_table$description <- paste0("The DNA voucher specimen ", NCBI_BioSample_Voucher_table$specimen_voucher," is deposited at ", NCBI_BioSample_Voucher_table$biomaterial_provider," with data available on www.antweb.org. The voucher is the same specimen as the one used for extraction, or is by default from the same nest/series/locality as the sequenced specimen.")

### Save/Export the Voucher table
saveRDS(object = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
write.xlsx(x = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")



## Info
# dev_stage = lifestage # Need to be manually curated (or removed)
# `infra specific rank` = to tag morphospecies
# `infra specific name` = to name subspecies