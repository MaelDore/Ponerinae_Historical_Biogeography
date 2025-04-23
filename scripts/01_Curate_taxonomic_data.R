##### Script 01: Curate taxonomic information #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Create taxa-level database for taxonomic, morphologic, ecologic and biogeographic information
# Aggregate data from AntWeb.org and AntMaps.org
# Extract list of taxa from phylogeny
# Aggregate current trait measurement data

###

### Inputs

# Ponerinae phylogeny: 793-taxa
# Trait measurements dataset
# AntWeb specimen database accessed on 30 May 2024
# AntMap: GABI Database Release 1.0 from 18 January 2020

###

### Sources

# AntWeb. Version 8.101. California Academy of Science, online at https://www.antweb.org. Accessed on 30 May 2024.
# GABI: Data release 1.0 of 18 January 2020.
    # Guénard, B., Weiser, M., Gomez, K., Narula, N., Economo, E.P. (2017) The Global Ant Biodiversity Informatics (GABI) database: synthesizing data on the geographic distributions of ant species. Myrmecological News 24: 83-89.

###

### Outputs

# Taxa-level database for taxonomic, ecologic and biogeographic information
# Curated AntWeb database
# Curate GABI database
# Biogeographic database merging the two latter (still need to be curated for geographic errors. See Script 2)
# Curated trait measurement database

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####
library(MASS)
library(tidyverse)
library(readxl)
library(xlsx)      # Need the Java Development Kit (JDK) installed
library(openxlsx)  # Use Rccp. No need of Java
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(phytools)
library(ape)

# Use utils::install.packages("package_name") to install package if error

### 1.2/ Load databases ####

AntWeb_database_old <- read_excel("./input_data/AntWeb_data/AntWeb_database_Ponerinae_2024_01_23.xlsx")
AntWeb_database <- read_excel("./input_data/AntWeb_data/AntWeb_database_Ponerinae_2024_02_06.xlsx")

setdiff(names(AntWeb_database_old), names(AntWeb_database))

GABI_database <- read.csv(file = "./input_data/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020.csv", header = T, sep = ",", na.strings = "", quote = '"')

## GABI database was initially broken !!! Need to fix that by detecting wrong access numbers

# All errors comes from misplaced quotation " that needed to be removed or added manually

# which.min(str_sub(GABI_database$gabi_acc_number, start = 0, end = 5) == "GABI_")
# GABI_database$gabi_acc_number[which.min(str_sub(GABI_database$gabi_acc_number, start = 0, end = 5) == "GABI_") - 1]

## Save curated GABI database
saveRDS(GABI_database, file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database.rds")

GABI_database <- readRDS(file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database.rds")

### 1.3/ Load phylogeny ####

# Ponerinae_phylogeny_805t <- read.iqtree(file = "./input_data/Phylogenies/ponerinae-805t-90p.phylip.contree")
# Ponerinae_phylogeny_805t@phylo$tip.label
# 
# Ponerinae_phylogeny_793t <- read.iqtree(file = "./input_data/Phylogenies/ponerinae-793t-75p-iqtree-nopart-gtrg.tre")
# Ponerinae_phylogeny_793t@phylo$tip.label

Ponerinae_phylogeny_792t <- read.tree(file = "./input_data/Phylogenies/ponerinae-792t-0p-spruce83-iqtree-prelim.tre")
Ponerinae_phylogeny_792t$tip.label

# Prune outgroups
outgroups <- c("Amblyopone_australis_D0872", "Paraponera_clavata_EX1573", "Proceratium_google_MAMI0434")
Ponerinae_phylogeny_789t <- tidytree::drop.tip(object = Ponerinae_phylogeny_792t, tip = outgroups)

# Save pruned phylogeny with only ingroups
saveRDS(Ponerinae_phylogeny_789t, file = "./input_data/Phylogenies/Ponerinae_phylogeny_789t.rds")


### 1.4/ Load molecular sample database to retrieve Extraction codes and Antweb codes ####

Phylogeny_samples_data <- read_excel("input_data/Phylogenies/ponerinae-dataset-sample-list-current.xlsx") %>%
  filter(`Ponerinae-dataset` == "KEEP") %>%
  select(`New sample name`, `Extraction code`, `Specimen code`, country, bioregion) %>% 
  rename("Extraction_Code" = `Extraction code`) %>%
  rename("Specimen_Code" = `Specimen code`) %>%
  rename("Phylo_label" = `New sample name`) %>%
  rename('Specimen_phylogeny_Country' = country) %>%
  rename('Specimen_phylogeny_Bioregion' = bioregion) %>%
  filter(!Phylo_label %in% outgroups) %>% # Remove outliers
  filter(!Phylo_label == "Hypoponera_molesta_EX2761") # Remove duplicate from 793t tree

# Check matches

setdiff(Phylogeny_samples_data$Phylo_label, Ponerinae_phylogeny_789t$tip.label)
setdiff(Ponerinae_phylogeny_789t$tip.label, Phylogeny_samples_data$Phylo_label)

# Adjust names with "-"
Phylogeny_samples_data$Phylo_label <- str_replace_all(string = Phylogeny_samples_data$Phylo_label, pattern = "-", replacement = "_")

# Extract the 'name' part of the label
Phylogeny_samples_data$Phylo_name <- NA
for (i in 1:nrow(Phylogeny_samples_data))
{
  # i <- 1
  
  # Extract information
  focal_label <- Phylogeny_samples_data$Phylo_label[i]
  focal_Extraction_code <- Phylogeny_samples_data$Extraction_Code[i]
  focal_Specimen_code <- Phylogeny_samples_data$Specimen_Code[i]
  
  # Remove Extraction code
  focal_name <- str_remove(string = focal_label, pattern = paste0("_", focal_Extraction_code))
  
  # Remove Specimen code
  focal_name <- str_remove(string = focal_name, pattern = paste0("_", focal_Specimen_code))
  
  # Return Phylo_name
  Phylogeny_samples_data$Phylo_name[i] <- focal_name
}

# Save metadata of 789 ingroup taxa in the phylogeny
saveRDS(Phylogeny_samples_data, file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")


##### 2/ Curate databases for taxonomy #####

### 2.1/ Obtain list of valid taxa ####

## Obtain valid Ponerinae Genus names from AntWeb
Ponerinae_Genus_AntWeb_df <- read_excel("input_data/AntWeb_data/Ponerinae_Genus_list.xlsx")
Ponerinae_Genus_AntWeb <- Ponerinae_Genus_AntWeb_df$`Taxon Name`[(Ponerinae_Genus_AntWeb_df$Type == "Valid name") & (Ponerinae_Genus_AntWeb_df$Status == "Extant")]

## Obtain valid Ponerinae Genus names from AntCat
Ponerinae_Genus_AntCat_df <- read_excel("input_data/AntWeb_data/AntCat_names_Ponerinae_2023_11_23.xlsx", sheet = "Ponerinae_Genus")
Ponerinae_Genus_AntCat <- Ponerinae_Genus_AntCat_df$genus

setdiff(Ponerinae_Genus_AntWeb, Ponerinae_Genus_AntCat)
setdiff(Ponerinae_Genus_AntCat, Ponerinae_Genus_AntWeb)

## Obtain valid Ponerinae species names from AntWeb
AntWeb_database$AntWeb_Genus_species <- paste0(str_to_title(AntWeb_database$genus), "_", AntWeb_database$species)
Ponerinae_Species_AntWeb <- unique(AntWeb_database$AntWeb_Genus_species[(AntWeb_database$status == "valid") & !AntWeb_database$is_introduced])

## Obtain valid Ponerinae species names from AntCat
Ponerinae_Species_AntCat_df <- read_excel("./input_data/AntWeb_data/AntCat_names_Ponerinae_2023_11_23.xlsx", sheet = "Ponerinae_Species")
Ponerinae_Species_AntCat_df$Genus_species <- paste0(Ponerinae_Species_AntCat_df$genus, "_", Ponerinae_Species_AntCat_df$species)
Ponerinae_Species_AntCat <- Ponerinae_Species_AntCat_df$Genus_species

# Remove synonym species names
Ponerinae_Species_AntCat_df <- Ponerinae_Species_AntCat_df[!(Ponerinae_Species_AntCat_df$Genus_species == "Leptogenys_foraminosa"), ] # Synonym of Leptogenys_volcanica
Ponerinae_Species_AntCat_df <- Ponerinae_Species_AntCat_df[!(Ponerinae_Species_AntCat_df$Genus_species == "Brachyponera_croceicornis"), ] # Synonym of Brachyponera_obscurans

# Remove subspecies
Ponerinae_Species_AntCat_df <- Ponerinae_Species_AntCat_df[paste0(Ponerinae_Species_AntCat_df$genus, " ", Ponerinae_Species_AntCat_df$species) != Ponerinae_Species_AntCat_df$`current valid parent`, ]

# Save the df for Valid Ponerinae species from AntCat
saveRDS(object = Ponerinae_Species_AntCat_df, file = "./input_data/AntWeb_data/Ponerinae_Species_AntCat_df.rds")

setdiff(Ponerinae_Species_AntWeb, Ponerinae_Species_AntCat) # All taxa in AntWeb that are not in the AntCat list are fossils. Need to remove them from AntWeb by filtering to keep only AntCat taxa
setdiff(Ponerinae_Species_AntCat, Ponerinae_Species_AntWeb) # 165 species without records in AntWeb (but it may still be possible to provide a Bioregion so keep them in the list)


### 2.2/ Extract Ponerinae only from GABI database ####

## 2.2.1/ Filter for valid Genera ####

## Extract 'valid' Genus and genus names from GABI
GABI_database_Genus_names_pub <- unique(GABI_database$genus_name_pub)
GABI_database$Genus_valid <- str_split(GABI_database$valid_species_name, pattern = "\\.", simplify = T)[, 1]
GABI_database$Species_valid <- str_split(GABI_database$valid_species_name, pattern = "\\.", simplify = T)[, 2]
GABI_database_Genus_names_valid <- unique(GABI_database$Genus_valid)

# Convert Genus.species to Genus_species to match names from AntCat
GABI_database$valid_species_name <- str_replace(GABI_database$valid_species_name, pattern = "\\.", replacement = "_")

## Flag GABI dataset for matches of original publication or curated names with Ponerinae Genera from AntCat
# Flag against AntCat names + everything starting with Poner*

GABI_database$Genus_pub_Ponerinae <- GABI_database$genus_name_pub %in% Ponerinae_Genus_AntCat
GABI_database$Genus_pub_Ponerinae_suspected <- !GABI_database$Genus_pub_Ponerinae & (str_sub(string = GABI_database$genus_name_pub, start = 0, end = 5) == "Poner")

table(GABI_database$Genus_pub_Ponerinae)
table(GABI_database$Genus_pub_Ponerinae_suspected)

GABI_database$Genus_GABI_Ponerinae <- GABI_database$Genus_valid %in% Ponerinae_Genus_AntCat
GABI_database$Genus_GABI_Ponerinae_suspected <- !GABI_database$Genus_GABI_Ponerinae & (str_sub(string = GABI_database$Genus_valid, start = 0, end = 5) == "Poner")

table(GABI_database$Genus_GABI_Ponerinae)
table(GABI_database$Genus_GABI_Ponerinae_suspected) # No curated Genus names starting with Poner* is not already matching a valid Ponerinae name

## Validate entry based on matching flags for Genus

# Entry with Genus names from curated GABI database matching AntCat list are considered valid
# GABI_database$Ponerinae_entry <- GABI_database$Genus_GABI_Ponerinae

# Entry with Genus names from original publication matching AntCat list, but not listed as Ponerinae by GABI need to be inspected.
# See if it reflects a change of taxonomy or if there is a mistake
View(GABI_database[GABI_database$Genus_pub_Ponerinae & !GABI_database$Genus_GABI_Ponerinae, ])
# Seems invalid Ponerinae if we have to trust GABI's curation

# Entry with Genus names from original publication starting with Poner* need to be inspected
# Check if misspelled Genus or true Genus outside of Ponerinae
View(GABI_database[GABI_database$Genus_pub_Ponerinae_suspected, ])

GABI_database[GABI_database$Genus_pub_Ponerinae_suspected, ]$genus_name_pub
# Wrong space included in the Genus pub name that prevents matching. But still detected as valid from the curated name so not an issue.

## 2.2.2/ Filter for valid Species ####

## Flag GABI dataset for matches of original publication or curated names with Ponerinae species from AntCat

# Original publication
GABI_database$pub_Genus_species_name <- paste0(GABI_database$genus_name_pub, "_", GABI_database$species_name_pub)
GABI_database$Species_pub_Ponerinae <- GABI_database$pub_Genus_species_name %in% Ponerinae_Species_AntCat

table(GABI_database$Species_pub_Ponerinae)

# Curated name
GABI_database$Species_GABI_Ponerinae <- GABI_database$valid_species_name %in% Ponerinae_Species_AntCat

table(GABI_database$Genus_GABI_Ponerinae)
table(GABI_database$Species_GABI_Ponerinae)

## Detect and correct curated names that do not match AntCat anymore, but were matching in original pub

# Entry with Species names from original publication matching AntCat list, but not listed as Ponerinae by GABI need to be inspected.
# See if it reflects a change of taxonomy or if there is a mistake
View(GABI_database[GABI_database$Species_pub_Ponerinae & !GABI_database$Species_GABI_Ponerinae, ])
# Neoponera_procidua needs to be changed back to Pachycondyla_procidua
# Cryptopone_guianensis needs to be changed back to Wadeura_guianensis
# Leptogenys_indigatrix needs to be changed back to Leptogenys_indagatrix
# Leptogenys_longiscapa needs to be changed back to Leptogenys_longiscapus
# Heteroponera_mayri needs to be changed back to Anochetus_mayri
# Iridomyrmex_mayri needs to be changed back to Anochetus_mayri
# Myopias_maligna punctigera needs to be changed for Myopias_maligna; punctigera is the subspecies
# Other are not Ponerinae, or the former name does not match the location, so better remove the entry anyway.

## 2.2.3/ Correct mistakes to match AntCat names ####

# Create a new field to keep tracks of initial GABI's curated names
GABI_database$AntCat_Genus_species_name <- GABI_database$valid_species_name
GABI_database$AntCat_Genus_species_name[!(GABI_database$AntCat_Genus_species_name %in% Ponerinae_Species_AntCat)] <- NA

# Neoponera_procidua needs to be changed back to Pachycondyla_procidua
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Neoponera_procidua") & (GABI_database$pub_Genus_species_name == "Pachycondyla_procidua")] <- "Pachycondyla_procidua"
# Cryptopone_guianensis needs to be changed back to Wadeura_guianensis
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Cryptopone_guianensis") & (GABI_database$pub_Genus_species_name == "Wadeura_guianensis")] <- "Wadeura_guianensis"
# Leptogenys_indigatrix needs to be changed back to Leptogenys_indagatrix
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Leptogenys_indigatrix") & (GABI_database$pub_Genus_species_name == "Leptogenys_indagatrix")] <- "Leptogenys_indagatrix"
# Leptogenys_longiscapa needs to be changed back to Leptogenys_longiscapus
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Leptogenys_longiscapa") & (GABI_database$pub_Genus_species_name == "Leptogenys_longiscapus")] <- "Leptogenys_longiscapus"
# Heteroponera_mayri needs to be changed back to Anochetus_mayri
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Heteroponera_mayri") & (GABI_database$pub_Genus_species_name == "Anochetus_mayri")] <- "Anochetus_mayri"
# Iridomyrmex_mayri needs to be changed back to Anochetus_mayri
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Iridomyrmex_mayri") & (GABI_database$pub_Genus_species_name == "Anochetus_mayri")] <- "Anochetus_mayri"
# Myopias_maligna punctigera needs to be changed for Myopias_maligna; punctigera is the subspecies
GABI_database$AntCat_Genus_species_name[(GABI_database$valid_species_name == "Myopias_maligna punctigera") & (GABI_database$pub_Genus_species_name == "Myopias_papua")] <- "Myopias_maligna"

## 2.2.4/ Validate entry based on matching flags for Species ####

# Rerun flagging on AntCat names to validate corrected entries
GABI_database$Ponerinae_entry <- GABI_database$AntCat_Genus_species_name %in% Ponerinae_Species_AntCat

## Extract valid Ponerinae entries
GABI_database_Ponerinae <- GABI_database[GABI_database$Ponerinae_entry, ]

## Remove dubious entries
table(GABI_database_Ponerinae$dubious)

GABI_database_Ponerinae <- GABI_database_Ponerinae[is.na(GABI_database_Ponerinae$dubious), ]

## Save curated GABI Ponerinae database
saveRDS(GABI_database_Ponerinae, file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database_Ponerinae.rds")

### 2.3/ Curate AntWeb database for valid names and dubious data ####

AntWeb_database_curated <- AntWeb_database

## 2.3.1/ Create new fields to record changes and keep track of AntWeb names ####

AntWeb_database_curated$Current_name <- AntWeb_database_curated$AntWeb_Genus_species
AntWeb_database_curated$Current_genus <- AntWeb_database_curated$genus
AntWeb_database_curated$Current_species <- AntWeb_database_curated$species
AntWeb_database_curated$Current_status <- AntWeb_database_curated$status
AntWeb_database_curated$Current_subspecies <- F


## 2.3.2/ Update names of identified indet to prevent removing them ####

# CASENT0007651 (EX2728): Hagensia_indet  =>  Hagensia_havilandi_marleyi
AntWeb_database_curated$Current_name[(AntWeb_database_curated$SpecimenCode == "casent0007651")] <- "Hagensia_havilandi_marleyi"
AntWeb_database_curated$Current_species[(AntWeb_database_curated$SpecimenCode == "casent0007651")] <- "havilandi_marleyi"
AntWeb_database_curated$Current_status[(AntWeb_database_curated$SpecimenCode == "casent0007651")] <- "valid"
AntWeb_database_curated$Current_subspecies[(AntWeb_database_curated$SpecimenCode == "casent0007651")] <- T

# CASENT0649902 (EX2442): Neoponera_(indet) => Neoponera_bucki
AntWeb_database_curated$Current_name[(AntWeb_database_curated$SpecimenCode == "casent0649902")] <- "Neoponera_bucki"
AntWeb_database_curated$Current_species[(AntWeb_database_curated$SpecimenCode == "casent0649902")] <- "bucki"
AntWeb_database_curated$Current_status[(AntWeb_database_curated$SpecimenCode == "casent0649902")] <- "valid"

# # CASENT0295217 (EX2727): Fisheropone_(indet)  => Fisheropone_ambigua
# AntWeb_database_curated$Current_name[(AntWeb_database_curated$SpecimenCode == "casent0295217")] <- "Fisheropone_ambigua"
# AntWeb_database_curated$Current_species[(AntWeb_database_curated$SpecimenCode == "casent0295217")] <- "ambigua"
# AntWeb_database_curated$Current_status[(AntWeb_database_curated$SpecimenCode == "casent0295217")] <- "valid"

# #	CASENT0776822 (EX3642): Parvaponera_(indet) =>  Parvaponera_suspecta
# AntWeb_database_curated$Current_name[(AntWeb_database_curated$SpecimenCode == "casent0776822")] <- "Parvaponera_suspecta"
# AntWeb_database_curated$Current_species[(AntWeb_database_curated$SpecimenCode == "casent0776822")] <- "suspecta"
# AntWeb_database_curated$Current_status[(AntWeb_database_curated$SpecimenCode == "casent0776822")] <- "valid"


## 2.3.3/ Filter out fossils, indet and introduced entries ####

## Keep only valid names and morphotaxa
table(AntWeb_database_curated$Current_status)

AntWeb_database_curated <- AntWeb_database_curated[AntWeb_database_curated$Current_status %in% c("valid", "morphotaxon", "new"), ]

## Filter out fossils
# All valid taxa in AntWeb that are not in the AntCat list are fossils. Need to remove them from AntWeb by filtering to keep only AntCat taxa
fossils_names <- setdiff(Ponerinae_Species_AntWeb, Ponerinae_Species_AntCat)
fossils_names <- fossils_names[!is.na(fossils_names)]; fossils_names 

AntWeb_database_curated <- AntWeb_database_curated[!(AntWeb_database_curated$Current_name %in% fossils_names), ]

## Remove entries for introduced species
table(AntWeb_database_curated$is_introduced)
AntWeb_database_curated <- AntWeb_database_curated[AntWeb_database_curated$is_introduced == F, ]

## Adjust status of valid names from AntWeb to match valid names from AntCat
Ponerinae_valid_Species_AntWeb <- unique(AntWeb_database_curated$Current_name[AntWeb_database_curated$Current_status %in% c("valid")])
Ponerinae_morphotaxa_Species_AntWeb <- unique(AntWeb_database_curated$Current_name[AntWeb_database_curated$Current_status %in% c("morphotaxon")])

Ponerinae_valid_Species_AntCat_in_AntWeb <- Ponerinae_valid_Species_AntWeb[Ponerinae_valid_Species_AntWeb %in% Ponerinae_Species_AntCat]
Ponerinae_not_valid_Species_AntCat_valid_in_AntWeb <- Ponerinae_valid_Species_AntWeb[!(Ponerinae_valid_Species_AntWeb %in% Ponerinae_Species_AntCat)]
# All valid species in Antweb are valid in AntCat. Only differences are subspecies

Ponerinae_morphotaxa_AntWeb_valid_in_AntCat <- Ponerinae_morphotaxa_Species_AntWeb[Ponerinae_morphotaxa_Species_AntWeb %in% Ponerinae_Species_AntCat]
# One morphotaxon that is valid in AntCat: "Bothroponera_pumicosa sculpturata_nr". It is a subspecies of a valid species. Not a mistake.
Ponerinae_morphotaxa_AntWeb_not_valid_AntCat <- Ponerinae_morphotaxa_Species_AntWeb[!(Ponerinae_morphotaxa_Species_AntWeb %in% Ponerinae_Species_AntCat)]


## 2.3.4/ Curate names in Antweb based on recent taxonomic changes ####

# Add field to record case where changes in AntWeb is not recommended
AntWeb_database_curated$AntWeb_changes_required <- NA

## Cases with wrong name in AntWeb. Must change all specimens

# Ponera_sinensis_nr => Ponera guangxiensis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Ponera_sinensis_nr"] <- "Parvaponera_suspecta"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Ponera_sinensis_nr"] <- "guangxiensis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Ponera_sinensis_nr"] <- "valid"

# Mesoponera_afrc-gh01  => Mesoponera_afr05
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-gh01"] <- "Mesoponera_afr05"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-gh01"] <- "afr05"

# # Neoponera_bucki => NewGenus_bucki
# AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Neoponera_bucki"] <- "NewGenus_bucki"
# AntWeb_database_curated$Current_genus[AntWeb_database_curated$AntWeb_Genus_species == "Neoponera_bucki"] <- "NewGenus"
# AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Neoponera_bucki"] <- "new"

# Brachyponera_croceicornis => Brachyponera_obscurans
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Brachyponera_croceicornis"] <- "Brachyponera_obscurans"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Brachyponera_croceicornis"] <- "obscurans"

# Hypoponera_casc-mz26 => Hypoponera_dulcis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz26"] <- "Hypoponera_dulcis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz26"] <- "dulcis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz26"] <- "valid"

# Hypoponera_ao03 => Hypoponera_dulcis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao03"] <- "Hypoponera_dulcis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao03"] <- "dulcis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao03"] <- "valid"

# Hypoponera_ug01 => Hypoponera_dulcis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug01"] <- "Hypoponera_dulcis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug01"] <- "dulcis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug01"] <- "valid"

# Hypoponera_afrc-za04 => Hypoponera_spei
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_afrc-za04"] <- "Hypoponera_spei"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_afrc-za04"] <- "spei"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_afrc-za04"] <- "Hypoponera_spei"

# Hypoponera_casc-mz31 => Hypoponera_spei
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz31"] <- "Hypoponera_spei"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz31"] <- "spei"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz31"] <- "valid"

# Hypoponera_ug03 => Hypoponera_tristis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug03"] <- "Hypoponera_tristis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug03"] <- "tristis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug03"] <- "valid"

# Hypoponera_ug04 => Hypoponera_tristis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug04"] <- "Hypoponera_tristis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug04"] <- "tristis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug04"] <- "valid"

# Hypoponera_ug05 => Hypoponera_tristis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug05"] <- "Hypoponera_tristis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug05"] <- "tristis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug05"] <- "valid"

# Hypoponera_ug06 => Hypoponera_tristis
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug06"] <- "Hypoponera_tristis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug06"] <- "tristis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug06"] <- "valid"

#	Leptogenys_casc-mz04 => Leptogenys_castanea
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_casc-mz04"] <- "Leptogenys_castanea"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_casc-mz04"] <- "castanea"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_casc-mz04"] <- "valid"

#	Mesoponera_afrc-ug01 => Mesoponera_ambigua
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-ug01"] <- "Mesoponera_ambigua"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-ug01"] <- "ambigua"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-ug01"] <- "valid"

# Mesoponera_afrc-zm01 => Mesoponera_ambigua
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-zm01"] <- "Mesoponera_ambigua"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-zm01"] <- "ambigua"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_afrc-zm01"] <- "valid"

#	Odontoponera_ph01 => Odontoponera_denticulata
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_ph01"] <- "Odontoponera_denticulata"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_ph01"] <- "denticulata"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_ph01"] <- "valid"

# Odontoponera_th01 => Odontoponera_denticulata
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_th01"] <- "Odontoponera_denticulata"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_th01"] <- "denticulata"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_th01"] <- "valid"

# Odontoponera_np01 => Odontoponera_transversa
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_np01"] <- "Odontoponera_transversa"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_np01"] <- "transversa"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Odontoponera_np01"] <- "valid"

#	Hypoponera_casc-mz06 => Hypoponera_angustata
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz06"] <- "Hypoponera_angustata"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz06"] <- "angustata"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_casc-mz06"] <- "valid"

# Hypoponera_ug12 => Hypoponera_blanda
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug12"] <- "Hypoponera_blanda"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug12"] <- "blanda"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug12"] <- "valid"

# Hypoponera_ug07 => Hypoponera_fatiga
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug07"] <- "Hypoponera_fatiga"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug07"] <- "fatiga"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug07"] <- "valid"

# Hypoponera_ug08 => Hypoponera_fatiga
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug08"] <- "Hypoponera_fatiga"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug08"] <- "fatiga"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ug08"] <- "valid"

# Hypoponera_ao02 => Hypoponera_jeanneli
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao02"] <- "Hypoponera_jeanneli"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao02"] <- "jeanneli"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_ao02"] <- "valid"

# Thaumatomyrmex_ferox_complex => Thaumatomyrmex_ferox  
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Thaumatomyrmex_ferox_complex"] <- "Thaumatomyrmex_ferox"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Thaumatomyrmex_ferox_complex"] <- "ferox"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Thaumatomyrmex_ferox_complex"] <- "valid"
AntWeb_database_curated$AntWeb_changes_required[AntWeb_database_curated$AntWeb_Genus_species == "Thaumatomyrmex_ferox_complex"] <- "No"

# Leptogenys_pubiceps_complex => Leptogenys_pubiceps
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_pubiceps_complex"] <- "Leptogenys_pubiceps"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_pubiceps_complex"] <- "pubiceps"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_pubiceps_complex"] <- "valid"
AntWeb_database_curated$AntWeb_changes_required[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_pubiceps_complex"] <- "No"

#	Odontomachus_tuneri_cf  =>  Odontomachus_turneri  
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Odontomachus_tuneri_cf"] <- "Odontomachus_turneri"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Odontomachus_tuneri_cf"] <- "turneri"
AntWeb_database_curated$Current_status[AntWeb_database_curated$AntWeb_Genus_species == "Odontomachus_tuneri_cf"] <- "valid"

# Hypoponera_sc-mora => Hypoponera_sc_mora
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-mora"] <- "Hypoponera_sc_mora"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-mora"] <- "sc_mora"

# Hypoponera_sc-nosy => Hypoponera_sc_nosy
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-nosy"] <- "Hypoponera_sc_nosy"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-nosy"] <- "sc_nosy"

# Hypoponera_sc-rano => Hypoponera_sc_rano
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-rano"] <- "Hypoponera_sc_rano"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-rano"] <- "sc_rano"

# Hypoponera_sc-tamp => Hypoponera_sc_tamp
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-tamp"] <- "Hypoponera_sc_tamp"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-tamp"] <- "sc_tamp"

#	Hypoponera_sc-zaha => Hypoponera_sc_zaha
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-zaha"] <- "Hypoponera_sc_zaha"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_sc-zaha"] <- "sc_zaha"

# Leptogenys_foraminosa => Leptogenys_volcanica  
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_foraminosa"] <- "Leptogenys_volcanica"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_foraminosa"] <- "volcanica"

# Diacamma_timor_01 => Diacamma_mj_timor01 
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_timor_01"] <- "Diacamma_mj_timor01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_timor_01"] <- "mj_timor01"

# Diacamma_janda_sp1 => Diacamma_mj_sp1 
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_janda_sp1"] <- "Diacamma_mj_sp1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_janda_sp1"] <- "mj_sp1"

# Ectomomyrmex_janda_sp7 => Ectomomyrmex_mj _sp7
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_janda_sp1"] <- "Diacamma_mj_sp1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Diacamma_janda_sp1"] <- "mj_sp1"

#	Hypoponera_fogo07 => Hypoponera_mj_fogo07 (Check for case non-sensitive matches)
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_fogo07"] <- "Hypoponera_mj_fogo07"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_fogo07"] <- "mj_fogo07"

# Hypoponera_weam01 => Hypoponera_mj_weam01
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_weam01"] <- "Hypoponera_mj_weam01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_weam01"] <- "mj_weam01"

# Mesoponera_janda_sp6 => Mesoponera_mj_sp6 
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp6"] <- "Mesoponera_mj_sp6"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp6"] <- "mj_sp6"

# Myopias_janda_sp6 => Myopias_mj_sp6
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_janda_sp6"] <- "Myopias_mj_sp6"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_janda_sp6"] <- "mj_sp6"

# Myopias_utai_1 => Myopias_mj_utai_1 
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_utai_1"] <- "Myopias_mj_utai_1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_utai_1"] <- "mj_utai_1"

# Myopias_janda_sp9 => Myopias_mj_sp9 
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_janda_sp9"] <- "Myopias_mj_sp9"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Myopias_janda_sp9"] <- "mj_sp9"

# Leptogenys_janda_sp1 => Leptogenys_mj_sp1
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_janda_sp1"] <- "Leptogenys_mj_sp1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Leptogenys_janda_sp1"] <- "mj_sp1"

# Hypoponera_mad02 => Hypoponera_mj_mad02  
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_mad02"] <- "Hypoponera_mj_mad02"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_mad02"] <- "mj_mad02"

# Ponera_psw-my01 => Ponera_psw_my01
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Ponera_psw-my01"] <- "Ponera_psw_my01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Ponera_psw-my01"] <- "psw_my01"

# Mesoponera_janda_sp2 => Mesoponera_mj_sp2
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp2"] <- "Mesoponera_mj_sp2"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp2"] <- "mj_sp2"

# Mesoponera_janda_sp3 => Mesoponera_mj_sp3
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp3"] <- "Mesoponera_mj_sp3"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Mesoponera_janda_sp3"] <- "mj_sp3"

# Hypoponera_jtl033 => Hypoponera_dias19_2
AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_jtl033"] <- "Hypoponera_dias19_2"
AntWeb_database_curated$Current_species[AntWeb_database_curated$AntWeb_Genus_species == "Hypoponera_jtl033"] <- "mj_dias19_2"

## Cases with systematic changes of '-' to '_' in morphotaxa (not in valid species names)
table((AntWeb_database_curated$Current_status == "morphotaxon") & str_detect(string = AntWeb_database_curated$Current_name, pattern = "-"))
AntWeb_database_curated$AntWeb_changes_required[(AntWeb_database_curated$Current_status == "morphotaxon") & str_detect(string = AntWeb_database_curated$Current_name, pattern = "-")] <- "dash"
AntWeb_database_curated$Current_name[AntWeb_database_curated$Current_status == "morphotaxon"] <- str_replace_all(string = AntWeb_database_curated$Current_name[AntWeb_database_curated$Current_status == "morphotaxon"], pattern = "-", replacement = "_")
AntWeb_database_curated$Current_species[AntWeb_database_curated$Current_status == "morphotaxon"] <- str_replace_all(string = AntWeb_database_curated$Current_species[AntWeb_database_curated$Current_status == "morphotaxon"], pattern = "-", replacement = "_")
table(AntWeb_database_curated$AntWeb_changes_required)

## Cases with systematic changes of '_janda_' to '_mj_'
table(str_detect(string = AntWeb_database_curated$Current_name, pattern = "janda"))
# No more taxa with 'janda'


## Cases with wrong name in AntWeb. Change only the voucher specimen (not sure about synonymy of the other ones)

#	CASENT0820826 (D2423): Anochetus_punctaticeps_cf => Anochetus_katonae_nr1
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0820826"] <- "Anochetus_katonae_nr1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0820826"] <- "katonae_nr1"

#	CASENT0066876 (EX2927):  Anochetus_siphneus  =>  Anochetus_katonae_nr2 
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0066876"] <- "Anochetus_katonae_nr2"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0066876"] <- "katonae_nr2"

# CASENT0886661: Hypoponera_opacior => Hypoponera_opacior_nr1
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0886661"] <- "Hypoponera_opacior_nr1"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0886661"] <- "opacior_nr1"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0886661"] <- "morphotaxon"

# CASENT0886861: Hypoponera_opacior => Hypoponera_opacior_nr2
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0886861"] <- "Hypoponera_opacior_nr2"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0886861"] <- "opacior_nr2"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0886861"] <- "morphotaxon"

# CASENT0886662: Hypoponera_opacior => Hypoponera_opacior_nr3
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0886662"] <- "Hypoponera_opacior_nr3"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0886662"] <- "opacior_nr3"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0886662"] <- "morphotaxon"

# CASENT0649078 (EX2340): Odontomachus_brunneus  =>  Odontomachus_brunneus_nr 
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "Odontomachus_brunneus_nr"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "brunneus_nr"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "morphotaxon"

# CASENT0649078 (EX2340): Odontomachus_brunneus  =>  Odontomachus_brunneus_nr 
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "Odontomachus_brunneus_nr"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "brunneus_nr"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0649078"] <- "morphotaxon"

# CASENT0817421 (D2422): Anochetus_afrc-tz01 => Anochetus_talpa  
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0817421"] <- "Anochetus_talpa"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0817421"] <- "talpa"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0817421"] <- "valid"
AntWeb_database_curated$AntWeb_changes_required[AntWeb_database_curated$SpecimenCode == "casent0817421"] <- "No"  # According to Peter Hawkes's will

#	CASENT0650371: Brachyponera_weam01 => Brachyponera_lutea
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650371"] <- "Brachyponera_lutea"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650371"] <- "lutea"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650371"] <- "valid"

# CASENT0650282 (EX3044) => Hypoponera_aliena
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650282"] <- "Hypoponera_aliena"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650282"] <- "aliena"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650282"] <- "valid"

# CASENT0923397 (D2588) => Hypoponera_aliena
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "Hypoponera_aliena"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "aliena"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "valid"

# CASENT0923397 (EX2810) => Hypoponera_aliena
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "Hypoponera_aliena"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "aliena"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0923397"] <- "valid"

#	CASENT0803835: Mesoponera_ambigua  => Mesoponera_elisae
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0803835"] <- "Mesoponera_elisae"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0803835"] <- "elisae"

# CASENT0816544: Mesoponera_ingesta  => Mesoponera_caffraria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0816544"] <- "Mesoponera_caffraria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0816544"] <- "caffraria"

# CASENT0374097: Cryptopone_th01 (AntWeb) / Ectomomyrmex_th01 (Phylogeny) => Ectomomyrmex_th05
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0374097"] <- "Ectomomyrmex_th05"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0374097"] <- "th05"

# CASENT0631798 (EX2282):  Neoponera_crenata  =>  Neoponera_moesta
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0631798"] <- "Neoponera_moesta"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0631798"] <- "moesta"

# CASENT0649900 (EX2451): Neoponera_moesta  =>  Neoponera_crenata  
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0649900"] <- "Neoponera_crenata"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0649900"] <- "crenata"

#	ZRC_HYM_0000557 (EX2685): Brachyponera_jerdonii  =>  Brachyponera_nigrita
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_hym_0000557"] <- "Brachyponera_nigrita"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "zrc_hym_0000557"] <- "nigrita"

# ZRC_ENT00007921 (EX2689): Ectomomyrmex_overbecki_cf  =>  Ectomomyrmex_overbecki  
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_ent00007921"] <- "Ectomomyrmex_overbecki"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "zrc_ent00007921"] <- "overbecki"

# CASENT0650114 (EX2777): Leptogenys_peninsularis (AntWeb) / Leptogenys_sonora (phylogeny) =>  Leptogenys_peninsularis_nr
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650114"] <- "Leptogenys_peninsularis_nr"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650114"] <- "peninsularis_nr"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650114"] <- "morphotaxon"

# CASENT0650188 (EX3010): Odontomachus_opaculus_c  => Odontomachus_imperator
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650188"] <- "Odontomachus_imperator"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650188"] <- "imperator"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650188"] <- "valid"

# CASENT0650198 (EX3020): Ectomomyrmex_scobinus  =>  Ectomomyrmex_aciculatus
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650198"] <- "Ectomomyrmex_aciculatus"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650198"] <- "aciculatus"

#	ZRC_ENT00047817 (EX2683): Leptogenys_kraepelini => Leptogenys_harmsi
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_ent00047817"] <- "Leptogenys_harmsi"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "zrc_ent00047817"] <- "harmsi"

# CASENT0284148 (EX2853): Leptogenys_th01 => Leptogenys_harmsi
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0284148"] <- "Leptogenys_harmsi"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0284148"] <- "harmsi"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0284148"] <- "valid"

# CASENT0650200 (EX3022): Ectomomyrmex_janda_sp10  =>  Ectomomyrmex_simillimus
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650200"] <- "Ectomomyrmex_simillimus"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650200"] <- "simillimus"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650200"] <- "valid"

#	CASENT0387827 (EX2916): Anochetus_my07  =>  Anochetus_rugosus
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0387827"] <- "Anochetus_rugosus"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0387827"] <- "rugosus"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0387827"] <- "valid"

# CASENT0354210 (EX2939): Anochetus_afr05 (Anochetus_ug01 in Phylogeny) => Anochetus_afr06 (Better to use afr06 than ug01 because CASENT0923412 (EX2937) in the same clade is in Gabon!)
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0354210"] <- "Anochetus_afr06"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0354210"] <- "afr06"

# CASENT0821692 (D2446): Anochetus_afrc_tz05  =>  Anochetus_afr05 (Different clade from above)
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0821692"] <- "Anochetus_afr05"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0821692"] <- "afr05"

#	CASENT0774577 (EX3640): Parvaponera_casc_mz02  => Fisheropone_ambigua
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0774577"] <- "Fisheropone_ambigua"
AntWeb_database_curated$Current_genus[AntWeb_database_curated$SpecimenCode == "casent0774577"] <- "Fisheropone"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0774577"] <- "ambigua"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0774577"] <- "valid"

# CASENT0821449 (D2436): Bothroponera_afrc-za03  =>  Bothroponera_strigulosa
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0821449"] <- "Bothroponera_strigulosa"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0821449"] <- "strigulosa"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0821449"] <- "valid"

# CASENT0356474 (EX2786): Hypoponera_ug02 =>  Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0356474"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0356474"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0356474"] <- "valid"

# CASENT0351350 (EX2797): Hypoponera_ug15 =>  Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0351350"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0351350"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0351350"] <- "valid"

# CASENT0888109 (D2459): Hypoponera_afrc_za06 =>  Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0888109"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0888109"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0888109"] <- "valid"

# CASENT0774326 (EX2802): Hypoponera_casc-mz11 =>  Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0774326"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0774326"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0774326"] <- "valid"

# CASENT0160001 (MAMI1267):  Hypoponera_mg116 =>  Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0160001"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0160001"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0160001"] <- "valid"

# CASENT0261090 (MAMI0742): Hypoponera_punctatissima => Hypoponera_ergatandria
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0261090"] <- "Hypoponera_ergatandria"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0261090"] <- "ergatandria"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0261090"] <- "valid"

#	CASENT0361603 (EX2798): Hypoponera_ug16  =>  Hypoponera_molesta
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0361603"] <- "Hypoponera_molesta"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0361603"] <- "molesta"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0361603"] <- "valid"

# CASENT0391901 (EX2748): Cryptopone_my09  => Cryptopone_typhlos
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0391901"] <- "Cryptopone_typhlos"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0391901"] <- "typhlos"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0391901"] <- "valid"

# CASENT0386062 (EX2743): Cryptopone_my04 => Cryptopone_typhlos
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0386062"] <- "Cryptopone_typhlos"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0386062"] <- "typhlos"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0386062"] <- "valid"

# CASENT0923391 (EX2738): Cryptopone_id01 => Cryptopone_typhlos
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0923391"] <- "Cryptopone_typhlos"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0923391"] <- "typhlos"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0923391"] <- "valid"

# CASENT0923392 (EX2739): Cryptopone_my03 => Cryptopone_butteli
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0923392"] <- "Cryptopone_butteli"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0923392"] <- "butteli"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0923392"] <- "valid"

# CASENT0390060 (EX2745): Cryptopone_my06 => Cryptopone_butteli
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0390060"] <- "Cryptopone_butteli"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0390060"] <- "butteli"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0390060"] <- "valid"

# CASENT0278864 (EX2963): Ectomomyrmex_th03 => Ectomomyrmex_punctatus
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0278864"] <- "Ectomomyrmex_punctatus"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0278864"] <- "punctatus"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0278864"] <- "valid"

#	CASENT0114965 (EX2954): Odontomachus_tuneri_cf  =>  Odontomachus_turneri  
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0114965"] <- "Odontomachus_turneri"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0114965"] <- "turneri"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0114965"] <- "valid"

# CASENT0386602 (EX2996): Myopias_my07  =>  Myopias_my01
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0386602"] <- "Myopias_my01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0386602"] <- "my01"

#	CASENT0250487 (D2479): Leptogenys_afrc_tz10  =>  Leptogenys_comajojo
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0250487"] <- "Leptogenys_comajojo"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0250487"] <- "comajojo"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0250487"] <- "valid"

#	CASENT0777597 (EX2890):  Leptogenys_casc_mz05  =>  Leptogenys_comajojo
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0777597"] <- "Leptogenys_comajojo"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0777597"] <- "comajojo"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0777597"] <- "valid"

#	CASENT0131922 (EX2978): Emeryopone_TH02  =>  Emeryopone_buttelreepeni
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0131922"] <- "Emeryopone_buttelreepeni"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0131922"] <- "buttelreepeni"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0131922"] <- "valid"

#	CASENT0389113 (EX2981): Harpegnathos_my02  => Harpegnathos_my01
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0389113"] <- "Harpegnathos_my01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0389113"] <- "my01"

#	CASENT0379431 (EX2913):  Anochetus_ag01  =>  Anochetus_neglectus
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0379431"] <- "Anochetus_neglectus"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0379431"] <- "neglectus"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0379431"] <- "valid"

# CASENT0887835 (EX2581): Pseudoneoponera_au04 => Pseudoneoponera_sublaevis
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0887835"] <- "Pseudoneoponera_sublaevis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0887835"] <- "sublaevis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0887835"] <- "valid"

# CASENT0887834 (EX2582): Pseudoneoponera_au07 => Pseudoneoponera_sublaevis
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0887834"] <- "Pseudoneoponera_sublaevis"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0887834"] <- "sublaevis"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0887834"] <- "valid"

#	CASENT0722234 (EX2993): Myopias_my04  =>  Myopias_maligna
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0722234"] <- "Myopias_maligna"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0722234"] <- "maligna"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0722234"] <- "valid"

#	CASENT0634948 (EX2368): Myopias_bidens_cf  =>  Myopias_concava_nr
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0634948"] <- "Myopias_concava_nr"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0634948"] <- "concava_nr"

#	CASENT0650203 (EX3025): Myopias_bgc33  =>  Myopias_bg02
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650203"] <- "Myopias_bg02"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650203"] <- "bg02"

# ZRC_ENT00007263 (EX2687) Myopias_mayri =>  Myopias_breviloba
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_ent00007263"] <- "Myopias_breviloba"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "zrc_ent00007263"] <- "breviloba"

# CASENT0923423 (EX2999):  Myopias_my10 =>  Myopias_breviloba
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0923423"] <- "Myopias_breviloba"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0923423"] <- "breviloba"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0923423"] <- "valid"

# CASENT0783788 (EX2892): Leptogenys_casc_mz07  => Leptogenys_afrc_tz02
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0783788"] <- "Leptogenys_afrc_tz02"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0783788"] <- "afrc_tz02"

# CASENT0775384 (EX2887): Leptogenys_casc_mz01  =>  Leptogenys_afrc_tz07
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0775384"] <- "Leptogenys_afrc_tz07"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0775384"] <- "afrc_tz07"

#	CASENT0906223 (EX3637): Parvaponera_afr04  =>  Parvaponera_afr01
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0906223"] <- "Parvaponera_afr01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0906223"] <- "afr01"

# CASENT0356584 (EX3641): Parvaponera_ug02  =>  Parvaponera_afr01
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0356584"] <- "Parvaponera_afr01"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0356584"] <- "afr01"

#	JTLC000000507 (EX3636): Bothroponera_silvestrii  =>  Bothroponera_silvestrii_nr
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "jtlc000000507"] <- "Bothroponera_silvestrii_nr"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "jtlc000000507"] <- "silvestrii_nr"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "jtlc000000507"] <- "morphotaxon"

# CASENT0650197 (EX3019): Myopias_janda_sp4 => Myopias_papua
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650197"] <- "Myopias_papua"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0650197"] <- "papua"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0650197"] <- "valid"

# CASENT0888244 (D2452): Hagensia_havilandi => Hagensia_havilandi_marleyi  (Only in our database so we record it differently than Hagensia_havilandi)
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0888244"] <- "Hagensia_havilandi_marleyi"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0888244"] <- "havilandi_marleyi"
AntWeb_database_curated$Current_subspecies[AntWeb_database_curated$SpecimenCode == "casent0888244"] <- TRUE

# CASENT0787946 (EX2729): Hagensia_peringueyi => Hagensia_peringueyi_saldanhae (Only in our database so we record it differently than Hagensia_peringueyi)
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0787946"] <- "Hagensia_peringueyi_saldanhae"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0787946"] <- "peringueyi_saldanhae"
AntWeb_database_curated$Current_subspecies[AntWeb_database_curated$SpecimenCode == "casent0787946"] <- TRUE

# CASENT0888216 (D2429): Bothroponera_afrc_mz02 => Bothroponera_silvestrii
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0888216"] <- "Bothroponera_silvestrii"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0888216"] <- "silvestrii"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0888216"] <- "valid"

# CASENT0872952 (EX2988): Neoponera_fiebrigi_cf => Neoponera_fiebrigi 
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0872952"] <- "Neoponera_fiebrigi"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "casent0872952"] <- "fiebrigi"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "casent0872952"] <- "valid"

# UFV-LABECOL-000500 (EX2438): Neoponera_schultzi1 => Neoponera_schultzi 
AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "ufv-labecol-000500"] <- "Neoponera_schultzi"
AntWeb_database_curated$Current_species[AntWeb_database_curated$SpecimenCode == "ufv-labecol-000500"] <- "schultzi"
AntWeb_database_curated$Current_status[AntWeb_database_curated$SpecimenCode == "ufv-labecol-000500"] <- "valid"


## 2.3.5/ Curate names in the phylogeny based on curated AntWeb data ####

## Check presence of specimens in AntWeb database
Phylogeny_sample_data_789t <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")
Introduced_vouchers_Codes <- Phylogeny_sample_data_789t$Specimen_Code[!(str_to_lower(Phylogeny_sample_data_789t$Specimen_Code) %in% AntWeb_database_curated$SpecimenCode)]
Introduced_vouchers_Codes
Phylogeny_sample_data_789t$Specimen_Code[!(str_to_lower(Phylogeny_sample_data_789t$Specimen_Code) %in% AntWeb_database$SpecimenCode)]

# "CASENT0649763" "CASENT0158872" "CASENT0261090" "CASENT0818943" "CASENT0071744" "CASENT0136413" "CASENT0159701" "CASENT0376164"

# Missing specimens from the curated AntWeb database are in introduced locations
# Retrieve taxonomic information using the non-curated database only for those specimens
# CASENT0261090 need to be renamed: Hypoponera_punctatissima => Hypoponera_ergatandria

# Extract voucher specimen absent from the curated database because they are in introduced locations
Introduced_vouchers_database <- AntWeb_database %>% 
  mutate(SpecimenCode = str_to_upper(SpecimenCode)) %>% 
  filter(SpecimenCode %in% Introduced_vouchers_Codes) %>% 
  mutate(Current_name_introduced = AntWeb_Genus_species)
  
# CASENT0261090 need to be renamed: Hypoponera_punctatissima => Hypoponera_ergatandria
Introduced_vouchers_database$Current_name_introduced[Introduced_vouchers_database$SpecimenCode == "CASENT0261090"] <- "Hypoponera_ergatandria"

# Merge data with the taxa-level df for samples in the phylogeny
Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  left_join(y = Introduced_vouchers_database[, c("SpecimenCode", "Current_name_introduced")], by = c("Specimen_Code" = "SpecimenCode"))

# Add in the taxa-level df the information on occurrence location of vouchers: introduced or not
Phylogeny_sample_data_789t$Specimen_phylogeny_Is_introduced <- FALSE
Phylogeny_sample_data_789t$Specimen_phylogeny_Is_introduced[(Phylogeny_sample_data_789t$Specimen_Code %in% Introduced_vouchers_Codes)] <- TRUE

## Match names using the curated AntWeb database

Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  mutate(Specimen_Code_lower = str_to_lower(Specimen_Code)) %>% 
  left_join(y = AntWeb_database_curated[, c("SpecimenCode", "Current_name")], by = c("Specimen_Code_lower" = "SpecimenCode")) %>% 
  select(-Specimen_Code_lower)

# Fill Current_name for introduced vouchers
Phylogeny_sample_data_789t$Current_name[Phylogeny_sample_data_789t$Specimen_phylogeny_Is_introduced] <- Phylogeny_sample_data_789t$Current_name_introduced[Phylogeny_sample_data_789t$Specimen_phylogeny_Is_introduced]
Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  select(-Current_name_introduced)

## Compare names and tag mismatches

Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  mutate(Specimen_phylogeny_Mismatched_name = Phylo_name != Current_name) %>% 
  mutate(Specimen_phylogeny_Mismatched_name_CAPS = Specimen_phylogeny_Mismatched_name & (str_to_title(Phylo_name) == Current_name)) %>% 
  mutate(Specimen_phylogeny_Mismatched_name_Content = Specimen_phylogeny_Mismatched_name & !Specimen_phylogeny_Mismatched_name_CAPS)

table(Phylogeny_sample_data_789t$Specimen_phylogeny_Mismatched_name)
table(Phylogeny_sample_data_789t$Specimen_phylogeny_Mismatched_name_CAPS)
table(Phylogeny_sample_data_789t$Specimen_phylogeny_Mismatched_name_Content)

View(Phylogeny_sample_data_789t[Phylogeny_sample_data_789t$Specimen_phylogeny_Mismatched_name, ])

# AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "casent0650201"]
# AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_ent00007921"]
# AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "zrc_hym_0000557"]
# AntWeb_database_curated$Current_name[AntWeb_database_curated$SpecimenCode == "atpfor2006"]


# Look for duplicates that highlight possible taxonomic update
Phylogeny_sample_data_789t$Current_name[duplicated(Phylogeny_sample_data_789t$Current_name)]
# No duplicates


# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_789t, file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")


### 2.3.7/ Flag for presence or not in the phylogeny ####

# In the taxa-level summary table
Phylogeny_sample_data_789t$In_phylogeny <- T

# In the curated AntWeb database
AntWeb_database_curated$In_phylogeny <- F
AntWeb_database_curated$In_phylogeny[(AntWeb_database_curated$Current_name %in% Phylogeny_sample_data_789t$Current_name)] <- T

table(AntWeb_database_curated$In_phylogeny)

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")


### 2.3.8/ Flag for specimens that need names updates in AntWeb ####

AntWeb_database_curated$AntWeb_Mismatch_name <- AntWeb_database_curated$AntWeb_Genus_species != AntWeb_database_curated$Current_name
AntWeb_database_curated$AntWeb_changes_required[is.na(AntWeb_database_curated$AntWeb_changes_required) & AntWeb_database_curated$AntWeb_Mismatch_name] <- "Yes"

table(AntWeb_database_curated$AntWeb_Mismatch_name)
table(AntWeb_database_curated$AntWeb_changes_required)

table(AntWeb_database_curated$Current_name[AntWeb_database_curated$AntWeb_changes_required == "No"])

View(AntWeb_database_curated[AntWeb_database_curated$AntWeb_Mismatch_name, ])

## Extract list of specimens requiring changes in AntWeb

AntWeb_database_Need_updates <- AntWeb_database_curated %>% 
  filter(AntWeb_Mismatch_name) %>% 
  filter(AntWeb_changes_required %in% c("dash", "Yes")) %>%
  # filter(AntWeb_Genus_species != "Neoponera_bucki") %>%
  select(SpecimenCode, AntWeb_Genus_species, Current_name, In_phylogeny, AntWeb_changes_required, ownedby, access_group) %>% 
  arrange(AntWeb_changes_required)

table(AntWeb_database_Need_updates$AntWeb_changes_required)

# Export list of specimens requiring changes in AntWeb
saveRDS(object = AntWeb_database_Need_updates, file = "./input_data/AntWeb_data/AntWeb_database_Need_updates.rds")
write.xlsx(x = AntWeb_database_Need_updates, file = "./input_data/AntWeb_data/AntWeb_database_Need_updates.xlsx")

# Aggregate cases to count number of changes per taxa
AntWeb_database_Need_updates_aggregated <- AntWeb_database_Need_updates %>% 
  group_by(AntWeb_Genus_species, Current_name) %>% 
  select(-SpecimenCode) %>%
  summarise_all(first) %>%
  ungroup() %>%
  arrange(AntWeb_changes_required)

table(AntWeb_database_Need_updates_aggregated$AntWeb_changes_required)

# Export list of specimens requiring changes in AntWeb, aggregated to taxa level
saveRDS(object = AntWeb_database_Need_updates_aggregated, file = "./input_data/AntWeb_data/AntWeb_database_Need_updates_aggregated.rds")
write.xlsx(x = AntWeb_database_Need_updates_aggregated, file = "./input_data/AntWeb_data/AntWeb_database_Need_updates_aggregated.xlsx")

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")


## 2.3.9/ Inform Taxa-level summary df for Changes in AntWeb names and Current status ####

## Inform from changes in AntWeb and current name status

# # Remove columns if already present
# Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>%
#   select(-Current_status, Current_subspecies, AntWeb_Genus_species, AntWeb_Mismatch_name, AntWeb_changes_required)

# Join with AntWeb data
Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
    mutate(Specimen_Code_lower = str_to_lower(Specimen_Code)) %>% 
    left_join(y = AntWeb_database_curated[, c("SpecimenCode", "Current_status", "Current_subspecies", "AntWeb_Genus_species", "AntWeb_Mismatch_name", "AntWeb_changes_required")], by = c("Specimen_Code_lower" = "SpecimenCode"))

## Complete data for specimen that are introduced, thus not in curated database
Introduced_vouchers_database_curated <- Introduced_vouchers_database %>% 
  rename(Current_status_introduced = status) %>% 
  mutate(Current_subspecies_introduced = FALSE) %>%
  rename(AntWeb_Genus_species_introduced = AntWeb_Genus_species) %>%
  mutate(AntWeb_Mismatch_name_introduced = AntWeb_Genus_species_introduced != Current_name_introduced)
  
Introduced_vouchers_database_curated$AntWeb_changes_required_introduced <- NA
Introduced_vouchers_database_curated$AntWeb_changes_required_introduced[Introduced_vouchers_database_curated$AntWeb_Mismatch_name_introduced] <- "Yes"

Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>%
  left_join(y = Introduced_vouchers_database_curated[, c("SpecimenCode", "Current_status_introduced", "Current_subspecies_introduced", "AntWeb_Genus_species_introduced", "AntWeb_Mismatch_name_introduced", "AntWeb_changes_required_introduced")], by = c("Specimen_Code" = "SpecimenCode"))

Phylogeny_sample_data_789t$Current_status[is.na(Phylogeny_sample_data_789t$Current_status)] <- Phylogeny_sample_data_789t$Current_status_introduced[is.na(Phylogeny_sample_data_789t$Current_status)] 
Phylogeny_sample_data_789t$Current_subspecies[is.na(Phylogeny_sample_data_789t$Current_subspecies)] <- Phylogeny_sample_data_789t$Current_subspecies_introduced[is.na(Phylogeny_sample_data_789t$Current_subspecies)] 
Phylogeny_sample_data_789t$AntWeb_Genus_species[is.na(Phylogeny_sample_data_789t$AntWeb_Genus_species)] <- Phylogeny_sample_data_789t$AntWeb_Genus_species_introduced[is.na(Phylogeny_sample_data_789t$AntWeb_Genus_species)] 
Phylogeny_sample_data_789t$AntWeb_Mismatch_name[is.na(Phylogeny_sample_data_789t$AntWeb_Mismatch_name)] <- Phylogeny_sample_data_789t$AntWeb_Mismatch_name_introduced[is.na(Phylogeny_sample_data_789t$AntWeb_Mismatch_name)] 
Phylogeny_sample_data_789t$AntWeb_changes_required[is.na(Phylogeny_sample_data_789t$AntWeb_changes_required)] <- Phylogeny_sample_data_789t$AntWeb_changes_required_introduced[is.na(Phylogeny_sample_data_789t$AntWeb_changes_required)]

Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  select(-Current_status_introduced, -Current_subspecies_introduced, -AntWeb_Genus_species_introduced, -AntWeb_Mismatch_name_introduced, -AntWeb_changes_required_introduced)

## Rename and select adequate fields
Phylogeny_sample_data_789t <- Phylogeny_sample_data_789t %>% 
  rename(Subspecies = Current_subspecies) %>% 
  rename(Specimen_phylogeny_AntWeb_name = AntWeb_Genus_species) %>% 
  rename(Specimen_phylogeny_AntWeb_name_Mismatch = AntWeb_Mismatch_name) %>% 
  rename(Specimen_phylogeny_AntWeb_name_Changes_required = AntWeb_changes_required) %>% 
  select(Specimen_Code, Extraction_Code, Phylo_name, Phylo_label, Current_name, Current_status, Subspecies, In_phylogeny,
         Specimen_phylogeny_Country, Specimen_phylogeny_Bioregion, Specimen_phylogeny_Is_introduced,
         Specimen_phylogeny_Mismatched_name, Specimen_phylogeny_Mismatched_name_CAPS, Specimen_phylogeny_Mismatched_name_Content,
         Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_AntWeb_name_Mismatch, Specimen_phylogeny_AntWeb_name_Changes_required)
  
# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_789t, file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

## 2.3.10/ Add previous comments from ponerinae-dataset-sample-list-current ####

## Merge previous comments from ponerinae-dataset-sample-list-current

Former_summary_df <- read_excel("./input_data/Phylogenies/ponerinae-dataset-sample-list-current.xlsx")
Former_comments <- Former_summary_df %>% 
  rename(Notes = `Notes - 2022-04-26`) %>% 
  rename(Specimen_Code = `Specimen code`) %>%
  filter(Former_summary_df$`Ponerinae-dataset` == "KEEP") %>% 
  select(Specimen_Code, Notes)
  
Phylogeny_sample_data_789t <- left_join(x = Phylogeny_sample_data_789t, y = Former_comments, by = c("Specimen_Code" = "Specimen_Code"))

# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_789t, file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

## 2.3.11/ Extract changes in Phylo_labels to update phylogenetic data and Gabi-current file ####

Phylogeny_sample_data_Need_updates <- Phylogeny_sample_data_789t %>% 
  mutate(Phylo_label_Updated = paste0(Current_name, "_", Extraction_Code, "_", Specimen_Code)) %>%
  mutate(Phylo_label_To_update = Phylo_label != Phylo_label_Updated) %>%
  mutate(Phylo_label_Codes_mismatch = Phylo_label_To_update & !Specimen_phylogeny_Mismatched_name) %>%
  # filter(Phylo_label_To_update) %>% 
  select(Specimen_Code, Phylo_label, Phylo_label_Updated, Phylo_name, Current_name, Phylo_label_To_update, Phylo_label_Codes_mismatch, Specimen_phylogeny_Mismatched_name, Specimen_phylogeny_Mismatched_name_CAPS, Specimen_phylogeny_Mismatched_name_Content, Specimen_phylogeny_AntWeb_name_Changes_required)

# Detect reason for required change
Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required <- NA
for (i in 1:nrow(Phylogeny_sample_data_Need_updates))
{
  if(Phylogeny_sample_data_Need_updates$Phylo_label_Codes_mismatch[i]) { Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required[i] <- "Codes" }
  if(Phylogeny_sample_data_Need_updates$Specimen_phylogeny_Mismatched_name_CAPS[i]) { Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required[i] <- "CAPS" }
  if(is.na(Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required[i])) { Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required[i] <- "Name" }
  if(!Phylogeny_sample_data_Need_updates$Phylo_label_To_update[i]) { Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required[i] <- NA }
}
table(Phylogeny_sample_data_Need_updates$Phylo_label_Changes_required)

# Order taxa by reasons
Phylogeny_sample_data_Need_updates <- Phylogeny_sample_data_Need_updates %>% 
  select(Specimen_Code, Phylo_label, Phylo_label_Updated, Phylo_name, Current_name, Phylo_label_Changes_required) %>%
  # filter(Phylo_label_To_update) %>% 
  arrange(desc(Phylo_label_Changes_required)) 

## Export summary table for changes
saveRDS(object = Phylogeny_sample_data_Need_updates, file = "./input_data/Phylogenies/Phylogeny_sample_data_Need_updates.rds")
write.xlsx(x = Phylogeny_sample_data_Need_updates, file = "./input_data/Phylogenies/Phylogeny_sample_data_Need_updates.xlsx")


## 2.3.12/ Flag discrepancy with ponerinae-dataset-sample-list-current and update names ####


Former_summary_df <- read_excel("./input_data/Phylogenies/ponerinae-dataset-sample-list-current.xlsx")

Former_summary_df$`Specimen code`[!(Former_summary_df$`Specimen code` %in% Phylogeny_sample_data_Need_updates$Specimen_Code)]
Phylogeny_sample_data_Need_updates$Specimen_Code[!(Phylogeny_sample_data_Need_updates$Specimen_Code %in% Former_summary_df$`Specimen code`)]

Phylogeny_sample_data_Need_updates$Phylo_label[!(Phylogeny_sample_data_Need_updates$Phylo_label %in% Former_summary_df$`New sample name`)]

# Discrepancy are similar to Phylogeny label + 5 cases of 'sc-' taxa that needs to be changed for 'sc_'

# Extract useful information to add to Gabi's current file
Updates_for_former_summary_df <- AntWeb_database_curated %>% 
  mutate(Specimen_Code = str_to_upper(SpecimenCode)) %>%
  select(Specimen_Code, Current_name, AntWeb_Genus_species)

# Add new labels/names to Gabi's current file
Former_summary_df <- dplyr::left_join(x = Former_summary_df, y = Updates_for_former_summary_df, by = c("Specimen code" = "Specimen_Code")) %>% 
   rename(Specimen_Code = `Specimen code`) %>%
   rename(Extraction_Code = `Extraction code`) %>%
   mutate(Phylo_label_updated = paste0(Current_name, "_", Extraction_Code, "_", Specimen_Code))

# Export results to apply to Gabi's current file on DropBox
write.xlsx(x = Former_summary_df, file = "./input_data/Phylogenies/Former_summary_df.xlsx")

# Also add manually the 3 outgroups, and the specimens that are not in the curated AntWeb database...


## 2.3.13/ Export curated AntWeb database with different levels of inclusion ####

### Remove morphospecies outside of the phylogeny

# Morphospecies may be useful to compute sampling fraction across bioregions, even if we may use only the ones in the phylogeny.
# Valid species without molecular data may be useful to compute sampling fraction across bioregions.

# Save curated AntWeb database with all morphospecies
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")
nrow(AntWeb_database_curated)
# 62,374 specimen data

AntWeb_database_curated_without_morphospecies_no_molecular <- AntWeb_database_curated[(AntWeb_database_curated$In_phylogeny | AntWeb_database_curated$Current_status %in% c("valid", "new")), ]
nrow(AntWeb_database_curated_without_morphospecies_no_molecular)
# 51,678 specimen data

# Save curated AntWeb database without morphospecies without molecular data
saveRDS(AntWeb_database_curated_without_morphospecies_no_molecular, file = "./input_data/AntWeb_data/AntWeb_database_curated_without_morphospecies_no_molecular.rds")

### Keep only taxa in the phylogeny

AntWeb_database_curated_only_phylo <- AntWeb_database_curated[(AntWeb_database_curated$In_phylogeny), ]
nrow(AntWeb_database_curated_only_phylo)
# 45,278 specimen data

## Save curated AntWeb database with only taxa in the phylogeny
saveRDS(AntWeb_database_curated_only_phylo, file = "./input_data/AntWeb_data/AntWeb_database_curated_only_phylo.rds")


### 2.3.14/ Go back to GABI database to flag names using the current curated names for terminals in the phylogeny ####

### Add fields for Names status (Current_status), presence in phylogeny (In_phylogeny), and Name update (Current_name), 

# All GABI entries have valid AntCat species names
GABI_database_Ponerinae$Current_status <- "valid"

unique(GABI_database_Ponerinae$AntCat_Genus_species_name[!(GABI_database_Ponerinae$AntCat_Genus_species_name %in% Phylogeny_sample_data_789t$Current_name)])

# Flag for match with a taxa in the phylogeny
GABI_database_Ponerinae$In_phylogeny <- F
GABI_database_Ponerinae$In_phylogeny[(GABI_database_Ponerinae$AntCat_Genus_species_name %in% Phylogeny_sample_data_789t$Current_name)] <- T

table(GABI_database_Ponerinae$In_phylogeny)

# Flag for a name update between the curated GABI name ("valid_species_name") and the current one curated by myself ("AntCat_Genus_species_name")

GABI_database_Ponerinae$Name_updated <- !(GABI_database_Ponerinae$valid_species_name == GABI_database_Ponerinae$AntCat_Genus_species_name)

table(GABI_database_Ponerinae$Name_updated)

## Save curated GABI Ponerinae database
saveRDS(GABI_database_Ponerinae, file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database_Ponerinae.rds")


### 2.4/ Merge both database ####

# Load the curated database
GABI_database_Ponerinae <- readRDS(file = "./input_data/GABI_Data_Release1.0_18012020/GABI_database_Ponerinae.rds")
AntWeb_database_curated <- readRDS(file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")

# Add a Source field for GABI or AntWeb
GABI_database_Ponerinae$Source <- "GABI"
AntWeb_database_curated$Source <- "AntWeb"

# Adapt databases for merging: remove useless fields and rename fields to be matched
GABI_database_Ponerinae_for_merging <- GABI_database_Ponerinae %>% 
  rename(GABI_accession_ID = gabi_acc_number) %>% 
  rename(Accession_ID = accession_number) %>% 
  rename(Initial_name = valid_species_name) %>% 
  rename(Current_name = AntCat_Genus_species_name) %>% 
  rename(Locality = locality) %>% 
  rename(Country = country) %>% 
  rename(Latitude_dec = dec_lat) %>%  
  rename(Longitude_dec = dec_long) %>% 
  select(GABI_accession_ID, Accession_ID, Current_name, Initial_name, Name_updated, Current_status, In_phylogeny, Latitude_dec, Longitude_dec, Country, adm1, adm2, Locality, bentity2_name, Source)
  
AntWeb_database_for_merging <- AntWeb_database_curated %>% 
  rename(Specimen_code = SpecimenCode) %>% 
  rename(Country = country) %>% 
  rename(Locality = localityname) %>%
  rename(Locality_code = localitycode) %>%
  rename(Locality_notes = localitynotes) %>% 
  rename(Latitude_dec = decimal_latitude) %>%  
  rename(Longitude_dec = decimal_longitude) %>% 
  rename(Elevation = elevation) %>% 
  rename(Initial_name = AntWeb_Genus_species) %>% 
  rename(Name_updated = AntWeb_Mismatch_name) %>% 
  select(Specimen_code, Current_name, Initial_name, Name_updated, Current_status, In_phylogeny, Latitude_dec, Longitude_dec, Country, adm1, adm2, Locality, Locality_code, Locality_notes, Elevation, Source)

# Merge databases
Biogeographic_database_Ponerinae <- full_join(GABI_database_Ponerinae_for_merging, AntWeb_database_for_merging)

nrow(GABI_database_Ponerinae_for_merging) + nrow(AntWeb_database_for_merging)
dim(Biogeographic_database_Ponerinae)
# 155,942 geolocalized specimen records before biogeographic curation

# Flag duplicates based on taxa name and geographic information 
Biogeographic_database_Ponerinae$duplicate <- duplicated(Biogeographic_database_Ponerinae[, c("Current_name", "Latitude_dec", "Longitude_dec", "adm1", "adm2")])

table(Biogeographic_database_Ponerinae$duplicate)
# 63,862 unique localized records before biogeographic curation

## Save merged specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")


##### 3/ Build taxa-level summary df #####

## Extent Phylogeny_sample_data_789t to include all valid species.
## Build several df with different levels of inclusion: all, phylo + valid.

### 3.1/ Load and prepare Phylogeny summary table and AntCat Species catalog for merging ####

## Load df for Metadata of samples in the phylogeny
Phylogeny_sample_data_789t <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")
Phylogeny_sample_data_789t_for_merging <- Phylogeny_sample_data_789t %>% 
  # rename(Specimen_phylogeny_AntWeb_name = AntWeb_Genus_species) %>%
  rename(Specimen_phylogeny_Specimen_code = "Specimen_Code") %>% 
  rename(Specimen_phylogeny_Extraction_code = "Extraction_Code") %>% 
  rename(Specimen_phylogeny_Phylo_name_Mismatch = Specimen_phylogeny_Mismatched_name)
  
# Load AntCat list for valid Ponerinae species
Ponerinae_Species_AntCat_df <- readRDS("./input_data/AntWeb_data/Ponerinae_Species_AntCat_df.rds")

# Remove subspecies
Ponerinae_Species_AntCat_df <- Ponerinae_Species_AntCat_df[paste0(Ponerinae_Species_AntCat_df$genus, " ", Ponerinae_Species_AntCat_df$species) != Ponerinae_Species_AntCat_df$`current valid parent`, ]

Ponerinae_Species_AntCat_for_merging <- Ponerinae_Species_AntCat_df %>% 
  rename(Current_status = status) %>% 
  rename(Type_Country = country) %>% 
  rename(Type_Bioregion = bioregion) %>% 
  rename(Current_name = Genus_species) %>% 
  select(Current_name, Current_status, Type_Country, Type_Bioregion)

### 3.2/ Merge summary df ####

## Merge both summary df to add other valid species in the df
Ponerinae_Macroevolution_taxa_database <- full_join(Phylogeny_sample_data_789t_for_merging, Ponerinae_Species_AntCat_for_merging) %>% 
  mutate(Available_occurrences = NA) %>%  # Create field to record availability of occurrence data in GABI + AntWeb
  mutate(To_include_in_analyses = NA) %>%  # Create field to use in analyses
  mutate(Conservative_clade_Node_ID = NA) %>%  # Create fields to record description, node ID and identity of terminal to use to recover the MRCA for a conservative clade we are sure the terminal belongs to
  mutate(Conservative_clade_Terminals_with_MRCA = NA) %>%
  mutate(Conservative_clade_Notes = NA) %>%
  mutate(Conservative_clade_Source = NA) %>%
  # mutate(Exclusive_subclades_Node_ID = NA) %>%  # Create fields to record description, node ID and identity of terminal to use to recover the MRCA for the least inclusive clade we think the terminal may belong to
  # mutate(Exclusive_subclades_Terminals_with_MRCA = NA) %>%
  # mutate(Exclusive_subclades_Notes = NA) %>%
  select(Current_name, Current_status, In_phylogeny, Specimen_phylogeny_Specimen_code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_AntWeb_name_Mismatch, Specimen_phylogeny_Extraction_code, Phylo_label, Phylo_name, Specimen_phylogeny_Phylo_name_Mismatch, Subspecies, Specimen_phylogeny_Is_introduced,
         Conservative_clade_Node_ID, Conservative_clade_Terminals_with_MRCA, Conservative_clade_Notes, Conservative_clade_Source,
         # Exclusive_subclades_Node_ID, Exclusive_subclades_Terminals_with_MRCA, Exclusive_subclades_Notes,
         Specimen_phylogeny_Country, Type_Country, Specimen_phylogeny_Bioregion, Type_Bioregion, Available_occurrences, To_include_in_analyses, Notes)

table(duplicated(Ponerinae_Macroevolution_taxa_database$Current_name))

View(Ponerinae_Macroevolution_taxa_database[duplicated(Ponerinae_Macroevolution_taxa_database$Current_name), ])

### 3.3/ Fill new fields ####

## Fill current fields for valid species outside of the phylogeny
Ponerinae_Macroevolution_taxa_database$In_phylogeny[is.na(Ponerinae_Macroevolution_taxa_database$In_phylogeny)] <- F
Ponerinae_Macroevolution_taxa_database$Subspecies[is.na(Ponerinae_Macroevolution_taxa_database$Subspecies)] <- F

## Format as Logical
Ponerinae_Macroevolution_taxa_database$In_phylogeny <- as.logical(Ponerinae_Macroevolution_taxa_database$In_phylogeny)
Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_AntWeb_name_Mismatch <- as.logical(Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_AntWeb_name_Mismatch)
Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Phylo_name_Mismatch <- as.logical(Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Phylo_name_Mismatch)
# Ponerinae_Macroevolution_taxa_database$Current_duplicate <- as.logical(Ponerinae_Macroevolution_taxa_database$Current_duplicate)
Ponerinae_Macroevolution_taxa_database$Subspecies <- as.logical(Ponerinae_Macroevolution_taxa_database$Subspecies)

table(Ponerinae_Macroevolution_taxa_database$Subspecies)
table(Ponerinae_Macroevolution_taxa_database$In_phylogeny)
# table(Ponerinae_Macroevolution_taxa_database$Current_duplicate)

## Fill fields for clade inclusion for terminals in the phylogeny
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Terminals_with_MRCA[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- Ponerinae_Macroevolution_taxa_database$Phylo_label[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] 
# Ponerinae_Macroevolution_taxa_database$Exclusive_subclades_Terminals_with_MRCA[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- Ponerinae_Macroevolution_taxa_database$Phylo_label[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] 
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Node_ID <- tidytree::nodeid(tree = Ponerinae_phylogeny_789t, label = Ponerinae_Macroevolution_taxa_database$Phylo_label)
# Ponerinae_Macroevolution_taxa_database$Exclusive_subclades_Node_ID <- tidytree::nodeid(tree = Ponerinae_phylogeny_789t, label = Ponerinae_Macroevolution_taxa_database$Phylo_label)
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- NA
# Ponerinae_Macroevolution_taxa_database$Exclusive_subclades_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- NA
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Source <- NA
# Ponerinae_Macroevolution_taxa_database$Conservative_clade_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- "Terminal in the current phylogeny. Conservative clade is the terminal itself. No need for grafting."
# Ponerinae_Macroevolution_taxa_database$Exclusive_subclades_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- "Terminal in the current phylogeny. Least inclusive clade is the termial itself. No need for grafting."

## Fill field for presence/absence of occurrences

AntWeb_database_taxa_names <- unique(AntWeb_database_curated$Current_name)
GABI_database_taxa_names <- unique(GABI_database_Ponerinae$AntCat_Genus_species_name)

Ponerinae_Macroevolution_taxa_database$Available_occurrences <- (Ponerinae_Macroevolution_taxa_database$Current_name %in% c(AntWeb_database_taxa_names, GABI_database_taxa_names))

table(Ponerinae_Macroevolution_taxa_database$Available_occurrences)
# Only 11 taxa without occurrences (but bioregions can still be retrieved from type descriptions)

## Add field for inclusion or not in the analyses

# Only flag as FALSE the taxa with no geographic information in the database (GABI + AntWeb) and from AntCat (type)
No_geographic_info <- is.na(Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Bioregion) & is.na(Ponerinae_Macroevolution_taxa_database$Type_Bioregion) & !(Ponerinae_Macroevolution_taxa_database$Available_occurrences)
table(No_geographic_info) # All taxa have geographic information!

# Flag as TRUE taxa in the phylogeny
Ponerinae_Macroevolution_taxa_database$To_include_in_analyses[Ponerinae_Macroevolution_taxa_database$In_phylogeny] <- T

## Save taxa-level summary df
saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")


### 3.4/ Merge with trait dataset ####

# Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")
# Ponerinae_Macroevolution_taxa_database <- read.xlsx(xlsxFile = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")

# Load the curated AntWeb database
AntWeb_database_curated <- readRDS(file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")

# Load trait dataset

# Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2024_01_12.xlsx")
# Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2024_02_05.xlsx")
# Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2024_11_14.xlsx")
# Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2025_01_08.xlsx")
# Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2025_01_13.xlsx")
Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2025_03_10.xlsx")

Trait_database <- Trait_database %>% 
  rename(Specimen_measured_Code = MeasSpecimen) %>% 
  rename(Extraction_Code = EXcode) %>%
  rename(Taxa_measurement_status = EX_status) %>% 
  mutate(Is_Ponerinae = TRUE) %>% # Add new field for Ponerinae membership
  # dplyr::select(-c("PmSP", "PrSP", "PtSP")) %>%   # Remove spinescence measurements that have been abandoned
  dplyr::select(Date, Lab, Source, Extraction_Code, Specimen_measured_Code, Current_name, Taxa_measurement_status, Is_Ponerinae, HW, HL, SL, ED, WL, PW, MtFL, PmSP, PrSP, PtSP, Notes) %>%
  filter((Extraction_Code != "EX3093") | is.na(Extraction_Code)) %>%   # Specimen is a male
  filter(Specimen_measured_Code != "INBIOCRI0012278883"| is.na(Specimen_measured_Code)) # Remove duplicate with wrong Specimen code

# Inspect duplicates based on Specimen_measured_Code
table(duplicated(Trait_database$Specimen_measured_Code, incomparables = NA))
View(Trait_database[duplicated(Trait_database$Specimen_measured_Code), ])
View(Trait_database[duplicated(Trait_database$Specimen_measured_Code, incomparables = NA), ])

# Remove duplicates based on Extraction code

table(duplicated(Trait_database$Extraction_Code, incomparables = NA))
View(Trait_database[duplicated(Trait_database$Extraction_Code), ])
View(Trait_database[duplicated(Trait_database$Extraction_Code, incomparables = NA), ])

Trait_database <- Trait_database[!duplicated(Trait_database$Extraction_Code, incomparables = NA), ]
# Trait_database <- Trait_database[!duplicated(Trait_database$Extraction_Code), ] # Remove NA too

# # Remove entries with no measurements data
# no_measurements <- apply(X = Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")], MARGIN = 1, FUN = function (x) { all(is.na(x))} )
# Trait_database <- Trait_database[!no_measurements, ]

### 3.4.1/ Retrieve name of specimens from the phylogeny ####

# Load metadata from phylogeny
Phylogeny_sample_data_789t <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

# Extract useful fields
Phylogeny_sample_data_Specimen_code <- Phylogeny_sample_data_789t[, c("Specimen_Code", "Current_name")] %>% 
  rename(Current_name_from_Specimen_code = Current_name)
Phylogeny_sample_data_Extraction_code <- Phylogeny_sample_data_789t[, c("Extraction_Code", "Current_name")] %>% 
  rename(Current_name_from_Extraction_code = Current_name)

# Match on specimen code
Trait_database <- left_join(Trait_database, Phylogeny_sample_data_Extraction_code, by = c("Extraction_Code" = "Extraction_Code"))
# Match on extraction code
Trait_database <- left_join(Trait_database, Phylogeny_sample_data_Specimen_code, by = c("Specimen_measured_Code" = "Specimen_Code"))

# Merge both matches
Trait_database <- Trait_database %>% 
  rename(Current_name_Check = Current_name_from_Extraction_code)
Trait_database$Current_name_Check[is.na(Trait_database$Current_name_Check)] <- Trait_database$Current_name_from_Specimen_code[is.na(Trait_database$Current_name_Check)]
Trait_database <- Trait_database %>% 
  dplyr::select(-Current_name_from_Specimen_code)

### 3.4.2/ Retrieve missing name of specimens from AntWeb databases ####

# # Remove columns if already present
# Trait_database <- Trait_database %>%
#   select(-Current_name, -Specimen_measured_AntWeb_name)
#   # select(-Current_name.x, -Specimen_measured_AntWeb_name.x, -Current_name.y, -Specimen_measured_AntWeb_name.y)

AntWeb_database_extract <- AntWeb_database_curated[, c("SpecimenCode", "Current_name", "AntWeb_Genus_species", "Current_status", "AntWeb_Mismatch_name")] %>% 
  rename(Current_name_AntWeb_database = Current_name) %>% 
  rename(Specimen_measured_AntWeb_name = AntWeb_Genus_species) %>% 
  rename(Specimen_measured_AntWeb_name_Mismatch = AntWeb_Mismatch_name) %>% 
  mutate(SpecimenCode = str_to_upper(SpecimenCode))
Trait_database <- left_join(Trait_database, AntWeb_database_extract, by = c("Specimen_measured_Code" = "SpecimenCode"))

# Update Current names with new matches
Trait_database$Current_name_Check[is.na(Trait_database$Current_name_Check)] <- Trait_database$Current_name_AntWeb_database[is.na(Trait_database$Current_name_Check)]
Trait_database <- Trait_database %>% 
  dplyr::select(-Current_name_AntWeb_database)

View(Trait_database)

## Complete data for specimen that are introduced, thus not in curated database, but still have a specimen code and can be found in AntWeb
Trait_data_Missing_specimens_Codes <- Trait_database$Specimen_measured_Code[is.na(Trait_database$Current_name_Check)]

Missing_vouchers_database <- AntWeb_database %>% # Non curated AntWeb database
  filter(SpecimenCode %in% str_to_lower(Trait_data_Missing_specimens_Codes))

Missing_vouchers_database_curated <- Missing_vouchers_database %>% 
  mutate(Current_name_Missing = paste0(genus, "_", species)) %>%
  mutate(Specimen_measured_AntWeb_name_Missing = Current_name_Missing) %>%
  rename(Current_status_Missing = status) %>% 
  mutate(Specimen_measured_AntWeb_name_Mismatch_Missing = Current_name_Missing != Specimen_measured_AntWeb_name_Missing) %>% 
  mutate(SpecimenCode = str_to_upper(SpecimenCode))

Trait_database <- Trait_database %>% 
  left_join(y = Missing_vouchers_database_curated[, c("SpecimenCode", "Current_name_Missing", "Specimen_measured_AntWeb_name_Missing", "Current_status_Missing", "Specimen_measured_AntWeb_name_Mismatch_Missing")], by = c("Specimen_measured_Code" = "SpecimenCode"))

Trait_database$Current_name_Check[is.na(Trait_database$Current_name_Check)] <- Trait_database$Current_name_Missing[is.na(Trait_database$Current_name_Check)] 
Trait_database$Specimen_measured_AntWeb_name[is.na(Trait_database$Specimen_measured_AntWeb_name)] <- Trait_database$Specimen_measured_AntWeb_name_Missing[is.na(Trait_database$Specimen_measured_AntWeb_name)] 
Trait_database$Current_status[is.na(Trait_database$Current_status)] <- Trait_database$Current_status_Missing[is.na(Trait_database$Current_status)] 
Trait_database$Specimen_measured_AntWeb_name_Mismatch[is.na(Trait_database$Specimen_measured_AntWeb_name_Mismatch)] <- Trait_database$Specimen_measured_AntWeb_name_Mismatch_Missing[is.na(Trait_database$Specimen_measured_AntWeb_name_Mismatch)] 

Trait_database <- Trait_database %>% 
  dplyr::select(-Current_name_Missing) %>% 
  dplyr::select(-Specimen_measured_AntWeb_name_Missing) %>%
  dplyr::select(-Current_status_Missing) %>%
  dplyr::select(-Specimen_measured_AntWeb_name_Mismatch_Missing) 

# Save curated Trait database
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")

# Load curated Trait database
Trait_database <- readRDS(file = "./input_data/Traits_data/Trait_database.rds")

### 3.4.3/ Complete data for specimen that are not in AntWeb (ex: holotypes) ####

# Update missing Current_name with Current_name_Check
View(Trait_database[is.na(Trait_database$Current_name_Check), ])
View(Trait_database[is.na(Trait_database$Current_name), ])
table(Trait_database$Current_name != Trait_database$Current_name_Check)

Trait_database$Current_name[is.na(Trait_database$Current_name)] <- Trait_database$Current_name_Check[is.na(Trait_database$Current_name)] 


# Get Current name and current status from Macroevolution taxa database
No_AntWeb_vouchers_database <- Trait_database[is.na(Trait_database$Current_name_Check), ]

# # Create current names from recorded Genus and species
# No_AntWeb_vouchers_database$Current_name <- paste0(No_AntWeb_vouchers_database$Genus, "_", No_AntWeb_vouchers_database$Species)

# Check if found in the Macroevolution taxa database
No_AntWeb_vouchers_database <- No_AntWeb_vouchers_database %>% 
  dplyr::select(-Current_status) %>%
  left_join(y = Ponerinae_Macroevolution_taxa_database[, c("Current_name", "Current_status")], join_by(Current_name))

# If no Current status is retrieved, the taxa is not in the Macroevolutionary database!
# If a Current status is retrieved, the current name can be used
table(is.na(No_AntWeb_vouchers_database$Current_status))
View(No_AntWeb_vouchers_database)

No_AntWeb_vouchers_database <- No_AntWeb_vouchers_database %>% 
  filter(!is.na(Current_status))

# Merge with Trait database
Trait_database$Current_name_Check[match(x = No_AntWeb_vouchers_database$Current_name, table = Trait_database$Current_name)] <- No_AntWeb_vouchers_database$Current_name
Trait_database$Current_status[match(x = No_AntWeb_vouchers_database$Current_name, table = Trait_database$Current_name)] <- No_AntWeb_vouchers_database$Current_status

View(Trait_database[is.na(Trait_database$Current_name_Check), ])

### 3.4.4/ Fix last errors on case by case ####

# ## Case of Specimens in Introduced location (Not in the curated AntWeb dataset)
# # CASENT0637780 # Ponera_swezeyi in AntWeb # Ponera_swezeyi in Phylogeny # Ponera_swezeyi in Current_name
# Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "Ponera_swezeyi"
# Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "Ponera_swezeyi"
# Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "valid"
# # CASENT0649763 # Brachyponera_chinensis in AntWeb # Brachyponera_chinensis in Phylogeny # Brachyponera_chinensis in Current_name 
# Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "Brachyponera_chinensis"
# Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "Brachyponera_chinensis"
# Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "valid"
# 
# ## Case of indeterminate not in the phylogeny, so need to be removed
# # CASENT0650354 # Odontomachus_(indet) not used in the phylogeny
# Trait_database <- Trait_database[!Trait_database$Specimen_measured_Code == "CASENT0650354", ]
# # CASENT0872814 # Euponera_(indet) not used in the phylogeny
# Trait_database <- Trait_database[!Trait_database$Specimen_measured_Code == "CASENT0872814", ]

# ## Case of mistake in the Specimen Code: CASENT02706126 does not exist. It is CASENT0270616
# # CASENT0270616 is Odontomachus_animosus
# Trait_database$Specimen_measured_Code[Trait_database$Genus_species == "Odontomachus_animosus"] <- "CASENT0270616"
# Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0270616"] <- "Odontomachus_animosus"
# Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT02706126"] <- "Odontomachus_animosus"
# Trait_database$Current_status[Trait_database$Specimen_measured_Code == "CASENT0270616"] <- "valid"
# Trait_database$Specimen_measured_AntWeb_name_Mismatch[Trait_database$Specimen_measured_Code == "CASENT0270616"] <- FALSE

## Case of Specimen that is not a Ponerinae

# D0872 ; CASENT0106229 # Amblyopone_australis is an Outgroup in the phylogeny. It is an Amblyoponinae.
Trait_database$Current_name_Check[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- "Amblyopone_australis"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- "Amblyopone_australis"
Trait_database$Current_status[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- "valid"
Trait_database$Specimen_measured_AntWeb_name_Mismatch[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- FALSE
Trait_database$Is_Ponerinae[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- F
Trait_database$Taxa_measurement_status[Trait_database$Specimen_measured_Code == "CASENT0106229"] <- "outgroup"

# EX1573 # Paraponera_clavata is an Outgroup in the phylogeny. It is a Paraponerinae.
Trait_database$Current_name_Check[Trait_database$Extraction_Code == "EX1573"] <- "Paraponera_clavata"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Extraction_Code == "EX1573"] <- "Paraponera_clavata"
Trait_database$Current_status[Trait_database$Extraction_Code == "EX1573"] <- "valid"
Trait_database$Specimen_measured_AntWeb_name_Mismatch[Trait_database$Extraction_Code == "EX1573"] <- FALSE
Trait_database$Is_Ponerinae[Trait_database$Extraction_Code == "EX1573"] <- F
Trait_database$Taxa_measurement_status[Trait_database$Extraction_Code == "EX1573"] <- "outgroup"

# MAMI0434 # Proceratium_google is an Outgroup in the phylogeny. It is a Proceratiinae.
Trait_database$Current_name_Check[Trait_database$Extraction_Code == "MAMI0434"] <- "Proceratium_google"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Extraction_Code == "MAMI0434"] <- "Proceratium_google"
Trait_database$Current_status[Trait_database$Extraction_Code == "MAMI0434"] <- "valid"
Trait_database$Specimen_measured_AntWeb_name_Mismatch[Trait_database$Extraction_Code == "MAMI0434"] <- FALSE
Trait_database$Is_Ponerinae[Trait_database$Extraction_Code == "MAMI0434"] <- F
Trait_database$Taxa_measurement_status[Trait_database$Extraction_Code == "MAMI0434"] <- "outgroup"

# Should be no more errors
View(Trait_database[is.na(Trait_database$Current_name_Check), ])
View(Trait_database[is.na(Trait_database$Current_name), ])

table(Trait_database$Current_name == Trait_database$Current_name_Check)

Trait_database <- Trait_database %>% 
  dplyr::select(-Current_name_Check)

### 3.4.5/ Match entries based on metadata in the phylogeny and GTL file ####

# Load updated df for Metadata of samples in the phylogeny
Phylogeny_sample_data_789t <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_789t.rds")

# Load complete list of voucher specimens
GTL_database <- openxlsx::read.xlsx(xlsxFile = "./input_data/Traits_data/Global_Taxon_List_2025_01_13.xlsx", sheet = "taxa")


## Check for match of taxa name
Trait_database$Taxa_in_phylogeny <- Trait_database$Current_name %in% Phylogeny_sample_data_789t$Current_name
table(Trait_database$Taxa_in_phylogeny)

# Add the 3 outgroups
Trait_database$Taxa_in_phylogeny[Trait_database$Taxa_measurement_status == "outgroup"] <- TRUE
table(Trait_database$Taxa_in_phylogeny)

## Check for match of Extraction code
Trait_database$Extraction_Code_In_phylogeny <- Trait_database$Extraction_Code %in% Phylogeny_sample_data_789t$Extraction_Code
table(Trait_database$Extraction_Code_In_phylogeny)

# Add the 3 outgroups
Trait_database$Extraction_Code_In_phylogeny[Trait_database$Taxa_measurement_status == "outgroup"] <- TRUE
table(Trait_database$Extraction_Code_In_phylogeny)

## Extract Specimen_Code of the extract/phylogeny voucher and Current status from phylogeny metadata
Phylogeny_metadata_to_merge <- Phylogeny_sample_data_789t %>% 
  dplyr::select(Extraction_Code, Specimen_Code, Current_status) %>% 
  rename(Current_status_to_update = Current_status) %>%
  rename(Specimen_Code_Extract_voucher = Specimen_Code) 

Trait_database <- Trait_database %>% 
  left_join(y = Phylogeny_metadata_to_merge, join_by(Extraction_Code))

# Add the 3 outgroups
Trait_database$Specimen_Code_Extract_voucher[Trait_database$Extraction_Code == "D0872"] <- "CASENT0106229" # Amblyopone_australis
Trait_database$Specimen_Code_Extract_voucher[Trait_database$Extraction_Code == "MAMI0434"] <- "CASENT0035028" # Proceratium_google
Trait_database$Specimen_Code_Extract_voucher[Trait_database$Extraction_Code == "EX1573"] <- "CASENT0633292" # Paraponera_clavata

table(is.na(Trait_database$Specimen_Code_Extract_voucher)) # 792 voucher specimens with codes

# Update Current status
table(Trait_database$Current_status == Trait_database$Current_status_to_update)
View(Trait_database[Trait_database$Current_status != Trait_database$Current_status_to_update, ])

Trait_database$Current_status[!is.na(Trait_database$Current_status_to_update)] <- Trait_database$Current_status_to_update[!is.na(Trait_database$Current_status_to_update)] 
Trait_database <- Trait_database %>% 
  dplyr::select(-Current_status_to_update)

# Add the 3 outgroups
Trait_database$Current_status[Trait_database$Taxa_measurement_status == "outgroup"] <- "valid"
table(Trait_database$Current_status)

## Extract Specimen_Code of the extract/phylogeny voucher from GTL file

Extract_metadata_to_merge <- GTL_database %>% 
  dplyr::select(ExtractionCode, voucher.specimen.code) %>% 
  rename(Extraction_Code = ExtractionCode) %>%
  rename(Specimen_Code_Extract_voucher_to_update = voucher.specimen.code) 

Trait_database$Entry_ID <- 1:nrow(Trait_database)
Trait_database <- Trait_database %>% 
  left_join(y = Extract_metadata_to_merge, join_by(Extraction_Code), na_matches = "never") %>% 
  filter(!(Specimen_Code_Extract_voucher_to_update == "CASENT0060381") | is.na(Specimen_Code_Extract_voucher_to_update)) # Error in GTL files with duplicated Extraction code

# Remove space in Extraction_code
Trait_database$Specimen_Code_Extract_voucher_to_update <- str_remove_all(string = Trait_database$Specimen_Code_Extract_voucher_to_update, pattern = " ")

Trait_database <- Trait_database %>% 
  select(-Entry_ID)

# Update Specimen_Code_Extract_voucher
table(Trait_database$Specimen_Code_Extract_voucher == Trait_database$Specimen_Code_Extract_voucher_to_update)
View(Trait_database[Trait_database$Specimen_Code_Extract_voucher != Trait_database$Specimen_Code_Extract_voucher_to_update, ])

table(is.na(Trait_database$Specimen_Code_Extract_voucher))

Trait_database$Specimen_Code_Extract_voucher[is.na(Trait_database$Specimen_Code_Extract_voucher)] <- Trait_database$Specimen_Code_Extract_voucher_to_update[is.na(Trait_database$Specimen_Code_Extract_voucher)] 
Trait_database <- Trait_database %>% 
  dplyr::select(-Specimen_Code_Extract_voucher_to_update)

table(is.na(Trait_database$Specimen_Code_Extract_voucher))

# test$EX_voucher[!(test$EX_voucher %in% Trait_database$Specimen_Code_Extract_voucher)]

## Check for match between measured specimen and voucher specimen based on Specimen code
Trait_database$Specimen_measured_Is_phylogeny_voucher <- Trait_database$Specimen_measured_Code == Trait_database$Specimen_Code_Extract_voucher
table(Trait_database$Specimen_measured_Is_phylogeny_voucher)

table(Trait_database$Specimen_measured_Is_phylogeny_voucher, Trait_database$Taxa_in_phylogeny)
View(Trait_database[Trait_database$Taxa_in_phylogeny & !Trait_database$Specimen_measured_Is_phylogeny_voucher, ])

# # # Remove columns if already present
# # Trait_database <- Trait_database %>%
# #   # select(-Specimen_phylogeny_Code, -Specimen_phylogeny_AntWeb_name, -Phylo_name)
# #   dplyr::select(-Specimen_phylogeny_Specimen_code, -Specimen_phylogeny_AntWeb_name, -Phylo_name)
# 
# # Add info on Specimen in phylogeny
# Macroevol_database_extract <- Ponerinae_Macroevolution_taxa_database[, c("Current_name", "Specimen_phylogeny_Extraction_code", "Specimen_phylogeny_Specimen_code", "Specimen_phylogeny_AntWeb_name", "Phylo_name")]
# Trait_database <- left_join(Trait_database, Macroevol_database_extract, by = c("Current_name" = "Current_name"))


### 3.4.6/ Detect conflicts of "measured vouchers" with "extract voucher"

## Detect measured specimen which are extract vouchers of another extract that the one they are attributed too

# Load complete list of voucher specimens
GTL_database <- openxlsx::read.xlsx(xlsxFile = "./input_data/Traits_data/Global_Taxon_List_2025_01_13.xlsx", sheet = "taxa")

## Based on taxa in the phylogeny
# Without outgroups
Trait_database$Specimen_measured_Is_any_phylogeny_voucher <- Trait_database$Specimen_measured_Code %in% Phylogeny_sample_data_789t$Specimen_Code
# With outgroups
Trait_database$Specimen_measured_Is_any_phylogeny_voucher <- Trait_database$Specimen_measured_Code %in% c(Phylogeny_sample_data_789t$Specimen_Code, "CASENT0106229", "CASENT0035028", "CASENT0633292")

## Based on all taxa
Trait_database$Specimen_measured_Is_any_phylogeny_voucher <- Trait_database$Specimen_measured_Code %in% c(GTL_database$voucher.specimen.code, Phylogeny_sample_data_789t$Specimen_Code, "CASENT0106229", "CASENT0035028", "CASENT0633292")

table(Trait_database$Specimen_measured_Is_any_phylogeny_voucher, Trait_database$Specimen_measured_Is_phylogeny_voucher)

Trait_database$Specimen_measured_Is_voucher_for_another_extract <- (Trait_database$Specimen_measured_Is_any_phylogeny_voucher & !Trait_database$Specimen_measured_Is_phylogeny_voucher)

View(Trait_database[Trait_database$Specimen_measured_Is_voucher_for_another_extract, ])


### 3.4.7/ Add flags for partial and complete trait measurements ####

Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")] <- round(apply(X = Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")], MARGIN = 2, FUN = as.numeric), 3)
Trait_matrix <- Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")]

Trait_database$Complete_trait_measurements <- complete.cases(Trait_matrix)
Trait_matrix_binary <- is.na(Trait_matrix)
Trait_database$Partial_trait_measurements <- apply(X = Trait_matrix_binary, MARGIN = 1, FUN = function (x) {(sum(x) > 0) & (sum(x) < length(x))})

table(Trait_database$Taxa_in_phylogeny)
table(Trait_database$Complete_trait_measurements)
table(Trait_database$Partial_trait_measurements)

## Save curated Trait database
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")


### 3.4.8/ Detect backup for extract in the phylogeny that are lacking measurements ####

table(Trait_database$Extraction_Code_In_phylogeny, (is.na(Trait_database$Complete_trait_measurements) | !Trait_database$Complete_trait_measurements & !(Trait_database$Taxa_measurement_status %in% c("male based taxon, no workers", "queen based taxon, no workers"))))

Trait_database$Missing_measurements_In_phylogeny <- Trait_database$Extraction_Code_In_phylogeny & (is.na(Trait_database$Complete_trait_measurements) | !Trait_database$Complete_trait_measurements)
Trait_database$To_measure <- Trait_database$Extraction_Code_In_phylogeny & (is.na(Trait_database$Complete_trait_measurements) | !Trait_database$Complete_trait_measurements & !(Trait_database$Taxa_measurement_status %in% c("male based taxon, no workers", "queen based taxon, no workers")))

View(Trait_database[Trait_database$Missing_measurements_In_phylogeny, ])
View(Trait_database[Trait_database$To_measure, ])

## Reorder as Taxa in phylogeny > Taxa name C Extract in Phylogeny > Extract
# Easier to detect available backup for extract missing measurements

Trait_database <- Trait_database %>% 
  dplyr::arrange(desc(Taxa_in_phylogeny), Current_name, desc(Extraction_Code_In_phylogeny), Extraction_Code)

View(Trait_database[Trait_database$To_measure, ])

### 3.4.9/ Discriminate between "Description" and "Notes" ####

## Record comments in the Trait dataset

# All entries where the measured specimen is not the extract voucher specimen should have a description

# Export as Excel file. Write comments. Reimport. Save as .rds
openxlsx::write.xlsx(x = Trait_database, file = "./input_data/Traits_data/Trait_database.xlsx", overwrite = T)
Trait_database <- read_excel("./input_data/Traits_data/Trait_database.xlsx")

# Save updated df for Trait dataset
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")

### 3.4.10/ Reorder columns and save ####

Trait_database <- Trait_database %>% 
  dplyr::select(Date, Lab, Source, Extraction_Code, Specimen_measured_Code, Current_name, Current_status,  
                Taxa_measurement_status, Is_Ponerinae, Taxa_in_phylogeny, Extraction_Code_In_phylogeny,
                HW, HL, SL, ED, WL, PW, MtFL, PmSP, PrSP, PtSP,
                Complete_trait_measurements, Partial_trait_measurements, Missing_measurements_In_phylogeny, To_measure,
                Specimen_measured_AntWeb_name, Specimen_Code_Extract_voucher, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_any_phylogeny_voucher, Specimen_measured_Is_voucher_for_another_extract,
                Description, Notes)

# Save updated df for Trait dataset
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")
# Excel version to copy-paste in DropBox
openxlsx::write.xlsx(x = Trait_database, file = "./input_data/Traits_data/Trait_database.xlsx", overwrite = T)

Trait_database_new <- Trait_database

### 3.4.6/ Merge the curated trait dataset with the taxa-level summary df ####

names(Trait_database)

Trait_database_for_merging <- Trait_database %>% 
  filter(Is_Ponerinae & Extraction_Code_In_phylogeny) %>% 
  rename(Notes_measurements = Notes) %>% 
  rename(Description_measurements = Description) %>% 
  rename(Source_measurements = Source) %>% 
  rename(Specimen_measured_Specimen_code = Specimen_measured_Code) %>%
  rename(Specimen_measured_Extraction_code = Extraction_Code) %>%
  dplyr::select(Current_name, Specimen_measured_Extraction_code, Specimen_measured_Specimen_code, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, PmSP, PrSP, PtSP, Complete_trait_measurements, Partial_trait_measurements, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_any_phylogeny_voucher, Specimen_measured_Is_voucher_for_another_extract, Source_measurements, Taxa_measurement_status, Notes_measurements, Description_measurements)

## Merge both summary df to add other valid species in the df

# Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")
# Ponerinae_Macroevolution_taxa_database <- openxlsx::read.xlsx(xlsxFile = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx", sheet = 1)

# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>%
#   rename(`Western Palearctic` = Western.Palearctic) %>%
#   rename(`Eastern Palearctic` = Eastern.Palearctic)

### Check what fields are missing after updating the measurements data!

# Remove fields that are already there only because I am rerunning the script.
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>%
  dplyr::select(-Specimen_measured_Specimen_code, -Specimen_measured_AntWeb_name, -HW, -HL, -SL, -ED, -WL, -PW, -MtFL, -PmSP, -PrSP, -PtSP, -Complete_trait_measurements, -Partial_trait_measurements, -Specimen_measured_Is_phylogeny_voucher, -Specimen_measured_Is_any_phylogeny_voucher, -Specimen_measured_Is_voucher_for_another_extract, -Source_measurements, -Taxa_measurement_status, -Notes_measurements, -Description_measurements)
  # dplyr::select(-Specimen_measured_Extraction_code, -Specimen_measured_Specimen_code, -Specimen_measured_Name, -Specimen_measured_AntWeb_name, -Specimen_measured_AntWeb_name_Mismatch, -Specimen_measured_Is_matching_phylogeny_extract, -Specimen_measured_Is_phylogeny_voucher, -HW, -HL, -SL, -ED, -WL, -PW, -MtFL, -Complete_trait_measurements, -Partial_trait_measurements, -Source_measurements, -Notes_measurements)

nrow(Ponerinae_Macroevolution_taxa_database) # 1534 rows in the Ponerinae_Macroevolution_taxa_database including all taxa in phylogeny + other valid species name
Ponerinae_Macroevolution_taxa_database <- left_join(x = Ponerinae_Macroevolution_taxa_database, y = Trait_database_for_merging, by = c("Current_name" = "Current_name", "Specimen_phylogeny_Extraction_code" = "Specimen_measured_Extraction_code")) %>% 
  # mutate(Nesting_habits = NA) %>%
  # rename(Conservative_clade_Terminals_with_MRCA = Conservative_clade_Terminal_with_MRCA) %>%
  dplyr::select(Current_name, Current_status, In_phylogeny, 
                Specimen_phylogeny_Specimen_code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_AntWeb_name_Mismatch, Specimen_phylogeny_Extraction_code,
                Phylo_label, Phylo_name, Specimen_phylogeny_Phylo_name_Mismatch, Specimen_phylogeny_Is_introduced, Subspecies,
                Conservative_clade_Node_ID, Conservative_clade_Terminals_with_MRCA, Conservative_clade_Notes, 
                # Conservative_clade_Source,
                # Exclusive_subclades_Node_ID, 
                Exclusive_subclades_Terminals_with_MRCA, Exclusive_subclades_Notes,
                Specimen_phylogeny_Country, Type_Country, Specimen_phylogeny_Bioregion, Type_Bioregion, Available_occurrences,
                To_include_in_analyses, Notes, 
                Specimen_measured_Specimen_code, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, PmSP, PrSP, PtSP, Complete_trait_measurements, Partial_trait_measurements, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_any_phylogeny_voucher, Specimen_measured_Is_voucher_for_another_extract, Source_measurements, Taxa_measurement_status, Notes_measurements, Description_measurements,
                # Specimen_measured_Is_matching_phylogeny_extract, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Specimen_code, Specimen_measured_Extraction_code, Specimen_measured_Name, Specimen_measured_AntWeb_name, Specimen_measured_AntWeb_name_Mismatch, HW, HL, SL, ED, WL, PW, MtFL, Complete_trait_measurements, Partial_trait_measurements, Source_measurements, Measurement_status, Notes_measurements,
                Nesting_habits,
                # To_measure, Specimen_to_measure_Code, Specimen_to_measure_Location, Specimen_to_measure_Only_queens_or_males, Specimen_to_measure_Notes,
                Occurrences_nb, Occurrences_nb_with_duplicates, Afrotropics, Malagasy, Indomalaya, Australasia, Neotropics, `Western Palearctic`, `Eastern Palearctic`, Nearctic, Palearctic, Afrotropics_Malagasy)

nrow(Ponerinae_Macroevolution_taxa_database) # 1534 rows after merging


# ### 3.4.7/ Check for duplicates after curation ####
# 
# # Reorder rows per taxa and keep entry using four hierarchical criteria: 1/ Specimen measured is phylogeny voucher, 2/ Match the extract used in the phylogeny, 3/ Data obtained from specimens rather than images, 4/ At random.
# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
#   arrange(Current_name, desc(Specimen_measured_Is_phylogeny_voucher), desc(Specimen_measured_Is_matching_phylogeny_extract), desc(Source_measurements)) %>% 
#   # select(Current_name, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_matching_phylogeny_extract, Source_measurements, everything()) %>%
#   group_by(Current_name) %>% 
#   mutate(Duplicates_counter = row_number(Current_name)) %>% 
#   ungroup()
# 
# # Check visually the duplicates
# Taxa_with_duplicates <- Ponerinae_Macroevolution_taxa_database$Current_name[Ponerinae_Macroevolution_taxa_database$Duplicates_counter > 1]
# Trait_database_duplicates <- Trait_database_for_merging[Trait_database_for_merging$Current_name %in% Taxa_with_duplicates, ]
# Trait_database_duplicates <- Trait_database_duplicates %>% 
#   arrange(Current_name, desc(Specimen_measured_Is_phylogeny_voucher), desc(Specimen_measured_Is_matching_phylogeny_extract), desc(Source_measurements)) %>% 
#   # select(Current_name, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_matching_phylogeny_extract, Source_measurements, everything()) %>%
#   group_by(Current_name) %>% 
#   mutate(Duplicates_counter = row_number(Current_name)) %>% 
#   ungroup()
# View(Trait_database_duplicates)
# write.xlsx(x = Trait_database_duplicates, file = "./input_data/Traits_data/Trait_database_duplicates.xlsx")
# 
# # Remove duplicates
# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
#   # filter(Duplicates_counter == 1 | Current_duplicate == T) %>%  # Remove duplicates from the trait dataset, but not duplicates from the phylogeny
#   filter(Duplicates_counter == 1) %>%  # Remove duplicates from the trait dataset
#   dplyr::select(-Duplicates_counter)
# 
# nrow(Ponerinae_Macroevolution_taxa_database)
# # 1534 rows again after cleaning of duplicates # If difference, it is because new taxa in the trait measurements need to be updated from the current AntWeb database I am using...

## Save taxa-level summary df
# saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

## Export in Excel to provide in AoW folder
# openxlsx::write.xlsx(x = Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")

## Load from Excel if manual modifications 
# Ponerinae_Macroevolution_taxa_database <- openxlsx::read.xlsx(xlsxFile = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")


### 3.5/ Retrieve grafting information from DropBox file ####

Ponerinae_Macroevolution_taxa_database_DropBox_version <- read_excel("input_data/Ponerinae_Macroevolution_taxa_database_DropBox_version.xlsx")

View(Ponerinae_Macroevolution_taxa_database_DropBox_version[which(!Ponerinae_Macroevolution_taxa_database_DropBox_version$Current_name %in% Ponerinae_Macroevolution_taxa_database$Current_name), ])
View(Ponerinae_Macroevolution_taxa_database[which(!Ponerinae_Macroevolution_taxa_database$Current_name %in% Ponerinae_Macroevolution_taxa_database_DropBox_version$Current_name), ])

# Extract grafting info from Dropbox file
Grafting_info <- Ponerinae_Macroevolution_taxa_database_DropBox_version %>% 
  select(Current_name, Conservative_clade_Node_ID, Conservative_clade_Terminals_with_MRCA, Conservative_clade_Notes, Conservative_clade_Source,
         Exclusive_subclades_Node_ID, Exclusive_subclades_Terminals_with_MRCA, Exclusive_subclades_Notes)

# Remove previous field to be replaced by the update
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  select(-Conservative_clade_Node_ID, -Conservative_clade_Terminals_with_MRCA, -Conservative_clade_Notes, -Conservative_clade_Source,
         -Exclusive_subclades_Node_ID, -Exclusive_subclades_Terminals_with_MRCA, -Exclusive_subclades_Notes)

# Merge update
Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, Grafting_info)

# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
#   filter(Duplicates_counter == 1) %>% 
#   select(-Duplicates_counter)


### 3.6/ Retrieve "To_measure" information from DropBox file ####

# ## From extra sheet
# 
# To_measure_info <- read_excel("input_data/Ponerinae_Macroevolution_taxa_database_DropBox_version.xlsx", sheet = "to measure")
# 
# # Rename columns adequately
# To_measure_info <- To_measure_info %>% 
#   rename(Specimen_to_measure_Code = `specimen to measure`) %>%
#   rename(Specimen_to_measure_Location = "where") %>%
#   rename(Specimen_to_measure_Notes = "note") %>%
#   select(Current_name, Specimen_to_measure_Code, Specimen_to_measure_Location, Specimen_to_measure_Notes)
# 
# # Merge with Macroevolution database
# Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, To_measure_info)
# 
# ## Add field to detect that measurements are still needed
# 
# Ponerinae_Macroevolution_taxa_database$To_measure <- is.na(Ponerinae_Macroevolution_taxa_database$Specimen_measured_Code)
# 
# table(Ponerinae_Macroevolution_taxa_database$To_measure)
# 
# ## Add field to notify if only 'queens or males' are available, disqualifying them from being measured
# 
# # Based on Jack's saying, there should be all specimens that have no measurements, but that are not listed in its subsheet
# 
# Ponerinae_Macroevolution_taxa_database$Specimen_to_measure_Only_queens_or_males <- NA
# 
# In_Jack_list <- Ponerinae_Macroevolution_taxa_database$Current_name %in% To_measure_info$Current_name
# 
# test <- !In_Jack_list & Ponerinae_Macroevolution_taxa_database$To_measure
# table(test) # 80 cases of taxa that are yet to be measured and not in Jack's list ! 
# View(Ponerinae_Macroevolution_taxa_database[test, ])

## From DropBox info entered manually

table(Ponerinae_Macroevolution_taxa_database$To_measure)
table(Ponerinae_Macroevolution_taxa_database$To_measure, Ponerinae_Macroevolution_taxa_database$In_phylogeny)

# Ponerinae_Macroevolution_taxa_database$To_measure <- is.na(Ponerinae_Macroevolution_taxa_database$Specimen_measured_Specimen_code)
Ponerinae_Macroevolution_taxa_database$To_measure <- is.na(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements) | !Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements

# Extract grafting info from Dropbox file
To_measure_info <- Ponerinae_Macroevolution_taxa_database_DropBox_version %>% 
  select(Current_name, To_measure, Specimen_to_measure_Code,	Specimen_to_measure_Location,	Specimen_to_measure_Only_queens_or_males,	Specimen_to_measure_Notes)

# Remove previous field to be replaced by the update
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  select(-Specimen_to_measure_Code,	-Specimen_to_measure_Location, -Specimen_to_measure_Only_queens_or_males,	-Specimen_to_measure_Notes)

# Merge update
Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, To_measure_info)

### 3.7/ Correct the Neoponera_bucki issue ####

table(Ponerinae_Macroevolution_taxa_database$Current_name == "NewGenus_bucki")
table(Ponerinae_Macroevolution_taxa_database$Current_name == "Neoponera_bucki")

View(Ponerinae_Macroevolution_taxa_database[Ponerinae_Macroevolution_taxa_database$Current_name %in% c("NewGenus_bucki", "Neoponera_bucki"), ])

# Remove Neoponera_bucki entry
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database[!(Ponerinae_Macroevolution_taxa_database$Current_name == "Neoponera_bucki"), ]

# Rename NewGenus_bucki as Neoponera_bucki
Ponerinae_Macroevolution_taxa_database$Current_name[Ponerinae_Macroevolution_taxa_database$Current_name == "NewGenus_bucki"] <- "Neoponera_bucki"

### 3.8/ Export Macroevolution database = taxa-level summary df ####

# Sort columns in a logical order
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  select(Current_name, Current_status, In_phylogeny, Specimen_phylogeny_Specimen_code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_AntWeb_name_Mismatch, Specimen_phylogeny_Extraction_code,
         Phylo_label, Phylo_name, Specimen_phylogeny_Phylo_name_Mismatch, Specimen_phylogeny_Is_introduced, Subspecies,
         Conservative_clade_Node_ID, Conservative_clade_Terminals_with_MRCA, Conservative_clade_Notes,
         Exclusive_subclades_Node_ID, Exclusive_subclades_Terminals_with_MRCA, Exclusive_subclades_Notes,
         Specimen_phylogeny_Country, Type_Country, Specimen_phylogeny_Bioregion, Type_Bioregion, Available_occurrences, To_include_in_analyses, Notes,
         Specimen_measured_Specimen_code, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, PmSP, PrSP, PtSP, Complete_trait_measurements, Partial_trait_measurements,
         Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Is_any_phylogeny_voucher, Specimen_measured_Is_voucher_for_another_extract, Source_measurements, Taxa_measurement_status, Notes_measurements, Description_measurements,
         # Specimen_measured_Is_matching_phylogeny_extract, Specimen_measured_Is_phylogeny_voucher, Specimen_measured_Specimen_code, Specimen_measured_Extraction_code,  Specimen_measured_Name, Specimen_measured_AntWeb_name, Specimen_measured_AntWeb_name_Mismatch,
         # HW, HL, SL, ED, WL, PW, MtFL, Complete_trait_measurements, Partial_trait_measurements, Source_measurements, Notes_measurements,
         # To_measure, Specimen_to_measure_Code, Specimen_to_measure_Location, Specimen_to_measure_Only_queens_or_males, Specimen_to_measure_Notes
         Occurrences_nb, Occurrences_nb_with_duplicates, Afrotropics, Malagasy, Indomalaya, Australasia, Neotropics, `Western Palearctic`, `Eastern Palearctic`, Nearctic, Palearctic, Afrotropics_Malagasy
         )

## Save taxa-level summary df
saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

## Export in Excel to provide in AoW folder
openxlsx::write.xlsx(x = Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")


table(Ponerinae_Macroevolution_taxa_database$In_phylogeny, Ponerinae_Macroevolution_taxa_database$To_measure)




##### 4/ Generate clean version for Supplementary Data 3 #####

## Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

Ponerinae_Macroevolution_taxa_database_cleaned <- Ponerinae_Macroevolution_taxa_database %>% 
  mutate(Taxa_ID = 1:nrow(Ponerinae_Macroevolution_taxa_database)) %>%
  mutate(Grafted = !In_phylogeny) %>%
  rename(Taxa_name = Current_name,
         Specimen_code_Phylogeny_voucher = Specimen_phylogeny_Specimen_code,
         Extraction_code_Phylogeny_voucher = Specimen_phylogeny_Extraction_code,
         Grafting_clade_Terminals_with_MRCA = Conservative_clade_Terminals_with_MRCA,
         Grafting_clade_Notes = Conservative_clade_Notes,
         Excluded_subclades_Terminals_with_MRCA = Exclusive_subclades_Terminals_with_MRCA,
         Excluded_subclades_Notes = Exclusive_subclades_Notes) %>%
  dplyr::select(Taxa_ID, Taxa_name, Current_status,
                Occurrences_nb, Occurrences_nb_with_duplicates,
                Grafted, Specimen_code_Phylogeny_voucher, Extraction_code_Phylogeny_voucher,
                Grafting_clade_Terminals_with_MRCA, Grafting_clade_Notes,
                Excluded_subclades_Terminals_with_MRCA, Excluded_subclades_Notes)

## Save Supplementary Data 3
saveRDS(Ponerinae_Macroevolution_taxa_database_cleaned, file = "./input_data/Ponerinae_Macroevolution_taxa_database_cleaned.rds")
## Export in Excel
openxlsx::write.xlsx(x = Ponerinae_Macroevolution_taxa_database_cleaned, file = "./input_data/Ponerinae_Macroevolution_taxa_database_cleaned.xlsx")


##### 5/ Evaluate database #####

Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")


# Evaluate trait measurement collection
Ponerinae_Macroevolution_taxa_database[replace_na(data = !Ponerinae_Macroevolution_taxa_database$In_phylogeny & Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements, replace = FALSE), ]

table(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements)
table(Ponerinae_Macroevolution_taxa_database$In_phylogeny)
table(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements, Ponerinae_Macroevolution_taxa_database$In_phylogeny)
table(replace_na(data = Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements[Ponerinae_Macroevolution_taxa_database$In_phylogeny], replace = F)) 
table(replace_na(data = Ponerinae_Macroevolution_taxa_database$Partial_trait_measurements[Ponerinae_Macroevolution_taxa_database$In_phylogeny], replace = F)) 

# Evaluate grafting information

table(!is.na(Ponerinae_Macroevolution_taxa_database$Conservative_clade_Terminals_with_MRCA))
table(!is.na(Ponerinae_Macroevolution_taxa_database$Conservative_clade_Terminals_with_MRCA), Ponerinae_Macroevolution_taxa_database$In_phylogeny)

