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
# AntWeb specimen database accessed on 04 February 2024
# AntMap: GABI Database Release 1.0 from 18 January 2020

###

### Sources

# AntWeb. Version 8.101. California Academy of Science, online at https://www.antweb.org. Accessed on 04 February 2024.
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
AntWeb_database <- read_excel("./input_data/AntWeb_data/AntWeb_database_Ponerinae_2024_02_04.xlsx")

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

# Extract labels from the phylogeny
taxa_names_phylogeny_789t_label <- Ponerinae_phylogeny_789t$tip.label

# Extract data for specimens used in the phylogeny, including name in AntWeb database (may have changed)
Phylogeny_sample_data <- read_excel("input_data/Phylogenies/ponerinae-dataset-sample-list-current.xlsx")
Phylogeny_sample_data <- Phylogeny_sample_data %>% 
  filter(`Ponerinae-dataset` == "KEEP") %>%
  select(`Sample name`, `Extraction code`, `Specimen code`, `New sample name`, `no change?`, country, bioregion) %>% 
  rename("AntWeb_name" = `Sample name`) %>% 
  mutate("AntWeb_name" = str_to_title(AntWeb_name)) %>% 
  rename("Phylo_label" = `New sample name`) %>%
  rename("Name_update" = `no change?`) %>% 
  mutate("Name_update" = !(Name_update))

# Remove outgroups
Phylogeny_sample_data <- Phylogeny_sample_data[!(Phylogeny_sample_data$AntWeb_name %in% c("Amblyopone_australis", "Paraponera_clavata", "Proceratium_google")), ]

# Correct mistakes in AntWeb Names
Phylogeny_sample_data$AntWeb_name[Phylogeny_sample_data$AntWeb_name == "Cryptopone_guatemalensis_large"] <- "Cryptopone_guatemalensis"

## 1.3.1/ Extract phylogeny name from label ####

# Extract taxa names from the phylogeny
taxa_in_phylogeny <- str_split(string = Phylogeny_sample_data$Phylo_label, pattern = "_", simplify = F)

taxa_in_phylogeny <- lapply(taxa_in_phylogeny, FUN = function (x) { y <- x[-length(x)] ; return(y) } )
taxa_in_phylogeny <- unlist(lapply(taxa_in_phylogeny, FUN = paste0, collapse = "_"))

taxa_in_phylogeny <- str_remove(taxa_in_phylogeny, pattern = "_MAMI.*")
taxa_in_phylogeny <- str_remove(taxa_in_phylogeny, pattern = "MAMI.*")
taxa_in_phylogeny <- str_remove(taxa_in_phylogeny, pattern = "_BBX.*")

Phylogeny_sample_data$Phylo_name <- str_to_title(taxa_in_phylogeny)

## Transform format to match AntWeb names

Phylo_name_formated <- Phylogeny_sample_data$Phylo_name
AntWeb_name_formated <- Phylogeny_sample_data$AntWeb_name

# Adjust Phylo_names
# Change underscore in specific epithet name for "-"
Phylo_name_formated_genus <- str_split(string = Phylo_name_formated, pattern = "_", simplify = T)[,1]
Phylo_name_formated_sp_epithet <- str_split(string = Phylo_name_formated, pattern = "_", simplify = F)
Phylo_name_formated_sp_epithet <- lapply(X = Phylo_name_formated_sp_epithet, FUN = function (x) { y <- x[-1] ; return(y) })
Phylo_name_formated_sp_epithet <- unlist(lapply(X = Phylo_name_formated_sp_epithet, paste0, collapse = "-"))

View(str_split(string = Phylo_name_formated, pattern = "_", simplify = T))
# Exception for "_01", "_1", "_a", "_b", "_cf", "_cf2", "_group", "_madecassa", "_marleyi", "_nr*", "_saldanhae", "_sp*".
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-01", replacement = "_01")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-1", replacement = "_1")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-a$", replacement = "_a")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-b$", replacement = "_b")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-cf", replacement = "_cf")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-group", replacement = "_group")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-madecassa", replacement = "_madecassa")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-marleyi", replacement = "_marleyi")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-nr", replacement = "_nr")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-saldanhae", replacement = "_saldanhae")
Phylo_name_formated_sp_epithet <- str_replace(string = Phylo_name_formated_sp_epithet, pattern = "-sp", replacement = "_sp")

# Merge Genus and specific epithets
Phylo_name_formated <- paste0(Phylo_name_formated_genus, "_", Phylo_name_formated_sp_epithet)
Phylogeny_sample_data$Phylo_name <- Phylo_name_formated

# Adjust AntWeb_names
# Change underscore in specific epithet name for "-"
AntWeb_name_formated_genus <- str_split(string = AntWeb_name_formated, pattern = "_", simplify = T)[,1]
AntWeb_name_formated_sp_epithet <- str_split(string = AntWeb_name_formated, pattern = "_", simplify = F)
AntWeb_name_formated_sp_epithet <- lapply(X = AntWeb_name_formated_sp_epithet, FUN = function (x) { y <- x[-1] ; return(y) })
AntWeb_name_formated_sp_epithet <- unlist(lapply(X = AntWeb_name_formated_sp_epithet, paste0, collapse = "-"))

View(str_split(string = AntWeb_name_formated, pattern = "_", simplify = T))
# Exception for "_1", "_a", "_b", "_cf", "_cf2", "_group", "_madecassa",  "_marleyi", "_nr*", "_saldanhae", "_sp*".
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-1", replacement = "_1")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-a$", replacement = "_a")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-b$", replacement = "_b")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-cf", replacement = "_cf")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-group", replacement = "_group")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-madecassa", replacement = "_madecassa")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-marleyi", replacement = "_marleyi")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-nr", replacement = "_nr")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-saldanhae", replacement = "_saldanhae")
AntWeb_name_formated_sp_epithet <- str_replace(string = AntWeb_name_formated_sp_epithet, pattern = "-sp", replacement = "_sp")
# Merge Genus and specific epithets
AntWeb_name_formated <- paste0(AntWeb_name_formated_genus, "_", AntWeb_name_formated_sp_epithet)
Phylogeny_sample_data$AntWeb_name <- AntWeb_name_formated


# Save metadata of 802 ingroup taxa in the phylogeny
saveRDS(Phylogeny_sample_data, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")

# Save pruned phylogeny with only ingroups
saveRDS(Ponerinae_phylogeny_789t, file = "./input_data/Phylogenies/Ponerinae_phylogeny_789t_treedata.rds")


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
AntWeb_database$Genus_species <- paste0(str_to_title(AntWeb_database$genus), "_", AntWeb_database$species)
Ponerinae_Species_AntWeb <- unique(AntWeb_database$Genus_species[(AntWeb_database$status == "valid") & !AntWeb_database$is_introduced])

## Obtain valid Ponerinae species names from AntCat
Ponerinae_Species_AntCat_df <- read_excel("./input_data/AntWeb_data/AntCat_names_Ponerinae_2023_11_23.xlsx", sheet = "Ponerinae_Species")
Ponerinae_Species_AntCat_df$Genus_species <- paste0(Ponerinae_Species_AntCat_df$genus, "_", Ponerinae_Species_AntCat_df$species)
Ponerinae_Species_AntCat <- Ponerinae_Species_AntCat_df$Genus_species

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
# Flag against AntWeb names + everything starting with Poner*

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

## 2.3.1/ Filter out fossils, indet and introduced entries ####

## Create Genus_species field
AntWeb_database_curated$Genus_species <- paste0(AntWeb_database_curated$genus, "_", AntWeb_database_curated$species)

## Keep only valid names and morphotaxa
table(AntWeb_database$status)

AntWeb_database_curated <- AntWeb_database[AntWeb_database$status %in% c("valid", "morphotaxon"), ]

## Filter out fossils
# All valid taxa in AntWeb that are not in the AntCat list are fossils. Need to remove them from AntWeb by filtering to keep only AntCat taxa
fossils_names <- setdiff(Ponerinae_Species_AntWeb, Ponerinae_Species_AntCat)
fossils_names <- fossils_names[!is.na(fossils_names)]; fossils_names 

AntWeb_database_curated <- AntWeb_database_curated[!(AntWeb_database_curated$Genus_species %in% fossils_names), ]

## Remove entries for introduced species
table(AntWeb_database_curated$is_introduced)
AntWeb_database_curated <- AntWeb_database_curated[AntWeb_database_curated$is_introduced == F, ]

## Adjust status of valid names from AntWeb to match valid names from AntCat
Ponerinae_valid_Species_AntWeb <- unique(AntWeb_database_curated$Genus_species[AntWeb_database_curated$status %in% c("valid")])
Ponerinae_morphotaxa_Species_AntWeb <- unique(AntWeb_database_curated$Genus_species[AntWeb_database_curated$status %in% c("morphotaxon")])

Ponerinae_valid_Species_AntCat_in_AntWeb <- Ponerinae_valid_Species_AntWeb[Ponerinae_valid_Species_AntWeb %in% Ponerinae_Species_AntCat]
Ponerinae_not_valid_Species_AntCat_valid_in_AntWeb <- Ponerinae_valid_Species_AntWeb[!(Ponerinae_valid_Species_AntWeb %in% Ponerinae_Species_AntCat)]
# All valid species in Antweb are valid in AntCat

Ponerinae_morphotaxa_AntWeb_valid_in_AntCat <- Ponerinae_morphotaxa_Species_AntWeb[Ponerinae_morphotaxa_Species_AntWeb %in% Ponerinae_Species_AntCat]
# One morphotaxon that is valid in AntCat: "Bothroponera_pumicosa sculpturata_nr". It is a subspecies of a valid species. Not a mistake.
Ponerinae_morphotaxa_AntWeb_not_valid_AntCat <- Ponerinae_morphotaxa_Species_AntWeb[!(Ponerinae_morphotaxa_Species_AntWeb %in% Ponerinae_Species_AntCat)]

## 2.3.2/ Curate names so they match between AntWeb, the phylogeny, and this database ####

Phylogeny_sample_data_802 <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")
taxa_names_phylogeny_789t_AntWeb_format <- Phylogeny_sample_data_802$AntWeb_name

AntWeb_database_taxa_names <- unique(AntWeb_database_curated$Genus_species)

# Detect Phylo-taxa missing from AntWeb
taxa_names_phylogeny_789t_AntWeb_format[!(taxa_names_phylogeny_789t_AntWeb_format %in% AntWeb_database_taxa_names)]

## Need to adjust case-by-case. Use Specimen code to track them back.

## Cases with wrong name in Ponerinae-dataset-sample-list-current => Adjust to right name to match the AntWeb data
# Specimen CASENT0381057 is Anochetus_altisquamis as labeled in the phylogeny. Not Anochetus_ag02
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_ag02"] <- "Anochetus_altisquamis"
# Specimen CASENT0906625 is Anochetus_ghilianii as in the phylogeny. Not "Anochetus_jlrl-nard". Not sure if other specimens should be changed too. In doubt, since they are all from Morocco, do no change. 
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_jlrl-nard"] <- "Anochetus_ghilianii"
# Specimen CASENT0645962 is Corrieopone_nouragues as labeled in the phylogeny. Not an indeterminate "Ponerinae"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponerinae_"] <- "Corrieopone_nouragues"
# Specimen CASENT0882091 is Euponera_sjostedti as labeled in the phylogeny. Not Bothroponera_picardi_cf
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Bothroponera_picardi_cf"] <- "Euponera_sjostedti"
# Specimen CASENT0217039 is Fisheropone_hartwigi as labeled in the phylogeny. Not Cryptopone_hartwigi
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_hartwigi"] <- "Fisheropone_hartwigi"
# Specimen CASENT0650367 is Hypoponera_pruinosa as labeled in the phylogeny. Not Hypoponera_newbrit01
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_newbrit01"] <- "Hypoponera_pruinosa"
# Specimen INB0003695621 is Mayaponera_pergandei as labeled in the phylogeny. Not Rasopone_pergandei
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_pergandei"] <- "Mayaponera_pergandei"
# Specimen CASENT0649888 is Neoponera_magnifica-4 in Antweb. Labeled as Neoponera_magnifica4 in the phylogeny.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_magnifica4"] <- "Neoponera_magnifica-4"
# Specimen CASENT0650355 is Odontomachus_saevissimus in Antweb. Wrongly labeled as Odontomachus_saevissmus in the phylogeny.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Odontomachus_saevissmus"] <- "Odontomachus_saevissimus"
# Specimen CASENT0646046 is Pachycondyla_procidua in Antweb and in the phylogeny. Not Neoponera_procidua (former name).
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_procidua"] <- "Pachycondyla_procidua"
# Specimen CASENT0761219 is Ponera_sinensis_nr in Antweb. Labeled as Ponera_chinensis_cf in the phylogeny.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponera_chinensis_cf"] <- "Ponera_sinensis_nr"
# Specimen CASENT0611162 is Ponera_exotica as labeled in the phylogeny. Not Ponera_jtl002
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponera_jtl002"] <- "Ponera_exotica"
# Specimen INB0003660648 is Rasopone_cryptergates as labeled in the phylogeny. Not Rasopone_jtl017
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_jtl017"] <- "Rasopone_cryptergates"
# Specimen CASENT0633282 is Rasopone_cubitalis as labeled in the phylogeny. Not Rasopone_jtl018
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_jtl018"] <- "Rasopone_cubitalis"
# Specimen JTLC000015360 is Rasopone_mesoamericana as labeled in the phylogeny. Not Rasopone_jtl022
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "JTLC000015360"] <- "Rasopone_mesoamericana"
# Specimen CASENT0633053 is Rasopone_jtl029 almost as labeled in the phylogeny. Not Rasopone_jtl022
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0633053"] <- "Rasopone_jtl029"
# Specimen CASENT0633075 is Rasopone_jtl030 almost as labeled in the phylogeny. Not Rasopone_jtl022
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0633075"] <- "Rasopone_jtl030"
# Specimen INB0004099727 is Rasopone_pluviselva as labeled in the phylogeny. Not Rasopone_jtl016
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_jtl016"] <- "Rasopone_pluviselva"
# Specimen CASENT0637779 is Wadeura_holmgrenita. Also wrongly labeled as Wadeura_holmgreni in the phylogeny. Not Cryptopone_holmgreni
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_holmgreni"] <- "Wadeura_holmgrenita"
# Specimen INB0003694616 is Wadeura_guianensis as labeled in the phylogeny. Not Cryptopone_guianensis
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_guianensis"] <- "Wadeura_guianensis"
# Specimen CASENT0637806 is Wadeura_pauli as labeled in the phylogeny. Not Rasopone_jtl016
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_jtl001"] <- "Wadeura_pauli"
# Specimen CASENT0845817 is Anochetus_bequaerti as labeled in the phylogeny. Not Anochetus_cm05
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_cm05"] <- "Anochetus_bequaerti"
# Specimen CASENT0287021 is Anochetus_graeffei as labeled in the phylogeny. Not Anochetus_th02
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_th02"] <- "Anochetus_graeffei"
# Specimen CASENT0387827 is Anochetus_rugosus as labeled in the phylogeny. Not Anochetus_my07
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_my07"] <- "Anochetus_rugosus"
# Specimen CASENT0372252 is Anochetus_targionii as labeled in the phylogeny. Not Anochetus_PE03
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_pe03"] <- "Anochetus_targionii"
# Specimen CASENT0641046 is Cryptopone_gilvatumida as labeled in the phylogeny. Not Cryptopone_guatemalensis-large
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_guatemalensis-large"] <- "Cryptopone_gilvatumida"
# Specimen CASENT0650356 is Diacamma_timor_01 almost as labeled in the phylogeny: Diacamma_Timor_01. Not Diacamma_timor_1
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Diacamma_timor_1"] <- "Diacamma_timor_01"
# Specimen CASENT0650200 is Ectomomyrmex_simillimus. Wrongly labeled in the phylogeny as Ectomomyrmex_Janda_sp10.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ectomomyrmex_janda_sp10"] <- "Ectomomyrmex_simillimus"
# Specimen CASENT0650200 is Hypoponera_dias19_2. Wrongly labeled in the phylogeny as Hypoponera_jtl033 because synonimized with another morphotaxon from molecular results
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_jtl033"] <- "Hypoponera_dias19_2"
# Specimen CASENT0046628 is Hypoponera_sc-mora as labeled in the phylogeny. Not Hypoponera_moraganga. May be updated later if new valid species name.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_moraganga"] <- "Hypoponera_sc-mora"
# Specimen CASENT0779195 is Hypoponera_sc-nosy as labeled in the phylogeny. Not Hypoponera_nosyankao.  May be updated later if new valid species name.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_nosyankao"] <- "Hypoponera_sc-nosy"
# Specimen CASENT0497401 is Hypoponera_sc-rano as labeled in the phylogeny. Not Hypoponera_ranomafana.  May be updated later if new valid species name.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ranomafana"] <- "Hypoponera_sc-rano"
# Specimen CASENT0872732 is Hypoponera_sc-tamp as labeled in the phylogeny. Not Hypoponera_tampolo.  May be updated later if new valid species name.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_tampolo"] <- "Hypoponera_sc-tamp"
# Specimen CASENT0152601 is Hypoponera_sc-zaha as labeled in the phylogeny. Not Hypoponera_zahamena.  May be updated later if new valid species name.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_zahamena"] <- "Hypoponera_sc-zaha"
# Specimen CASENT0373353 is Leptogenys_unistimulosa as labeled in the phylogeny. Not Leptogenys_pe01
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_pe01"] <- "Leptogenys_unistimulosa"
# Specimen CASENT0650367 is Mayaponera_pergandei as labeled in the phylogeny. Not Rasopone_pergandei
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_pergandei"] <- "Mayaponera_pergandei"
# Specimen INBIOCRI001241794 is Mayaponera_arhuaca as labeled in the phylogeny. Not Rasopone_arhuaca
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_arhuaca"] <- "Mayaponera_arhuaca"
# Specimen CASENT0633229 is Mayaponera_becculata as labeled in the phylogeny. Not Rasopone_becculata
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_becculata"] <- "Mayaponera_becculata"
# Specimen CASENT0650197 is Myopias_papua. Wrongly labeled in the phylogeny as Myopias_janda_sp4.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Myopias_janda_sp4"] <- "Myopias_papua"
# Specimen MEPNINV2923 is Neoponera_ecu2923. Wrongly labeled in the phylogeny as Neoponera_ecu2323.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_ecu2323"] <- "Neoponera_ecu2923"
# Specimen DZUP549444 is Neoponera_gojira as labeled in the phylogeny. Not Neoponera_bra549444.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_bra549444"] <- "Neoponera_gojira"
# Specimen ATPFOR1962 is Neoponera_magnifica-1. Wrongly labeled in the phylogeny as Neoponera_magnifica1.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_magnifica1"] <- "Neoponera_magnifica-1"
# Specimen ATPFOR1965 is Neoponera_magnifica-2. Wrongly labeled in the phylogeny as Neoponera_magnifica2.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_magnifica2"] <- "Neoponera_magnifica-2"
# Specimen UFV-LABECOL-000625 is Neoponera_metanotalis. Wrongly labeled in the phylogeny as Neoponera_metanotalis1.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_metanotalis1"] <- "Neoponera_metanotalis"
# Specimen ATPFOR1997 is Neoponera_metanotalis-2. Wrongly labeled in the phylogeny as Neoponera_metanotalis2.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_metanotalis2"] <- "Neoponera_metanotalis-2"
# Specimen UFV-LABECOL-000500 is Neoponera_schultzi. Wrongly labeled in the phylogeny as Neoponera_schultzi1.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_schultzi1"] <- "Neoponera_schultzi"
# Specimen CASENT0636865 is Odontomachus_erythrocephalus. Wrongly labeled in the phylogeny as Odontomachus_erythrocephala.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Odontomachus_erythrocephala"] <- "Odontomachus_erythrocephalus"
# Specimen CASENT0650350 is Odontomachus_papuanus as labeled in the phylogeny. Not Odontomachus_indet.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650350"] <- "Odontomachus_papuanus"
# Specimen CASENT0159375 is Ponera_adumbrans as labeled in the phylogeny. Not Ponera_sc02.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponera_sc02"] <- "Ponera_adumbrans"
# Specimen CASENT0249221 is Ponera_colaensis as labeled in the phylogeny. Not Ponera_coalensis.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponera_coalensis"] <- "Ponera_colaensis"
# Specimen INB0003621410 is Ponera_jtl001 as labeled in the phylogeny. Not Ponerinae_jtl001.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ponerinae_jtl001"] <- "Ponera_jtl001"
# Specimen CASENT0887836 is Pseudoneoponera_au_a as labeled in the phylogeny. Not Pseudoponera_au_a.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Pseudoponera_au_a"] <- "Pseudoneoponera_au_a"
# Specimen CASENT0887835 is Pseudoneoponera_au04 as labeled in the phylogeny. Not Pseudoponera_au04.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Pseudoponera_au04"] <- "Pseudoneoponera_au04"
# Specimen CASENT0887834 is Pseudoneoponera_au07 as labeled in the phylogeny. Not Pseudoponera_au07.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Pseudoponera_au07"] <- "Pseudoneoponera_au07"
# Specimen CASENT0887837 is Pseudoneoponera_au08 as labeled in the phylogeny. Not Pseudoponera_au08.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Pseudoponera_au08"] <- "Pseudoneoponera_au08"
# Specimen CASENT0644253 is Rasopone_costaricensis as labeled in the phylogeny. Not Rasopone_jtl014.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0644253"] <- "Rasopone_costaricensis"
# Specimen CASENT0633224 is Rasopone_guatemalensis as labeled in the phylogeny. Not Rasopone_jtl014.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0633224"] <- "Rasopone_guatemalensis"
# Specimen CASENT0612523 is Rasopone_politognatha as labeled in the phylogeny. Not Rasopone_jtl014.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0612523"] <- "Rasopone_politognatha"
# Specimen CASENT0640257 is Rasopone_ferruginea as labeled in the phylogeny. Not Rasopone_jtl033.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name ==  "Rasopone_jtl033"] <- "Rasopone_ferruginea"
# Specimen CASENT0611483 is Rasopone_minuta as labeled in the phylogeny. Not Rasopone_jtl025.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name ==  "Rasopone_jtl025"] <- "Rasopone_minuta"
# Specimen CASENT0614335 is Rasopone_subcubitalis as labeled in the phylogeny. Not Rasopone_jtl024.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name ==  "Rasopone_jtl024"] <- "Rasopone_subcubitalis"
# Specimen CASENT0914502 is Simopelta_vieirai in Antweb. Wrongly labeled as Simopelta_vierirai in the phylogeny.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Simopelta_vierirai"] <- "Simopelta_vieirai"
# Specimen CASENT0639794 is Leptogenys_ug01 as labeled in the phylogeny. Not Leptogenys_elegans.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_elegans"] <- "Leptogenys_ug01"
# Specimen CASENT0779894 is Anochetus_afr01 as labeled in the phylogeny. Not Anochetus_casc-mz03.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_casc-mz03"] <- "Anochetus_afr01"
# Specimen CASENT0919893 is Hypoponera_javana as labeled in the phylogeny. Not Hypoponera_confinis.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_confinis"] <- "Hypoponera_javana"
# Specimen CASENT0270584 is Leptogenys_mjobergi_nr as labeled in the phylogeny. Not Leptogenys_ebenina_nr.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_ebenina_nr"] <- "Leptogenys_mjobergi_nr"
# Specimen CASENT0384495 is Leptogenys_myops as labeled in the phylogeny. Not Leptogenys_my01.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_my01"] <- "Leptogenys_myops"
# Specimen CASENT0650362 is Ectomomyrmex_acutus as labeled in the phylogeny. Not Ectomomyrmex_bg01
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650362"] <- "Ectomomyrmex_acutus"
# Specimen CASENT0650199 is Ectomomyrmex_janda_sp7 as labeled in the phylogeny. Not Ectomomyrmex_acutus
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650199"] <- "Ectomomyrmex_janda_sp7"
Phylogeny_sample_data_802$Name_status[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650199"] <- "morphotaxon"
# Specimen MCZ-ENT00759860 is Ectomomyrmex_zhengi in AntWeb. Found in Yunnan, China. Possibly wrongly labeled in the phylogeny as Ectomomyrmex_obtusus which is found on a small island in Indonesia.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "MCZ-ENT00759860"] <- "Ectomomyrmex_zhengi"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "MCZ-ENT00759860"] <- "Ectomomyrmex_zhengi"
Phylogeny_sample_data_802$Name_changed_from_phylogeny[Phylogeny_sample_data_802$`Specimen code` == "MCZ-ENT00759860"] <- T



## Cases of subspecies. Need to change the name in AntWeb to match the subspecies
Phylogeny_sample_data_802$Subspecies <- F # Create new field to record specific cases of Subspecies
# Specimen CASENT0888244 is Hagensia_havilandi_marleyi. Change in AntWeb and flag for special case of subspecies
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hagensia_havilandi" & AntWeb_database_curated$subspecies == "marleyi"] <- "havilandi_marleyi"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hagensia_havilandi" & AntWeb_database_curated$subspecies == "marleyi"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hagensia_havilandi" & AntWeb_database_curated$subspecies == "marleyi"] <- "Hagensia_havilandi_marleyi"
Phylogeny_sample_data_802$Subspecies[Phylogeny_sample_data_802$AntWeb_name == "Hagensia_havilandi_marleyi"] <- T
# Specimen CASENT0300040 is Parvaponera_darwinii_madecassa. Change in AntWeb and flag for special case of subspecies
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Parvaponera_darwinii" & AntWeb_database_curated$subspecies == "madecassa"] <- "darwinii_madecassa"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Parvaponera_darwinii" & AntWeb_database_curated$subspecies == "madecassa"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Parvaponera_darwinii" & AntWeb_database_curated$subspecies == "madecassa"] <- "Parvaponera_darwinii_madecassa"
Phylogeny_sample_data_802$Subspecies[Phylogeny_sample_data_802$AntWeb_name == "Parvaponera_darwinii_madecassa"] <- T
# Specimen CASENT0787946 is Hagensia_peringueyi_saldanhae. Keep like this, but flag for special case of subspecies
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hagensia_peringueyi" & AntWeb_database_curated$subspecies == "saldanhae"] <- "peringueyi_saldanhae"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hagensia_peringueyi" & AntWeb_database_curated$subspecies == "saldanhae"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hagensia_peringueyi" & AntWeb_database_curated$subspecies == "saldanhae"] <- "Hagensia_peringueyi_saldanhae"
Phylogeny_sample_data_802$Subspecies[Phylogeny_sample_data_802$AntWeb_name == "Hagensia_peringueyi_saldanhae"] <- T


## Cases with wrong name in AntWeb => Taxonomy update from morphotaxon to valid name. Need correction for all specimens in the curated AntWeb database to avoid duplicates and allow merging with GABI.
# Specimen CASENT0782827 is Parvaponera_casc-mz01 in AntWeb. Labeled as Parvaponera_suspecta in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Parvaponera_casc-mz01"] <- "suspecta"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Parvaponera_casc-mz01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Parvaponera_casc-mz01"] <- "Parvaponera_suspecta"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Cryptopone_casc-mz01"] <- "Parvaponera_suspecta"
# Specimen CASENT0817421 is Anochetus_afrc-tz01 in AntWeb. Labeled as Anochetus_talpa in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Anochetus_afrc-tz01"] <- "talpa"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Anochetus_afrc-tz01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Anochetus_afrc-tz01"] <- "Anochetus_talpa"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_afrc-tz01"] <- "Anochetus_talpa"
# Specimen CASENT0888216 is Bothroponera_afrc-mz02 in AntWeb. Labeled as Bothroponera_silvestrii in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Bothroponera_afrc-mz02"] <- "silvestrii"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Bothroponera_afrc-mz02"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Bothroponera_afrc-mz02"] <- "Bothroponera_silvestrii"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Bothroponera_afrc-mz02"] <- "Bothroponera_silvestrii"
# Specimen CASENT0650372 is Brachyponera_mad01 in AntWeb. Labeled as Brachyponera_croceicornis in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Brachyponera_mad01"] <- "croceicornis"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Brachyponera_mad01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Brachyponera_mad01"] <- "Brachyponera_croceicornis"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Brachyponera_mad01"] <- "Brachyponera_croceicornis"
# Specimen CASENT0650371 is Brachyponera_weam01 in AntWeb. Labeled as Brachyponera_lutea in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Brachyponera_weam01"] <- "lutea"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Brachyponera_weam01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Brachyponera_weam01"] <- "Brachyponera_lutea"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Brachyponera_weam01"] <- "Brachyponera_lutea"
# Specimen CASENT0650282 is Hypoponera_jtl027 in AntWeb. Labeled as Hypoponera_aliena in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hypoponera_jtl027"] <- "aliena"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hypoponera_jtl027"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hypoponera_jtl027"] <- "Hypoponera_aliena"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_jtl027"] <- "Hypoponera_aliena"
# Specimen CASENT0783280 is Hypoponera_casc-mz26 in AntWeb. Labeled as Hypoponera_dulcis in the phylogeny. Update all specimens in AntWeb data, including Hypoponera_ao03, and Hypoponera_ug01 according to Brian's comments
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz26", "Hypoponera_ao03", "Hypoponera_ug01")] <- "dulcis"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz26", "Hypoponera_ao03", "Hypoponera_ug01")] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz26", "Hypoponera_ao03", "Hypoponera_ug01")] <- "Hypoponera_dulcis"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_casc-mz26"] <- "Hypoponera_dulcis"
# Specimen CASENT0785139 is Hypoponera_casc-mz31 in AntWeb. Labeled as Hypoponera_spei in the phylogeny. Update all specimens in AntWeb data, including Hypoponera_afrc-za04 according to Brian's comments
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz31", "Hypoponera_afrc-za04")] <- "spei"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz31", "Hypoponera_afrc-za04")] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_casc-mz31", "Hypoponera_afrc-za04")] <- "Hypoponera_spei"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_casc-mz31"] <- "Hypoponera_spei"
# Specimen CASENT0361652 is Hypoponera_ug03 in AntWeb. Labeled as Hypoponera_tristis in the phylogeny. Update all specimens in AntWeb data, including Hypoponera_ug04, Hypoponera_ug05, Hypoponera_ug06 according to Brian's comments
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug03", "Hypoponera_ug04", "Hypoponera_ug05", "Hypoponera_ug06")] <- "tristis"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug03", "Hypoponera_ug04", "Hypoponera_ug05", "Hypoponera_ug06")] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug03", "Hypoponera_ug04", "Hypoponera_ug05", "Hypoponera_ug06")] <- "Hypoponera_tristis"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ug03"] <- "Hypoponera_tristis"
# Specimen CASENT0774569 is Leptogenys_castanea in AntWeb as labeled in the phylogeny. Not Leptogenys_casc-mz04. Update all other Leptogenys_casc-mz04 specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Leptogenys_casc-mz04"] <- "castanea"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Leptogenys_casc-mz04"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Leptogenys_casc-mz04"] <- "Leptogenys_castanea"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_casc-mz04"] <- "Leptogenys_castanea"
# Specimen CASENT0254413 is Mesoponera_afrc-zm01 in AntWeb. Labeled as Mesoponera_ambigua in the phylogeny. Update all specimens in AntWeb data, including Mesoponera_afrc-ug01 according to Phil's comments
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species %in% c("Mesoponera_afrc-zm01", "Mesoponera_afrc-ug01")] <- "ambigua"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species %in% c("Mesoponera_afrc-zm01", "Mesoponera_afrc-ug01")] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species %in% c("Mesoponera_afrc-zm01", "Mesoponera_afrc-ug01")] <- "Mesoponera_ambigua"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Mesoponera_afrc-zm01"] <- "Mesoponera_ambigua"
# Specimen CASENT0267032 is Odontoponera_ph01 in AntWeb. Labeled as Odontoponera_denticulata in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Odontoponera_ph01"] <- "denticulata"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Odontoponera_ph01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Odontoponera_ph01"] <- "Odontoponera_denticulata"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Odontoponera_ph01"] <- "Odontoponera_denticulata"
# Specimen CASENT0923419 is Odontoponera_np01 in AntWeb. Labeled as Odontoponera_transversa in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Odontoponera_np01"] <- "transversa"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Odontoponera_np01"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Odontoponera_np01"] <- "Odontoponera_transversa"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Odontoponera_np01"] <- "Odontoponera_transversa"
# Specimen CASENT0354210 is Anochetus_afr05 in AntWeb as labeled in the phylogeny. Not Anochetus_ug01. Update all other Anochetus_ug01 specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Anochetus_ug01"] <- "afr05"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Anochetus_ug01"] <- "Anochetus_afr05"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_ug01"] <- "Anochetus_afr05"
# Specimen CASENT0820826 is Anochetus_punctaticeps_cf in AntWeb. Labeled as another morphotaxon Anochetus_katonae_nr1 in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Anochetus_punctaticeps_cf"] <- "katonae_nr1"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Anochetus_punctaticeps_cf"] <- "Anochetus_katonae_nr1"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Anochetus_punctaticeps_cf"] <- "Anochetus_katonae_nr1"
# Specimen CASENT0374097 is Cryptopone_th01 in AntWeb. Labeled as another morphotaxon Ectomomyrmex_th01 in the phylogeny, but would be a synonym of another terminal so update to Ectomomyrmex_th05. Update this unique specimen in AntWeb data.
AntWeb_database_curated$genus[AntWeb_database_curated$code == "casent0374097"] <- "Ectomomyrmex"
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0374097"] <- "th05"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0374097"] <- "Ectomomyrmex_th05"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0374097"] <- "Ectomomyrmex_th05"
# Specimen CASENT0784819 is Hypoponera_casc-mz06 in AntWeb. Labeled as Hypoponera_angustata in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hypoponera_casc-mz06"] <- "angustata"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hypoponera_casc-mz06"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hypoponera_casc-mz06"] <- "Hypoponera_angustata"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_casc-mz06"] <- "Hypoponera_angustata"
# Specimen CASENT0356319 is Hypoponera_ug12 in AntWeb. Labeled as Hypoponera_blanda in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hypoponera_ug12"] <- "blanda"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hypoponera_ug12"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hypoponera_ug12"] <- "Hypoponera_blanda"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ug12"] <- "Hypoponera_blanda"
# Specimen CASENT0356532 is Hypoponera_ug07 in AntWeb. Labeled as Hypoponera_fatiga in the phylogeny. Update all specimens in AntWeb data, including Hypoponera_ug08 according to Brian's comments
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug07", "Hypoponera_ug08")] <- "fatiga"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug07", "Hypoponera_ug08")] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species %in% c("Hypoponera_ug07", "Hypoponera_ug08")] <- "Hypoponera_fatiga"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ug07"] <- "Hypoponera_fatiga"
# Specimen CASENT0128290 is Hypoponera_ao02 in AntWeb. Labeled as Hypoponera_jeanneli in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "jeanneli"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "Hypoponera_jeanneli"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ao02"] <- "Hypoponera_jeanneli"
# Specimen CASENT0128290 is Hypoponera_ao02 in AntWeb. Labeled as Hypoponera_jeanneli in the phylogeny. Update this unique specimen in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "jeanneli"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Hypoponera_ao02"] <- "Hypoponera_jeanneli"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hypoponera_ao02"] <- "Hypoponera_jeanneli"
# Specimen CASENT0644557 is Rasopone_jtl049 in AntWeb. Labeled as Rasopone_lunaris in the phylogeny. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Rasopone_jtl049"] <- "lunaris"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Rasopone_jtl049"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Rasopone_jtl049"] <- "Rasopone_lunaris"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Rasopone_jtl049"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_lunaris"] <- "Rasopone_jtl049"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_jtl049"] <- "Rasopone_lunaris"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$AntWeb_name == "Rasopone_jtl049"] <- T
# Specimen CASENT0649093 is Thaumatomyrmex_ferox_complex in AntWeb. Labeled as Thaumatomyrmex_ferox in the phylogeny. We disregard the complex issue and label all specimens as a unique valid species. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Thaumatomyrmex_ferox_complex"] <- "ferox"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Thaumatomyrmex_ferox_complex"] <- T
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Thaumatomyrmex_ferox_complex"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Thaumatomyrmex_ferox_complex"] <- "Thaumatomyrmex_ferox"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Thaumatomyrmex_ferox"] <- "Thaumatomyrmex_ferox_complex"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$AntWeb_name == "Thaumatomyrmex_ferox_complex"] <- "Thaumatomyrmex_ferox"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$AntWeb_name == "Thaumatomyrmex_ferox_complex"] <- T
# Specimen CASENT0628711 is Leptogenys_pubiceps_complex in AntWeb. Labeled as Leptogenys_pubiceps in the phylogeny. We disregard the complex issue and label all specimens as a unique valid species. Update all specimens in AntWeb data.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Leptogenys_pubiceps_complex"] <- "pubiceps"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Leptogenys_pubiceps_complex"] <- T
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Leptogenys_pubiceps_complex"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Leptogenys_pubiceps_complex"] <- "Leptogenys_pubiceps"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_pubiceps"] <- "Leptogenys_pubiceps_complex"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_pubiceps_complex"] <- "Leptogenys_pubiceps"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$AntWeb_name == "Leptogenys_pubiceps_complex"] <- T
# Specimen ZRC_ENT00007921 is Ectomomyrmex_overbecki_cf in AntWeb. Labeled as Ectomomyrmex_overbecki in the phylogeny. Update this unique specimen in AntWeb.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Ectomomyrmex_overbecki_cf"] <- "overbecki"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Ectomomyrmex_overbecki_cf"] <- T
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Ectomomyrmex_overbecki_cf"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Ectomomyrmex_overbecki_cf"] <- "Ectomomyrmex_overbecki"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Ectomomyrmex_overbecki"] <- "Ectomomyrmex_overbecki_cf"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$AntWeb_name == "Ectomomyrmex_overbecki_cf"] <- "Ectomomyrmex_overbecki"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$AntWeb_name == "Ectomomyrmex_overbecki_cf"] <- T
# Specimen CASENT0650188 is Odontomachus_opaculus_c in AntWeb. Labeled as Odontomachus_imperator in the phylogeny. Update this unique specimen in AntWeb.
AntWeb_database_curated$species[AntWeb_database_curated$Genus_species == "Odontomachus_opaculus_c"] <- "imperator"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Odontomachus_opaculus_c"] <- T
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$Genus_species == "Odontomachus_opaculus_c"] <- "valid"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$Genus_species == "Odontomachus_opaculus_c"] <- "Odontomachus_imperator"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Odontomachus_imperator"] <- "Odontomachus_opaculus_c"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$AntWeb_name == "Odontomachus_opaculus_c"] <- "Odontomachus_imperator"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$AntWeb_name == "Odontomachus_opaculus_c"] <- T


## Cases with wrong name in AntWeb => Taxonomy update from valid name to another valid name. Need correction in the curated database to avoid duplicates and allow merging with GABI. May change only the given specimen if other specimens are not dubious.
# Specimen CASENT0816349 is Bothroponera_ilgii in AntWeb. Labeled as Bothroponera_ancilla in the phylogeny. Update only this specimen AntWeb as Bothroponera_ilgii has a larger distribution across Africa. 
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0816349"] <- "ancilla"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0816349"] <- "Bothroponera_ancilla"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Bothroponera_ilgii"] <- "Bothroponera_ancilla"
# Specimen CASENT0803835 is Mesoponera_ambigua in AntWeb. Labeled as Mesoponera_elisae in the phylogeny because considered different from Mainland's Mesoponera_ambigua. Update this unique specimen in AntWeb data since not sure this change applies to all Mesoponera_ambigua from Madagascar. Note: some specimens of Mesoponera_elisae are found in Africa's mainland too.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0803835"] <- "elisae"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0803835"] <- "Mesoponera_elisae"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0803835"] <- "Mesoponera_elisae"
# Specimen CASENT0816544 is Mesoponera_ingesta in AntWeb. Labeled as Mesoponera_caffaria in the phylogeny. Both species have large scale distribution, so in the absence of rationale to justify this change, only update this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0816544"] <- "caffaria"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0816544"] <- "Mesoponera_caffaria"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816544"] <- "Mesoponera_caffaria"
# Specimen CASENT0631798 is Neoponera_moesta in AntWeb. Labeled as Neoponera_crenata in the phylogeny. Both species have large scale distribution across Neotropics, so in the absence of rationale to justify this change, only update this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0631798"] <- "crenata"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0631798"] <- "Neoponera_crenata"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0631798"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0631798"] <- "Neoponera_moesta"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0631798"] <- T
# Specimen CASENT0649900 is Neoponera_crenata in AntWeb. Labeled as Neoponera_moesta in the phylogeny. Both species have large scale distribution across Neotropics, so in the absence of rationale to justify this change, only update this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0649900"] <- "moesta"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0649900"] <- "Neoponera_moesta"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0649900"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0649900"] <- "Neoponera_crenata"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0649900"] <- T
# Specimen ZRC_HYM_0000557 is Brachyponera_jerdonii in AntWeb. Labeled as Brachyponera_nigrita in the phylogeny. Both species have large scale distribution across SouthEast Asia, so in the absence of rationale to justify this change, only update this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "zrc_hym_0000557"] <- "nigrita"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "zrc_hym_0000557"] <- "Brachyponera_nigrita"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "zrc_hym_0000557"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "ZRC_HYM_0000557"] <- "Brachyponera_jerdonii"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "ZRC_HYM_0000557"] <- T
# Specimen CASENT0650114 is Leptogenys_peninsularis in AntWeb. Labeled as Leptogenys_sonora in the phylogeny. Other specimens of Leptogenys_peninsularis are found in Baja California (a peninsula) with one outlier in Guerrero state, South of Mexico City. By default, I renamed all five specimens in Sonora as Leptogenys_sonora. Not sure what to do with the outlier.
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Leptogenys_peninsularis") & (AntWeb_database_curated$adm1 == "Sonora")] <- "sonora"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Leptogenys_peninsularis") & (AntWeb_database_curated$adm1 == "Sonora")] <- "Leptogenys_sonora"
AntWeb_database_curated$Name_updated_from_Antweb[(AntWeb_database_curated$Genus_species == "Leptogenys_peninsularis") & (AntWeb_database_curated$adm1 == "Sonora")] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650114"] <- "Leptogenys_peninsularis"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650114"] <- T
# Specimen CASENT0650198 is Ectomomyrmex_scobinus in AntWeb. Labeled as Ectomomyrmex_aciculatus in the phylogeny. Both species are found in PNG. In absence of rationale to justify the change, I only renamed this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0650198"] <- "aciculatus"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0650198"] <- "Ectomomyrmex_aciculatus"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0650198"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650198"] <- "Ectomomyrmex_scobinus"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650198"] <- T
# Specimen CASENT0882084 is Leptogenys_kraepelini in AntWeb. Labeled as Leptogenys_chinensis in the phylogeny. Both species are found scattered in Southeast Asia. In absence of rationale to justify the change, I only renamed this unique specimen.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0882084"] <- "chinensis"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0882084"] <- "Leptogenys_chinensis"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0882084"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0882084"] <- "Leptogenys_kraepelini"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0882084"] <- T


View(AntWeb_database_curated[AntWeb_database_curated$Genus_species == "Hypoponera_ug04", ])
View(AntWeb_database_curated[AntWeb_database_curated$Genus_species == "Parvaponera_suspecta", ])


## Cases with indet. name in AntWeb => Retrieve specimen data (and associated specimens) from complete (non-curated) AntWeb database, correct the name, and add it to the curated database
# Specimen CASENT0649902 is Neoponera (indet.) in AntWeb. Labeled as Neoponera bucki in the phylogeny as it is "probably a male Neoponera bucki" (Jack Longino). Need a new Genus because not close to other Neoponera
Neoponera_indet_data <- AntWeb_database[AntWeb_database$code == "casent0649902", ]
Neoponera_indet_data$Genus_species <- paste0(Neoponera_indet_data$genus, "_", Neoponera_indet_data$species)
Neoponera_indet_data$Status_AntCat <-  Neoponera_indet_data$status
AntWeb_database_curated <- rbind(AntWeb_database_curated, Neoponera_indet_data)
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Neoponera_indet"] <- Neoponera_indet_data$Genus_species
# Specimen CASENT0295217 is Fisheropone (indet) in AntWeb. Labeled as Fisheropone_ambigua but also labeled as new so better keep aside.
Fisheropone_indet_data <- AntWeb_database[AntWeb_database$code == "casent0295217", ]
Fisheropone_indet_data$Genus_species <- paste0(Fisheropone_indet_data$genus, "_", Fisheropone_indet_data$species)
Fisheropone_indet_data$Status_AntCat <- Fisheropone_indet_data$status
AntWeb_database_curated <- rbind(AntWeb_database_curated, Fisheropone_indet_data)
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Fisheropone_indet"] <- Fisheropone_indet_data$Genus_species
# Specimen CASENT0007651 is Hagensia_(indet) in AntWeb. Labeled as Hagensia_indet in the phylogeny.
Hagensia_indet_data <- AntWeb_database[AntWeb_database$code == "casent0007651", ]
Hagensia_indet_data$Genus_species <- paste0(Hagensia_indet_data$genus, "_", Hagensia_indet_data$species)
Hagensia_indet_data$Status_AntCat <- Hagensia_indet_data$status
AntWeb_database_curated <- rbind(AntWeb_database_curated, Hagensia_indet_data)
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$AntWeb_name == "Hagensia_indet"] <- Hagensia_indet_data$Genus_species
# Specimen CASENT0650351 is Odontomachus_(indet) in AntWeb. Labeled as Odontomachus_indet in the phylogeny.
Odontomachus_indet_data <- AntWeb_database[AntWeb_database$code == "casent0650351", ]
Odontomachus_indet_data$Genus_species <- paste0(Odontomachus_indet_data$genus, "_", Odontomachus_indet_data$species)
Odontomachus_indet_data$Status_AntCat <- Odontomachus_indet_data$status
AntWeb_database_curated <- rbind(AntWeb_database_curated, Odontomachus_indet_data)
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650351"] <- Odontomachus_indet_data$Genus_species
# Specimen CASENT0635390 is Ponera_(indet) in AntWeb. Labeled as Ponera_sp in the phylogeny.
Ponera_indet_data <- AntWeb_database[AntWeb_database$code == "casent0635390", ]
Ponera_indet_data$Genus_species <- paste0(Ponera_indet_data$genus, "_", Ponera_indet_data$species)
Ponera_indet_data$Status_AntCat <- Ponera_indet_data$status
AntWeb_database_curated <- rbind(AntWeb_database_curated, Ponera_indet_data)
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0635390"] <- Ponera_indet_data$Genus_species


## Case of specimens with updated names that were valid initially in AntWeb, and are now attributed to another name (valid or not).
# Need to change the specimen name in AntWeb database, and see if this change extend to other specimens

# Look for duplicates that highlight possible taxonomic update
Phylogeny_sample_data_802$AntWeb_name[duplicated(Phylogeny_sample_data_802$AntWeb_name)]

# Cryptopone_gilva # CASENT0749266 # Cryptopone_gilva # Stay the same
# Cryptopone_gilva # CASENT0614525 # Cryptopone_gilvagrande # Valid to another valid name
# Cryptopone_gilva # CASENT0623898 # Cryptopone_guatemalensis # Valid to another valid name
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0614525"] <- "Cryptopone_gilvagrande"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0623898"] <- "Cryptopone_guatemalensis"

# Neoponera_laevigata # ATPFOR2006 # Neoponera_laevigata # Stay the same
# Neoponera_laevigata # CASENT0637827 # Neoponera_mashpi # Valid to another valid name
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0637827"] <- "Neoponera_mashpi"

# Anochetus_siphneus # CASENT0815716 # Anochetus_siphneus # Stay the same
# Anochetus_siphneus # CASENT0066876 # Anochetus_katonae_nr2 # Valid to morphotaxon # Need change in AntWeb database.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0066876"] <- "katonae_nr2"
AntWeb_database_curated$status[AntWeb_database_curated$code == "casent0066876"] <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$code == "casent0066876"] <- "morphotaxon"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0066876"] <- "Anochetus_katonae_nr2"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0066876"] <- "Anochetus_katonae_nr2"
# Also need to change all specimens records from Zambia!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "katonae_nr2"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Zambia")] <- "Anochetus_katonae_nr2"
# Rename specimen from DRC for a new unique name (Anochetus_katonae_nr3) not to use in our database!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "katonae_nr3"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Anochetus_siphneus") & (AntWeb_database_curated$country == "Democratic Republic of Congo")] <- "Anochetus_katonae_nr3"

# Anochetus_rugosus # CASENT0635363 # Anochetus_rugosus # Stay the same
# Anochetus_MY07 # CASENT0387827 # Anochetus_rugosus # True duplicate. May need to be removed from phylogeny

# Ectomomyrmex_simillimus # CASENT0923418 # Ectomomyrmex_simillimus # Stay the same
# Ectomomyrmex_Janda_sp10 # CASENT0650200 # Ectomomyrmex_simillimus # Wrongly labeled in the phylogeny as Ectomomyrmex_Janda_sp10. # True duplicate. May need to be removed from phylogeny.

# Hypoponera_opacior_nr # CASENT0886661 # Hypoponera_opacior_nr1 # Morphotaxon to morphotaxon. Need change in AntWeb database.
# Hypoponera_us_ca01 # CASENT0886861 # Hypoponera_opacior_nr2 # Morphotaxon to morphotaxon. Need change in AntWeb database.
# Hypoponera_opacior_nr # CASENT0886662 # Hypoponera_opacior_nr3 # Morphotaxon to morphotaxon
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0886661"] <- "opacior_nr1"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0886661"] <- "Hypoponera_opacior_nr1"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886661"] <- "Hypoponera_opacior_nr1"
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0886861"] <- "opacior_nr2"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0886861"] <- "Hypoponera_opacior_nr2"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886861"] <- "Hypoponera_opacior_nr2"
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0886662"] <- "opacior_nr3"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0886662"] <- "Hypoponera_opacior_nr3"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886662"] <- "Hypoponera_opacior_nr3"

# Hypoponera_opacior_nr # CASENT0882410 # Hypoponera_opacior_nr1 # Morphotaxon to morphotaxon. Not in the phylogeny but used for trait measurements. Need change in AntWeb database.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0882410"] <- "opacior_nr1"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0882410"] <- "Hypoponera_opacior_nr1"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0882410"] <- T

# Hypoponera_trigona # CASENT0635113 # Hypoponera_trigona # Stay the same
# Hypoponera_trigona # CASENT0650231 # Hypoponera_JTL038 # Valid to morphotaxon
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650231"] <- "Hypoponera_jtl038"

# Odontomachus_brunneus # CASENT0649761 # Odontomachus_brunneus # Stay the same
# Odontomachus_brunneus # CASENT0649078 # Odontomachus_brunneus_nr # Valid to morphotaxon. Need change in AntWeb database.
AntWeb_database_curated$species[AntWeb_database_curated$code == "casent0649078"] <- "brunneus_nr"
AntWeb_database_curated$status[AntWeb_database_curated$code == "casent0649078"] <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$code == "casent0649078"] <- "morphotaxon"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "casent0649078"] <- "Odontomachus_brunneus_nr"
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0649078"] <- "Odontomachus_brunneus_nr"
# Also need to change all specimens records from South America/Outside South East USA!!!
AntWeb_database_curated$species[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "brunneus_nr"
AntWeb_database_curated$status[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")]  <- "morphotaxon"
AntWeb_database_curated$Status_AntCat[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "morphotaxon"
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "Odontomachus_brunneus_nr"

# Ponera_incerta # CASENT0882069 # Ponera_incerta # Stay the same
# Ponera_incerta # CASENT0637782 # Ponera_petila # Valid to another valid name
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0637782"] <- "Ponera_petila"

# Specimen DZUP549431 is Neoponera_bucki in Antweb and as labeled in the phylogeny. Not Relict_Neoponera. 
# However, not even a Neoponera and need its own new genus. # Named as RelictNeoponera_sp for the moment.
# Not sure if all specimens of Neoponera bucki in Antweb are in this case, but since they are all localized in the same area in Brazil, it is not an issue for Biogeography.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "DZUP549431"] <- "RelictNeoponera_(indet)"
AntWeb_database_curated$genus[AntWeb_database_curated$code == "dzup549431"] <- "RelictNeoponera"
AntWeb_database_curated$species[AntWeb_database_curated$code == "dzup549431"] <- "(indet)"
AntWeb_database_curated$status[AntWeb_database_curated$code == "dzup549431"] <- "indetermined"
AntWeb_database_curated$Status_AntCat[AntWeb_database_curated$code == "dzup549431"] <- "indetermined"
AntWeb_database_curated$Genus_species[AntWeb_database_curated$code == "dzup549431"] <- "RelictNeoponera_(indet)"

# Specimen CASENT0650201 is Mesoponera_papuana as labeled in the phylogeny. Not Mesoponera_manni.
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650201"] <- "Mesoponera_papuana"


# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")


### 2.3.3/ Flag for presence or not in the phylogeny ####

# In the taxa-level summary table
Phylogeny_sample_data_802$In_phylogeny <- T

# In the curated AntWeb database
AntWeb_database_curated$In_phylogeny <- F
AntWeb_database_curated$In_phylogeny[(AntWeb_database_curated$Genus_species %in% Phylogeny_sample_data_802$AntWeb_name)] <- T
# AntWeb_database_curated$In_phylogeny[(AntWeb_database_curated$Genus_species %in% Phylogeny_sample_data_802$Current_name)] <- T

### 2.3.4/ Flag for specimens where name from AntWeb is "wrong" (different from the one we are using) ####

# Need to add a current name and an AntWeb name field in both specimen and taxa-level datasets (which I should have kept in the first place...)

# In the taxa-level summary table
Phylogeny_sample_data_802$Current_name <- Phylogeny_sample_data_802$AntWeb_name
Phylogeny_sample_data_802$Name_updated_from_Antweb <- F
Phylogeny_sample_data_802$Name_changed_from_phylogeny <- NA

# In the curated AntWeb database
AntWeb_database_curated$AntWeb_name <- AntWeb_database_curated$Genus_species
AntWeb_database_curated$Name_updated_from_Antweb <- F

## Reupdate AntWeb name with their initial name and flag them for this changes

## Updated species names (all specimens with this name)
# CASENT0782827 # Parvaponera_casc-mz01 # "Parvaponera_suspecta"
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Parvaponera_casc-mz01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Parvaponera_casc-mz01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0782827"] <- "Parvaponera_casc-mz01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0782827"] <- T
# CASENT0817421 # Anochetus_afrc-tz01 # Anochetus_talpa
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Anochetus_afrc-tz01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Anochetus_afrc-tz01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0817421"] <- "Anochetus_afrc-tz01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0817421"] <- T
# CASENT0888216 # Bothroponera_afrc-mz02 # Bothroponera_silvestrii
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Bothroponera_afrc-mz02"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Bothroponera_afrc-mz02"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0888216"] <- "Bothroponera_afrc-mz02"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0888216"] <- T
# CASENT0650372 # Brachyponera_mad01 # Brachyponera_croceicornis
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Brachyponera_mad01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Brachyponera_mad01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650372"] <- "Brachyponera_mad01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650372"] <- T
# CASENT0650371 # Brachyponera_weam01 # Brachyponera_lutea
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Brachyponera_weam01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Brachyponera_weam01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650371"] <- "Brachyponera_weam01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650371"] <- T
# CASENT0650282 # Hypoponera_jtl027 # Hypoponera_aliena
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_jtl027"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_jtl027"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650282"] <- "Hypoponera_jtl027"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0650282"] <- T
# CASENT0783280 # Hypoponera_casc-mz26, Hypoponera_ao03, and Hypoponera_ug01 # Hypoponera_dulcis
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_casc-mz26"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_casc-mz26"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ao03"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ao03"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0783280"] <- "Hypoponera_casc-mz26"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0783280"] <- T
# CASENT0785139 # Hypoponera_casc-mz31 # Hypoponera_spei
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_casc-mz31"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_casc-mz31"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0785139"] <- "Hypoponera_casc-mz31"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0785139"] <- T
# CASENT0361652 # Hypoponera_ug03 # Hypoponera_tristis
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug03"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug03"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug04"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug04"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug05"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug05"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug06"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug06"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0361652"] <- "Hypoponera_ug03"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0361652"] <- T
# CASENT0774569 # Leptogenys_casc-mz04 # Leptogenys_castanea # No AntWeb update for the specimen with molecular data
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Leptogenys_casc-mz04"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Leptogenys_casc-mz04"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
# CASENT0254413 # Mesoponera_afrc-zm01, Mesoponera_afrc-ug01 # Mesoponera_ambigua
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Mesoponera_afrc-zm01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Mesoponera_afrc-zm01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Mesoponera_afrc-ug01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Mesoponera_afrc-ug01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0254413"] <- "Mesoponera_afrc-zm01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0254413"] <- T
# CASENT0267032 # Odontoponera_ph01 # Odontoponera_denticulata
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Odontoponera_ph01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Odontoponera_ph01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0267032"] <- "Odontoponera_ph01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0267032"] <- T
# CASENT0923419 # Odontoponera_np01 # Odontoponera_transversa
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Odontoponera_np01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Odontoponera_np01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0923419"] <- "Odontoponera_np01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0923419"] <- T
# CASENT0354210 # Anochetus_ug01 # Anochetus_afr05 # No AntWeb update for the specimen with molecular data 
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Anochetus_ug01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Anochetus_ug01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
# CASENT0820826 # Anochetus_punctaticeps_cf # Anochetus_katonae_nr1
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Anochetus_punctaticeps_cf"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Anochetus_punctaticeps_cf"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0820826"] <- "Anochetus_punctaticeps_cf"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0820826"] <- T
# CASENT0374097 # Cryptopone_th01 # Ectomomyrmex_th05
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Cryptopone_th01"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Cryptopone_th01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0374097"] <- "Cryptopone_th01"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0374097"] <- T
# CASENT0784819 # Hypoponera_casc-mz06 # Hypoponera_angustata
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_casc-mz06"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_casc-mz06"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0784819"] <- "Hypoponera_casc-mz06"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0784819"] <- T
# CASENT0356319 # Hypoponera_ug12 # Hypoponera_blanda
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug12"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug12"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0356319"] <- "Hypoponera_ug12"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0356319"] <- T
# CASENT0356532 # Hypoponera_ug07, Hypoponera_ug08 # Hypoponera_fatiga
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug07"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug07"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ug08"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ug08"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0356532"] <- "Hypoponera_ug07"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0356532"] <- T
# CASENT0128290 # Hypoponera_ao02 # Hypoponera_jeanneli
specimen_codes <- AntWeb_database$code[paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Hypoponera_ao02"]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Hypoponera_ao02"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0128290"] <- "Hypoponera_ao02"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0128290"] <- T

## Update specimens (single specimen)

# CASENT0816349 # Bothroponera_ilgii # Bothroponera_ancilla
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0816349"] <- "Bothroponera_ilgii"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0816349"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816349"] <- "Bothroponera_ilgii"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816349"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816349"] <- T
# CASENT0803835 # Mesoponera_ambigua # Mesoponera_elisae
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0803835"] <- "Mesoponera_ambigua"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0803835"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0803835"] <- "Mesoponera_ambigua"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0803835"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0803835"] <- T
# CASENT0816544 # Mesoponera_ingesta # Mesoponera_caffaria
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0816544"] <- "Mesoponera_ingesta"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0816544"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816544"] <- "Mesoponera_ingesta"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816544"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0816544"] <- T
# CASENT0886661 # Hypoponera_opacior_nr # Hypoponera_opacior_nr1
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0886661"] <- "Hypoponera_opacior_nr"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0886661"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886661"] <- "Hypoponera_opacior_nr"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886661"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886661"] <- T
# CASENT0886861 # Hypoponera_us_ca01 # Hypoponera_opacior_nr2
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0886861"] <- "Hypoponera_us_ca01"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0886861"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886861"] <- "Hypoponera_us_ca01"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886861"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886861"] <- T
# CASENT0886662 # Hypoponera_opacior_nr # Hypoponera_opacior_nr3
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "casent0886662"] <- "Hypoponera_opacior_nr"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "casent0886662"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886662"] <- "Hypoponera_opacior_nr"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886662"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0886662"] <- T
# DZUP549431 # Neoponera_bucki # RelictNeoponera_(indet)
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code == "dzup549431"] <- "Neoponera_bucki"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code == "dzup549431"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "DZUP549431"] <- "Neoponera_bucki"
Phylogeny_sample_data_802$Current_name[Phylogeny_sample_data_802$`Specimen code` == "DZUP549431"]
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "DZUP549431"] <- T

## Special cases
# CASENT0066876 # Anochetus_siphneus # Anochetus_katonae_nr2 # All specimens from Zambia
specimen_codes <- AntWeb_database$code[(paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Anochetus_siphneus") & (AntWeb_database$country == "Zambia")]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Anochetus_siphneus"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0066876"] <- "Anochetus_siphneus"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0066876"] <- T
# Rename specimen from DRC for a new unique name (Anochetus_katonae_nr3) not to use in our database!!!
specimen_codes <- AntWeb_database$code[(paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Anochetus_siphneus") & (AntWeb_database$country == "Democratic Republic of Congo")]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Anochetus_siphneus"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
# CASENT0649078 # Odontomachus_brunneus # Odontomachus_brunneus_nr # All specimens records from South America/Outside South East USA.
specimen_codes <- AntWeb_database$code[(paste0(AntWeb_database$genus, "_", AntWeb_database$species) == "Odontomachus_brunneus") & (AntWeb_database$biogeographicregion == "Neotropical")]
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$code %in% specimen_codes] <- "Odontomachus_brunneus"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$code %in% specimen_codes] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$`Specimen code` == "CASENT0649078"] <- "Odontomachus_brunneus"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$`Specimen code` == "CASENT0649078"] <- T
AntWeb_database_curated$Genus_species[(AntWeb_database_curated$Genus_species == "Odontomachus_brunneus") & (AntWeb_database_curated$biogeographicregion == "Neotropical")] <- "Odontomachus_brunneus_nr"

## Updated subspecies
# "Hagensia_havilandi" # "Hagensia_havilandi_marleyi"
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$Genus_species == "Hagensia_havilandi_marleyi"] <- "Hagensia_havilandi"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Hagensia_havilandi_marleyi"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$Current_name == "Hagensia_havilandi_marleyi"] <- "Hagensia_havilandi"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$Current_name == "Hagensia_havilandi_marleyi"] <- T
# "Parvaponera_darwinii" # "Parvaponera_darwinii_madecassa"
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$Genus_species == "Parvaponera_darwinii_madecassa"] <- "Parvaponera_darwinii"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Parvaponera_darwinii_madecassa"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$Current_name == "Parvaponera_darwinii_madecassa"] <- "Parvaponera_darwinii"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$Current_name == "Parvaponera_darwinii_madecassa"] <- T
# "Hagensia_peringueyi" # "Hagensia_peringueyi_saldanhae"
AntWeb_database_curated$AntWeb_name[AntWeb_database_curated$Genus_species == "Hagensia_peringueyi_saldanhae"] <- "Hagensia_peringueyi"
AntWeb_database_curated$Name_updated_from_Antweb[AntWeb_database_curated$Genus_species == "Hagensia_peringueyi_saldanhae"] <- T
Phylogeny_sample_data_802$AntWeb_name[Phylogeny_sample_data_802$Current_name == "Hagensia_peringueyi_saldanhae"] <- "Hagensia_peringueyi"
Phylogeny_sample_data_802$Name_updated_from_Antweb[Phylogeny_sample_data_802$Current_name == "Hagensia_peringueyi_saldanhae"] <- T

# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")


## Adjust "change name" fields by comparing current names with AntWeb and Phylo

# Remove Name_update
Phylogeny_sample_data_802 <- Phylogeny_sample_data_802 %>% 
  select(-Name_update)

# Check comparison of current names with AntWeb names in Summary taxa-level df
Updated_AntWeb_names_match <- !(Phylogeny_sample_data_802$Current_name == Phylogeny_sample_data_802$AntWeb_name)
View(Phylogeny_sample_data_802[!(Updated_AntWeb_names_match == Phylogeny_sample_data_802$Name_updated_from_Antweb), ])

# Check comparison of current names with AntWeb names in currated AntWeb database
Updated_AntWeb_names_match <- !(AntWeb_database_curated$Genus_species == AntWeb_database_curated$AntWeb_name)
View(AntWeb_database_curated[!(Updated_AntWeb_names_match == AntWeb_database_curated$Name_updated_from_Antweb), ])

# Record change between current names and phylogey names
Phylogeny_sample_data_802$Name_changed_from_phylogeny <- !(Phylogeny_sample_data_802$Current_name == Phylogeny_sample_data_802$Phylo_name)


### 2.3.5/ Flag true duplicates to eventually remove ####

# Look for duplicates that highlight possible taxonomic update
Phylogeny_sample_data_802$Current_name[duplicated(Phylogeny_sample_data_802$Current_name)]

Phylogeny_sample_data_802$Current_duplicate <- F
Phylogeny_sample_data_802$Current_duplicate[duplicated(Phylogeny_sample_data_802$Current_name)] <- T

# Add comments on duplicate to justify selection based on UCE results

# Anochetus_rugosus: CASENT0635363 vs. CASENT0387827
  # CASENT0635363 should be removed.
  # CASENT0387827 has higher number of contigs, median lenght and max length. Should be kept.
  
# Ectomomyrmex_simillimus: CASENT0650200 vs. CASENT0923418
  # CASENT0650200 has higher number of contigs, median lenght and max length. Should be kept.
  # CASENT0923418 should be removed.


### 2.3.6/ Add name status in metadata table for phylogeny data: valid, morphotaxon, indeterminate ####

# Add new field to record name status
Phylogeny_sample_data_802$Name_status <- NA

# Extract name status of specimen entries from taxa in the phylogeny
AntWeb_database_names_curated_phylo_only <- AntWeb_database_curated[AntWeb_database_curated$In_phylogeny, c("Genus_species", "Status_AntCat")] %>% 
  distinct() 

Phylogeny_sample_data_802$Name_status <- AntWeb_database_names_curated_phylo_only$Status_AntCat[match(Phylogeny_sample_data_802$Current_name, AntWeb_database_names_curated_phylo_only$Genus_species)]

# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")

# Save curated AntWeb database
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")

## 2.3.7/ Add comments ####

## Merge previous comments from ponerinae-dataset-sample-list-current

Former_summary_df <- read_excel("./input_data/Phylogenies/ponerinae-dataset-sample-list-current.xlsx")
Former_comments <- Former_summary_df %>% 
  rename(Notes_2022_04_26 = `Notes - 2022-04-26`) %>% 
  rename(Specimen_code = `Specimen code`) %>%
  filter(Former_summary_df$`Ponerinae-dataset` == "KEEP") %>% 
  select(Specimen_code, Notes_2022_04_26)
  
Phylogeny_sample_data_802 <- left_join(x = Phylogeny_sample_data_802, y = Former_comments, by = c("Specimen code" = "Specimen_code"))

## Add manual comments on the name changes

# Export as Excel file. Write comments. Reimport. Save as .rds
# xlsx::write.xlsx(x = Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.xlsx", showNA = F, row.names = F)
openxlsx::write.xlsx(x = Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.xlsx", overwrite = T)
Phylogeny_sample_data_802 <- read_excel("./input_data/Phylogenies/Phylogeny_sample_data_802.xlsx")

# Save updated df for Metadata of samples in the phylogeny
saveRDS(Phylogeny_sample_data_802, file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")

## 2.3.8/ Export curated AntWeb database with different levels of inclusion ####

### Remove morphospecies outside of the phylogeny

# Morphospecies may be useful to compute sampling fraction across bioregions, even if we may use only the ones in the phylogeny.
# Valid species without molecular data may be useful to compute sampling fraction across bioregions.

# Save curated AntWeb database with all morphospecies
saveRDS(AntWeb_database_curated, file = "./input_data/AntWeb_data/AntWeb_database_curated.rds")
nrow(AntWeb_database_curated)
# 62,292 specimen data

AntWeb_database_curated_without_morphospecies_no_molecular <- AntWeb_database_curated[(AntWeb_database_curated$In_phylogeny | AntWeb_database_curated$Status_AntCat == "valid"), ]
nrow(AntWeb_database_curated_without_morphospecies_no_molecular)
# 51,605 specimen data

# Save curated AntWeb database without morphospecies without molecular data
saveRDS(AntWeb_database_curated_without_morphospecies_no_molecular, file = "./input_data/AntWeb_data/AntWeb_database_curated_without_morphospecies_no_molecular.rds")

### Keep only taxa in the phylogeny

AntWeb_database_curated_only_phylo <- AntWeb_database_curated[(AntWeb_database_curated$In_phylogeny), ]
nrow(AntWeb_database_curated_only_phylo)
# 43,665 specimen data

## Save curated AntWeb database with only taxa in the phylogeny
saveRDS(AntWeb_database_curated_only_phylo, file = "./input_data/AntWeb_data/AntWeb_database_curated_only_phylo.rds")


### 2.3.9/ Go back to GABI database to flag names using the current curated names for terminals in the phylogeny ####

### Add fields for Names status (Status_AntCat), presence in phylogeny (In_phylogeny), and Name update (Name_updated), 

# All GABI entries have valid AntCat species names
GABI_database_Ponerinae$Status_AntCat <- "valid"

unique(GABI_database_Ponerinae$AntCat_Genus_species_name[!(GABI_database_Ponerinae$AntCat_Genus_species_name %in% Phylogeny_sample_data_802$Current_name)])

# Flag for match with a taxa in the phylogeny
GABI_database_Ponerinae$In_phylogeny <- F
GABI_database_Ponerinae$In_phylogeny[(GABI_database_Ponerinae$AntCat_Genus_species_name %in% Phylogeny_sample_data_802$Current_name)] <- T

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
  select(GABI_accession_ID, Accession_ID, Current_name, Initial_name, Name_updated, Status_AntCat, In_phylogeny, Latitude_dec, Longitude_dec, Country, adm1, adm2, Locality, bentity2_name, Source)
  
AntWeb_database_for_merging <- AntWeb_database_curated %>% 
  rename(Specimen_code = code) %>% 
  rename(Country = country) %>% 
  rename(Locality = localityname) %>%
  rename(Locality_code = localitycode) %>%
  rename(Locality_notes = localitynotes) %>% 
  rename(Latitude_dec = decimal_latitude) %>%  
  rename(Longitude_dec = decimal_longitude) %>% 
  rename(Elevation = elevation) %>% 
  rename(Initial_name = AntWeb_name) %>% 
  rename(Current_name = Genus_species) %>% 
  rename(Name_updated = Name_updated_from_Antweb) %>% 
  select(Specimen_code, Current_name, Initial_name, Name_updated, Status_AntCat, In_phylogeny, Latitude_dec, Longitude_dec, Country, adm1, adm2, Locality, Locality_code, Locality_notes, Elevation, Source)

# Merge databases
Biogeographic_database_Ponerinae <- full_join(GABI_database_Ponerinae_for_merging, AntWeb_database_for_merging)

nrow(GABI_database_Ponerinae_for_merging) + nrow(AntWeb_database_for_merging)
dim(Biogeographic_database_Ponerinae)
# 155,860 geolocalized specimen records before biogeographic curation

# Flag duplicates based on taxa name and geographic information 
Biogeographic_database_Ponerinae$duplicate <- duplicated(Biogeographic_database_Ponerinae[, c("Current_name", "Latitude_dec", "Longitude_dec", "adm1", "adm2")])

table(Biogeographic_database_Ponerinae$duplicate)
# 63,807 unique localized records before biogeographic curation

## Save merged specimen database for Biogeographic data on Ponerinae
saveRDS(Biogeographic_database_Ponerinae, file = "./input_data/Biogeographic_data/Biogeographic_database_Ponerinae.rds")


##### 3/ Build taxa-level summary df #####

## Extent Phylogeny_sample_data_802 to include all valid species.
## Build several df with different levels of inclusion: all, phylo + valid.

### 3.1/ Load and prepare Phylogeny summary table and AntCat Species catalog for merging ####

## Load df for Metadata of samples in the phylogeny
Phylogeny_sample_data_802 <- readRDS(file = "./input_data/Phylogenies/Phylogeny_sample_data_802.rds")
Phylogeny_sample_data_802_for_merging <- Phylogeny_sample_data_802 %>% 
  rename(Specimen_phylogeny_AntWeb_name = AntWeb_name) %>%
  rename(Specimen_phylogeny_Code = "Specimen code") %>% 
  rename(Specimen_phylogeny_Extraction_code = "Extraction code") %>% 
  rename(Specimen_phylogeny_Name_updated_from_AntWeb = Name_updated_from_Antweb) %>% 
  rename(Specimen_phylogeny_Name_changed_from_phylogeny = Name_changed_from_phylogeny) %>% 
  rename(Status_AntCat = Name_status) %>% 
  rename(Country_Specimen_Phylogeny = country) %>% 
  rename(Bioregion_Specimen_Phylogeny = bioregion)
  
# Load AntCat list for valid Ponerinae species
Ponerinae_Species_AntCat_df <- readRDS("./input_data/AntWeb_data/Ponerinae_Species_AntCat_df.rds")

# Remove subspecies
Ponerinae_Species_AntCat_df <- Ponerinae_Species_AntCat_df[paste0(Ponerinae_Species_AntCat_df$genus, " ", Ponerinae_Species_AntCat_df$species) != Ponerinae_Species_AntCat_df$`current valid parent`, ]

Ponerinae_Species_AntCat_for_merging <- Ponerinae_Species_AntCat_df %>% 
  rename(Status_AntCat = status) %>% 
  rename(Country_type = country) %>% 
  rename(Bioregion_type = bioregion) %>% 
  rename(Current_name = Genus_species) %>% 
  select(Current_name, Status_AntCat, Country_type, Bioregion_type)

### 3.2/ Merge summary df ####

## Merge both summary df to add other valid species in the df
Ponerinae_Macroevolution_taxa_database <- full_join(Phylogeny_sample_data_802_for_merging, Ponerinae_Species_AntCat_for_merging) %>% 
  mutate(Available_occurrences = NA) %>%  # Create field to record availability of occurrence data in GABI + AntWeb
  mutate(To_include_in_analyses = NA) %>%  # Create field to record availability of occurrence data in GABI + AntWeb
  mutate(Conservative_clade_Node_ID = NA) %>%  # Create fields to record description, node ID and identity of terminal to use to recover the MRCA for a conservative clade we are sure the terminal belongs to
  mutate(Conservative_clade_Terminal_with_MRCA = NA) %>%
  mutate(Conservative_clade_Notes = NA) %>%
  mutate(Least_inclusive_clade_Node_ID = NA) %>%  # Create fields to record description, node ID and identity of terminal to use to recover the MRCA for the least inclusive clade we think the terminal may belong to
  mutate(Least_inclusive_clade_Terminal_with_MRCA = NA) %>%
  mutate(Least_inclusive_clade_Notes = NA) %>%
  select(Current_name, Status_AntCat, In_phylogeny, Specimen_phylogeny_Code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_Name_updated_from_AntWeb, Specimen_phylogeny_Extraction_code, Phylo_label, Phylo_name, Specimen_phylogeny_Name_changed_from_phylogeny, Current_duplicate, Subspecies, Conservative_clade_Node_ID, Conservative_clade_Terminal_with_MRCA, Conservative_clade_Notes, Least_inclusive_clade_Node_ID, Least_inclusive_clade_Terminal_with_MRCA, Least_inclusive_clade_Notes, Country_Specimen_Phylogeny, Country_type, Bioregion_Specimen_Phylogeny, Bioregion_type, Available_occurrences, To_include_in_analyses, Notes_2022_04_26, Notes_2023_11_28)

table(duplicated(Ponerinae_Macroevolution_taxa_database$Current_name))

Ponerinae_Macroevolution_taxa_database[duplicated(Ponerinae_Macroevolution_taxa_database$Current_name), ]

### 3.3/ Fill new fields ####

## Fill current fields for valid species outside of the phylogeny
Ponerinae_Macroevolution_taxa_database$In_phylogeny[is.na(Ponerinae_Macroevolution_taxa_database$In_phylogeny)] <- F
Ponerinae_Macroevolution_taxa_database$Subspecies[is.na(Ponerinae_Macroevolution_taxa_database$Subspecies)] <- F

## Format as Logical
Ponerinae_Macroevolution_taxa_database$In_phylogeny <- as.logical(Ponerinae_Macroevolution_taxa_database$In_phylogeny)
Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Name_updated_from_AntWeb <- as.logical(Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Name_updated_from_AntWeb)
Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Name_changed_from_phylogeny <- as.logical(Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Name_changed_from_phylogeny)
Ponerinae_Macroevolution_taxa_database$Current_duplicate <- as.logical(Ponerinae_Macroevolution_taxa_database$Current_duplicate)
Ponerinae_Macroevolution_taxa_database$Subspecies <- as.logical(Ponerinae_Macroevolution_taxa_database$Subspecies)

table(Ponerinae_Macroevolution_taxa_database$Subspecies)
table(Ponerinae_Macroevolution_taxa_database$In_phylogeny)
table(Ponerinae_Macroevolution_taxa_database$Current_duplicate)

## Fill fields for clade inclusion for terminals in the phylogeny
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Terminal_with_MRCA[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- Ponerinae_Macroevolution_taxa_database$Phylo_label[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] 
Ponerinae_Macroevolution_taxa_database$Least_inclusive_clade_Terminal_with_MRCA[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- Ponerinae_Macroevolution_taxa_database$Phylo_label[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] 
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Node_ID <- tidytree::nodeid(tree = Ponerinae_phylogeny_789t, label = Ponerinae_Macroevolution_taxa_database$Phylo_label)
Ponerinae_Macroevolution_taxa_database$Least_inclusive_clade_Node_ID <- tidytree::nodeid(tree = Ponerinae_phylogeny_789t, label = Ponerinae_Macroevolution_taxa_database$Phylo_label)
Ponerinae_Macroevolution_taxa_database$Conservative_clade_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- "Terminal in the current phylogeny. Conservative clade is the terminal itself. No need for grafting."
Ponerinae_Macroevolution_taxa_database$Least_inclusive_clade_Notes[(Ponerinae_Macroevolution_taxa_database$In_phylogeny == T)] <- "Terminal in the current phylogeny. Least inclusive clade is the termial itself. No need for grafting."

## Fill field for presence/absence of occurrences

AntWeb_database_taxa_names <- unique(AntWeb_database_curated$Genus_species)
GABI_database_taxa_names <- unique(GABI_database_Ponerinae$AntCat_Genus_species_name)

Ponerinae_Macroevolution_taxa_database$Available_occurrences <- (Ponerinae_Macroevolution_taxa_database$Current_name %in% c(AntWeb_database_taxa_names, GABI_database_taxa_names))

table(Ponerinae_Macroevolution_taxa_database$Available_occurrences)
# Only 11 taxa without occurrences (but bioregions can still be retrieved from type descriptions)

## Add field for inclusion or not in the analyses

# Only flag as FALSE the taxa with no geographic information in the database (GABI + AntWeb) and from AntCat (type)
No_geographic_info <- is.na(Ponerinae_Macroevolution_taxa_database$Bioregion_Specimen_Phylogeny) & is.na(Ponerinae_Macroevolution_taxa_database$Bioregion_type) & !(Ponerinae_Macroevolution_taxa_database$Available_occurrences)
table(No_geographic_info) # All taxa have geographic information!

# Flag as TRUE taxa in the phylogeny
Ponerinae_Macroevolution_taxa_database$To_include_in_analyses[Ponerinae_Macroevolution_taxa_database$In_phylogeny] <- T

## Save taxa-level summary df
saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")


### 3.4/ Merge with trait dataset ####

# Load taxa-level summary df
Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

# Load trait dataset
Trait_database <- read_excel("input_data/Traits_data/AoW_Ponerinae_measurements_2024_01_12.xlsx")
Trait_database <- Trait_database %>% 
  rename(Specimen_measured_Code = MeasSpecimen) %>% 
  select(-c("PmSP", "PrSP", "PtSP")) %>%   # Remove spinescence measurements that have been abandoned
  filter(GenSp != "Centromyrmex alfaroi") %>%  # Remove entry with no data
  filter(EXcode != "EX3093") %>%   
  filter(Specimen_measured_Code != "INBIOCRI0012278883") # Remove duplicate with wrong Specimen code


# Add new field to match names
Trait_database$Genus_species <- paste0(Trait_database$Genus, "_", Trait_database$Species) 

# Remove duplicates

table(duplicated(Trait_database$Specimen_measured_Code))
View(Trait_database[duplicated(Trait_database$Specimen_measured_Code), ])

Trait_database <- Trait_database[!duplicated(Trait_database$Specimen_measured_Code), ]

### 3.4.1/ Retrieve name of specimens from AntWeb database ####

AntWeb_database_extract <- AntWeb_database_curated[, c("code", "Genus_species", "AntWeb_name", "Status_AntCat")] %>% 
  rename(Current_name = Genus_species) %>% 
  mutate(code = str_to_upper(code))
Trait_database <- left_join(Trait_database, AntWeb_database_extract, by = c("Specimen_measured_Code" = "code"))
Trait_database$AntWeb_match <- !is.na(Trait_database$AntWeb_name) 

table(Trait_database$AntWeb_match)

# Trait_database <- Trait_database %>%
#   rename(Current_name = Current_name.y) %>%
#   select(-Current_name.x)
# Trait_database <- Trait_database %>%
#   rename(Status_AntCat = Status_AntCat.y) %>%
#   select(-Status_AntCat.x)

# Save curated Trait database
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")

# Load curated Trait database
Trait_database <- readRDS(file = "./input_data/Traits_data/Trait_database.rds")

### 3.4.2/ Fix error on case by case ####

View(Trait_database[!Trait_database$AntWeb_match, ])

## Case of Specimens in Introduced location (Not in the curated AntWeb dataset)
# CASENT0637780 # Ponera_swezeyi in AntWeb # Ponera_swezeyi in Phylogeny # Ponera_swezeyi in Current_name
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "Ponera_swezeyi"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "Ponera_swezeyi"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0637780"] <- "valid"
# CASENT0649763 # Brachyponera_chinensis in AntWeb # Brachyponera_chinensis in Phylogeny # Brachyponera_chinensis in Current_name 
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "Brachyponera_chinensis"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "Brachyponera_chinensis"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0649763"] <- "valid"

## Case of Specimens with error in the Specimen code
# # Specimen LACMENT140941. Not LACM ENT 140941.
# Trait_database$Specimen_measured_Code[Trait_database$Specimen_measured_Code == "LACM ENT 140941"] <- "LACMENT140941"
# Trait_database$Current_name[Trait_database$Specimen_measured_Code == "LACMENT140941"] <- "Ponera_coarctata"
# Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "LACMENT140941"] <- "Ponera_coarctata"
# Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "LACMENT140941"] <- "valid"
# # Specimen of Anochetus_hohenbergiae is indicated as the "holotype" which corresponds to Specimen CASENT0919827 found in AntWeb and used as voucher for molecular extraction.
# Trait_database$Specimen_measured_Code[Trait_database$Specimen_measured_Code == "holotype"] <- "CASENT0919827"
# Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0919827"] <- "Anochetus_hohenbergiae"
# Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0919827"] <- "Anochetus_hohenbergiae"
# Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0919827"] <- "valid"

## Case of Specimens with a Specimen code different from the extracted voucher, and not found in AntWeb
# Specimen ZRC_ENT00028294. Not in AntWeb but used as a substitute for voucher ZRC_ENT00047844 which is a Centromyrmex_hamulatus.
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "ZRC_ENT00028294"] <- "Centromyrmex_hamulatus"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "ZRC_ENT00028294"] <- "Centromyrmex_hamulatus"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "ZRC_ENT00028294"] <- "valid"

## Case of indeterminate not in the phylogeny, so need to be removed
# CASENT0650354 # Odontomachus_(indet) not used in the phylogeny
Trait_database <- Trait_database[!Trait_database$Specimen_measured_Code == "CASENT0650354", ]
# CASENT0872814 # Euponera_(indet) not used in the phylogeny
Trait_database <- Trait_database[!Trait_database$Specimen_measured_Code == "CASENT0872814", ]

## Case of Specimen that is not a Ponerinae
# CASENT0106229 # Amblyopone_australis is an Outgroup in the phylogeny. It is an Amblyoponinae.
Trait_database <- Trait_database[!Trait_database$Specimen_measured_Code == "CASENT0106229", ]

## New specimens in AntWeb. Should be present in the new updated AntWeb database
# CASENT0012708 si Odontomachus_simillimus
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "Odontomachus_simillimus"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "Odontomachus_simillimus"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "valid"
Trait_database$Specimen_phylogeny_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "Odontomachus_simillimus"
Trait_database$Phylo_name[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "Odontomachus_simillimus"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0012708"] <- "Odontomachus_simillimus"
# CASENT0818943 is Hypoponera_ragusai
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "Hypoponera_ragusai"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "Hypoponera_ragusai"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "valid"
Trait_database$Specimen_phylogeny_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "Hypoponera_ragusai"
Trait_database$Phylo_name[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "Hypoponera_ragusai"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0818943"] <- "Hypoponera_ragusai"
# CASENT0898545 is Promyopias_silvestrii
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT0898545"] <- "Promyopias_silvestrii"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0898545"] <- "Promyopias_silvestrii"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT0898545"] <- "valid"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT0898545"] <- "Promyopias_silvestrii"
# CASENT02706126 is Odontomachus_animosus
Trait_database$Current_name[Trait_database$Specimen_measured_Code == "CASENT02706126"] <- "Odontomachus_animosus"
Trait_database$Specimen_measured_AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT02706126"] <- "Odontomachus_animosus"
Trait_database$Status_AntCat[Trait_database$Specimen_measured_Code == "CASENT02706126"] <- "valid"
Trait_database$AntWeb_name[Trait_database$Specimen_measured_Code == "CASENT02706126"] <- "Odontomachus_animosus"


# Should be no more errors
View(Trait_database[!Trait_database$AntWeb_match, ])


### 3.4.3/ Match specimen based on Extraction code to detect presence in the phylogeny ####

# (Not Specimen code since they can use different specimens)

Trait_database$Taxa_in_phylogeny <- Trait_database$EXcode %in% Ponerinae_Macroevolution_taxa_database$Specimen_phylogeny_Extraction_code
table(Trait_database$Taxa_in_phylogeny)

# Compare names
Macroevol_database_extract <- Ponerinae_Macroevolution_taxa_database[, c("Current_name", "Specimen_phylogeny_Extraction_code", "Specimen_phylogeny_Code", "Specimen_phylogeny_AntWeb_name", "Phylo_name")]
Trait_database <- left_join(Trait_database, Macroevol_database_extract, by = c("EXcode" = "Specimen_phylogeny_Extraction_code"))

table((Trait_database$Current_name.x == Trait_database$Current_name.y))
test <- Trait_database[!(Trait_database$Current_name.x == Trait_database$Current_name.y),]
test <- test[!is.na(test$Current_name.x), ]
View(test) # Should be no entries with different Current names

# Clean the duplicate of Current-name fields
Trait_database <- Trait_database %>% 
  rename(Current_name = Current_name.x) %>% 
  select(-Current_name.y)
# Trait_database <- Trait_database %>%
#   rename(Specimen_phylogeny_Code = Specimen_phylogeny_Code.y) %>%
#   select(-Specimen_phylogeny_Code.x)
# Trait_database <- Trait_database %>%
#   rename(Specimen_phylogeny_AntWeb_name = Specimen_phylogeny_AntWeb_name.y) %>%
#   select(-Specimen_phylogeny_AntWeb_name.x)
# Trait_database <- Trait_database %>%
#   rename(Phylo_name = Phylo_name.y) %>%
#   select(-Phylo_name.x)

# Rename AntWeb names to distinguish between Specimen in the phylogeny and specimen measured
Trait_database <- Trait_database %>% 
  # select(-Specimen_measured_AntWeb_name) %>% 
  rename(Specimen_measured_AntWeb_name = AntWeb_name) %>% 
  rename(Specimen_measured_Name = Genus_species)


## Flag for discrepancy between specimen measured and specimen in phylogeny

Trait_database$Specimen_measured_is_in_phylogeny <- Trait_database$Specimen_measured_Code == Trait_database$Specimen_phylogeny_Code

### 3.4.4/ Add flags for partial and complete trait measurements ####

Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")] <- round(apply(X = Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")], MARGIN = 2, FUN = as.numeric), 3)
Trait_matrix <- Trait_database[, c("HW", "HL", "SL", "ED", "WL", "PW", "MtFL")]

Trait_database$Complete_trait_measurements <- complete.cases(Trait_matrix)
Trait_matrix_binary <- is.na(Trait_matrix)
Trait_database$Partial_trait_measurements <- apply(X = Trait_matrix_binary, MARGIN = 1, FUN = function (x) {(sum(x) > 0) & (sum(x) < length(x))})

table(Trait_database$Complete_trait_measurements)
table(Trait_database$Partial_trait_measurements)

## Save curated Trait database
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")


### 3.4.5/ Merge the curated trait dataset with the taxa-level summary df ####

names(Trait_database)

Trait_database_for_merging <- Trait_database %>% 
  rename(Notes_measurements = Notes) %>% 
  rename(Source_measurements = Source) %>% 
  rename(Specimen_measured_Extraction_code = EXcode) %>%
  # rename(Specimen_measured_Code = MeasSpecimen) %>%
  select(Current_name, Specimen_measured_Extraction_code, Specimen_measured_Code, Specimen_measured_Name, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, Complete_trait_measurements, Partial_trait_measurements, Source_measurements, Notes_measurements)

## Merge both summary df to add other valid species in the df

# # Remove fields that are already there only because I am rerunning the script.
# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
#   select(-Specimen_measured_Extraction_code, -Specimen_measured_Code, -Specimen_measured_Name, -Specimen_measured_AntWeb_name, -HW, -HL, -SL, -ED, -WL, -PW, -MtFL, -Complete_trait_measurements, -Partial_trait_measurements, -Source_measurements, -Notes_measurements)

nrow(Ponerinae_Macroevolution_taxa_database) # 1569 rows in the Ponerinae_Macroevolution_taxa_database including all taxa in phylogeny + other valid species name
Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, Trait_database_for_merging, by = "Current_name") %>% 
  mutate(Specimen_measured_used_in_phylogeny = Specimen_measured_Code == Specimen_phylogeny_Code) %>% 
  select(Current_name, Status_AntCat, In_phylogeny, Specimen_phylogeny_Code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_Name_updated_from_AntWeb, Specimen_phylogeny_Extraction_code, Phylo_label, Phylo_name, Specimen_phylogeny_Name_changed_from_phylogeny, Current_duplicate, Subspecies, Conservative_clade_Node_ID, Conservative_clade_Terminal_with_MRCA, Conservative_clade_Notes, Least_inclusive_clade_Node_ID, Least_inclusive_clade_Terminal_with_MRCA, Least_inclusive_clade_Notes, Country_Specimen_Phylogeny, Country_type, Bioregion_Specimen_Phylogeny, Bioregion_type, Available_occurrences, To_include_in_analyses, Notes_2022_04_26, Notes_2023_11_28, Specimen_measured_used_in_phylogeny, Specimen_measured_Extraction_code, Specimen_measured_Code, Specimen_measured_Name, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, Complete_trait_measurements, Partial_trait_measurements, Source_measurements, Notes_measurements)

nrow(Ponerinae_Macroevolution_taxa_database) # 1616 rows after merging
# Number of row is not equivalent because there are duplicates in the trait dataset. Need to remove them.


## Inspect case of specimen different in phylogeny and in trait measurement AND name in trait measurement disagree with AntWeb (so Antweb records may need to be changed)
# May need to go back and change the name of the Measured specimen in AntWeb and Trait dataset if it is a mistake...
test <- test[!test$Specimen_measured_used_in_phylogeny, ]
test <- test[!is.na(test$Current_name), ]
test <- test %>% 
  filter(Specimen_measured_Name != Specimen_measured_AntWeb_name) %>% 
  select(Current_name, Specimen_phylogeny_Code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_Extraction_code, Notes_2022_04_26, Notes_2023_11_28, Specimen_measured_Extraction_code, Specimen_measured_Code, Specimen_measured_Name, Specimen_measured_AntWeb_name, Notes_measurements, everything())
View(test) # Inspect for possible mistakes

## Case of possible mistake in the Trait dataset
# Specimen CASENT0842143, found in North of Mozambique is labeled as Anochetus_punctaticeps but recorded as Anochetus_talpa in AntWeb. Other Anochetus_punctaticeps are found n South Africa, except a dubious outlier in Cameroon. Other Anochetus_talpa can be found in Mozambique, so I assumed CASENT0842143 was an Anochetus_talpa as recorded in AntWeb.
# Specimen CASENT0634243, found in Sabah, Borneo, Malaysia is labeled as Leptogenys_iridescens but recorded as Leptogenys_borneensis in AntWeb. Leptogenys_borneensis has only records in Sabah while Leptogenys_iridescens can be found in Sarawak and Sulawesi, but not in Sabah.  Thus, I assumed CASENT0634243 was a Leptogenys_borneensis as recorded in AntWeb.
# Specimen CASENT0629612, found in Uganda is labeled as Fisheropone_ambigua but recorded as Mesoponera_ambigua in AntWeb. Both have large distribution across Africa so I believed AntWeb and kept it as a Mesoponera_ambigua.
# Specimen JTLC000008642, found in Queensland, Australia is labeled as Anochetus_armstrongi but recorded as Anochetus_rectangularis in AntWeb. Both species have large distribution across Australia, including Queensland, so I assumed AntWeb was right and kept it as a Anochetus_rectangularis.
# Specimen CASENT0650262 is labeled as Hypoponera_JTL030 but have been changed for Hypoponera_vc01 according to Jack’s comment in the other file.

## Record comments in the Trait dataset
# Export as Excel file. Write comments. Reimport. Save as .rds
openxlsx::write.xlsx(x = Trait_database, file = "./input_data/Traits_data/Trait_database.xlsx", overwrite = T)
Trait_database <- read_excel("./input_data/Traits_data/Trait_database.xlsx")

# Save updated df for Trait dataset
saveRDS(Trait_database, file = "./input_data/Traits_data/Trait_database.rds")


### 3.4.6/ Check for duplicates after curation ####

# Reorder rows per taxa and keep entry using three hierarchical criteria: 1/ the Specimen used in the phylogeny, 2/ Data obtained from specimens, 3/ At random.
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  arrange(Current_name, desc(Specimen_measured_used_in_phylogeny), desc(Source_measurements)) %>% 
  # select(Current_name, Specimen_measured_used_in_phylogeny, Source_measurements, everything()) %>%
  group_by(Current_name) %>% 
  mutate(Duplicates_counter = row_number(Current_name)) %>% 
  ungroup()

# Remove duplicates
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  filter(Duplicates_counter == 1 | Current_duplicate == T) %>%  # Remove duplicates from the trait dataset, but not duplicates from the phylogeny
  select(-Duplicates_counter)

nrow(Ponerinae_Macroevolution_taxa_database)
# 1571 rows again after cleaning of duplicates # If difference, it is because new taxa in the trait measurements need to be updated from the current AntWeb database I am using...

## Save taxa-level summary df
saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

## Export in Excel to provide in AoW folder
openxlsx::write.xlsx(x = Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")

### 3.5/ Retrieve grafting information from DropBox file ####

Ponerinae_Macroevolution_taxa_database_DropBox_version <- read_excel("input_data/Ponerinae_Macroevolution_taxa_database_DropBox_version.xlsx")

View(Ponerinae_Macroevolution_taxa_database_DropBox_version[which(!Ponerinae_Macroevolution_taxa_database_DropBox_version$Current_name %in% Ponerinae_Macroevolution_taxa_database$Current_name), ])
View(Ponerinae_Macroevolution_taxa_database[which(!Ponerinae_Macroevolution_taxa_database$Current_name %in% Ponerinae_Macroevolution_taxa_database_DropBox_version$Current_name), ])

# Extract grafting info from Dropbox file
Grafting_info <- Ponerinae_Macroevolution_taxa_database_DropBox_version %>% 
  select(Current_name, Conservative_clade_Node_ID, Conservative_clade_Terminal_with_MRCA, Conservative_clade_Notes, Least_inclusive_clade_Node_ID, Least_inclusive_clade_Terminal_with_MRCA, Least_inclusive_clade_Notes)

# Remove previous field to be replaced by the update
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  select(-Conservative_clade_Node_ID, -Conservative_clade_Terminal_with_MRCA, -Conservative_clade_Notes, -Least_inclusive_clade_Node_ID, -Least_inclusive_clade_Terminal_with_MRCA, -Least_inclusive_clade_Notes)

# Merge update
Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, Grafting_info)

# Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
#   filter(Duplicates_counter == 1) %>% 
#   select(-Duplicates_counter)
                                                    
### 3.6/ Retrieve "To_measure" information from DropBox file ####

To_measure_info <- read_excel("input_data/Ponerinae_Macroevolution_taxa_database_DropBox_version.xlsx", sheet = "to measure")

# Rename columns adequately
To_measure_info <- To_measure_info %>% 
  rename(Specimen_to_measure_Code = `specimen to measure`) %>%
  rename(Specimen_to_measure_Location = "where") %>%
  rename(Specimen_to_measure_Notes = "note") %>%
  select(Current_name, Specimen_to_measure_Code, Specimen_to_measure_Location, Specimen_to_measure_Notes)

# Merge with Macroevolution database
Ponerinae_Macroevolution_taxa_database <- left_join(Ponerinae_Macroevolution_taxa_database, To_measure_info)

## Add field to detect that measurements are still needed

Ponerinae_Macroevolution_taxa_database$To_measure <- is.na(Ponerinae_Macroevolution_taxa_database$Specimen_measured_Code)

table(Ponerinae_Macroevolution_taxa_database$To_measure)

## Add field to notify if only 'queens or males' are available, disqualifying them from being measured

# Based on Jack's saying, there should be all specimens that have no measurements, but that are not listed in its subsheet

Ponerinae_Macroevolution_taxa_database$Specimen_to_measure_Only_queens_or_males <- NA

In_Jack_list <- Ponerinae_Macroevolution_taxa_database$Current_name %in% To_measure_info$Current_name

test <- !In_Jack_list & Ponerinae_Macroevolution_taxa_database$To_measure
table(test) # 80 cases of taxa that are yet to be measured and not in Jack's list ! 
View(Ponerinae_Macroevolution_taxa_database[test, ])

## Export Macroevolution database = taxa-level summary df

# Sort columns in a logical order
Ponerinae_Macroevolution_taxa_database <- Ponerinae_Macroevolution_taxa_database %>% 
  select(Current_name, Status_AntCat, In_phylogeny, Specimen_phylogeny_Code, Specimen_phylogeny_AntWeb_name, Specimen_phylogeny_Name_updated_from_AntWeb, Specimen_phylogeny_Extraction_code, Phylo_label, Phylo_name, Specimen_phylogeny_Name_changed_from_phylogeny, Current_duplicate, Subspecies,
         Country_Specimen_Phylogeny, Country_type, Bioregion_Specimen_Phylogeny, Bioregion_type, Available_occurrences, To_include_in_analyses, Notes_2022_04_26, Notes_2023_11_28,
         Conservative_clade_Node_ID, Conservative_clade_Terminal_with_MRCA, Conservative_clade_Notes, Least_inclusive_clade_Node_ID, Least_inclusive_clade_Terminal_with_MRCA, Least_inclusive_clade_Notes, Specimen_measured_used_in_phylogeny, Specimen_measured_Extraction_code, Specimen_measured_Code, Specimen_measured_Name, Specimen_measured_AntWeb_name, HW, HL, SL, ED, WL, PW, MtFL, Complete_trait_measurements, Partial_trait_measurements, Source_measurements, Notes_measurements, To_measure, Specimen_to_measure_Code, Specimen_to_measure_Location, Specimen_to_measure_Only_queens_or_males, Specimen_to_measure_Notes)

## Save taxa-level summary df
saveRDS(Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

## Export in Excel to provide in AoW folder
openxlsx::write.xlsx(x = Ponerinae_Macroevolution_taxa_database, file = "./input_data/Ponerinae_Macroevolution_taxa_database.xlsx")



Ponerinae_Macroevolution_taxa_database <- readRDS(file = "./input_data/Ponerinae_Macroevolution_taxa_database.rds")

table(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements)
table(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements & Ponerinae_Macroevolution_taxa_database$In_phylogeny)
table(Ponerinae_Macroevolution_taxa_database$Complete_trait_measurements & !Ponerinae_Macroevolution_taxa_database$In_phylogeny)
