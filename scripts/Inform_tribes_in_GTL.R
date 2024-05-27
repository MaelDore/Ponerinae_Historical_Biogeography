##### Script 00: Inform tribe level in Global Taxa List file #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Use Catalogue of Life to retrieve tribe level of each taxa in the Global Taxa List file

### Inputs

# Global Taxa List file

### Sources

# Bánki, O., et al. (2024). Catalogue of Life Checklist (Version 2024-02-22). Catalogue of Life. https://doi.org/10.48580/dfvll

### Outputs

# Global Taxa List file with tribe level assigned

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load libraries ####

library(taxize)
library(tidyverse)
library(openxlsx)

### 1.2/ Load files ####

## Load GTL

Global_Taxon_List <- read.xlsx(xlsxFile = "./input_data/Global_Taxon_List_v3.xlsx", sheet = "taxa")

## Load AntCat taxa summary

AntCat_names <- read.xlsx(xlsxFile = "./input_data/AntWeb_data/AntCat_names_2023_11_23.xlsx")

### 1.3/ Set key for NCBI ###

usethis::edit_r_environ()
# Paste this in the .Renviron 
# ENTREZ_KEY='d6dd1d73965cd2f3978cffc4863247ba7308'

##### 2/ Find tribes of all Genera ####

?taxize::col_classification # CoL API is deprecated
?taxize::classification

# Extract Genera

genera_list <- unique(Global_Taxon_List$genus)

### 2.1/ Get tribe from NCBI ####

tribe_NCBI_raw_output_1_100 <- taxize::classification(genera_list[1:100], db = "ncbi")
tribe_NCBI_raw_output_101_200 <- taxize::classification(genera_list[101:200], db = "ncbi")
tribe_NCBI_raw_output_201_293 <- taxize::classification(genera_list[201:293], db = "ncbi")
tribe_NCBI_raw_output_294_300 <- taxize::classification(genera_list[294:300], db = "ncbi")
tribe_NCBI_raw_output_301_359 <- taxize::classification(genera_list[301:359], db = "ncbi")

tribe_NCBI_raw_output <- append(x = tribe_NCBI_raw_output_1_100, values = tribe_NCBI_raw_output_101_200)
tribe_NCBI_raw_output <- append(x = tribe_NCBI_raw_output, values = tribe_NCBI_raw_output_201_293)
tribe_NCBI_raw_output <- append(x = tribe_NCBI_raw_output, values = tribe_NCBI_raw_output_294_300)
tribe_NCBI_raw_output <- append(x = tribe_NCBI_raw_output, values = tribe_NCBI_raw_output_301_359)

names(tribe_NCBI_raw_output)

# Summarize results in a df
summary_genera_classif_df <- data.frame(Genus = genera_list, Genus_NCBI = names(tribe_NCBI_raw_output), Match_NCBI = NA, Tribe_NCBI = NA)

# Record results
for (i in 1:length(tribe_NCBI_raw_output))
{
  # i <- 2
  
  taxa_df <- tribe_NCBI_raw_output[[i]]
  
  summary_genera_classif_df$Match_NCBI[i] <- FALSE
  
  # Does it have a match
  if (length(taxa_df) > 1)
  {
    # Does it have an identified Family?
    if ("family" %in% taxa_df$rank)
    {
      family_index <- which(taxa_df$rank == "family")
      family_name <- taxa_df$name[family_index]
      
      # Is it an ant?
      if (family_name == "Formicidae")
      {
        summary_genera_classif_df$Match_NCBI[i] <- TRUE
        
        # Does it have a tribe?
        if ("tribe" %in% taxa_df$rank)
        {
          tribe_index <- which(taxa_df$rank == "tribe")
          
          # Record the tribe
          summary_genera_classif_df$Tribe_NCBI[i] <- taxa_df$name[tribe_index]
        }
      }
    }
  }
}

View(summary_genera_classif_df)


### 2.2/ Get tribe from ITIS ####

tribe_ITIS_raw_output <- taxize::classification(genera_list, db = "itis")

# Initialize summary variables
summary_genera_classif_df$Genus_ITIS <- names(tribe_ITIS_raw_output)
summary_genera_classif_df$Match_ITIS <- NA
summary_genera_classif_df$Tribe_ITIS <- NA

# Record results
for (i in 1:length(tribe_ITIS_raw_output))
{
  # i <- 2
  
  taxa_df <- tribe_ITIS_raw_output[[i]]
  
  summary_genera_classif_df$Match_ITIS[i] <- FALSE
  
  # Does it have a match
  if (length(taxa_df) > 1)
  {
    # Does it have an identified Family?
    if ("family" %in% taxa_df$rank)
    {
      family_index <- which(taxa_df$rank == "family")
      family_name <- taxa_df$name[family_index]
      
      # Is it an ant?
      if (family_name == "Formicidae")
      {
        summary_genera_classif_df$Match_ITIS[i] <- TRUE
        
        # Does it have a tribe?
        if ("tribe" %in% taxa_df$rank)
        {
          tribe_index <- which(taxa_df$rank == "tribe")
          
          # Record the tribe
          summary_genera_classif_df$Tribe_ITIS[i] <- taxa_df$name[tribe_index]
        }
      }
    }
  }
}

View(summary_genera_classif_df)

### 2.3/ Get tribe from AntCat ####

# Keep only unique genus entry in AntCat
AntCat_names_Genera <- AntCat_names %>% 
  select(genus, tribe) %>%
  distinct(genus, tribe, .keep_all = TRUE) %>%
  rename(Tribe_AntCat = tribe) %>% 
  drop_na()

# Match with summary df
summary_genera_classif_df <- left_join(x = summary_genera_classif_df, y = AntCat_names_Genera, by = c("Genus" = "genus"))

summary_genera_classif_df$Match_AntCat <- !is.na(summary_genera_classif_df$Tribe_AntCat)


#### 3/ Correct errors ####

### 3.1/ Retrieve missing AntCat from ITIS ####

View(summary_genera_classif_df[is.na(summary_genera_classif_df$Tribe_AntCat), ])


summary_genera_classif_df$Current_tribe <- summary_genera_classif_df$Tribe_AntCat
summary_genera_classif_df$Current_tribe[is.na(summary_genera_classif_df$Tribe_AntCat)] <- summary_genera_classif_df$Tribe_ITIS[is.na(summary_genera_classif_df$Tribe_AntCat)]
  
### 3.2/ Correct missing ####

summary_genera_classif_df$Genus[is.na(summary_genera_classif_df$Current_tribe)]
View(summary_genera_classif_df[is.na(summary_genera_classif_df$Current_tribe), ])

#### 4/ Fill initial GTL file ####

# Join
Global_Taxon_List_with_tribes <- left_join(x = Global_Taxon_List, y = summary_genera_classif_df[, c("Genus", "Tribe_ITIS", "Tribe_AntCat", "Current_tribe")], by = c("genus" = "Genus"))

# Summary
table(Global_Taxon_List_with_tribes$Current_tribe)
table(is.na(Global_Taxon_List_with_tribes$Current_tribe))

table(Global_Taxon_List_with_tribes$subfamily[is.na(Global_Taxon_List_with_tribes$Current_tribe)])
# Most entry missing tribes are Dorylinae

# Export
write.xlsx(Global_Taxon_List_with_tribes, file = "./input_data/Global_Taxon_List_with_tribes.xlsx")
