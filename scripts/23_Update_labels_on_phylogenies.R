
##### Script 23: Update taxonmy on phylogenies  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Update names of species on phylogenies according to latest revision in Fisher et al., 2025

###

### Inputs

# ML phylogram - 792t
# Time-calibrated MCC phylogeny - 789t
# 1000 posterior for Time-calibrated phylogenies - 789t
# Grafted time-calibrated MCC phylogeny - 1534t
# 1000 posterior for grafted time-calibrated phylogenies - 1534t

# List of taxonomic changes from Fisher et al., 2025 (DOI: 10.3897/zookeys.1264.173399)

###

### Outputs

## Supplementary Data S9
  # Same as above with updated labels
    # .rds
    # .tree
    # treedata objects
  # Associated supplementary figures (PDF)

###

# Clean environment
rm(list = ls())

library(tidyverse)
library(ape)
library(phytools)
library(openxlsx)
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(qpdf)      # To merge PDF

##### 1/ Load input files ####

### 1.1/ Load UCE-based ‘792-taxon ML tree’ ####

ML_phylogram_792t <- readRDS(file = "./input_data/Phylogenies/Ponerinae_uncalibrated_UCE_phylogeny_792t.rds")

### 1.2/ Load UCE-based 789t time-calibrated tree ####

MCC_phylogeny_789t_treedata <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")
MCC_phylogeny_789t <- MCC_phylogeny_789t_treedata@phylo

### 1.3/ Load all 1000 posterior 789t time-calibrated phylogenies ####

all_posteriors_789t <- readRDS(file = "input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.rds")

### 1.4/ Load the 'grafted time-calibrated MCC 1534t tree’ ####

# With full names
MCC_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")

# Update Neoponera_bucki label
MCC_phylogeny_1534t_treedata@phylo$tip.label[str_detect(string = MCC_phylogeny_1534t_treedata@phylo$tip.label, pattern = "bucki")] <- "Neoponera_bucki_EX2455_DZUP549431"

MCC_phylogeny_1534t <- MCC_phylogeny_1534t_treedata@phylo

# With short names
MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")

### 1.5/ Load all 1000 posterior grafted time-calibrated 1534t phylogenies ####

all_posteriors_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.rds")

### 1.6/ Load list of taxonomic changes from Fisher et al., 2025

Taxonomic_changes_df <- read.xlsx(xlsxFile = "./input_data/Phylogenies/Taxonomic_changes_Fisher_2025-s001.xlsx")


##### 2/ Build matching table for labels #####

### 2.1/ Extract labels from grafted phylogeny ####

All_labels_df <- data.frame(Full_labels = MCC_phylogeny_1534t$tip.label,
                            Short_labels = MCC_phylogeny_1534t_short_names$tip.label,
                            Status = "Grafted_ingroup")

All_labels_df$Status[All_labels_df$Full_labels != All_labels_df$Short_labels] <- "UCE_ingroup"

### 2.2/ Add outgroups ####

outgroups_full_labels <- ML_phylogram_792t$tip.label[!(ML_phylogram_792t$tip.label %in% All_labels_df$Full_labels)]
outgroups_short_labels <- str_split(string = outgroups_full_labels, pattern = "_", simplify = T)[, 1:2]
outgroups_short_labels <- paste0(outgroups_short_labels[, 1], "_", outgroups_short_labels[, 2])

outgroups_df <- data.frame(Full_labels = outgroups_full_labels, Short_labels = outgroups_short_labels, Status = "UCE_outgroup")

All_labels_df <- rbind(All_labels_df, outgroups_df)

### 2.3/ Assign updated names ####

## Build new names
Taxonomic_changes_df$Prior_name <- str_replace_all(string = Taxonomic_changes_df$Prior.placement, pattern = " ", replacement = "_")
Taxonomic_changes_df$Prior_Genus <- str_split(string = Taxonomic_changes_df$Prior_name, pattern = "_", simplify = T)[, 1]
Taxonomic_changes_df$Specific_epithet <- apply(X = str_split(string = Taxonomic_changes_df$Prior_name, pattern = "_", simplify = T)[, -1], MARGIN = 1, FUN = paste, collapse = "_")
Taxonomic_changes_df$Specific_epithet <- str_remove(string = Taxonomic_changes_df$Specific_epithet, pattern = "_$")

Taxonomic_changes_df$New_name <- paste0(Taxonomic_changes_df$Revised.genus, "_", Taxonomic_changes_df$Specific_epithet)
  
## Join into label df
All_labels_df <- left_join(x = All_labels_df, y = Taxonomic_changes_df[, c("Prior_name", "New_name")], by = join_by(Short_labels == Prior_name))

## Complete manually for morphospecies that are not listed in Fisher et al., 2025
write.xlsx(x = All_labels_df, file = "./input_data/Phylogenies/All_labels_df.xlsx")

SD9_Update_names_Fisher_2025_matching_table <- read.xlsx(xlsxFile = "./input_data/Phylogenies/SD9_Update_names_Fisher_2025_matching_table.xlsx")
saveRDS(SD9_Update_names_Fisher_2025_matching_table, file = "./input_data/Phylogenies/SD9_Update_names_Fisher_2025_matching_table.xlsx")


##### 3/ Update and export phylogenies #####

### 3.1/ Update UCE-based ‘792-taxon ML tree’ ####

Current_labels <- ML_phylogram_792t$tip.label
Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_full_label[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_label)]
Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]

cbind(Current_labels, Updated_labels)

ML_phylogram_792t_updated <- ML_phylogram_792t
ML_phylogram_792t_updated$tip.label <- Updated_labels

# Export
saveRDS(object = ML_phylogram_792t_updated, file = "./input_data/Phylogenies/SD_9.1_ML_phylogram_792t.rds")
write.tree(phy = ML_phylogram_792t_updated, file = "./input_data/Phylogenies/SD_9.1_ML_phylogram_792t.tree")

### 3.2/ Update UCE-based 789t time-calibrated tree ####

Current_labels <- MCC_phylogeny_789t$tip.label
Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_full_label[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_label)]
Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]

cbind(Current_labels, Updated_labels)

MCC_phylogeny_789t_updated <- MCC_phylogeny_789t
MCC_phylogeny_789t_updated$tip.label <- Updated_labels

MCC_phylogeny_789t_treedata_updated <- MCC_phylogeny_789t_treedata
MCC_phylogeny_789t_treedata_updated@phylo <- MCC_phylogeny_789t_updated

# Export phylogeny
saveRDS(object = MCC_phylogeny_789t_updated, file = "./input_data/Phylogenies/SD_9.2_Time-calibrated_MCC_phylogeny_789t.rds")
write.tree(phy = MCC_phylogeny_789t_updated, file = "./input_data/Phylogenies/SD_9.2_Time-calibrated_MCC_phylogeny_789t.tree")

# Export treedata
saveRDS(object = MCC_phylogeny_789t_treedata_updated, file = "./input_data/Phylogenies/SD_9.2_Time-calibrated_MCC_phylogeny_789t_treedata.rds")
write.beast.newick(treedata = MCC_phylogeny_789t_treedata_updated, file = "./input_data/Phylogenies/SD_9.2_Time-calibrated_MCC_phylogeny_789t_treedata.tree")

### 3.3/ Update all 1000 posterior 789t time-calibrated phylogenies ####

all_posteriors_789t_updated <- list()

for (i in seq_along(all_posteriors_789t))
{
  # i <- 1
  
  tree_i <- all_posteriors_789t[[i]]
  
  Current_labels <- tree_i$tip.label
  Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_full_label[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_label)]
  Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]
  
  # cbind(Current_labels, Updated_labels)
  
  tree_i$tip.label <- Updated_labels

  # Append list
  # all_posteriors_789t_updated <- c(all_posteriors_789t_updated, tree_i)
  all_posteriors_789t_updated <- append(x = all_posteriors_789t_updated, values = list(tree_i))
}
all_posteriors_789t_updated <- TreeTools::as.multiPhylo(all_posteriors_789t_updated)

# Export phylogenies
saveRDS(object = all_posteriors_789t_updated, file = "./input_data/Phylogenies/SD_9.3_1000_posteriors_time-calibrated_phylogenies_789t.rds")
write.tree(phy = all_posteriors_789t_updated, file = "./input_data/Phylogenies/SD_9.3_1000_posteriors_time-calibrated_phylogenies_789t.tree")

### 3.4/ Update the 'grafted time-calibrated MCC 1534t tree’ ####

Current_labels <- MCC_phylogeny_1534t$tip.label
Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_full_label[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_label)]
Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]

cbind(Current_labels, Updated_labels)

MCC_phylogeny_1534t_updated <- MCC_phylogeny_1534t
MCC_phylogeny_1534t_updated$tip.label <- Updated_labels

MCC_phylogeny_1534t_treedata_updated <- MCC_phylogeny_1534t_treedata
MCC_phylogeny_1534t_treedata_updated@phylo <- MCC_phylogeny_1534t_updated

# Export phylogeny
saveRDS(object = MCC_phylogeny_789t_updated, file = "./input_data/Phylogenies/SD_9.4_Grafted_time-calibrated_MCC_phylogeny_1534t.rds")
write.tree(phy = MCC_phylogeny_789t_updated, file = "./input_data/Phylogenies/SD_9.4_Grafted_time-calibrated_MCC_phylogeny_1534t.tree")

# Export treedata
saveRDS(object = MCC_phylogeny_789t_treedata_updated, file = "./input_data/Phylogenies/SD_9.4_Grafted_time-calibrated_MCC_phylogeny_1534t_treedata.rds")
write.beast.newick(treedata = MCC_phylogeny_789t_treedata_updated, file = "./input_data/Phylogenies/SD_9.4_Grafted_time-calibrated_MCC_phylogeny_1534t_treedata.tree")

## With short names

Current_labels <- MCC_phylogeny_1534t_short_names$tip.label
Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_name[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_name)]
Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]

cbind(Current_labels, Updated_labels)

MCC_phylogeny_1534t_short_names_updated <- MCC_phylogeny_1534t_short_names
MCC_phylogeny_1534t_short_names_updated$tip.label <- Updated_labels

# Export phylogeny with short names
saveRDS(object = MCC_phylogeny_1534t_short_names_updated, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names_updated.rds")


### 3.5/ Update all 1000 posterior grafted time-calibrated 1534t phylogenies ####

all_posteriors_1534t_updated <- list()

for (i in seq_along(all_posteriors_1534t))
{
  # i <- 1
  
  tree_i <- all_posteriors_1534t[[i]]
  
  Current_labels <- tree_i$tip.label
  Updated_labels <- SD9_Update_names_Fisher_2025_matching_table$New_full_label[match(x = Current_labels, table = SD9_Update_names_Fisher_2025_matching_table$Current_label)]
  Updated_labels[is.na(Updated_labels)] <- Current_labels[is.na(Updated_labels)]
  
  # cbind(Current_labels, Updated_labels)
  
  tree_i$tip.label <- Updated_labels
  
  # Append list
  # all_posteriors_1534t_updated <- c(all_posteriors_1534t_updated, tree_i)
  all_posteriors_1534t_updated <- append(x = all_posteriors_1534t_updated, values = list(tree_i))
}
all_posteriors_1534t_updated <- TreeTools::as.multiPhylo(all_posteriors_1534t_updated)

all_posteriors_1534t_updated[[105]]$tip.label

# Export phylogenies
saveRDS(object = all_posteriors_1534t_updated, file = "./input_data/Phylogenies/SD_9.5_1000_posteriors_grafted_time-calibrated_phylogenies_1534t.rds")
write.tree(phy = all_posteriors_1534t_updated, file = "./input_data/Phylogenies/SD_9.5_1000_posteriors_grafted_time-calibrated_phylogenies_1534t.tree")
