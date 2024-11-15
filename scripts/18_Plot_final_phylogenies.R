##### Script 18: Plot nice phylogenies for publication  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Plot UCE-based tree from IQ-TREE ML = ‘792-taxon ML tree’
  # Rectangular non-grafted, non-calibrated topology with node support values
     # Tip labels = with Specimen codes
     # Node colors = UFBS scores (> 99%, 99-95%, < 95%) ; Node shape = SH-like test success/fail (95% threshold)
     # Should we also display the non-partitioned, and loci-partitioned ones besides the SWSC one?
     # Include genus-groups and ant images but no geological epochs as it is not dated.
     # Add a substitutions/site scale
  # Rectangular non-grafted, but time-calibrated (MCMCTree + TreePL) = show the summarized MCC tree 
     # Rectangular phylo with with 95% HPD for node divergence time estimates
     # Include genus-groups and ant images
     # Include geological epochs and vertical time lines

# Plot nice time-calibrated phylogenies of publication
  # Rectangular phylo with 95% HPD and grafted edge colors
  # Circular (Fig. 1A) phylo with density mapping of bioregions + RISR crown nodes
  # Circular (Fig. 1B) phylo with BAMM rates and regime shifts 
  # All phylos include genus-groups, geological epochs, ant images

###

### Inputs

# Time-calibrated phylogeny from the Maximum Credibility Clade
# Node support metadata
# Genus-group metadata table

###

### Outputs

# Rectangular UCE-based phylos
   # with node support values (UFBS, SH-like test)
# Include genus-groups as encompassing color rectangles and ant images, but no geological epochs as it is not dated.

# Rectangular and Circular time-calibrated phylos
   # Rectangular with with 95% HPD and grafted edge colors
   # Circular with density mapping of bioregions (Figure 1A)
   # Circular with BAMM rates and regime shifts (Figure 1B)
# All phylos include genus-groups, geological epochs, ant images

###

# Clean environment
rm(list = ls())

library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)
library(xlsx)
library(ape)
library(phytools)
library(dispRity)
library(deeptime)  # Library to display geological scale on plot
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(qpdf)      # To merge PDF


##### 1/ Load input files ####

### 1.1/ Load UCE-based ‘792-taxon ML tree’ ####

# With the complete name, including extraction and specimen codes
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata <- read.beast.newick(file = "./input_data/Phylogenies/ponerinae-792t-spruce-75p-iqtree-swscmerge-mfp_v2_v2.tre")
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label

# Update Neoponera_bucki label
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label[str_detect(string = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label, pattern = "bucki")] <- "Neoponera_bucki_EX2455_DZUP549431"

# Save phylo in RDS
saveRDS(object = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo, file = "./input_data/Phylogenies/Ponerinae_uncalibrated_UCE_phylogeny_792t.rds")

### 1.2/ Load UCE-based 789t time-calibrated tree ####

## Load non-grafted MCC phylogeny with node metadata (and short names)
Ponerinae_MCC_phylogeny_789t_treedata <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")

Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label
View(Ponerinae_MCC_phylogeny_789t_treedata@data)

# Demultiplex the HPD range values
Ponerinae_MCC_phylogeny_789t_treedata@data$age_HPD_0.025 <- unlist(lapply(X = Ponerinae_MCC_phylogeny_789t_treedata@data$height_0.95_HPD, FUN = function (x) {x[1]} ))
Ponerinae_MCC_phylogeny_789t_treedata@data$age_HPD_0.975 <- unlist(lapply(X = Ponerinae_MCC_phylogeny_789t_treedata@data$height_0.95_HPD, FUN = function (x) {x[2]} ))

## Save non-grafted MCC phylogeny with node metadata
saveRDS(Ponerinae_MCC_phylogeny_789t_treedata, file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")
write.beast.newick(treedata = Ponerinae_MCC_phylogeny_789t_treedata, file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.tree")

# # Match with previous long names including Specimen codes and Extraction codes
# Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label
# 
# # Extract long names from the non-calibrated tree
# long_names <- Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label
# # Remove outliers
# long_names <- long_names[!long_names %in% c("Paraponera_clavata_EX1573_CASENT0633292", "Proceratium_google_MAMI0434_CASENT0035028", "Amblyopone_australis_D0872_CASENT0106229")]
# # Order alphabetically
# long_names <- long_names[order(long_names)]
# 
# # Extract short names
# short_names <- Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label
# short_names <- short_names[order(short_names)]
# 
# # Associate names and ID
# Ponerinae_789t_names_df <- data.frame(short_names = short_names, long_names = long_names)
# View(Ponerinae_789t_names_df)
# Ponerinae_789t_names_df$tip_ID <- match(x = short_names, table = Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label)
# Ponerinae_789t_names_df$inverse_ID <- match(x = Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label, table = short_names)
# 
# # Replace short names with long names
# Ponerinae_MCC_phylogeny_789t_treedata <- Ponerinae_MCC_phylogeny_789t_treedata_short_names
# Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label <- Ponerinae_789t_names_df$long_names[Ponerinae_789t_names_df$inverse_ID]


### 1.3/ Load the 'grafted MCC' 1534t tree’ ####

## Load grafted MCC phylogeny with node metadata
Ponerinae_MCC_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")

Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label
View(Ponerinae_MCC_phylogeny_1534t_treedata@data)

# Demultiplex the HPD range values
Ponerinae_MCC_phylogeny_1534t_treedata@data$age_HPD_0.025_backbone <- unlist(lapply(X = Ponerinae_MCC_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[1]} } ))
Ponerinae_MCC_phylogeny_1534t_treedata@data$age_HPD_0.975_backbone <- unlist(lapply(X = Ponerinae_MCC_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[2]} } ))

## Save non-grafted MCC phylogeny with node metadata
saveRDS(Ponerinae_MCC_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")
write.beast.newick(treedata = Ponerinae_MCC_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.tree")


## Load grafted MCC tree with short names, to make it cleaner with the grafted taxa
Ponerinae_MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")

# Remove the underscore in taxa names
Ponerinae_MCC_phylogeny_1534t_short_names$tip.label <- str_replace(string = Ponerinae_MCC_phylogeny_1534t_short_names$tip.label, pattern = "_", replacement = " ")

# Replace phylo with the short name version
Ponerinae_MCC_phylogeny_1534t_treedata@phylo <- Ponerinae_MCC_phylogeny_1534t_short_names

## Load the grafted tree with ARE metadata (with short names)
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE.rds")


### 1.4/ Prepare geological epochs metadata ####

# Use seven time slices associated to geological periods as in Kawahara et al., 2023: 
  # 0-5.33 My = Holocene + Pleistocene + Pliocene
  # 5.33-23.03 = Miocene 
  # 23.03-33.9 = Oligocene
  # 33.9-56.0 = Eocene
  # 56.0-66.0 = Paleocene
  # 66.0-100.5 = Late Cretaceous
  # 100.5-145.0 = Early Cretaceous (Root of Ponerinae = 113 My ?)

# Extract root age
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_789t_treedata@phylo))
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_treedata@phylo))
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo))

# Create custom epochs
epochs_custom_df <- data.frame(
  name = c("PPHo", "Miocene", "Oligo.", "Eocene", "Paleo.", "Late Cretaceous", "Early Cret."),
  true_name = c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cretaceous", "Early Cretaceous"),
  max_age = c(5.3, 23.0, 33.9, 55.8, 66.0, 100.5, root_age),
  min_age = c(0, 5.3, 23.0, 33.9, 55.8, 66.0, 100.5),
  color = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"),
  grey_color = c("white", "grey90", "white", "grey90", "white", "grey90", "white"))

# Convert age in time to root
epochs_custom_df$time_since_root_min <- root_age - epochs_custom_df$min_age
epochs_custom_df$time_since_root_max <- root_age - epochs_custom_df$max_age

##### 2/ Plot UCE-based 792t ML tree #####

### 2.1/ Prepare support values metadata ####

## Prepare template for metadata df 

nb_tips <- length(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$tip.label)
nb_nodes <- length(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$node.label)

Ponerinae_uncalibrated_UCE_phylogeny_792t_node_metadata_df <- as_tibble(data.frame(node = 1:(nb_tips + nb_nodes)))
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data <- Ponerinae_uncalibrated_UCE_phylogeny_792t_node_metadata_df

## Extract support values from node labels
UFBS_values <- str_remove(string = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$node.label, pattern = "/.*")
SH_aLRT_values <- str_remove(string = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$node.label, pattern = ".*/")

plot(UFBS_values, SH_aLRT_values)

# Add support values to metadata df
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS[(nb_tips+1):(nb_tips + nb_nodes)] <- UFBS_values
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT[(nb_tips+1):(nb_tips + nb_nodes)] <- SH_aLRT_values
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS <- as.numeric(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS)
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT <- as.numeric(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT)

## Explore UFBS values

hist(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS)
table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS)

## Use 95 and 100 as thresholds to define color scheme

Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status <- NA
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status[Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS > 0] <- "low"
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status[Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS >= 95] <- "high"
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status[Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS == 100] <- "max"

table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status)

## Explore SH-aLRT values

hist(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT)
table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT)

## Use 95 as threshold to define shape scheme

Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT_status <- NA
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT_status[Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS < 95] <- "failed"
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT_status[Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS >= 95] <- "success"

table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT_status)

table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$UFBS_status, Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$SH_aLRT_status)
# Redundant. Only shows UFBS values?


### 2.2/ Prepare genus-groups metadata #####

# 7 genus-groups from Schmidt & Shattuck, 2014
  # Platythyrea
  # Pachycondyla: Pachycondyla + Simopelta + Thaumatomyrmex + Belonopelta + Mayaponera + Raspone + Dinoponera + Neoponera
  # Ponera: Ponera + Diacamma + Emeryopone + Austroponera + Pseudoponera + Parvaponera + Wadeura + Cryptoponera + Iroponera + Ectomomyrmex
  # Harpegnathos
  # Hypoponera 
  # Plectroctena: Plectoponera + Centromyrmex + Psalidomyrmex + Loboponera + Boloponera
  # Odontomachus: A lot of stuff... (ca. 20 genera)

## 2.2.1/ Initiate metadata df ####

# List Genus-groups from Schmidt & Shattuck, 2014
genus_groups_list <- c("Platythyrea", "Pachycondyla", "Ponera", "Harpegnathos", "Hypoponera", "Plectroctena", "Odontomachus")

# List associated taxa to use to retrieve MRCA
genus_groups_MRCA_taxa_list <- list(Platythyrea = c("Platythyrea_arnoldi_D2118_CASENT0842276", "Platythyrea_schultzei_D2414a_CASENT0815127"),
                                    Pachycondyla = c("Simopelta_minima_EX2354_ANTWEB1038199", "Belonopelta_attenuata_EX2501_ICN100255"),
                                    Ponera = c("Diacamma_magdalenae_EX2692_ZRC_ENT00007536", "Austroponera_castanea_EX2322_CASENT0618835"),
                                    Harpegnathos = c("Harpegnathos_saltator_D0887_CASENT0179535", "Harpegnathos_my01_EX2981_CASENT0389113"),
                                    Hypoponera = c("Hypoponera_inexorata_EX2317_CASENT0647598", "Hypoponera_comis_D2470_CASENT0250654"),
                                    Plectroctena =  c("Centromyrmex_fugator_D2431_CASENT0235867", "Boloponera_ikemkha_EX2977_CASENT0254323"),
                                    Odontomachus = c("Neoponera_bucki_EX2455_DZUP549431", "Anochetus_talpa_D2422_CASENT0817421"))

# Initiate metadata df
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_node_ID <- NA
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_edge_ID <- NA
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$crown_height <- NA
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$stem_height <- NA

# Extract edge heights (not age, as not time-calibrated. Unit = substitutions/site)
all_edge_heights <- phytools::nodeHeights(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo)
root_height <- max(all_edge_heights[, 2])
all_edge_heights <- round(-1 * all_edge_heights + root_height, 5)

Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$Genus_group <- NA
for (i in 1:nrow(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- genus_groups_MRCA_taxa_list[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant nodes and tips
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo, node = MRCA_node_ID_i)
  
  # Store information of genus-group membership in the tree metadata
  Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$Genus_group[c(MRCA_node_ID_i, all_descendants_ID_i)] <- genus_group_i
  
  # Extract crown and stem heights
  MRCA_edge_ID_i <- which(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_height_i <- all_edge_heights[MRCA_edge_ID_i, 2]
  stem_height_i <- all_edge_heights[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and heights
  Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$crown_height[i] <- round(crown_height_i, 4)
  Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$stem_height[i] <- round(stem_height_i, 4)
  
}

View(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df)
table(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata@data$Genus_group)

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- genus_groups_list
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$color <- genus_groups_colors

## Save Genus-group metadata for 792t tree
saveRDS(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df.rds")

## Save treedata object 
saveRDS(Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata.rds")


### 2.3/ Rectangular plot ####

## Load treedata object for uncalibrated UCE 792t phylogeny
Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata.rds")

## Load Genus-group metadata for 792t tree
Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df.rds")

# Extract MRCA to flip position of Ponera and Pachycondyla
Ponera_MRCA_node_ID <- Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$group_name == "Ponera"]
Pachycondyla_MRCA_node_ID <- Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$group_name == "Pachycondyla"]


## Plot tree as background
plot_792t_phylo_rect <- ggtree(tr = Ponerinae_uncalibrated_UCE_phylogeny_792t_treedata,
                               layout = "rectangular",
                               ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                               size = 1.2)
                               # aes(color = (pp1 > 0.5))

## Add Genus-group cladelabels/colored rectangles
  
# Adjust colored rectangles
extend_rect_to <- 0.25

# Loop per genus-group
for (i in 1:nrow(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df))
{
  plot_792t_phylo_rect <- plot_792t_phylo_rect +
    
    #  Add Background color on the whole clade
    geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
                 fill = Ponerinae_uncalibrated_UCE_phylogeny_792t_all_genus_groups_metadata_df$color[i],
                 # type = "gradient", # gradient.direction = 'rt',
                 # align = "right",
                 extendto = extend_rect_to,
                 alpha = .5) # +
  
  # # Add a label to the median tip. Add a bar alongside all the clade tips
  # geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
  #                 # label = all_genus_groups_metadata_df$group_name[i], 
  #                 label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
  #                 align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
  #                 geom = 'text',     # Tip of label, with or without rectangle
  #                 fontsize = 7,      # Size of the text
  #                 fill = NA,         # Color of the label background
  #                 color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
  #                 angle = 270,
  #                 hjust = 'center',
  #                 offset = 5.0,
  #                 offset.text = 0.4,
  #                 barsize = 4)    # Width of the bar
}

## Add tree over genus-group background rectangles
plot_792t_phylo_rect <- plot_792t_phylo_rect +
  
  geom_tree(layout = "rectangular",
            size = 1.2) +
            # aes(color = (pp1 > 0.5))
  
  # ## Add nod support scores as side text
  # geom_text(aes(label = UFBS), hjust = -0.2, size = 3) +
  
  # ## Add nod support scores as label on the nod
  # geom_label(aes(label = UFBS, fill = UFBS_status),
  #            label.padding = unit(0.15, "lines"),
  #            hjust = 0.5, size = 3) +
  
  ## Add symbols as nod support
  geom_nodepoint(aes(fill = UFBS_status),
                 shape = 22,
                 size = 5) +
  
  ### Add species tip labels
  geom_tiplab(aes(label = label),
              color = "black",
              align = F,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 4,
              offset = 0.00,    # To move to the exterior. Useful if you wish to put multiple layer of labels with different offsets.
              hjust = 0)  +   # To center, left-align, right-aligned text
  
  ### Set title
  labs(title = paste0("IQTREE ML tree - 792t")) +
  
  ### Set legend fill scheme
  scale_fill_manual("UFBS support",
                    values = c("orange", "lightblue", "limegreen"),
                    labels = c("< 95%", "95% - 99%", "100%"),
                    breaks = c("low", "high", "max")
  ) +
  
  
  ### Add scale in substitutions per sites
  geom_treescale(x = 0.01, y = 780,
                 width = 0.02, # In substitutions per sites
                 label = "Substitutions / sites",
                 color = "black",
                 offset = 3.0,
                 offset.label = -3.0,
                 linesize = 4.0,
                 fontsize = 20) +
  
  # Avoid labels to be cutted
  coord_cartesian(clip = "off") +
  
  # Set aesthetics
  theme(plot.margin = unit(c(20, 120, 5, 20), "mm"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        # axis.line.x = element_line(color = "black", size = 2, vjust = 3),
        legend.title = element_text(size = 48, vjust = 3, hjust = 0, face = "bold"),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 40, face = "plain"),
        legend.key.size = unit(5.0,"lines"),
        # legend.key.width = unit(2.0, 'lines'),
        # legend.key.height = unit(5.0, 'lines'),
        legend.position = c(0.10, 0.92),
        # legend.position = unit(c(1.20, 0.50), "npc"),
        # legend.key = element_rect(fill = NA))
)

# Flip position of Ponera and Pachycondyla
plot_792t_phylo_rect <- ggtree::flip(tree_view = plot_792t_phylo_rect, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

## Plot final tree
pdf(file = "./outputs/Final_phylogenies/Ponerinae_IQTREE_phylogeny_792t_rect.pdf", height = 120, width = 40)

print(plot_792t_phylo_rect)

dev.off()


### To do in Postprod in Illustrator

# Add genus-group names as label
# Add cladelabel for Genera in alternating colors (beware of polyphyletic groups (Use numbers))
# Replace color rectangles with background fading rectangles
# Add nicer legend for nod symbols
# Add ant images with shadow effect


##### 3/ Plot UCE-based 789t time-calibrated MCC tree #####

## Load Genus-group metadata for 789t tree
Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_789t_all_genus_groups_UCE_only_metadata_df.rds")

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$group_name
Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$color <- genus_groups_colors


### 3.1/ Rectangular plot ####

# Extract MRCA to flip position of Ponera and Pachycondyla
Ponera_MRCA_node_ID <- Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$group_name == "Ponera"]
Pachycondyla_MRCA_node_ID <- Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$group_name == "Pachycondyla"]

## Plot tree as background
plot_789t_phylo_rect <- ggtree(tr = Ponerinae_MCC_phylogeny_789t_treedata,
                               layout = "rectangular",
                               ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                               size = 1.2)
# aes(color = (pp1 > 0.5))

## Add Genus-group cladelabels/colored rectangles

# # Adjust colored rectangles
# extend_rect_to <- 0.25

# Loop per genus-group to add labels/rectangles
for (i in 1:nrow(Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df))
{
  plot_789t_phylo_rect <- plot_789t_phylo_rect +
    
    # #  Add Background color on the whole clade
    # geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
    #              fill = Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$color[i],
    #              # type = "gradient", # gradient.direction = 'rt',
    #              # align = "right",
    #              extendto = extend_rect_to,
    #              alpha = .5) # +
  
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$group_name[i],
                    label = paste0('bold(',Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 20,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$color[i], Ponerinae_MCC_phylogeny_789t_all_genus_groups_metadata_df$color[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = 30.0,
                    offset.text = c(5.0, 5.0, 5.0, 11.0, 5.0, 5.0, 5.0)[i],
                    barsize = 20)    # Width of the bar
}

## Loop per geological epochs to add polygons
plot_789t_phylo_rect <- plot_789t_phylo_rect +
  
  geom_rect(data = epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
                xmin = time_since_root_min, xmax = time_since_root_max, 
                # fill = color),
                fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

## Add tree over geological epochs background rectangles
plot_789t_phylo_rect <- plot_789t_phylo_rect +
  
  geom_tree(layout = "rectangular",
            size = 1.2) +
  # aes(color = (pp1 > 0.5))
  
  # Add HPD 95% intervals on backbone nodes
  geom_range(range = "height_0.95_HPD", center = "height",
             color = "steelblue", size = 2, alpha = 0.5) +
  
  ### Add species tip labels
  geom_tiplab(aes(label = label),
              color = "black",
              align = F,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 4,
              offset = 0.00,    # To move to the exterior. Useful if you wish to put multiple layer of labels with different offsets.
              hjust = 0)  +   # To center, left-align, right-aligned text
  
  ### Set title
  labs(title = paste0("MCMCTree/TreePL MCC tree - 789t")) +
  
  # ### Add scale in My
  # geom_treescale(x = 0.01, y = 780,
  #                width = 20, # In substitutions per sites
  #                label = "My",
  #                color = "black",
  #                offset = 3.0,
  #                offset.label = -3.0,
  #                linesize = 4.0,
  #                fontsize = 20) +
  
  # Avoid labels to be cutted
  coord_cartesian(clip = "off") +
  
  # Set aesthetics
  theme(plot.margin = unit(c(20, 120, 5, 20), "mm"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        # axis.line.x = element_line(color = "black", size = 2, vjust = 3),
        legend.title = element_text(size = 48, vjust = 3, hjust = 0, face = "bold"),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 40, face = "plain"),
        legend.key.size = unit(5.0,"lines"),
        # legend.key.width = unit(2.0, 'lines'),
        # legend.key.height = unit(5.0, 'lines'),
        legend.position = c(0.10, 0.92),
        # legend.position = unit(c(1.20, 0.50), "npc"),
        # legend.key = element_rect(fill = NA))
)

## Add age lines
age_step <- 20
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_789t_treedata@phylo))
age <- age_step
while (age < root_age)
{ 
  plot_789t_phylo_rect <- plot_789t_phylo_rect +
    geom_vline(xintercept = age, 
               # linetype = 92,
               linetype = "dotted",
               linewidth = 1, color = 'grey20')
  age  <- age + age_step
}

## Add the geoscale
epochs_custom_df$min_age <- epochs_custom_df$time_since_root_min
epochs_custom_df$max_age <- epochs_custom_df$time_since_root_max

plot_789t_phylo_rect <- plot_789t_phylo_rect +
  
  # scale_x_reverse("Age (Mya)") +

  deeptime::coord_geo(dat = list(epochs_custom_df), lwd = 0.3,
                      size = 12, fontface = "bold",
                      height = unit(5.0, "line"),
                      expand = T,
                      # xlim = c(root_age, 0), 
                      # ylim = c(0, 789),
                      clip = "off")

# Flip position of Ponera and Pachycondyla
plot_789t_phylo_rect <- ggtree::flip(tree_view = plot_789t_phylo_rect, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

## Plot final tree
pdf(file = "./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_789t_rect_with_HPD.pdf", height = 120, width = 40)

print(plot_789t_phylo_rect)

dev.off()


### To do in Postprod in Illustrator

# Add genus-group names as label
# Add cladelabel for Genera in alternating colors (beware of polyphyletic groups (Use numbers))
# Replace color rectangles with background fading rectangles
# Add nice legend for 95% HPD of divergence time estimates
# Adjust geotime scale
# Add ant images with shadow effect


##### 4/ Plot the 1543t time-calibrated MCC tree #####

## Load Genus-group metadata for 1534t tree
Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df.rds")

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name
Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$color <- genus_groups_colors

### 4.1/ Rectangular plot for 95% HPD and UCE/grafted identity ####

Ponerinae_MCC_phylogeny_1534t_treedata@data

# Set the color scheme for UCE data vs. Grafted
# UCE_grafted_colors <- c("black", "red")
# UCE_grafted_colors <- c("black", "darkgoldenrod")
UCE_grafted_colors <- c("black", "grey70")

# Set color scheme for edges (Not used because use aes mapping to get a legend)
# branches_color <- UCE_grafted_colors[Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$missing_node + 1]
branches_color <- UCE_grafted_colors[Ponerinae_MCC_phylogeny_1534t_treedata@data$new_parental_edge + 1]

# Set color scheme for tips
# nb_tips <- length(Ponerinae_phylogeny_1534t_treedata_for_ARE@phylo$tip.label)
nb_tips <- length(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$tip.label)
# tips_color <- UCE_grafted_colors[Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$missing_node + 1]
tips_color <- UCE_grafted_colors[Ponerinae_MCC_phylogeny_1534t_treedata@data$new_parental_edge + 1]
tips_color <- tips_color[1:nb_tips]


## Plot tree as background
Ponerinae_phylogeny_plot <- ggtree(tr = Ponerinae_MCC_phylogeny_1534t_treedata,
                                   # tr = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE,
                                   layout = "rectangular",
                                   ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                                   linewidth = NA)

## Add Genus-group cladelabels/colored rectangles

# # Adjust colored rectangles
# extend_rect_to <- 0.25

# Loop per genus-group to add labels/rectangles
for (i in 1:nrow(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # #  Add Background color on the whole clade
    # geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
    #              fill = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$color[i],
    #              # type = "gradient", # gradient.direction = 'rt',
    #              # align = "right",
    #              extendto = extend_rect_to,
    #              alpha = .5) # +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],
                    label = paste0('bold(',Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 20,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$color[i], Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$color[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    # offset = 8.0,
                    offset = 2.0,
                    offset.text = c(5.0, 5.0, 5.0, 11.0, 5.0, 5.0, 5.0)[i],
                    barsize = 60)    # Width of the bar
                    # barsize = 20)    # Width of the bar
}

## Loop per geological epochs to add polygons
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_rect(data = epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
              xmin = time_since_root_min, xmax = time_since_root_max, 
              # fill = color),
              fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

## Plot tree, tip labels, node HPD ranges
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_tree(mapping = aes(colour = new_parental_edge), # Color as branch status (UCE or grafted)
            # mapping = aes(colour = missing_node), # For ARE treedata
            # color = branches_color,
            linewidth = 1.2) +
  
  # Add tip labels
  geom_tiplab(# mapping = aes(colour = missing_node),
    color = tips_color, # Color as tip status (UCE or grafted) 
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = UCE_grafted_colors, labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Add HPD 95% intervals on backbone nodes
  geom_range(range = "height_0.95_HPD_backbone", center = "height",
             color = "steelblue", size = 2, alpha = 0.5) +
  
  ### Set title
  labs(title = paste0("MCC grafted tree with 95% HPD - 1534t")) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.01,0.05,0,0), "npc"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        legend.title = element_text(size = 48, vjust = 0, hjust = 0, face = "bold"),
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(5.0,"lines"),
        legend.position = c(0.10, 0.92)
  )

## Add age lines
age_step <- 20
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_treedata@phylo))
age <- root_age - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, 
               # linetype = 92,
               linetype = "dotted",
               linewidth = 1, color = 'grey20')
  age  <- age - age_step
}

## Add the geoscale
epochs_custom_df$min_age <- epochs_custom_df$time_since_root_min
epochs_custom_df$max_age <- epochs_custom_df$time_since_root_max

Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  # scale_x_reverse("Age (Mya)") +
  
  deeptime::coord_geo(dat = list(epochs_custom_df), lwd = 0.3,
                      size = 12, fontface = "bold",
                      height = unit(5.0, "line"),
                      expand = T,
                      # xlim = c(root_age, 0), 
                      # ylim = c(0, 1534),
                      clip = "off")


# Export PDF
pdf(file = "./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_1534t_rect_HPD_and_grafting_identities.pdf", height = 120, width = 40)

print(Ponerinae_phylogeny_plot)

dev.off()


### To do in Postprod in Illustrator

# Add genus-group names as label
# Add cladelabel for Genera in alternating colors (beware of polyphyletic groups (Use numbers))
# Replace color rectangles with background fading rectangles
# Adjust geotime scale
# Add ant images with shadow effect


### 4.2/ Circular plot for ARE and RISR nodes ####

## Figure 1A

Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$tip.label

## Load backbone ARE pie lists
# backbone_ARE_pies_list <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_backbone_ARE_pies_list.rds")
backbone_ARE_pies_list <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_backbone_ARE_pies_list.rds")

# Extract only the MRCA of S&S genus-groups
genus_groups_ARE_pies_list <- backbone_ARE_pies_list[names(backbone_ARE_pies_list) %in% Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID]

## Load RISR metadata
RISR_clades_for_Bioregions_0.75_df <- readRDS(file = "./outputs/RISR/Ponerinae_MCC_phylogeny_1534t/RISR_clades_for_Bioregions_0.75_df.rds")

# Merge RISR data with treedata

Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE <- Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE %>% 
  left_join(RISR_clades_for_Bioregions_0.75_df)
View(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
areas_list <- c("A", "U", "E", "I", "R", "N", "W")
colors_list_for_areas <- colors_list_for_areas[areas_list]

bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name


# Initiate plot
Ponerinae_phylogeny_plot <- ggtree(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE,
                                   color = NA,
                                   open.angle = 8, # Doe snot work anti-clockwise
                                   layout = "fan") +
                                   # layout = "circular") + # Do use to get the circular closed background rectangles of geological epochs
  # Add title
  ggtitle(label = paste0("MCC grafted tree with ARE and RISR - 1534t"))

## Loop per geological epochs to add polygons
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_rect(data = epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
              xmin = time_since_root_min, xmax = time_since_root_max, 
              # fill = color),
              fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = "dotted",
               linewidth = 0.5, color = 'grey20')
  age  <- age - age_step
}

## Loop per bioregions
for (i in seq_along(bioregion_names))
{
  # i <- 1
  
  # Extract bioregion
  bioregion_i <- bioregion_names[i]
  
  # Find index in metadata df
  bioregion_ID_i <- which(names(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo) == paste0("mean_state_",bioregion_i))
  
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    geom_tree(color = colors_list_for_areas[i],
              alpha = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo[, bioregion_ID_i, drop = T],
              linewidth = 0.5)
}

## Add Genus-group cladelabels
for (j in 1:nrow(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[j],
                    # label = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[j], 
                    label = paste0('bold(',str_to_upper(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[j]),')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[j], genus_groups_colors[j]),   # Color of the bar and the text
                    angle = genus_groups_angles[j],
                    hjust = 'center',
                    offset = 6.8,
                    offset.text = genus_groups_offsets[j],
                    barsize = 25)    # Width of the bar
}

# Add RISR nodes
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  ggnewscale::new_scale_fill() +
  
  # geom_nodepoint(# data = RISR_clades_for_Bioregions_0.75_df, # Should work, but does not. So need need to be merged with treedata
  #                mapping = aes(#node = node,
  #                              fill = Region),
  #                # color = "black",
  #                size = 10, alpha = 1.0) +
  
  geom_point2(# data = RISR_clades_for_Bioregions_0.75_df, # Should work, but does not. So need need to be merged with treedata
              mapping = aes(subset = (node %in% RISR_clades_for_Bioregions_0.75_df$node), node = node, fill = Region),
              color = "black", shape = 22,
              size = 10, alpha = 1.0) +
  
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_areas), labels = names(colors_list_for_areas), values = colors_list_for_areas)

# Add tip labels
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  geom_tiplab(# mapping = aes(colour = missing_node),
    # color = tips_color,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0,
    size = 0.8,
    offset = 0.2,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0)      # To center, left-align, right-aligned text

## Does not work on circular/fan layouts...
# # Add ARE pie charts on backbone nodes
# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   geom_inset(insets = genus_groups_ARE_pies_list,
#              width = 0.2, height = 0.2,
#              hjust = 0, vjust = 0,
#              x = "node")

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(-0.05,-0.2,-0.1,-0.25), "npc"), # trbl
        plot.title = element_text(color = "black", size = 30, face = "bold", hjust = 0.5, vjust = -10, angle = 0),
        legend.title = element_text(size = 24, vjust = 3, hjust = 0, face = "bold"),
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        legend.background = element_rect(fill = NA),
        # legend.key.size = unit(5.0,"lines"),
        legend.position = c(0.90, 0.20),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(color = "transparent", fill = NA),
        panel.background = element_rect(color = NA, fill = "transparent"))

# Plot PDF
# pdf(file = paste0("./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_1534t_circular_Density_ARE_and_RISR_circular_template.pdf"), height = 16, width = 17)
pdf(file = paste0("./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_1534t_circular_Density_ARE_and_RISR.pdf"), height = 16, width = 17)

print(Ponerinae_phylogeny_plot)

dev.off()


### To do in Postprod in Illustrator

# Add ARE for genus-group from the rectangular layout
# Extend geological epochs rectangles from the circular layout
# Add full legend for bioregions
# Add legend for RISR (rectangles), ARE (circles), both (circles within rectangles)
# Add genus-group names as label
# Replace color rectangles with background fading rectangles
# Adjust geotime scale
# Add ant images with shadow effect


### 4.3/ Circular plot for BAMM rates and regime shifts ####

## Figure 1B

## 4.3.1/ Load data ####

Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$tip.label

# Use the MSC configuration for plotting the shifts
# Load metadata on rate shifts in the MSC tree
# MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_shifts_df.rds")
MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_shifts_df.rds")

# Merge MSC shifts data with treedata object
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE <- Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE %>% 
  left_join(MSC_shifts_df, by = join_by(node == tipward_node))
View(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo)


## Plot ‘core shifts’ with OR > 50? rather than MCS configuration? May fit better the color pattern

# Load phylogeny scaled according to the Odd-Ratio of marginal posterior probability / prior probability of shift (proportional to branch length)
branch_marginal_odd_ratios <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/branch_marginal_odd_ratios.rds")

# Check that both phylogenies are ordered similarly
attr(x = branch_marginal_odd_ratios, which = "order")
attr(x = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo, which = "order")
View(cbind(branch_marginal_odd_ratios$edge, Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$edge))

# Load edge rate shifts proba/OR
edge_rate_shift_probs_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/edge_rate_shift_probs_df.rds")

dim(edge_rate_shift_probs_df)
dim(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$edge)
dim(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo)

# Issue with using core shifts is that some shift have high probability to happen on successive branches but likely represent a unique shift in all configurations
edge_rate_shift_probs_df[(edge_rate_shift_probs_df$marginal_odd_ratios > 50), ]

table(edge_rate_shift_probs_df$marginal_odd_ratios > 50)
table(edge_rate_shift_probs_df$marg_posterior_probs > 0.5)

# Load BAMM output object
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

# Check BAMM output data are in the same order than the initial phylogeny
View(cbind(BAMM_output$edge, Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$edge))
table(BAMM_output$edge[, 1] == Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$edge[, 1])
table(BAMM_output$edge[, 2] == Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$edge[, 2])

## 4.3.2/ Compute mean/median node rates ####

## Easy version = use mean branch rates as the tipward node rates and smooth rates using the 'continuous' option in ggtree
 # Issue = will tend to displace true rate values toward the tips
## Complex version = design a function to extract mean rates at nodes, and not at edges!

# BAMMtools:::getMarginalBranchRateMatrix
# Branch-specific rates are the mean rates computed by integrating the relevant rate-through-time function along each branch, then dividing by the length of the branch.

source("./functions/get_BAMM_nodes_rates_summary.R")

# Compute node rates summary stats from BAMM output
BAMM_nodes_rates_summary_array <- get_BAMM_nodes_rates_summary(x = BAMM_posterior_samples_data, verbose = TRUE)

# Save node rates summary stats
saveRDS(object = BAMM_nodes_rates_summary_array, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_nodes_rates_summary_array.rds")

## Load node rates summary stats
BAMM_nodes_rates_summary_array <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_nodes_rates_summary_array.rds")

# Extract mean and median net diversification rates
BAMM_nodes_net_div_rates_df <- as.data.frame(BAMM_nodes_rates_summary_array[ , "net diversification", c("mean", "median")])
BAMM_nodes_net_div_rates_df$node_ID <- as.numeric(row.names(BAMM_nodes_net_div_rates_df))

# Inform node mean and median net diversification rates in the treedata object
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE <- left_join(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE, BAMM_nodes_net_div_rates_df, by = join_by(node == node_ID))
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo <- Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo %>% 
  rename(mean_net_div_rates = mean,
         median_net_div_rates = median)


## 4.3.3/ Plot circular phylogeny with BAMM rates ####

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
areas_list <- c("A", "U", "E", "I", "R", "N", "W")
colors_list_for_areas <- colors_list_for_areas[areas_list]

bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name

## Prepare color gradient
nb_col <- 101
# Get breaks from Jenks method
# scale_breaks <- BAMMtools::getJenksBreaks(var = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$mean_net_div_rates, k = nb_col)
scale_breaks <- BAMMtools::getJenksBreaks(var = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$median_net_div_rates, k = nb_col)
# Get breaks from quantile method
# scale_breaks <- as.numeric(quantile(x = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$mean_net_div_rates, probs = seq(from = 0, to = 1, length.out = nb_col)))
# Rescale values between 0 and 1
rescale_breaks <- scales::rescale(x = scale_breaks)
# Get colors by interpolating on the RColorBrewer palette
scale_colors <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"), space = 'Lab')(nb_col)
scale_colors <- rev(scale_colors)

## Initiate plot
Ponerinae_phylogeny_plot <- ggtree(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE,
                                   color = NA,
                                   open.angle = 8, # Doe snot work anti-clockwise
                                   layout = "fan") +
  # layout = "circular") + # Do use to get the circular closed background rectangles of geological epochs
  # Add title
  ggtitle(label = paste0("MCC grafted tree with BAMM rates and regime shifts - 1534t"))

## Loop per geological epochs to add polygons
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_rect(data = epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
              xmin = time_since_root_min, xmax = time_since_root_max, 
              # fill = color),
              fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = "dotted",
               linewidth = 0.5, color = 'grey20')
  age  <- age - age_step
}

# Plot tree with edge colored according to BAMM rates
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_tree(# mapping = aes(color = mean_net_div_rates),
            mapping = aes(color = median_net_div_rates),
            continuous = 'colour',
            linewidth = 0.5) +
  
  scale_color_gradientn(name = "Net div. rates",
                        values = rescale_breaks, # Should use the "values" arguments and not "breaks" (for discrete scales)
                        colours = scale_colors)
  
## Add Genus-group cladelabels
for (j in 1:nrow(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[j],
                    # label = Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[j], 
                    label = paste0('bold(',str_to_upper(Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name[j]),')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[j], genus_groups_colors[j]),   # Color of the bar and the text
                    angle = genus_groups_angles[j],
                    hjust = 'center',
                    offset = 6.8,
                    offset.text = genus_groups_offsets[j],
                    barsize = 25)    # Width of the bar
}

# # Add RISR nodes
# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   
#   ggnewscale::new_scale_fill() +
#   
#   # geom_nodepoint(# data = RISR_clades_for_Bioregions_0.75_df, # Should work, but does not. So need need to be merged with treedata
#   #                mapping = aes(#node = node,
#   #                              fill = Region),
#   #                # color = "black",
#   #                size = 10, alpha = 1.0) +
#   
#   geom_point2(# data = RISR_clades_for_Bioregions_0.75_df, # Should work, but does not. So need need to be merged with treedata
#     mapping = aes(subset = (node %in% RISR_clades_for_Bioregions_0.75_df$node), node = node, fill = Region),
#     color = "black", shape = 22,
#     size = 10, alpha = 1.0) +
#   
#   scale_fill_manual("Bioregion", breaks = names(colors_list_for_areas), labels = names(colors_list_for_areas), values = colors_list_for_areas)

# Add regime shifts from the Maximum Shift Credibility (MSC) configuration

Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  ggnewscale::new_scale_fill() +
  
  geom_point2(# data = MSC_shifts_df, # Should work, but does not. So need need to be merged with treedata
              mapping = aes(subset = (node %in% MSC_shifts_df$tipward_node[-1]), 
                            node = node, size = branch_marg_posterior_prob,
                            fill = rate_shift_type), # Color according to the absolute rate shift values
                            # fill = alpha_shift_type), # Color according to the change in time-variation trends 
              color = "black", shape = 21,
              alpha = 0.7) +
  
  scale_fill_manual("Regime shifts", breaks = c("decrease", "increase"),
                    labels = c("Decrease", "Increase"),
                    # values = c("dodgerblue", "brown1")) +
                    values = c(scale_colors[1], scale_colors[length(scale_colors)])) +
  

  scale_size_continuous("Marginal posterior probability\nof regime shift", 
                        range = c(8, 15))
  
  # guides(colour = guide_legend(order = 1), 
  #        size = guide_legend(order = 2),
  #                            # override.aes = list(size = 10)),
  #        fill = guide_legend(order = 3))

# BAMM rate shift: size = Marginal PP, color = increase (red) or decrease (blue)
# May add letter on shifts to describe them. Linked to table in SM. (Gave letters only to the shift in the ??? tree, but plot all 'core shifts'?)


# Add tip labels
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  geom_tiplab(# mapping = aes(colour = missing_node),
    # color = tips_color,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0,
    size = 0.8,
    offset = 0.2,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0)      # To center, left-align, right-aligned text

## Does not work on circular/fan layouts...
# # Add ARE pie charts on backbone nodes
# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   geom_inset(insets = genus_groups_ARE_pies_list,
#              width = 0.2, height = 0.2,
#              hjust = 0, vjust = 0,
#              x = "node")

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  # Adjust legend for taxa source
  theme(# plot.margin = unit(c(-0.05,-0.2,-0.1,-0.25), "npc"), # trbl # To preserve matching layout
        plot.margin = unit(c(-0.05,0.25,-0.1,-0.25), "npc"), # trbl # To include legend nicely
        plot.title = element_text(color = "black", size = 30, face = "bold", hjust = 0.5, vjust = -10, angle = 0),
        legend.title = element_text(size = 24, vjust = 3, hjust = 0, face = "bold"),
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(3.0,"lines"),
        # legend.position = c(0.25, 0.60), # To preserve matching layout
        legend.position = c(1.05, 0.50), # To include legend nicely
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(color = "transparent", fill = NA),
        panel.background = element_rect(color = NA, fill = "transparent"))

# Plot PDF
# pdf(file = paste0("./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_1534t_circular_BAMM_rates_and_regime_shifts.pdf"), height = 16, width = 17)
pdf(file = paste0("./outputs/Final_phylogenies/Ponerinae_MCC_phylogeny_1534t_circular_BAMM_rates_and_regime_shifts_with_legend.pdf"), height = 16, width = 17)

print(Ponerinae_phylogeny_plot)

dev.off()


# BAMM rate shift: size = Marginal PP, color = increase (red) or decrease (blue)
# May add letter on shifts to describe them. Linked to table in SM. (Gave letters only to the shift in the ??? tree, but plot all 'core shifts'?)


### To do in Postprod in Illustrator

# Extend geological epochs rectangles from the circular layout
# Add full legend for Regime shifts?
# Add genus-group names as label
# Replace color rectangles with background fading rectangles
# Adjust geotime scale
# Add ant images with shadow effect
# Add BAMM rates color scale above the geological line for circular phylo 
# Adjust position of regime shifts on the parental edge rather than the tipward node





