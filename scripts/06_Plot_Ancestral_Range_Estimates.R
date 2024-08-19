##### Script 06: Plot ancestral ranges #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Plot complete, UCE-only, and genus-group collapsed phylogenies
# Plot ancestral range estimates
   # At nodes = marginal likelihoods from BioGeoBEARS models
   # Along branches = posterior probabilities of bioregion membership from Biogeographic Stochastic Maps
# Collapse at genus (group) level

###

### Inputs

# ARE results from BioGeoBEARS models
# Simmaps of Biogeographic Stochastic Maps

###

### Outputs

# Plot complete, UCE-only, and genus-group collapsed phylogenies
## Plot ancestral range estimates at nodes
  # Overall and collapsed at Genus-level
## Plot ancestral range estimates along branches
  # Posterior probability per bioregions + Aggregated
## Plot nodes and branches together

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(phytools)
library(ape)
library(BioGeoBEARS)
library(qpdf)
library(magick)   # For animated GIF
library(pdftools)  # To read PDF
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(ggimage)   # To plot pies as insets on ggtree

### 1.2/ Load BioGeoBEARS results ####

# Load time-stratified DEC+J model output
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")

# Inspect marginal likelihoods of ancestral states
dim(DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each branch tipward end = each node/tip, given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start/rootward end of each branch, just after cladogenesis/speciation (and eventually cladogenetic transition).


### 1.3/ Load phylogeny ####

# Load imputed phylogeny with short names
Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")

# Reorder in "cladewise" order to match with the Simmap data
Ponerinae_phylogeny_1534t_short_names_cladewise <- reorder(Ponerinae_phylogeny_1534t_short_names, order = "cladewise")

# Load list of Missing taxa that were imputed on the backbone phylogeny
Missing_taxa_list <- readRDS(file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")
Missing_taxa_list <- Missing_taxa_list$Taxa_name

##### 2/ Detect tips and branches of imputed taxa ####

### 2.1/ Prepare df of species (tips) data

tips_metadata_df <- data.frame(label = Ponerinae_phylogeny_1534t_short_names_cladewise$tip.label)
tips_metadata_df <- tips_metadata_df %>% 
  mutate(missing_taxa = label %in% Missing_taxa_list)

### 2.2/ Transform into a treedata object (treeio package)
Ponerinae_phylogeny_1534t_treedata <- Ponerinae_phylogeny_1534t_short_names_cladewise %>%
  as.treedata() %>%
  full_join(x = ., y = tips_metadata_df, by = "label") # Associate the list_sp table data to each tip

# Node metadata is stored in @extraInfo
View(Ponerinae_phylogeny_1534t_treedata@extraInfo)

# Add status of each nodes based on its index
Ponerinae_phylogeny_1534t_treedata@extraInfo$status <- NA
Ponerinae_phylogeny_1534t_treedata@extraInfo$status[1:length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)] <- "tip"
Ponerinae_phylogeny_1534t_treedata@extraInfo$status[(length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label) + 2):(length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label) + Ponerinae_phylogeny_1534t_treedata@phylo$Nnode)] <- "internal_node"
Ponerinae_phylogeny_1534t_treedata@extraInfo$status[(length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label) + 1)] <- "root"
table(Ponerinae_phylogeny_1534t_treedata@extraInfo$status)

# Add parental edge ID
Ponerinae_phylogeny_1534t_treedata@extraInfo$parental_edge_ID <- match(Ponerinae_phylogeny_1534t_treedata@extraInfo$node, table = Ponerinae_phylogeny_1534t_treedata@phylo$edge[, 2])

# Ponerinae_phylogeny_1534t_treedata@extraInfo[1530:1550,]

# Save treedata for grafted phylogeny
# saveRDS(Ponerinae_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")


### 2.3/ Prune tree from grafted taxa

Ponerinae_phylogeny_789t_treedata <- tidytree::drop.tip(object = Ponerinae_phylogeny_1534t_treedata, tip = Missing_taxa_list)

# Save treedata for UCE only phylogeny
saveRDS(Ponerinae_phylogeny_789t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_789t_treedata.rds")


### 2.4/ Detect origin of each node (and associated descending edge) ####

all_missing_taxa_ID <- which(tips_metadata_df$missing_taxa)
Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node <- NA

for (i in 1:nrow(Ponerinae_phylogeny_1534t_treedata@extraInfo))
{
  # i <- 1
  
  # Extract node ID
  node_ID_i <- Ponerinae_phylogeny_1534t_treedata@extraInfo$node[i]
    
  # Extract all current descendants
  descendants_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo, node = node_ID_i)
  
  # Check if all current descendants are missing taxa
  missing_i <- all(descendants_i %in% all_missing_taxa_ID)
  
  # Record status
  Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node[i] <- missing_i
}
table(Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node)


# Save treedata for complete phylogeny
saveRDS(Ponerinae_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")


##### 3/ Prepare metadata for Genera #####

### 3.1/ Initiate metadata df ####

# Extract all genera from tip labels
all_genera_list <- str_split(string = Ponerinae_phylogeny_1534t_treedata@phylo$tip.label, pattern = "_", simplify = T)[,1]
unique_genera_list <- unique(all_genera_list)
unique_genera_list <- unique_genera_list[order(unique_genera_list)]

# Initiate metadata df
all_genera_metadata_df <- data.frame(genus_name = unique_genera_list)
all_genera_metadata_df$MRCA_node_ID <- NA
all_genera_metadata_df$MRCA_edge_ID <- NA
all_genera_metadata_df$crown_age <- NA
all_genera_metadata_df$stem_age <- NA
all_genera_metadata_df$current_richness <- NA
all_genera_metadata_df$missing_taxa_nb <- NA
all_genera_metadata_df$missing_taxa_perc <- NA

### 3.2/ Extract genera metadata from tree ####

all_edge_ages <- phytools::nodeHeights(Ponerinae_phylogeny_1534t_treedata@phylo)
root_age <- max(all_edge_ages[, 2])
all_edge_ages <- round(-1 * all_edge_ages + root_age, 5)


for (i in 1:nrow(all_genera_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus
  genus_i <- all_genera_metadata_df$genus_name[i]
  
  # Get current species list
  genus_match_i <- which(all_genera_list == genus_i)
  sp_match_i <- Ponerinae_phylogeny_1534t_treedata@phylo$tip.label[genus_match_i]
  missing_sp_match_i <- Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_taxa[genus_match_i]
  
  # Record current richness data
  all_genera_metadata_df$current_richness[i] <- length(sp_match_i)
  all_genera_metadata_df$missing_taxa_nb[i] <- sum(missing_sp_match_i)
  all_genera_metadata_df$missing_taxa_perc[i] <- round(all_genera_metadata_df$missing_taxa_nb[i] / all_genera_metadata_df$current_richness[i] * 100, 1)
  
  # Get MRCA node 
  if (length(sp_match_i) > 1) 
  {
    # If genus has multiple current species, identify MRCA has the oldest MRCA of all pairs
    sp_pairs_i <- expand.grid(sp_match_i, sp_match_i)
    MRCA_all_pairs_i <- apply(X = sp_pairs_i, MARGIN = 1, FUN = ape::getMRCA, phy = Ponerinae_phylogeny_1534t_treedata@phylo)
    all_MRCA_nodes_pairs_i <- unique(MRCA_all_pairs_i)
    all_MRCA_edges_i <- match(x = all_MRCA_nodes_pairs_i, table = Ponerinae_phylogeny_1534t_treedata@phylo$edge[, 2])
    if (length(all_MRCA_edges_i) > 1)
    {
      MRCA_all_pairs_ages_i <- all_edge_ages[all_MRCA_edges_i, ]
      MRCA_edge_ID_i <- all_MRCA_edges_i[which.max(MRCA_all_pairs_ages_i[,2])]
      MRCA_node_ID_i <- Ponerinae_phylogeny_1534t_treedata@phylo$edge[MRCA_edge_ID_i, 2]
    } else {
      MRCA_node_ID_i <- all_MRCA_nodes_pairs_i
    }
  } else {
    # If genus is monospecific, MRCA is the species tip
    MRCA_node_ID_i <- genus_match_i

  }
  
  # Extract crown and stem ages
  MRCA_edge_ID_i <- which(Ponerinae_phylogeny_1534t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_age_i <- all_edge_ages[MRCA_edge_ID_i, 2]
  stem_age_i <- all_edge_ages[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and ages
  all_genera_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  all_genera_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  all_genera_metadata_df$crown_age[i] <- round(crown_age_i, 1)
  all_genera_metadata_df$stem_age[i] <- round(stem_age_i, 1)
  
}

View(all_genera_metadata_df)

# Save Genera metadata
saveRDS(all_genera_metadata_df, file = "./outputs/Grafting_missing_taxa/all_genera_metadata_df.rds")


##### 4/ Prepare metadata for Genus-groups from Schmidt & Shattuck, 2014 #####

# 7 genus-groups from Schmidt & Shattuck, 2014
  # Platythyrea
  # Pachycondyla: Pachycondyla + Simopelta + Thaumatomyrmex + Belonopelta + Mayaponera + Raspone + Dinoponera + Neoponera
  # Ponera: Ponera + Diacamma + Emeryopone + Austroponera + Pseudoponera + Parvaponera + Wadeura + Cryptoponera + Iroponera + Ectomomyrmex
  # Harpegnathos
  # Hypoponera 
  # Plectroctena: Plectoponera + Centromyrmex + Psalidomyrmex + Loboponera + Boloponera
  # Odontomachus: A lot of stuff... (ca. 20 genera)
  
 


### 4.1/ Initiate metadata df ####

# List Genus-groups from Schmidt & Shattuck, 2014
genus_groups_list <- c("Platythyrea", "Pachycondyla", "Ponera", "Harpegnathos", "Hypoponera", "Plectroctena", "Odontomachus")

# List associated taxa to use to retrieve MRCA
genus_groups_MRCA_taxa_list <- list(Platythyrea = c("Platythyrea_arnoldi", "Platythyrea_schultzei"),
                                    Pachycondyla = c("Simopelta_minima", "Belonopelta_attenuata"),
                                    Ponera = c("Diacamma_magdalenae", "Austroponera_castanea"),
                                    Harpegnathos = c("Harpegnathos_saltator", "Harpegnathos_my01"),
                                    Hypoponera = c("Hypoponera_inexorata", "Hypoponera_comis"),
                                    Plectroctena =  c("Centromyrmex_fugator", "Boloponera_ikemkha"),
                                    Odontomachus = c("Neoponera_bucki", "Anochetus_talpa"))

# Initiate metadata df
all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
all_genus_groups_metadata_df$MRCA_node_ID <- NA
all_genus_groups_metadata_df$MRCA_edge_ID <- NA
all_genus_groups_metadata_df$crown_age <- NA
all_genus_groups_metadata_df$stem_age <- NA
all_genus_groups_metadata_df$current_richness <- NA
all_genus_groups_metadata_df$missing_taxa_nb <- NA
all_genus_groups_metadata_df$missing_taxa_perc <- NA


### 4.2/ Extract Genus-groups metadata from complete tree ####

all_edge_ages <- phytools::nodeHeights(Ponerinae_phylogeny_1534t_treedata@phylo)
root_age <- max(all_edge_ages[, 2])
all_edge_ages <- round(-1 * all_edge_ages + root_age, 5)

for (i in 1:nrow(all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- genus_groups_MRCA_taxa_list[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_phylogeny_1534t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant taxa
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo, node = MRCA_node_ID_i)
  all_descendants_ID_i <- all_descendants_ID_i[all_descendants_ID_i <= length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)]
  all_descendants_names_i <- Ponerinae_phylogeny_1534t_treedata@phylo$tip.label[all_descendants_ID_i]
  missing_sp_match_i <- Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_taxa[all_descendants_ID_i]
  
  # Record current richness data
  all_genus_groups_metadata_df$current_richness[i] <- length(all_descendants_names_i)
  all_genus_groups_metadata_df$missing_taxa_nb[i] <- sum(missing_sp_match_i)
  all_genus_groups_metadata_df$missing_taxa_perc[i] <- round(all_genus_groups_metadata_df$missing_taxa_nb[i] / all_genus_groups_metadata_df$current_richness[i] * 100, 1)
  
  # Extract crown and stem ages
  MRCA_edge_ID_i <- which(Ponerinae_phylogeny_1534t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_age_i <- all_edge_ages[MRCA_edge_ID_i, 2]
  stem_age_i <- all_edge_ages[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and ages
  all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  all_genus_groups_metadata_df$crown_age[i] <- round(crown_age_i, 1)
  all_genus_groups_metadata_df$stem_age[i] <- round(stem_age_i, 1)
  
}

View(all_genus_groups_metadata_df)

# Save Genus-groups metadata
saveRDS(all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/all_genus_groups_metadata_df.rds")

### 4.3/ Extract Genus-groups metadata from UCE only tree ####

# Initiate metadata df for UCE only taxa
all_genus_groups_UCE_only_metadata_df <- data.frame(group_name = genus_groups_list)
all_genus_groups_UCE_only_metadata_df$MRCA_node_ID <- NA
all_genus_groups_UCE_only_metadata_df$MRCA_edge_ID <- NA
all_genus_groups_UCE_only_metadata_df$crown_age <- NA
all_genus_groups_UCE_only_metadata_df$stem_age <- NA
all_genus_groups_UCE_only_metadata_df$sampled_richness <- NA

all_edge_ages <- phytools::nodeHeights(Ponerinae_phylogeny_789t_treedata@phylo)
root_age <- max(all_edge_ages[, 2])
all_edge_ages <- round(-1 * all_edge_ages + root_age, 5)

for (i in 1:nrow(all_genus_groups_UCE_only_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- all_genus_groups_UCE_only_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- genus_groups_MRCA_taxa_list[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_phylogeny_789t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant taxa
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_789t_treedata@phylo, node = MRCA_node_ID_i)
  all_descendants_ID_i <- all_descendants_ID_i[all_descendants_ID_i <= length(Ponerinae_phylogeny_789t_treedata@phylo$tip.label)]
  all_descendants_names_i <- Ponerinae_phylogeny_789t_treedata@phylo$tip.label[all_descendants_ID_i]

  # Record current richness data
  all_genus_groups_UCE_only_metadata_df$sampled_richness[i] <- length(all_descendants_names_i)

  # Extract crown and stem ages
  MRCA_edge_ID_i <- which(Ponerinae_phylogeny_789t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_age_i <- all_edge_ages[MRCA_edge_ID_i, 2]
  stem_age_i <- all_edge_ages[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and ages
  all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  all_genus_groups_UCE_only_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  all_genus_groups_UCE_only_metadata_df$crown_age[i] <- round(crown_age_i, 1)
  all_genus_groups_UCE_only_metadata_df$stem_age[i] <- round(stem_age_i, 1)
  
}

View(all_genus_groups_UCE_only_metadata_df)

# Save Genus-groups metadata
saveRDS(all_genus_groups_UCE_only_metadata_df, file = "./outputs/Grafting_missing_taxa/all_genus_groups_UCE_only_metadata_df.rds")


##### 5/ Prepare edge/node pie color scheme for bioregion membership ####

# ggtree only allows one color per edge
# See phytools::densityMap (in Script 5) for a continuous output
# Compute mean membership per edge per bioregions from Simmap

### 5.1/ Compute mean and final membership per edges ####

# Load treedata for complete phylogeny
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")

# Load density maps with bioregion membership probabilities
DEC_J_density_map_all_areas <- readRDS(file = paste0("./outputs/Density_maps/DEC_J_density_map_all_areas.rds"))

names(DEC_J_density_map_all_areas)
str(DEC_J_density_map_all_areas[[i]])


##### Need to REDO all with split areas, not random attribution !!! #####

# Initiate edge df
egde_bioregion_membership_df <- data.frame(edge_ID = 1:length(DEC_J_density_map_all_areas[[i]]$tree$maps))
mean_egde_bioregion_membership_df <- final_egde_bioregion_membership_df <- egde_bioregion_membership_df

for (i in seq_along(DEC_J_density_map_all_areas))
{
  # i <- 1
  
  # Extract bioregion name
  bioregion_i <- str_remove(string = names(DEC_J_density_map_all_areas)[i], pattern = "Density_map_")
  
  # Extract state maps
  state_maps_i <- DEC_J_density_map_all_areas[[i]]$tree$maps
  
  # Initiate mean probability vector
  mean_prob_per_edge_i <- c()
  final_prob_per_edge_i <- c()
  edge_length <- c()
  
  # Loop per edge
  for (j in seq_along(state_maps_i))
  {
    # j <- 1
    
    # Compute mean membership probability per edge
    states_j <- state_maps_i[[j]]
    binary_states <- as.numeric(names(states_j)) / 1000 # Weighted states
    weights <- states_j # Residence times
    mean_prob_j <- weighted.mean(x = binary_states, w = weights)
    edge_length_j <- sum(weights)
    
    # Extract final probability
    final_prob_j <- binary_states[length(binary_states)]
    
    # Store mean probability
    mean_prob_per_edge_i <- c(mean_prob_per_edge_i, mean_prob_j)
    # Store final probability
    final_prob_per_edge_i <- c(final_prob_per_edge_i, final_prob_j)
    # Store edge length
    edge_length <- c(edge_length, edge_length_j)
    
  }
  
  # Store mean probability in final df
  old_names <- names(mean_egde_bioregion_membership_df)
  mean_egde_bioregion_membership_df <- cbind(mean_egde_bioregion_membership_df, mean_prob_per_edge_i)
  names(mean_egde_bioregion_membership_df) <- c(old_names, paste0("mean_state_",bioregion_i))
  
  # Store mean probability in final df
  old_names <- names(final_egde_bioregion_membership_df)
  final_egde_bioregion_membership_df <- cbind(final_egde_bioregion_membership_df, final_prob_per_edge_i)
  names(final_egde_bioregion_membership_df) <- c(old_names, paste0("final_state_",bioregion_i))
  
  # Add edge length
  egde_bioregion_membership_df$edge_length <- edge_length
  
}

# Merge data for mean and final edge states
egde_bioregion_membership_df <- egde_bioregion_membership_df %>% 
  left_join(mean_egde_bioregion_membership_df) %>% 
  left_join(final_egde_bioregion_membership_df)
  
## Save df of mean membership probability per edges
saveRDS(egde_bioregion_membership_df, file = "./outputs/Density_maps/egde_bioregion_membership_df.rds")


### 5.2/ Deal with the special case of the root ####

# Does not have a parental edge
# Compute the mean of first states of descending edges

# Get ID of descending edges of the root
root_node_ID <- Ponerinae_phylogeny_1534t_treedata@extraInfo$node[Ponerinae_phylogeny_1534t_treedata@extraInfo$status == "root"]
root_descendants_ID <- which(Ponerinae_phylogeny_1534t_treedata@phylo$edge[,1] == root_node_ID)

# Get first states of those descending edges
root_descendants_first_states <- matrix(nrow = length(root_descendants_ID), ncol = length(DEC_J_density_map_all_areas))
colnames(root_descendants_first_states) <- str_remove(string = names(DEC_J_density_map_all_areas), pattern = "Density_map_")
row.names(root_descendants_first_states) <- root_descendants_ID

# Loop per bioregions
for (j in seq_along(DEC_J_density_map_all_areas))
{
  # j <- 1
  
  # Extract state maps
  state_maps_j <- DEC_J_density_map_all_areas[[j]]$tree$maps
  
  for (i in seq_along(root_descendants_ID))
  {
    # i <- 1
    
    # Extract descending edge ID
    edge_i <- root_descendants_ID[i]
    
    # Extract state map
    state_maps_ij <- state_maps_j[[edge_i]]
    # Extract first state probability
    first_prob <- as.numeric(names(state_maps_ij))[length(state_maps_ij)] / 1000
    
    # Store first probability
    root_descendants_first_states[i,j] <- first_prob
  }
}

# Compute mean from starting states of both descending edges
root_states <- apply(X = root_descendants_first_states, MARGIN = 2, FUN = mean)


### 5.3/ Join bioregion membership probability data in treedata object ####

# Check agreement in edge length
plot(egde_bioregion_membership_df$edge_length, Ponerinae_phylogeny_1534t_treedata@phylo$edge.length)

# Both order should be the same
attr(Ponerinae_phylogeny_1534t_treedata@phylo, "order") # Postorder/Cladewise
attr(DEC_J_density_map_all_areas[[1]]$tree, "order") # Cladewise

Ponerinae_phylogeny_1534t_treedata <- Ponerinae_phylogeny_1534t_treedata %>%
  full_join(x = ., y = egde_bioregion_membership_df, by = join_by("parental_edge_ID" == "edge_ID"))

# Add root states info
var_names <- paste0("final_state_", names(root_states))
Ponerinae_phylogeny_1534t_treedata@extraInfo[root_node_ID, var_names] <- as.data.frame(t(root_states))

View(Ponerinae_phylogeny_1534t_treedata@extraInfo)
View(Ponerinae_phylogeny_1534t_treedata@extraInfo[1530:1550,])

# Save treedata for complete phylogeny
saveRDS(Ponerinae_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")


### 5.4/ Prepare node metadata df for likelihood pies (including multiarea ranges) #####

## 5.4.1/ Extract marginal likelihoods from model outputs ####

# Load time-stratified DEC+J model output
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")

# Inspect marginal likelihoods of ancestral states
dim(Ponerinae_phylogeny_1534t_treedata@extraInfo)
dim(DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each branch tipward end = each node/tip, given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start/rootward end of each parental branch of each node/tip, just after cladogenesis/speciation (and eventually cladogenetic transition).
# Values of ARE can be found for the root, but not "corner values" because it does not have a parental edge
DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node[root_node_ID,]
DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node[root_node_ID,]

# Compute list of all possible states
generate_list_ranges <- function (areas_list, max_range_size, include_null_range = TRUE)
{
  states_list_0based = cladoRcpp::rcpp_areas_list_to_states_list(areas = areas_list, maxareas = max_range_size, include_null_range = include_null_range)
  ranges_list <- NULL
  for (i in 1:length(states_list_0based))
  {    
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
    {
      tmprange = "_"
    } else {
      tmprange = paste(areas_list[states_list_0based[[i]]+1], collapse="")
    }
    ranges_list = c(ranges_list, tmprange)
  }
  return(ranges_list)
}

# Update paths with current working directory if needed
DEC_J_fit$inputs$geogfn <- np(paste0(getwd(), "/input_data/BioGeoBEARS_setup/lagrange_area_data_file_7_regions_PaleA.data"))
DEC_J_fit$inputs$trfn <- np(paste0(getwd(), "/outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.tree"))
DEC_J_fit$inputs$timesfn <- np(paste0(getwd(), "/input_data/BioGeoBEARS_setup/time_boundaries.txt"))
  
# Extract full range list
returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
# ranges_list <- returned_mats$ranges_list
reduced_ranges_list = returned_mats$ranges_list
max_range_size <- max(nchar(reduced_ranges_list))
ranges_list <- generate_list_ranges(areas_list = returned_mats$areanames, max_range_size = max_range_size, include_null_range = DEC_J_fit$inputs$include_null_range)

# Extract posterior probabilities based on marginal likelihood of all ranges
node_marginal_likelihoods_df <- DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node
colnames(node_marginal_likelihoods_df) <- ranges_list 

# All probabilities should sum to one 
table(apply(X = node_marginal_likelihoods_df, MARGIN = 1, FUN = sum))

# Remove unused ranges
range_usages <- apply(X = node_marginal_likelihoods_df, MARGIN = 2, FUN = sum)
range_usages <- range_usages / sum(range_usages) * 100
range_usages[order(range_usages, decreasing = T)]

node_marginal_likelihoods_df_only_used_ranges <- node_marginal_likelihoods_df[, range_usages > 0]
dim(node_marginal_likelihoods_df_only_used_ranges) # Only 58 states used
node_marginal_likelihoods_df_reduced_ranges <- node_marginal_likelihoods_df[, range_usages > 0.01]
dim(node_marginal_likelihoods_df_reduced_ranges) # Only 19 states with 0.01% usage

root_states_PP <- node_marginal_likelihoods_df_only_used_ranges[root_node_ID,]
root_states_PP <- node_marginal_likelihoods_df_reduced_ranges[root_node_ID,]
root_states_PP

## Save edge marginal likelihoods df
saveRDS(object = node_marginal_likelihoods_df, file = "./outputs/Ancestral_range_estimates_maps/DEC_J_node_marginal_likelihoods_df.rds")
saveRDS(object = node_marginal_likelihoods_df_only_used_ranges, file = "./outputs/Ancestral_range_estimates_maps/DEC_J_node_marginal_likelihoods_df_only_used_ranges.rds")
saveRDS(object = node_marginal_likelihoods_df_reduced_ranges, file = "./outputs/Ancestral_range_estimates_maps/DEC_J_node_marginal_likelihoods_df_reduced_ranges.rds")

## 5.4.2/ Create ggplot pies ####

## Load edge marginal likelihoods df
node_marginal_likelihoods_df <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/DEC_J_node_marginal_likelihoods_df.rds")
node_marginal_likelihoods_df_reduced_ranges <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/DEC_J_node_marginal_likelihoods_df_reduced_ranges.rds")

# Load color scheme
ranges_list <- colnames(node_marginal_likelihoods_df)
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
reduced_ranges_list <- colnames(node_marginal_likelihoods_df_reduced_ranges) 
colors_list_for_reduced_ranges <- colors_list_for_states[reduced_ranges_list]

# Melt edge marginal likelihood df
node_marginal_likelihoods_df <- as.data.frame(node_marginal_likelihoods_df)
node_marginal_likelihoods_df$node_ID <- 1:nrow(node_marginal_likelihoods_df)
node_marginal_likelihoods_df_melted <- reshape2::melt(node_marginal_likelihoods_df, id.var = "node_ID")
names(node_marginal_likelihoods_df_melted) <- c("node_ID", "range", "PP")

# # Melt edge marginal likelihood df for reduced ranges
# node_marginal_likelihoods_df <- as.data.frame(node_marginal_likelihoods_df_reduced_ranges)
# node_marginal_likelihoods_df_reduced_ranges$node_ID <- 1:nrow(node_marginal_likelihoods_df_reduced_ranges)
# node_marginal_likelihoods_df_reduced_ranges_melted <- reshape2::melt(node_marginal_likelihoods_df_reduced_ranges, id.var = "node_ID")
# names(node_marginal_likelihoods_df_reduced_ranges_melted) <- c("node_ID", "range", "PP")


## Create list of pie charts
ARE_pies_list <- list()
for (i in 1:nrow(node_marginal_likelihoods_df))
# for (i in 1:nrow(node_marginal_likelihoods_df_reduced_ranges)) 
    
{
  # i <- 2155
  
  # Extract data for a given node
  ARE_data_i <- node_marginal_likelihoods_df_melted %>% 
  # ARE_data_i <- node_marginal_likelihoods_df_reduced_ranges_melted %>% 
    filter(node_ID == i) %>% 
    select(range, PP)
  
  # Create a ggplot object for each pie chart
  ARE_pie_i <- ggplot(data = ARE_data_i,
                      aes(y = PP, fill = range, x = "")) + 
    geom_bar(stat = "identity", color = "black", show.legend = F) +
    coord_polar("y", start = 0) +
    
    scale_fill_manual("Ranges", breaks = ranges_list, values = colors_list_for_states, guide = F) +
    # scale_fill_manual("Ranges", breaks = reduced_ranges_list, values = colors_list_for_reduced_ranges, guide = F) +
    
    theme_inset()
    # theme_void()
    
  # plot(ARE_pie_i)  
    
  ARE_pies_list[[i]] <- ARE_pie_i
  
  ## Print progress
  if (i %% 1000 == 0)
  {
    cat(paste0(Sys.time(), " - Created pie chart for edge n°", i, "/", nrow(node_marginal_likelihoods_df),"\n"))
  }
}

# Assign edge ID to pies
names(ARE_pies_list) <- 1:nrow(node_marginal_likelihoods_df)
# names(ARE_pies_list) <- 1:nrow(node_marginal_likelihoods_df_reduced_ranges)

## Save ARE pie lists
saveRDS(ARE_pies_list, file = "./outputs/Ancestral_range_estimates_maps/DEC_J_ARE_pies_list.rds")

print(ARE_pies_list[[root_node_ID]])
print(ARE_pies_list[[3000]])

## Make list of pie charts with ggtree

?ggtree::nodepie

test <- node_marginal_likelihoods_df
test <- test %>% 
  rename(node = node_ID)
ARE_pies_ggtree <- ggtree::nodepie(data = test,
                                   cols = 1:(ncol(test)-1),
                                   color = colors_list_for_states,
                                   alpha = 1.0,
                                   outline.color = "black",
                                   outline.size = 0.2)

## Test plots

# Ponerinae_phylogeny_plot <- ggtree::inset(tree_view = Ponerinae_phylogeny_plot,
#                                           insets = ARE_pies_list[3000:3060],
#                                           # insets = ARE_pies_ggtree[3000:3100],
#                                           hjust = 0, 
#                                           width = 0.1, height = 0.1)
# 
# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   geom_inset(insets = ARE_pies_list[3000:3060],
#              # insets = ARE_pies_ggtree[3000:3100],
#              width = 0.1, height = 0.1,
#              hjust = 0, vjust = 0,
#              x = "node")

### 5.5/ Plot color scheme of ranges based on root states ####

# Extract data for a given node
ARE_data_i <- node_marginal_likelihoods_df_melted %>% 
  # ARE_data_i <- node_marginal_likelihoods_df_reduced_ranges_melted %>% 
  filter(node_ID == i) %>% 
  select(range, PP)

# Get ID of the root node
root_node_ID <- Ponerinae_phylogeny_1534t_treedata@extraInfo$node[Ponerinae_phylogeny_1534t_treedata@extraInfo$status == "root"]
root_node_ID

# Extract data for a root node
ARE_data_root_node <- node_marginal_likelihoods_df_melted %>% 
  filter(node_ID == root_node_ID) %>% 
  select(range, PP) %>% 
  mutate(PP = round(PP, 3) * 100) %>% 
  filter(PP > 1)

# Adjust color scheme for root states
root_states_reduced <- ARE_data_root_node$range
root_states_reduced <- c("I", "AI", "A", "AN", "N", "AIN") # Manual rearrangement
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_root_reduced_ranges <- colors_list_for_states[root_states_reduced]

# Adjust order of states
ARE_data_root_node$range <- factor(ARE_data_root_node$range, levels = root_states_reduced, labels = root_states_reduced)


# Create a pie chart for the root data
ARE_pie_root <- ggplot(data = ARE_data_root_node,
                    aes(y = PP, fill = range, x = "")) + 
  geom_bar(stat = "identity", color = "black", show.legend = T) +
  coord_polar("y", direction = -1, start = 0) +
  
  scale_fill_manual("Ranges", breaks = root_states_reduced, values = colors_list_for_root_reduced_ranges) +
  
  theme_void() +
  theme(plot.margin = unit(c(0.0,0.1,0.0,0), "npc"), # trbl
        legend.title = element_text(size  = 20, face = "bold", margin = margin(b = 15)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"))

# Export backbone ARE color scheme

pdf(file = "./outputs/Ancestral_range_estimates_maps/ARE_pie_root.pdf", width = 5, height = 5)

print(ARE_pie_root)

dev.off()


### 5.6/ Arrange state/color order in backbone pies as in root ARE ####

## 5.6.1/ Subset the nodes occurring before the genus-groups MRCA = backbone ####

# Get all descendant nodes of genus-groups
all_descendant_nodes <- c()
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  # Get descendant nodes for genus-group i
  MRCA_node_ID_i <- all_genus_groups_metadata_df$MRCA_node_ID[i]
  all_descendant_nodes_i <- getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                           node = MRCA_node_ID_i)
  # Store descendants
  all_descendant_nodes <- c(all_descendant_nodes, all_descendant_nodes_i)
}

# Extract backbone nodes by difference
backbone_nodes_ID <- setdiff(1:length(ARE_pies_list), all_descendant_nodes)


## Create list of pie charts
backbone_ARE_pies_list <- list()
for (i in seq_along(backbone_nodes_ID))
{
  # i <- 1
  
  # Extract node ID
  node_ID_i <- backbone_nodes_ID[i]
  
  # Extract data for a given node
  ARE_data_i <- node_marginal_likelihoods_df_melted %>% 
    # ARE_data_i <- node_marginal_likelihoods_df_reduced_ranges_melted %>% 
    filter(node_ID == node_ID_i) %>% 
    filter(range %in% root_states_reduced) %>%
    select(range, PP)
  
  # Adjust order of states
  ARE_data_i$range <- factor(ARE_data_i$range, levels = root_states_reduced, labels = root_states_reduced)
    
  # Create a ggplot object for each pie chart
  ARE_pie_i <- ggplot(data = ARE_data_i,
                      aes(y = PP, fill = range, x = "")) + 
    geom_bar(stat = "identity", color = "black", show.legend = F) +
    coord_polar("y", start = 0) +
    
    scale_fill_manual("Ranges", breaks = root_states_reduced, values = colors_list_for_root_reduced_ranges) +

    theme_inset()
  # theme_void()
  
  # print(ARE_pie_i)  
  
  backbone_ARE_pies_list[[i]] <- ARE_pie_i
  
}

# Assign backbone node ID to pies
names(backbone_ARE_pies_list) <- backbone_nodes_ID

## Save backbone ARE pie lists
saveRDS(backbone_ARE_pies_list, file = "./outputs/Ancestral_range_estimates_maps/DEC_J_backbone_ARE_pies_list.rds")



##### 6/ Plot complete phylogenies ####

# Load Genus-groups metadata
all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/all_genus_groups_metadata_df.rds")

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]

# Set color scheme for edges (Not used because use aes mapping to get a legend)
branches_color <- c("black", "red")[Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node + 1]

# Set color scheme for tips
nb_tips <- length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)
tips_color <- c("black", "red")[Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node + 1]
tips_color <- tips_color[1:nb_tips]

## 6.1/ Rectangular with all nod/branch labels ####

pdf(file = "./outputs/Grafting_missing_taxa/Full_phylo_rectangular_all_labels.pdf", height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5) +
  
  geom_tiplab(# mapping = aes(colour = missing_node),
              color = tips_color,
              align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 2, 
              offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
              hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
                                 
  # Add ID labels to nodes
  geom_label(aes(label = node),
             fill = 'steelblue', size = 1) +
  
  # Add ID labels to edges
  geom_label(aes(x = branch, label = parental_edge_ID),
             fill = 'lightgreen', size = 1) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.002,0.05,0.002,0), "npc"), # trbl
        legend.position = c(1.0, 0.5),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
        )
  

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 7,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = 5.0,
                    offset.text = 0.4,
                    barsize = 4)    # Width of the bar
}

print(Ponerinae_phylogeny_plot)

dev.off()

# Look at the circular version
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  layout_circular()
print(Ponerinae_phylogeny_plot)


## 6.2/ Circular layout with no edge/node and tip labels but genus-group labels ####

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/Full_phylo_circular_all_genus_groups_no_tiplabels.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "circular") +

  geom_tree(mapping = aes(colour = missing_node),
          # color = branches_color,
          linewidth = 0.5) +
  
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   # color = tips_color,
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +

  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )
  
# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = genus_groups_angles[i],
                    hjust = 'center',
                    offset = 8.0,
                    offset.text = genus_groups_offsets[i],
                    barsize = 25)    # Width of the bar
}

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
        plot.margin = unit(c(-0.08,-0.2,-0.1,-0.25), "npc"), # trbl
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(color = "transparent", fill = NA),
        panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()



## 6.3/ Circular layout with no edge/node but tip labels and genus-group labels ####

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/Full_phylo_circular_all_genus_groups_with_tiplabels.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "circular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5) +
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = genus_groups_angles[i],
                    hjust = 'center',
                    offset = 5.5,
                    offset.text = genus_groups_offsets[i],
                    barsize = 25)    # Width of the bar
}


# Add tip labels
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  geom_tiplab(# mapping = aes(colour = missing_node),
    # color = tips_color,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0,
    size = 0.8,
    offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0)      # To center, left-align, right-aligned text
  
# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(-0.08,-0.2,-0.1,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()


### 6.4/ Circular layout of edge bioregion membership: Per bioregions ####

# Loop per bioregion on different PDFs
# Make color gradient from focal color to grey80

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
areas_list <- c("A", "U", "E", "I", "R", "N", "W")
colors_list_for_areas <- colors_list_for_states[areas_list]
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

## 6.4.1/ Loop per bioregions ####
for (i in seq_along(bioregion_names))
{
  # i <- 1
  
  # Extract bioregion
  bioregion_i <- bioregion_names[i]
  
  # Find index in metadata df
  bioregion_ID_i <- which(names(Ponerinae_phylogeny_1534t_treedata@extraInfo) == bioregion_i)
  
  # Plot PDF for bioregion i
  pdf(file = paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_",bioregion_i,".pdf"), height = 16, width = 17)
  
  Ponerinae_phylogeny_plot_i <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                       # color = NA,
                                       color = "grey90",
                                       layout = "circular") +
    
    geom_tree(color = colors_list_for_areas[i],
              alpha = Ponerinae_phylogeny_1534t_treedata@extraInfo[, bioregion_ID_i, drop = T],
              linewidth = 0.5) +
  
  # Add title
  ggtitle(label = paste0("Posterior probability of range including ", bioregion_i)) 
  
  # Add Genus-group cladelabels
  for (j in 1:nrow(all_genus_groups_metadata_df))
  {
    Ponerinae_phylogeny_plot_i <- Ponerinae_phylogeny_plot_i +
      
      # Add a bar alongside all the clade tips. Add a label to the median tip.
      geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[j],
                      # label = all_genus_groups_metadata_df$group_name[j], 
                      label = paste0('bold(',all_genus_groups_metadata_df$group_name[j],')'), parse = T,
                      align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                      geom = 'text',     # Tip of label, with or without rectangle
                      fontsize = 10,      # Size of the text
                      fill = NA,         # Color of the label background
                      color = c(genus_groups_colors[j], genus_groups_colors[j]),   # Color of the bar and the text
                      angle = genus_groups_angles[j],
                      hjust = 'center',
                      offset = 8.0,
                      offset.text = genus_groups_offsets[j],
                      barsize = 25)    # Width of the bar
  }
  
  # Add age lines
  age_step <- 20
  age <- max(Ponerinae_phylogeny_plot_i$data$x) - age_step
  while (age > 0)
  { 
    Ponerinae_phylogeny_plot_i <- Ponerinae_phylogeny_plot_i +
      geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
    age  <- age - age_step
  }
 
  # Adapt margins
  Ponerinae_phylogeny_plot_i <- Ponerinae_phylogeny_plot_i +
    theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
          plot.title = element_text(face = "bold", colour = colors_list_for_areas[i], size = 30, hjust = 0.5,
                                    angle = 0, vjust = -10),
          plot.margin = unit(c(-0.05,-0.2,-0.1,-0.25), "npc"), # trbl
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(color = "transparent", fill = NA),
          panel.background = element_rect(color = NA, fill = "transparent"))
  
  print(Ponerinae_phylogeny_plot_i)
  
  dev.off()
  
  
}

## 6.4.2/ Make a single multipages PDF ####

# Recreate paths to density maps in predefined order
all_density_maps_paths <- paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_",bioregion_names,".pdf")

# Combine PDFs in a unique PDF
qpdf::pdf_combine(input = all_density_maps_paths, output = paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_all_bioregions_ordered.pdf"))


## 6.4.3/ Make a GIF ####

source("./functions/image_resize_and_write_gif.R")

pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_all_bioregions_ordered.pdf"),
                                      pages = NULL, density = 100)
magick::image_info(pdf_pointer)

image_resize_and_write_gif(image = pdf_pointer,
                           path =  paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_all_bioregions_ordered.gif"),
                           delay = 1, # Time between frames in seconds
                           width = 1700, height = 1600,
                           loop = FALSE,
                           progress = TRUE)


### 6.5/ Circular layout of edge bioregion membership: Aggregated ####

# Loop over bioregion to overlay colors on a single PDF
# Make an alpha gradient, not a color gradient !

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
areas_list <- c("A", "U", "E", "I", "R", "N", "W")
colors_list_for_areas <- colors_list_for_states[areas_list]
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name


## 6.5.1/ Generate PDF for all bioregions aggregated ####

# Plot PDF to aggregate all bioregions
pdf(file = paste0("./outputs/Density_maps/Full_phylo_circular_bioregion_membership_all_bioregions_aggregated.pdf"), height = 16, width = 17)

# Initiate plot
Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   color = NA,
                                   # color = "grey90",
                                   layout = "circular") +
  # Add title
  ggtitle(label = paste0("Posterior probability of lineage ranges")) 

## Loop per bioregions
for (i in seq_along(bioregion_names))
{
  # i <- 1
  
  # Extract bioregion
  bioregion_i <- bioregion_names[i]
  
  # Find index in metadata df
  bioregion_ID_i <- which(names(Ponerinae_phylogeny_1534t_treedata@extraInfo) == bioregion_i)
  
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
    geom_tree(color = colors_list_for_areas[i],
              alpha = Ponerinae_phylogeny_1534t_treedata@extraInfo[, bioregion_ID_i, drop = T],
              linewidth = 0.5)
}

## Add Genus-group cladelabels
for (j in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[j],
                    # label = all_genus_groups_metadata_df$group_name[j], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[j],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[j], genus_groups_colors[j]),   # Color of the bar and the text
                    angle = genus_groups_angles[j],
                    hjust = 'center',
                    offset = 8.0,
                    offset.text = genus_groups_offsets[j],
                    barsize = 25)    # Width of the bar
}

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.title = element_text(face = "bold", colour = "black", size = 30, hjust = 0.5,
                              angle = 0, vjust = -10),
    plot.margin = unit(c(-0.05,-0.2,-0.1,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()


### 6.6/ Plot rectangular layout with ARE pie charts ####

pdf(file = "./outputs/Ancestral_range_estimates_maps/Full_phylo_rectangular_ARE_pie_charts.pdf", height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5) +
  
  geom_tiplab(# mapping = aes(colour = missing_node),
    color = tips_color,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.002,0.05,0.002,0), "npc"), # trbl
        legend.position = c(1.0, 0.5),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )


# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 7,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = 5.0,
                    offset.text = 0.4,
                    barsize = 4)    # Width of the bar
}

## Add ARE pie charts
Ponerinae_phylogeny_plot <- ggtree::inset(tree_view = Ponerinae_phylogeny_plot,
                                          insets = ARE_pies_list,
                                          # insets = ARE_pies_ggtree[3000:3100],
                                          hjust = 0, 
                                          width = 0.01, height = 0.01)

# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   geom_inset(insets = ARE_pies_list[3000:3060],
#              # insets = ARE_pies_ggtree[3000:3100],
#              width = 0.1, height = 0.1,
#              hjust = 0, vjust = 0,
#              x = "node")

# print(Ponerinae_phylogeny_plot)

dev.off()



### 6.7/ Plot circular layout with ARE pie charts ####

## Load ARE pie lists
ARE_pies_list <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/DEC_J_ARE_pies_list.rds")

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

### Add reduced legend of fill color scheme including multiarea ranges #####


# Plot PDF
pdf(file = "./outputs/Ancestral_range_estimates_maps/Full_phylo_circular_ARE_pie_charts.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "circular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5) +
  
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   # color = tips_color,
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = genus_groups_angles[i],
                    hjust = 'center',
                    offset = 8.0,
                    offset.text = genus_groups_offsets[i],
                    barsize = 25)    # Width of the bar
}

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(-0.08,-0.2,-0.1,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

## Add ARE pie charts
Ponerinae_phylogeny_plot <- ggtree::inset(tree_view = Ponerinae_phylogeny_plot,
                                          insets = ARE_pies_list[3000:3060],
                                          # insets = ARE_pies_ggtree[3000:3100],
                                          hjust = 0, 
                                          width = 0.1, height = 0.1)

# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   geom_inset(insets = ARE_pies_list[3000:3060],
#              # insets = ARE_pies_ggtree[3000:3100],
#              width = 0.1, height = 0.1,
#              hjust = 0, vjust = 0,
#              x = "node")

print(Ponerinae_phylogeny_plot)

dev.off()


# Does not work on circular tree
# See fix in https://github.com/YuLab-SMU/ggtree/issues/419


##### 7/ Plot UCE data only phylogenies #####

Ponerinae_phylogeny_789t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_789t_treedata.rds")

# To flip position of Ponera and Pachycondyla
Ponera_MRCA_node_ID <- all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[all_genus_groups_UCE_only_metadata_df$group_name == "Ponera"]
Pachycondyla_MRCA_node_ID <- all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[all_genus_groups_UCE_only_metadata_df$group_name == "Pachycondyla"]

## 7.1/ Rectangular with all nod/branch labels ####

# Set genus-group offset to avoid overlap
genus_groups_offsets <- c(0.4, 0.4, 0.4, 2.0, 0.4, 0.4, 0.4)
names(genus_groups_offsets) <- all_genus_groups_UCE_only_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/UCE_only_phylo_rectangular_all_labels.pdf", height = 50, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_789t_treedata,
                                   layout = "rectangular") +
  
  geom_tree(linewidth = 0.5) +
  
  geom_tiplab(color = "black",
              align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 2, 
              offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
              hjust = 0) +    # To center, left-align, right-aligned text
  
  # Add ID labels to nodes
  geom_label(aes(label = node),
             fill = 'steelblue', size = 1) +
  
  # Add ID labels to edges
  geom_label(aes(x = branch, label = node),
             fill = 'lightgreen', size = 1) +
  
  # Adjust margins
  theme(plot.margin = unit(c(0.002,0,0.002,0), "npc"), # trbl
       )

# Flip position of Ponera and Pachycondyla
Ponerinae_phylogeny_plot <- ggtree::flip(tree_view = Ponerinae_phylogeny_plot, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)
  
# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_UCE_only_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_UCE_only_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_UCE_only_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 8,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = 5.0,
                    offset.text = genus_groups_offsets[i],
                    barsize = 4)    # Width of the bar
}

print(Ponerinae_phylogeny_plot)

dev.off()

# Look at the circular version
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  layout_circular()
print(Ponerinae_phylogeny_plot)


## 7.2/ Circular layout with no edge/node and tip labels but genus-group labels ####

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 310, 358, 15, 58, 21, 2)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_UCE_only_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_UCE_only_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/UCE_only_phylo_circular_all_genus_groups_no_tiplabels.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_789t_treedata,
                                   layout = "circular") +
  
  geom_tree(linewidth = 0.5)
  
  # geom_tiplab(
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
# Flip position of Ponera and Pachycondyla
Ponerinae_phylogeny_plot <- ggtree::flip(tree_view = Ponerinae_phylogeny_plot, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_UCE_only_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_UCE_only_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_UCE_only_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,     # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    # angle = 0,
                    angle = genus_groups_angles[i],
                    hjust = 'center',
                    offset = 8.0,
                    # offset.text = 0.2,
                    offset.text = genus_groups_offsets[i],
                    barsize = 25)    # Width of the bar
}

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(-0.08,-0.25,-0.10,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()



## 7.3/ Circular layout with no edge/node but tip labels and genus-group labels ####

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 310, 358, 15, 58, 21, 2)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_UCE_only_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_UCE_only_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/UCE_only_circular_all_genus_groups_with_tiplabels.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_789t_treedata,
                                   layout = "circular") +
  
  geom_tree(linewidth = 0.5)
  
# Flip position of Ponera and Pachycondyla
Ponerinae_phylogeny_plot <- ggtree::flip(tree_view = Ponerinae_phylogeny_plot, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_UCE_only_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a bar alongside all the clade tips. Add a label to the median tip.
    geom_cladelabel(node = all_genus_groups_UCE_only_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_UCE_only_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 10,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = genus_groups_angles[i],
                    hjust = 'center',
                    offset = 5.5,
                    offset.text = genus_groups_offsets[i],
                    barsize = 25)    # Width of the bar
}


# Add tip labels
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  geom_tiplab(align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0,
              size = 0.8,
              offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
              hjust = 0)      # To center, left-align, right-aligned text

# Add age lines
age_step <- 20
age <- max(Ponerinae_phylogeny_plot$data$x) - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(-0.08,-0.25,-0.10,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()






##### 8/ Plot phylogenies collapsed to Genus-groups #####

?ggtree:collapse

### 8.1/ Plot circular phylogeny with genus-groups labels ####

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

pdf(file = "./outputs/Grafting_missing_taxa/Genus_group_phylo_circular_no_tiplabels.pdf", height = 16, width = 17)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "circular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5, show.legend = F) +
  
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   # color = tips_color,
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )

# Collapse Genus-groups
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- ggtree::collapse(x = Ponerinae_phylogeny_plot,
                                               node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                                               mode = "max", 
                                               clade_name = all_genus_groups_metadata_df$group_name[i],
                                               alpha = 1.0, 
                                               color = "black",
                                               fill = genus_groups_colors[i])
  
}

# Add age lines
age_step <- 20
age <- 113 - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(-0.08,-0.2,-0.1,-0.25), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()


### 8.2/ Plot circular phylogeny with ARE pies and genus-groups labels ####

# Does not work on circular tree
# See fix in https://github.com/YuLab-SMU/ggtree/issues/419


### 8.3/ Plot rectangular phylogeny with genus-groups labels ####

pdf(file = "./outputs/Grafting_missing_taxa/Genus_group_phylo_rectangular_no_tiplabels.pdf",
    height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5, show.legend = F) +
  
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   # color = tips_color,
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )


# Collapse Genus-groups
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- ggtree::collapse(x = Ponerinae_phylogeny_plot,
                                               node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                                               mode = "max", 
                                               clade_name = all_genus_groups_metadata_df$group_name[i],
                                               alpha = 1.0, 
                                               color = "black",
                                               fill = genus_groups_colors[i])
  
}

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 15,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = min(all_genus_groups_metadata_df$crown_age) + 3.0,
                    offset.text = 0.4,
                    barsize = 0)    # Width of the bar
}



# Add age lines
age_step <- 20
age <- 113 - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

# # Adapt margins
# Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
#   theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
#     plot.margin = unit(c(-0.08,-0.2,-0.1,-0.25), "npc"), # trbl
#     plot.background = element_rect(fill = "transparent",colour = NA),
#     panel.border = element_rect(color = "transparent", fill = NA),
#     panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()


### 8.4/ Plot rectangular phylogeny with ARE pie charts and genus-groups labels ####

## Load ARE pie lists
ARE_pies_list <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/DEC_J_ARE_pies_list.rds")

## Load backbone ARE pie lists
backbone_ARE_pies_list <- readRDS(file = "./outputs/Ancestral_range_estimates_maps/DEC_J_backbone_ARE_pies_list.rds")


## 8.4.1/ Subset the nodes occurring before the genus-groups MRCA = backbone ####

# Get all descendant nodes of genus-groups
all_descendant_nodes <- c()
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  # Get descendant nodes for genus-group i
  MRCA_node_ID_i <- all_genus_groups_metadata_df$MRCA_node_ID[i]
  all_descendant_nodes_i <- getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                           node = MRCA_node_ID_i)
  # Store descendants
  all_descendant_nodes <- c(all_descendant_nodes, all_descendant_nodes_i)
}

# Extract backbone nodes by difference
backbone_nodes_ID <- setdiff(1:length(ARE_pies_list), all_descendant_nodes)

## 8.4.2/ Plot with backbone ARE pie charts ####

pdf(file = "./outputs/Ancestral_range_estimates_maps/Genus_group_phylo_rectangular_ARE_pie_charts.pdf",
    height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  geom_tree(mapping = aes(colour = missing_node),
            # color = branches_color,
            linewidth = 0.5, show.legend = F) +
  
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   # color = tips_color,
  #   align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0,
  #   size = 2,
  #   offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Adjust legend for taxa source
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )


# Collapse Genus-groups
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- ggtree::collapse(x = Ponerinae_phylogeny_plot,
                                               node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                                               mode = "max", 
                                               clade_name = all_genus_groups_metadata_df$group_name[i],
                                               alpha = 1.0, 
                                               color = "black",
                                               fill = genus_groups_colors[i])
  
}

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 15,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    offset = min(all_genus_groups_metadata_df$crown_age) + 3.0,
                    offset.text = 0.4,
                    barsize = 0)    # Width of the bar
}

# Add age lines
age_step <- 20
age <- 113 - age_step
while (age > 0)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, linetype = 92, linewidth = 1, color = 'grey80')
  age  <- age - age_step
}

## Add ARE pie charts
# Ponerinae_phylogeny_plot <- ggtree::inset(tree_view = Ponerinae_phylogeny_plot,
#                                           insets = ARE_pies_list[backbone_nodes_ID],
#                                           # insets = ARE_pies_ggtree[backbone_nodes_ID],
#                                           hjust = 0, 
#                                           width = 0.2, height = 0.2)

Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  geom_inset(# insets = ARE_pies_list[backbone_nodes_ID],
             insets = backbone_ARE_pies_list,
             width = 0.2, height = 0.2,
             hjust = 0, vjust = 0,
             x = "node")

# Adapt margins
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  theme(# plot.margin = unit(c(-80, -20, -100, -20), "mm"),
    plot.margin = unit(c(0.0,0.0,0.02,0.0), "npc"), # trbl
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(color = "transparent", fill = NA),
    panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_plot)

dev.off()


### Make a zoom in to show the ARE of the genus-groups ???


