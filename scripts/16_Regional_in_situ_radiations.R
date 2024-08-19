##### Script 16: Investigate regional in situ radiations  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Define ‘Regional In Situ Radiations’ (RISR) as independent clades including at least X% of residence times along descendant edges within the focal region
# Fit BD time-dependent models (RPANDA) per RISR as an exponential function of time and fixed extinction rate
  # Extract current richenss, crown-age, net diversification rate at the crown (lambda_0), net diversification rate at present, time-variation (alpha).
# Explore frequency of EB fits: evidence for burst of diversification following disperal to new bioregion
  # Compare fits of time-dependent (decaying) vs. constant BD
  # Count relative frequency of decaying model selection as the best fit
# Compare Bioregions, OW vs. NW, Tropics vs. Temperate
  # Compare crown-ages ~ Bioregions: evidence for “Time-for-accumulation hypothesis”
  # Compare net_div_0 ~ Bioregions => evidence for "Niche opportunity"
  # Compare current net_div ~ Bioregions => current source/sinks of biodiversity
# Test which parameter explains the best the current richness: multi-factor ANOVA with model selection based on AICc
  # If crown age: "Time-for-accumulation hypothesis"
  # If lambda_0: Cradle hypothesis (higher initial net div?)
  # If alpha: Museum hypothesis? (less decline in rates with time?)

###

### Inputs

# (Set of) Time-calibrated phylogeny(ies)
# Probability of bioregion membership per branch

###

### Outputs

# BD time-dependent models per RISR
# EB vs. constant BD model fits per RISR
  # Count relative frequency of EB selection as the best fit
# Plot comparing Bioregions, OW vs. NW, Tropics vs. Temperate
  # Plot current richness ~ Bioregions
  # Plot crown-ages ~ Bioregions: evidence for “Time-for-accumulation hypothesis”
  # Plot net_div_0 ~ Bioregions => evidence for "Niche opportunity"
  # Plot alpha ~ Bioregions
  # Plot current net_div ~ Bioregions => current source/sinks of biodiversity
  # Plot current richness ~ crown age x Bioregions (fit regression lines per Bioregions)
# Current richness ~ Multi-factor ANOVA with model selection based on AICc
  # If crown age: "Time-for-accumulation hypothesis"
  # If lambda_0: Cradle hypothesis (higher initial net div?)
  # If alpha: Museum hypothesis? (less decline in rates with time?)

###

### Bonus in other scripts (?)

# For a more smooth model with many small cladogenetic changes, see CLaDS
# GeoSSE, HiSSE, BiSSE for Old World vs. New world

# Clean environment
rm(list = ls())

library(ape)
library(phytools)
library(BAMMtools)
library(tidyverse)
library(diversitree)  # For time-dependent BD models
library(ggpubr) # For paiwise significance bars on ggplots
library(reghelper)


##### 1/ Load input files ####

### 1.1/ Load the phylogeny ####

# Load phylogeny
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")

### 1.2/ Prepare node metadata ####

# Get node ages per branch (no root edge)
all_edges_df <- phytools::nodeHeights(Ponerinae_phylogeny_1534t_treedata@phylo)
root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_treedata@phylo)[,2])
all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
all_edges_df$edge_ID <- row.names(all_edges_df)
# Join node_ID
all_edges_ID_df <- Ponerinae_phylogeny_1534t_treedata@phylo$edge
colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
all_edges_df <- cbind(all_edges_df, all_edges_ID_df) %>% 
  dplyr::select("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age") %>% 
  arrange(tipward_node_ID)

# Extract binary bioregion membership per node
RISR_nodes_metadata_df <- Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("node", "status", "edge_length", "Bioregion",
                                                                         "mean_state_Neotropics", "mean_state_Afrotropics", "mean_state_Eastern Palearctic",
                                                                         "mean_state_Indomalaya", "mean_state_Western Palearctic", "mean_state_Australasia", "mean_state_Nearctic")]  %>%
  rename(PP_Afrotropics = mean_state_Afrotropics,
         PP_Australasia = mean_state_Australasia,
         PP_Indomalaya = mean_state_Indomalaya,
         PP_Nearctic = mean_state_Nearctic,
         PP_Neotropics = mean_state_Neotropics,
         PP_Eastern_Palearctic = `mean_state_Eastern Palearctic`,
         PP_Western_Palearctic = `mean_state_Western Palearctic`)

# Assign node ages
RISR_nodes_metadata_df <- left_join(x = RISR_nodes_metadata_df, y = all_edges_df[c("tipward_node_ID", "tipward_node_age")], by = join_by(node == tipward_node_ID)) %>% 
  rename(node_age = tipward_node_age)
RISR_nodes_metadata_df$node_age[RISR_nodes_metadata_df$status == "root"] <- root_age

# Assign Old World vs New Word binary and PP status 
RISR_nodes_metadata_df$OW_vs_NW <- "Old_World"
RISR_nodes_metadata_df$OW_vs_NW[RISR_nodes_metadata_df$Bioregion %in% c("Neotropics", "Nearctic")] <- "New_World"
RISR_nodes_metadata_df$PP_New_World <- RISR_nodes_metadata_df$PP_Neotropics + RISR_nodes_metadata_df$PP_Nearctic
RISR_nodes_metadata_df$PP_Old_World <- 1 - RISR_nodes_metadata_df$PP_New_World

# Assign Tropics vs Temperate binary and PP status 
RISR_nodes_metadata_df$Trop_vs_Temp <- "Tropics"
RISR_nodes_metadata_df$Trop_vs_Temp[RISR_nodes_metadata_df$Bioregion %in% c("Nearctic", "Eastern Palearctic", "Western Palearctic")] <- "Temperate"
RISR_nodes_metadata_df$PP_Temperate <- RISR_nodes_metadata_df$PP_Nearctic + RISR_nodes_metadata_df$PP_Eastern_Palearctic + RISR_nodes_metadata_df$PP_Western_Palearctic
RISR_nodes_metadata_df$PP_Tropics <- 1 - RISR_nodes_metadata_df$PP_Temperate

# Save nodes metadata for RISR
saveRDS(RISR_nodes_metadata_df, "./outputs/RISR/RISR_nodes_metadata_df.rds")

### 1.3/ AICc function for model fits ####

# Home-made function to compute AICc from AIC, number of observations, and number of parameters (k)
compute_AICc <- function (AIC, nobs, k)
{
  AICc <- AIC + (2 * k) * (k - 1) / (nobs - k - 1)
}


##### 2/ Identify ‘Regional In Situ Radiations’ (RISR) ####

# Rules to identify RISR per Bioregions
  # Compute proportion of residence times along descendant edges in focal node Bioregion for all nodes
      # Use binary membership for focal node
      # Use PP membership for descending nodes
  # Identify nodes with at least X% of endemicity of residence time
  # Remove clades with less than 5 tips
  # Select RISR in stepwise fashion
    # Select nodes from oldest to most recent
    # Remove all descending nodes when a node is selected as RISR
  # Plot RISR per Bioregions: Highlight nodes and color descending edges


### 2.1/ Compute proportion of endemicity in descendant nodes ####

## Loop per nodes
RISR_nodes_metadata_df$nb_descendant_tips <- NA
RISR_nodes_metadata_df$residence_time_desc_sum_Bioregion <- NA
RISR_nodes_metadata_df$residence_time_desc_prop_Bioregion <- NA
RISR_nodes_metadata_df$residence_time_desc_sum_OW_NW <- NA
RISR_nodes_metadata_df$residence_time_desc_prop_OW_NW <- NA
RISR_nodes_metadata_df$residence_time_desc_sum_Trop_Temp <- NA
RISR_nodes_metadata_df$residence_time_desc_prop_Trop_Temp <- NA

for (i in 1:nrow(RISR_nodes_metadata_df))
{
  # i <- 3000
  
  # Extract tipward node ID
  node_ID_i <- RISR_nodes_metadata_df$node[i]
  
  # Identify focal regions
  bioregion_i <- RISR_nodes_metadata_df$Bioregion[i]
  OW_NW_i <- RISR_nodes_metadata_df$OW_vs_NW[i]
  Trop_Temp_i <- RISR_nodes_metadata_df$Trop_vs_Temp[i]
  
  # Extract all descendant nodes
  descendant_nodes_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                 node = node_ID_i)
  descendant_nodes_i <- setdiff(descendant_nodes_i, node_ID_i) # Remove focal node itself
  
  # Extract tips
  descendant_tips_nb_i <- sum(RISR_nodes_metadata_df$status[RISR_nodes_metadata_df$node %in% descendant_nodes_i] == "tip")
  # Store nb of descendant tips
  RISR_nodes_metadata_df$nb_descendant_tips[i] <- descendant_tips_nb_i
  
  if (length(descendant_nodes_i) > 0)
  {
    # Identify descendant nodes rows
    descendant_nodes_indices_i <- which(RISR_nodes_metadata_df$node %in% descendant_nodes_i)
    
    # Extract PP of descendants in focal regions
    PP_var_bioregion_ID <- which(names(RISR_nodes_metadata_df) == str_replace(string = paste0("PP_", bioregion_i), pattern = " ", replacement = "_"))
    PP_focal_bioregion_i <- RISR_nodes_metadata_df[descendant_nodes_indices_i, PP_var_bioregion_ID, drop = T]
    PP_var_OW_NW_ID <- which(names(RISR_nodes_metadata_df) == str_replace(string = paste0("PP_", OW_NW_i), pattern = " ", replacement = "_"))
    PP_focal_OW_NW_i <- RISR_nodes_metadata_df[descendant_nodes_indices_i, PP_var_OW_NW_ID, drop = T]
    PP_var_Trop_Temp_ID <- which(names(RISR_nodes_metadata_df) == str_replace(string = paste0("PP_", Trop_Temp_i), pattern = " ", replacement = "_"))
    PP_focal_Trop_Temp_i <- RISR_nodes_metadata_df[descendant_nodes_indices_i, PP_var_Trop_Temp_ID, drop = T]
    
    # Extract length of descendant edges
    edge_length_i <- RISR_nodes_metadata_df$edge_length[descendant_nodes_indices_i]
    
    # Compute residence times in focal region
    # Residence times = length x PP in focal region
    RISR_nodes_metadata_df$residence_time_desc_sum_Bioregion[i] <- sum(PP_focal_bioregion_i * edge_length_i)
    RISR_nodes_metadata_df$residence_time_desc_sum_OW_NW[i] <- sum(PP_focal_OW_NW_i * edge_length_i)
    RISR_nodes_metadata_df$residence_time_desc_sum_Trop_Temp[i] <- sum(PP_focal_Trop_Temp_i * edge_length_i)
    
    # Convert in proportion
    RISR_nodes_metadata_df$residence_time_desc_prop_Bioregion[i] <- sum(PP_focal_bioregion_i * edge_length_i) / sum(edge_length_i)
    RISR_nodes_metadata_df$residence_time_desc_prop_OW_NW[i] <- sum(PP_focal_OW_NW_i * edge_length_i) / sum(edge_length_i)
    RISR_nodes_metadata_df$residence_time_desc_prop_Trop_Temp[i] <- sum(PP_focal_Trop_Temp_i * edge_length_i) / sum(edge_length_i)
  }
}
  
hist(RISR_nodes_metadata_df$residence_time_desc_prop_Bioregion)
hist(RISR_nodes_metadata_df$residence_time_desc_prop_OW_NW)
hist(RISR_nodes_metadata_df$residence_time_desc_prop_Trop_Temp)

# Save nodes metadata for RISR
saveRDS(RISR_nodes_metadata_df, "./outputs/RISR/RISR_nodes_metadata_df.rds")


### 2.2/ Identify RISR per Bioregions ~ endemicity threshold ####

## 2.2.1/ Identify RISR per Bioregions ~ endemicity threshold

# Load nodes metadata for RISR
RISR_nodes_metadata_df <- readRDS(file = "./outputs/RISR/RISR_nodes_metadata_df.rds")


bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

endemicity_threshold_list <- seq(from = 0.70, to = 1.00, by = 0.05)

# Initiate RISR status variables
RISR_nodes_metadata_df$RISR_Bioregion_0.70 <- FALSE
RISR_nodes_metadata_df$RISR_Bioregion_0.75 <- FALSE
RISR_nodes_metadata_df$RISR_Bioregion_0.80 <- FALSE
RISR_nodes_metadata_df$RISR_Bioregion_0.85 <- FALSE
RISR_nodes_metadata_df$RISR_Bioregion_0.90 <- FALSE 
RISR_nodes_metadata_df$RISR_Bioregion_0.95 <- FALSE
RISR_nodes_metadata_df$RISR_Bioregion_1.00 <- FALSE

## Loop per endemicity threshold
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_Bioregion_", format(threshold_i, nsmall = 2)))
  
  # Detect candidate edges with endemicity higher than the threshold
  candidate_nodes_df_i <- RISR_nodes_metadata_df %>% 
    filter(residence_time_desc_prop_Bioregion >= threshold_i) %>% 
    filter(nb_descendant_tips >= 5) # Remove clades with less than 5 tips
  
  # table(candidate_nodes_df_i$Bioregion)
  
  # ## Loop per Bioregions (may lead to non-independant RISR clades across Bioregions) 
  # for (j in seq_along(bioregion_names))
  # {
  #   # j <- 4
  #   
  #   # Extract focal bioregion
  #   bioregion_j <- bioregion_names[j]
  #   
  #   # Extract candidate nodes for the focal bioregion
  #   candidate_nodes_df_ij <- candidate_nodes_df_i %>% 
  #     filter(Bioregion == bioregion_j) %>%
  #     arrange(desc(node_age)) %>%
  #     dplyr::select(node, node_age, residence_time_desc_prop_Bioregion, nb_descendant_tips)
  #   
  #   RISR_nodes_df_ij <- c()
  #   candidate_nodes_df_ij_temp <- candidate_nodes_df_ij
  #     
  #   while (nrow(candidate_nodes_df_ij_temp) > 0)
  #   {
  #     # Extract oldest candidate node
  #     candidate_node_k <- candidate_nodes_df_ij_temp$node[1]
  #     # Identify descendants
  #     descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
  #                                                    node = candidate_node_k)
  #     descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
  #     # Remove descendant nodes from the identified clade
  #     candidate_nodes_df_ij_temp <- candidate_nodes_df_ij_temp %>% 
  #       filter(!(node %in% descendant_nodes_k))
  #     # Store candidate node
  #     RISR_nodes_df_ij <- c(RISR_nodes_df_ij, candidate_node_k)
  #   }
  #   
  #   # Store identified RISR nodes
  #   RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_ij, RISR_status_var_index] <- TRUE
  # }
  
  # Extract candidate nodes
  candidate_nodes_df_i <- candidate_nodes_df_i %>%
    arrange(desc(node_age)) %>%
    dplyr::select(node, node_age, residence_time_desc_prop_Bioregion, nb_descendant_tips)
  
  RISR_nodes_df_i <- c()
  candidate_nodes_df_i_temp <- candidate_nodes_df_i
  
  while (nrow(candidate_nodes_df_i_temp) > 0)
  {
    # Extract oldest candidate node
    candidate_node_k <- candidate_nodes_df_i_temp$node[1]
    # Identify descendants
    descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                   node = candidate_node_k)
    descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
    # Remove descendant nodes from the identified clade
    candidate_nodes_df_i_temp <- candidate_nodes_df_i_temp %>%
      filter(!(node %in% descendant_nodes_k))
    # Store candidate node
    RISR_nodes_df_i <- c(RISR_nodes_df_i, candidate_node_k)
  }
  
  # Store identified RISR nodes
  RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_i, RISR_status_var_index] <- TRUE
    
  # Print progress
  cat(paste0(Sys.time(), " - RISR nodes identified for Bioregions for endemicity threshold = ", format(threshold_i, nsmall = 2)," - n° ", i, "/", length(endemicity_threshold_list),"\n"))
}

## Check that the identified RISR are independent

RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_Bioregion_0.80]
RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_Bioregion_1.00]

RISR_descendant_nodes <- c()
for (i in seq_along(RISR_nodes))
{
  # i <- 1
  
  RISR_node_i <- RISR_nodes[i]
  
  # Identify descendants
  descendant_nodes_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                 node = RISR_node_i)
  descendant_nodes_i <- unique(c(descendant_nodes_i, RISR_node_i)) # Add focal node itself
  
  # Store all descendant nodes
  RISR_descendant_nodes <- c(RISR_descendant_nodes, descendant_nodes_i)
  
}

# Only one occurrence per descending nodes across the RISR clades
table(table(RISR_descendant_nodes))

## Save nodes metadata for RISR updated with RISR nodes identification for Bioregions
saveRDS(RISR_nodes_metadata_df, "./outputs/RISR/RISR_nodes_metadata_df.rds")

## 2.2.2/ Compute metadata of RISR clades per endemicity threshold 

RISR_clades_for_Bioregions_df <- data.frame()
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_Bioregion_", format(threshold_i, nsmall = 2)))
  RISR_status_var <- RISR_nodes_metadata_df[, RISR_status_var_index, drop = T]
  
  # Extract RISR nodes metadata
  RISR_clades_for_Bioregions_df_i <- RISR_nodes_metadata_df %>% 
    filter(RISR_status_var) %>%
    dplyr::select(node, node_age, Bioregion, residence_time_desc_prop_Bioregion, nb_descendant_tips) %>% 
    mutate(endemicity_threshold = threshold_i) %>% 
    rename(Region = Bioregion)
  
  # Store RISR nodes metadata
  RISR_clades_for_Bioregions_df <- rbind(RISR_clades_for_Bioregions_df, RISR_clades_for_Bioregions_df_i)
}

# No RISR in Temperate bioregions as they are always nested in Tropical RISR!
View(RISR_clades_for_Bioregions_df)
table(RISR_clades_for_Bioregions_df$Region, RISR_clades_for_Bioregions_df$endemicity_threshold)

# Save metadata df of RISR clades per endemicity threshold 
saveRDS(RISR_clades_for_Bioregions_df, file = "./outputs/RISR/RISR_clades_for_Bioregions_df.rds")


### 2.3/ Identify RISR per Old World vs. New World ####

## 2.3.1/ Identify RISR per Old World vs. New World ~ endemicity threshold

# Load nodes metadata for RISR
RISR_nodes_metadata_df <- readRDS(file = "./outputs/RISR/RISR_nodes_metadata_df.rds")


OW_NW_names <- c("Old World", "New World")

endemicity_threshold_list <- seq(from = 0.70, to = 1.00, by = 0.05)

# Initiate RISR status variables
RISR_nodes_metadata_df$RISR_OW_NW_0.70 <- FALSE
RISR_nodes_metadata_df$RISR_OW_NW_0.75 <- FALSE
RISR_nodes_metadata_df$RISR_OW_NW_0.80 <- FALSE
RISR_nodes_metadata_df$RISR_OW_NW_0.85 <- FALSE
RISR_nodes_metadata_df$RISR_OW_NW_0.90 <- FALSE 
RISR_nodes_metadata_df$RISR_OW_NW_0.95 <- FALSE
RISR_nodes_metadata_df$RISR_OW_NW_1.00 <- FALSE

## Loop per endemicity threshold
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_OW_NW_", format(threshold_i, nsmall = 2)))
  
  # Detect candidate edges with endemicity higher than the threshold
  candidate_nodes_df_i <- RISR_nodes_metadata_df %>% 
    filter(residence_time_desc_prop_OW_NW >= threshold_i) %>% 
    filter(nb_descendant_tips >= 5) # Remove clades with less than 5 tips
  
  # table(candidate_nodes_df_i$Bioregion)
  
  # ## Loop per Bioregions (may lead to non-independant RISR clades across Bioregions) 
  # for (j in seq_along(bioregion_names))
  # {
  #   # j <- 4
  #   
  #   # Extract focal bioregion
  #   bioregion_j <- bioregion_names[j]
  #   
  #   # Extract candidate nodes for the focal bioregion
  #   candidate_nodes_df_ij <- candidate_nodes_df_i %>% 
  #     filter(Bioregion == bioregion_j) %>%
  #     arrange(desc(node_age)) %>%
  #     dplyr::select(node, node_age, residence_time_desc_prop_OW_NW, nb_descendant_tips)
  #   
  #   RISR_nodes_df_ij <- c()
  #   candidate_nodes_df_ij_temp <- candidate_nodes_df_ij
  #     
  #   while (nrow(candidate_nodes_df_ij_temp) > 0)
  #   {
  #     # Extract oldest candidate node
  #     candidate_node_k <- candidate_nodes_df_ij_temp$node[1]
  #     # Identify descendants
  #     descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
  #                                                    node = candidate_node_k)
  #     descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
  #     # Remove descendant nodes from the identified clade
  #     candidate_nodes_df_ij_temp <- candidate_nodes_df_ij_temp %>% 
  #       filter(!(node %in% descendant_nodes_k))
  #     # Store candidate node
  #     RISR_nodes_df_ij <- c(RISR_nodes_df_ij, candidate_node_k)
  #   }
  #   
  #   # Store identified RISR nodes
  #   RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_ij, RISR_status_var_index] <- TRUE
  # }
  
  # Extract candidate nodes
  candidate_nodes_df_i <- candidate_nodes_df_i %>%
    arrange(desc(node_age)) %>%
    dplyr::select(node, node_age, residence_time_desc_prop_OW_NW, nb_descendant_tips)
  
  RISR_nodes_df_i <- c()
  candidate_nodes_df_i_temp <- candidate_nodes_df_i
  
  while (nrow(candidate_nodes_df_i_temp) > 0)
  {
    # Extract oldest candidate node
    candidate_node_k <- candidate_nodes_df_i_temp$node[1]
    # Identify descendants
    descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                   node = candidate_node_k)
    descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
    # Remove descendant nodes from the identified clade
    candidate_nodes_df_i_temp <- candidate_nodes_df_i_temp %>%
      filter(!(node %in% descendant_nodes_k))
    # Store candidate node
    RISR_nodes_df_i <- c(RISR_nodes_df_i, candidate_node_k)
  }
  
  # Store identified RISR nodes
  RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_i, RISR_status_var_index] <- TRUE
  
  # Print progress
  cat(paste0(Sys.time(), " - RISR nodes identified for OW vs. NW for endemicity threshold = ", format(threshold_i, nsmall = 2)," - n° ", i, "/", length(endemicity_threshold_list),"\n"))
}

## Check that the identified RISR are independent

RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_OW_NW_0.80]
RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_OW_NW_1.00]

RISR_descendant_nodes <- c()
for (i in seq_along(RISR_nodes))
{
  # i <- 1
  
  RISR_node_i <- RISR_nodes[i]
  
  # Identify descendants
  descendant_nodes_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                 node = RISR_node_i)
  descendant_nodes_i <- unique(c(descendant_nodes_i, RISR_node_i)) # Add focal node itself
  
  # Store all descendant nodes
  RISR_descendant_nodes <- c(RISR_descendant_nodes, descendant_nodes_i)
  
}

# Only one occurrence per descending nodes across the RISR clades
table(table(RISR_descendant_nodes))

## Save nodes metadata for RISR updated with RISR nodes identification for Bioregions
saveRDS(RISR_nodes_metadata_df, "./outputs/RISR/RISR_nodes_metadata_df.rds")

## 2.3.2/ Compute metadata of RISR clades per endemicity threshold 

RISR_clades_for_OW_NW_df <- data.frame()
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_OW_NW_", format(threshold_i, nsmall = 2)))
  RISR_status_var <- RISR_nodes_metadata_df[, RISR_status_var_index, drop = T]
  
  # Extract RISR nodes metadata
  RISR_clades_for_OW_NW_df_i <- RISR_nodes_metadata_df %>% 
    filter(RISR_status_var) %>%
    dplyr::select(node, node_age, OW_vs_NW, residence_time_desc_prop_OW_NW, nb_descendant_tips) %>% 
    rename(Region = OW_vs_NW) %>%
    mutate(endemicity_threshold = threshold_i)
  
  # Store RISR nodes metadata
  RISR_clades_for_OW_NW_df <- rbind(RISR_clades_for_OW_NW_df, RISR_clades_for_OW_NW_df_i)
}

# No RISR in Temperate bioregions as they are always nested in Tropical RISR!
View(RISR_clades_for_OW_NW_df)
table(RISR_clades_for_OW_NW_df$Region, RISR_clades_for_OW_NW_df$endemicity_threshold)

# Save metadata df of RISR clades per endemicity threshold 
saveRDS(RISR_clades_for_OW_NW_df, file = "./outputs/RISR/RISR_clades_for_OW_NW_df.rds")


### 2.4/ Identify RISR per Tropics vs. Temperate ####

## 2.4.1/ Identify RISR per Tropics vs. Temperate ~ endemicity threshold

# Load nodes metadata for RISR
RISR_nodes_metadata_df <- readRDS(file = "./outputs/RISR/RISR_nodes_metadata_df.rds")

Trop_Temp_names <- c("Tropics", "Temperate")

endemicity_threshold_list <- seq(from = 0.70, to = 1.00, by = 0.05)

# Initiate RISR status variables
RISR_nodes_metadata_df$RISR_Trop_Temp_0.70 <- FALSE
RISR_nodes_metadata_df$RISR_Trop_Temp_0.75 <- FALSE
RISR_nodes_metadata_df$RISR_Trop_Temp_0.80 <- FALSE
RISR_nodes_metadata_df$RISR_Trop_Temp_0.85 <- FALSE
RISR_nodes_metadata_df$RISR_Trop_Temp_0.90 <- FALSE 
RISR_nodes_metadata_df$RISR_Trop_Temp_0.95 <- FALSE
RISR_nodes_metadata_df$RISR_Trop_Temp_1.00 <- FALSE

## Loop per endemicity threshold
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_Trop_Temp_", format(threshold_i, nsmall = 2)))
  
  # Detect candidate edges with endemicity higher than the threshold
  candidate_nodes_df_i <- RISR_nodes_metadata_df %>% 
    filter(residence_time_desc_prop_Trop_Temp >= threshold_i) %>% 
    filter(nb_descendant_tips >= 5) # Remove clades with less than 5 tips
  
  # table(candidate_nodes_df_i$Bioregion)
  
  # ## Loop per Bioregions (may lead to non-independant RISR clades across Bioregions) 
  # for (j in seq_along(bioregion_names))
  # {
  #   # j <- 4
  #   
  #   # Extract focal bioregion
  #   bioregion_j <- bioregion_names[j]
  #   
  #   # Extract candidate nodes for the focal bioregion
  #   candidate_nodes_df_ij <- candidate_nodes_df_i %>% 
  #     filter(Bioregion == bioregion_j) %>%
  #     arrange(desc(node_age)) %>%
  #     dplyr::select(node, node_age, residence_time_desc_prop_Trop_Temp, nb_descendant_tips)
  #   
  #   RISR_nodes_df_ij <- c()
  #   candidate_nodes_df_ij_temp <- candidate_nodes_df_ij
  #     
  #   while (nrow(candidate_nodes_df_ij_temp) > 0)
  #   {
  #     # Extract oldest candidate node
  #     candidate_node_k <- candidate_nodes_df_ij_temp$node[1]
  #     # Identify descendants
  #     descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
  #                                                    node = candidate_node_k)
  #     descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
  #     # Remove descendant nodes from the identified clade
  #     candidate_nodes_df_ij_temp <- candidate_nodes_df_ij_temp %>% 
  #       filter(!(node %in% descendant_nodes_k))
  #     # Store candidate node
  #     RISR_nodes_df_ij <- c(RISR_nodes_df_ij, candidate_node_k)
  #   }
  #   
  #   # Store identified RISR nodes
  #   RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_ij, RISR_status_var_index] <- TRUE
  # }
  
  # Extract candidate nodes
  candidate_nodes_df_i <- candidate_nodes_df_i %>%
    arrange(desc(node_age)) %>%
    dplyr::select(node, node_age, residence_time_desc_prop_Trop_Temp, nb_descendant_tips)
  
  RISR_nodes_df_i <- c()
  candidate_nodes_df_i_temp <- candidate_nodes_df_i
  
  while (nrow(candidate_nodes_df_i_temp) > 0)
  {
    # Extract oldest candidate node
    candidate_node_k <- candidate_nodes_df_i_temp$node[1]
    # Identify descendants
    descendant_nodes_k <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                   node = candidate_node_k)
    descendant_nodes_k <- unique(c(descendant_nodes_k, candidate_node_k)) # Add focal node itself
    # Remove descendant nodes from the identified clade
    candidate_nodes_df_i_temp <- candidate_nodes_df_i_temp %>%
      filter(!(node %in% descendant_nodes_k))
    # Store candidate node
    RISR_nodes_df_i <- c(RISR_nodes_df_i, candidate_node_k)
  }
  
  # Store identified RISR nodes
  RISR_nodes_metadata_df[RISR_nodes_metadata_df$node %in% RISR_nodes_df_i, RISR_status_var_index] <- TRUE
  
  # Print progress
  cat(paste0(Sys.time(), " - RISR nodes identified for OW vs. NW for endemicity threshold = ", format(threshold_i, nsmall = 2)," - n° ", i, "/", length(endemicity_threshold_list),"\n"))
}

## Check that the identified RISR are independent

RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_Trop_Temp_0.80]
RISR_nodes <- RISR_nodes_metadata_df$node[RISR_nodes_metadata_df$RISR_Trop_Temp_1.00]

RISR_descendant_nodes <- c()
for (i in seq_along(RISR_nodes))
{
  # i <- 1
  
  RISR_node_i <- RISR_nodes[i]
  
  # Identify descendants
  descendant_nodes_i <- phytools::getDescendants(tree = Ponerinae_phylogeny_1534t_treedata@phylo,
                                                 node = RISR_node_i)
  descendant_nodes_i <- unique(c(descendant_nodes_i, RISR_node_i)) # Add focal node itself
  
  # Store all descendant nodes
  RISR_descendant_nodes <- c(RISR_descendant_nodes, descendant_nodes_i)
  
}

# Only one occurrence per descending nodes across the RISR clades
table(table(RISR_descendant_nodes))

## Save nodes metadata for RISR updated with RISR nodes identification for Bioregions
saveRDS(RISR_nodes_metadata_df, "./outputs/RISR/RISR_nodes_metadata_df.rds")

## 2.4.2/ Compute metadata of RISR clades per endemicity threshold 

RISR_clades_for_Trop_Temp_df <- data.frame()
for (i in seq_along(endemicity_threshold_list))
{
  # i <- 1
  
  # Extract endemicity threshold
  threshold_i <- endemicity_threshold_list[i]
  
  # Identify RISR status variable
  RISR_status_var_index <- which(names(RISR_nodes_metadata_df) == paste0("RISR_Trop_Temp_", format(threshold_i, nsmall = 2)))
  RISR_status_var <- RISR_nodes_metadata_df[, RISR_status_var_index, drop = T]
  
  # Extract RISR nodes metadata
  RISR_clades_for_Trop_Temp_df_i <- RISR_nodes_metadata_df %>% 
    filter(RISR_status_var) %>%
    dplyr::select(node, node_age, Trop_vs_Temp, residence_time_desc_prop_Trop_Temp, nb_descendant_tips) %>% 
    rename(Region = Trop_vs_Temp) %>%
    mutate(endemicity_threshold = threshold_i)
  
  # Store RISR nodes metadata
  RISR_clades_for_Trop_Temp_df <- rbind(RISR_clades_for_Trop_Temp_df, RISR_clades_for_Trop_Temp_df_i)
}

# No RISR in Temperate bioregions as they are always nested in Tropical RISR!
View(RISR_clades_for_Trop_Temp_df)
table(RISR_clades_for_Trop_Temp_df$Region, RISR_clades_for_Trop_Temp_df$endemicity_threshold)

# Save metadata df of RISR clades per endemicity threshold 
saveRDS(RISR_clades_for_Trop_Temp_df, file = "./outputs/RISR/RISR_clades_for_Trop_Temp_df.rds")

# Not a single Temperate RISR clade! They are all nested in Tropics RISR clades, or do not have enough tips...


### 2.5/ Aggregate RISR clades metadata across region schemes ####

# Load metadata df of RISR clades per endemicity threshold per Bioregions
RISR_clades_for_Bioregions_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_df.rds")

# Load metadata df of RISR clades per endemicity threshold for Old World vs. New World
RISR_clades_for_OW_NW_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_OW_NW_df.rds")


# Aggregate nb of RISR clades and RISR descendants per endemicity threshold per Bioregion scheme
RISR_clades_for_all_Bioregions_df_per_threshold <- RISR_clades_for_Bioregions_df %>% 
  group_by(endemicity_threshold) %>%
  summarize(RISR_clades_nb = n(),
            RISR_descendant_tips_nb = sum(nb_descendant_tips),
            RISR_descendant_tips_prop = sum(nb_descendant_tips) / length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)) %>%
  ungroup() %>%
  mutate(Region = "all_Bioregions",
         Region_scheme = "Bioregions")

RISR_clades_for_all_Bioregions_df_per_threshold

# Aggregate nb of RISR clades and RISR descendants per endemicity threshold for Old World vs. New World scheme
RISR_clades_for_all_OW_NW_df_per_threshold <- RISR_clades_for_OW_NW_df %>% 
  group_by(endemicity_threshold) %>%
  summarize(RISR_clades_nb = n(),
            RISR_descendant_tips_nb = sum(nb_descendant_tips),
            RISR_descendant_tips_prop = sum(nb_descendant_tips) / length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)) %>%
  ungroup() %>%
  mutate(Region = "all_Hemispheres",
         Region_scheme = "OW_vs_NW")

RISR_clades_for_all_OW_NW_df_per_threshold

# Aggregate nb of RISR clades and RISR descendants per endemicity threshold per Bioregions
RISR_clades_for_Bioregions_df_per_threshold <- RISR_clades_for_Bioregions_df %>% 
  group_by(endemicity_threshold, Region) %>%
  summarize(RISR_clades_nb = n(),
            RISR_descendant_tips_nb = sum(nb_descendant_tips),
            RISR_descendant_tips_prop = sum(nb_descendant_tips) / length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)) %>%
  ungroup() %>%
  mutate(Region_scheme = "Bioregions")

RISR_clades_for_Bioregions_df_per_threshold

# Aggregate nb of RISR clades and RISR descendants per endemicity threshold per Old World vs. New World
RISR_clades_for_OW_NW_df_per_threshold <- RISR_clades_for_OW_NW_df %>% 
  group_by(endemicity_threshold, Region) %>%
  summarize(RISR_clades_nb = n(),
            RISR_descendant_tips_nb = sum(nb_descendant_tips),
            RISR_descendant_tips_prop = sum(nb_descendant_tips) / length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)) %>%
  ungroup() %>%
  mutate(Region_scheme = "OW_vs_NW")

RISR_clades_for_OW_NW_df_per_threshold

# Merge all summary tables
RISR_clades_per_threshold_df <- rbind(RISR_clades_for_Bioregions_df_per_threshold, RISR_clades_for_all_Bioregions_df_per_threshold,
                                      RISR_clades_for_OW_NW_df_per_threshold, RISR_clades_for_all_OW_NW_df_per_threshold)

View(RISR_clades_per_threshold_df)

# Save RISR clades metadata per endemicity threshold
saveRDS(RISR_clades_per_threshold_df, file = "./outputs/RISR/RISR_clades_per_threshold_df.rds")


### 2.6/  Plot nb of RISR per endemicity threshold ####

# Set color scheme for regions
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_regions <- colors_list_for_areas[c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")]
colors_list_for_regions <- c(colors_list_for_regions, "black", "peachpuff2", "mediumpurple2", "grey")
names(colors_list_for_regions)[5:8] <- c("all Bioregions", "New World", "Old World", "all Hemispheres")

region_breaks <- str_replace(string = names(colors_list_for_regions), pattern = " ", replacement = "_")

# GGplot
pdf(file = "./outputs/RISR/RISR_clades_nb_per_threshold_per_Regions.pdf", height = 8, width = 12)

RISR_clades_nbs_per_threshold_per_Regions_plot <- ggplot(data = RISR_clades_per_threshold_df) +
  
  geom_line(mapping = aes(y = RISR_clades_nb, x = endemicity_threshold*100, col = Region, lty = Region_scheme),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add vertical lines for selected thresholds
  geom_vline(xintercept = 80, linewidth = 1.5, linetype = "solid", col = "red") +
  geom_vline(xintercept = 95, linewidth = 1.5, linetype = "dashed", col = "red") +

  # Adjust color scheme and legend
  scale_color_manual("Region", breaks = region_breaks, labels = names(colors_list_for_regions), values = unname(colors_list_for_regions)) +
  
  # Adjust linetype scheme and legend
  scale_linetype_manual("Region scheme", labels = c("Bioregions", "Hemispheres"), values = c("solid", 22)) +
  
  # Set plot title +
  ggtitle(label = paste0("Nb of RISR clades ~ Endemicity threshold")) +
  
  # Set axes labels
  xlab("Endemicity threshold [%]") +
  ylab("Nb of RISR clades") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_clades_nbs_per_threshold_per_Regions_plot)

dev.off()

### 2.7/ Plot proportion of tips included per endemicity threshold ####

# Set color scheme for regions
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_regions <- colors_list_for_areas[c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")]
colors_list_for_regions <- c(colors_list_for_regions, "black", "peachpuff2", "mediumpurple2", "grey")
names(colors_list_for_regions)[5:8] <- c("all Bioregions", "New World", "Old World", "all Hemispheres")

region_breaks <- str_replace(string = names(colors_list_for_regions), pattern = " ", replacement = "_")

# GGplot
pdf(file = "./outputs/RISR/RISR_tips_prop_per_threshold_per_Regions.pdf", height = 8, width = 12)

RISR_tips_prop_per_threshold_per_Regions_plot <- ggplot(data = RISR_clades_per_threshold_df) +
  
  geom_line(mapping = aes(y = RISR_descendant_tips_prop*100, x = endemicity_threshold*100, col = Region, lty = Region_scheme),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add vertical lines for selected thresholds
  geom_vline(xintercept = 80, linewidth = 1.5, linetype = "solid", col = "red") +
  geom_vline(xintercept = 95, linewidth = 1.5, linetype = "dashed", col = "red") +
  
  # Adjust color scheme and legend
  scale_color_manual("Region", breaks = region_breaks, labels = names(colors_list_for_regions), values = unname(colors_list_for_regions)) +
  
  # Adjust linetype scheme and legend
  scale_linetype_manual("Region scheme", labels = c("Bioregions", "Hemispheres"), values = c("solid", 22)) +
  
  # Set plot title +
  ggtitle(label = paste0("Tips included in RISR clades ~ Endemicity threshold")) +
  
  # Set axes labels
  xlab("Endemicity threshold  [%]") +
  ylab("Tips included in RISR clades  [%]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_tips_prop_per_threshold_per_Regions_plot)

dev.off()


## Conclusion: Use a 80% endemicity threshold!

RISR_clades_per_threshold_df[as.character(RISR_clades_per_threshold_df$endemicity_threshold) == "0.8", ]

## Total bioregions: 47
# 16 Afrotropics
# 8 Australasia
# 15 Indomalaya
# 8 Neotropics

## Total independent OW vs. NW: 16
# 11 Old World
# 5 New World

## Total aggregated OW vs. NW: 47
# 41 Old World
# 8 New World


##### 3/ Fit time-dependent BD models #####

?diversitree::make.bd.t # For time-dependent BD models

# Load metadata df of RISR clades per endemicity threshold per Bioregions
RISR_clades_for_Bioregions_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_df.rds")

# Load metadata df of RISR clades per endemicity threshold for Old World vs. New World
RISR_clades_for_OW_NW_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_OW_NW_df.rds")

# Define selected threshold
selected_endemicity_threshold_Bioregions <- 0.80 # For Bioregions (and aggregated OW/NW)
selected_endemicity_threshold_OW_NW <- 0.95 # For independent OW/NW

# Extract RISR for the selected threshold for Bioregions
RISR_clades_for_Bioregions_0.80_df <- RISR_clades_for_Bioregions_df %>% 
  filter(as.character(endemicity_threshold) == as.character(selected_endemicity_threshold_Bioregions))

# Extract RISR for the selected threshold for Old World vs. New World
RISR_clades_for_OW_NW_0.95_df <- RISR_clades_for_OW_NW_df %>% 
  filter(as.character(endemicity_threshold) == as.character(selected_endemicity_threshold_OW_NW))


### 3.1/ Fit time-dependent BD models for RISR clades in Bioregions ####

## With diversitree in ML

# Initiate variables for model parameters
RISR_clades_for_Bioregions_0.80_df$lambda_0 <- NA
RISR_clades_for_Bioregions_0.80_df$alpha <- NA
RISR_clades_for_Bioregions_0.80_df$mu <- NA
RISR_clades_for_Bioregions_0.80_df$net_div_0 <- NA
RISR_clades_for_Bioregions_0.80_df$lambda_current <- NA
RISR_clades_for_Bioregions_0.80_df$net_div_current <- NA
RISR_clades_for_Bioregions_0.80_df$BD_AICc <- NA
RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_AICc <- NA
RISR_clades_for_Bioregions_0.80_df$BD_Akaike_w <- NA
RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_Akaike_w <- NA

## Loop per RISR clades

for (i in 1:nrow(RISR_clades_for_Bioregions_0.80_df))
{
  # i <- 1
  
  # Extract the RISR root nodes
  node_i <- RISR_clades_for_Bioregions_0.80_df$node[i]
  bioregion_i <- RISR_clades_for_Bioregions_0.80_df$Region[i]
  
  # Extract RISR clade from the phylogeny
  phylo_i <- ape::extract.clade(phy = Ponerinae_phylogeny_1534t_treedata@phylo, node = node_i)
  # plot(phylo_i)
  
  # Step 1: Make constant birth-death likelihood function
  BD_Lk_fn <- make.bd(tree = phylo_i,
                      sampling.f = 1.00)
  
  # Step 2: Find MLE for the constant birth-death model using optimization function
  starting_values <- c(0.1,0.01)
  BD_MLE <- find.mle(func = BD_Lk_fn, x.init = starting_values)
  
  # BD_MLE
  
  # Save constant BD model output
  saveRDS(object = BD_MLE, file = paste0("./outputs/RISR/model_fits/BD_MLE_Bioregions_node_",node_i,".rds"))
  
  # Step 3: Make time-dependent BD likelihood function
  BD_lambda_var_Lk_fn <- make.bd.t(tree = phylo_i,
                                   sampling.f = 1.00,
                                   functions = c("exp.t","constant.t")) # set the time of relationship of rates with time
  
  # Step 4: Find MLE for the time-dependent BD model using optimization function
  # Start optimization with the MLE parameters from the constant BD
  starting_values <- c(BD_MLE$par[1], 0.00, BD_MLE$par[2]) # Baseline speciation (lambda.l), Decaying parameter for speciation (lambda.a), extinction (mu)
  BD_lambda_var_MLE <- find.mle(func = BD_lambda_var_Lk_fn, x.init = starting_values)
  
  # BD_lambda_var_MLE
  
  # Save model output
  saveRDS(object = BD_lambda_var_MLE, file = paste0("./outputs/RISR/model_fits/BD_lambda_var_MLE_Bioregions_node_",node_i,".rds"))
  
  # Extract MLE of parameters from our fitted model
  lambda_0 <- BD_lambda_var_MLE$par[1]
  alpha <- BD_lambda_var_MLE$par[2]
  mu <- BD_lambda_var_MLE$par[3]
  net_div_0 <- lambda_0 - mu
    
  # Compute current rates
  clade_age <- RISR_clades_for_Bioregions_0.80_df$node_age[i]
  lambda_current <- lambda_0 * exp(-alpha * clade_age)
  net_div_current <- lambda_current - mu
  
  # Store parameters
  RISR_clades_for_Bioregions_0.80_df$lambda_0[i] <- lambda_0
  RISR_clades_for_Bioregions_0.80_df$alpha[i] <- -alpha # Record as positive alpha = increasing rates ; negative alpha = decaying rates
  RISR_clades_for_Bioregions_0.80_df$mu[i] <- mu
  RISR_clades_for_Bioregions_0.80_df$net_div_0[i] <- net_div_0
  RISR_clades_for_Bioregions_0.80_df$lambda_current[i] <- lambda_current
  RISR_clades_for_Bioregions_0.80_df$net_div_current[i] <- net_div_current
  
  ## Plot estimated lambda parameter though time
  
  # Extract time range to create time scale of the RISR clade
  time_scale <- seq(from = clade_age, to = 0, length.out = 100)
  
  # Create exponential or linear functions of the parameters through time
  lambda_t <- lambda_0 * exp(-alpha * time_scale) # Computation use time since root
  mu_t <- rep(mu, length(time_scale))
  max_rates <- max(lambda_t, mu_t)
  
  pdf(file = paste0("./outputs/RISR/model_fits/BD_lambda_var_MLE_Bioregions_node_",node_i,".pdf"), width = 8, height = 6)
  # Set plotting parameters
  par(mar = c(5.1, 5.1, 2.6, 2.1))
  # Plot curves from the fitted model
  plot(x = clade_age - time_scale, # Transform in age rather than time since root
       y = lambda_t,
       main = paste0("RISR ",bioregion_i, " - Node n°", node_i),
       xlim = c(clade_age, 0), ylim = c(0, max_rates*1.2),
       type = "l" , col = "darkgreen", lwd = 3,
       bty = "n", xlab = "time", las = 1, cex.axis = 0.8,
       ylab = expression(paste("rate (",lambda," or ",mu,")")))
  lines(x = time_scale, y = mu_t, col = "red", lwd = 3)
  legend(x = "topright", inset = c(0.05, 0.05), lwd = 3, col = c("darkgreen","red"),
         legend = c(expression(paste("speciation (",lambda,")")),
                    expression(paste("extinction (",mu,")"))),
         cex = 0.8, bty = "n")
  dev.off()
  
  ## Compare model fits
  
  BD_MLE
  BD_lambda_var_MLE
  
  # Extract number of observations for all models
  nobs <- length(phylo_i$tip.label)
  
  # Combine models' fit n a list
  models_list <- list(BD_MLE, BD_lambda_var_MLE)
  
  # Extract ln-likelihood, number of parameters, and AIC from models' list
  models_comparison <- data.frame(model = c("Constant", "Lambda var"),
                                  logL = sapply(X = models_list, FUN = logLik),
                                  k = sapply(X = models_list, FUN = function(x) length(x$par)),
                                  AIC = sapply(X = models_list, FUN = AIC))
  # Compute AICc from AIC, k and nobs
  models_comparison$AICc <- unlist(purrr::map2(.x = models_comparison$AIC, nobs = nobs, .y = models_comparison$k, .f = compute_AICc))
  
  # Compute Akaike's weights from AIC
  models_comparison$Akaike_weights <- round(phytools::aic.w(models_comparison$AICc), 3) * 100
  
  # Display result
  models_comparison
  
  # Store model comparison results
  RISR_clades_for_Bioregions_0.80_df$BD_AICc[i] <- models_comparison$AICc[1]
  RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_AICc[i] <- models_comparison$AICc[2]
  RISR_clades_for_Bioregions_0.80_df$BD_Akaike_w[i] <- models_comparison$Akaike_weights[1]
  RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_Akaike_w[i] <- models_comparison$Akaike_weights[2]
  
  # Print progress
  cat(paste0(Sys.time(), " - Time-dependent BD models fit for RISR clades identified for Bioregions - n° ", i, "/", nrow(RISR_clades_for_Bioregions_0.80_df),"\n"))
  
}

# Save metadata df of RISR clades per Bioregions for the selected endemicity threshold
saveRDS(object = RISR_clades_for_Bioregions_0.80_df, file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")


### 3.2/ Fit time-dependent BD models for RISR clades for Old World vs. New World ####

## With diversitree in ML

# Initiate variables for model parameters
RISR_clades_for_OW_NW_0.95_df$lambda_0 <- NA
RISR_clades_for_OW_NW_0.95_df$alpha <- NA
RISR_clades_for_OW_NW_0.95_df$mu <- NA
RISR_clades_for_OW_NW_0.95_df$net_div_0 <- NA
RISR_clades_for_OW_NW_0.95_df$lambda_current <- NA
RISR_clades_for_OW_NW_0.95_df$net_div_current <- NA
RISR_clades_for_OW_NW_0.95_df$BD_AICc <- NA
RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_AICc <- NA
RISR_clades_for_OW_NW_0.95_df$BD_Akaike_w <- NA
RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_Akaike_w <- NA

## Loop per RISR clades

for (i in 1:nrow(RISR_clades_for_OW_NW_0.95_df))
{
  # i <- 1
  
  # Extract the RISR root nodes
  node_i <- RISR_clades_for_OW_NW_0.95_df$node[i]
  bioregion_i <- RISR_clades_for_OW_NW_0.95_df$Region[i]
  
  # Extract RISR clade from the phylogeny
  phylo_i <- ape::extract.clade(phy = Ponerinae_phylogeny_1534t_treedata@phylo, node = node_i)
  # plot(phylo_i)
  
  # Step 1: Make constant birth-death likelihood function
  BD_Lk_fn <- make.bd(tree = phylo_i,
                      sampling.f = 1.00)
  
  # Step 2: Find MLE for the constant birth-death model using optimization function
  starting_values <- c(0.1,0.01)
  BD_MLE <- find.mle(func = BD_Lk_fn, x.init = starting_values)
  
  # BD_MLE
  
  # Save constant BD model output
  saveRDS(object = BD_MLE, file = paste0("./outputs/RISR/model_fits/BD_MLE_OW_NW_node_",node_i,".rds"))
  
  # Step 3: Make time-dependent BD likelihood function
  BD_lambda_var_Lk_fn <- make.bd.t(tree = phylo_i,
                                   sampling.f = 1.00,
                                   functions = c("exp.t","constant.t")) # set the time of relationship of rates with time
  
  # Step 4: Find MLE for the time-dependent BD model using optimization function
  # Start optimization with the MLE parameters from the constant BD
  starting_values <- c(BD_MLE$par[1], 0.00, BD_MLE$par[2]) # Baseline speciation (lambda.l), Decaying parameter for speciation (lambda.a), extinction (mu)
  BD_lambda_var_MLE <- find.mle(func = BD_lambda_var_Lk_fn, x.init = starting_values)
  
  # BD_lambda_var_MLE
  
  # Save model output
  saveRDS(object = BD_lambda_var_MLE, file = paste0("./outputs/RISR/model_fits/BD_lambda_var_MLE_OW_NW_node_",node_i,".rds"))
  
  # Extract MLE of parameters from our fitted model
  lambda_0 <- BD_lambda_var_MLE$par[1]
  alpha <- BD_lambda_var_MLE$par[2]
  mu <- BD_lambda_var_MLE$par[3]
  net_div_0 <- lambda_0 - mu
  
  # Compute current rates
  clade_age <- RISR_clades_for_OW_NW_0.95_df$node_age[i]
  lambda_current <- lambda_0 * exp(-alpha * clade_age)
  net_div_current <- lambda_current - mu
  
  # Store parameters
  RISR_clades_for_OW_NW_0.95_df$lambda_0[i] <- lambda_0
  RISR_clades_for_OW_NW_0.95_df$alpha[i] <- -alpha # Record as positive alpha = increasing rates ; negative alpha = decaying rates
  RISR_clades_for_OW_NW_0.95_df$mu[i] <- mu
  RISR_clades_for_OW_NW_0.95_df$net_div_0[i] <- net_div_0
  RISR_clades_for_OW_NW_0.95_df$lambda_current[i] <- lambda_current
  RISR_clades_for_OW_NW_0.95_df$net_div_current[i] <- net_div_current
  
  ## Plot estimated lambda parameter though time
  
  # Extract time range to create time scale of the RISR clade
  time_scale <- seq(from = clade_age, to = 0, length.out = 100)
  
  # Create exponential or linear functions of the parameters through time
  lambda_t <- lambda_0 * exp(-alpha * time_scale) # Computation use time since root
  mu_t <- rep(mu, length(time_scale))
  max_rates <- max(lambda_t, mu_t)
  
  pdf(file = paste0("./outputs/RISR/model_fits/BD_lambda_var_MLE_OW_NW_node_",node_i,".pdf"), width = 8, height = 6)
  # Set plotting parameters
  par(mar = c(5.1, 5.1, 2.6, 2.1))
  # Plot curves from the fitted model
  plot(x = clade_age - time_scale, # Transform in age rather than time since root
       y = lambda_t,
       main = paste0("RISR ",bioregion_i, " - Node n°", node_i),
       xlim = c(clade_age, 0), ylim = c(0, max_rates*1.2),
       type = "l" , col = "darkgreen", lwd = 3,
       bty = "n", xlab = "time", las = 1, cex.axis = 0.8,
       ylab = expression(paste("rate (",lambda," or ",mu,")")))
  lines(x = time_scale, y = mu_t, col = "red", lwd = 3)
  legend(x = "topright", inset = c(0.05, 0.05), lwd = 3, col = c("darkgreen","red"),
         legend = c(expression(paste("speciation (",lambda,")")),
                    expression(paste("extinction (",mu,")"))),
         cex = 0.8, bty = "n")
  dev.off()
  
  ## Compare model fits
  
  BD_MLE
  BD_lambda_var_MLE
  
  # Extract number of observations for all models
  nobs <- length(phylo_i$tip.label)
  
  # Combine models' fit n a list
  models_list <- list(BD_MLE, BD_lambda_var_MLE)
  
  # Extract ln-likelihood, number of parameters, and AIC from models' list
  models_comparison <- data.frame(model = c("Constant", "Lambda var"),
                                  logL = sapply(X = models_list, FUN = logLik),
                                  k = sapply(X = models_list, FUN = function(x) length(x$par)),
                                  AIC = sapply(X = models_list, FUN = AIC))
  # Compute AICc from AIC, k and nobs
  models_comparison$AICc <- unlist(purrr::map2(.x = models_comparison$AIC, nobs = nobs, .y = models_comparison$k, .f = compute_AICc))
  
  # Compute Akaike's weights from AIC
  models_comparison$Akaike_weights <- round(phytools::aic.w(models_comparison$AICc), 3) * 100
  
  # Display result
  models_comparison
  
  # Store model comparison results
  RISR_clades_for_OW_NW_0.95_df$BD_AICc[i] <- models_comparison$AICc[1]
  RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_AICc[i] <- models_comparison$AICc[2]
  RISR_clades_for_OW_NW_0.95_df$BD_Akaike_w[i] <- models_comparison$Akaike_weights[1]
  RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_Akaike_w[i] <- models_comparison$Akaike_weights[2]
  
  # Print progress
  cat(paste0(Sys.time(), " - Time-dependent BD models fit for RISR clades identified for Bioregions - n° ", i, "/", nrow(RISR_clades_for_OW_NW_0.95_df),"\n"))
  
}

# Save metadata df of RISR clades for Old World vs. New World for the selected endemicity threshold
saveRDS(object = RISR_clades_for_OW_NW_0.95_df, file = "./outputs/RISR/RISR_clades_for_OW_NW_0.95_df.rds")


##### 4/ Compare fits of models ####

### 4.1/ Compare fits of models for RISR clades in Bioregions ####

## Count relative frequency of EB selection as the best fit

Constant_count <- sum(RISR_clades_for_Bioregions_0.80_df$BD_Akaike_w > RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_Akaike_w)
Decaying_count <- sum(RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$BD_Akaike_w <= RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_Akaike_w] <= 0)
Growth_count <- sum(RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$BD_Akaike_w <= RISR_clades_for_Bioregions_0.80_df$BD_lambda_var_Akaike_w] > 0)

RISR_clades_for_Bioregions_model_comparison_df <- data.frame(model_type = c("Constant", "Decaying", "Growth"),
                                                             counts = c(Constant_count, Decaying_count, Growth_count))
RISR_clades_for_Bioregions_model_comparison_df

## Khi² goodness of fit test

RISR_clades_for_Bioregions_model_comparison_chisq_test <- chisq.test(x = RISR_clades_for_Bioregions_model_comparison_df$counts)

RISR_clades_for_Bioregions_model_comparison_chisq_test
RISR_clades_for_Bioregions_model_comparison_chisq_test$statistic # Khi² = 66.4
RISR_clades_for_Bioregions_model_comparison_chisq_test$p.value # p < 0.001
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_model_comparison_chisq_test$parameter) ; khi_Q95  # Q95% = 5.99

## Bar plot

# Reorder model types
model_types <- c("Decaying", "Constant", "Growth")
RISR_clades_for_Bioregions_model_comparison_df$model_type <- factor(x = RISR_clades_for_Bioregions_model_comparison_df$model_type, levels = model_types, labels = model_types)

# GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_model_comparison_barplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_model_comparison_barplot <- ggplot(data = RISR_clades_for_Bioregions_model_comparison_df) +
  
  geom_col(mapping = aes(y = counts, x = model_type, fill = model_type),
           alpha = 1.0) +
  
  geom_text(mapping = aes(y = counts/2, x = model_type, label = counts),
            col = "black", size = 5, fontface = "bold",
            alpha = 1.0) +
  
  # # Add density of alpha on background ????
  # geom_density(data = RISR_clades_for_Bioregions_0.80_df,
  #              mapping = aes(y = alpha),
  #              fill = NA, col = "black") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Model type", breaks = model_types, labels = model_types, values = c("royalblue", "grey90", "tomato")) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 40, label = "Khi² test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 33,
           label = paste0("Khi² = ", round(RISR_clades_for_Bioregions_model_comparison_chisq_test$statistic, 1), "\n",
                          "Q95 = ",  round(khi_Q95, 2), "\n",
                          "p < 0.001"),
           hjust = 0, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Best model fits of RISR clades for Bioregions")) +
  
  # Set axes labels
  xlab("Model type") +
  ylab("Counts") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 15)),
        legend.position = c(0.85, 0.6),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_model_comparison_barplot)

dev.off()


### 4.2/ Compare fits of models for RISR clades for Old World vs. New World ####

## Count relative frequency of EB selection as the best fit

Constant_count <- sum(RISR_clades_for_OW_NW_0.95_df$BD_Akaike_w > RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_Akaike_w)
Decaying_count <- sum(RISR_clades_for_OW_NW_0.95_df$alpha[RISR_clades_for_OW_NW_0.95_df$BD_Akaike_w <= RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_Akaike_w] <= 0)
Growth_count <- sum(RISR_clades_for_OW_NW_0.95_df$alpha[RISR_clades_for_OW_NW_0.95_df$BD_Akaike_w <= RISR_clades_for_OW_NW_0.95_df$BD_lambda_var_Akaike_w] > 0)

RISR_clades_for_OW_NW_model_comparison_df <- data.frame(model_type = c("Constant", "Decaying", "Growth"),
                                                             counts = c(Constant_count, Decaying_count, Growth_count))
RISR_clades_for_OW_NW_model_comparison_df

## Khi² goodness of fit test

RISR_clades_for_OW_NW_model_comparison_chisq_test <- chisq.test(x = RISR_clades_for_OW_NW_model_comparison_df$counts)

RISR_clades_for_OW_NW_model_comparison_chisq_test
RISR_clades_for_OW_NW_model_comparison_chisq_test$statistic # Khi² = 34.16
RISR_clades_for_OW_NW_model_comparison_chisq_test$p.value # p < 0.001
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_OW_NW_model_comparison_chisq_test$parameter) ; khi_Q95  # Q95% = 5.99

## Bar plot

# Reorder model types
model_types <- c("Decaying", "Constant", "Growth")
RISR_clades_for_OW_NW_model_comparison_df$model_type <- factor(x = RISR_clades_for_OW_NW_model_comparison_df$model_type, levels = model_types, labels = model_types)

# GGplot
pdf(file = "./outputs/RISR/RISR_for_OW_NW_model_comparison_barplot.pdf", height = 6, width = 8)

RISR_for_OW_NW_model_comparison_barplot <- ggplot(data = RISR_clades_for_OW_NW_model_comparison_df) +
  
  geom_col(mapping = aes(y = counts, x = model_type, fill = model_type),
           alpha = 1.0) +
  
  geom_text(mapping = aes(y = counts/2, x = model_type, label = counts),
            col = "black", size = 5, fontface = "bold",
            alpha = 1.0) +
  
  # # Add density of alpha on background ????
  # geom_density(data = RISR_clades_for_OW_NW_0.95_df,
  #              mapping = aes(y = alpha),
  #              fill = NA, col = "black") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Model type", breaks = model_types, labels = model_types, values = c("royalblue", "grey90", "tomato")) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 21, label = "Khi² test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 18,
           label = paste0("Khi² = ", round(RISR_clades_for_OW_NW_model_comparison_chisq_test$statistic, 1), "\n",
                          "Q95 = ",  round(khi_Q95, 2), "\n",
                          "p < 0.001"),
           hjust = 0, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Best model fits of RISR clades\nOld World vs. New World")) +
  
  # Set axes labels
  xlab("Model type") +
  ylab("Counts") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 15)),
        legend.position = c(0.85, 0.7),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_OW_NW_model_comparison_barplot)

dev.off()


##### 5/ Compare model parameters of RISR per Bioregions #####

# Load metadata df of RISR clades per Bioregions for the selected endemicity threshold
RISR_clades_for_Bioregions_0.80_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")

# Set color scheme for Bioregions
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_bioregions <- colors_list_for_areas[c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")]

# Reorder Bioregions
RISR_clades_for_Bioregions_0.80_df$Region <- factor(x = RISR_clades_for_Bioregions_0.80_df$Region, levels = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions))

### 5.1/ Plot current richness ~ Bioregions ####

## Overall Kruskal-Wallis test

RISR_clades_for_Bioregions_kruskal_test <- kruskal.test(x = RISR_clades_for_Bioregions_0.80_df$nb_descendant_tips,
                                                        g = RISR_clades_for_Bioregions_0.80_df$Region)
RISR_clades_for_Bioregions_kruskal_test
RISR_clades_for_Bioregions_kruskal_test$statistic # Khi² = 1.20
RISR_clades_for_Bioregions_kruskal_test$p.value # p = 0.753
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_kruskal_test$parameter) ; khi_Q95  # Q95% = 7.81

## Pairwise Mann-Whitney tests

bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)
pairwise_bioregion_levels <- t(combn(x = bioregion_levels, m = 2))
pairwise_bioregion_indices <- t(combn(x = 1:length(bioregion_levels), m = 2))
pairwise_bioregion_list <- split(pairwise_bioregion_indices, row(pairwise_bioregion_indices))
names(pairwise_bioregion_list) <- NULL

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests <- list()
for (i in 1:nrow(pairwise_bioregion_levels))
{
  # i <- 1
  
  # Extract pair of Bioregions
  bioregion_1 <- pairwise_bioregion_levels[i, 1]
  bioregion_2 <- pairwise_bioregion_levels[i, 2]
  
  wilcox_test_i <- wilcox.test(x = RISR_clades_for_Bioregions_0.80_df$nb_descendant_tips[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_1],
                               y = RISR_clades_for_Bioregions_0.80_df$nb_descendant_tips[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_2],
                               alternative = "two.sided")
  
  RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests[[i]] <- wilcox_test_i
  names(RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests)[i] <- paste0("Wilcox_test_", bioregion_1, "_", bioregion_2)
}

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests

# Detect significant pairwise tests
pairwise_signif <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value < 0.1 } ))
pairwise_signif
# No significant pairwise differences
  
## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_current_richness_boxplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_current_richness_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                       mapping = aes(y = nb_descendant_tips, x = Region)) +
  
  geom_boxplot(mapping = aes(y = nb_descendant_tips, x = Region, fill = Region),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Add pairwise comparisons
  ggpubr::stat_compare_means(method = "wilcox.test", p.adjust.method = "none",
                             comparisons = pairwise_bioregion_list[pairwise_signif],
                             hide.ns = TRUE,
                             label = NULL,  label.x = NULL, label.y = NULL) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 155, label = "Kruskal-Wallis test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 145,
           label = paste0("Khi² = ", format(x = round(RISR_clades_for_Bioregions_kruskal_test$statistic, 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p = ", format(round(RISR_clades_for_Bioregions_kruskal_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness of RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Current richness") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_current_richness_boxplot)

dev.off()


### 5.2/ Plot crown ages ~ Bioregions: evidence for “Time-for-accumulation hypothesis” ####

## Overall Kruskal-Wallis test

RISR_clades_for_Bioregions_kruskal_test <- kruskal.test(x = RISR_clades_for_Bioregions_0.80_df$node_age,
                                                        g = RISR_clades_for_Bioregions_0.80_df$Region)
RISR_clades_for_Bioregions_kruskal_test
RISR_clades_for_Bioregions_kruskal_test$statistic # Khi² = 5.38
RISR_clades_for_Bioregions_kruskal_test$p.value # p = 0.146
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_kruskal_test$parameter) ; khi_Q95  # Q95% = 7.81

## Pairwise Mann-Whitney tests

bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)
pairwise_bioregion_levels <- t(combn(x = bioregion_levels, m = 2))

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests <- list()
for (i in 1:nrow(pairwise_bioregion_levels))
{
  # i <- 1
  
  # Extract pair of Bioregions
  bioregion_1 <- pairwise_bioregion_levels[i, 1]
  bioregion_2 <- pairwise_bioregion_levels[i, 2]
  
  wilcox_test_i <- wilcox.test(x = RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_1],
                               y = RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_2],
                               alternative = "two.sided")
  
  RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests[[i]] <- wilcox_test_i
  names(RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests)[i] <- paste0("Wilcox_test_", bioregion_1, "_", bioregion_2)
}

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests

# Detect significant pairwise tests
pairwise_p_values <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value } ))
pairwise_signif <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value < 0.1 } ))
pairwise_p_values_signif <- round(pairwise_p_values[pairwise_signif], 3)

# Two slightly significant pairwise differences


## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_crown_age_boxplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_crown_age_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                       aes(y = node_age, x = Region)) +
  
  geom_boxplot(mapping = aes(y = node_age, x = Region, fill = Region),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # # Add pairwise comparisons automatically
  # ggpubr::stat_compare_means(method = "wilcox.test", p.adjust.method = "none",
  #                            comparisons = pairwise_bioregion_list[pairwise_signif],
  #                            bracket.size = 0.5, # label = "p.signif",
  #                            hide.ns = TRUE) +
  
  # Add pairwise comparisons manually
  geom_bracket(xmin = c("Afrotropics", "Australasia"), xmax = c("Australasia", "Neotropics"),
               y.position = c(80, 120),
               label.size = 5, size = 0.3,
               vjust = -0.5,
               label = pairwise_p_values_signif,
               # label = c(".", "*"),
               tip.length = 0.02) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 125, label = "Kruskal-Wallis test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 115,
           label = paste0("Khi² = ", format(round(RISR_clades_for_Bioregions_kruskal_test$statistic, 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p = ", format(round(RISR_clades_for_Bioregions_kruskal_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown age of RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Crown age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_crown_age_boxplot)

dev.off()


### 5.3/ Plot net diversification rates at crown (netDiv crown) ~ Bioregions => evidence for "Niche opportunity" ####

## Overall Kruskal-Wallis test

RISR_clades_for_Bioregions_kruskal_test <- kruskal.test(x = RISR_clades_for_Bioregions_0.80_df$net_div_0,
                                                        g = RISR_clades_for_Bioregions_0.80_df$Region)
RISR_clades_for_Bioregions_kruskal_test
RISR_clades_for_Bioregions_kruskal_test$statistic # Khi² = 5.53
RISR_clades_for_Bioregions_kruskal_test$p.value # p = 0.137
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_kruskal_test$parameter) ; khi_Q95  # Q95% = 7.81

## Pairwise Mann-Whitney tests

bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)
pairwise_bioregion_levels <- t(combn(x = bioregion_levels, m = 2))

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests <- list()
for (i in 1:nrow(pairwise_bioregion_levels))
{
  # i <- 1
  
  # Extract pair of Bioregions
  bioregion_1 <- pairwise_bioregion_levels[i, 1]
  bioregion_2 <- pairwise_bioregion_levels[i, 2]
  
  wilcox_test_i <- wilcox.test(x = RISR_clades_for_Bioregions_0.80_df$net_div_0[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_1],
                               y = RISR_clades_for_Bioregions_0.80_df$net_div_0[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_2],
                               alternative = "two.sided")
  
  RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests[[i]] <- wilcox_test_i
  names(RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests)[i] <- paste0("Wilcox_test_", bioregion_1, "_", bioregion_2)
}

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests

# Detect significant pairwise tests
pairwise_p_values <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value } ))
pairwise_signif <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value < 0.1 } ))
pairwise_p_values_signif <- round(pairwise_p_values[pairwise_signif], 3)

# Two slighly significant pairwise differences


## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_net_div_crown_boxplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_net_div_crown_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                aes(y = net_div_0, x = Region)) +
  
  geom_boxplot(mapping = aes(y = net_div_0, x = Region, fill = Region),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # # Add pairwise comparisons automatically
  # ggpubr::stat_compare_means(method = "wilcox.test", p.adjust.method = "none",
  #                            comparisons = pairwise_bioregion_list[pairwise_signif],
  #                            bracket.size = 0.5, # label = "p.signif",
  #                            hide.ns = TRUE) +
  
  # Add pairwise comparisons manually
  geom_bracket(xmin = c("Afrotropics", "Australasia"), xmax = c("Australasia", "Indomalaya"),
               y.position = c(0.38, 0.42),
               label.size = 5, size = 0.3,
               vjust = -0.5,
               label = pairwise_p_values_signif,
               # label = c(".", "*"),
               tip.length = 0.02) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.65, label = "Kruskal-Wallis test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.60,
           label = paste0("Khi² = ", format(round(RISR_clades_for_Bioregions_kruskal_test$statistic, 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p = ", format(round(RISR_clades_for_Bioregions_kruskal_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown net div. rates of RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Crown net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_net_div_crown_boxplot)

dev.off()

### 5.4/ Plot time variation (alpha) ~ Bioregions ####

## Overall Kruskal-Wallis test

RISR_clades_for_Bioregions_kruskal_test <- kruskal.test(x = RISR_clades_for_Bioregions_0.80_df$alpha,
                                                        g = RISR_clades_for_Bioregions_0.80_df$Region)
RISR_clades_for_Bioregions_kruskal_test
RISR_clades_for_Bioregions_kruskal_test$statistic # Khi² = 2.82
RISR_clades_for_Bioregions_kruskal_test$p.value # p = 0.421
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_kruskal_test$parameter) ; khi_Q95  # Q95% = 7.81

## Pairwise Mann-Whitney tests

bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)
pairwise_bioregion_levels <- t(combn(x = bioregion_levels, m = 2))

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests <- list()
for (i in 1:nrow(pairwise_bioregion_levels))
{
  # i <- 1
  
  # Extract pair of Bioregions
  bioregion_1 <- pairwise_bioregion_levels[i, 1]
  bioregion_2 <- pairwise_bioregion_levels[i, 2]
  
  wilcox_test_i <- wilcox.test(x = RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_1],
                               y = RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_2],
                               alternative = "two.sided")
  
  RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests[[i]] <- wilcox_test_i
  names(RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests)[i] <- paste0("Wilcox_test_", bioregion_1, "_", bioregion_2)
}

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests

# Detect significant pairwise tests
pairwise_p_values <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value } ))
pairwise_signif <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value < 0.1 } ))
pairwise_p_values_signif <- round(pairwise_p_values[pairwise_signif], 3)

# No significant pairwise differences


## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_alpha_boxplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_alpha_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                    aes(y = alpha, x = Region)) +
  
  geom_boxplot(mapping = aes(y = alpha, x = Region, fill = Region),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # # Add pairwise comparisons automatically
  # ggpubr::stat_compare_means(method = "wilcox.test", p.adjust.method = "none",
  #                            comparisons = pairwise_bioregion_list[pairwise_signif],
  #                            bracket.size = 0.5, # label = "p.signif",
  #                            hide.ns = TRUE) +
  
  # # Add pairwise comparisons manually
  # geom_bracket(xmin = c("Afrotropics", "Australasia"), xmax = c("Australasia", "Indomalaya"),
  #              y.position = c(0.38, 0.42),
  #              label.size = 5, size = 0.3,
  #              vjust = -0.5,
  #              label = pairwise_p_values_signif,
  #              # label = c(".", "*"),
  #              tip.length = 0.02) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.12, label = "Kruskal-Wallis test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.105,
           label = paste0("Khi² = ", format(round(RISR_clades_for_Bioregions_kruskal_test$statistic, 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p = ", format(round(RISR_clades_for_Bioregions_kruskal_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Time variation trends of RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste("Time variation (",alpha,")"))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_alpha_boxplot)

dev.off()

### 5.5/ Plot current net div. rates ~ Bioregions => current source/sinks of biodiversity ####

## Overall Kruskal-Wallis test

RISR_clades_for_Bioregions_kruskal_test <- kruskal.test(x = RISR_clades_for_Bioregions_0.80_df$net_div_current,
                                                        g = RISR_clades_for_Bioregions_0.80_df$Region)
RISR_clades_for_Bioregions_kruskal_test
RISR_clades_for_Bioregions_kruskal_test$statistic # Khi² = 2.29
RISR_clades_for_Bioregions_kruskal_test$p.value # p = 0.515
khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_kruskal_test$parameter) ; khi_Q95  # Q95% = 7.81

## Pairwise Mann-Whitney tests

bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)
pairwise_bioregion_levels <- t(combn(x = bioregion_levels, m = 2))

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests <- list()
for (i in 1:nrow(pairwise_bioregion_levels))
{
  # i <- 1
  
  # Extract pair of Bioregions
  bioregion_1 <- pairwise_bioregion_levels[i, 1]
  bioregion_2 <- pairwise_bioregion_levels[i, 2]
  
  wilcox_test_i <- wilcox.test(x = RISR_clades_for_Bioregions_0.80_df$net_div_current[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_1],
                               y = RISR_clades_for_Bioregions_0.80_df$net_div_current[RISR_clades_for_Bioregions_0.80_df$Region == bioregion_2],
                               alternative = "two.sided")
  
  RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests[[i]] <- wilcox_test_i
  names(RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests)[i] <- paste0("Wilcox_test_", bioregion_1, "_", bioregion_2)
}

RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests

# Detect significant pairwise tests
pairwise_p_values <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value } ))
pairwise_signif <- unlist(lapply(X = RISR_clades_for_Bioregions_pairwise_Mann_Whitney_tests, FUN = function (x) { x$p.value < 0.1 } ))
pairwise_p_values_signif <- round(pairwise_p_values[pairwise_signif], 3)

# No significant pairwise differences


## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_net_div_current_boxplot.pdf", height = 6, width = 8)

RISR_for_Bioregions_net_div_current_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                            aes(y = net_div_current, x = Region)) +
  
  geom_boxplot(mapping = aes(y = net_div_current, x = Region, fill = Region),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # # Add pairwise comparisons automatically
  # ggpubr::stat_compare_means(method = "wilcox.test", p.adjust.method = "none",
  #                            comparisons = pairwise_bioregion_list[pairwise_signif],
  #                            bracket.size = 0.5, # label = "p.signif",
  #                            hide.ns = TRUE) +
  
  # # Add pairwise comparisons manually
  # geom_bracket(xmin = c("Afrotropics", "Australasia"), xmax = c("Australasia", "Indomalaya"),
  #              y.position = c(0.38, 0.42),
  #              label.size = 5, size = 0.3,
  #              vjust = -0.5,
  #              label = pairwise_p_values_signif,
  #              # label = c(".", "*"),
  #              tip.length = 0.02) +

  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.35, label = "Kruskal-Wallis test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.31,
           label = paste0("Khi² = ", format(round(RISR_clades_for_Bioregions_kruskal_test$statistic, 2), nsmall = 2), "\n",
                          "Q95 = ", format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p = ", format(round(RISR_clades_for_Bioregions_kruskal_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current net div. rates of RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Current net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_net_div_current_boxplot)

dev.off()


### 5.6/ Plot current richness ~ crown age x Bioregions ####

## Fit regression lines per Bioregions with a GLM Poisson

RISR_clades_for_Bioregions_richness_age_GLM <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                   formula = nb_descendant_tips ~ node_age * Region,
                                                   family = "poisson")

summary(RISR_clades_for_Bioregions_richness_age_GLM) ; RISR_clades_for_Bioregions_richness_age_GLM_summary_table
RISR_clades_for_Bioregions_richness_age_GLM_coefs <- coefficients(RISR_clades_for_Bioregions_richness_age_GLM)
RISR_clades_for_Bioregions_richness_age_GLM_anova_table <- anova(RISR_clades_for_Bioregions_richness_age_GLM, test = 'LRT') ; RISR_clades_for_Bioregions_richness_age_GLM_anova_table

Q95 <- qchisq(p = 0.95, df = RISR_clades_for_Bioregions_richness_age_GLM_anova_table$Df[3])

## Compute segments for regression lines
# For Afrotropics (base level)
Afrotropics_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Afrotropics"]) * 0.9
Afrotropics_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Afrotropics"]) * 1.1
Afrotropics_ymin <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + Afrotropics_xmin * RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"]
Afrotropics_ymax <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + Afrotropics_xmax * RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"]
# For Australasia
Australasia_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Australasia"]) * 0.9
Australasia_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Australasia"]) * 1.1
Australasia_ymin <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionAustralasia"] + Australasia_xmin * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionAustralasia"])
Australasia_ymax <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionAustralasia"] + Australasia_xmax * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionAustralasia"])
# For Indomalaya
Indomalaya_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Indomalaya"]) * 0.9
Indomalaya_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Indomalaya"]) * 1.1
Indomalaya_ymin <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionIndomalaya"] + Indomalaya_xmin * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionIndomalaya"])
Indomalaya_ymax <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionIndomalaya"] + Indomalaya_xmax * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionIndomalaya"])
# For Neotropics
Neotropics_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Neotropics"]) * 0.9
Neotropics_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Region == "Neotropics"]) * 1.1
Neotropics_ymin <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionNeotropics"] + Neotropics_xmin * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionNeotropics"])
Neotropics_ymax <- RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionNeotropics"] + Neotropics_xmax * (RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionNeotropics"])


## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_current_richness_per_crown_ages_plot.pdf", height = 6, width = 8)

RISR_for_Bioregions_current_richness_per_crown_ages_plot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                                   aes(y = nb_descendant_tips, x = node_age)) +
  
  geom_point(mapping = aes(y = nb_descendant_tips, x = node_age, col = Region),
             size = 3,
             alpha = 1.0, show.legend = T) +
  
  # Adjust col scheme and legend
  scale_color_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Convert Y-scale to log
  scale_y_continuous(transform = "log",
                     breaks = c(7.389, 20.085, 54.598, 148.41),
                     labels = c("7", "20", "55", "150")) +
  
  # Reverse X-scale
  # scale_x_continuous(transform = "reverse") +
  
  # # Add regression lines
  # # For Afrotropics (base level)
  # geom_abline(slope = RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"],
  #             intercept = RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"],
  #             col = colors_list_for_bioregions["Afrotropics"],
  #             linewidth  = 1.5) +
  # # For Australasia
  # geom_abline(slope = RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionAustralasia"],
  #             intercept = RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionAustralasia"],
  #             col = colors_list_for_bioregions["Australasia"],
  #             linewidth  = 1.5) +
  # # For Indomalaya
  # geom_abline(slope = RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionIndomalaya"],
  #             intercept = RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionIndomalaya"],
  #             col = colors_list_for_bioregions["Indomalaya"],
  #             linewidth  = 1.5) +
  # # For Neotropics
  # geom_abline(slope = RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["node_age:RegionNeotropics"],
  #             intercept = RISR_clades_for_Bioregions_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_Bioregions_richness_age_GLM_coefs["RegionNeotropics"],
  #             col = colors_list_for_bioregions["Neotropics"],
  #             linewidth  = 1.5) +
  
  # Add regression lines as segments
  # For Afrotropics (base level)
  geom_segment(x = Afrotropics_xmin, y = Afrotropics_ymin, xend = Afrotropics_xmax, yend = Afrotropics_ymax,
               col = colors_list_for_bioregions["Afrotropics"], linewidth  = 1.5) +
  # For Australasia
  geom_segment(x = Australasia_xmin, y = Australasia_ymin, xend = Australasia_xmax, yend = Australasia_ymax,
               col = colors_list_for_bioregions["Australasia"], linewidth  = 1.5) +
  # For Indomalaya
  geom_segment(x = Indomalaya_xmin, y = Indomalaya_ymin, xend = Indomalaya_xmax, yend = Indomalaya_ymax,
               col = colors_list_for_bioregions["Indomalaya"], linewidth  = 1.5) +
  # For Neotropics
  geom_segment(x = Neotropics_xmin, y = Neotropics_ymin, xend = Neotropics_xmax, yend = Neotropics_ymax,
               col = colors_list_for_bioregions["Neotropics"], linewidth  = 1.5) +
  
  # Add test results
  annotate(geom = "label", x = 125, y = 80, label = "LRT Bioregions",
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"), 
           hjust = 1, fontface = "bold", size = 4.5) +
  # Add test results
  annotate(geom = "label", x = 125, y = 65,
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"),
           label = paste0("Khi² = ", format(round(RISR_clades_for_Bioregions_richness_age_GLM_anova_table$Deviance[3], 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p < 0.001"),
           hjust = 1, vjust = 1, fontface = "plain", size = 4.5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness ~ Crown ages\nRISR clades per Bioregions")) +
  
  # Set axes labels
  xlab("Crown ages  [My]") +
  ylab("Current richness (log)") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size  = 14, margin = margin(b = 10), hjust = 0.8), 
        legend.text = element_text(size = 12),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.5, "line"),
        legend.spacing.y = unit(0.2, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_current_richness_per_crown_ages_plot)

dev.off()


### 5.7/ Multifaceted plot ####

## List all ggplot

RISR_for_Bioregions_all_plots_list <- list(RISR_for_Bioregions_current_richness_boxplot, # A/ Current richness
                                           RISR_for_Bioregions_crown_age_boxplot, # B/ Crown age
                                           RISR_for_Bioregions_current_richness_per_crown_ages_plot, # C/ Current richness ~ Crown age
                                           RISR_for_Bioregions_net_div_current_boxplot, # D/ Current net div.
                                           RISR_for_Bioregions_net_div_crown_boxplot, # E/ Crown net div.
                                           RISR_for_Bioregions_alpha_boxplot) # F/ Time variation
## Multifaceted plot
  # Rows = type of events
  # Columns = Time-strata

pdf(file = paste0("./outputs/RISR/RISR_for_Bioregions_all_plots.pdf"),
    height = 6*2, width = 8*3)

gridExtra::grid.arrange(
  grobs = RISR_for_Bioregions_all_plots_list, # List of ggplots
  widths = c(1, 1, 1),  # Width of columns
  heights = c(1, 1),
  nrow = 2,
  ncol = 3)

dev.off()


##### 6/ Compare model parameters of RISR for aggregated Old World vs. New World #####

# Load metadata df of RISR clades per Bioregions for the selected endemicity threshold
RISR_clades_for_Bioregions_0.80_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")

# Assign Old World vs. New World groups
RISR_clades_for_Bioregions_0.80_df$Hemisphere <- "Old World"
RISR_clades_for_Bioregions_0.80_df$Hemisphere[RISR_clades_for_Bioregions_0.80_df$Region %in% c("Nearctic", "Neotropics")] <- "New World"
table(RISR_clades_for_Bioregions_0.80_df$Hemisphere)

# Set color scheme for Old World vs. New World
colors_list_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_names <- c("Old World", "New World")
names(colors_list_for_OW_NW) <- OW_NW_names

# Reorder Hemispheres
RISR_clades_for_Bioregions_0.80_df$Hemisphere <- factor(x = RISR_clades_for_Bioregions_0.80_df$Hemisphere, levels = OW_NW_names, labels = OW_NW_names)

# Save metadata df of RISR clades per Bioregions for the selected endemicity threshold
saveRDS(RISR_clades_for_Bioregions_0.80_df, file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")

### 6.1/ Plot current richness ~ aggregated OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_Bioregions_0.80_df$nb_descendant_tips[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_Bioregions_0.80_df$nb_descendant_tips[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]

RISR_clades_for_aggregated_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                          y = NW_data,
                                                          alternative = "two.sided")
RISR_clades_for_aggregated_OW_NW_wilcox_test
RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic # W = 127.5
RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value # p = 0.427
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12))


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_current_richness_boxplot.pdf", height = 6, width = 6)

RISR_for_aggregated_OW_NW_current_richness_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                       mapping = aes(y = nb_descendant_tips, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = nb_descendant_tips, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 155, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 145,
           label = paste0("W = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness of RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Current richness") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_current_richness_boxplot)

dev.off()

### 6.2/ Plot crown ages ~ aggregated OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]

RISR_clades_for_aggregated_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                          y = NW_data,
                                                          alternative = "two.sided")
RISR_clades_for_aggregated_OW_NW_wilcox_test
RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic # W = 127.5
RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value # p = 0.427
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_crown_age_boxplot.pdf", height = 6, width = 6)

RISR_for_aggregated_OW_NW_crown_age_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                             mapping = aes(y = node_age, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = node_age, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 110, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 100,
           label = paste0("W = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown age of RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Crown age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_crown_age_boxplot)

dev.off()

### 6.3/ Plot net diversification rates at crown (netDiv crown) ~ aggregated OW vs. NW => evidence for "Niche opportunity" ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_Bioregions_0.80_df$net_div_0[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_Bioregions_0.80_df$net_div_0[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]

RISR_clades_for_aggregated_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                          y = NW_data,
                                                          alternative = "two.sided")
RISR_clades_for_aggregated_OW_NW_wilcox_test
RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic # W = 162.0
RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value # p = 0.87
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_net_div_crown_boxplot.pdf", height = 6, width = 6)

RISR_for_aggregated_OW_NW_net_div_crown_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                      mapping = aes(y = net_div_0, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = net_div_0, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.50, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.46,
           label = paste0("W = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown net div. rates of RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Crown net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_net_div_crown_boxplot)

dev.off()

### 6.4/ Plot time variation (alpha) ~ aggregated OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_Bioregions_0.80_df$alpha[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]

RISR_clades_for_aggregated_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                          y = NW_data,
                                                          alternative = "two.sided")
RISR_clades_for_aggregated_OW_NW_wilcox_test
RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic # W = 194.0
RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value # p = 0.294
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_alpha_boxplot.pdf", height = 6, width = 6)

RISR_for_aggregated_OW_NW_alpha_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                  mapping = aes(y = alpha, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = alpha, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.12, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.10,
           label = paste0("W = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Time variation trends of RISR clades\nper aggregated OW vs. NW")) +

  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste("Time variation (",alpha,")"))) +

  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_alpha_boxplot)

dev.off()


### 6.5/ Plot current net div. rates ~ aggregated OW vs. NW => current source/sinks of biodiversity? ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_Bioregions_0.80_df$net_div_current[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_Bioregions_0.80_df$net_div_current[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]

RISR_clades_for_aggregated_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                          y = NW_data,
                                                          alternative = "two.sided")
RISR_clades_for_aggregated_OW_NW_wilcox_test
RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic # W = 152.0
RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value # p = 0.922
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_net_div_current_boxplot.pdf", height = 6, width = 6)

RISR_for_aggregated_OW_NW_net_div_current_boxplot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                  mapping = aes(y = net_div_current, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = net_div_current, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.32, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.29,
           label = paste0("W = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_aggregated_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current net div. rates of RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Current net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_net_div_current_boxplot)

dev.off()


### 6.6/ Plot current richness ~ crown age x aggregated OW vs. NW ####

## Fit regression lines per aggregated OW vs. NW with a GLM Poisson

RISR_clades_for_aggregated_OW_NW_richness_age_GLM <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                   formula = nb_descendant_tips ~ node_age * Hemisphere,
                                                   family = "poisson")

summary(RISR_clades_for_aggregated_OW_NW_richness_age_GLM)
RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs <- coefficients(RISR_clades_for_aggregated_OW_NW_richness_age_GLM)
RISR_clades_for_aggregated_OW_NW_richness_age_GLM_anova_table <- anova(RISR_clades_for_aggregated_OW_NW_richness_age_GLM, test = 'LRT') ; RISR_clades_for_aggregated_OW_NW_richness_age_GLM_anova_table

khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_aggregated_OW_NW_richness_age_GLM_anova_table$Df[3])

## Compute segments for regression lines
# For Old World (base level)
Old_World_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]) * 0.9
Old_World_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "Old World"]) * 1.1
Old_World_ymin <- RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["(Intercept)"] + Old_World_xmin * RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age"]
Old_World_ymax <- RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["(Intercept)"] + Old_World_xmax * RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age"]
# For New World
New_World_xmin <- min(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]) * 0.9
New_World_xmax <- max(RISR_clades_for_Bioregions_0.80_df$node_age[RISR_clades_for_Bioregions_0.80_df$Hemisphere == "New World"]) * 1.1
New_World_ymin <- RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["HemisphereNew World"] + New_World_xmin * (RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age"] + RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age:HemisphereNew World"])
New_World_ymax <- RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["HemisphereNew World"] + New_World_xmax * (RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age"] + RISR_clades_for_aggregated_OW_NW_richness_age_GLM_coefs["node_age:HemisphereNew World"])


## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_current_richness_per_crown_ages_plot.pdf", height = 6, width = 7)

RISR_for_aggregated_OW_NW_current_richness_per_crown_ages_plot <- ggplot(data = RISR_clades_for_Bioregions_0.80_df,
                                                                   aes(y = nb_descendant_tips, x = node_age)) +
  
  geom_point(mapping = aes(y = nb_descendant_tips, x = node_age, col = Hemisphere),
             size = 3, alpha = 1.0, show.legend = T) +
  
  # Adjust col scheme and legend
  scale_color_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Convert Y-scale to log
  scale_y_continuous(transform = "log",
                     breaks = c(7.389, 20.085, 54.598, 148.41),
                     labels = c("7", "20", "55", "150")) +
  
  # Reverse X-scale
  # scale_x_continuous(transform = "reverse") +
  
  # Add regression lines as segments
  # For Old World (base level)
  geom_segment(x = Old_World_xmin, y = Old_World_ymin, xend = Old_World_xmax, yend = Old_World_ymax,
               col = colors_list_for_OW_NW["Old World"], linewidth  = 1.5) +
  # For New World
  geom_segment(x = New_World_xmin, y = New_World_ymin, xend = New_World_xmax, yend = New_World_ymax,
               col = colors_list_for_OW_NW["New World"], linewidth  = 1.5) +
  
  # Add test results
  annotate(geom = "label", x = 125, y = 60, label = "LRT Bioregions",
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"), 
           hjust = 1, fontface = "bold", size = 4.5) +
  # Add test results
  annotate(geom = "label", x = 125, y = 45,
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"),
           label = paste0("Khi² = ", format(round(RISR_clades_for_aggregated_OW_NW_richness_age_GLM_anova_table$Deviance[3], 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p < 0.001"),
           hjust = 1, vjust = 1, fontface = "plain", size = 4.5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness ~ Crown ages\nRISR clades per aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Crown ages  [My]") +
  ylab("Current richness (log)") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.position = c(0.85, 0.20),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size  = 14, margin = margin(b = 10), hjust = 0.8), 
        legend.text = element_text(size = 12),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.5, "line"),
        legend.spacing.y = unit(0.2, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_current_richness_per_crown_ages_plot)

dev.off()


### 6.7/ Multifaceted plot ####

## List all ggplot

RISR_for_aggregated_OW_NW_all_plots_list <- list(RISR_for_aggregated_OW_NW_current_richness_boxplot, # A/ Current richness
                                                 RISR_for_aggregated_OW_NW_crown_age_boxplot, # B/ Crown age
                                                 RISR_for_aggregated_OW_NW_current_richness_per_crown_ages_plot, # C/ Current richness ~ Crown age
                                                 RISR_for_aggregated_OW_NW_net_div_current_boxplot, # D/ Current net div.
                                                 RISR_for_aggregated_OW_NW_net_div_crown_boxplot, # E/ Crown net div.
                                                 RISR_for_aggregated_OW_NW_alpha_boxplot) # F/ Time variation
## Multifaceted plot
# Rows = type of events
# Columns = Time-strata

pdf(file = paste0("./outputs/RISR/RISR_for_aggregated_OW_NW_all_plots.pdf"),
    height = 6*2, width = 6.5*3)

gridExtra::grid.arrange(
  grobs = RISR_for_aggregated_OW_NW_all_plots_list, # List of ggplots
  widths = c(1, 1, 1),  # Width of columns
  heights = c(1, 1),
  nrow = 2,
  ncol = 3)

dev.off()


##### 7/ Compare model parameters of RISR for independent Old World vs. New World #####

# Load metadata df of RISR clades per independent for Old World vs. New World for the selected endemicity threshold
RISR_clades_for_OW_NW_0.95_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_OW_NW_0.95_df.rds")


# Assign Old World vs. New World groups
RISR_clades_for_OW_NW_0.95_df$Hemisphere <- str_replace(string = RISR_clades_for_OW_NW_0.95_df$Region, pattern = "_", replacement = " ")
table(RISR_clades_for_OW_NW_0.95_df$Hemisphere)

# Set color scheme for Old World vs. New World
colors_list_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_names <- c("Old World", "New World")
names(colors_list_for_OW_NW) <- OW_NW_names

# Reorder Hemispheres
RISR_clades_for_OW_NW_0.95_df$Hemisphere <- factor(x = RISR_clades_for_OW_NW_0.95_df$Hemisphere, levels = OW_NW_names, labels = OW_NW_names)


### 7.1/ Plot current richness ~ independent OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_OW_NW_0.95_df$nb_descendant_tips[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_OW_NW_0.95_df$nb_descendant_tips[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]

RISR_clades_for_independent_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                            y = NW_data,
                                                            alternative = "two.sided")
RISR_clades_for_independent_OW_NW_wilcox_test
RISR_clades_for_independent_OW_NW_wilcox_test$statistic # W = 92.5
RISR_clades_for_independent_OW_NW_wilcox_test$p.value # p = 0.161
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) #  96.2


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_current_richness_boxplot.pdf", height = 6, width = 6)

RISR_for_independent_OW_NW_current_richness_boxplot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                             mapping = aes(y = nb_descendant_tips, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = nb_descendant_tips, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Convert Y-scale to log
  scale_y_continuous(transform = "log",
                     breaks = c(7.389, 20.085, 54.598, 148.41, 403.43),
                     labels = c("7", "20", "55", "150", "400")
                     ) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 800, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 550,
           label = paste0("W = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness of RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Current richness") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_current_richness_boxplot)

dev.off()

### 7.2/ Plot crown ages ~ independent OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]

RISR_clades_for_independent_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                            y = NW_data,
                                                            alternative = "two.sided")
RISR_clades_for_independent_OW_NW_wilcox_test
RISR_clades_for_independent_OW_NW_wilcox_test$statistic # W = 101
RISR_clades_for_independent_OW_NW_wilcox_test$p.value # p = 0.057
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_crown_age_boxplot.pdf", height = 6, width = 6)

RISR_for_independent_OW_NW_crown_age_boxplot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                      mapping = aes(y = node_age, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = node_age, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 110, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 100,
           label = paste0("W = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown age of RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Crown age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_crown_age_boxplot)

dev.off()


### 7.3/ Plot net diversification rates at crown (netDiv crown) ~ independent OW vs. NW => evidence for "Niche opportunity" ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_OW_NW_0.95_df$net_div_0[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_OW_NW_0.95_df$net_div_0[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]

RISR_clades_for_independent_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                            y = NW_data,
                                                            alternative = "two.sided")
RISR_clades_for_independent_OW_NW_wilcox_test
RISR_clades_for_independent_OW_NW_wilcox_test$statistic # W = 44.0
RISR_clades_for_independent_OW_NW_wilcox_test$p.value # p = 0.174
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_net_div_crown_boxplot.pdf", height = 6, width = 6)

RISR_for_independent_OW_NW_net_div_crown_boxplot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                          mapping = aes(y = net_div_0, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = net_div_0, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.17, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.155,
           label = paste0("W = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Crown net div. rates of RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Hemispheres") +
  ylab("Crown net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_net_div_crown_boxplot)

dev.off()


### 7.4/ Plot time variation (alpha) ~ independent OW vs. NW ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_OW_NW_0.95_df$alpha[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_OW_NW_0.95_df$alpha[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]

RISR_clades_for_independent_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                            y = NW_data,
                                                            alternative = "two.sided")
RISR_clades_for_independent_OW_NW_wilcox_test
RISR_clades_for_independent_OW_NW_wilcox_test$statistic # W = 74.0
RISR_clades_for_independent_OW_NW_wilcox_test$p.value # p = 0.754
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_alpha_boxplot.pdf", height = 6, width = 6)

RISR_for_independent_OW_NW_alpha_boxplot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                  mapping = aes(y = alpha, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = alpha, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.10, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.09,
           label = paste0("W = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ",  format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Time variation trends of RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste("Time variation (",alpha,")"))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_alpha_boxplot)

dev.off()


### 7.5/ Plot current net div. rates ~ independent OW vs. NW => current source/sinks of biodiversity? ####

## Overall Wilcox-Mann-Whitney test

OW_data <- RISR_clades_for_OW_NW_0.95_df$net_div_current[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]
NW_data <- RISR_clades_for_OW_NW_0.95_df$net_div_current[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]

RISR_clades_for_independent_OW_NW_wilcox_test <- wilcox.test(x = OW_data,
                                                            y = NW_data,
                                                            alternative = "two.sided")
RISR_clades_for_independent_OW_NW_wilcox_test
RISR_clades_for_independent_OW_NW_wilcox_test$statistic # W = 71.0
RISR_clades_for_independent_OW_NW_wilcox_test$p.value # p = 0.887
W_Q95 <- qnorm(p = 0.95, mean = length(OW_data)*length(NW_data)/2, sd = sqrt(length(OW_data)*length(NW_data)*(length(OW_data)+length(NW_data)+1) / 12)) ; W_Q95


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_net_div_current_boxplot.pdf", height = 6, width = 6)

RISR_for_independent_OW_NW_net_div_current_boxplot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                            mapping = aes(y = net_div_current, x = Hemisphere)) +
  
  geom_boxplot(mapping = aes(y = net_div_current, x = Hemisphere, fill = Hemisphere),
               col = "black", alpha = 1.0, show.legend = F) +
  
  # Horizontal line for Rate = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.40, label = "Mann-Whitney test",
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "text", x = 0.55, y = 0.37,
           label = paste0("W = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$statistic, 1), nsmall = 1), "\n",
                          "Q95 = ", format(round(W_Q95, 1), nsmall = 1), "\n",
                          "p = ", format(round(RISR_clades_for_independent_OW_NW_wilcox_test$p.value, 3), nsmall = 3)),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current net div. rates of RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Current net div. rates") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_net_div_current_boxplot)

dev.off()


### 7.6/ Plot current richness ~ crown age x independent OW vs. NW ####

## Fit regression lines per independent OW vs. NW with a GLM Poisson

RISR_clades_for_independent_OW_NW_richness_age_GLM <- glm(data = RISR_clades_for_OW_NW_0.95_df,
                                                         formula = nb_descendant_tips ~ node_age * Hemisphere,
                                                         family = "poisson")

summary(RISR_clades_for_independent_OW_NW_richness_age_GLM)
RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs <- coefficients(RISR_clades_for_independent_OW_NW_richness_age_GLM)
RISR_clades_for_independent_OW_NW_richness_age_GLM_anova_table <- anova(RISR_clades_for_independent_OW_NW_richness_age_GLM, test = 'LRT') ; RISR_clades_for_independent_OW_NW_richness_age_GLM_anova_table

khi_Q95 <- qchisq(p = 0.95, df = RISR_clades_for_independent_OW_NW_richness_age_GLM_anova_table$Df[3])

## Compute segments for regression lines
# For Old World (base level)
Old_World_xmin <- min(RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]) * 0.9
Old_World_xmax <- max(RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "Old World"]) * 1.06
Old_World_ymin <- RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["(Intercept)"] + Old_World_xmin * RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age"]
Old_World_ymax <- RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["(Intercept)"] + Old_World_xmax * RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age"]
# For New World
New_World_xmin <- min(RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]) * 0.9
New_World_xmax <- max(RISR_clades_for_OW_NW_0.95_df$node_age[RISR_clades_for_OW_NW_0.95_df$Hemisphere == "New World"]) * 1.03
New_World_ymin <- RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["HemisphereNew World"] + New_World_xmin * (RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age"] + RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age:HemisphereNew World"])
New_World_ymax <- RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["(Intercept)"] + RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["HemisphereNew World"] + New_World_xmax * (RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age"] + RISR_clades_for_independent_OW_NW_richness_age_GLM_coefs["node_age:HemisphereNew World"])


## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_current_richness_per_crown_ages_plot.pdf", height = 6, width = 7)

RISR_for_independent_OW_NW_current_richness_per_crown_ages_plot <- ggplot(data = RISR_clades_for_OW_NW_0.95_df,
                                                                         aes(y = nb_descendant_tips, x = node_age)) +
  
  geom_point(mapping = aes(y = nb_descendant_tips, x = node_age, col = Hemisphere),
             size = 3, alpha = 1.0, show.legend = T) +
  
  # Adjust col scheme and legend
  scale_color_manual("Hemisphere", breaks = names(colors_list_for_OW_NW), labels = names(colors_list_for_OW_NW), values = colors_list_for_OW_NW) +
  
  # Convert Y-scale to log
  scale_y_continuous(transform = "log",
                     breaks = c(7.389, 20.085, 54.598, 148.41, 403.43),
                     labels = c("7", "20", "55", "150", "400")) +
  
  # Reverse X-scale
  # scale_x_continuous(transform = "reverse") +
  
  # Add regression lines as segments
  # For Old World (base level)
  geom_segment(x = Old_World_xmin, y = Old_World_ymin, xend = Old_World_xmax, yend = Old_World_ymax,
               col = colors_list_for_OW_NW["Old World"], linewidth  = 1.5) +
  # For New World
  geom_segment(x = New_World_xmin, y = New_World_ymin, xend = New_World_xmax, yend = New_World_ymax,
               col = colors_list_for_OW_NW["New World"], linewidth  = 1.5) +
  
  # Add test results
  annotate(geom = "label", x = 10, y = 500, label = "LRT Bioregions",
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"), 
           hjust = 0, fontface = "bold", size = 5) +
  # Add test results
  annotate(geom = "label", x = 10, y = 380,
           fill = "white", col = "black", label.size = 0, label.padding = unit(0.0, "lines"),
           label = paste0("Khi² = ", format(round(RISR_clades_for_independent_OW_NW_richness_age_GLM_anova_table$Deviance[3], 2), nsmall = 2), "\n",
                          "Q95 = ",  format(round(khi_Q95, 2), nsmall = 2), "\n",
                          "p < 0.001"),
           hjust = 0, vjust = 1, fontface = "plain", size = 5) +
  
  # Set plot title +
  ggtitle(label = paste0("Current richness ~ Crown ages\nRISR clades per independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Crown ages  [My]") +
  ylab("Current richness (log)") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.position = c(0.85, 0.20),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size  = 16, margin = margin(b = 10), hjust = 0.8), 
        legend.text = element_text(size = 14),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.5, "line"),
        legend.spacing.y = unit(0.2, "line"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_current_richness_per_crown_ages_plot)

dev.off()


### 7.7/ Multifaceted plot ####

## List all ggplot

RISR_for_independent_OW_NW_all_plots_list <- list(RISR_for_independent_OW_NW_current_richness_boxplot, # A/ Current richness
                                                 RISR_for_independent_OW_NW_crown_age_boxplot, # B/ Crown age
                                                 RISR_for_independent_OW_NW_current_richness_per_crown_ages_plot, # C/ Current richness ~ Crown age
                                                 RISR_for_independent_OW_NW_net_div_current_boxplot, # D/ Current net div.
                                                 RISR_for_independent_OW_NW_net_div_crown_boxplot, # E/ Crown net div.
                                                 RISR_for_independent_OW_NW_alpha_boxplot) # F/ Time variation
## Multifaceted plot
# Rows = type of events
# Columns = Time-strata

pdf(file = paste0("./outputs/RISR/RISR_for_independent_OW_NW_all_plots.pdf"),
    height = 6*2, width = 6.5*3)

gridExtra::grid.arrange(
  grobs = RISR_for_independent_OW_NW_all_plots_list, # List of ggplots
  widths = c(1, 1, 1),  # Width of columns
  heights = c(1, 1),
  nrow = 2,
  ncol = 3)

dev.off()


##### 8/ Identify factors explaining current richness, for Bioregions #####

# Load metadata df of RISR clades per Bioregions for the selected endemicity threshold
RISR_clades_for_Bioregions_0.80_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")

# Current richness ~ Binomial GLM with model selection based on AICc
  # If crown age: "Time-for-accumulation hypothesis"
  # If lambda_0: Cradle hypothesis (higher initial net div?)
  # If alpha: Museum hypothesis? (less decline in rates with time?)

### 8.1/ Run unstandardized model ####

RISR_clades_for_Bioregions_richness_model_GLM <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                     formula = nb_descendant_tips ~ (node_age + net_div_0 + alpha + 0) * (Region),
                                                     family = "poisson")
# Add + 0 to remove the Intercept and compute coefficients for each levels => RegionAfrotropics = (Intercept) ; RegionAustralasia = RegionAustralasia + (Intercept)
# Yet, still have to add the base-level slope coefficient to the other slope-levels coefficients

summary(RISR_clades_for_Bioregions_richness_model_GLM)
RISR_clades_for_Bioregions_richness_model_GLM_coefs <- coefficients(RISR_clades_for_Bioregions_richness_model_GLM)

# Type I Anova (sequential)
RISR_clades_for_Bioregions_richness_model_GLM_anova_table <- anova(RISR_clades_for_Bioregions_richness_model_GLM, test = 'LRT') ; RISR_clades_for_Bioregions_richness_model_GLM_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_Bioregions_richness_model_GLM, type = "II", test.statistic = "LR")

### 8.2/ Run standardized model ####

## Standardized version (only continuous predictors are scaled)

# reghelper::beta(model = RISR_clades_for_Bioregions_richness_model_GLM, skip = "Region")

RISR_clades_for_Bioregions_richness_model_GLM_scaled <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                            formula = nb_descendant_tips ~ (scale(node_age) + scale(net_div_0) + scale(alpha) + 0) * (Region),
                                                            family = "poisson")

summary(RISR_clades_for_Bioregions_richness_model_GLM_scaled)
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std <- coefficients(RISR_clades_for_Bioregions_richness_model_GLM_scaled)

# Type I Anova (sequential)
RISR_clades_for_Bioregions_richness_model_GLM_scaled_anova_table <- anova(RISR_clades_for_Bioregions_richness_model_GLM_scaled, test = 'LRT') ; RISR_clades_for_Bioregions_richness_model_GLM_scaled_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_Bioregions_richness_model_GLM_scaled, type = "II", test.statistic = "LR")

# Compare design matrices
model.matrix(RISR_clades_for_Bioregions_richness_model_GLM)
model.matrix(RISR_clades_for_Bioregions_richness_model_GLM_scaled)

### 8.3/ Extract Beta coefficients ####

# Extract standardized coefficients data
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_summary_table <- summary(RISR_clades_for_Bioregions_richness_model_GLM_scaled)$coefficients
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_summary_table

RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df <- data.frame(par_type = c("Crown age", "Crown net div. rate", "Time variation", rep("Bioregion", 4), rep("Crown age", 3), rep("Crown net div. rate", 3), rep("Time variation", 3)),
                                                                         bioregion = c(rep("Afrotropics", 4), rep(c("Australasia", "Indomalaya", "Neotropics"), 4)),
                                                                         mean = RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_summary_table[,1], 
                                                                         se = RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_summary_table[,2])

# Sort in logical order for plotting
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df %>% 
  arrange(par_type, bioregion) %>% 
  # Add par_name
  mutate(par_name = paste0(par_type, "_", bioregion)) %>% 
  dplyr::select(par_name, par_type, bioregion, mean, se)


# Add base level values to other Bioregion coefficients (except the Bioregion = Intercept)
# For Crown age
base_coef <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion == "Afrotropics"]
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] + base_coef
# For Crown net div. rate
base_coef <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion == "Afrotropics"]
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] + base_coef
# For Time variation
base_coef <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion == "Afrotropics"]
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion != "Afrotropics"] + base_coef

# Use Neotropics as the base level for Bioregion effects (only New World bioregion)
# Substract Neotropics Beta-coefficient to all Bioregion coefficents
base_coef <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Bioregion" & RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$bioregion == "Neotropics"]
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Bioregion"] <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df$par_type == "Bioregion"] - base_coef


# Get multiplier for Z-scores CI
Q2.5_mult <- qnorm(p = 0.025, mean = 0, sd = 1)
Q97.5_mult <- qnorm(p = 0.975, mean = 0, sd = 1)

# Add CI95% intervals
RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df <- RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df %>%
  mutate(CI2.5 = mean + Q2.5_mult*se,
         CI97.5 = mean + Q97.5_mult*se) %>%
  dplyr::select(par_name, par_type, bioregion, mean, se, CI2.5, CI97.5)

# Save std coeffs df
saveRDS(RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df, file = "./outputs/RISR/RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df.rds")


### 8.4/ Plot Beta coefficients ####

# Set color scheme for Bioregions
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_bioregions <- colors_list_for_areas[c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")]
bioregion_labels <- c("Afrotrp", "Austr", "IndoM", "Neotrp")

## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_richness_model_coefs_plot.pdf", height = 8, width = 8)

RISR_for_Bioregions_richness_model_coefs_plot <- ggplot(data = RISR_clades_for_Bioregions_richness_model_GLM_coefs_std_df) +
  
  # Add horizontal line: Beta-coef = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Add error bars
  geom_errorbar(mapping = aes(x = bioregion, y = mean, ymin = CI2.5, ymax = CI97.5, col = bioregion),
                width = 0.2, linewidth = 1.2,
                alpha = 1.0, show.legend = F) +
  
  # Add mean points
  geom_point(mapping = aes(x = bioregion, y = mean, fill = bioregion),
             shape = 21, size = 5, col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = bioregion_labels, values = colors_list_for_bioregions) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = bioregion_labels, values = colors_list_for_bioregions) +
  
  # Adjust x-axis labels
  scale_x_discrete(labels = c("Afrotropics" = "Afrotrp", "Australasia" = "Austr", "Indomalaya" = "IndoM", "Neotropics" = "Neotrp")) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ par_type,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
             ) +
  
  # Add an asterisk beside Neotropics bioregion coefficient to signify it is used as the reference level
  geom_text(data = data.frame(bioregion = "Neotropics", mean = 0.3, lab = "*",
                              par_type = factor("Bioregion", levels = c("Bioregion", "Crown age", "Crown net div. rate", "Time variation"))),
            aes(x = bioregion, y = mean, col = bioregion),
            label = "*", size = 10, nudge_x = 0.3, show.legend = F) +
  
  # Flip axes
  # coord_flip() + 
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste(beta,"-coefficients"))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_richness_model_coefs_plot)

dev.off()


### 8.5/ Compare % of variance (deviance) explained by each factor in each bioregion ####

# One model per Bioregion
bioregion_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Region)

## 8.5.1/ Run GLM and extract Anova Type II outputs ####

RISR_clades_for_Bioregions_0.80_Anova_df <- data.frame()
for (i in seq_along(bioregion_levels))
{
  # i <- 1
  
  # Extract bioregion
  bioregion_i <- bioregion_levels[i]
  
  # Extract RISR data
  RISR_clades_for_Bioregions_0.80_df_i <- RISR_clades_for_Bioregions_0.80_df %>% 
    filter(Region == bioregion_i)
  
  # Run GLM
  RISR_clades_for_Bioregions_richness_model_GLM_scaled_i <- glm(data = RISR_clades_for_Bioregions_0.80_df_i,
                                                              formula = nb_descendant_tips ~ scale(node_age) + scale(net_div_0) + scale(alpha),
                                                              family = "poisson")
  
  summary(RISR_clades_for_Bioregions_richness_model_GLM_scaled_i)
  
  # Type I Anova (sequential)
  # anova(RISR_clades_for_Bioregions_richness_model_GLM_scaled_i, test = 'LRT')
  
  # Use Type II Anova to obtain marginal deviance
  # Type II Anova (marginal)
  RISR_clades_for_Bioregions_richness_model_GLM_scaled_anova_table_i <- car::Anova(RISR_clades_for_Bioregions_richness_model_GLM_scaled_i, type = "II", test.statistic = "LR")
  RISR_clades_for_Bioregions_richness_model_GLM_scaled_anova_table_i
  
  # Extract predictor names
  predictors <- str_remove(string = names(RISR_clades_for_Bioregions_richness_model_GLM_scaled_i$model)[-1], pattern = "scale\\(")
  predictors <- str_remove(string = predictors, pattern = "\\)")
  
  # Extract % of marginal deviance explain by each factor
  # Convert to % using the Null Deviance
  RISR_clades_for_Bioregions_richness_model_Anofa_df_i <- data.frame(bioregion = bioregion_i,
                                                                     predictor = predictors,
                                                                     deviance_perc = round(RISR_clades_for_Bioregions_richness_model_GLM_scaled_anova_table_i$`LR Chisq` / RISR_clades_for_Bioregions_richness_model_GLM_scaled_i$null.deviance * 100, 1))
  
  # Store in final df
  RISR_clades_for_Bioregions_0.80_Anova_df <- rbind(RISR_clades_for_Bioregions_0.80_Anova_df, RISR_clades_for_Bioregions_richness_model_Anofa_df_i)
}

# Adjust predictor levels
RISR_clades_for_Bioregions_0.80_Anova_df$predictor <- factor(x = RISR_clades_for_Bioregions_0.80_Anova_df$predictor,
                                                             levels = c("node_age", "net_div_0", "alpha"),
                                                             labels = c("Crown age", "Crown net div. rate", "Time variation"))

# Save Anova type II summary df for Bioregions
saveRDS(object = RISR_clades_for_Bioregions_0.80_Anova_df, file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_Anova_df.rds")


## 8.5.2/ Plot % of marginal deviance ~ predictors for Bioregions ####

## GGplot
pdf(file = "./outputs/RISR/RISR_for_Bioregions_richness_model_Anova_predictors_plot.pdf", height = 9, width = 8)

RISR_for_Bioregions_richness_model_Anova_predictors_plot <- ggplot(data = RISR_clades_for_Bioregions_0.80_Anova_df) +
  
  # Add column bars
  geom_col(mapping = aes(x = predictor, y = deviance_perc, fill = bioregion),
           col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_bioregions), labels = names(colors_list_for_bioregions), values = colors_list_for_bioregions) +
  
  # Adjust x-axis labels, , 
  scale_x_discrete(labels = c("Crown age" = bquote(Age[crown]), "Crown net div. rate" = bquote(Rate[crown]), "Time variation" = bquote(TV~(alpha)))) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ bioregion,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
  ) +
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper Bioregions")) +
  
  # Set axes labels
  xlab("Predictors") +
  ylab("Marginal deviance  [%]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed", linewidth = 0.5),
        panel.grid.major.x = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5), angle = 0),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_Bioregions_richness_model_Anova_predictors_plot)

dev.off()



##### 9/ Identify factors explaining current richness, aggregated Old World vs. New World #####

# Load metadata df of RISR clades per aggregated OW vs. NW for the selected endemicity threshold
RISR_clades_for_Bioregions_0.80_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_Bioregions_0.80_df.rds")

# Current richness ~ Binomial GLM with model selection based on AICc
# If crown age: "Time-for-accumulation hypothesis"
# If lambda_0: Cradle hypothesis (higher initial net div?)
# If alpha: Museum hypothesis? (less decline in rates with time?)

### 9.1/ Run unstandardized model ####

RISR_clades_for_aggregated_OW_NW_richness_model_GLM <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                           formula = nb_descendant_tips ~ (node_age + net_div_0 + alpha + 0) * (Hemisphere),
                                                           family = "poisson")
# Add + 0 to remove the Intercept and compute coefficients for each levels => HemisphereOld World = (Intercept) ; HemisphereNew World = HemisphereNew World + (Intercept)
# Yet, still have to add the base-level slope coefficient to the other slope-levels coefficients

summary(RISR_clades_for_aggregated_OW_NW_richness_model_GLM)
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs <- coefficients(RISR_clades_for_aggregated_OW_NW_richness_model_GLM)

# Type I Anova (sequential)
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_anova_table <- anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM, test = 'LRT') ; RISR_clades_for_aggregated_OW_NW_richness_model_GLM_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM, type = "II", test.statistic = "LR")

### 9.2/ Run standardized model ####

## Standardized version (only continuous predictors are scaled)

# reghelper::beta(model = RISR_clades_for_aggregated_OW_NW_richness_model_GLM, skip = "Hemisphere")

RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled <- glm(data = RISR_clades_for_Bioregions_0.80_df,
                                                            formula = nb_descendant_tips ~ (scale(node_age) + scale(net_div_0) + scale(alpha) + 0) * (Hemisphere),
                                                            family = "poisson")

summary(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled)
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std <- coefficients(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled)

# Type I Anova (sequential)
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_anova_table <- anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled, test = 'LRT') ; RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled, type = "II", test.statistic = "LR")

# Compare design matrices
model.matrix(RISR_clades_for_aggregated_OW_NW_richness_model_GLM)
model.matrix(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled)

### 9.3/ Extract Beta coefficients ####

# Extract standardized coefficients data
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_summary_table <- summary(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled)$coefficients
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_summary_table

RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df <- data.frame(par_type = c("Crown age", "Crown net div. rate", "Time variation", rep("Hemisphere", 2), "Crown age", "Crown net div. rate", "Time variation"),
                                                                         hemisphere = c(rep("Old World", 4), rep("New World", 4)),
                                                                         mean = RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_summary_table[,1], 
                                                                         se = RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_summary_table[,2])

# Sort in logical order for plotting
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df %>% 
  arrange(par_type, hemisphere) %>% 
  # Add par_name
  mutate(par_name = paste0(par_type, "_", hemisphere)) %>% 
  dplyr::select(par_name, par_type, hemisphere, mean, se)

# Add base level values to other Bioregion coefficients (except the Bioregion = Intercept)
# For Crown age
base_coef <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"]
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] + base_coef
# For Crown net div. rate
base_coef <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"]
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] + base_coef
# For Time variation
base_coef <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"]
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"] + base_coef

# Use New World as the base level for Hemisphere effects (to be consistent with using Neotropics for Bioregions)
# Substract New World Beta-coefficient to the Hemisphere coefficents
base_coef <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere" & RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"]
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere"] <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere"] - base_coef

# Get multiplier for Z-scores CI
Q2.5_mult <- qnorm(p = 0.025, mean = 0, sd = 1)
Q97.5_mult <- qnorm(p = 0.975, mean = 0, sd = 1)

# Add CI95% intervals
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df <- RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df %>%
  mutate(CI2.5 = mean + Q2.5_mult*se,
         CI97.5 = mean + Q97.5_mult*se) %>%
  dplyr::select(par_name, par_type, hemisphere, mean, se, CI2.5, CI97.5)

# Save std coeffs df
saveRDS(object = RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df, file = "./outputs/RISR/RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df.rds")


### 9.4/ Plot Beta coefficients ####

# Set color scheme for Old World vs. New World
colors_list_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_names <- c("Old World", "New World")
names(colors_list_for_OW_NW) <- OW_NW_names
OW_NW_labels <- c("OW", "NW")

# Order par_type so facets are properly ordered
RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type <- factor(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df$par_type,
                                                                                    levels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"),
                                                                                    labels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"))



## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_richness_model_coefs_plot.pdf", height = 8, width = 8)

RISR_for_aggregated_OW_NW_richness_model_coefs_plot <- ggplot(data = RISR_clades_for_aggregated_OW_NW_richness_model_GLM_coefs_std_df) +
  
  # Add horizontal line: Beta-coef = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Add error bars
  geom_errorbar(mapping = aes(x = hemisphere, y = mean, ymin = CI2.5, ymax = CI97.5, col = hemisphere),
                width = 0.2, linewidth = 1.2,
                alpha = 1.0, show.legend = F) +
  
  # Add mean points
  geom_point(mapping = aes(x = hemisphere, y = mean, fill = hemisphere),
             shape = 21, size = 5, col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust x-axis labels
  scale_x_discrete(labels = c("Afrotropics" = "Afrotrp", "Australasia" = "Austr", "Indomalaya" = "IndoM", "Neotropics" = "Neotrp")) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ par_type,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
  ) +
  
  # Add an asterisk beside Neotropics bioregion coefficient to signify it is used as the reference level
  geom_text(data = data.frame(hemisphere = "New World", mean = 0.22, lab = "*",
                              par_type = factor("Hemisphere", levels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"))),
            aes(x = hemisphere, y = mean, col = hemisphere),
            label = "*", size = 10, nudge_x = 0.22, show.legend = F) +
  
  # Flip axes
  # coord_flip() + 
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste(beta,"-coefficients"))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_richness_model_coefs_plot)

dev.off()


### 9.5/ Compare % of variance (deviance) explained by each factor in each aggregated Old World vs. New World ####

# One model per Hemisphere
OW_NW_levels <- levels(RISR_clades_for_Bioregions_0.80_df$Hemisphere)

## 9.5.1/ Run GLM and extract Anova Type II outputs ####

RISR_clades_for_aggregated_OW_NW_0.80_Anova_df <- data.frame()
for (i in seq_along(OW_NW_levels))
{
  # i <- 1
  
  # Extract hemisphere
  hemisphere_i <- OW_NW_levels[i]
  
  # Extract RISR data
  RISR_clades_for_aggregated_OW_NW_0.80_df_i <- RISR_clades_for_Bioregions_0.80_df %>% 
    filter(Hemisphere == hemisphere_i)
  
  # Run GLM
  RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i <- glm(data = RISR_clades_for_aggregated_OW_NW_0.80_df_i,
                                                                formula = nb_descendant_tips ~ scale(node_age) + scale(net_div_0) + scale(alpha),
                                                                family = "poisson")
  
  summary(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i)
  
  # Type I Anova (sequential)
  # anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i, test = 'LRT')
  
  # Use Type II Anova to obtain marginal deviance
  # Type II Anova (marginal)
  RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_anova_table_i <- car::Anova(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i, type = "II", test.statistic = "LR")
  RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_anova_table_i
  
  # Extract predictor names
  predictors <- str_remove(string = names(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i$model)[-1], pattern = "scale\\(")
  predictors <- str_remove(string = predictors, pattern = "\\)")
  
  # Extract % of marginal deviance explain by each factor
  # Convert to % using the Null Deviance
  RISR_clades_for_aggregated_OW_NW_richness_model_Anofa_df_i <- data.frame(hemisphere = hemisphere_i,
                                                                     predictor = predictors,
                                                                     deviance_perc = round(RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_anova_table_i$`LR Chisq` / RISR_clades_for_aggregated_OW_NW_richness_model_GLM_scaled_i$null.deviance * 100, 1))
  
  # Store in final df
  RISR_clades_for_aggregated_OW_NW_0.80_Anova_df <- rbind(RISR_clades_for_aggregated_OW_NW_0.80_Anova_df, RISR_clades_for_aggregated_OW_NW_richness_model_Anofa_df_i)
}

# Adjust predictor levels
RISR_clades_for_aggregated_OW_NW_0.80_Anova_df$predictor <- factor(x = RISR_clades_for_aggregated_OW_NW_0.80_Anova_df$predictor,
                                                             levels = c("node_age", "net_div_0", "alpha"),
                                                             labels = c("Crown age", "Crown net div. rate", "Time variation"))

# Save Anova type II summary df for aggregated Old World vs. New World
saveRDS(object = RISR_clades_for_aggregated_OW_NW_0.80_Anova_df, file = "./outputs/RISR/RISR_clades_for_aggregated_OW_NW_0.80_Anova_df.rds")


## 9.5.2/ Plot % of marginal deviance ~ predictors for aggregated Old World vs. New World ####

## GGplot
pdf(file = "./outputs/RISR/RISR_for_aggregated_OW_NW_richness_model_Anova_predictors_plot.pdf", height = 6, width = 8)

RISR_for_aggregated_OW_NW_richness_model_Anova_predictors_plot <- ggplot(data = RISR_clades_for_aggregated_OW_NW_0.80_Anova_df) +
  
  # Add column bars
  geom_col(mapping = aes(x = predictor, y = deviance_perc, fill = hemisphere),
           col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust x-axis labels, , 
  scale_x_discrete(labels = c("Crown age" = bquote(Age[crown]), "Crown net div. rate" = bquote(Rate[crown]), "Time variation" = bquote(TV~(alpha)))) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ hemisphere,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
  ) +
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper aggregated OW vs. NW")) +
  
  # Set axes labels
  xlab("Predictors") +
  ylab("Marginal deviance  [%]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed", linewidth = 0.5),
        panel.grid.major.x = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5), angle = 0),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_aggregated_OW_NW_richness_model_Anova_predictors_plot)

dev.off()


##### 10/ Identify factors explaining current richness, for independent Old World vs. New World #####


# Load metadata df of RISR clades per independent OW vs. NW for the selected endemicity threshold
RISR_clades_for_OW_NW_0.95_df <- readRDS(file = "./outputs/RISR/RISR_clades_for_OW_NW_0.95_df.rds")
RISR_clades_for_OW_NW_0.95_df$Hemisphere <- as.factor(str_replace(string = RISR_clades_for_OW_NW_0.95_df$Region, pattern = "_", replacement = " "))

# Current richness ~ Binomial GLM with model selection based on AICc
# If crown age: "Time-for-accumulation hypothesis"
# If lambda_0: Cradle hypothesis (higher initial net div?)
# If alpha: Museum hypothesis? (less decline in rates with time?)

### 10.1/ Run unstandardized model ####

RISR_clades_for_independent_OW_NW_richness_model_GLM <- glm(data = RISR_clades_for_OW_NW_0.95_df,
                                                           formula = nb_descendant_tips ~ (node_age + net_div_0 + alpha + 0) * (Hemisphere),
                                                           family = "poisson")
# Add + 0 to remove the Intercept and compute coefficients for each levels => HemisphereOld World = (Intercept) ; HemisphereNew World = HemisphereNew World + (Intercept)
# Yet, still have to add the base-level slope coefficient to the other slope-levels coefficients

summary(RISR_clades_for_independent_OW_NW_richness_model_GLM)
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs <- coefficients(RISR_clades_for_independent_OW_NW_richness_model_GLM)

# Type I Anova (sequential)
RISR_clades_for_independent_OW_NW_richness_model_GLM_anova_table <- anova(RISR_clades_for_independent_OW_NW_richness_model_GLM, test = 'LRT') ; RISR_clades_for_independent_OW_NW_richness_model_GLM_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_independent_OW_NW_richness_model_GLM, type = "II", test.statistic = "LR")

### 10.2/ Run standardized model ####

## Standardized version (only continuous predictors are scaled)

# reghelper::beta(model = RISR_clades_for_independent_OW_NW_richness_model_GLM, skip = "Hemisphere")

RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled <- glm(data = RISR_clades_for_OW_NW_0.95_df,
                                                                  formula = nb_descendant_tips ~ (scale(node_age) + scale(net_div_0) + scale(alpha) + 0) * (Hemisphere),
                                                                  family = "poisson")

summary(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled)
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std <- coefficients(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled)

# Type I Anova (sequential)
RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_anova_table <- anova(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled, test = 'LRT') ; RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_anova_table
# Type II Anova (marginal)
car::Anova(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled, type = "II", test.statistic = "LR")

# Compare design matrices
model.matrix(RISR_clades_for_independent_OW_NW_richness_model_GLM)
model.matrix(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled)

### 10.3/ Extract Beta coefficients ####

# Extract standardized coefficients data
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_summary_table <- summary(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled)$coefficients
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_summary_table

RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df <- data.frame(par_type = c("Crown age", "Crown net div. rate", "Time variation", rep("Hemisphere", 2), "Crown age", "Crown net div. rate", "Time variation"),
                                                                               hemisphere = c(rep("New World", 4), rep("Old World", 4)),
                                                                               mean = RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_summary_table[,1], 
                                                                               se = RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_summary_table[,2])

# Sort in logical order for plotting
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df %>% 
  arrange(par_type, hemisphere) %>% 
  # Add par_name
  mutate(par_name = paste0(par_type, "_", hemisphere)) %>% 
  dplyr::select(par_name, par_type, hemisphere, mean, se)

# Add base level values to other Bioregion coefficients (except the Bioregion = Intercept)
# For Crown age
base_coef <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"]
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown age" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] + base_coef
# For Crown net div. rate
base_coef <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"]
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Crown net div. rate" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] + base_coef
# For Time variation
base_coef <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"]
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Time variation" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "Old World"] + base_coef

# Use New World as the base level for Hemisphere effects (to be consistent with using Neotropics for Bioregions)
# Substract New World Beta-coefficient to the Hemisphere coefficents
base_coef <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere" & RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$hemisphere == "New World"]
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere"] <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$mean[RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type == "Hemisphere"] - base_coef

# Get multiplier for Z-scores CI
Q2.5_mult <- qnorm(p = 0.025, mean = 0, sd = 1)
Q97.5_mult <- qnorm(p = 0.975, mean = 0, sd = 1)

# Add CI95% intervals
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df <- RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df %>%
  mutate(CI2.5 = mean + Q2.5_mult*se,
         CI97.5 = mean + Q97.5_mult*se) %>%
  dplyr::select(par_name, par_type, hemisphere, mean, se, CI2.5, CI97.5)

# Save std coeffs df
saveRDS(object = RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df, file = "./outputs/RISR/RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df.rds")


### 10.4/ Plot Beta coefficients ####

# Set color scheme for Old World vs. New World
colors_list_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_names <- c("Old World", "New World")
names(colors_list_for_OW_NW) <- OW_NW_names
OW_NW_labels <- c("OW", "NW")

# Order par_type so facets are properly ordered
RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type <- factor(RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df$par_type,
                                                                                    levels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"),
                                                                                    labels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"))



## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_richness_model_coefs_plot.pdf", height = 8, width = 8)

RISR_for_independent_OW_NW_richness_model_coefs_plot <- ggplot(data = RISR_clades_for_independent_OW_NW_richness_model_GLM_coefs_std_df) +
  
  # Add horizontal line: Beta-coef = 0
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
  
  # Add error bars
  geom_errorbar(mapping = aes(x = hemisphere, y = mean, ymin = CI2.5, ymax = CI97.5, col = hemisphere),
                width = 0.2, linewidth = 1.2,
                alpha = 1.0, show.legend = F) +
  
  # Add mean points
  geom_point(mapping = aes(x = hemisphere, y = mean, fill = hemisphere),
             shape = 21, size = 5, col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust x-axis labels
  scale_x_discrete(labels = c("Afrotropics" = "Afrotrp", "Australasia" = "Austr", "Indomalaya" = "IndoM", "Neotropics" = "Neotrp")) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ par_type,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
  ) +
  
  # Add an asterisk beside Neotropics bioregion coefficient to signify it is used as the reference level
  geom_text(data = data.frame(hemisphere = "New World", mean = 0.22, lab = "*",
                              par_type = factor("Hemisphere", levels = c("Hemisphere", "Crown age", "Crown net div. rate", "Time variation"))),
            aes(x = hemisphere, y = mean, col = hemisphere),
            label = "*", size = 10, nudge_x = 0.22, show.legend = F) +
  
  # Flip axes
  # coord_flip() + 
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab(expression(paste(beta,"-coefficients"))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_richness_model_coefs_plot)

dev.off()


### 10.5/ Compare % of variance (deviance) explained by each factor in each independent Old World vs. New World ####

# One model per Hemisphere
OW_NW_levels <- levels(RISR_clades_for_OW_NW_0.95_df$Hemisphere)

## 10.5.1/ Run GLM and extract Anova Type II outputs ####

RISR_clades_for_independent_OW_NW_0.80_Anova_df <- data.frame()
for (i in seq_along(OW_NW_levels))
{
  # i <- 1
  
  # Extract hemisphere
  hemisphere_i <- OW_NW_levels[i]
  
  # Extract RISR data
  RISR_clades_for_independent_OW_NW_0.80_df_i <- RISR_clades_for_OW_NW_0.95_df %>% 
    filter(Hemisphere == hemisphere_i)
  
  # Run GLM
  RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i <- glm(data = RISR_clades_for_independent_OW_NW_0.80_df_i,
                                                                      formula = nb_descendant_tips ~ scale(node_age) + scale(net_div_0) + scale(alpha),
                                                                      family = "poisson")
  
  summary(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i)
  
  # Type I Anova (sequential)
  # anova(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i, test = 'LRT')
  
  # Use Type II Anova to obtain marginal deviance
  # Type II Anova (marginal)
  RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_anova_table_i <- car::Anova(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i, type = "II", test.statistic = "LR")
  RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_anova_table_i
  
  # Extract predictor names
  predictors <- str_remove(string = names(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i$model)[-1], pattern = "scale\\(")
  predictors <- str_remove(string = predictors, pattern = "\\)")
  
  # Extract % of marginal deviance explain by each factor
  # Convert to % using the Null Deviance
  RISR_clades_for_independent_OW_NW_richness_model_Anofa_df_i <- data.frame(hemisphere = hemisphere_i,
                                                                           predictor = predictors,
                                                                           deviance_perc = round(RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_anova_table_i$`LR Chisq` / RISR_clades_for_independent_OW_NW_richness_model_GLM_scaled_i$null.deviance * 100, 1))
  
  # Store in final df
  RISR_clades_for_independent_OW_NW_0.80_Anova_df <- rbind(RISR_clades_for_independent_OW_NW_0.80_Anova_df, RISR_clades_for_independent_OW_NW_richness_model_Anofa_df_i)
}

# Adjust predictor levels
RISR_clades_for_independent_OW_NW_0.80_Anova_df$predictor <- factor(x = RISR_clades_for_independent_OW_NW_0.80_Anova_df$predictor,
                                                                   levels = c("node_age", "net_div_0", "alpha"),
                                                                   labels = c("Crown age", "Crown net div. rate", "Time variation"))

# Save Anova type II summary df for independent Old World vs. New World
saveRDS(object = RISR_clades_for_independent_OW_NW_0.80_Anova_df, file = "./outputs/RISR/RISR_clades_for_independent_OW_NW_0.80_Anova_df.rds")


## 10.5.2/ Plot % of marginal deviance ~ predictors for independent Old World vs. New World ####

## GGplot
pdf(file = "./outputs/RISR/RISR_for_independent_OW_NW_richness_model_Anova_predictors_plot.pdf", height = 6, width = 8)

RISR_for_independent_OW_NW_richness_model_Anova_predictors_plot <- ggplot(data = RISR_clades_for_independent_OW_NW_0.80_Anova_df) +
  
  # Add column bars
  geom_col(mapping = aes(x = predictor, y = deviance_perc, fill = hemisphere),
           col = "black", alpha = 1.0, show.legend = F) +
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregion", breaks = names(colors_list_for_OW_NW), labels = OW_NW_labels, values = colors_list_for_OW_NW) +
  
  # Adjust x-axis labels, , 
  scale_x_discrete(labels = c("Crown age" = bquote(Age[crown]), "Crown net div. rate" = bquote(Rate[crown]), "Time variation" = bquote(TV~(alpha)))) +
  
  # Make a facetted plot
  facet_wrap(facets = ~ hemisphere,
             scales = "fixed",
             # scales = "free",
             nrow = 2, ncol = 2,
             # nrow = 4, ncol = 1
  ) +
  
  # Set plot title +
  ggtitle(label = paste0("Predictors of current richness in RISR clades\nper independent OW vs. NW")) +
  
  # Set axes labels
  xlab("Predictors") +
  ylab("Marginal deviance  [%]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed", linewidth = 0.5),
        panel.grid.major.x = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        strip.text.x = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5), angle = 0),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(RISR_for_independent_OW_NW_richness_model_Anova_predictors_plot)

dev.off()



### Bonus in other scripts (?)

# For a more smooth model with many small cladogenetic changes, see CLaDS
# GeoSSE, GeoHiSSE for Old World vs. New world. See Kawahara ?




