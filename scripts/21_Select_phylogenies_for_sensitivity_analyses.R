##### Script 21: Select alternative dating hypothesis for Sensitivity Analyses  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Select the two most extreme hypotheses for divergence dating to use for Sensitivity Analyses
# Those two grafted time-calibrated phylogenies will be used to run the whoel analytic pipeline and
# show robustness of our results to dating uncertainty

###

### Inputs

# Set of non-grafted time-calibrated backbone phylogenies
# Associated set of grafted time-calibrated phylogenies

###

### Outputs

# Youngest_phylogeny = the phylogeny with the youngest divergence dating hypothesis for backbone nodes
# Oldest_phylogeny = the phylogeny with the oldest divergence dating hypothesis for backbone nodes

###

# Clean environment
rm(list = ls())


##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(phytools)
library(deeptime)  # Library to display geological scale on plot
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges


### 1.2/ Load non-grafted time-calibrated backbone phylogenies ####

Ponerinae_all_posteriors_phylogeny_789t <- readRDS(file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.rds")

Ponerinae_MCC_phylogeny_789t_treedata <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")

### 1.3/ Load associated grafted time-calibrated phylogenies ####

Ponerinae_all_posteriors_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.rds")


##### 2/ Identify most extreme hypotheses #####

### 2.1/ Extract sum of ages ###

root_ages <- c()
ages_sums <- c()
for (i in seq_along(Ponerinae_all_posteriors_phylogeny_1534t))
{
  # i <- 1
  
  # Extract phylo
  phy_i <- Ponerinae_all_posteriors_phylogeny_1534t[[i]]
  # Get root age
  root_ages[i] <- max(phytools::nodeHeights(phy_i)[,2])
  # Get node ages
  node_ages_df_i <- round(root_ages[i] - phytools::nodeHeights(phy_i), 5)
  # Get sum of internal node ages
  ages_sums[i] <- sum(node_ages_df_i[,1])
}

hist(root_ages)  
summary(root_ages)
quantile(root_ages, probs = c(0.025, 0.975))
hist(ages_sums)
summary(ages_sums)

### 2.2/ Identify the most extreme posterior trees ####

Youngest_ID <- which(rank(ages_sums) == 1) # Youngest
Young_Q2.5_ID <- which(rank(ages_sums) == 25) # Young Q2.5% 
Oldest_ID <- which(rank(ages_sums) == 1000) # Oldest
Old_Q97.5_ID <-which(rank(ages_sums) == 975) # Old Q97.5% 

root_ages[c(Youngest_ID, Young_Q2.5_ID, Old_Q97.5_ID, Oldest_ID)]
ages_sums[c(Youngest_ID, Young_Q2.5_ID, Old_Q97.5_ID, Oldest_ID)]

##### 3/ Export the most extreme hypotheses #####

### 3.1/ Extract from backbone trees ###

Ponerinae_Youngest_phylogeny_789t <- Ponerinae_all_posteriors_phylogeny_789t[[Youngest_ID]]
Ponerinae_Oldest_phylogeny_789t <- Ponerinae_all_posteriors_phylogeny_789t[[Oldest_ID]]

### 3.2/ Extract from grafted trees ###

Ponerinae_Youngest_phylogeny_1534t <- Ponerinae_all_posteriors_phylogeny_1534t[[Youngest_ID]]
Ponerinae_Oldest_phylogeny_1534t <- Ponerinae_all_posteriors_phylogeny_1534t[[Oldest_ID]]

### 3.3/ Add node metadata ####

## Load grafted MCC phylogeny with node metadata
Ponerinae_MCC_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")

View(Ponerinae_MCC_phylogeny_1534t_treedata@data)
str(Ponerinae_MCC_phylogeny_1534t_treedata@data$node, 1)

## Need to find correspondence between previous backbone nodes and current nodes using matching descendant tips and age

## 3.3.1/ For Ponerinae_Youngest_phylogeny_1534t ####

## Initiate new treedata object for the grafted tree
Ponerinae_Youngest_phylogeny_1534t_treedata <- Ponerinae_MCC_phylogeny_789t_treedata
Ponerinae_Youngest_phylogeny_1534t_treedata@phylo <- Ponerinae_Youngest_phylogeny_1534t # Replace phylogeny
Ponerinae_Youngest_phylogeny_1534t_treedata@data$backbone_node <- Ponerinae_Youngest_phylogeny_1534t_treedata@data$node # Save backbone node ID
Ponerinae_Youngest_phylogeny_1534t_treedata@data$node <- NA # Initiate new node ID

# Get node ages for backbone tree
root_age_789t <- max(phytools::nodeHeights(Ponerinae_Youngest_phylogeny_789t))
Ponerinae_Youngest_phylogeny_789t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_Youngest_phylogeny_789t) + root_age_789t, 5))
names(Ponerinae_Youngest_phylogeny_789t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_Youngest_phylogeny_789t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_Youngest_phylogeny_789t_edges_ages_df)
Ponerinae_Youngest_phylogeny_789t_edges_ages_df <- Ponerinae_Youngest_phylogeny_789t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_Youngest_phylogeny_789t_edges_ages_df <- cbind(Ponerinae_Youngest_phylogeny_789t_edges_ages_df, Ponerinae_Youngest_phylogeny_789t$edge)
names(Ponerinae_Youngest_phylogeny_789t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_Youngest_phylogeny_789t_nodes_ages_df <- Ponerinae_Youngest_phylogeny_789t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_Youngest_phylogeny_789t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_Youngest_phylogeny_789t$tip.label),
                                                        node_age = 0)
Ponerinae_Youngest_phylogeny_789t_nodes_ages_df <- rbind(Ponerinae_Youngest_phylogeny_789t_tips_ages_df, Ponerinae_Youngest_phylogeny_789t_nodes_ages_df)

# Get node ages for grafted tree
root_age_1534t <- max(phytools::nodeHeights(Ponerinae_Youngest_phylogeny_1534t))
Ponerinae_Youngest_phylogeny_1534t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_Youngest_phylogeny_1534t) + root_age_1534t, 5))
names(Ponerinae_Youngest_phylogeny_1534t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_Youngest_phylogeny_1534t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_Youngest_phylogeny_1534t_edges_ages_df)
Ponerinae_Youngest_phylogeny_1534t_edges_ages_df <- Ponerinae_Youngest_phylogeny_1534t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_Youngest_phylogeny_1534t_edges_ages_df <- cbind(Ponerinae_Youngest_phylogeny_1534t_edges_ages_df, Ponerinae_Youngest_phylogeny_1534t$edge)
names(Ponerinae_Youngest_phylogeny_1534t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_Youngest_phylogeny_1534t_nodes_ages_df <- Ponerinae_Youngest_phylogeny_1534t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_Youngest_phylogeny_1534t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_Youngest_phylogeny_1534t$tip.label),
                                                         node_age = 0)
Ponerinae_Youngest_phylogeny_1534t_nodes_ages_df <- rbind(Ponerinae_Youngest_phylogeny_1534t_tips_ages_df, Ponerinae_Youngest_phylogeny_1534t_nodes_ages_df)

# Extract descendants for all nodes in backbone tree
Ponerinae_Youngest_phylogeny_789t_descendants_list <- list()
nb_nodes_789t <- dim(Ponerinae_MCC_phylogeny_789t_treedata@data)[1]
for (i in 1:nb_nodes_789t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_Youngest_phylogeny_789t, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Youngest_phylogeny_789t$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_Youngest_phylogeny_789t$tip.label[all_descendants_i]
  
  # Store in list
  Ponerinae_Youngest_phylogeny_789t_descendants_list[[i]] <- all_descendants_i
}

# Extract descendants for all nodes in grafted tree
Ponerinae_Youngest_phylogeny_1534t_descendants_list <- list()
nb_nodes_1534t <- Ponerinae_Youngest_phylogeny_1534t$Nnode + length(Ponerinae_Youngest_phylogeny_1534t$tip.label)
for (i in 1:nb_nodes_1534t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_Youngest_phylogeny_1534t, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Youngest_phylogeny_1534t$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_Youngest_phylogeny_1534t$tip.label[all_descendants_i]
  # Keep only tips from the backbone
  all_descendants_i <- all_descendants_i[all_descendants_i %in% Ponerinae_Youngest_phylogeny_789t$tip.label]
  # Reorder alphabetically
  all_descendants_i <- all_descendants_i[order(all_descendants_i)]
  
  # Store in list
  Ponerinae_Youngest_phylogeny_1534t_descendants_list[[i]] <- all_descendants_i
}

# Find corresponding node in the grafted tree based on list of descending tips and node age
for (i in 1:nb_nodes_789t)
{
  # Extract descending tips from backbone node
  all_descendants_backbone_i <- Ponerinae_Youngest_phylogeny_789t_descendants_list[[i]]
  # Reorder alphabetically
  all_descendants_backbone_i <- all_descendants_backbone_i[order(all_descendants_backbone_i)]
  # Find the matches among list of descendants from grafted tree
  descendant_match_test_i <- unlist(lapply(X = Ponerinae_Youngest_phylogeny_1534t_descendants_list, FUN = function (x) { identical(x, all_descendants_backbone_i) }))
  descendant_match_nodes_ID <- which(descendant_match_test_i)
  
  # Extract node age in the backbone tree
  focal_node_age <- Ponerinae_Youngest_phylogeny_789t_nodes_ages_df$node_age[i]
  # Extract node ages of the matching nodes in the grafted tree
  descendant_match_nodes_ages <- Ponerinae_Youngest_phylogeny_1534t_nodes_ages_df$node_age[descendant_match_nodes_ID]
  # Find the matching node age
  age_matching_node_ID <- descendant_match_nodes_ID[abs(focal_node_age - descendant_match_nodes_ages) < 0.001]
  
  # Store matching ID in treedata object
  Ponerinae_Youngest_phylogeny_1534t_treedata@data$node[i] <- age_matching_node_ID
}

# Create template for all nodes (height = age, length of parental edge, posterior (tips = NA, non-tips = 1))
Ponerinae_Youngest_phylogeny_1534t_nodes_data_df <- Ponerinae_Youngest_phylogeny_1534t_nodes_ages_df %>% 
  rename(node = node_ID,
         height = node_age)
Ponerinae_Youngest_phylogeny_1534t_nodes_data_df$posterior <- NA
Ponerinae_Youngest_phylogeny_1534t_nodes_data_df$posterior[Ponerinae_Youngest_phylogeny_1534t_nodes_data_df$node > length(Ponerinae_Youngest_phylogeny_1534t$tip.label)] <- 1
Ponerinae_Youngest_phylogeny_1534t_nodes_data_df$length <- Ponerinae_Youngest_phylogeny_1534t$edge.length[match(x = Ponerinae_Youngest_phylogeny_1534t_nodes_data_df$node, table = Ponerinae_Youngest_phylogeny_1534t$edge[,2])]

# Join backbone data with template data
Ponerinae_Youngest_phylogeny_1534t_treedata@data <- left_join(Ponerinae_Youngest_phylogeny_1534t_nodes_data_df, Ponerinae_Youngest_phylogeny_1534t_treedata@data, by = "node") %>% 
  rename(height = height.x) %>% # Keep the height/age as measured in the grafted tree
  select(-height.y) %>% # Remove height as measured in the backbone tree (should be equal modulo numerical precision)
  rename(posterior = posterior.x) %>% # Keep info on non-tips from the grafted tree
  select(-posterior.y) %>% # Remove info on non-tips from the backbone tree as some nodes are missing
  rename(length = length.x) %>% # Keep parental branch length from the grafted tree as the current length
  rename(length_backbone = length.y) %>% # Record parental branch length from the backbone tree (may be different due to grafting)
  rename(length_median_backbone = length_median,
         length_0.95_HPD_backbone = length_0.95_HPD,
         length_range_backbone = length_range,
         height_median_backbone = height_median,
         height_0.95_HPD_backbone = height_0.95_HPD,
         height_range_backbone = height_range) %>% 
  mutate(new_node = is.na(backbone_node)) %>%
  select(node, backbone_node, posterior, height, height_median_backbone, height_0.95_HPD_backbone, height_range_backbone, length, length_backbone, length_median_backbone, length_0.95_HPD_backbone, length_range_backbone, new_node) %>% 
  as_tibble()

# Record nodes with new parental edge
all_missing_taxa_ID <- which(Ponerinae_Youngest_phylogeny_1534t_treedata@data$new_node)
all_missing_taxa_ID <- all_missing_taxa_ID[all_missing_taxa_ID <= length(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label)]
Ponerinae_Youngest_phylogeny_1534t_treedata@data$new_parental_edge <- FALSE

for (i in 1:nrow(Ponerinae_Youngest_phylogeny_1534t_treedata@data))
{
  # i <- 1
  
  # Extract node ID
  node_ID_i <- Ponerinae_Youngest_phylogeny_1534t_treedata@data$node[i]
  
  # Extract all current descendants
  all_descendants_i <- phytools::getDescendants(tree = Ponerinae_Youngest_phylogeny_1534t_treedata@phylo, node = node_ID_i)
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label)]
  
  # Check if all current descendants are missing taxa
  missing_i <- all(all_descendants_i %in% all_missing_taxa_ID)
  
  # Record status
  Ponerinae_Youngest_phylogeny_1534t_treedata@data$new_parental_edge[i] <- missing_i
}
table(Ponerinae_Youngest_phylogeny_1534t_treedata@data$new_parental_edge)

# Demultiplex the range values for node ages
Ponerinae_Youngest_phylogeny_1534t_treedata@data$age_min_backbone <- unlist(lapply(X = Ponerinae_Youngest_phylogeny_1534t_treedata@data$height_range_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[1]} } ))
Ponerinae_Youngest_phylogeny_1534t_treedata@data$age_max_backbone <- unlist(lapply(X = Ponerinae_Youngest_phylogeny_1534t_treedata@data$height_range_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[2]} } ))

# Demultiplex the 95% HPD range values for node ages
Ponerinae_Youngest_phylogeny_1534t_treedata@data$age_HPD_0.025_backbone <- unlist(lapply(X = Ponerinae_Youngest_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[1]} } ))
Ponerinae_Youngest_phylogeny_1534t_treedata@data$age_HPD_0.975_backbone <- unlist(lapply(X = Ponerinae_Youngest_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[2]} } ))

View(Ponerinae_Youngest_phylogeny_1534t_treedata@data)
View(Ponerinae_MCC_phylogeny_789t_treedata@data)


## 3.3.2/ For Ponerinae_Oldest_phylogeny_1534t ####

## Initiate new treedata object for the grafted tree
Ponerinae_Oldest_phylogeny_1534t_treedata <- Ponerinae_MCC_phylogeny_789t_treedata
Ponerinae_Oldest_phylogeny_1534t_treedata@phylo <- Ponerinae_Oldest_phylogeny_1534t # Replace phylogeny
Ponerinae_Oldest_phylogeny_1534t_treedata@data$backbone_node <- Ponerinae_Oldest_phylogeny_1534t_treedata@data$node # Save backbone node ID
Ponerinae_Oldest_phylogeny_1534t_treedata@data$node <- NA # Initiate new node ID

# Get node ages for backbone tree
root_age_789t <- max(phytools::nodeHeights(Ponerinae_Oldest_phylogeny_789t))
Ponerinae_Oldest_phylogeny_789t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_Oldest_phylogeny_789t) + root_age_789t, 5))
names(Ponerinae_Oldest_phylogeny_789t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_Oldest_phylogeny_789t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_Oldest_phylogeny_789t_edges_ages_df)
Ponerinae_Oldest_phylogeny_789t_edges_ages_df <- Ponerinae_Oldest_phylogeny_789t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_Oldest_phylogeny_789t_edges_ages_df <- cbind(Ponerinae_Oldest_phylogeny_789t_edges_ages_df, Ponerinae_Oldest_phylogeny_789t$edge)
names(Ponerinae_Oldest_phylogeny_789t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_Oldest_phylogeny_789t_nodes_ages_df <- Ponerinae_Oldest_phylogeny_789t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_Oldest_phylogeny_789t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_Oldest_phylogeny_789t$tip.label),
                                                          node_age = 0)
Ponerinae_Oldest_phylogeny_789t_nodes_ages_df <- rbind(Ponerinae_Oldest_phylogeny_789t_tips_ages_df, Ponerinae_Oldest_phylogeny_789t_nodes_ages_df)

# Get node ages for grafted tree
root_age_1534t <- max(phytools::nodeHeights(Ponerinae_Oldest_phylogeny_1534t))
Ponerinae_Oldest_phylogeny_1534t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_Oldest_phylogeny_1534t) + root_age_1534t, 5))
names(Ponerinae_Oldest_phylogeny_1534t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_Oldest_phylogeny_1534t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_Oldest_phylogeny_1534t_edges_ages_df)
Ponerinae_Oldest_phylogeny_1534t_edges_ages_df <- Ponerinae_Oldest_phylogeny_1534t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_Oldest_phylogeny_1534t_edges_ages_df <- cbind(Ponerinae_Oldest_phylogeny_1534t_edges_ages_df, Ponerinae_Oldest_phylogeny_1534t$edge)
names(Ponerinae_Oldest_phylogeny_1534t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_Oldest_phylogeny_1534t_nodes_ages_df <- Ponerinae_Oldest_phylogeny_1534t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_Oldest_phylogeny_1534t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_Oldest_phylogeny_1534t$tip.label),
                                                           node_age = 0)
Ponerinae_Oldest_phylogeny_1534t_nodes_ages_df <- rbind(Ponerinae_Oldest_phylogeny_1534t_tips_ages_df, Ponerinae_Oldest_phylogeny_1534t_nodes_ages_df)

# Extract descendants for all nodes in backbone tree
Ponerinae_Oldest_phylogeny_789t_descendants_list <- list()
nb_nodes_789t <- dim(Ponerinae_MCC_phylogeny_789t_treedata@data)[1]
for (i in 1:nb_nodes_789t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_Oldest_phylogeny_789t, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Oldest_phylogeny_789t$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_Oldest_phylogeny_789t$tip.label[all_descendants_i]
  
  # Store in list
  Ponerinae_Oldest_phylogeny_789t_descendants_list[[i]] <- all_descendants_i
}

# Extract descendants for all nodes in grafted tree
Ponerinae_Oldest_phylogeny_1534t_descendants_list <- list()
nb_nodes_1534t <- Ponerinae_Oldest_phylogeny_1534t$Nnode + length(Ponerinae_Oldest_phylogeny_1534t$tip.label)
for (i in 1:nb_nodes_1534t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_Oldest_phylogeny_1534t, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Oldest_phylogeny_1534t$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_Oldest_phylogeny_1534t$tip.label[all_descendants_i]
  # Keep only tips from the backbone
  all_descendants_i <- all_descendants_i[all_descendants_i %in% Ponerinae_Oldest_phylogeny_789t$tip.label]
  # Reorder alphabetically
  all_descendants_i <- all_descendants_i[order(all_descendants_i)]
  
  # Store in list
  Ponerinae_Oldest_phylogeny_1534t_descendants_list[[i]] <- all_descendants_i
}

# Find corresponding node in the grafted tree based on list of descending tips and node age
for (i in 1:nb_nodes_789t)
{
  # Extract descending tips from backbone node
  all_descendants_backbone_i <- Ponerinae_Oldest_phylogeny_789t_descendants_list[[i]]
  # Reorder alphabetically
  all_descendants_backbone_i <- all_descendants_backbone_i[order(all_descendants_backbone_i)]
  # Find the matches among list of descendants from grafted tree
  descendant_match_test_i <- unlist(lapply(X = Ponerinae_Oldest_phylogeny_1534t_descendants_list, FUN = function (x) { identical(x, all_descendants_backbone_i) }))
  descendant_match_nodes_ID <- which(descendant_match_test_i)
  
  # Extract node age in the backbone tree
  focal_node_age <- Ponerinae_Oldest_phylogeny_789t_nodes_ages_df$node_age[i]
  # Extract node ages of the matching nodes in the grafted tree
  descendant_match_nodes_ages <- Ponerinae_Oldest_phylogeny_1534t_nodes_ages_df$node_age[descendant_match_nodes_ID]
  # Find the matching node age
  age_matching_node_ID <- descendant_match_nodes_ID[abs(focal_node_age - descendant_match_nodes_ages) < 0.001]
  
  # Store matching ID in treedata object
  Ponerinae_Oldest_phylogeny_1534t_treedata@data$node[i] <- age_matching_node_ID
}

# Create template for all nodes (height = age, length of parental edge, posterior (tips = NA, non-tips = 1))
Ponerinae_Oldest_phylogeny_1534t_nodes_data_df <- Ponerinae_Oldest_phylogeny_1534t_nodes_ages_df %>% 
  rename(node = node_ID,
         height = node_age)
Ponerinae_Oldest_phylogeny_1534t_nodes_data_df$posterior <- NA
Ponerinae_Oldest_phylogeny_1534t_nodes_data_df$posterior[Ponerinae_Oldest_phylogeny_1534t_nodes_data_df$node > length(Ponerinae_Oldest_phylogeny_1534t$tip.label)] <- 1
Ponerinae_Oldest_phylogeny_1534t_nodes_data_df$length <- Ponerinae_Oldest_phylogeny_1534t$edge.length[match(x = Ponerinae_Oldest_phylogeny_1534t_nodes_data_df$node, table = Ponerinae_Oldest_phylogeny_1534t$edge[,2])]

# Join backbone data with template data
Ponerinae_Oldest_phylogeny_1534t_treedata@data <- left_join(Ponerinae_Oldest_phylogeny_1534t_nodes_data_df, Ponerinae_Oldest_phylogeny_1534t_treedata@data, by = "node") %>% 
  rename(height = height.x) %>% # Keep the height/age as measured in the grafted tree
  select(-height.y) %>% # Remove height as measured in the backbone tree (should be equal modulo numerical precision)
  rename(posterior = posterior.x) %>% # Keep info on non-tips from the grafted tree
  select(-posterior.y) %>% # Remove info on non-tips from the backbone tree as some nodes are missing
  rename(length = length.x) %>% # Keep parental branch length from the grafted tree as the current length
  rename(length_backbone = length.y) %>% # Record parental branch length from the backbone tree (may be different due to grafting)
  rename(length_median_backbone = length_median,
         length_0.95_HPD_backbone = length_0.95_HPD,
         length_range_backbone = length_range,
         height_median_backbone = height_median,
         height_0.95_HPD_backbone = height_0.95_HPD,
         height_range_backbone = height_range) %>% 
  mutate(new_node = is.na(backbone_node)) %>%
  select(node, backbone_node, posterior, height, height_median_backbone, height_0.95_HPD_backbone, height_range_backbone, length, length_backbone, length_median_backbone, length_0.95_HPD_backbone, length_range_backbone, new_node) %>% 
  as_tibble()

# Record nodes with new parental edge
all_missing_taxa_ID <- which(Ponerinae_Oldest_phylogeny_1534t_treedata@data$new_node)
all_missing_taxa_ID <- all_missing_taxa_ID[all_missing_taxa_ID <= length(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label)]
Ponerinae_Oldest_phylogeny_1534t_treedata@data$new_parental_edge <- FALSE

for (i in 1:nrow(Ponerinae_Oldest_phylogeny_1534t_treedata@data))
{
  # i <- 1
  
  # Extract node ID
  node_ID_i <- Ponerinae_Oldest_phylogeny_1534t_treedata@data$node[i]
  
  # Extract all current descendants
  all_descendants_i <- phytools::getDescendants(tree = Ponerinae_Oldest_phylogeny_1534t_treedata@phylo, node = node_ID_i)
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label)]
  
  # Check if all current descendants are missing taxa
  missing_i <- all(all_descendants_i %in% all_missing_taxa_ID)
  
  # Record status
  Ponerinae_Oldest_phylogeny_1534t_treedata@data$new_parental_edge[i] <- missing_i
}
table(Ponerinae_Oldest_phylogeny_1534t_treedata@data$new_parental_edge)

# Demultiplex the range values for node ages
Ponerinae_Oldest_phylogeny_1534t_treedata@data$age_min_backbone <- unlist(lapply(X = Ponerinae_Oldest_phylogeny_1534t_treedata@data$height_range_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[1]} } ))
Ponerinae_Oldest_phylogeny_1534t_treedata@data$age_max_backbone <- unlist(lapply(X = Ponerinae_Oldest_phylogeny_1534t_treedata@data$height_range_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[2]} } ))

# Demultiplex the 95% HPD range values for node ages
Ponerinae_Oldest_phylogeny_1534t_treedata@data$age_HPD_0.025_backbone <- unlist(lapply(X = Ponerinae_Oldest_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[1]} } ))
Ponerinae_Oldest_phylogeny_1534t_treedata@data$age_HPD_0.975_backbone <- unlist(lapply(X = Ponerinae_Oldest_phylogeny_1534t_treedata@data$height_0.95_HPD_backbone, FUN = function (x) { if (is.null(x)) {x <- NA} else {x[2]} } ))

View(Ponerinae_Oldest_phylogeny_1534t_treedata@data)
View(Ponerinae_MCC_phylogeny_789t_treedata@data)

### 3.4/ Convert labels to short names without underscores ####

Ponerinae_MCC_phylogeny_789t_treedata_short_names <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata_short_names.rds")
Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label

# Tips are ordered similarly as in MCC tree
View(cbind(Ponerinae_Youngest_phylogeny_789t$tip.label, Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label))

# Get short names without underscores
short_tip_names_no_underscore_Genera <- str_split(Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label, pattern = "_", simplify = T)[,1]
short_tip_names_no_underscore_df <- str_split(Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label, pattern = "_", simplify = T)[,-1]
short_tip_names_no_underscore_epithets <- apply(X = short_tip_names_no_underscore_df, MARGIN = 1, FUN = paste0, collapse = "_")
short_tip_names_no_underscore_epithets <- str_remove(string = short_tip_names_no_underscore_epithets, pattern = "_$")
short_tip_names_no_underscore_epithets <- str_remove(string = short_tip_names_no_underscore_epithets, pattern = "_$")
short_tip_names_no_underscore_789t <- paste(short_tip_names_no_underscore_Genera, short_tip_names_no_underscore_epithets)

short_tip_names_no_underscore_789t[short_tip_names_no_underscore_789t == "NewGenus bucki"] <- "Neoponera bucki"

## Update tip labels for non-grafted phylogenies
Ponerinae_Youngest_phylogeny_789t$tip.label <- short_tip_names_no_underscore_789t
Ponerinae_Oldest_phylogeny_789t$tip.label <- short_tip_names_no_underscore_789t

## Update tip labels for grafted phylogenies

# Remove underscore from full names
Youngest_short_tip_names_no_underscore_1534t <- str_replace(string = Ponerinae_Youngest_phylogeny_1534t$tip.label, pattern = "_", replacement = " ")
# Replace backbone tips with short names
Youngest_short_tip_names_no_underscore_1534t[Ponerinae_Youngest_phylogeny_1534t$tip.label %in% Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label] <- short_tip_names_no_underscore

# Remove underscore from full names
Oldest_short_tip_names_no_underscore_1534t <- str_replace(string = Ponerinae_Oldest_phylogeny_1534t$tip.label, pattern = "_", replacement = " ")
# Replace backbone tips with short names
Oldest_short_tip_names_no_underscore_1534t[Ponerinae_Oldest_phylogeny_1534t$tip.label %in% Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label] <- short_tip_names_no_underscore

cbind(Ponerinae_Youngest_phylogeny_1534t$tip.label, Youngest_short_tip_names_no_underscore_1534t)
cbind(Ponerinae_Oldest_phylogeny_1534t$tip.label, Oldest_short_tip_names_no_underscore_1534t)

Ponerinae_Youngest_phylogeny_1534t$tip.label <- Youngest_short_tip_names_no_underscore_1534t
Ponerinae_Oldest_phylogeny_1534t$tip.label <- Oldest_short_tip_names_no_underscore_1534t

Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label <- Youngest_short_tip_names_no_underscore_1534t
Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label <- Oldest_short_tip_names_no_underscore_1534t

### 3.5/ Export the trees ####

## Save non-grafted phylogenies
saveRDS(object = Ponerinae_Youngest_phylogeny_789t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Youngest_789t.rds")
saveRDS(object = Ponerinae_Oldest_phylogeny_789t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Oldest_789t.rds")

## Save grafted phylogenies
saveRDS(object = Ponerinae_Youngest_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Youngest_1534t.rds")
saveRDS(object = Ponerinae_Oldest_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Oldest_1534t.rds")

## Save grafted phylogenies with node metadata
saveRDS(object = Ponerinae_Youngest_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Youngest_1534t_treedata.rds")
saveRDS(object = Ponerinae_Oldest_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Oldest_1534t_treedata.rds")


### 3.6/ Plot the trees ####

# Plot nice time-calibrated phylogenies of publication
# Rectangular phylo with 95% HPD and grafted edge colors
# Include genus-groups, geological epochs, ant images

## 3.6.1/ Prepare geological epochs metadata ####

# Use seven time slices associated to geological periods as in Kawahara et al., 2023: 
# 0-5.33 My = Holocene + Pleistocene + Pliocene
# 5.33-23.03 = Miocene 
# 23.03-33.9 = Oligocene
# 33.9-56.0 = Eocene
# 56.0-66.0 = Paleocene
# 66.0-100.5 = Late Cretaceous
# 100.5-145.0 = Early Cretaceous
# 145.0-201.3 = Jurassic

# Create custom epochs
epochs_custom_df <- data.frame(
  name = c("PPHo", "Miocene", "Oligo.", "Eocene", "Paleo.", "Late Cretaceous", "Early Cret.", "Jurassic"),
  true_name = c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cretaceous", "Early Cretaceous", "Jurassic"),
  max_age = c(5.3, 23.0, 33.9, 55.8, 66.0, 100.5, 145.0, 201.3),
  min_age = c(0, 5.3, 23.0, 33.9, 55.8, 66.0, 100.5, 145.0),
  color = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e", "#4eb2ca"),
  grey_color = c("white", "grey90", "white", "grey90", "white", "grey90", "white", "grey90"))

# Extract root age
Youngest_root_age <- max(phytools::nodeHeights(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo))
Oldest_root_age <- max(phytools::nodeHeights(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo))

# Convert age in time to root
Youngest_epochs_custom_df <- Oldest_epochs_custom_df <- epochs_custom_df
Youngest_epochs_custom_df$time_since_root_min <- Youngest_root_age - epochs_custom_df$min_age
Youngest_epochs_custom_df$time_since_root_max <- Youngest_root_age - epochs_custom_df$max_age
Oldest_epochs_custom_df$time_since_root_min <- Oldest_root_age - epochs_custom_df$min_age
Oldest_epochs_custom_df$time_since_root_max <- Oldest_root_age - epochs_custom_df$max_age

### 3.6.2/ Prepare Genus-group metadata for 1534t phylogenies

# 7 genus-groups from Schmidt & Shattuck, 2014
  # Platythyrea
  # Pachycondyla: Pachycondyla + Simopelta + Thaumatomyrmex + Belonopelta + Mayaponera + Raspone + Dinoponera + Neoponera
  # Ponera: Ponera + Diacamma + Emeryopone + Austroponera + Pseudoponera + Parvaponera + Wadeura + Cryptoponera + Iroponera + Ectomomyrmex
  # Harpegnathos
  # Hypoponera 
  # Plectroctena: Plectoponera + Centromyrmex + Psalidomyrmex + Loboponera + Boloponera
  # Odontomachus: A lot of stuff... (ca. 20 genera)

# List Genus-groups from Schmidt & Shattuck, 2014
genus_groups_list <- c("Platythyrea", "Pachycondyla", "Ponera", "Harpegnathos", "Hypoponera", "Plectroctena", "Odontomachus")


## For Youngest_phylogeny

# List associated taxa to use to retrieve MRCA from the grafted 1534t tree
Youngest_genus_groups_MRCA_taxa_list_1534t <- list(Platythyrea = c("Platythyrea arnoldi", "Platythyrea schultzei"),
                                                Pachycondyla = c("Simopelta minima", "Belonopelta attenuata"),
                                                Ponera = c("Diacamma magdalenae", "Austroponera castanea"),
                                                # Harpegnathos = c("Harpegnathos saltator", "Harpegnathos my01"), # For UCE tree
                                                Harpegnathos = c("Harpegnathos honestoi", "Harpegnathos my01"), # For MCC grafted tree
                                                Hypoponera = c("Hypoponera inexorata", "Hypoponera comis"),
                                                Plectroctena =  c("Centromyrmex fugator", "Boloponera ikemkha"),
                                                Odontomachus = c("Neoponera bucki", "Anochetus talpa"))

# Initiate metadata df
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID <- NA
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_edge_ID <- NA
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$crown_age <- NA
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$stem_age <- NA

all_edge_ages <- phytools::nodeHeights(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo)

root_age <- max(all_edge_ages[, 2])
all_edge_ages <- round(-1 * all_edge_ages + root_age, 5)

for (i in 1:nrow(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- Youngest_genus_groups_MRCA_taxa_list_1534t[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_Youngest_phylogeny_1534t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant taxa
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_Youngest_phylogeny_1534t_treedata@phylo, node = MRCA_node_ID_i)
  all_descendants_ID_i <- all_descendants_ID_i[all_descendants_ID_i <= length(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label)]
  all_descendants_names_i <- Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label[all_descendants_ID_i]
  
  # Record current richness data
  Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$current_richness[i] <- length(all_descendants_names_i)
 
  # Extract crown and stem ages
  MRCA_edge_ID_i <- which(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_age_i <- all_edge_ages[MRCA_edge_ID_i, 2]
  stem_age_i <- all_edge_ages[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and ages
  Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$crown_age[i] <- round(crown_age_i, 1)
  Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$stem_age[i] <- round(stem_age_i, 1)
  
}

View(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df)

# Save Genus-groups metadata
saveRDS(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df.rds")


## For Oldest_phylogeny

# List associated taxa to use to retrieve MRCA from the grafted 1534t tree
Oldest_genus_groups_MRCA_taxa_list_1534t <- list(Platythyrea = c("Platythyrea arnoldi", "Platythyrea schultzei"),
                                               Pachycondyla = c("Simopelta minima", "Belonopelta attenuata"),
                                               Ponera = c("Diacamma magdalenae", "Austroponera castanea"),
                                               # Harpegnathos = c("Harpegnathos saltator", "Harpegnathos my01"), # For UCE tree
                                               Harpegnathos = c("Harpegnathos venator", "Harpegnathos my01"), # For MCC grafted tree
                                               Hypoponera = c("Hypoponera inexorata", "Hypoponera comis"),
                                               Plectroctena =  c("Centromyrmex fugator", "Boloponera ikemkha"),
                                               Odontomachus = c("Neoponera bucki", "Anochetus talpa"))

# Initiate metadata df
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID <- NA
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_edge_ID <- NA
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$crown_age <- NA
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$stem_age <- NA

all_edge_ages <- phytools::nodeHeights(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo)

root_age <- max(all_edge_ages[, 2])
all_edge_ages <- round(-1 * all_edge_ages + root_age, 5)

for (i in 1:nrow(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- Oldest_genus_groups_MRCA_taxa_list_1534t[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_Oldest_phylogeny_1534t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant taxa
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_Oldest_phylogeny_1534t_treedata@phylo, node = MRCA_node_ID_i)
  all_descendants_ID_i <- all_descendants_ID_i[all_descendants_ID_i <= length(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label)]
  all_descendants_names_i <- Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label[all_descendants_ID_i]
  
  # Record current richness data
  Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$current_richness[i] <- length(all_descendants_names_i)
  
  # Extract crown and stem ages
  MRCA_edge_ID_i <- which(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_age_i <- all_edge_ages[MRCA_edge_ID_i, 2]
  stem_age_i <- all_edge_ages[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and ages
  Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$crown_age[i] <- round(crown_age_i, 1)
  Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$stem_age[i] <- round(stem_age_i, 1)
  
}

View(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df)

# Save Genus-groups metadata
saveRDS(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df.rds")


### 3.6.2/ Rectangular plot for Youngest_phylogeny with node age range ####

## Load Genus-group metadata for 1534t tree
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df.rds")

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name
Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$color <- genus_groups_colors

Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df
Ponerinae_Youngest_phylogeny_1534t_treedata@data

## Plot tree as background
Ponerinae_phylogeny_plot <- ggtree(tr = Ponerinae_Youngest_phylogeny_1534t_treedata,
                                   layout = "rectangular",
                                   ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                                   linewidth = NA)

## Add Genus-group cladelabels/colored rectangles

# # Adjust colored rectangles
# extend_rect_to <- 0.25

# Loop per genus-group to add labels/rectangles
for (i in 1:nrow(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # #  Add Background color on the whole clade
    # geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
    #              fill = Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$color[i],
    #              # type = "gradient", # gradient.direction = 'rt',
    #              # align = "right",
    #              extendto = extend_rect_to,
    #              alpha = .5) # +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],
                    label = paste0('bold(',Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 20,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$color[i], Ponerinae_Youngest_phylogeny_1534t_all_genus_groups_metadata_df$color[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    # offset = 8.0,
                    offset = 2.0,
                    offset.text = c(10.0, 10.0, 10.0, 16.0, 10.0, 10.0, 10.0)[i],
                    barsize = 60)    # Width of the bar
  # barsize = 20)    # Width of the bar
}

## Loop per geological epochs to add polygons
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_rect(data = Youngest_epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
              xmin = time_since_root_min, xmax = time_since_root_max, 
              # fill = color),
              fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

## Plot tree, tip labels, node age ranges
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_tree(colour = "black",
            # mapping = aes(colour = new_parental_edge), # Color as branch status (UCE or grafted)
            # mapping = aes(colour = missing_node), # For ARE treedata
            # color = branches_color,
            linewidth = 1.2) +
  
  # Add tip labels
  geom_tiplab(# mapping = aes(colour = missing_node),
    # color = tips_color, # Color as tip status (UCE or grafted) 
    color = "black",
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.15,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # # Adjust tip/edge color scheme
  # scale_colour_manual("Taxa source", breaks = c(F, T), values = UCE_grafted_colors, labels = c("UCE data", "Grafted")) +
  # guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Add age range intervals on backbone nodes
  geom_range(range = "height_range_backbone", center = "height",
             color = "steelblue", size = 2, alpha = 0.5) +
  
  ### Set title
  labs(title = paste0("Youngest ages - Grafted tree with age ranges - 1534t")) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.01,0.05,0,0), "npc"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        legend.title = element_text(size = 48, vjust = 0, hjust = 0, face = "bold"),
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(5.0,"lines"),
        legend.position.inside = c(0.10, 0.92)
  )

## Add age lines
age_step <- 20
root_age <- max(phytools::nodeHeights(Ponerinae_Youngest_phylogeny_1534t_treedata@phylo))
age <- root_age - age_step
while (age > -100)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, 
               # linetype = 92,
               linetype = "dotted",
               linewidth = 1, color = 'grey20')
  age  <- age - age_step
}

## Add the geoscale
Youngest_epochs_custom_df$min_age <- Youngest_epochs_custom_df$time_since_root_min
Youngest_epochs_custom_df$max_age <- Youngest_epochs_custom_df$time_since_root_max

Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  # scale_x_reverse("Age (Mya)") +
  
  deeptime::coord_geo(dat = list(Youngest_epochs_custom_df), lwd = 0.3,
                      size = 12, fontface = "bold",
                      height = unit(5.0, "line"),
                      expand = T,
                      # xlim = c(root_age, 0), 
                      # ylim = c(0, 1534),
                      clip = "off")


# Export PDF
pdf(file = "./outputs/Final_phylogenies/Ponerinae_Youngest_phylogeny_1534t_rect_with_age_ranges.pdf", height = 120, width = 40)

print(Ponerinae_phylogeny_plot)

dev.off()


### 3.6.3/ Rectangular plot for Oldest_phylogeny with node age range ####

## Load Genus-group metadata for 1534t tree
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df.rds")

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df$group_name
Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$color <- genus_groups_colors

Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df
Ponerinae_Oldest_phylogeny_1534t_treedata@data

## Plot tree as background
Ponerinae_phylogeny_plot <- ggtree(tr = Ponerinae_Oldest_phylogeny_1534t_treedata,
                                   layout = "rectangular",
                                   ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                                   linewidth = NA)

## Add Genus-group cladelabels/colored rectangles

# # Adjust colored rectangles
# extend_rect_to <- 0.25

# Loop per genus-group to add labels/rectangles
for (i in 1:nrow(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    
    # #  Add Background color on the whole clade
    # geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
    #              fill = Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$color[i],
    #              # type = "gradient", # gradient.direction = 'rt',
    #              # align = "right",
    #              extendto = extend_rect_to,
    #              alpha = .5) # +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],
                    label = paste0('bold(',Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 20,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$color[i], Ponerinae_Oldest_phylogeny_1534t_all_genus_groups_metadata_df$color[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    # offset = 8.0,
                    offset = 2.0,
                    offset.text = c(10.0, 10.0, 10.0, 16.0, 10.0, 10.0, 10.0)[i],
                    barsize = 60)    # Width of the bar
  # barsize = 20)    # Width of the bar
}

## Loop per geological epochs to add polygons
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_rect(data = Oldest_epochs_custom_df, inherit.aes = FALSE,
            aes(# xmin = min_age, xmax = max_age, 
              xmin = time_since_root_min, xmax = time_since_root_max, 
              # fill = color),
              fill = grey_color),
            ymin = 0, ymax = length(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo$tip.label),
            show.legend = F) +
  
  scale_fill_manual("Geological epochs",
                    values = c("white", "grey90"),
                    breaks = c("white", "grey90"))

## Plot tree, tip labels, node age ranges
Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  geom_tree(colour = "black",
            # mapping = aes(colour = new_parental_edge), # Color as branch status (UCE or grafted)
            # mapping = aes(colour = missing_node), # For ARE treedata
            # color = branches_color,
            linewidth = 1.2) +
  
  # Add tip labels
  geom_tiplab(# mapping = aes(colour = missing_node),
    # color = tips_color, # Color as tip status (UCE or grafted) 
    color = "black",
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.15,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # # Adjust tip/edge color scheme
  # scale_colour_manual("Taxa source", breaks = c(F, T), values = UCE_grafted_colors, labels = c("UCE data", "Grafted")) +
  # guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Add age range intervals on backbone nodes
  geom_range(range = "height_range_backbone", center = "height",
             color = "steelblue", size = 2, alpha = 0.5) +
  
  ### Set title
  labs(title = paste0("Oldest ages - Grafted tree with age ranges - 1534t"))
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.01,0.05,0,0), "npc"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        legend.title = element_text(size = 48, vjust = 0, hjust = 0, face = "bold"),
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(5.0,"lines"),
        legend.position.inside = c(0.10, 0.92)
  )

## Add age lines
age_step <- 20
root_age <- max(phytools::nodeHeights(Ponerinae_Oldest_phylogeny_1534t_treedata@phylo))
age <- root_age - age_step
while (age > -40)
{ 
  Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
    geom_vline(xintercept = age, 
               # linetype = 92,
               linetype = "dotted",
               linewidth = 1, color = 'grey20')
  age  <- age - age_step
}

## Add the geoscale
Oldest_epochs_custom_df$min_age <- Oldest_epochs_custom_df$time_since_root_min
Oldest_epochs_custom_df$max_age <- Oldest_epochs_custom_df$time_since_root_max

Ponerinae_phylogeny_plot <- Ponerinae_phylogeny_plot +
  
  # scale_x_reverse("Age (Mya)") +
  
  deeptime::coord_geo(dat = list(Oldest_epochs_custom_df), lwd = 0.3,
                      size = 12, fontface = "bold",
                      height = unit(5.0, "line"),
                      expand = T,
                      # xlim = c(root_age, 0), 
                      # ylim = c(0, 1534),
                      clip = "off")


# Export PDF
pdf(file = "./outputs/Final_phylogenies/Ponerinae_Oldest_phylogeny_1534t_rect_with_age_ranges.pdf", height = 120, width = 40)

print(Ponerinae_phylogeny_plot)

dev.off()


### To do in Postprod in Illustrator

# Add genus-group names as label
# Add cladelabel for Genera in alternating colors (beware of polyphyletic groups (Use numbers))
# Replace color rectangles with background fading rectangles
# Adjust geotime scale


