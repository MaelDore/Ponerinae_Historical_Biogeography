
##### Script 00: Run full analyses on Platythyrea #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Obtain a BSM GIF to illustrate the method on a subclade

###

### Inputs

# Ponerinae phyogeny
# Bioregion assignments

###

### Outputs

# Pruned phyogeny
# DEC+J model fit
# ARE map
# BSM simulations
# BSM GIF

###


# Clean environment
rm(list = ls())


##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(openxlsx)  # Use Rccp. No need of Java
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(phytools)
library(geiger)
library(ape)
library(BioGeoBEARS)
library(parallel)
library(qpdf)

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Time-calibrated phylogeny(ies) ####

Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")

# Detect Platythyrea species

Genus_all_sp_list <- str_split(string = Ponerinae_phylogeny_1534t_short_names$tip.label, pattern = "_", simplify = T)[,1]
Platythyrea_sp_ID <- which(Genus_all_sp_list == "Platythyrea")
Platythyrea_sp_names <- Ponerinae_phylogeny_1534t_short_names$tip.label[Platythyrea_sp_ID]
length(Platythyrea_sp_names) # 44 Platythyrea species

# Prune phylogeny to keep only Platythyrea species
Platythyrea_phylogeny_44t <- ape::keep.tip(phy = Ponerinae_phylogeny_1534t_short_names, tip = Platythyrea_sp_names)

plot(Platythyrea_phylogeny_44t)

# Save Platythyrea phylogeny
saveRDS(object = Platythyrea_phylogeny_44t, file = "./outputs/Grafting_missing_taxa/Platythyrea_phylogeny_44t.rds")


### 1.3/ Load biogeographic ranges per species ####

Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Keep only Platythyrea data
Platythyrea_bioregions_binary_table <- Taxa_bioregions_binary_table %>% 
  filter(Current_name %in% Platythyrea_sp_names)

# Save biogeographic range table
saveRDS(object = Platythyrea_bioregions_binary_table, file = "./input_data/Biogeographic_data/Platythyrea_bioregions_binary_table.rds")


##### 2/ Format data for BioGeoBEARS #####

### 2.1/ Phylogeny ####

# Load Platythyrea phylogeny
Platythyrea_phylogeny_44t <- readRDS(file = "./outputs/Grafting_missing_taxa/Platythyrea_phylogeny_44t.rds")

# Export under newick format
write.tree(Platythyrea_phylogeny_44t, file = "./outputs/Grafting_missing_taxa/Platythyrea_phylogeny_44t.tree")

# Set path to Newick tree
path_to_tree <- "./outputs/Grafting_missing_taxa/Platythyrea_phylogeny_44t.tree"
path_to_tree <- BioGeoBEARS::np(path_to_tree)


### 2.2/ Biogeographic range table ####

# Load biogeographic range table
Platythyrea_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Platythyrea_bioregions_binary_table.rds")

# Check match in taxa names
table(Platythyrea_bioregions_binary_table$Current_name %in% Platythyrea_phylogeny_44t$tip.label)
table(Platythyrea_phylogeny_44t$tip.label %in% Platythyrea_bioregions_binary_table$Current_name)

# Extract binary df of presence/absence
binary_df <- Platythyrea_bioregions_binary_table[, -1]
# Convert character strings into numerical factors
binary_df_num <- as.data.frame(apply(X = binary_df, MARGIN = 2, FUN = as.numeric))
binary_df_factors <- binary_df_num
binary_df_factors[] <- lapply(binary_df_factors, factor) 
row.names(binary_df_factors) <- Platythyrea_bioregions_binary_table$Current_name
row.names(binary_df_num) <- Platythyrea_bioregions_binary_table$Current_name

dim(binary_df_factors)

# Convert area names in unique letters (needed for Biogeographic stochastic mapping)
names(binary_df_num) <- c("A", "U", "I", "R", "N", "E", "W")

# Produce tipranges object from numeric df
Platythyrea_bioregions_tipranges_obj <- define_tipranges_object(tmpdf = binary_df_num)

# Produce Lagrange PHYLIP biogeographic data
save_tipranges_to_LagrangePHYLIP(tipranges_object = Platythyrea_bioregions_tipranges_obj,
                                 lgdata_fn = "./input_data/BioGeoBEARS_setup/Platythyrea_lagrange_area_data_file_7_regions_PaleA.data",
                                 areanames = colnames(Platythyrea_bioregions_tipranges_obj@df))

# Set path to Lagrange PHYLIP biogeographic data
path_to_geo_data <- "./input_data/BioGeoBEARS_setup/Platythyrea_lagrange_area_data_file_7_regions_PaleA.data"
path_to_geo_data <- BioGeoBEARS::np(path_to_geo_data)


##### 3/ Set hyperparameters for all runs #####

### 3.1/ Set the maximum number of bioregions occupied by a lineage at any time ####

# Extract the range size (nb of bioregion occupied by terminals)
range_size <- rowSums(dfnums_to_numeric(Platythyrea_bioregions_tipranges_obj@df))
current_max_range_size <- max(range_size) # Current max range size

# Extract number of areas (bioregions)
nb_areas <- ncol(Platythyrea_bioregions_tipranges_obj@df)

# Cannot be lower than the observed max range size at the tips (3)
current_max_range_size
# Cannot be higher than the number of bioregions (7)
nb_areas

# max_range_size <- 3 # Smallest possible
max_range_size <- 5 # As in Kawahara et al., 2023
# max_range_size <- 7 # Most extensive

# Estimate the number of states needed to represent all possible combinations of bioregions
cladoRcpp::numstates_from_numareas(numareas = nb_areas, maxareas = max_range_size)

# Advice on number of ranges/states
#	less than 1000 if you want the analysis to run in under a day
#	less than 1500 if you want the analysis to run in under a week
#	less than 2500 if you want it to run at all


### 3.2/ Set time boundaries of strata ####

# Use seven time slices associated to geological periods as in Kawahara et al., 2023: 
# 0-5.33 My = Holocene + Pleistocene + Pliocene
# 5.33-23.03 = Miocene 
# 23.03-33.9 = Oligocene
# 33.9-56.0 = Eocene
# 56.0-66.0 = Paleocene
# 66.0-100.5 = Late Cretaceous (Root of Platythyrea = 96.7 My ?)

# Oldest boundary should be any value older than the root
time_boundaries <- c(5.33, 23.03, 33.9, 56.0, 66.0, 100.5)

# Write output file for time_boundaries
writeLines(text = as.character(time_boundaries), con = "./input_data/BioGeoBEARS_setup/Platythyrea_time_boundaries.txt")

# Set path to boundaries for the time-periods
path_to_time_strata_boundaries <- "./input_data/BioGeoBEARS_setup/Platythyrea_time_boundaries.txt"
path_to_time_strata_boundaries <- BioGeoBEARS::np(path_to_time_strata_boundaries)

### 3.3/ Set areas (monomorphic ranges) allowed per time-strata ####

# Provide list of areas (monomorphic ranges) allowed by time-frame as binary matrices (should be just lists, but are matrices due to historical contingency in programming)
writeLines(readLines("./input_data/BioGeoBEARS_setup/Platythyrea_areas_allowed.txt"))
# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by removing all ranges/states that include the non-allowed areas
# Better option is to use the manual list of allowed states as $lists_of_states_lists_0based

# Optional if using adjacency matrices and $lists_of_states_lists_0based

# Set path to areas allowed
path_to_areas_allowed <- "./input_data/BioGeoBEARS_setup/Platythyrea_areas_allowed.txt"
path_to_areas_allowed <- BioGeoBEARS::np(path_to_areas_allowed)


### 3.4/ Set time-stratified adjacency matrices ####

# Provide list of bimorphic states allowed by time-frame as binary adjacency matrices
writeLines(readLines("./input_data/BioGeoBEARS_setup/Platythyrea_Adjacency_matrix_7areas_6TS.txt"))
# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by removing all ranges/states that include combination of non-adjacent areas
# Better option is to use the manual list of allowed states as $lists_of_states_lists_0based

# Set path to time-stratified adjacency matrices
path_to_adjacency_matrices <- "./input_data/BioGeoBEARS_setup/Platythyrea_Adjacency_matrix_7areas_6TS.txt"
path_to_adjacency_matrices <- BioGeoBEARS::np(path_to_adjacency_matrices)

# Import as list of df the time-stratified adjacency matrices
list_of_adjacency_matrices <- BioGeoBEARS:::read_areas_adjacency_fn(areas_adjacency_fn = path_to_adjacency_matrices)


## 3.5/ Check/Adjust lists of allowed states per time-strata using manual modification ####

## 3.5.1/ Get the list of monomorphic states = bioregion = areas from the tip data ####
areas_list <- getareas_from_tipranges_object(Platythyrea_bioregions_tipranges_obj)
areas_list

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

ranges_list <- generate_list_ranges(areas_list, max_range_size, include_null_range = TRUE)

# How many states/ranges, by default: 120 (119 + null state)
length(ranges_list)

# Get the 0based version of all ranges
ranges_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(areas = areas_list, maxareas = max_range_size, include_null_range = TRUE)


## 3.5.2/ Remove automatically the non-valid ranges using information from the adjacency matrices ####

# Function to convert binary df (adjacency matrix) into list of bimorphic states
convert_binary_df_to_2states <- function (binary_df, areas_list = NULL, include_null_range = TRUE, non_adjacent = F)
{
  # Initiate vectors for 0based bimorphic states
  bimorphic_states_list <- c()
  
  # Add row.names 
  row.names(binary_df) <- names(binary_df)
  
  # Reorder areas to match the biogeographic table
  if (is.null(areas_list))
  {
    areas_list <- names(binary_df)
  } else {
    binary_df <- binary_df[areas_list, areas_list]
  }
  
  # Extract bimorphic states
  for (i in (1:(ncol(binary_df)-1)))
  {
    for (j in ((i+1):ncol(binary_df)))
    {
      binary_test <- binary_df[i,j]
      
      # If asking for non-adjacent, use the opposite of the binary test
      if (non_adjacent) { binary_test <- !binary_test}
      
      if (binary_test)
      {
        new_state <- paste0(areas_list[i], areas_list[j])
        bimorphic_states_list <- c(bimorphic_states_list, new_state)
      }
    }
  }
  
  return(bimorphic_states_list)
}

# Function to detect all allowed adjacent ranges (0based states) from an adjacent matrix
detect_valid_ranges_from_adjacent_df <- function (adjacent_df, areas_list = NULL, ranges_list = NULL, max_range_size = NULL, include_null_range = TRUE)
{
  # If no areas_list is provided, extract all areas from the adjacent_df
  if (is.null(areas_list))
  {
    areas_list <- names(adjacent_df)
  } else {
    # If areas_list is provided, reorder adjacent_df so it fits the order in areas_list
    row.names(adjacent_df) <- names(adjacent_df)
    adjacent_df <- adjacent_df[areas_list,areas_list]
  }
  
  # If no ranges_list is provided, generate all possible ranges from the areas in the adjacent df
  if (is.null(ranges_list))
  {
    ranges_list <- generate_list_ranges(areas_list = areas_list, max_range_size = max_range_size, include_null_range = include_null_range)
  }
  
  # Initiate vector of valid ranges
  valid_ranges <- c()
  
  # Include null range
  if (include_null_range)
  {
    valid_ranges <- c("_")
  }
  
  # Include monomorphic ranges = areas
  valid_ranges <- c(valid_ranges, areas_list)
  
  ## Generate list of bimorphic states (ranges2) that are allowed in the adjacency matrix
  # valid_ranges2_list <- ranges_list[nchar(ranges_list)] # Do not use the full list of ranges2 as some are not allowed in the adjacency matrix
  valid_ranges2_list <- convert_binary_df_to_2states(binary_df = adjacent_df, areas_list = areas_list, include_null_range = include_null_range)
  
  # Add bimorphic states in the list of valid ranges
  valid_ranges <- c(valid_ranges, valid_ranges2_list)
  
  ## Extract higher ranges (size > 2)
  higher_ranges_list <- ranges_list[nchar(ranges_list) > 2]
  
  if (length(higher_ranges_list) > 0)
  {
    ## Check validity of higher range
    # Higher range is valid if it exists a combination of valid bimorphic states that encompass all areas in the higher range
    # As such, it means it exists a path to connect all areas in the higher range using valid adjacency transitions
    
    for (i in seq_along(higher_ranges_list))
    {
      # i <- 90
      focal_range <- higher_ranges_list[i]
      
      # Extract all areas in the higher range
      focal_areas <- unlist(str_split(string = focal_range, pattern = ""))
      # Reorder so they fits the order in areas_list
      focal_areas <- focal_areas[order(match(x = focal_areas, table = areas_list))]
      
      # Generate all nested bimorphic states (ranges2)
      all_nested_ranges2 <- unlist(apply(X = combn(x = focal_areas, m = 2), MARGIN = 2, FUN = paste, collapse = ""))
      # Extract valid nested bimorphic states (ranges2)
      valid_nested_ranges2 <- all_nested_ranges2[all_nested_ranges2 %in% valid_ranges2_list]
      # Identify areas encompasses by the set of valid nested bimorphic states (ranges2)
      valid_areas_in_nested_ranges2 <- unique(unlist(str_split(string = paste(valid_nested_ranges2, collapse = ""), pattern = "")))
      # Check if valid areas encompass all areas in the higher range
      test_validity <- all(focal_areas %in% valid_areas_in_nested_ranges2)
      
      # Add focal higher range only if valid
      if (test_validity) { valid_ranges <- c(valid_ranges, focal_range) }
    }
  }
  
  # Order unique valid ranges as found in ranges_list
  valid_ranges <- valid_ranges[order(match(x = valid_ranges, table = ranges_list))]
  
  # Return valid ranges
  return(valid_ranges)
}

## Initiate final list of lists of time-stratified valid 0based states
lists_of_states_lists_0based <- list()

## STRATUM 1 (0-5.33 My = Holocene + Pleistocene + Pliocene)

# Extract adjacency matrix for Stratum 1
adjacency_df_1 <- list_of_adjacency_matrices[[1]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_1 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_1,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_1) # Now 36 states
# Make a selection vector to keep only valid ranges  
keep_stratum_1 <- (ranges_list %in% valid_ranges_stratum_1)

# Inform the list of valid states
lists_of_states_lists_0based[[1]] <- ranges_list_0based[keep_stratum_1]

## STRATUM 2 (5.33-23.03 = Miocene)

# Extract adjacency matrix for Stratum 2
adjacency_df_2 <- list_of_adjacency_matrices[[2]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_2 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_2,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_2) # Now 36 states
setdiff(valid_ranges_stratum_2, valid_ranges_stratum_1) # Gain no states
setdiff(valid_ranges_stratum_1, valid_ranges_stratum_2) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_2 <- (ranges_list %in% valid_ranges_stratum_2)

# Inform the list of valid states
lists_of_states_lists_0based[[2]] <- ranges_list_0based[keep_stratum_2]

## STRATUM 3 (23.03-33.9 = Oligocene)

# Extract adjacency matrix for Stratum 3
adjacency_df_3 <- list_of_adjacency_matrices[[3]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_3 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_3,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_3) # Now 27 states
setdiff(valid_ranges_stratum_3, valid_ranges_stratum_2) # Gain no states
setdiff(valid_ranges_stratum_2, valid_ranges_stratum_3) # Lose "UI" states
# Make a selection vector to keep only valid ranges  
keep_stratum_3 <- (ranges_list %in% valid_ranges_stratum_3)

# Inform the list of valid states
lists_of_states_lists_0based[[3]] <- ranges_list_0based[keep_stratum_3]

## STRATUM 4 (33.9-56.0 = Eocene)

# Extract adjacency matrix for Stratum 4
adjacency_df_4 <- list_of_adjacency_matrices[[4]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_4 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_4,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_4) # Now 20 states
setdiff(valid_ranges_stratum_4, valid_ranges_stratum_3) # Gain no states
setdiff(valid_ranges_stratum_3, valid_ranges_stratum_4) # Lose "AW" states
# Make a selection vector to keep only valid ranges  
keep_stratum_4 <- (ranges_list %in% valid_ranges_stratum_4)

# Inform the list of valid states
lists_of_states_lists_0based[[4]] <- ranges_list_0based[keep_stratum_4]

## STRATUM 5 (56.0-66.0 = Paleocene)

# Extract adjacency matrix for Stratum 5
adjacency_df_5 <- list_of_adjacency_matrices[[5]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_5 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_5,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_5) # Now 24 states
setdiff(valid_ranges_stratum_5, valid_ranges_stratum_4) # Gain "UN" and "RW" states
setdiff(valid_ranges_stratum_4, valid_ranges_stratum_5) # Lose "RN" states
# Make a selection vector to keep only valid ranges  
keep_stratum_5 <- (ranges_list %in% valid_ranges_stratum_5)

# Inform the list of valid states
lists_of_states_lists_0based[[5]] <- ranges_list_0based[keep_stratum_5]

## STRATUM 6 (66.0-100.5 = Late Cretaceous)

# Extract adjacency matrix for Stratum 6
adjacency_df_6 <- list_of_adjacency_matrices[[6]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_6 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_6,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_6) # Now 20 states
setdiff(valid_ranges_stratum_6, valid_ranges_stratum_5) # Gain no states
setdiff(valid_ranges_stratum_5, valid_ranges_stratum_6) # Lose "EW" states
# Make a selection vector to keep only valid ranges  
keep_stratum_6 <- (ranges_list %in% valid_ranges_stratum_6)

# Inform the list of valid states
lists_of_states_lists_0based[[6]] <- ranges_list_0based[keep_stratum_6]

# Final list of valid states/ranges per time-strata
str(lists_of_states_lists_0based)


## 3.5.3/ Ensure all current tip range are valid ranges in STRATUM 1 ####

## Extract all current tip ranges
current_tip_ranges <- c()
for (i in 1:nrow(Platythyrea_bioregions_tipranges_obj@df))
{
  # i <- 36
  
  binary_range_tip_i <- as.logical(Platythyrea_bioregions_tipranges_obj@df[i, ])
  current_range_tip_i <- names(Platythyrea_bioregions_tipranges_obj@df)[which(binary_range_tip_i)]
  current_range_tip_i <- paste(current_range_tip_i, collapse = "")
  
  # Store tip ranges
  current_tip_ranges <- c(current_tip_ranges, current_range_tip_i)
}
# Name tips
names(current_tip_ranges) <- row.names(Platythyrea_bioregions_tipranges_obj@df)
# Extract unique tip ranges
unique_tip_ranges <- unique(current_tip_ranges)
# Order current_tip_ranges as found in ranges_list
unique_tip_ranges <- unique_tip_ranges[order(match(x = unique_tip_ranges, table = ranges_list))]

## Check that all current tip states are found in the list of allowed states for STRATUM 1
abnormal_current_ranges <- unique_tip_ranges[!(unique_tip_ranges %in% valid_ranges_stratum_1)]
abnormal_current_ranges


##### 4/ Set BioGeoBEARS DEC+J model #####

## 4.1/ Create run object for Time-stratified DEC model ####

# Detect available number of cores
# nb_cores <- parallel::detectCores()
# nb_cores <- round(nb_cores*0.75)

# No need for parallelization for such a small model
nb_cores <- 1

Platythyrea_DEC_run <- define_BioGeoBEARS_run(
  num_cores_to_use = nb_cores, # To set parallel processing
  cluster_already_open = FALSE, # To ask BioGeoBEARS to open (and close) the cluster itself
  max_range_size = max_range_size, # To set the maximum number of bioregion encompassed by a lineage range at any time
  states_list = NULL, # To provide list of allowed states manually when using a single time stratum. Very useful to avoid overparametrization when authorizing a few carefully selected 'polymorphic' ranges
  trfn = path_to_tree, # To provide path to the input tree file
  geogfn = path_to_geo_data, # To provide path to the LagrangePHYLIP file with binary ranges
  timesfn = path_to_time_strata_boundaries, # To provide path to optional file defining the stratified time frames (typically, geographic epochs)
  distsfn = NA, # To provide path to optional file defining the changes in distances between bioregions across time-strates (+X models)
  dispersal_multipliers_fn = NA, # To provide path to optional file defining changes in dispersal multipliers across time-strates (+W models)
  area_of_areas_fn = NA, # To provide path to optional file defining the area of each bioregion (used to weight local extinction rates)
  # areas_allowed_fn = path_to_areas_allowed, # To provide path to file defining the allowed areas at a given time-frame.
  areas_allowed_fn = NA, # No need to provide them if using $lists_of_states_lists_0based
  # areas_adjacency_fn = path_to_adjacency_matrices, # To provide path to optional file defining the allowed bimorphic states/connections between bioregions at a given time-frame.
  areas_adjacency_fn = NA, # No need to provide them if using $lists_of_states_lists_0based
  return_condlikes_table = TRUE, # To ask to obtain all marginal likelihoods computed by the model and used to display ancestral states
  # fixnode = NA, # To provide fixed states for specific nodes
  fixlikes = NA) # To provide weighting scheme for tip likelihoods if you want to account for uncertainty in tip ranges, instead of using 1 for the observed state and 0 for the others.

# Add manually the list of allowed states (per time-strata). Very useful to avoid overparametrization when authorizing a few carefully selected 'polymorphic' ranges
Platythyrea_DEC_run$lists_of_states_lists_0based <- lists_of_states_lists_0based

## 4.2/ Create Time-strata ####

# Loads the time-periods from the text files into the model object.
Platythyrea_DEC_run <- readfiles_BioGeoBEARS_run(Platythyrea_DEC_run)

# Divide the tree up by timeperiods/strata
Platythyrea_DEC_run <- BioGeoBEARS::section_the_tree(inputs = Platythyrea_DEC_run,
                                         make_master_table = TRUE,
                                         plot_pieces = FALSE,
                                         fossils_older_than = 0.001,
                                         cut_fossils = FALSE)
# The stratified tree is summarized in this table describing stratum membership per branches
Platythyrea_DEC_run$master_table

# Check information provided as input to BioGeoBEARS for our run
str(Platythyrea_DEC_run, max.level = 2)

### 4.3/ Inspect the parameter table ####

## Properties of DEC model

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# No jump dispersal (j = 0)
# Here, relative weight of y,v,s are fixed to 1.
# Model as only 2 free parameters = d (range extension) and e (range contraction)

## Parameters

# d = dispersal rate = anagenetic range extension. Ex: A -> AB
# e = extinction rate = anagenetic range contraction. Ex: AB -> A
# a = range-switching rate = anagenetic jump dispersal. Ex: A -> B

# x = used in DEC+X to account for effect of geographic distances between bioregions
# n = used to account for effect of environmental distances between bioregions
# w = used in DEC+W to account for effect of manual dispersal multipliers
# u = used to account for bioregion-specific risk of extinction based on areas of bioregions

# j = Jump-dispersal relative weight = cladogenetic founder-event. Ex: A -> (A),(B)
# y = Non-transitional speciation relative weight = cladogenetic inheritence. Ex: A -> (A),(A)
# s = Subset speciation relative weight = cladogenetic sympatric speciation. Ex: AB -> (AB),(A)
# v = Vicariance relative weight = cladogenetic vicariance. Ex: AB -> (A),(B)

Platythyrea_DEC_run$BioGeoBEARS_model_object@params_table

# Check that starting parameter values are inside the min/max
Platythyrea_DEC_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = Platythyrea_DEC_run)
# Check validity of set-up before run
check_BioGeoBEARS_run(Platythyrea_DEC_run)

# Save run settings
saveRDS(object = Platythyrea_DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/Platythyrea_DEC_run.rds")

## 4.4/ Adjust parameters for DEC+J models ####

# Allows (A) -> (A),(B) # Cladogenetic transition as founder-event speciation = jump-speciation

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# Allows jump dispersal (j)
# Here, relative weight of y,v,s are fixed to 1.
# Model as only 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

Platythyrea_DEC_J_run <- Platythyrea_DEC_run

# Update status of jump speciation parameter to be estimated
Platythyrea_DEC_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
# Set initial value of J for optimization to an arbitrary low non-null value
j_start <- 0.0001
Platythyrea_DEC_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
Platythyrea_DEC_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Check validity of set-up before run
Platythyrea_DEC_J_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = Platythyrea_DEC_J_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(Platythyrea_DEC_J_run)

# Save run settings
saveRDS(object = Platythyrea_DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Platythyrea_DEC_J_run.rds")


##### 5/ Run BioGeoBEARS model ####

### 5.1/ Run DEC+J model ####

## Load run settings
Platythyrea_DEC_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Platythyrea_DEC_J_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DEC_J_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Platythyrea_DEC_J_run_log.txt", open = "wt")
sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DEC+J model run\n"))
Start_time <- Sys.time()
Platythyrea_DEC_J_fit <- bears_optim_run(Platythyrea_DEC_J_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DEC+J model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

print(Platythyrea_DEC_J_fit)

# Save model fit
saveRDS(object = Platythyrea_DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Platythyrea_DEC_J_fit.rds")

### 5.2/ Inspect results ####

compute_parameter_CI_from_Hessian <- function (BioGeoBEARS_output, subset_parameter_indices = NULL, alpha = 0.05, plot = F, model_name = NULL)
{
  # Use a Hessian matrix to compute secondary derivative of the likelihood function in each parameter dimension to obtain a quantification of the curvature
  #	SEs are the square root of the opposite of the diagonal of the inverse of the Hessian matrix
  
  # # Check if the Hessian matrix is positive-definite
  # if (!BioGeoBEARS_output$optim_result$kkt2)
  # {
  #   stop("The Hessian matrix is not positive-definite. Is is not inversible and parameter SD cannot be computed.\n")
  # }
  
  # Detect free parameters
  params_table <- BioGeoBEARS_output$outputs@params_table
  all_pars_names <- row.names(params_table)[params_table$type == "free"]
  
  # Set parameters subsetting if needed
  if (is.null(subset_parameter_indices))
  {
    subset_parameter_indices <- 1:length(all_pars_names)
  }
  
  # Extract parameters MLE
  pars_MLE <- as.numeric(BioGeoBEARS_output$optim_result[1:length(all_pars_names)])[subset_parameter_indices]
  
  # Extract the Hessian matrix
  optim_details <- attr(x = BioGeoBEARS_output$optim_result, which = "details")
  Hessian_matrix <- optim_details[[3]][subset_parameter_indices, subset_parameter_indices]
  
  # Extract parameter variance/sd
  pars_var <- diag(-1 * solve(Hessian_matrix))
  pars_sd <- sqrt(pars_var)
  
  # Transform in CI assuming normal distribution and using the alpha value
  CI_min <- qnorm(p = alpha/2, mean = pars_MLE, sd = pars_sd)
  CI_max <- qnorm(p = 1 - (alpha/2), mean = pars_MLE, sd = pars_sd)
  
  # Store results in data.frame
  parameters_df <- data.frame(pars = all_pars_names[subset_parameter_indices], MLE = pars_MLE, sd = pars_sd, CI_min = CI_min, CI_max = CI_max)
  
  # Plot parameter estimates
  if (plot)
  {
    ggplot_CI <- ggplot(data = parameters_df) +
      geom_point(mapping = aes(x = pars, y = MLE), size = 8) +
      geom_errorbar(mapping = aes(x = pars, ymin = CI_min, ymax = CI_max), size = 2, width = 0.2) +
      geom_hline(yintercept = 0, col = "red", linetype = "dotted", size = 2) +
      
      # Titles
      ggtitle(paste0("Parameter estimates for ",model_name," model", "\n",
                     "CI ",(1-alpha)*100,"%")) +
      xlab("Parameters") +
      ylab("MLE") +
      
      # Aesthetics
      theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", size = 0.5),
            # panel.grid.major.x = element_blank(),
            panel.background = element_rect(fill = NA, color = NA),
            plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
            axis.title = element_text(size = 20, color = "black"),
            axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 12)),
            axis.line = element_line(linewidth = 1.5),
            axis.ticks = element_line(linewidth = 1.5),
            axis.ticks.length = unit(8, "pt"),
            axis.text = element_text(size = 18, color = "black"),
            axis.text.x = element_text(margin = margin(t = 5)),
            axis.text.y = element_text(margin = margin(r = 5)))
    
    print(ggplot_CI)
    
  }
  
  # Return parameter table
  return(parameters_df)
}

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

Platythyrea_DEC_J_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

# Inspect uncertainty in parameter estimates
Platythyrea_DEC_J_pars_table <- compute_parameter_CI_from_Hessian(BioGeoBEARS_output = Platythyrea_DEC_J_fit, subset_parameter_indices = c(1, 3), plot = T, model_name = "DEC+J") ; Platythyrea_DEC_J_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
Platythyrea_DEC_J_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = Platythyrea_DEC_J_fit$inputs, include_null_range = Platythyrea_DEC_J_fit$inputs$include_null_range)
Platythyrea_DEC_J_Qmat <- Platythyrea_DEC_J_mat$Qmat
colnames(Platythyrea_DEC_J_Qmat) <- row.names(Platythyrea_DEC_J_Qmat) <- Platythyrea_DEC_J_mat$ranges_list
Platythyrea_DEC_J_Qmat

# Inspect final parameter table
Platythyrea_DEC_J_fit$outputs

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# With jump dispersal (j)
# Here, relative weights of y,v,s are fixed to 1 - j/3.
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

# Inspect marginal likelihoods of ancestral states
Platythyrea_DEC_J_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(Platythyrea_DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
Platythyrea_DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
Platythyrea_DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


##### 6/ Plot ancestral range estimates #####

# Load DEC+J model output
Platythyrea_DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Platythyrea_DEC_J_fit.rds")

pdf(file = "./outputs/Ancestral_range_estimates_maps/Platythyrea_DEC_J_ARE_map.pdf", width = 20, height = 30)
# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = Platythyrea_DEC_J_fit,
                         analysis_titletxt = "DEC+J - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()


##### 7/ Run Biogeographic Stochastic Mapping ####

?BioGeoBEARS::get_inputs_for_stochastic_mapping
?BioGeoBEARS::runBSM

### 7.1/ Extract inputs needed for BSM from model fit object ####

# Load DEC+J model output
Platythyrea_DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Platythyrea_DEC_J_fit.rds")

Platythyrea_DEC_J_BSM_inputs <- get_inputs_for_stochastic_mapping(res = Platythyrea_DEC_J_fit)

# Save BSM inputs
saveRDS(object = Platythyrea_DEC_J_BSM_inputs, file = "./outputs/BSM/Platythyrea_DEC_J_BSM_inputs.rds")

# Load BSM inputs
Platythyrea_DEC_J_BSM_inputs <- readRDS(file = "./outputs/BSM/Platythyrea_DEC_J_BSM_inputs.rds")

# BSM input object as one list per time-stratum
str(Platythyrea_DEC_J_BSM_inputs, max.level = 1)
str(Platythyrea_DEC_J_BSM_inputs[[1]], max.level = 2)

# Inspect inputs
names(Platythyrea_DEC_J_BSM_inputs[[1]])
unlist(Platythyrea_DEC_J_BSM_inputs[[1]]$ranges_list) # All states/ranges list (All included, even the actually non-allowed ones)
Platythyrea_DEC_J_BSM_inputs[[1]]$numstates # Total number of states/ranges (All included, even the actually non-allowed ones)
Platythyrea_DEC_J_BSM_inputs[[1]]$states_allowed_this_timeperiod_TF # T/F to select allowed states for this time-period
unlist(Platythyrea_DEC_J_BSM_inputs[[1]]$ranges_list)[Platythyrea_DEC_J_BSM_inputs[[1]]$states_allowed_this_timeperiod_TF] # List of allowed states/ranges for this time-period
Platythyrea_DEC_J_BSM_inputs[[1]]$state_indices_0based_this_timeperiod # 0Based states/ranges. Only the allowed ones
Platythyrea_DEC_J_BSM_inputs[[1]]$tipranges # Tip geographic data (binary ranges)
Platythyrea_DEC_J_BSM_inputs[[1]]$Qmat # Transition Q matrix with d/e parameters. Only the allowed states ordered as in States/Ranges list
length(Platythyrea_DEC_J_BSM_inputs[[1]]$independent_likelihoods_by_tree_piece_for_timeperiod_i) # Independent likelihood of each branch for each possible/allowed parental node state

# Transition Q matrix with d/e parameters. 
# States are ordered as in States/Ranges list
# They are different across time-strata as the trait spaces are different!
Platythyrea_DEC_J_BSM_inputs[[1]]$Qmat 
Platythyrea_DEC_J_BSM_inputs[[3]]$Qmat
Platythyrea_DEC_J_BSM_inputs[[5]]$Qmat

### 7.2/ Run BSM  ####

nb_BSM_maps <- 50

Platythyrea_DEC_J_BSM_output <- runBSM(res = Platythyrea_DEC_J_fit,  # Model fit object
                           stochastic_mapping_inputs_list = Platythyrea_DEC_J_BSM_inputs,
                           maxnum_maps_to_try = 200, # Maximum number of stochastic maps to simulate
                           nummaps_goal = nb_BSM_maps, # Number of accepted stochastic maps to target
                           maxtries_per_branch = 40000,
                           save_after_every_try = TRUE,
                           savedir = "./outputs/BSM/",
                           seedval = 12345,
                           wait_before_save = 0.01,
                           master_nodenum_toPrint = 0)

## Save BSM outputs
saveRDS(object = Platythyrea_DEC_J_BSM_output, file = "./outputs/BSM/Platythyrea_DEC_J_BSM_output.rds")

### 7.3/ Inspect BSM outputs = cladogenetic and anagenetic event tables ####

## Load BSM outputs
Platythyrea_DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Platythyrea_DEC_J_BSM_output.rds")

# Extract records of events
Platythyrea_DEC_J_clado_events_tables <- Platythyrea_DEC_J_BSM_output$RES_clado_events_tables
Platythyrea_DEC_J_ana_events_tables <- Platythyrea_DEC_J_BSM_output$RES_ana_events_tables

length(Platythyrea_DEC_J_clado_events_tables) # One summary table per stochastic history
length(Platythyrea_DEC_J_ana_events_tables) # One summary table per stochastic history

## Cladogenetic event tables (one table per simulated biogeographic history)

# One row = 1 node/tip
names(Platythyrea_DEC_J_clado_events_tables[[1]])
nrow(Platythyrea_DEC_J_clado_events_tables[[1]]) # Have more events than nodes/tips because some branches are split across time-strata and one "event" (or absence of event) per branch per strata is recorded
Platythyrea_phylogeny_44t$Nnode + length(Platythyrea_phylogeny_44t$tip.label)
View(Platythyrea_DEC_J_clado_events_tables[[1]])

# Node metadata
Platythyrea_DEC_J_clado_events_tables[[1]]$node # Node ID as in ape
Platythyrea_DEC_J_clado_events_tables[[1]]$label # Node label
table(Platythyrea_DEC_J_clado_events_tables[[1]]$node.type) # Tip, root or internal
Platythyrea_DEC_J_clado_events_tables[[1]]$parent_br # ID of the edge leading to the node = gives you the row of tree$edge and tree$edge.length)
Platythyrea_DEC_J_clado_events_tables[[1]]$edge.length # Length of the branch leading to the node
Platythyrea_DEC_J_clado_events_tables[[1]]$ancestor # ID of the node at the bottom of the branch leading to the node
Platythyrea_DEC_J_clado_events_tables[[1]]$daughter_nds # ID of the nodes descending from the two descending lineages of the node
Platythyrea_DEC_J_clado_events_tables[[1]]$tipnames # Comma separated list of all descending terminals
Platythyrea_DEC_J_clado_events_tables[[1]]$node_ht # Time since the root
Platythyrea_DEC_J_clado_events_tables[[1]]$time_bp # Time before present of the node

# States description. States are sampled using the MLE scaled marginal likelihoods
Platythyrea_DEC_J_clado_events_tables[[1]]$sampled_states_AT_nodes # Node state in this particular history / stochastic map
Platythyrea_DEC_J_clado_events_tables[[1]]$sampled_states_AT_brbots # State at the start of the branch leading to the node in this particular history / stochastic map
Platythyrea_DEC_J_clado_events_tables[[1]]$samp_LEFT_dcorner # State at the start of the left descendant edge of the node in this particular history / stochastic map
Platythyrea_DEC_J_clado_events_tables[[1]]$samp_RIGHT_dcorner # State at the start of the right descendant edge of the node in this particular history / stochastic map

# Event description
Platythyrea_DEC_J_clado_events_tables[[1]]$clado_event_type # Type of cladogenetic event
Platythyrea_DEC_J_clado_events_tables[[1]]$clado_event_txt # Text format describing the event using states
Platythyrea_DEC_J_clado_events_tables[[1]]$clado_dispersal_to # In case of cladogenetic dispersal (i.e, founder event, j), record the area of destination

# Stratum metadata
Platythyrea_DEC_J_clado_events_tables[[1]]$stratum # ID of the time-stratum where is the node
Platythyrea_DEC_J_clado_events_tables[[1]]$time_top # Tipward time boundary of the stratum where is the node
Platythyrea_DEC_J_clado_events_tables[[1]]$time_bot # Rootward time boundary of the stratum where is the node
Platythyrea_DEC_J_clado_events_tables[[1]]$reltimept # Time-width of the time-stratum where is the node
Platythyrea_DEC_J_clado_events_tables[[1]]$piecenum # ID of the tree piece where the branch is within the time-stratum
Platythyrea_DEC_J_clado_events_tables[[1]]$piececlass # Type of tree piece where the branch is within the time-stratum. Can be a segment of a single branch ('subbranch') or a 'subtree' if speciation occurs within the time-stratum

# Subtree pieces metadata (part of the tree within the time-stratum where the branch where the anagenetic event happened can be found)
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBnode
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBnode # ID of the descending node of the branch within the subpiece
table(Platythyrea_DEC_J_clado_events_tables[[1]]$SUBnode.type) # Tip, root or internal branch within the subpiece
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBparent_br # ID of the parent branch within the subpiece (NA for roots)
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBedge.length # Length of the branch within the subpiece
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBancestor # ID of the node at the bottom of the branch within the subpiece
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBdaughter_nds # ID of the nodes descending from the two descending lineages of the node descending from the branch, within the subpiece
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBlabel # Label of the descending node within the subpiece
Platythyrea_DEC_J_clado_events_tables[[1]]$SUBtime_bp # Time of the descending node before the end of the time-stratum 

## Anagenetic event tables (one table per simulated biogeographic history)

names(Platythyrea_DEC_J_ana_events_tables[[1]])
nrow(Platythyrea_DEC_J_ana_events_tables[[1]]) # Number of rows depend on the number of anagenetic events
View(Platythyrea_DEC_J_ana_events_tables[[1]])

# Have similar columns for the metadata copied from the cladogenetic table based on the descending node of the branch where the event is happening

# Branch metadata
Platythyrea_DEC_J_ana_events_tables[[1]]$parent_br # ID of the branch where the event is happening
Platythyrea_DEC_J_ana_events_tables[[1]]$edge.length # Length of the branch where the event is happening
Platythyrea_DEC_J_ana_events_tables[[1]]$ancestor # ID of the node at the bottom of the branch where the event is happening
Platythyrea_DEC_J_ana_events_tables[[1]]$tipnames # Comma separated list of all descending terminals

# Descending node metadata
Platythyrea_DEC_J_ana_events_tables[[1]]$node # ID of the descending node as in ape
Platythyrea_DEC_J_ana_events_tables[[1]]$nodenum_at_top_of_branch # Same = ID of the descending node as in ape
Platythyrea_DEC_J_ana_events_tables[[1]]$label # Label of the descending node
table(Platythyrea_DEC_J_ana_events_tables[[1]]$node.type) # Type of the descending node: Tip or internal (No anagenetic event can happen on the root since it has no branch leading to it)
Platythyrea_DEC_J_ana_events_tables[[1]]$daughter_nds # ID of the nodes descending from the two descending lineages of the descending node
Platythyrea_DEC_J_ana_events_tables[[1]]$node_ht # Time since the root
Platythyrea_DEC_J_ana_events_tables[[1]]$time_bp # Time before present of the node

# Have new columns to describe the anagenetic events
Platythyrea_DEC_J_ana_events_tables[[1]]$trynum # Number of trials needed to generate/simulate an evolutionary history compatible with the starting and ending states of the branch
Platythyrea_DEC_J_ana_events_tables[[1]]$brlen # Length of the branch where the event is happening OR the segment within the time-stratum, when stratified
Platythyrea_DEC_J_ana_events_tables[[1]]$current_rangenum_1based # State before event in numerical format
Platythyrea_DEC_J_ana_events_tables[[1]]$new_rangenum_1based # State after event in numerical format
Platythyrea_DEC_J_ana_events_tables[[1]]$current_rangetxt # State before event in text format
Platythyrea_DEC_J_ana_events_tables[[1]]$new_rangetxt # State after event in text format
Platythyrea_DEC_J_ana_events_tables[[1]]$abs_event_time # Time before present of the event
Platythyrea_DEC_J_ana_events_tables[[1]]$event_time # Relative time of the event on the branch (from the tipward node, or the end of the segment within the time-stratum)
Platythyrea_DEC_J_ana_events_tables[[1]]$abs_event_time + Platythyrea_DEC_J_ana_events_tables[[1]]$event_time # Should equal branch length if not time-stratified
Platythyrea_DEC_J_ana_events_tables[[1]]$event_type # Type of event: range extension (d), range-switching (a), or range contraction (e)
Platythyrea_DEC_J_ana_events_tables[[1]]$event_txt # Text format describing the event using states
Platythyrea_DEC_J_ana_events_tables[[1]]$new_area_num_1based # New destination areas in case of d/a event, in numerical format
Platythyrea_DEC_J_ana_events_tables[[1]]$lost_area_num_1based # Lost source areas in case of e event
Platythyrea_DEC_J_ana_events_tables[[1]]$ana_dispersal_from # Source area selected probabilistically as the source using simulate_source_areas_ana_clado() in case the source range encompasses multiple area
Platythyrea_DEC_J_ana_events_tables[[1]]$dispersal_to # Destination area in d/a event, in text format
Platythyrea_DEC_J_ana_events_tables[[1]]$dispersal_to # Lost source area in case of e event, in text format


##### 8/ Plot all BSM histories ####

# Extract number of stochastic maps
nb_maps <- length(Platythyrea_DEC_J_clado_events_tables)

# Extract complete list of areas used across time-strata
all_areas <- lapply(X = Platythyrea_DEC_J_BSM_inputs, FUN = function (x) x$areas)
all_areas <- unique(unlist(all_areas))

# Extract complete list of "character" states used across time-strata
all_states <- unlist(Platythyrea_DEC_J_BSM_inputs[[1]]$ranges_list)

# Extract complete list of 0based states used across time-strata
all_0based_states <- Platythyrea_DEC_J_BSM_inputs[[1]]$state_indices_0based_all_timeperiods

# Extract max number of areas combined in a range/state
max_range_size <- max(unlist(lapply(X = all_0based_states, FUN = length)))

# Is there a null range?
null_range_check <- any(unlist(lapply(X = all_0based_states, FUN = function (x) any(is.na(x)))))

# Extract list of colors for each possible state across all Time-strata
colors_list_for_states <- get_colors_for_states_list_0based(areanames = all_areas,
                                                            states_list_0based = all_0based_states,
                                                            max_range_size = max_range_size,
                                                            plot_null_range = null_range_check)
names(colors_list_for_states) <- all_states

# Save set of colors used in BSM
saveRDS(object = colors_list_for_states, file = "./outputs/BSM/BSM_maps/Platythyrea_colors_list_for_states.rds")

# Load path to additional script
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

### 8.1/ Loop through the maps and plot to a single PDF (typically fail) ####

pdf(file = paste0("./outputs/BSM/BSM_maps/Platythyrea_DEC_J_BSMs_",nb_maps,"_maps.pdf"),
    height = 8, width = 8)

for (i in 1:nb_maps)
  # for (i in 1:100)
{
  # i <- 1
  
  # Extract the master tables of events summary for map n°i
  Platythyrea_DEC_J_clado_events_table_i <- Platythyrea_DEC_J_clado_events_tables[[i]]
  Platythyrea_DEC_J_ana_events_table_i <- Platythyrea_DEC_J_ana_events_tables[[i]]
  
  ## Extend the master table for cladogenetic events to account for branches that are split across different time-strata
  cols_to_get <- names(Platythyrea_DEC_J_clado_events_table_i[,-ncol(Platythyrea_DEC_J_clado_events_table_i)])
  colnums <- match(cols_to_get, names(Platythyrea_DEC_J_ana_events_table_i))
  ana_events_table_cols_to_add <- Platythyrea_DEC_J_ana_events_table_i[,colnums]
  anagenetic_events_txt_below_node <- rep("none", nrow(ana_events_table_cols_to_add))
  ana_events_table_cols_to_add <- cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  rows_to_get_TF <- ana_events_table_cols_to_add$node <= length(Ponerinae_phylogeny_1534t_short_names$tip.label)
  Platythyrea_DEC_J_clado_events_master_table_i <- rbind(ana_events_table_cols_to_add[rows_to_get_TF,], Platythyrea_DEC_J_clado_events_table_i)
  
  # View(Platythyrea_DEC_J_clado_events_table_i)
  # View(Platythyrea_DEC_J_clado_events_master_table_i)
  
  # Convert the BSM into a modified fit object
  Platythyrea_DEC_J_BSM_map_i <- stochastic_map_states_into_res(res = Platythyrea_DEC_J_fit,  # Model fit object
                                                    master_table_cladogenetic_events = Platythyrea_DEC_J_clado_events_master_table_i,
                                                    stratified = TRUE)
  # Set title
  main_i <- paste0("\n\nDEC + J - Stochastic Map #", i, "/", nb_maps)
  
  # Plot stochastic map with node states
  plot_BioGeoBEARS_results(results_object = Platythyrea_DEC_J_BSM_map_i,
                           analysis_titletxt = main_i,
                           addl_params = list("j"), 
                           label.offset = 0.5,
                           plotwhat = "text",
                           tipcex = 0.7, 
                           titlecex = 1.0,
                           statecex = 0.3, # Works only on pies
                           splitcex = 0.3, # Works only on pies
                           cornercoords_loc = scriptdir,
                           root.edge = TRUE, 
                           colors_list_for_states = colors_list_for_states,
                           skiptree = FALSE,
                           show.tip.label = TRUE)
  # Paint on the branch states
  paint_stochastic_map_branches(res = Platythyrea_DEC_J_BSM_map_i,
                                master_table_cladogenetic_events = Platythyrea_DEC_J_clado_events_master_table_i,
                                colors_list_for_states = colors_list_for_states,
                                lwd = 5, lty = par("lty"),
                                root.edge = TRUE, stratified = TRUE)
  # Replot the nodes states
  plot_BioGeoBEARS_results(results_object = Platythyrea_DEC_J_BSM_map_i,
                           analysis_titletxt = main_i,
                           addl_params = list("j"), 
                           label.offset = 0.5,
                           plotwhat = "text",
                           tipcex = 0.7, 
                           titlecex = 3.0,
                           statecex = 0.3, # Works only on pies
                           splitcex = 0.3, # Works only on pies
                           cornercoords_loc = scriptdir,
                           root.edge = TRUE, 
                           colors_list_for_states = colors_list_for_states,
                           skiptree = TRUE, # DO not plot the tree, only the node states
                           show.tip.label = TRUE)
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Stochastic map n°", i, "/", nb_maps,"\n"))
  }
}

dev.off()


### 8.2/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

# fps <- 5
fps <- 10

# Load pdf as image frames
pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/BSM/BSM_maps/Platythyrea_DEC_J_BSMs_",nb_maps,"_maps.pdf"),
                                                 pages = NULL, density = 100)
magick::image_info(pdf_pointer)

image_resize_and_write_gif(image = pdf_pointer,
                           path =  paste0("./outputs/BSM/BSM_maps/Platythyrea_DEC_J_BSMs_",nb_maps,"_maps_",fps,"fps.gif"),
                           delay = 1/fps, # Time between frames in seconds
                           width = 800, height = 800,
                           loop = TRUE,
                           progress = TRUE)


