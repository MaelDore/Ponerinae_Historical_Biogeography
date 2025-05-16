##### Script 04: Fit models in BioGeoBEARS #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Fit different biogeographic models in BioGeoBEARS
# Compare model fits and select the best one
  # DEC, DEC + J, DIVALIKE, DIVALIKE + J
  # + W versions of the above (?)

###

### Inputs

# (Set of) Time-calibrated phylogeny/ies with imputed missing taxa
# Biogeographic ranges per species
# Hyperparameters for Biogeographic models:
  # Geological time boundaries
  # Time-stratified adjacency matrices (affects range space by canceling some combined areas from the range space)
  # Time-stratified dispersal multiplier matrices (for +W models) (modulates dispersal events (d and j))

###

### Outputs

# Biogeographic model fits
# Comparison of data likelihood according to several models (only possible if using similar biogeographic coding scheme)
# Comparison of posterior probabilities (scaled marginal likelihoods) of ancestral range estimates of important nodes (See table S30 in Kawahara et al., 2023)

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

### 1.2/ Load biogeographic ranges per species ####

Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Odontomachus_simillimus", ]
# View(Biogeographic_database_Ponerinae_curated[Biogeographic_database_Ponerinae_curated$Current_name == "Ponera_swezeyi", ])

# Temporary fixes waiting for occurrence correction

View(Taxa_bioregions_binary_table)

range_size <- apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 1, FUN = sum)
max(range_size)

# Wide-range taxa
Taxa_bioregions_binary_table$Current_name[range_size >= 3]
# Overall counts per bioregions
apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 2, FUN = sum)

# Save updated range table
saveRDS(object = Taxa_bioregions_binary_table, file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")


### 1.3/ Time-calibrated phylogeny(ies) ####

# # Roughly calibrated phylogeny
# Ponerinae_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t.rds")

# MCC phylogeny
Ponerinae_MCC_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_MCC_1534t.rds")

# Youngest phylogeny
Ponerinae_Youngest_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Youngest_1534t.rds")

# Oldest phylogeny
Ponerinae_Oldest_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Oldest_1534t.rds")


# # Posterior phylogenies
# Ponerinae_all_posteriors_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.rds")

# Match names

Ponerinae_MCC_phylogeny_1534t_short_names <- Ponerinae_MCC_phylogeny_1534t
tips_short_names <- Ponerinae_MCC_phylogeny_1534t_short_names$tip.label

# Remove extraction and specimen codes
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_EX.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_CASENT.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_MAMI.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BBX\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BEB\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_D\\d{4}.*")

tips_short_names[!(tips_short_names %in% Taxa_bioregions_binary_table$Current_name)]
Taxa_bioregions_binary_table$Current_name[!(Taxa_bioregions_binary_table$Current_name %in% tips_short_names)]

tips_short_names[tips_short_names == "NewGenus_bucki"] <- "Neoponera_bucki"

Ponerinae_MCC_phylogeny_1534t_short_names$tip.label <- tips_short_names

# # Rescale so the root is at 113 My
# root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_short_names)) ; root_age
# Ponerinae_phylogeny_1534t_short_names <- rescale(x = Ponerinae_phylogeny_1534t_short_names, model = "depth", depth = 113)
# root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_short_names)) ; root_age

# Save imputed phylogeny with short names
# saveRDS(object = Ponerinae_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")
saveRDS(object = Ponerinae_MCC_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")


##### 2/ Format data for BioGeoBEARS #####

### 2.1/ Phylogeny ####

# Load imputed phylogeny with short names
# Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")
Ponerinae_MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")
Ponerinae_Youngest_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Youngest_1534t_short_names.rds")
Ponerinae_Oldest_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_Oldest_1534t_short_names.rds")

# Export under newick format
# write.tree(Ponerinae_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.tree")
write.tree(Ponerinae_MCC_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.tree")
write.tree(Ponerinae_Youngest_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_Youngest_phylogeny_1534t_short_names.tree")
write.tree(Ponerinae_Oldest_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_Oldest_phylogeny_1534t_short_names.tree")

# Set path to Newick tree
# path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.tree"
path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.tree"
# path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_Youngest_phylogeny_1534t_short_names.tree"
# path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_Oldest_phylogeny_1534t_short_names.tree"

path_to_tree <- BioGeoBEARS::np(path_to_tree)


### 2.2/ Biogeographic range table ####

# Load updated range table
Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Check match in taxa names
table(Taxa_bioregions_binary_table$Current_name %in% Ponerinae_phylogeny_1534t_short_names$tip.label)
table(Ponerinae_phylogeny_1534t_short_names$tip.label %in% Taxa_bioregions_binary_table$Current_name)

# Extract binary df of presence/absence
binary_df <- Taxa_bioregions_binary_table[, -1]
# Convert character strings into numerical factors
binary_df_num <- as.data.frame(apply(X = binary_df, MARGIN = 2, FUN = as.numeric))
binary_df_factors <- binary_df_num
binary_df_factors[] <- lapply(binary_df_factors, factor) 
row.names(binary_df_factors) <- Taxa_bioregions_binary_table$Current_name
row.names(binary_df_num) <- Taxa_bioregions_binary_table$Current_name

dim(binary_df_factors)

# Convert area names in unique letters (needed for Biogeographic stochastic mapping)
names(binary_df_num) <- c("A", "U", "I", "R", "N", "E", "W")

# Produce tipranges object from numeric df
Taxa_bioregions_tipranges_obj <- define_tipranges_object(tmpdf = binary_df_num)

# Produce Lagrange PHYLIP biogeographic data
save_tipranges_to_LagrangePHYLIP(tipranges_object = Taxa_bioregions_tipranges_obj,
                                 lgdata_fn = "./input_data/BioGeoBEARS_setup/lagrange_area_data_file_7_regions_PaleA.data",
                                 areanames = colnames(Taxa_bioregions_tipranges_obj@df))

# Set path to Lagrange PHYLIP biogeographic data
path_to_geo_data <- "./input_data/BioGeoBEARS_setup/lagrange_area_data_file_7_regions_PaleA.data"
path_to_geo_data <- BioGeoBEARS::np(path_to_geo_data)


##### 3/ Plot biogeographic tip data ####

# Set the colors we’ll use for plotting
colors <- setNames(object = replicate(ncol(binary_df_factors),
                                      setNames(c("white", "gray30"), 0:1),
                                      simplify = FALSE),
                   nm = colnames(binary_df_factors))

# pdf(file = "./outputs/Phylo_1534t_with_ranges.pdf", width = 20, height = 200)
# pdf(file = "./outputs/Phylo_MCC_1534t_with_ranges.pdf", width = 20, height = 200)
# pdf(file = "./outputs/Phylo_Youngest_1534t_with_ranges.pdf", width = 20, height = 200)
pdf(file = "./outputs/Phylo_Oldest_1534t_with_ranges.pdf", width = 20, height = 200)

# Set plotting parameters
par(mar = c(0.1,0.1,0.1,0.1), oma = c(0,0,0,0)) # bltr
# Graph presence/absence using plotTree.datamatrix
range_map <- phytools::plotTree.datamatrix(# tree = Ponerinae_MCC_phylogeny_1534t_short_names,
                                           # tree = Ponerinae_phylogeny_1534t_short_names,
                                           # tree = Ponerinae_Youngest_phylogeny_1534t_short_names,
                                           tree = Ponerinae_Oldest_phylogeny_1534t_short_names,
                                           X = binary_df_factors,
                                           fsize = 0.7, yexp = 1.1,
                                           header = TRUE, xexp = 1.25, colors = colors)

# Get plot info in "last_plot.phylo"
plot_info <- get("last_plot.phylo", envir=.PlotPhyloEnv)

## Add time line

# Extract root age
# root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_short_names))
# root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_short_names))
# root_age <- max(phytools::nodeHeights(Ponerinae_Youngest_phylogeny_1534t_short_names))
root_age <- max(phytools::nodeHeights(Ponerinae_Oldest_phylogeny_1534t_short_names))

# Define ticks
# ticks_labels <- seq(from = 0, to = 100, by = 20)
ticks_labels <- seq(from = 0, to = 120, by = 20)
axis(side = 1, pos = 0, at = (-1 * ticks_labels) + root_age, labels = ticks_labels, cex.axis = 1.5)
legend(x = root_age/2,
       y = 0 - 5, adj = 0,
       bty = "n", legend = "", title = "Time  [My]", title.cex = 1.5)

## Add a legend
legend(x = plot_info$x.lim[2] - 10,
       y = mean(plot_info$y.lim),
       # adj = c(0,0),
       # x = "topleft",
       legend = c("Absence", "Presence"),
       pch = 22, pt.bg = c("white","gray30"), pt.cex =  1.8,
       cex = 1.5, bty = "n")

dev.off()

# Check for weird patterns = potential mistakes
# Check corrections have been applied


##### 4/ Set hyperparameters for all runs #####

### 4.1/ Set the maximum number of bioregions occupied by a lineage at any time ####

# Extract the range size (nb of bioregion occupied by terminals)
range_size <- rowSums(dfnums_to_numeric(Taxa_bioregions_tipranges_obj@df))
current_max_range_size <- max(range_size) # Current max range size

# Extract number of areas (bioregions)
nb_areas <- ncol(Taxa_bioregions_tipranges_obj@df)

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


### 4.2/ Set time boundaries of strata ####

# Use seven time slices associated to geological periods as in Kawahara et al., 2023: 
  # 0-5.33 My = Holocene + Pleistocene + Pliocene
  # 5.33-23.03 = Miocene 
  # 23.03-33.9 = Oligocene
  # 33.9-56.0 = Eocene
  # 56.0-66.0 = Paleocene
  # 66.0-100.5 = Late Cretaceous
  # 100.5-145.0 = Early Cretaceous (Root of Ponerinae = 123.5 My in MCC)
  # 145.0-201.3 = Jurassic (Only for Oldest hypothesis)

# Need to ensure the oldest boundary is older than the oldest root among posteriors
# Thus, use 201.3 My as oldest boundary

# Oldest boundary should be any value older than the root
time_boundaries <- c(5.33, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0)
# time_boundaries_for_Oldest <- c(5.33, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0, 201.3) # Need Jurassic

# Write output file for time_boundaries
writeLines(text = as.character(time_boundaries), con = "./input_data/BioGeoBEARS_setup/time_boundaries.txt")
# writeLines(text = as.character(time_boundaries_for_Oldest), con = "./input_data/BioGeoBEARS_setup/time_boundaries_for_Oldest.txt")

# Set path to boundaries for the time-periods
path_to_time_strata_boundaries <- "./input_data/BioGeoBEARS_setup/time_boundaries.txt"
path_to_time_strata_boundaries <- BioGeoBEARS::np(path_to_time_strata_boundaries)
# path_to_time_strata_boundaries_for_Oldest <- "./input_data/BioGeoBEARS_setup/time_boundaries_for_Oldest.txt" # Include Stratum 8 (Jurassic)
# path_to_time_strata_boundaries_for_Oldest <- BioGeoBEARS::np(path_to_time_strata_boundaries_for_Oldest) # Include Stratum 8 (Jurassic)

### 4.3/ Set areas (monomorphic ranges) allowed per time-strata ####

# Provide list of areas (monomorphic ranges) allowed by time-frame as binary matrices (should be just lists, but are matrices due to historical contingency in programming)
writeLines(readLines("./input_data/BioGeoBEARS_setup/areas_allowed.txt"))
writeLines(readLines("./input_data/BioGeoBEARS_setup/areas_allowed_for_Oldest.txt")) # Include Stratum 8 (Jurassic)
# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by removing all ranges/states that include the non-allowed areas
# Better option is to use the manual list of allowed states as $lists_of_states_lists_0based

# Optional if using adjacency matrices and $lists_of_states_lists_0based

# Set path to areas allowed
path_to_areas_allowed <- "./input_data/BioGeoBEARS_setup/areas_allowed.txt"
path_to_areas_allowed <- BioGeoBEARS::np(path_to_areas_allowed)
# path_to_areas_allowed_for_Oldest <- "./input_data/BioGeoBEARS_setup/areas_allowed.txt" # Include Stratum 8 (Jurassic)
# path_to_areas_allowed_for_Oldest <- BioGeoBEARS::np(path_to_areas_allowed_for_Oldest)

### 4.4/ Set time-stratified adjacency matrices ####

# Provide list of bimorphic states allowed by time-frame as binary adjacency matrices
writeLines(readLines("./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Adjacency_matrix_7areas_7TS.txt"))
writeLines(readLines("./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Adjacency_matrix_7areas_8TS.txt")) # Include Stratum 8 (Jurassic)

# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by removing all ranges/states that include combination of non-adjacent areas
# Better option is to use the manual list of allowed states as $lists_of_states_lists_0based

# Set path to time-stratified adjacency matrices
path_to_adjacency_matrices <- "./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Adjacency_matrix_7areas_7TS.txt"
path_to_adjacency_matrices <- BioGeoBEARS::np(path_to_adjacency_matrices)
# path_to_adjacency_matrices_for_Oldest <- "./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Adjacency_matrix_7areas_8TS.txt" # Include Stratum 8 (Jurassic)
# path_to_adjacency_matrices_for_Oldest <- BioGeoBEARS::np(path_to_adjacency_matrices_for_Oldest) # Include Stratum 8 (Jurassic)

# Import as list of df the time-stratified adjacency matrices
list_of_adjacency_matrices <- BioGeoBEARS:::read_areas_adjacency_fn(areas_adjacency_fn = path_to_adjacency_matrices)
# list_of_adjacency_matrices_for_Oldest <- BioGeoBEARS:::read_areas_adjacency_fn(areas_adjacency_fn = path_to_adjacency_matrices_for_Oldest)


## 4.5/ Check/Adjust lists of allowed states per time-strata using manual modification ####

## 4.5.1/ Get the list of monomorphic states = bioregion = areas from the tip data ####
areas_list <- getareas_from_tipranges_object(Taxa_bioregions_tipranges_obj)
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


## 4.5.2/ Remove automatically the non-valid ranges using information from the adjacency matrices ####

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

## STRATUM 7 (100.5-145.0 = Early Cretaceous)

# Extract adjacency matrix for Stratum 7
adjacency_df_7 <- list_of_adjacency_matrices[[7]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_7 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_7,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_7) # Now 38 states
setdiff(valid_ranges_stratum_7, valid_ranges_stratum_6) # Gain "AU", "AI", "AN", "UI" states
setdiff(valid_ranges_stratum_6, valid_ranges_stratum_7) # Lose "RE" states
# Make a selection vector to keep only valid ranges  
keep_stratum_7 <- (ranges_list %in% valid_ranges_stratum_7)

# Inform the list of valid states
lists_of_states_lists_0based[[7]] <- ranges_list_0based[keep_stratum_7]

## STRATUM 8 (1145.0-201.3 = Jurassic)

## Only for Oldest_phylogenuy

# Extract adjacency matrix for Stratum 8
adjacency_df_8 <- list_of_adjacency_matrices_for_Oldest[[8]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_8 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_8,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_8) # Now 79 states
setdiff(valid_ranges_stratum_8, valid_ranges_stratum_7) # Gain "AR", "IN", "RN", "EW" states
setdiff(valid_ranges_stratum_7, valid_ranges_stratum_8) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_8 <- (ranges_list %in% valid_ranges_stratum_8)

# Inform the list of valid states
lists_of_states_lists_0based_for_Oldest <- lists_of_states_lists_0based
lists_of_states_lists_0based_for_Oldest[[8]] <- ranges_list_0based[keep_stratum_8]

## Final list of valid states/ranges per time-strata
str(lists_of_states_lists_0based, max.level = 1)
str(lists_of_states_lists_0based_for_Oldest, max.level = 1)

## 4.5.3/ Ensure all current tip range are valid ranges in STRATUM 1 ####

## Extract all current tip ranges
current_tip_ranges <- c()
for (i in 1:nrow(Taxa_bioregions_tipranges_obj@df))
{
  # i <- 36
  
  binary_range_tip_i <- as.logical(Taxa_bioregions_tipranges_obj@df[i, ])
  current_range_tip_i <- names(Taxa_bioregions_tipranges_obj@df)[which(binary_range_tip_i)]
  current_range_tip_i <- paste(current_range_tip_i, collapse = "")
  
  # Store tip ranges
  current_tip_ranges <- c(current_tip_ranges, current_range_tip_i)
}
# Name tips
names(current_tip_ranges) <- row.names(Taxa_bioregions_tipranges_obj@df)
# Extract unique tip ranges
unique_tip_ranges <- unique(current_tip_ranges)
# Order current_tip_ranges as found in ranges_list
unique_tip_ranges <- unique_tip_ranges[order(match(x = unique_tip_ranges, table = ranges_list))]

## Check that all current tip states are found in the list of allowed states for STRATUM 1
abnormal_current_ranges <- unique_tip_ranges[!(unique_tip_ranges %in% valid_ranges_stratum_1)]
abnormal_current_ranges

## Check those weird cases

if (length(abnormal_current_ranges) > 0)
{
  current_tip_ranges %in% abnormal_current_ranges
  
  Abnormal_current_range_df <- data.frame(Current_name = Taxa_bioregions_binary_table$Current_name[current_tip_ranges %in% abnormal_current_ranges],
                                          Current_range = current_tip_ranges[current_tip_ranges %in% abnormal_current_ranges],
                                          Evaluation = NA)
  View(Abnormal_current_range_df)
  write.xlsx(x = Abnormal_current_range_df, file = "./maps/Abnormal_current_range_df.xlsx")
  
  ## If valid. Allow those states in the STRATUM 1, even if they do not match the adjacency matrix
}


## 4.5.4/ Remove manually all unwanted N-areas ranges with N > 2 (not possible with adjacency matrices) ####

# Remove manually some non-valid ranges, per time-strata!
# Useful to remove specific N-areas ranges with N > 2 (not possible with adjacency matrices)

## Example with STRATUM 1 (0-5.33 My = Holocene + Pleistocene + Pliocene)

# # Extract adjacency matrix for Stratum 1
# adjacency_df_1 <- list_of_adjacency_matrices[[1]]
# 
# # Identify valid ranges (all sizes) from the adjacent df
# valid_ranges_stratum_1 <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_1,
#                                                                areas_list = areas_list, ranges_list = ranges_list)
# length(valid_ranges_stratum_1) # Now 36 states
# 
# # Add manual exclusion
# manual_invalid_ranges_stratum_1 <- c("ARNEW", "IRNEW")
# 
# # Make a selection vector to keep only valid ranges  
# keep_stratum_1 <- (ranges_list %in% valid_ranges_stratum_1) & !(ranges_list %in% manual_invalid_ranges_stratum_1)
# 
# sum(keep_stratum_1) # Now 34 states
# 
# # Inform the list of valid states
# lists_of_states_lists_0based[[1]] <- ranges_list_0based[keep_stratum_1]


### 4.6/ Set dispersal multiplier matrices for +W models ####

# Provide list of dispersal multiplier matrices used to weight dispersal (j and d parameters) by time-strata
writeLines(readLines("./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Dispersal_multipler_matrices_7areas_7TS.txt"))
# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by scaling initial dispersal parameter using the fixed multiplier powered by parameter W
   # d' = d ×〖fixed_modifier〗^ w
   # j' = j ×〖fixed_modifier〗^ w

# Set path to time-stratified dispersal multiplier matrices
path_to_dispersal_multiplier_matrices <- "./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Dispersal_multipler_matrices_7areas_7TS.txt"
path_to_dispersal_multiplier_matrices <- BioGeoBEARS::np(path_to_dispersal_multiplier_matrices)

# Import as list of df the time-stratified adjacency matrices
list_of_dispersal_multiplier_matrices <- BioGeoBEARS:::read_areas_adjacency_fn(areas_adjacency_fn = path_to_dispersal_multiplier_matrices)
list_of_dispersal_multiplier_matrices

# Max value = 1 so they relates to the base rate/weight
# Min value > 0 to avoid issue with Inf values for w = 0


## If using dispersal multipliers, need to allow all bimorphic and higher states for which dispersal is not zero!
# (Or it will only affect cladogenetic jump-dispersal)

### 4.7/ Adjust lists of allowed states per time-strata according to dispersal multiplier matrices ####

## Adjust adjacency matrices to fit dispersal multipliers

list_of_adjacency_matrices_for_W_models <- list_of_dispersal_multiplier_matrices

# Set all values >= 0.01 to 1 and all values < 0.01 to 0
list_of_adjacency_matrices_for_W_models <- lapply(X = list_of_adjacency_matrices_for_W_models, FUN = function (x) { as.data.frame(1* (x >= 0.01)) })
list_of_adjacency_matrices_for_W_models

## Export manually updated adjacency matrix in appropriate text file
print(list_of_adjacency_matrices_for_W_models)
print(list_of_adjacency_matrices)

# Set path to time-stratified adjacency matrices for +W models
path_to_adjacency_matrices_for_W_models <- "./input_data/BioGeoBEARS_setup/Dore_2024_BioGeoBears_Adjacency_matrix_7areas_7TS_for_W_models.txt"
path_to_adjacency_matrices_for_W_models <- BioGeoBEARS::np(path_to_adjacency_matrices_for_W_models)


## Initiate final list of lists of time-stratified valid 0based states for +W models
lists_of_states_lists_0based_for_W_models <- list()

areas_list <- getareas_from_tipranges_object(Taxa_bioregions_tipranges_obj)
ranges_list <- generate_list_ranges(areas_list, max_range_size, include_null_range = TRUE)
ranges_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(areas = areas_list, maxareas = max_range_size, include_null_range = TRUE)

## STRATUM 1 (0-5.33 My = Holocene + Pleistocene + Pliocene)

# Extract adjacency matrix for Stratum 1
adjacency_df_1_for_W_models <- list_of_adjacency_matrices_for_W_models[[1]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_1_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_1_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_1_for_W_models) # Now 117 states
# Make a selection vector to keep only valid ranges  
keep_stratum_1_for_W_models <- (ranges_list %in% valid_ranges_stratum_1_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[1]] <- ranges_list_0based[keep_stratum_1_for_W_models]

## STRATUM 2 (5.33-23.03 = Miocene)

# Extract adjacency matrix for Stratum 2
adjacency_df_2_for_W_models <- list_of_adjacency_matrices_for_W_models[[2]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_2_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_2_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_2_for_W_models) # Now 117 states
setdiff(valid_ranges_stratum_2_for_W_models, valid_ranges_stratum_1_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_1_for_W_models, valid_ranges_stratum_2_for_W_models) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_2_for_W_models <- (ranges_list %in% valid_ranges_stratum_2_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[2]] <- ranges_list_0based[keep_stratum_2_for_W_models]

## STRATUM 3 (23.03-33.9 = Oligocene)

# Extract adjacency matrix for Stratum 3
adjacency_df_3_for_W_models <- list_of_adjacency_matrices_for_W_models[[3]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_3_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_3_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_3_for_W_models) # Now 117 states
setdiff(valid_ranges_stratum_3_for_W_models, valid_ranges_stratum_2_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_2_for_W_models, valid_ranges_stratum_3_for_W_models) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_3_for_W_models <- (ranges_list %in% valid_ranges_stratum_3_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[3]] <- ranges_list_0based[keep_stratum_3_for_W_models]

## STRATUM 4 (33.9-56.0 = Eocene)

# Extract adjacency matrix for Stratum 4
adjacency_df_4_for_W_models <- list_of_adjacency_matrices_for_W_models[[4]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_4_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_4_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_4_for_W_models) # Now 117 states
setdiff(valid_ranges_stratum_4_for_W_models, valid_ranges_stratum_3_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_3_for_W_models, valid_ranges_stratum_4_for_W_models) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_4_for_W_models <- (ranges_list %in% valid_ranges_stratum_4_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[4]] <- ranges_list_0based[keep_stratum_4_for_W_models]

## STRATUM 5 (56.0-66.0 = Paleocene)

# Extract adjacency matrix for Stratum 5
adjacency_df_5_for_W_models <- list_of_adjacency_matrices_for_W_models[[5]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_5_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_5_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_5_for_W_models) # Now 117 states
setdiff(valid_ranges_stratum_5_for_W_models, valid_ranges_stratum_4_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_4_for_W_models, valid_ranges_stratum_5_for_W_models) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_5_for_W_models <- (ranges_list %in% valid_ranges_stratum_5_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[5]] <- ranges_list_0based[keep_stratum_5_for_W_models]

## STRATUM 6 (66.0-100.5 = Late Cretaceous)

# Extract adjacency matrix for Stratum 6
adjacency_df_6_for_W_models <- list_of_adjacency_matrices_for_W_models[[6]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_6_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_6_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_6_for_W_models) # Now 117 states
setdiff(valid_ranges_stratum_6_for_W_models, valid_ranges_stratum_5_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_5_for_W_models, valid_ranges_stratum_6_for_W_models) # Lose no states
# Make a selection vector to keep only valid ranges  
keep_stratum_6_for_W_models <- (ranges_list %in% valid_ranges_stratum_6_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[6]] <- ranges_list_0based[keep_stratum_6_for_W_models]

## STRATUM 7 (100.5-145.0 = Early Cretaceous)

# Extract adjacency matrix for Stratum 6
adjacency_df_7_for_W_models <- list_of_adjacency_matrices_for_W_models[[7]]

# Identify valid ranges (all sizes) from the adjacent df
valid_ranges_stratum_7_for_W_models <- detect_valid_ranges_from_adjacent_df(adjacent_df = adjacency_df_7_for_W_models,
                                                               areas_list = areas_list, ranges_list = ranges_list)
length(valid_ranges_stratum_7_for_W_models) # Now 113 states
setdiff(valid_ranges_stratum_7_for_W_models, valid_ranges_stratum_6_for_W_models) # Gain no states
setdiff(valid_ranges_stratum_6_for_W_models, valid_ranges_stratum_7_for_W_models) # Lose "NW" states
# Make a selection vector to keep only valid ranges  
keep_stratum_7_for_W_models <- (ranges_list %in% valid_ranges_stratum_7_for_W_models)

# Inform the list of valid states
lists_of_states_lists_0based_for_W_models[[7]] <- ranges_list_0based[keep_stratum_7_for_W_models]

# Final list of valid states/ranges per time-strata
str(lists_of_states_lists_0based_for_W_models)


##### 5/ Set BioGeoBEARS models #####

## 5.1/ Create run object for Time-stratified DEC model ####

# # Detect available number of cores
# nb_cores <- parallel::detectCores()
# nb_cores <- round(nb_cores*0.75)

# Set nb of cores to 1 as parallelization takes a lot of time to setup, thus is less efficient when number of states are not too big
nb_cores <- 1

# # Set a cluster
# ?parallel::makeCluster
# parallel::makeCluster(rep("localhost",nb_cores), type = "SOCK")

# Set default DEC run
DEC_run <- define_BioGeoBEARS_run(
  num_cores_to_use = nb_cores, # To set parallel processing
  cluster_already_open = FALSE, # To ask BioGeoBEARS to open (and close) the cluster itself
  max_range_size = max_range_size, # To set the maximum number of bioregion encompassed by a lineage range at any time
  states_list = NULL, # To provide list of allowed states manually when using a single time stratum. Very useful to avoid overparametrization when authorizing a few carefully selected 'polymorphic' ranges
  trfn = path_to_tree, # To provide path to the input tree file
  geogfn = path_to_geo_data, # To provide path to the LagrangePHYLIP file with binary ranges
  timesfn = path_to_time_strata_boundaries, # To provide path to optional file defining the stratified time frames (typically, geographic epochs)
  # timesfn = path_to_time_strata_boundaries_for_Oldest, # To provide path to optional file defining the stratified time frames (typically, geographic epochs)
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
DEC_run$lists_of_states_lists_0based <- lists_of_states_lists_0based
# DEC_run$lists_of_states_lists_0based <- lists_of_states_lists_0based_for_Oldest

## 5.2/ Create Time-strata ####

# Loads the time-periods from the text files into the model object.
DEC_run <- readfiles_BioGeoBEARS_run(DEC_run)

# Divide the tree up by timeperiods/strata
DEC_run <- BioGeoBEARS::section_the_tree(inputs = DEC_run,
                                 make_master_table = TRUE,
                                 plot_pieces = FALSE,
                                 fossils_older_than = 0.001,
                                 cut_fossils = FALSE)
# The stratified tree is summarized in this table describing stratum membership per branches
DEC_run$master_table

# Check information provided as input to BioGeoBEARS for our run
str(DEC_run, max.level = 2)

### 5.3/ Inspect the parameter table ####

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

DEC_run$BioGeoBEARS_model_object@params_table

# Check that starting parameter values are inside the min/max
DEC_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_run)
# Check validity of set-up before run
check_BioGeoBEARS_run(DEC_run)

# Save run settings
# saveRDS(object = DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_run.rds")
saveRDS(object = DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_run.rds")
# saveRDS(object = DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DEC_run.rds")
# saveRDS(object = DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DEC_run.rds")


## 5.4/ Adjust parameters for DIVALIKE models ####

# DIVALIKE is a likelihood interpretation of parsimony DIVA, and it is "like DIVA" 
# Similar to, but not identical to, parsimony DIVA.

# Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
# No subset sympatric speciation (s = 0)
# No jump dispersal (j = 0)
# Here, relative weight of y and v are fixed to 1.
# Model has 2 free parameters = d (range extension), e (range contraction)

# Use the DEC model as template
DIVA_run <- DEC_run 

# Remove subset sympatric speciation
DIVA_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVA_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVA_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# Adjust the fixed values of ysv, ys, y, and v to account for the absence of s
DIVA_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j" # Instead of 3-j in DEC
DIVA_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2" # Instead of ysv*2/3 in DEC
DIVA_run$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2" # Instead of ysv*1/3 in DEC
DIVA_run$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2" # Instead of ysv*1/3 in DEC

# Allow widespread vicariance; Narrow vicariance and Wide-range vicariance are set equiprobable
DIVA_run$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVA_run$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVA_run$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Check information provided as input to BioGeoBEARS for our run
str(DIVA_run, max.level = 2)

# Check that starting parameter values are inside the min/max
DIVA_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVA_run) # Check that starting parameter values are inside the min/max
# Check validity of set-up before run
check_BioGeoBEARS_run(DIVA_run)

# Save run settings
# saveRDS(object = DIVA_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_run.rds")
saveRDS(object = DIVA_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_run.rds")
# saveRDS(object = DIVA_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DIVA_run.rds")
# saveRDS(object = DIVA_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DIVA_run.rds")

## 5.5/ Adjust parameters for +J models ####

# Allows (A) -> (A),(B) # Cladogenetic transition as founder-event speciation = jump-speciation

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# Allows jump dispersal (j)
# Here, relative weight of y,v,s are fixed to 1.
# Model as only 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

## Update DEC + J run object with new parameter settings

DEC_J_run <- DEC_run

# Update status of jump speciation parameter to be estimated
DEC_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
# Set initial value of J for optimization to an arbitrary low non-null value
j_start <- 0.0001
DEC_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DEC_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Check validity of set-up before run
DEC_J_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_J_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DEC_J_run)

# Save run settings
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_J_run.rds")
saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_J_run.rds")
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DEC_J_run.rds")
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DEC_J_run.rds")


## Update DIVALIKE + J run object with new parameter settings

DIVA_J_run <- DIVA_run

# Update status of jump speciation parameter to be estimated
DIVA_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
# Set initial value of J for optimization to an arbitrary low non-null value
j_start <- 0.0001
DIVA_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DIVA_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Check validity of set-up before run
DIVA_J_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVA_J_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DIVA_J_run)

# Save run settings
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_J_run.rds")
saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_J_run.rds")
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_run.rds")
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_run.rds")


## 5.6/ Adjust dispersal multipliers and allowed states for +W models ####

# Works by scaling initial dispersal parameters using the fixed multiplier powered by parameter W
  # d' = d ×〖fixed_modifier〗^ w
  # j' = j ×〖fixed_modifier〗^ w

## Update DEC + W run object with new parameter settings

DEC_W_run <- DEC_run

# Inform of path to dispersal multiplier matrices
DEC_W_run$dispersal_multipliers_fn <- path_to_dispersal_multiplier_matrices
# Adjust adjacency matrices and list of allowed states so dispersal multipliers acts on both jump-dispersal (j) AND range extension (d)
# DEC_W_run$areas_adjacency_fn <- path_to_adjacency_matrices_for_W_models # Useless if providing $lists_of_states_lists_0based
# Adjust list of allowed states per time strata
DEC_W_run$lists_of_states_lists_0based <- lists_of_states_lists_0based_for_W_models

# Load the dispersal multiplier and modified adjacency matrices from the text files into the model object.
DEC_W_run <- readfiles_BioGeoBEARS_run(DEC_W_run)

# Update status of dispersal multiplier exponent (+W) to be estimated
DEC_W_run$BioGeoBEARS_model_object@params_table["w","type"] <- "free"

# Set initial value of W as 1 such as there is no modulation of the diserpsal multipliers
w_start <- 1

# # Set initial value of W for optimization to an arbitrary low non-null value
# w_start <- 0.001

# Inform initial value of W
DEC_W_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
DEC_W_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization

# Inspect the parameter table
DEC_W_run$BioGeoBEARS_model_object@params_table

#	ML search for DEC+W model can be hard to converge, 
# So it is recommended to try with different starting point,
# and with “speedup" = FALSE to increase number of maxiter in optimx search to 250 instead of 50*nb of free parameters = 150
DEC_W_run$speedup = FALSE

# Check validity of set-up before run
DEC_W_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_W_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DEC_W_run)

# Save run settings
# saveRDS(object = DEC_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_W_run.rds")
saveRDS(object = DEC_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_W_run.rds")


## Update DIVA + W run object with new parameter settings

DIVA_W_run <- DIVA_run

# Inform of path to dispersal multiplier matrices
DIVA_W_run$dispersal_multipliers_fn <- path_to_dispersal_multiplier_matrices
# Adjust adjacency matrices and list of allowed states so dispersal multipliers acts on both jump-dispersal (j) AND range extension (d)
# DIVA_W_run$areas_adjacency_fn <- path_to_adjacency_matrices_for_W_models # Useless if providing $lists_of_states_lists_0based
# Adjust list of allowed states per time strata
DIVA_W_run$lists_of_states_lists_0based <- lists_of_states_lists_0based_for_W_models

# Load the dispersal multiplier and modified adjacency matrices from the text files into the model object.
DIVA_W_run <- readfiles_BioGeoBEARS_run(DIVA_W_run)

# Update status of dispersal multiplier exponent (+W) to be estimated
DIVA_W_run$BioGeoBEARS_model_object@params_table["w","type"] <- "free"

# Set initial value of W as 1 such as there is no modulation of the diserpsal multipliers
w_start <- 1

# # Set initial value of W for optimization to an arbitrary low non-null value
# w_start <- 0.001

# Inform initial value of W
DIVA_W_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
DIVA_W_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization

# Inspect the parameter table
DIVA_W_run$BioGeoBEARS_model_object@params_table

#	ML search for DIVA+W model can be hard to converge, 
# So it is recommended to try with different starting point,
# and with “speedup" = FALSE to increase number of maxiter in optimx search to 250 instead of 50*nb of free parameters = 150
DIVA_W_run$speedup = FALSE

# Check validity of set-up before run
DIVA_W_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVA_W_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DIVA_W_run)

# Save run settings
# saveRDS(object = DIVA_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_W_run.rds")
saveRDS(object = DIVA_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_run.rds")


## 5.7/ Adjust parameters for +J+W models ####

## Use +W models as template to add the +J parameter
# (In practice, use MLE from +J models to set starting parameters of +J+W models)

## Update DEC + J + W run object with new parameter settings

DEC_JW_run <- DEC_W_run

# Update status of jump speciation parameter to be estimated
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
# Set initial value of J for optimization to an arbitrary low non-null value
j_start <- 0.0001
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Check validity of set-up before run
DEC_JW_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_JW_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DEC_JW_run)

# Save run settings
# saveRDS(object = DEC_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_JW_run.rds")
saveRDS(object = DEC_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_JW_run.rds")

## Update DIVALIKE + J run object with new parameter settings

DIVA_JW_run <- DIVA_W_run

# Update status of jump speciation parameter to be estimated
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
# Set initial value of J for optimization to an arbitrary low non-null value
j_start <- 0.0001
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Check validity of set-up before run
DIVA_JW_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVA_JW_run) # Check that starting parameter values are inside the min/max
check_BioGeoBEARS_run(DIVA_JW_run)

# Save run settings
# saveRDS(object = DIVA_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_JW_run.rds")
saveRDS(object = DIVA_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_run.rds")


##### 6/ Run models #####

# Make several runs with different starting parameter values to ensure optimal likelihood has been reached
# Use the stepwise approach by using MLE estimates of previous runs as staring parameter values when increasing complexity (adding parameters)

# Start with DIVALIKE and DEC
# Then DEC + J, DEC + W (and DIVALIKE + J, DIVALIKE + W)
# Then DEC + J + W (and DIVALIKE + J + W) using MLE from +J models 

# Compare LnLk such as bigger model always have better (higher) LnLk
# Use annealing if issues with optimization (slower but less prone to suboptimal pitfalls) BioGeoBEARS_run_object$use_optimx = "GenSA"

### 6.1/ Run DEC model ####

# Load run settings
# DEC_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DEC_run.rds")
DEC_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_run.rds")
# DEC_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DEC_run.rds")
# DEC_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DEC_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
  # Run time for DEC model, no parallel = 8.8 min
  # Run time for DEC model, 18 cores parallel = 36.9 min
# DEC_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DEC_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DEC_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Youngest_phylogeny_1534t/DEC_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Oldest_phylogeny_1534t/DEC_run_log.txt", open = "wt")

sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DEC model run\n"))
Start_time <- Sys.time()
DEC_fit <- bears_optim_run(DEC_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DEC model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

# Save model fit
# saveRDS(object = DEC_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_fit.rds")
saveRDS(object = DEC_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_fit.rds")
# saveRDS(object = DEC_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_fit.rds")
# saveRDS(object = DEC_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_fit.rds")

print(DEC_fit)

## Inspect optimization

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

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


## Inspect results
DEC_fit$optim_result
# Model as only 2 free parameters = d (range extension) and e (range contraction)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas

# Inspect uncertainty in parameter estimates
DEC_pars_table <- compute_parameter_CI_from_Hessian(DEC_fit, plot = T, model_name = "DEC") ; DEC_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DEC_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_fit$inputs, include_null_range = DEC_fit$inputs$include_null_range)
DEC_Qmat <- DEC_mat$Qmat
colnames(DEC_Qmat) <- row.names(DEC_Qmat) <- DEC_mat$ranges_list
DEC_Qmat

# Inspect final parameter table
DEC_fit$outputs

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# No jump dispersal (j = 0)
# Here, relative weights of y,v,s are fixed to 1.
# Model as only 2 free parameters = d (range extension) and e (range contraction)

# Inspect marginal likelihoods of ancestral states
str(DEC_fit)
DEC_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DEC_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).

# length(Ponerinae_phylogeny_1534t_short_names$tip.label) + Ponerinae_phylogeny_1534t_short_names$Nnode
length(Ponerinae_MCC_phylogeny_1534t_short_names$tip.label) + Ponerinae_MCC_phylogeny_1534t_short_names$Nnode
# length(Ponerinae_Youngest_phylogeny_1534t_short_names$tip.label) + Ponerinae_Youngest_phylogeny_1534t_short_names$Nnode
# length(Ponerinae_Oldest_phylogeny_1534t_short_names$tip.label) + Ponerinae_Oldest_phylogeny_1534t_short_names$Nnode


### 6.2/ Run DIVALIKE model ####

# Load run settings
# DIVA_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_run.rds")
DIVA_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_run.rds")
# DIVA_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DIVA_run.rds")
# DIVA_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DIVA_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
  # Run time for DIVALIKE model, no parallel = 8.9 min
  # Run time for DIVALIKE model, 18 cores parallel = 37.5 min
# DIVA_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DIVA_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DIVA_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Youngest_phylogeny_1534t/DIVA_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Oldest_phylogeny_1534t/DIVA_run_log.txt", open = "wt")

sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DIVALIKE model run\n"))
Start_time <- Sys.time()
DIVA_fit <- bears_optim_run(DIVA_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DIVALIKE model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

# # Fit model using ML
# DIVA_fit <- bears_optim_run(DIVA_run)

print(DIVA_fit)

# Save model fit
# saveRDS(object = DIVA_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")
saveRDS(object = DIVA_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_fit.rds")
# saveRDS(object = DIVA_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_fit.rds")
# saveRDS(object = DIVA_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_fit.rds")

## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DIVA_fit$optim_result
# Model as only 2 free parameters = d (range extension) and e (range contraction)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas

# Inspect uncertainty in parameter estimates
DIVA_pars_table <- compute_parameter_CI_from_Hessian(DIVA_fit, plot = T, model_name = "DIVALIKE") ; DIVA_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DIVA_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DIVA_fit$inputs, include_null_range = DIVA_fit$inputs$include_null_range)
DIVA_Qmat <- DIVA_mat$Qmat
colnames(DIVA_Qmat) <- row.names(DIVA_Qmat) <- DIVA_mat$ranges_list
DIVA_Qmat

# Inspect final parameter table
DIVA_fit$outputs

# Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
# No subset sympatric speciation (s = 0)
# No jump dispersal (j = 0)
# Here, relative weights of y and v are fixed to 1.
# Model has 2 free parameters = d (range extension), e (range contraction)

# Inspect marginal likelihoods of ancestral states
DIVA_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DIVA model with the MLE parameters
dim(DIVA_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DIVA_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DIVA model with the MLE parameters estimated globally
DIVA_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.3/ Run DEC+J model ####

## Load run settings
# DEC_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DEC_J_run.rds")
DEC_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_J_run.rds")
# DEC_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DEC_J_run.rds")
# DEC_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DEC_J_run.rds")

## Use MLE of DEC model as starting values
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_fit.rds")
DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_fit.rds")
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_fit.rds")
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_fit.rds")

## Set starting values for optimization based on MLE in the DEC model
d_start <- DEC_fit$outputs@params_table["d","est"]
e_start <- DEC_fit$outputs@params_table["e","est"]

DEC_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_J_run.rds")
saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_J_run.rds")
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DEC_J_run.rds")
# saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DEC_J_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DEC_J_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DEC_J_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DEC_J_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Youngest_phylogeny_1534t/DEC_J_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Oldest_phylogeny_1534t/DEC_J_run_log.txt", open = "wt")

sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DEC+J model run\n"))
Start_time <- Sys.time()
DEC_J_fit <- bears_optim_run(DEC_J_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DEC+J model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

print(DEC_J_fit)

# Save model fit
# saveRDS(object = DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
saveRDS(object = DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
# saveRDS(object = DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_J_fit.rds")
# saveRDS(object = DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_J_fit.rds")

## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DEC_J_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

# Inspect uncertainty in parameter estimates
DEC_J_pars_table <- compute_parameter_CI_from_Hessian(BioGeoBEARS_output = DEC_J_fit, subset_parameter_indices = c(1,3), plot = T, model_name = "DEC+J") ; DEC_J_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DEC_J_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
DEC_J_Qmat <- DEC_J_mat$Qmat
colnames(DEC_J_Qmat) <- row.names(DEC_J_Qmat) <- DEC_J_mat$ranges_list
DEC_J_Qmat

# Inspect final parameter table
DEC_J_fit$outputs

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# With jump dispersal (j)
# Here, relative weights of y,v,s are fixed to 1 - j/3.
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

# Inspect marginal likelihoods of ancestral states
DEC_J_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.4/ Run DIVALIKE+J model ####

## Load run settings
# DIVA_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_J_run.rds")
DIVA_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_J_run.rds")
# DIVA_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_run.rds")
# DIVA_J_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_run.rds")

## Use MLE of DEC model as starting values
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")
DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_fit.rds")

## Set starting values for optimization based on MLE in the DEC model
d_start <- DIVA_fit$outputs@params_table["d","est"]
e_start <- DIVA_fit$outputs@params_table["e","est"]

DIVA_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_J_run.rds")
saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_J_run.rds")
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_run.rds")
# saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DIVA_J_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DIVA_J_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DIVA_J_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_run_log.txt", open = "wt")

sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DIVALIKE+J model run\n"))
Start_time <- Sys.time()
DIVA_J_fit <- bears_optim_run(DIVA_J_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DIVALIKE+J model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

print(DIVA_J_fit)

# Save model fit
# saveRDS(object = DIVA_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_J_fit.rds")
saveRDS(object = DIVA_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_J_fit.rds")
# saveRDS(object = DIVA_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_fit.rds")
# saveRDS(object = DIVA_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_fit.rds")

## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DIVA_J_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

# Inspect uncertainty in parameter estimates
DIVA_J_pars_table <- compute_parameter_CI_from_Hessian(DIVA_J_fit, subset_parameter_indices = c(1,3), plot = T, model_name = "DIVALIKE+J") ; DIVA_J_pars_table


# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DIVA_J_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DIVA_J_fit$inputs, include_null_range = DIVA_J_fit$inputs$include_null_range)
DIVA_J_Qmat <- DIVA_J_mat$Qmat
colnames(DIVA_J_Qmat) <- row.names(DIVA_J_Qmat) <- DIVA_J_mat$ranges_list
DIVA_J_Qmat

# Inspect final parameter table
DIVA_J_fit$outputs

# Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
# No subset sympatric speciation (s = 0)
# With jump dispersal (j)
# Here, relative weights of y and v are fixed to 1 - j/2.
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

# Inspect marginal likelihoods of ancestral states
DIVA_J_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DIVA_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DIVA_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.5/ Run DEC+W model ####

## Load run settings
# DEC_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DEC_W_run.rds")
DEC_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_W_run.rds")

## Use MLE of DEC model as starting values
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_fit.rds")
DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_fit.rds")

## Set starting values for optimization based on MLE in the DEC model
d_start <- DEC_fit$outputs@params_table["d","est"]
e_start <- max(DEC_fit$outputs@params_table["e","est"], 0.01)

DEC_W_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_W_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_W_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_W_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

### Check what happen if using force_sparse == TRUE to save time in computation

# Save settings before run
# saveRDS(object = DEC_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_W_run.rds")
saveRDS(object = DEC_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_W_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DEC_W_run$num_cores_to_use <- 1

# # Fit model using ML while recording logs
# # log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DEC_W_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DEC_W_run_log.txt", open = "wt")
# sink(file = log_file, type = "output", split = TRUE)
# cat(paste0(Sys.time(), " - Start DEC+W model run\n"))
# Start_time <- Sys.time()
# DEC_W_fit <- bears_optim_run(DEC_W_run)
# End_time <- Sys.time()
# cat(paste0(Sys.time(), " - End of DEC+W model run\n"))
# cat(paste0("Total run time of:\n"))
# print(End_time - Start_time)
# sink()
# closeAllConnections()
# 
# print(DEC_W_fit)
# 
# # Save model fit
# # saveRDS(object = DEC_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds")
# saveRDS(object = DEC_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_fit.rds")

### Loop to try several starting value for W 
w_start_list <- seq(0, 1, 0.25)
w_start_list[1] <- 0.001

for (i in seq_along(w_start_list))
{
  ## Load run settings
  # DEC_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DEC_W_run.rds")
  DEC_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_W_run.rds")
  
  # No parallelization as it seems to be useful only when the number of states is big
  # DEC_W_run$num_cores_to_use <- 1
  
  w_start <- w_start_list[i]
  
  # Inform changes in initial value of W
  DEC_W_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
  DEC_W_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization
  
  # Fit model using ML while recording logs
  # log_file <- file(paste0("./outputs/BioGeoBEARS_models/model_logs/DEC_W_",w_start,"_run_log.txt"), open = "wt")
  log_file <- file(paste0("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DEC_W_",w_start,"_run_log.txt"), open = "wt")
  sink(file = log_file, type = "output", split = TRUE)
  cat(paste0(Sys.time(), " - Start DEC+W model run - W start = ",w_start,"\n"))
  Start_time <- Sys.time()
  DEC_W_fit <- bears_optim_run(DEC_W_run)
  End_time <- Sys.time()
  cat(paste0(Sys.time(), " - End of DEC+W model run - W start = ",w_start,"\n"))
  cat(paste0("Total run time of:\n"))
  print(End_time - Start_time)
  sink()
  closeAllConnections()
  
  print(DEC_W_fit)
  
  # Save model fit
  # saveRDS(object = DEC_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/DEC_W_",w_start,"_fit.rds"))
  saveRDS(object = DEC_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_",w_start,"_fit.rds"))
  
}

# Best run is with w_start = 1 for roughly calibrated tree
# DEC_W_fit <- readRDS(file = paste0("./outputs/BioGeoBEARS_models/model_fits/DEC_W_1_fit.rds"))

# Best run is with w_start = 0.25 for MCC tree
DEC_W_fit <- readRDS(file = paste0("./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_0.25_fit.rds"))

# Save model fit
# saveRDS(object = DEC_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds")
saveRDS(object = DEC_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_fit.rds")


## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DEC_W_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and w (exponent on dispersal multipliers)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w

# Inspect uncertainty in parameter estimates
DEC_W_pars_table <- compute_parameter_CI_from_Hessian(DEC_W_fit, plot = T, model_name = "DEC+W") ; DEC_W_pars_table


# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DEC_W_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_W_fit$inputs, include_null_range = DEC_W_fit$inputs$include_null_range)
DEC_W_Qmat <- DEC_W_mat$Qmat
colnames(DEC_W_Qmat) <- row.names(DEC_W_Qmat) <- DEC_W_mat$ranges_list
DEC_W_Qmat

# Inspect final parameter table
DEC_W_fit$outputs

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# No jump dispersal (j = 0)
# Here, relative weights of y,v,s are fixed to 1.
# Model as 3 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers)

# Inspect marginal likelihoods of ancestral states
DEC_W_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DEC_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_W_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.6/ Run DIVALIKE+W model ####

## Load run settings
# DIVA_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_W_run.rds")
DIVA_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_run.rds")

## Use MLE of DEC model as starting values
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")
DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_fit.rds")

## Set starting values for optimization based on MLE in the DEC model
d_start <- DIVA_fit$outputs@params_table["d","est"]
e_start <- max(DIVA_fit$outputs@params_table["e","est"], 0.01)

DIVA_W_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_W_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_W_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_W_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

### Check what happen if using force_sparse == TRUE to save time in computation

# Save settings before run
# saveRDS(object = DIVA_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_W_run.rds")
saveRDS(object = DIVA_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DIVA_W_run$num_cores_to_use <- 1

# # Fit model using ML while recording logs
# # log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DIVA_W_run_log.txt", open = "wt")
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_run_log.txt", open = "wt")
# sink(file = log_file, type = "output", split = TRUE)
# cat(paste0(Sys.time(), " - Start DIVALIKE+W model run\n"))
# Start_time <- Sys.time()
# DIVA_W_fit <- bears_optim_run(DIVA_W_run)
# End_time <- Sys.time()
# cat(paste0(Sys.time(), " - End of DIVALIKE+W model run\n"))
# cat(paste0("Total run time of:\n"))
# print(End_time - Start_time)
# sink()
# closeAllConnections()
# 
# print(DIVA_W_fit)
# 
# # Save model fit
# # saveRDS(object = DIVA_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds")
# saveRDS(object = DIVA_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_fit.rds")

### Loop to try several starting value for W 
w_start_list <- seq(0, 1, 0.25)
w_start_list[1] <- 0.001

for (i in seq_along(w_start_list))
{
  ## Load run settings
  # DIVA_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_W_run.rds")
  DIVA_W_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_run.rds")
  
  # No parallelization as it seems to be useful only when the number of states is big
  # DIVA_W_run$num_cores_to_use <- 1
  
  w_start <- w_start_list[i]
  
  # Inform changes in initial value of W
  DIVA_W_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
  DIVA_W_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization
  
  # Fit model using ML while recording logs
  # log_file <- file(paste0("./outputs/BioGeoBEARS_models/model_logs/DIVA_W_",w_start,"_run_log.txt"), open = "wt")
  log_file <- file(paste0("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DIVA_W_",w_start,"_run_log.txt"), open = "wt")
  sink(file = log_file, type = "output", split = TRUE)
  cat(paste0(Sys.time(), " - Start DIVALIKE+W model run - W start = ",w_start,"\n"))
  Start_time <- Sys.time()
  DIVA_W_fit <- bears_optim_run(DIVA_W_run)
  End_time <- Sys.time()
  cat(paste0(Sys.time(), " - End of DIVALIKE+W model run - W start = ",w_start,"\n"))
  cat(paste0("Total run time of:\n"))
  print(End_time - Start_time)
  sink()
  closeAllConnections()
  
  print(DIVA_W_fit)
  
  # Save model fit
  # saveRDS(object = DIVA_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/DIVA_W_",w_start,"_fit.rds"))
  saveRDS(object = DIVA_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_",w_start,"_fit.rds"))
  
}

# Load best model with w_start = 0.25
# DIVA_W_fit <- readRDS(file = paste0("./outputs/BioGeoBEARS_models/model_fits/DIVA_W_0.25_fit.rds"))
DIVA_W_fit <- readRDS(file = paste0("./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_0.25_fit.rds"))
DIVA_W_fit$total_loglikelihood

# Save best model for DIVA+W
# saveRDS(object = DIVA_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds"))
saveRDS(object = DIVA_W_fit, file = paste0("./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_fit.rds"))


## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DIVA_W_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and w (exponent on dispersal multipliers)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w

# Inspect uncertainty in parameter estimates
DIVA_W_pars_table <- compute_parameter_CI_from_Hessian(DIVA_W_fit, plot = T, model_name = "DIVALIKE+W") ; DIVA_W_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DIVA_W_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DIVA_W_fit$inputs, include_null_range = DIVA_W_fit$inputs$include_null_range)
DIVA_W_Qmat <- DIVA_W_mat$Qmat
colnames(DIVA_W_Qmat) <- row.names(DIVA_W_Qmat) <- DIVA_W_mat$ranges_list
DIVA_W_Qmat

# Inspect final parameter table
DIVA_W_fit$outputs

# Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
# No subset sympatric speciation (s = 0)
# No jump dispersal (j = 0)
# Here, relative weights of y and v are fixed to 1.
# Model as 3 free parameters = d (range extension), e (range contraction), and w (exponent on dispersal multipliers)

# Inspect marginal likelihoods of ancestral states
DIVA_W_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DIVA_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DIVA_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_W_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.7/ Run DEC+J+W model ####

## Load run settings
# DEC_JW_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DEC_JW_run.rds")
DEC_JW_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_JW_run.rds")

## Use MLE of DEC+J and DEC+W models as starting values
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
# DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds") # Best fit so far for w_start = 1
DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_fit.rds") # Best fit so far for w_start = 1

## Set starting values for optimization based on MLE in the DEC+J model
d_start <- DEC_J_fit$outputs@params_table["d","est"]
e_start <- max(DEC_J_fit$outputs@params_table["e","est"], 0.01)
j_start <- max(DEC_J_fit$outputs@params_table["j","est"], 0.0001)
w_start <- max(DEC_W_fit$outputs@params_table["w","est"], 0.0001)

DEC_JW_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_JW_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization
DEC_JW_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization

### Check what happen if using force_sparse == TRUE to save time in computation

# Save settings before run
# saveRDS(object = DEC_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_JW_run.rds")
saveRDS(object = DEC_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DEC_JW_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DEC_JW_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DEC_JW_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DEC_JW_run_log.txt", open = "wt")
sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DEC+J+W model run\n"))
Start_time <- Sys.time()
DEC_JW_fit <- bears_optim_run(DEC_JW_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DEC+J+W model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

print(DEC_JW_fit)

# Save model fit
# saveRDS(object = DEC_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_JW_fit.rds")
saveRDS(object = DEC_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_JW_fit.rds")


## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DEC_JW_fit$optim_result
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w
# p4 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

# Inspect uncertainty in parameter estimates
DEC_JW_pars_table <- compute_parameter_CI_from_Hessian(DEC_JW_fit, plot = T, model_name = "DEC+J+W") ; DEC_JW_pars_table


# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DEC_JW_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_JW_fit$inputs, include_null_range = DEC_JW_fit$inputs$include_null_range)
DEC_JW_Qmat <- DEC_JW_mat$Qmat
colnames(DEC_JW_Qmat) <- row.names(DEC_JW_Qmat) <- DEC_JW_mat$ranges_list
DEC_JW_Qmat

# Inspect final parameter table
DEC_JW_fit$outputs

# Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
# With jump dispersal (j)
# Here, relative weights of y,v,s are fixed to 1 - j/3.
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)

# Inspect marginal likelihoods of ancestral states
DEC_JW_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DEC_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DEC_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_JW_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.8/ Run DIVALIKE+J+W model ####

## Load run settings
#  DIVA_JW_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_JW_run.rds")
DIVA_JW_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_run.rds")

## Use MLE of DIVALIKE+J and DIVALIKE+W models as starting values
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_J_fit.rds")
DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_J_fit.rds")
# DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds") # Best fit so far for w_start = 0.25
DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_fit.rds")  # Best fit so far for w_start = 0.25

## Set starting values for optimization based on MLE in the DEC+J model
d_start <- DIVA_J_fit$outputs@params_table["d","est"]
e_start <- max(DIVA_J_fit$outputs@params_table["e","est"], 0.01)
j_start <- max(DIVA_J_fit$outputs@params_table["j","est"], 0.0001)
w_start <- max(DIVA_W_fit$outputs@params_table["w","est"], 0.0001)

DIVA_JW_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_JW_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization
DIVA_JW_run$BioGeoBEARS_model_object@params_table["w","init"] <- w_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["w","est"] <- w_start # MLE will evolved after optimization

### Check what happen if using force_sparse == TRUE to save time in computation

# Save settings before run
# saveRDS(object = DIVA_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_JW_run.rds")
saveRDS(object = DIVA_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_run.rds")

# No parallelization as it seems to be useful only when the number of states is big
# DIVA_JW_run$num_cores_to_use <- 1

# Fit model using ML while recording logs
# log_file <- file("./outputs/BioGeoBEARS_models/model_logs/DIVA_JW_run_log.txt", open = "wt")
log_file <- file("./outputs/BioGeoBEARS_models/model_logs/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_run_log.txt", open = "wt")
sink(file = log_file, type = "output", split = TRUE)
cat(paste0(Sys.time(), " - Start DIVALIKE+J+W model run\n"))
Start_time <- Sys.time()
DIVA_JW_fit <- bears_optim_run(DIVA_JW_run)
End_time <- Sys.time()
cat(paste0(Sys.time(), " - End of DIVALIKE+J+W model run\n"))
cat(paste0("Total run time of:\n"))
print(End_time - Start_time)
sink()
closeAllConnections()

print(DIVA_JW_fit)

# Save model fit
# saveRDS(object = DIVA_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_JW_fit.rds")
saveRDS(object = DIVA_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_fit.rds")


## Inspect results

?optimx
# fevals = number of iteration during optimization
# convcode = 0 => Successful convergence
# convcode = 1 => Hit maxit (50 x nb of free parameters) without convergence (Should adjust maxit as )

DIVA_JW_fit$optim_result
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w
# p4 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

# Inspect uncertainty in parameter estimates
DIVA_JW_pars_table <- compute_parameter_CI_from_Hessian(DEC_JW_fit, plot = T, model_name = "DIVALIKE+J+W") ; DIVA_JW_pars_table

# Extract Q matrix (only anagenetic transitions. Not affected by multipliers and distances...)
DIVA_JW_mat <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DIVA_JW_fit$inputs, include_null_range = DIVA_JW_fit$inputs$include_null_range)
DIVA_JW_Qmat <- DIVA_JW_mat$Qmat
colnames(DIVA_JW_Qmat) <- row.names(DIVA_JW_Qmat) <- DIVA_JW_mat$ranges_list
DIVA_JW_Qmat

# Inspect final parameter table
DIVA_JW_fit$outputs

# Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
# No subset sympatric speciation (s = 0)
# With jump dispersal (j)
# Here, relative weights of y and v are fixed to 1 - j/2.
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)

# Inspect marginal likelihoods of ancestral states
DIVA_JW_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
dim(DEC_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node) # Likelihood df: Rows = tips + nodes, Columns = states/ranges
DIVA_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_JW_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


##### 7/ Plot ancestral range estimates #####

?BioGeoBEARS::plot_BioGeoBEARS_results

### 7.1/ For DEC model ####

# Load DEC model output
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_fit.rds")
DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_fit.rds")
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_fit.rds")
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DEC_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DEC_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/DEC_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/DEC_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_fit,
                         analysis_titletxt = "DEC - 2p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.2/ For DIVALIKE model ####

# Load DIVALIKE model output
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")
DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DIVA_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DIVA_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/DIVA_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/DIVA_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_fit,
                         analysis_titletxt = "DIVALIKE - 2p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.3/ For DEC+J model ####

# Load DEC+J model output
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_J_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_J_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/DEC_J_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/DEC_J_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_J_fit,
                         analysis_titletxt = "DEC+J - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.4/ For DIVALIKE+J model ####

# Load DIVALIKE+J model output
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_J_fit.rds")
DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_J_fit.rds")
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_fit.rds")
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DIVA_J_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DIVA_J_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_ARE_map.pdf", width = 40, height = 150)
# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_J_fit,
                         analysis_titletxt = "DIVALIKE+J - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.5/ For DEC+W model ####

# Load DEC+W model output
# DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds")
DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DEC_W_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DEC_W_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_W_fit,
                         analysis_titletxt = "DEC+W - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.6/ For DIVALIKE+W model ####

# Load DIVALIKE+W model output
# DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds")
DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DIVA_W_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DIVA_W_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_W_fit,
                         analysis_titletxt = "DIVALIKE+W - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.7/ For DEC+J+W model ####

# Load DEC+J+W model output
# DEC_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_JW_fit.rds")
DEC_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_JW_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DEC_JW_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DEC_JW_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_JW_fit,
                         analysis_titletxt = "DEC+J+W - 4p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.8/ For DIVALIKE+J+W model ####

# Load DIVALIKE+J+W model output
# DIVA_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_JW_fit.rds")
DIVA_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_fit.rds")

# pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/DIVA_JW_ARE_map.pdf", width = 40, height = 150)
pdf(file = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_ARE_map.pdf", width = 40, height = 150)

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(6.1,0.1,6.1,0.1), cex = 0.8) # bltr

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_JW_fit,
                         analysis_titletxt = "DIVALIKE+J+W - 4p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotsplits = FALSE, # To not plot corner/split states
                         plotlegend = TRUE, tipcex = 0.7, 
                         titlecex = 3.0,
                         statecex = 0.3, splitcex = 0.3)
dev.off()

### 7.9/ Aggregate ARE pdfs ####

# For roughly calibrated phylogeny
# all_ARE_maps_path <- list.files(path = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/", pattern = "ARE_map.pdf", full.names = T)
# qpdf::pdf_combine(input = all_ARE_maps_path, output = "./outputs/Ancestral_range_estimates_maps/Ponerinae_rough_phylogeny_1534t/All_models_ARE_maps.pdf")

# For MCC tree
all_ARE_maps_path <- list.files(path = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/", pattern = "ARE_map.pdf", full.names = T)
qpdf::pdf_combine(input = all_ARE_maps_path, output = "./outputs/Ancestral_range_estimates_maps/Ponerinae_MCC_phylogeny_1534t/All_models_ARE_maps.pdf")

# # For Youngest tree
# all_ARE_maps_path <- list.files(path = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/", pattern = "ARE_map.pdf", full.names = T)
# qpdf::pdf_combine(input = all_ARE_maps_path, output = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Youngest_phylogeny_1534t/All_models_ARE_maps.pdf")
# 
# # For Oldest tree
# all_ARE_maps_path <- list.files(path = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/", pattern = "ARE_map.pdf", full.names = T)
# qpdf::pdf_combine(input = all_ARE_maps_path, output = "./outputs/Ancestral_range_estimates_maps/Ponerinae_Oldest_phylogeny_1534t/All_models_ARE_maps.pdf")


##### 8/ Compare models #####

?BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object
?BioGeoBEARS::AICstats_2models

## Make my own version with model name, k, display parameters MLE, logLK, AIC, AICc, Akaike's weights

### 8.1/ Home-made functions to extract info from BioGeoBEARS model outputs ####

extract_nobs_from_BioGeoBEARS_output <- function (x)
{
  nobs <- (length(x$computed_likelihoods_at_each_node) + 1)/2
  names(nobs) <- "nobs"
  return(nobs)
}

extract_k_from_BioGeoBEARS_output <- function (x)
{
  k <- sum(x$inputs$BioGeoBEARS_model_object@params_table$type == "free")
  names(k) <- "k"
  return(k)
}

extract_MLE_params_from_BioGeoBEARS_output <- function (x)
{
  params <- x$outputs@params_table[c("d", "e", "y", "v", "s", "j"), "est"]
  names(params) <- c("d", "e", "y", "v", "s", "j")
  return(params)
}

extract_MLE_params_with_XW_from_BioGeoBEARS_output <- function (x)
{
  params <- x$outputs@params_table[c("d", "e", "y", "v", "s", "j","x","w"), "est"]
  names(params) <- c("d", "e", "y", "v", "s", "j","x","w")
  return(params)
}


# Home-made function to compute AIC from logLk and number of parameters (k)
compute_AIC <- function (logLk, k)
{
  AIC <- - 2 * (logLk - k)
}

# Home-made function to compute AICc from AIC, number of observations, and number of parameters (k)
compute_AICc <- function (AIC, nobs, k)
{
  AICc <- AIC + (2 * k) * (k - 1) / (nobs - k - 1)
}

### 8.2/ Build summary table for all models ####

# ## Load model fits for roughly calibrated tree
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_J_fit.rds")
# DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds")
# DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds")
# DEC_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_JW_fit.rds")
# DIVA_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_JW_fit.rds")

## Load model fits for MCC tree
DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_fit.rds")
DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_J_fit.rds")
DEC_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_W_fit.rds")
DIVA_W_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_W_fit.rds")
DEC_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_JW_fit.rds")
DIVA_JW_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DIVA_JW_fit.rds")

# ## Load model fits for Youngest tree
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_J_fit.rds")
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DIVA_J_fit.rds")

# ## Load model fits for Oldest tree
# DEC_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_fit.rds")
# DIVA_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_J_fit.rds")
# DIVA_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DIVA_J_fit.rds")


# Extract number of observations (terminal taxa)
# nb_obs <- length(Ponerinae_phylogeny_1534t_short_names$tip.label)
nb_obs <- length(Ponerinae_MCC_phylogeny_1534t_short_names$tip.label)
# nb_obs <- length(Ponerinae_Youngest_phylogeny_1534t_short_names$tip.label)
# nb_obs <- length(Ponerinae_Oldest_phylogeny_1534t_short_names$tip.label)

# Combine models' fit n a list
models_list <- list(DEC_fit, DIVA_fit, DEC_J_fit, DIVA_J_fit, DEC_W_fit, DIVA_W_fit, DEC_JW_fit, DIVA_JW_fit)
# models_list <- list(DEC_fit, DIVA_fit, DEC_J_fit, DIVA_J_fit)

# Extract ln-likelihood, number of parameters, and AIC from models' list
BIOGeoBEARS_models_comparison <- data.frame(model = c("DEC", "DIVALIKE", "DEC+J", "DIVALIKE+J", "DEC+W", "DIVALIKE+W", "DEC+J+W", "DIVALIKE+J+W"),
                                logLk = sapply(X = models_list, FUN = get_LnL_from_BioGeoBEARS_results_object),
                                k = sapply(X = models_list, FUN = extract_k_from_BioGeoBEARS_output))
# BIOGeoBEARS_models_comparison <- data.frame(model = c("DEC", "DIVALIKE", "DEC+J", "DIVALIKE+J"),
#                                             logLk = sapply(X = models_list, FUN = get_LnL_from_BioGeoBEARS_results_object),
#                                             k = sapply(X = models_list, FUN = extract_k_from_BioGeoBEARS_output))

# Compute AIC from logLk and k
BIOGeoBEARS_models_comparison$AIC <- unlist(purrr::map2(.x = BIOGeoBEARS_models_comparison$logLk, .y = BIOGeoBEARS_models_comparison$k, .f = compute_AIC))

# Compute AICc from AIC, k and nobs
BIOGeoBEARS_models_comparison$AICc <- unlist(purrr::map2(.x = BIOGeoBEARS_models_comparison$AIC, nobs = nb_obs, .y = BIOGeoBEARS_models_comparison$k, .f = compute_AICc))

# Compute Delta AICc
best_AICc <- min(BIOGeoBEARS_models_comparison$AICc)
BIOGeoBEARS_models_comparison$delta_AICc <- BIOGeoBEARS_models_comparison$AICc - best_AICc

# Compute Akaike's weights from AIC
BIOGeoBEARS_models_comparison$Akaike_weights <- round(phytools::aic.w(BIOGeoBEARS_models_comparison$AICc), 3) * 100

# Extract MLE parameters
models_MLE_pars <- t(sapply(X = models_list, FUN = extract_MLE_params_with_XW_from_BioGeoBEARS_output))

# Merge with summary table
BIOGeoBEARS_models_comparison <- cbind(BIOGeoBEARS_models_comparison, models_MLE_pars)

# Display results
BIOGeoBEARS_models_comparison

# Save results
# saveRDS(object = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison.rds")
saveRDS(object = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_MCC_phylogeny_1534t.rds")
# saveRDS(object = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_Youngest_phylogeny_1534t.rds")
# saveRDS(object = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_Oldest_phylogeny_1534t.rds")

# Export results
# write.xlsx(x = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison.xlsx")
write.xlsx(x = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_MCC_phylogeny_1534t.xlsx")
# write.xlsx(x = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_Youngest_phylogeny_1534t.xlsx")
# write.xlsx(x = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison_for_Oldest_phylogeny_1534t.xlsx")


## Best model for 7/8 Time-Strata and 8 Bioregions is DEC+J



