##### Script 04bis: Reduced model in BioGeoBEARS #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Fit a reduced tree with DEC to figure out how obtain print information while running parallel on cluster

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

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Time-calibrated phylogeny(ies) ####

# Load imputed phylogeny with short names
Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")

# Extract Myopias subset
Myopias_MRCA_node <- ape::getMRCA(phy = Ponerinae_phylogeny_1534t_short_names, tip = c("Myopias_my09", "Myopias_breviloba"))
Myopias_phylogeny <- ape::extract.clade(phy = Ponerinae_phylogeny_1534t_short_names, node = Myopias_MRCA_node)

plot(Myopias_phylogeny)

# Extract list of terminal taxa
Myopias_taxa <- Myopias_phylogeny$tip.label

# Substract even more
Myopias_taxa <- c("Myopias_my09", "Myopias_shivalikensis", "Myopias_levigata", "Myopias_philippinensis", "Myopias_breviloba")

# Extract Myopias subset
Myopias_phylogeny <- ape::keep.tip(phy = Ponerinae_phylogeny_1534t_short_names, tip = Myopias_taxa)

plot(Myopias_phylogeny)

### 1.2/ Load biogeographic ranges per species ####

Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Extract only subset taxa
Myopias_bioregions_binary_table <- Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name %in% Myopias_taxa, ]

# Remove useless bioregions
Myopias_bioregions_binary_table <- Myopias_bioregions_binary_table[, c(TRUE, (apply(X = Myopias_bioregions_binary_table[,-1], MARGIN = 2, FUN = sum) != 0))]


##### 2/ Format data for BioGeoBEARS #####

### 2.1/ Phylogeny ####

# Export under newick format
write.tree(Myopias_phylogeny, file = "./outputs/Grafting_missing_taxa/Myopias_phylogeny.tree")

# Set path to Newick tree
path_to_tree <- "./outputs/Grafting_missing_taxa/Myopias_phylogeny.tree"
path_to_tree <- BioGeoBEARS::np(path_to_tree)


### 2.2/ Biogeographic range table ####

# Extract binary df of presence/absence
binary_df <- Myopias_bioregions_binary_table[, -1]

# Convert character strings into numerical factors
binary_df_num <- as.data.frame(apply(X = binary_df, MARGIN = 2, FUN = as.numeric))
binary_df_factors <- binary_df_num
binary_df_factors[] <- lapply(binary_df_factors, factor) 
row.names(binary_df_factors) <- Myopias_bioregions_binary_table$Current_name
row.names(binary_df_num) <- Myopias_bioregions_binary_table$Current_name

dim(binary_df_factors)

# Convert area names in unique letters (needed for Biogeographic stochastic mapping)
names(binary_df_num) <- c("A", "I", "E")

# Produce tipranges object from numeric df
Myopias_bioregions_tipranges_obj <- define_tipranges_object(tmpdf = binary_df_num)

# Produce Lagrange PHYLIP biogeographic data
save_tipranges_to_LagrangePHYLIP(tipranges_object = Myopias_bioregions_tipranges_obj,
                                 lgdata_fn = "./input_data/BioGeoBEARS_setup/Myopias_bioregions_tipranges.data",
                                 areanames = colnames(Myopias_bioregions_tipranges_obj@df))

# Set path to Lagrange PHYLIP biogeographic data
path_to_geo_data <- "./input_data/BioGeoBEARS_setup/Myopias_bioregions_tipranges.data"
path_to_geo_data <- BioGeoBEARS::np(path_to_geo_data)


##### 3/ Plot biogeographic tip data ####

# Set the colors we’ll use for plotting
colors <- setNames(object = replicate(ncol(binary_df_factors),
                                      setNames(c("white", "gray30"), 0:1),
                                      simplify = FALSE),
                   nm = colnames(binary_df_factors))

pdf(file = "./outputs/Myopias_with_ranges.pdf", width = 6, height = 4)

# Set plotting parameters
par(mar = c(0.1,0.1,0.1,0.1), oma = c(3,0,2,3)) # bltr
par(xpd = T)
# Graph presence/absence using plotTree.datamatrix
range_map <- phytools::plotTree.datamatrix(tree = Myopias_phylogeny,
                                           X = binary_df_factors,
                                           fsize = 0.8, yexp = 1.1,
                                           mar = c(0.1,0.1,0.1,8.1),
                                           header = TRUE, xexp = 1.25, colors = colors)
par(xpd = T)
# Get plot info in "last_plot.phylo"
plot_info <- get("last_plot.phylo", envir=.PlotPhyloEnv)

## Add time line
# Extract root age
root_age <- max(phytools::nodeHeights(Myopias_phylogeny))
# Define ticks
ticks_labels <- seq(from = 0, to = 10, by = 2)
axis(side = 1, pos = 0, at = (-1 * ticks_labels) + root_age, labels = ticks_labels, cex.axis = 1.5)
legend(x = root_age/2,
       y = -1, adj = 0,
       bty = "n", legend = "", title = "Time  [My]", title.cex = 1.5)

## Add a legend
legend(# x = plot_info$x.lim[2] - 10,
       # y = mean(plot_info$y.lim),
       inset = c(0, 0),
       x = "topleft",
       legend = c("Absence", "Presence"),
       pch = 22, pt.bg = c("white","gray30"), pt.cex =  1.8,
       cex = 1.2, bty = "n")
par(xpd = F)
dev.off()


##### 4/ Set hyperparameters for all runs #####

### 4.1/ Set the maximum number of bioregions occupied by a lineage at any time ####

# Extract the range size (nb of bioregion occupied by terminals)
range_size <- rowSums(dfnums_to_numeric(Myopias_bioregions_tipranges_obj@df))
current_max_range_size <- max(range_size) # Current max range size

# Extract number of areas (bioregions)
nb_areas <- ncol(Myopias_bioregions_tipranges_obj@df)

# Cannot be lower than the observed max range size at the tips (1)
current_max_range_size
# Cannot be higher than the number of bioregions (3)
nb_areas

# Set max range size
max_range_size <- 3

# Estimate the number of states needed to represent all possible combinations of bioregions
cladoRcpp::numstates_from_numareas(numareas = nb_areas, maxareas = max_range_size)


### 4.2/ Set time boundaries of strata ####

# Use seven time slices associated to geological periods as in Kawahara et al., 2023: 
# 0-5.33 My = Holocene + Pleistocene + Pliocene
# 5.33-23.03 = Miocene 
# 23.03-33.9 = Oligocene
# 33.9-56.0 = Eocene
# 56.0-66.0 = Paleocene
# 66.0-100.5 = Late Cretaceous
# 100.5-145.0 = Early Cretaceous (Root of Ponerinae = 113 My ?)

# Oldest boundary should be any value older than the root
time_boundaries <- c(5.33, 23.03)

# Write output file for time_boundaries
writeLines(text = as.character(time_boundaries), con = "./input_data/BioGeoBEARS_setup/time_boundaries_reduced.txt")

# Set path to boundaries for the time-periods
path_to_time_strata_boundaries <- "./input_data/BioGeoBEARS_setup/time_boundaries_reduced.txt"
path_to_time_strata_boundaries <- BioGeoBEARS::np(path_to_time_strata_boundaries)


### 4.3/ Set time-stratified adjacency matrices ####

# Provide list of bimorphic states allowed by time-frame as binary adjacency matrices
writeLines(readLines("./input_data/BioGeoBEARS_setup/Adjacency_matrix_reduced.txt"))
# Binary matrices are symmetrical
# Areas should be in the same order as in the Bioregion table!

# Works by removing all ranges/states that include combination of non-adjacent areas
# Better option is to use the manual list of allowed states as $lists_of_states_lists_0based

# Set path to time-stratified adjacency matrices
path_to_adjacency_matrices <- "./input_data/BioGeoBEARS_setup/Adjacency_matrix_reduced.txt"
path_to_adjacency_matrices <- BioGeoBEARS::np(path_to_adjacency_matrices)

# Import as list of df the time-stratified adjacency matrices
list_of_adjacency_matrices <- BioGeoBEARS:::read_areas_adjacency_fn(areas_adjacency_fn = path_to_adjacency_matrices)


## 4.5/ Check/Adjust lists of allowed states per time-strata using manual modification ####

## 4.5.1/ Get the list of monomorphic states = bioregion = areas from the tip data ####
areas_list <- getareas_from_tipranges_object(Myopias_bioregions_tipranges_obj)
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
length(valid_ranges_stratum_1) # Now 7 states
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
length(valid_ranges_stratum_2) # Now 5 states
setdiff(valid_ranges_stratum_2, valid_ranges_stratum_1) # Gain no states
setdiff(valid_ranges_stratum_1, valid_ranges_stratum_2) # Lose "AI" and "AIE" states
# Make a selection vector to keep only valid ranges  
keep_stratum_2 <- (ranges_list %in% valid_ranges_stratum_2)

# Inform the list of valid states
lists_of_states_lists_0based[[2]] <- ranges_list_0based[keep_stratum_2]


##### 5/ Set BioGeoBEARS models #####

## 5.1/ Create run object for Time-stratified DEC model ####

# Detect available number of cores
nb_cores <- parallel::detectCores()
nb_cores <- round(nb_cores*0.75)

# # Set a cluster
# ?parallel::makeCluster
# parallel::makeCluster(rep("localhost",nb_cores), type = "SOCK")


Myopias_DEC_run <- define_BioGeoBEARS_run(
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
  areas_allowed_fn = NA, # To provide path to file defining the allowed areas at a given time-frame.
  areas_adjacency_fn = path_to_adjacency_matrices, # To provide path to optional file defining the allowed bimorphic states/connections between bioregions at a given time-frame.
  return_condlikes_table = TRUE, # To ask to obtain all marginal likelihoods computed by the model and used to display ancestral states
  # fixnode = NA, # To provide fixed states for specific nodes
  fixlikes = NA) # To provide weighting scheme for tip likelihoods if you want to account for uncertainty in tip ranges, instead of using 1 for the observed state and 0 for the others.

# Add manually the list of allowed states (per time-strata). Very useful to avoid overparametrization when authorizing a few carefully selected 'polymorphic' ranges
Myopias_DEC_run$lists_of_states_lists_0based <- lists_of_states_lists_0based

## 5.2/ Create Time-strata ####

# Loads the time-periods from the text files into the model object.
Myopias_DEC_run <- readfiles_BioGeoBEARS_run(Myopias_DEC_run)

# Divide the tree up by timeperiods/strata
Myopias_DEC_run <- BioGeoBEARS::section_the_tree(inputs = Myopias_DEC_run,
                                         make_master_table = TRUE,
                                         plot_pieces = FALSE,
                                         fossils_older_than = 0.001,
                                         cut_fossils = FALSE)
# The stratified tree is summarized in this table describing stratum membership per branches
Myopias_DEC_run$master_table

# Check information provided as input to BioGeoBEARS for our run
str(Myopias_DEC_run, max.level = 2)

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

Myopias_DEC_run$BioGeoBEARS_model_object@params_table

# Check that starting parameter values are inside the min/max
Myopias_DEC_run <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = Myopias_DEC_run)
# Check validity of set-up before run
check_BioGeoBEARS_run(Myopias_DEC_run)

# Save settings before run
saveRDS(object = Myopias_DEC_run, file = "./outputs/BioGeoBEARS_models/model_runs/Myopias_DEC_run.rds")


##### 6/ Run models #####

# Make several runs with different starting parameter values to ensure optimal likelihood has been reached
# Use the stepwise approach by using MLE estimates of previous runs as staring parameter values when increasing complexity (adding parameters)

# Start with DIVALIKE and DEC
# Then DEC + J, DEC + W (and DIVALIKE + J, DIVALIKE + W)
# Then DEC + J + W (and DIVALIKE + J + W) using MLE from +J models 

# Compare LnLk such as bigger model always have better (higher) LnLk
# Use annealing if issues with optimization (slower but less prone to suboptimal pitfalls) BioGeoBEARS_run_object$use_optimx = "GenSA"

### 6.1/ Run DEC model ####

# Load model run
Myopias_DEC_run <- readRDS(file = "./outputs/BioGeoBEARS_models/model_runs/Myopias_DEC_run.rds")

# No parallelization
# Myopias_DEC_run$num_cores_to_use <- 18
# Myopias_DEC_run$num_cores_to_use <- 2

# Fit model using ML
Myopias_DEC_fit <- bears_optim_run(Myopias_DEC_run)

# Save model fit
saveRDS(object = Myopias_DEC_fit, file = "./outputs/BioGeoBEARS_models/model_fits/Myopias_DEC_fit.rds")


##### Look in the functions to find a way to make the likelihood optimization more verbose, even when using clusters for parallelization! #####

# Check what happen if using force_sparse == TRUE
# Should get a print to help estimating the time of computation (but may not work once the cluster is open)


# Inspect results
DEC_fit$optim_result
# Model as only 2 free parameters = d (range extension) and e (range contraction)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas

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
DEC_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.2/ Run DIVALIKE model ####

# Save settings before run
saveRDS(object = DIVA_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_run.rds")

# Fit model using ML
DIVA_fit <- bears_optim_run(DIVA_run)

# Save model fit
saveRDS(object = DIVA_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_fit.rds")


# Inspect results
DIVA_fit$optim_result
# Model has 2 free parameters = d (range extension), e (range contraction)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas

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
str(DIVA_fit)
DIVA_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DIVA model with the MLE parameters
DIVA_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DIVA model with the MLE parameters estimated globally
DIVA_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.3/ Run DEC+J model ####

## Use MLE of DEC model as starting values

## Set starting values for optimization based on MLE in the DEC model
d_start <- DEC_fit$outputs@params_table["d","est"]
e_start <- DEC_fit$outputs@params_table["e","est"]

DEC_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DEC_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_J_run.rds")

# Fit model using ML
DEC_J_fit <- bears_optim_run(DEC_J_run)

# Save model fit
saveRDS(object = DEC_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")


# Inspect results
DEC_J_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

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
str(DEC_J_fit)
DEC_J_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DEC_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.4/ Run DIVALIKE+J model ####

## Use MLE of DIVALIKE model as starting values

## Set starting values for optimization based on MLE in the DEC model
d_start <- DIVA_fit$outputs@params_table["d","est"]
e_start <- DIVA_fit$outputs@params_table["e","est"]

DIVA_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DIVA_J_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_J_run.rds")

# Fit model using ML
DIVA_J_fit <- bears_optim_run(DIVA_J_run)

# Save model fit
saveRDS(object = DIVA_J_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_J_fit.rds")


# Inspect results
DIVA_J_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

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
str(DIVA_J_fit)
DIVA_J_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DIVA_J_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_J_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.5/ Run DEC+W model ####

## Use MLE of DEC model as starting values

## Set starting values for optimization based on MLE in the DEC model
d_start <- DEC_fit$outputs@params_table["d","est"]
e_start <- DEC_fit$outputs@params_table["e","est"]

DEC_W_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_W_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_W_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_W_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DEC_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_W_run.rds")

# Fit model using ML
DEC_W_fit <- bears_optim_run(DEC_W_run)

# Save model fit
saveRDS(object = DEC_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_W_fit.rds")


# Inspect results
DEC_W_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and w (exponent on dispersal multipliers)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w

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
str(DEC_W_fit)
DEC_W_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DEC_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_W_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.6/ Run DIVALIKE+W model ####

## Use MLE of DIVALIKE model as starting values

## Set starting values for optimization based on MLE in the DEC model
d_start <- DIVA_fit$outputs@params_table["d","est"]
e_start <- DIVA_fit$outputs@params_table["e","est"]

DIVA_W_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_W_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_W_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_W_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DIVA_W_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_W_run.rds")

# Fit model using ML
DIVA_W_fit <- bears_optim_run(DIVA_W_run)

# Save model fit
saveRDS(object = DIVA_W_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_W_fit.rds")


# Inspect results
DIVA_W_fit$optim_result
# Model as 3 free parameters = d (range extension), e (range contraction), and w (exponent on dispersal multipliers)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w

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
str(DIVA_W_fit)
DIVA_W_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DIVA_W_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_W_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).


### 6.7/ Run DEC+J+W model ####

## Use MLE of DEC+J model as starting values

## Set starting values for optimization based on MLE in the DEC+J model
d_start <- DEC_J_fit$outputs@params_table["d","est"]
e_start <- DEC_J_fit$outputs@params_table["e","est"]
j_start <- DEC_J_fit$outputs@params_table["j","est"]

DEC_JW_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DEC_JW_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DEC_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DEC_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DEC_JW_run.rds")

# Fit model using ML
DEC_JW_fit <- bears_optim_run(DEC_JW_run)

# Save model fit
saveRDS(object = DEC_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DEC_JW_fit.rds")


# Inspect results
DEC_JW_fit$optim_result
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w
# p4 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

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
str(DEC_JW_fit)
DEC_JW_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DEC_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DEC_JW_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).

### 6.8/ Run DIVALIKE+J+W model ####

## Use MLE of DIVALIKE+J model as starting values

## Set starting values for optimization based on MLE in the DEC model
d_start <- DIVA_J_fit$outputs@params_table["d","est"]
e_start <- DIVA_J_fit$outputs@params_table["e","est"]
j_start <- DIVA_J_fit$outputs@params_table["j","est"]

DIVA_JW_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
DIVA_JW_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
DIVA_JW_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

# Save settings before run
saveRDS(object = DIVA_JW_run, file = "./outputs/BioGeoBEARS_models/model_runs/DIVA_JW_run.rds")

# Fit model using ML
DIVA_JW_fit <- bears_optim_run(DIVA_JW_run)

# Save model fit
saveRDS(object = DIVA_JW_fit, file = "./outputs/BioGeoBEARS_models/model_fits/DIVA_JW_fit.rds")


# Inspect results
DIVA_JW_fit$optim_result
# Model as 4 free parameters = d (range extension), e (range contraction), w (exponent on dispersal multipliers), and jump-dispersal (j)
# p1 = d = dispersal rate = rate of range expansion between areas along branches = transition from a range with N-areas to a range with N+1 areas
# p2 = e = extinction rate within each area along branches = transition from a range with N-areas to a range with N-1 areas
# p3 = w = exponent on dispersal multipliers (modifies d, j) = controls the relationship between d/j and dispersal multipliers such as d' = multiplier^w
# p4 = j = jump dispersal relative weight = relative weights of cladogenetic founder-event. Ex: A -> (A),(B)

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
str(DIVA_JW_fit)
DIVA_JW_fit$total_loglikelihood # log-likelihood of the geographic range data given the phylogenetic tree and DEC model with the MLE parameters
DIVA_JW_fit$ML_marginal_prob_each_state_at_branch_top_AT_node # Scaled marginal likelihood (probabilities) at each node/tip given the phylogenetic tree and DEC model with the MLE parameters estimated globally
DIVA_JW_fit$ML_marginal_prob_each_state_at_branch_bottom_below_node # Same but for at the 'corner' of internal node, meaning at the start of each descendant branch, after cladogenesis/speciation (and eventually cladogenetic transition).

##### 7/ Plot ancestral range estimates #####

?BioGeoBEARS::plot_BioGeoBEARS_results

### 7.1/ For DEC model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_fit,
                         analysis_titletxt = "DEC - 2p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.2/ For DIVALIKE model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_fit,
                         analysis_titletxt = "DIVALIKE - 2p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.3/ For DEC+J model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_J_fit,
                         analysis_titletxt = "DEC+J - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.4/ For DIVALIKE+J model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_J_fit,
                         analysis_titletxt = "DIVALIKE+J - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.5/ For DEC+W model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_W_fit,
                         analysis_titletxt = "DEC+W - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.6/ For DIVALIKE+W model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_W_fit,
                         analysis_titletxt = "DIVALIKE+W - 3p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.7/ For DEC+J+W model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DEC_JW_fit,
                         analysis_titletxt = "DEC+J+W - 4p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

### 7.8/ For DIVALIKE+J+W model ####

# Subdivide our plotting area using layout
layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1), cex = 0.8)

# Plot legend and tree
plot_BioGeoBEARS_results(results_object = DIVA_JW_fit,
                         analysis_titletxt = "DIVALIKE+J+W - 4p model",
                         plotwhat = "pie", # If setting "test", will only plot the states with the highest marginal likelihood
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

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

# Combine models' fit n a list
models_list <- list(DEC_fit, DIVA_fit, DEC_J_fit, DIVA_J_fit, DEC_W_fit, DIVA_W_fit, DEC_JW_fit, DIVA_JW_fit)

# Extract ln-likelihood, number of parameters, and AIC from models' list
BIOGeoBEARS_models_comparison <- data.frame(model = c("DEC", "DIVALIKE", "DEC+J", "DIVALIKE+J", "DEC+W", "DIVALIKE+W", "DEC+J+W", "DIVALIKE+J+W",),
                                            logLk = sapply(X = models_list, FUN = get_LnL_from_BioGeoBEARS_results_object),
                                            k = sapply(X = models_list, FUN = extract_k_from_BioGeoBEARS_output))
# Compute AIC from logLk and nobs
BIOGeoBEARS_models_comparison$AIC <- unlist(purrr::map2(.x = BIOGeoBEARS_models_comparison$logLk, .y = BIOGeoBEARS_models_comparison$k, .f = compute_AIC))

# Compute AICc from AIC, k and nobs
BIOGeoBEARS_models_comparison$AICc <- unlist(purrr::map2(.x = BIOGeoBEARS_models_comparison$AIC, nobs = nobs, .y = BIOGeoBEARS_models_comparison$k, .f = compute_AICc))

# Compute Akaike's weights from AIC
BIOGeoBEARS_models_comparison$Akaike_weights <- round(phytools::aic.w(BIOGeoBEARS_models_comparison$AICc), 3) * 100

# Extract MLE parameters
models_MLE_pars <- t(sapply(X = models_list, FUN = extract_MLE_params_with_XW_from_BioGeoBEARS_output))

# Merge with summary table
BIOGeoBEARS_models_comparison <- cbind(BIOGeoBEARS_models_comparison, models_MLE_pars)

# Display results
BIOGeoBEARS_models_comparison

# Save results
saveRDS(object = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison.rds")
# Export results
write.xlsx(x = BIOGeoBEARS_models_comparison, file = "./outputs/BioGeoBEARS_models/BIOGeoBEARS_models_comparison.xlsx")



