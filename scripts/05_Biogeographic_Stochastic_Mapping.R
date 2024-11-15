##### Script 05: Biogeographic Stochastic Mapping #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Simulate Biogeographic histories under given tree, tip data, model and parameter MLE, and node ACE.
# Aggregate results under biogeographic stochastic maps
  # Maps of continuous posterior probability of presence in each bioregion


## Other scripts

## Diversification dynamics
# Compare diversification rates across bioregions/clades (from rjMCMC / BAMM)
  # For bioregion membership, use the continuous probability from the aggregated BSM

###

### Inputs

# (Set of) Time-calibrated phylogeny/ies with imputed missing taxa
# Best fitted Biogeographic model

###

### Outputs

# 1000 Biogeographic stochastic maps
# Aggregated maps displaying the continuous probability of presence in each bioregion (1 per bioregion)

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(phytools)
library(geiger)
library(ape)
library(BioGeoBEARS)
library(qpdf)
library(MultinomialCI)    # For 95% CIs on BSM counts
library(magick)   # For animated GIF
library(pdftools)  # To read PDF
library(BayesTwin) # To compute HPD intervals

# devtools::install_github("nmatzke/BioGeoBEARS")


### 1.2/ Load data and path to data ####

# Load imputed phylogeny with short names
# Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")
Ponerinae_MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")

# Set path to Newick tree
# path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.tree"
path_to_tree <- "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.tree"
path_to_tree <- BioGeoBEARS::np(path_to_tree)


# Load updated range table
Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Set path to Lagrange PHYLIP biogeographic data
path_to_geo_data <- "./input_data/BioGeoBEARS_setup/lagrange_area_data_file_7_regions_PaleA.data"
path_to_geo_data <- BioGeoBEARS::np(path_to_geo_data)

# Read data from file
Taxa_bioregions_tipranges_obj <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = path_to_geo_data)


### 1.3/ Load modeling results from best model (DEC+J) #####

# Load time-stratified DEC+J model output
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")


##### 2/ Run Biogeographic Stochastic Mapping ####

?BioGeoBEARS::get_inputs_for_stochastic_mapping
?BioGeoBEARS::runBSM


### 2.1/ Extract inputs needed for BSM from model fit object ####

DEC_J_BSM_inputs <- get_inputs_for_stochastic_mapping(res = DEC_J_fit)

# Save BSM inputs
# saveRDS(object = DEC_J_BSM_inputs, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_inputs.rds")
saveRDS(object = DEC_J_BSM_inputs, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_inputs.rds")
# DEC_J_BSM_inputs <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_inputs.rds")
DEC_J_BSM_inputs <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_inputs.rds")

# BSM input object as one list per time-stratum
str(DEC_J_BSM_inputs, max.level = 1)
str(DEC_J_BSM_inputs[[1]], max.level = 2)

# Inspect inputs
names(DEC_J_BSM_inputs[[1]])
unlist(DEC_J_BSM_inputs[[1]]$ranges_list) # All states/ranges list (All included, even the actually non-allowed ones)
DEC_J_BSM_inputs[[1]]$numstates # Total number of states/ranges (All included, even the actually non-allowed ones)
DEC_J_BSM_inputs[[1]]$states_allowed_this_timeperiod_TF # T/F to select allowed states for this time-period
unlist(DEC_J_BSM_inputs[[1]]$ranges_list)[DEC_J_BSM_inputs[[1]]$states_allowed_this_timeperiod_TF] # List of allowed states/ranges for this time-period
DEC_J_BSM_inputs[[1]]$state_indices_0based_this_timeperiod # 0Based states/ranges. Only the allowed ones
DEC_J_BSM_inputs[[1]]$tipranges # Tip geographic data (binary ranges)
DEC_J_BSM_inputs[[1]]$Qmat # Transition Q matrix with d/e parameters. Only the allowed states ordered as in States/Ranges list
length(DEC_J_BSM_inputs[[1]]$independent_likelihoods_by_tree_piece_for_timeperiod_i) # Independent likelihood of each branch for each possible/allowed parental node state

# Transition Q matrix with d/e parameters. 
# States are ordered as in States/Ranges list
# They are different across time-strata as the trait spaces are different!
DEC_J_BSM_inputs[[1]]$Qmat 
DEC_J_BSM_inputs[[3]]$Qmat
DEC_J_BSM_inputs[[5]]$Qmat

### 2.2/ Run BSM  ####

# nb_BSM_maps <- 10
nb_BSM_maps <- 1000

DEC_J_BSM_output <- runBSM(res = DEC_J_fit,  # Model fit object
                             stochastic_mapping_inputs_list = DEC_J_BSM_inputs,
                             maxnum_maps_to_try = 200, # Maximum number of stochastic maps to simulate
                             nummaps_goal = nb_BSM_maps, # Number of accepted stochastic maps to target
                             maxtries_per_branch = 40000,
                             save_after_every_try = TRUE,
                             # savedir = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/",
                             savedir = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/",
                             seedval = 12345,
                             wait_before_save = 0.01,
                             master_nodenum_toPrint = 0)

## Save BSM outputs
# saveRDS(object = DEC_J_BSM_output, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_output.rds")
saveRDS(object = DEC_J_BSM_output, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_output.rds")

## Inspect outputs = cladogenetic and anagenetic event tables

# Load BSM outputs
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_output.rds")
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_output.rds")

# Extract records of events
DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
DEC_J_ana_events_tables <- DEC_J_BSM_output$RES_ana_events_tables

length(DEC_J_clado_events_tables) # One summary table per stochastic history
length(DEC_J_ana_events_tables) # One summary table per stochastic history

# Save independently the cladogenetic and anagenetic tables
# saveRDS(object = DEC_J_clado_events_tables, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_clado_events_tables.rds")
saveRDS(object = DEC_J_clado_events_tables, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_clado_events_tables.rds")
# saveRDS(object = DEC_J_ana_events_tables, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_ana_events_tables.rds")
saveRDS(object = DEC_J_ana_events_tables, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_ana_events_tables.rds")


## Cladogenetic event tables (one table per simulated biogeographic history)

# One row = 1 node/tip
names(DEC_J_clado_events_tables[[1]])
nrow(DEC_J_clado_events_tables[[1]]) # Have more events than nodes/tips because some branches are split across time-strata and one "event" (or absence of event) per branch per strata is recorded
# Ponerinae_phylogeny_1534t_short_names$Nnode + length(Ponerinae_phylogeny_1534t_short_names$tip.label)
Ponerinae_MCC_phylogeny_1534t_short_names$Nnode + length(Ponerinae_MCC_phylogeny_1534t_short_names$tip.label)

View(DEC_J_clado_events_tables[[1]])

# Node metadata
DEC_J_clado_events_tables[[1]]$node # Node ID as in ape
DEC_J_clado_events_tables[[1]]$label # Node label
table(DEC_J_clado_events_tables[[1]]$node.type) # Tip, root or internal
DEC_J_clado_events_tables[[1]]$parent_br # ID of the edge leading to the node = gives you the row of tree$edge and tree$edge.length)
DEC_J_clado_events_tables[[1]]$edge.length # Length of the branch leading to the node
DEC_J_clado_events_tables[[1]]$ancestor # ID of the node at the bottom of the branch leading to the node
DEC_J_clado_events_tables[[1]]$daughter_nds # ID of the nodes descending from the two descending lineages of the node
DEC_J_clado_events_tables[[1]]$tipnames # Comma separated list of all descending terminals
DEC_J_clado_events_tables[[1]]$node_ht # Time since the root
DEC_J_clado_events_tables[[1]]$time_bp # Time before present of the node

# States description. States are sampled using the MLE scaled marginal likelihoods
DEC_J_clado_events_tables[[1]]$sampled_states_AT_nodes # Node state in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$sampled_states_AT_brbots # State at the start of the branch leading to the node in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$samp_LEFT_dcorner # State at the start of the left descendant edge of the node in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$samp_RIGHT_dcorner # State at the start of the right descendant edge of the node in this particular history / stochastic map

# Event description
DEC_J_clado_events_tables[[1]]$clado_event_type # Type of cladogenetic event
DEC_J_clado_events_tables[[1]]$clado_event_txt # Text format describing the event using states
DEC_J_clado_events_tables[[1]]$clado_dispersal_to # In case of cladogenetic dispersal (i.e, founder event, j), record the area of destination

# Stratum metadata
DEC_J_clado_events_tables[[1]]$stratum # ID of the time-stratum where is the node
DEC_J_clado_events_tables[[1]]$time_top # Tipward time boundary of the stratum where is the node
DEC_J_clado_events_tables[[1]]$time_bot # Rootward time boundary of the stratum where is the node
DEC_J_clado_events_tables[[1]]$reltimept # Time-width of the time-stratum where is the node
DEC_J_clado_events_tables[[1]]$piecenum # ID of the tree piece where the branch is within the time-stratum
DEC_J_clado_events_tables[[1]]$piececlass # Type of tree piece where the branch is within the time-stratum. Can be a segment of a single branch ('subbranch') or a 'subtree' if speciation occurs within the time-stratum

# Subtree pieces metadata (part of the tree within the time-stratum where the branch where the anagenetic event happened can be found)
DEC_J_clado_events_tables[[1]]$SUBnode
DEC_J_clado_events_tables[[1]]$SUBnode # ID of the descending node of the branch within the subpiece
table(DEC_J_clado_events_tables[[1]]$SUBnode.type) # Tip, root or internal branch within the subpiece
DEC_J_clado_events_tables[[1]]$SUBparent_br # ID of the parent branch within the subpiece (NA for roots)
DEC_J_clado_events_tables[[1]]$SUBedge.length # Length of the branch within the subpiece
DEC_J_clado_events_tables[[1]]$SUBancestor # ID of the node at the bottom of the branch within the subpiece
DEC_J_clado_events_tables[[1]]$SUBdaughter_nds # ID of the nodes descending from the two descending lineages of the node descending from the branch, within the subpiece
DEC_J_clado_events_tables[[1]]$SUBlabel # Label of the descending node within the subpiece
DEC_J_clado_events_tables[[1]]$SUBtime_bp # Time of the descending node before the end of the time-stratum 

## Anagenetic event tables (one table per simulated biogeographic history)

names(DEC_J_ana_events_tables[[1]])
nrow(DEC_J_ana_events_tables[[1]]) # Number of rows depend on the number of anagenetic events
View(DEC_J_ana_events_tables[[1]])

# Have similar columns for the metadata copied from the cladogenetic table based on the descending node of the branch where the event is happening

# Branch metadata
DEC_J_ana_events_tables[[1]]$parent_br # ID of the branch where the event is happening
DEC_J_ana_events_tables[[1]]$edge.length # Length of the branch where the event is happening
DEC_J_ana_events_tables[[1]]$ancestor # ID of the node at the bottom of the branch where the event is happening
DEC_J_ana_events_tables[[1]]$tipnames # Comma separated list of all descending terminals

# Descending node metadata
DEC_J_ana_events_tables[[1]]$node # ID of the descending node as in ape
DEC_J_ana_events_tables[[1]]$nodenum_at_top_of_branch # Same = ID of the descending node as in ape
DEC_J_ana_events_tables[[1]]$label # Label of the descending node
table(DEC_J_ana_events_tables[[1]]$node.type) # Type of the descending node: Tip or internal (No anagenetic event can happen on the root since it has no branch leading to it)
DEC_J_ana_events_tables[[1]]$daughter_nds # ID of the nodes descending from the two descending lineages of the descending node
DEC_J_ana_events_tables[[1]]$node_ht # Time since the root
DEC_J_ana_events_tables[[1]]$time_bp # Time before present of the node

# Have new columns to describe the anagenetic events
DEC_J_ana_events_tables[[1]]$trynum # Number of trials needed to generate/simulate an evolutionary history compatible with the starting and ending states of the branch
DEC_J_ana_events_tables[[1]]$brlen # Length of the branch where the event is happening OR the segment within the time-stratum, when stratified
DEC_J_ana_events_tables[[1]]$current_rangenum_1based # State before event in numerical format
DEC_J_ana_events_tables[[1]]$new_rangenum_1based # State after event in numerical format
DEC_J_ana_events_tables[[1]]$current_rangetxt # State before event in text format
DEC_J_ana_events_tables[[1]]$new_rangetxt # State after event in text format
DEC_J_ana_events_tables[[1]]$abs_event_time # Time before present of the event
DEC_J_ana_events_tables[[1]]$event_time # Relative time of the event on the branch (from the tipward node, or the end of the segment within the time-stratum)
DEC_J_ana_events_tables[[1]]$abs_event_time + DEC_J_ana_events_tables[[1]]$event_time # Should equal branch length if not time-stratified
DEC_J_ana_events_tables[[1]]$event_type # Type of event: range extension (d), range-switching (a), or range contraction (e)
DEC_J_ana_events_tables[[1]]$event_txt # Text format describing the event using states
DEC_J_ana_events_tables[[1]]$new_area_num_1based # New destination areas in case of d/a event, in numerical format
DEC_J_ana_events_tables[[1]]$lost_area_num_1based # Lost source areas in case of e event
DEC_J_ana_events_tables[[1]]$ana_dispersal_from # Source area selected probabilistically as the source using simulate_source_areas_ana_clado() in case the source range encompasses multiple area
DEC_J_ana_events_tables[[1]]$dispersal_to # Destination area in d/a event, in text format
DEC_J_ana_events_tables[[1]]$dispersal_to # Lost source area in case of e event, in text format


##### 3/ Plot all Biogeographic Stochastic Maps ####

# Extract number of stochastic maps
nb_maps <- length(DEC_J_clado_events_tables)

# Extract complete list of areas used across time-strata
all_areas <- lapply(X = DEC_J_BSM_inputs, FUN = function (x) x$areas)
all_areas <- unique(unlist(all_areas))

# Extract complete list of "character" states used across time-strata
all_states <- unlist(DEC_J_BSM_inputs[[1]]$ranges_list)

# Extract complete list of 0based states used across time-strata
all_0based_states <- DEC_J_BSM_inputs[[1]]$state_indices_0based_all_timeperiods

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
saveRDS(object = colors_list_for_states, file = "./outputs/BSM/colors_list_for_states.rds")

# Load path to additional script
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

### 3.1/ Loop through the maps and plot to a single PDF (typically fail) ####

# pdf(file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/DEC_J_BSMs_",nb_maps,"_maps.pdf"),
#     height = 150, width = 40)
pdf(file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/DEC_J_BSMs_",nb_maps,"_maps.pdf"),
    height = 150, width = 40)

for (i in 1:nb_maps)
# for (i in 1:100)
{
  # i <- 1
  
  # Extract the master tables of events summary for map n°i
  DEC_J_clado_events_table_i <- DEC_J_clado_events_tables[[i]]
  DEC_J_ana_events_table_i <- DEC_J_ana_events_tables[[i]]
  
  ## Extend the master table for cladogenetic events to account for branches that are split across different time-strata
  cols_to_get <- names(DEC_J_clado_events_table_i[,-ncol(DEC_J_clado_events_table_i)])
  colnums <- match(cols_to_get, names(DEC_J_ana_events_table_i))
  ana_events_table_cols_to_add <- DEC_J_ana_events_table_i[,colnums]
  anagenetic_events_txt_below_node <- rep("none", nrow(ana_events_table_cols_to_add))
  ana_events_table_cols_to_add <- cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  rows_to_get_TF <- ana_events_table_cols_to_add$node <= length(Ponerinae_phylogeny_1534t_short_names$tip.label)
  DEC_J_clado_events_master_table_i <- rbind(ana_events_table_cols_to_add[rows_to_get_TF,], DEC_J_clado_events_table_i)
  
  # View(DEC_J_clado_events_table_i)
  # View(DEC_J_clado_events_master_table_i)
  
  # Convert the BSM into a modified fit object
  DEC_J_BSM_map_i <- stochastic_map_states_into_res(res = DEC_J_fit,  # Model fit object
                                                      master_table_cladogenetic_events = DEC_J_clado_events_master_table_i,
                                                      stratified = TRUE)
  # Set title
  main_i <- paste0("\n\nDEC + J - Stochastic Map #", i, "/", nb_maps)
  
  # Plot stochastic map with node states
  plot_BioGeoBEARS_results(results_object = DEC_J_BSM_map_i,
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
                           skiptree = FALSE,
                           show.tip.label = TRUE)
  # Paint on the branch states
  paint_stochastic_map_branches(res = DEC_J_BSM_map_i,
                                master_table_cladogenetic_events = DEC_J_clado_events_master_table_i,
                                colors_list_for_states = colors_list_for_states,
                                lwd = 5, lty = par("lty"),
                                root.edge = TRUE, stratified = TRUE)
  # Replot the nodes states
  plot_BioGeoBEARS_results(results_object = DEC_J_BSM_map_i,
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


### 3.2/ Loop through the maps and plot a PDF per map ####


for (i in 1:nb_maps)
  # for (i in 1:100)
{
  # i <- 3
  
  # Extract the master tables of events summary for map n°i
  DEC_J_clado_events_table_i <- DEC_J_clado_events_tables[[i]]
  DEC_J_ana_events_table_i <- DEC_J_ana_events_tables[[i]]
  
  ## Extend the master table for cladogenetic events to account for branches that are split across different time-strata
  cols_to_get <- names(DEC_J_clado_events_table_i[,-ncol(DEC_J_clado_events_table_i)])
  colnums <- match(cols_to_get, names(DEC_J_ana_events_table_i))
  ana_events_table_cols_to_add <- DEC_J_ana_events_table_i[,colnums]
  anagenetic_events_txt_below_node <- rep("none", nrow(ana_events_table_cols_to_add))
  ana_events_table_cols_to_add <- cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  # rows_to_get_TF <- ana_events_table_cols_to_add$node <= length(Ponerinae_phylogeny_1534t_short_names$tip.label)
  rows_to_get_TF <- ana_events_table_cols_to_add$node <= length(Ponerinae_MCC_phylogeny_1534t_short_names$tip.label)
  DEC_J_clado_events_master_table_i <- rbind(ana_events_table_cols_to_add[rows_to_get_TF,], DEC_J_clado_events_table_i)
  
  # View(DEC_J_clado_events_table_i)
  # View(DEC_J_clado_events_master_table_i)
  
  # Convert the BSM into a modified fit object
  DEC_J_BSM_map_i <- stochastic_map_states_into_res(res = DEC_J_fit,  # Model fit object
                                                    master_table_cladogenetic_events = DEC_J_clado_events_master_table_i,
                                                    stratified = TRUE)
  
  # pdf(file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/DEC_J_BSMs_map_",i,".pdf"),
  #     height = 150, width = 40)
  pdf(file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/DEC_J_BSMs_map_",i,".pdf"),
      height = 150, width = 40)
  
  # Set title
  main_i <- paste0("\n\nDEC + J - Stochastic Map #", i, "/", nb_maps)
  
  # Plot stochastic map with node states
  plot_BioGeoBEARS_results(results_object = DEC_J_BSM_map_i,
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
                           skiptree = FALSE,
                           show.tip.label = TRUE)
  # Paint on the branch states
  paint_stochastic_map_branches(res = DEC_J_BSM_map_i,
                                master_table_cladogenetic_events = DEC_J_clado_events_master_table_i,
                                colors_list_for_states = colors_list_for_states,
                                lwd = 5, lty = par("lty"),
                                root.edge = TRUE, stratified = TRUE)
  # Replot the nodes states
  plot_BioGeoBEARS_results(results_object = DEC_J_BSM_map_i,
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
  dev.off()
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Stochastic map n°", i, "/", nb_maps,"\n"))
  }
}

# ## Combine all PDF in a single one with N pages
# 
# # all_BSM_maps_path <- list.files(path = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/", pattern = "BSMs_map_", full.names = T)
# all_BSM_maps_path <- list.files(path = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/", pattern = "BSMs_map_", full.names = T)
# 
# nb_maps <- length(all_BSM_maps_path)
# 
# # Reorder in numerical rather than alphabetic order
# prefix <- str_remove(string = all_BSM_maps_path[1], pattern = "_\\d+.pdf")
# all_BSM_maps_path <- paste0(prefix,"_",1:nb_maps,".pdf")
# 
# # qpdf::pdf_combine(input = all_BSM_maps_path, output = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_",nb_maps,"_maps.pdf"))
# qpdf::pdf_combine(input = all_BSM_maps_path, output = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_",nb_maps,"_maps.pdf"))

## Do it by smaller batch to avoid saturating RAM

batch_size <- 100
start_i <- seq(from = 1, to = nb_maps, by = batch_size)

for (i in seq_along(start_i))
{
  subset_indices_i <- start_i[i]:(start_i[i] + batch_size - 1)
  
  all_BSM_maps_path_i <- all_BSM_maps_path[subset_indices_i]
  nb_maps_i <- length(all_BSM_maps_path_i)
  subset_label <- paste0(subset_indices_i[1],"_",subset_indices_i[nb_maps_i])
  
  # qpdf::pdf_combine(input = all_BSM_maps_path_i, output = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_",subset_label,"_maps.pdf"))
  qpdf::pdf_combine(input = all_BSM_maps_path_i, output = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_",subset_label,"_maps.pdf"))
  
  # Print progress
  cat(paste0(Sys.time(), " - Stochastic maps aggregated from n°", subset_indices_i[1], " to ", subset_indices_i[nb_maps_i],"\n"))
}

## Aggregate batches in a single PDF

# all_BSM_maps_batches_path <- list.files(path = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/", pattern = "All_DEC_J_BSMs_", full.names = T)
all_BSM_maps_batches_path <- list.files(path = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/", pattern = "All_DEC_J_BSMs_", full.names = T)

nb_batches <- length(all_BSM_maps_batches_path)

# # Reorder in numerical rather than alphabetic order
# prefix <- str_remove(string = all_BSM_maps_path[1], pattern = "_\\d+.pdf")
# all_BSM_maps_path <- paste0(prefix,"_",1:nb_maps,".pdf")

# qpdf::pdf_combine(input = all_BSM_maps_batches_path, output = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps,"_maps.pdf"))
qpdf::pdf_combine(input = all_BSM_maps_batches_path, output = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps,"_maps.pdf"))



##### 4/ Produce animated Biogeographic Stochastic Maps ####

# See magickR... to use the PDF as input and create a GIF

source("./functions/image_resize_and_write_gif.R")

# Limit number of maps to reduce GIF size and avoid crashing because of memory allocation limits
nb_maps_in_GIF <- nb_maps
nb_maps_in_GIF <- 100

# Select fps
# fps <- 1
fps <- 5
# fps <- 10

# Load pdf as image frames
# pdf_pointer_BSM_maps <- magick::image_read_pdf(path = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps_in_GIF,"_maps.pdf"),
#                                                pages = NULL, density = 75)
pdf_pointer_BSM_maps <- magick::image_read_pdf(path = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps_in_GIF,"_maps.pdf"),
                                                       pages = NULL, density = 75)
magick::image_info(pdf_pointer_BSM_maps)

image_resize_and_write_gif(image = pdf_pointer_BSM_maps,
                           # path =  paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps_in_GIF,"_maps.gif"),
                           path =  paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/BSM_maps/All_DEC_J_BSMs_1_",nb_maps_in_GIF,"_maps.gif"),
                           delay = 1/fps, # Time between frames in seconds
                           width = 300, height = 1250,
                           loop = FALSE,
                           progress = TRUE)


# If does not work because of memory allocation issue. Use online tools instead...


##### 5/ Convert BSM outputs in phytools.simmaps ####

source("./functions/BSM_to_phytools_SM_custom.R")
source("./functions/generate_list_ranges.R")

## Reorder phylogeny in "cladewise" order so edge ID match what is expected in BSM_to_phytools_SM
# str(Ponerinae_phylogeny_1534t_short_names)
# Ponerinae_phylogeny_1534t_cladewise <- reorder.phylo(x = Ponerinae_phylogeny_1534t_short_names, order = "cladewise")
str(Ponerinae_MCC_phylogeny_1534t_short_names)
Ponerinae_MCC_phylogeny_1534t_cladewise <- reorder.phylo(x = Ponerinae_MCC_phylogeny_1534t_short_names, order = "cladewise")

# Initiate list for simmap
DEC_J_simmaps <- list()

## Loop across all BSM maps 
for (i in 1:length(DEC_J_clado_events_tables))
{
  # i <- 1
  
  # Convert BSM output to a SIMMAP object for phytools while computing residence times
  DEC_J_simmap_times_i <- BSM_to_phytools_SM_custom(res = DEC_J_fit,
                                                    # tree = Ponerinae_phylogeny_1534t_cladewise,
                                                    tree = Ponerinae_MCC_phylogeny_1534t_cladewise,
                                                    clado_events_table = DEC_J_clado_events_tables[[i]],
                                                    ana_events_table = DEC_J_ana_events_tables[[i]])
  
  # Store simmap
  DEC_J_simmaps[[i]] <- DEC_J_simmap_times_i$simmap
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    # Save simmaps
    # saveRDS(object = DEC_J_simmaps, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps.rds")
    saveRDS(object = DEC_J_simmaps, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps.rds")
    
    cat(paste0(Sys.time(), " - BSM map output converted to simmap Stochastic map - n°", i, "/", length(DEC_J_clado_events_tables),"\n"))
  }
}
class(DEC_J_simmaps) <- c("list", "multiSimmap")

# Save simmaps
# saveRDS(object = DEC_J_simmaps, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps.rds")
saveRDS(object = DEC_J_simmaps, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps.rds")



##### 6/ Plot continuous maps of Bioregion membership probability ####

# Compute continuous probabilities of each branch to belong to a focal bioregion
# Probabilities computed across the 1000 BS maps.

?phytools::densityMap

# Load simmaps of BS maps
# DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps.rds")
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps.rds")

# Subset for testing
# DEC_J_simmaps <- DEC_J_simmaps[1:10]

class(DEC_J_simmaps) <- c("list", "multiSimmap", "multiPhylo")

### 6.1/ Extract states and areas list ####

# Extract all states
all_states <- unique(unlist(lapply(X = DEC_J_simmaps, FUN = function (x) { lapply(X = x$maps, FUN = names) })))

# Ignore ranges encompassing multiple areas
all_areas <- all_states[nchar(all_states) == 1]
all_areas

# all_areas_labels <- c("Neotropics", "Afrotropics", "Eastern Palearctic", "Indomalaya", "Western Palearctic", "Australasia", "Nearctic")
all_areas_labels <- c("Afrotropics", "Neotropics", "Indomalaya", "Australasia", "Eastern Palearctic",  "Western Palearctic",  "Nearctic")

### 6.2/ Modify simmaps to attribute randomly a unique area to multi-states ####

# Better approach would be to use equal weights to split multi-states rather than random attribution of a unique area

# Function to attribute randomly a unique area to multi-states
attribute_unique_areas_to_simmap <- function (simmap, all_areas = NULL, verbose = F, seed = 1234) 
{
  # Set seed
  if (!is.null(seed))
  {
    set.seed(seed = seed)
  }
  
  # Extract all states
  all_states <- unique(unlist(lapply(X = simmap$maps, FUN = names)))
  
  # If not provided, extract unique areas from states
  if (is.null(all_areas))
  {
    all_areas <- all_states[nchar(all_states) == 1]
  }
  
  # Loop per edge
  for (i in seq_along(simmap$maps))
  {
    # i <- 2536
    map_i <- simmap$maps[[i]]
    
    # Identify multi-states
    nb_areas_per_mapped_states_i <- nchar(names(map_i))
    
    if (any(nb_areas_per_mapped_states_i > 1))
    {
      multi_state_maps_i <- map_i[nb_areas_per_mapped_states_i > 1]
      
      # Loop per mapped multi-states
      for (j in 1:length(multi_state_maps_i))
      {
        # j <- 2
        
        # Extract multi-state
        state_j <- names(multi_state_maps_i)[j]
        # Extract unique ranges
        unique_areas_j <- unlist(str_split(string = state_j, pattern = ""))
        
        # Draw randomly a single area
        new_area <- sample(x = unique_areas_j, size = 1)
        # Replace initial multi-state with unique area
        names(multi_state_maps_i)[j] <- new_area
      }
      # Patch new states
      names(map_i)[nb_areas_per_mapped_states_i > 1] <- names(multi_state_maps_i)
      simmap$maps[[i]] <- map_i
      
      ## May wish to merge consecutive similar states
    }
    # Print progress every 1000 edges
    if ((i %% 1000 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Unique area attributed to edge n°", i, "/", length(simmap$maps),"\n"))
    }
  }
  # Return updated simmap
  return(simmap)
}

set.seed(seed = 1234) # Set seed outside of the lapply
DEC_J_simmaps_unique_areas <- lapply(X = DEC_J_simmaps, FUN = attribute_unique_areas_to_simmap, all_areas = all_areas)

# Check that all states are unique areas
all_states <- unique(unlist(lapply(X = DEC_J_simmaps_unique_areas, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
all_states

# Save simmaps with only unique areas
# saveRDS(object = DEC_J_simmaps_unique_areas, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")
saveRDS(object = DEC_J_simmaps_unique_areas, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")


### 6.3/ Compute and plot density maps of binary areas ####

# Load simmaps with only unique areas
# DEC_J_simmaps_unique_areas <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")
DEC_J_simmaps_unique_areas <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")

## Loop per unique area
DEC_J_density_map_all_areas <- list()
for (i in seq_along(all_areas))
{
  # i <- 1
  
  # Extract the focal area
  focal_area <- all_areas[i]
  focal_area_label <- all_areas_labels[i]
  
  # Define other states by contrast
  other_states <- all_states[all_states != focal_area]
  
  ## 6.3.1/ Label simmaps with binary states ####
  
  DEC_J_simmaps_binary <- DEC_J_simmaps_unique_areas
  DEC_J_simmaps_binary <- lapply(X = DEC_J_simmaps_binary, FUN = phytools::mergeMappedStates, old.states = other_states, new.state = "0")
  DEC_J_simmaps_binary <- lapply(X = DEC_J_simmaps_binary, FUN = phytools::mergeMappedStates, old.states = focal_area, new.state = "1")
  
  class(DEC_J_simmaps_binary) <- c("list", "multiSimmap", "multiPhylo")
  
  # Check that remaining states are all binary
  # unique(unlist(lapply(X = DEC_J_simmaps_binary, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
  # DEC_J_simmaps_binary[[1]]$maps
  
  ## 6.3.2/ Estimate the posterior probabilities of states along all branches (from the set of simulated maps) ####
  
  # Use a custom version to have better control
  source(file = "./functions/densityMap_custom.R")
  
  # # Create a "densityMap" object
  # DEC_J_densityMap <- phytools::densityMap(trees = DEC_J_simmaps_binary,
  #                                          plot = TRUE)
  
  # Create a "densityMap" object
  DEC_J_density_map <- densityMap_custom(trees = DEC_J_simmaps_binary,
                                        tol = 1e-5, verbose = T,
                                        col_scale = NULL,
                                        plot = FALSE)
  
  ## 6.3.3/ Update color gradient ####
  
  ##  Create a custom color scale
  # Use color scheme of BSM
  colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
  
  # Set colors for areas/bioregions
  focal_color <- colors_list_for_states[focal_area]
  col_fn <- colorRampPalette(colors = c("grey90", focal_color))
  col_scale <- col_fn(n = 1001)
  
  # Update color gradient
  DEC_J_density_map <- setMap(DEC_J_density_map, c("grey90", focal_color))
  
  ## Save the resulting Density map with updated color gradient
  # saveRDS(object = DEC_J_density_map, file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_",focal_area_label,".rds"))
  saveRDS(object = DEC_J_density_map, file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_",focal_area_label,".rds"))
  
  ## Store in final object with all density maps
  DEC_J_density_map_all_areas <- append(x = DEC_J_density_map_all_areas, values = list(DEC_J_density_map))
  
  ## Print progress for each bioregion
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Posterior probability computed for Bioregion = ",focal_area_label," - n°", i, "/", length(all_areas),"\n"))
  }
}
names(DEC_J_density_map_all_areas) <- paste0("Density_map_", all_areas_labels)

# Save density maps
# saveRDS(object = DEC_J_density_map_all_areas, file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_all_areas.rds"))
saveRDS(object = DEC_J_density_map_all_areas, file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_all_areas.rds"))


### 6.4/ Plot density maps of posterior probability to belong to focal area ####

# Load density maps
# DEC_J_density_map_all_areas <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_all_areas.rds"))
DEC_J_density_map_all_areas <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_all_areas.rds"))


## Loop per area
for (i in 1:length(DEC_J_density_map_all_areas))
{
  # Extract the focal area label from the name
  focal_area_label <- names(DEC_J_density_map_all_areas)[i]
  focal_area_label <- str_remove(string = focal_area_label, pattern = "Density_map_")
  
  # Extract density map
  DEC_J_density_map_i <- DEC_J_density_map_all_areas[[i]]
  
  # Extract focal color
  focal_color <- DEC_J_density_map_i$cols[length(DEC_J_density_map_i$cols)]
  
  ## Plot density map
  
  # pdf(file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_",focal_area_label,".pdf"), width = 40, height = 150)
  pdf(file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_",focal_area_label,".pdf"), width = 40, height = 150)
  
  plot(DEC_J_density_map_i, fsize = c(0.3,2.5), lwd = c(3,4))
  
  # Extract plot metadata
  plot_phylo_metadata <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x_mean <- mean(plot_phylo_metadata$x.lim)
  
  # Add title
  text(x = x_mean, y = plot_phylo_metadata$y.lim[2]*1.02, adj = 0.5,
       labels = paste0("Posterior probability of range including ", focal_area_label),
       bty = "n", col = "black", font = 2, cex = 6)
  # Add legend title
  text(x = 29, y = plot_phylo_metadata$y.lim[2]*-0.068, adj = 0.5,
       labels = paste0("PP of being in ", focal_area_label),
       bty = "n", col = focal_color, font = 2, cex = 4)
  
  dev.off()
  
  # Output looks continuous, but in practice it is the summary of discrete histories of trait evolution

  ## Print progress for each bioregion
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Posterior probability plotted for Bioregion = ",focal_area_label," - n°", i, "/", length(DEC_J_density_map_all_areas),"\n"))
  }
}

### 6.5/ Merge density map plots in a unique PDF ####

# Keep Bioregions in similar order than in the BioGeoBEARS model
returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
returned_mats$areanames

bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Recreate paths to density maps in predefined order
# all_density_maps_paths <- paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_",bioregion_names,".pdf")
all_density_maps_paths <- paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_",bioregion_names,".pdf")

# Combine PDFs in a unique PDF
# qpdf::pdf_combine(input = all_density_maps_paths, output = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/All_DEC_J_density_maps.pdf"))
qpdf::pdf_combine(input = all_density_maps_paths, output = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/All_DEC_J_density_maps.pdf"))


