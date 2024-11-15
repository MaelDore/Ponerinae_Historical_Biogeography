##### Script 08: Discrete events per time-strata #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Aggregate discrete event counts per time strata
   # Compute number/probabilities of each class of events (mean and sd), per time-strata
   # Plot the histogram of number of each class of events, per time-strata
   # Plot the stacked bars of event types, per time-strata (raw and percentage)

# Plot Emigration/Immigration counts between bioregions, per time-strata

## Type of events
   # Per type of transitions
       # Anagenetic:
           # Range extension (d)
           # Range extinction (e)
       # Cladogenetic
           # Inheritance (y)
           # Vicariance (v)
           # Subset speciation (s)
           # Jump-dispersal (j)
   # Per bioregion
       # Immigration: Source of d/j
       # Emigration: Destination of d/j

###

### Inputs

# Best fitted Biogeographic model
# Biogeographic Stochastic Maps outputs

###

### Outputs

## Summary arrays
  # Tables of number/probabilities of each class of events per time stratum
  # Summary statistics (mean, sd, 95% HPD) of each class of events per time stratum

## Event type plots
  # Density curves of number of each class of events across maps
  # Stacked bars of mean number of event types (raw and percentage)
      # Overall vs. per time-strata vs. per bioregions (source/destination) vs. per time-strata X bioregions (source/destination)

## Dispersal networks
  # Emigration/Immigration counts between bioregions
      # Overall vs. per time-strata

###


####### May wish to focus on specific clades or bioregions to compute the same series of plots #######

# Overall
# Across bioregions, overall time
# Across time-strata, overall bioregions
# Across time-strata x bioregions

# All can be segregated per clades!


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(phytools)
library(geiger)
library(ape)
library(BioGeoBEARS)
library(MultinomialCI)    # For 95% CIs on BSM counts
library(BayesTwin)  # To compute HPD intervals
library(gridExtra)  # For Multifaceted ggplot
library(ggnewscale) # To add multiple scale of the same type in ggplot
library(sf)
library(igraph)     # For network plots

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load modeling results from best model (DEC+J) #####

# Load time-stratified DEC+J model output
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")

# Load BSM outputs
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_output.rds")
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_output.rds")

# Extract records of events
DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
DEC_J_ana_events_tables <- DEC_J_BSM_output$RES_ana_events_tables

# Extract areas list
returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
all_areas <- returned_mats$areanames

### 1.3/ Function for ggplot aesthetics ####

# Aesthetics for density curve plots
add_aesthetics_density_curve <- function (ggplot)
{
  ggplot <- ggplot + 
    
    theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
          panel.background = element_rect(fill = NA, color = NA),
          legend.title = element_text(size  = 20, margin = margin(b = 8)), 
          legend.text = element_text(size = 15),
          legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
          legend.key.size = unit(1.8, "line"),
          legend.spacing.y = unit(1.0, "line"),
          plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
          axis.title = element_text(size = 20, color = "black"),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 12)),
          axis.line = element_line(linewidth = 1.5),
          axis.text = element_text(size = 18, color = "black"),
          axis.text.x = element_text(margin = margin(t = 5)),
          axis.text.y = element_text(margin = margin(r = 5)))
  
  return(ggplot)
}



##### 2/ Summarize events from stochastic maps ####

?BioGeoBEARS::get_dmat_times_from_res # To generate matrix of range expansion from source area to destination area
?BioGeoBEARS::simulate_source_areas_ana_clado # Use to select randomly a unique area source for transition from a multi-area state to a single area
?BioGeoBEARS::count_ana_clado_events # To count the number of events per type of events from BSM tables
?BioGeoBEARS::hist_event_counts # To plot histograms of event counts across BSM tables

### 2.1/ Function to reorganize event counts lists as arrays for clarity ####

## Reorganize as arrays for clarity

convert_BSM_counts_lists_to_arrays <- function (x)
{
  # Initiate output
  out <- list()
  
  # Extract area list
  areas_list <- colnames(x$all_dispersals_counts_fromto_means)
  
  # Extract number of stochastic maps
  nb_maps <- dim(x$a_counts_cube)[3]
  
  # List type of dispersal events
  dispersal_event_types <- c("range-switching (a)", "range extension (d)", "jump-dispersal (j)", "anagenetic dispersal (a + d)", "all dispersal (a + d + j)")
  # List all types of events
  all_event_types <- c("range-switching (a)", "range extension (d)", "range contraction (e)", "subset speciation (s)", "vicariance (v)", "range inheritance (y)", "jump-dispersal (j)", "anagenetic dispersal (a + d)", "all anagenetic (a + d + e)", "all cladogenetic (s + v + y + j)", "all dispersal (a + d + j)", "all events (a + d + e + s + v + y + j)")
  
  # Unique events stays as df
  out$unique_sub_counts <- x$unique_sub_counts # Counts of individual subset speciation (s) events per stochastic maps. Ex: AB -> (B),(AB)
  out$unique_vic_counts <- x$unique_vic_counts # Counts of individual vicariance (v) events per stochastic maps. Ex: AB -> (A),(B)
  out$unique_sym_counts <- x$unique_sym_counts # Counts of individual range-inheritance / narrow sympatry (y) events per stochastic maps. Ex: A -> (A),(A)
  
  # Source/dest dispersal count matrices are aggregated in an array across type of events
  out$dispersal_events_count_matrices <- array(data = NA,
                                               dim = c(length(areas_list), length(areas_list), length(dispersal_event_types), nb_maps),
                                               dimnames = list(areas_list, areas_list, dispersal_event_types, paste0("Map_", 1:nb_maps)))
  out$dispersal_events_count_matrices[ , , "range-switching (a)", ] <- x$a_counts_cube
  out$dispersal_events_count_matrices[ , , "range extension (d)", ] <- x$d_counts_cube
  out$dispersal_events_count_matrices[ , , "jump-dispersal (j)", ] <- x$founder_counts_cube
  out$dispersal_events_count_matrices[ , , "anagenetic dispersal (a + d)", ] <- x$anagenetic_dispersals_counts_cube
  out$dispersal_events_count_matrices[ , , "all dispersal (a + d + j)", ] <- x$all_dispersals_counts_cube
  
  # Source/dest summary matrices are aggregated in an array across type of events and stats
  out$dispersal_events_summary_matrices <- array(data = NA,
                                                 dim = c(length(areas_list), length(areas_list), length(dispersal_event_types), 5),
                                                 dimnames = list(areas_list, areas_list, dispersal_event_types, c("mean", "median", "sd", "2.5% HPD", "97.5% HPD")))
  out$dispersal_events_summary_matrices[ , , , "mean"] <- apply(X = out$dispersal_events_count_matrices, MARGIN = c(1,2,3), FUN = mean)
  out$dispersal_events_summary_matrices[ , , , "median"] <- apply(X = out$dispersal_events_count_matrices, MARGIN = c(1,2,3), FUN = median)
  out$dispersal_events_summary_matrices[ , , , "sd"] <- apply(X = out$dispersal_events_count_matrices, MARGIN = c(1,2,3), FUN = sd)
  HPD95 <- apply(X = out$dispersal_events_count_matrices, MARGIN = c(1,2,3), FUN = BayesTwin::HPD, cred_int = 0.95)
  out$dispersal_events_summary_matrices[ , , , c("2.5% HPD", "97.5% HPD")] <- aperm(HPD95, c(2,3,4,1))
  
  # Range extinction counts x states stays as df
  out$range_contraction_events_count <- x$e_counts_rectangle # Anagenetic range contraction (e) events. 
  colnames(out$range_contraction_events_count) <- areas_list
  row.names(out$range_contraction_events_count) <- paste0("Map_", 1:nb_maps)
  
  # Counts of events per maps are aggregated across all type of events
  out$all_events_count_per_maps <- array(data = NA,
                                         dim = c(length(all_event_types), nb_maps),
                                         dimnames = list(all_event_types, paste0("Map_", 1:nb_maps)))
  out$all_events_count_per_maps["range-switching (a)", ] <- x$a_totals_list 
  out$all_events_count_per_maps["range extension (d)", ] <- x$d_totals_list 
  out$all_events_count_per_maps["range contraction (e)", ] <- x$e_totals_list
  out$all_events_count_per_maps["subset speciation (s)", ] <- x$subsetSymp_totals_list
  out$all_events_count_per_maps["vicariance (v)", ] <- x$vicariance_totals_list
  out$all_events_count_per_maps["range inheritance (y)", ] <- x$clado_totals_list - (x$subsetSymp_totals_list + x$vicariance_totals_list + x$founder_totals_list)
  out$all_events_count_per_maps["jump-dispersal (j)", ] <- x$founder_totals_list
  out$all_events_count_per_maps["anagenetic dispersal (a + d)", ] <- x$anagenetic_dispersals_totals_list
  out$all_events_count_per_maps["all anagenetic (a + d + e)", ] <- x$ana_totals_list
  out$all_events_count_per_maps["all cladogenetic (s + v + y + j)", ] <- x$clado_totals_list
  out$all_events_count_per_maps["all dispersal (a + d + j)", ] <- x$all_dispersals_totals_list
  out$all_events_count_per_maps["all events (a + d + e + s + v + y + j)", ] <- x$all_totals_list
  
  # Recreate the summary table per event types with more stats
  out$all_events_summary_df <- array(data = NA,
                                     dim = c(length(all_event_types), 5),
                                     dimnames = list(all_event_types, c("mean", "median", "sd", "2.5% HPD", "97.5% HPD")))
  out$all_events_summary_df[ , "mean"] <- apply(X = out$all_events_count_per_maps, MARGIN = 1, FUN = mean)
  out$all_events_summary_df[ , "median"] <- apply(X = out$all_events_count_per_maps, MARGIN = 1, FUN = median)
  out$all_events_summary_df[ , "sd"] <- apply(X = out$all_events_count_per_maps, MARGIN = 1, FUN = sd)
  HPD95 <- apply(X = out$all_events_count_per_maps, MARGIN = 1, FUN = BayesTwin::HPD, cred_int = 0.95)
  out$all_events_summary_df[ , c("2.5% HPD", "97.5% HPD")] <- t(HPD95)
  
  # Return final list
  return(out)
}


### 2.2/ Summarize events on total time-frame ####

# Get the dmat (?) and time-strata boundaries (if any)
DEC_J_dmat_times <- get_dmat_times_from_res(res = DEC_J_fit, numstates = NULL)
DEC_J_dmat_times

# Simulate the areas chosen as the source for cases of dispersal events where the ancestral area occupied 2 or more areas
# Source areas are chosen probabilistically according to the dispersal constrains applied on transition (i.e., multipliers, distances, ...)
# Used to compute the number of events across areas using only single-range areas
DEC_J_BSM_sourceAreas <- simulate_source_areas_ana_clado(res = DEC_J_fit,
                                                         clado_events_tables = DEC_J_clado_events_tables, 
                                                         ana_events_tables = DEC_J_ana_events_tables,
                                                         areanames = all_areas)

DEC_J_clado_events_tables_with_unique_source <- DEC_J_BSM_sourceAreas$clado_events_tables
DEC_J_ana_events_tables_with_unique_source <- DEC_J_BSM_sourceAreas$ana_events_tables

# Save clado/ana events tables with unique source
# saveRDS(object = DEC_J_clado_events_tables_with_unique_source, "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/Ponerinae_MCC_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
# saveRDS(object = DEC_J_ana_events_tables_with_unique_source, "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")
saveRDS(object = DEC_J_clado_events_tables_with_unique_source, "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
saveRDS(object = DEC_J_ana_events_tables_with_unique_source, "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")

# Count all anagenetic and cladogenetic events
DEC_J_BSM_counts_list <- count_ana_clado_events(clado_events_tables = DEC_J_clado_events_tables_with_unique_source,
                                                ana_events_tables = DEC_J_ana_events_tables_with_unique_source,
                                                areanames = all_areas,
                                                actual_names = all_areas)
# Contains all information on counts of all possible cladogenetic and anagenetic events
# Ex: Number of transitions from A to C
# Contains also summary of events per class (i.e., range extension, vicariance, ...)
str(DEC_J_BSM_counts_list)
names(DEC_J_BSM_counts_list)

# Save summary lists of event counts
# saveRDS(object = DEC_J_BSM_counts_list, "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_list.rds")
saveRDS(object = DEC_J_BSM_counts_list, "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_list.rds")

# Summarize counts per type of events
DEC_J_summary_counts_BSMs <- DEC_J_BSM_counts_list$summary_counts_BSMs
print(conditional_format_table(DEC_J_summary_counts_BSMs))

# Get posterior probabilities of class of events at cladogenesis
DEC_J_summary_counts_BSMs["means", c("subset", "founder", "vicariance", "sympatry")] / DEC_J_summary_counts_BSMs["means", c("all_clado")]
# Different from the relative weights parameters!

# Plot histogram of event counts (generate PDF directly)
hist_event_counts(counts_list = DEC_J_BSM_counts_list, 
                  # pdffn = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_Histo_event_counts.pdf"))
                  pdffn = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_Histo_event_counts.pdf"))

# Convert to arrays for clarity
DEC_J_BSM_counts_arrays <- convert_BSM_counts_lists_to_arrays(DEC_J_BSM_counts_list)
str(DEC_J_BSM_counts_arrays, max.level = 2)

DEC_J_BSM_counts_arrays$all_events_summary_df

# Save summary arrays of event counts
# saveRDS(object = DEC_J_BSM_counts_arrays, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays.rds")
saveRDS(object = DEC_J_BSM_counts_arrays, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays.rds")

# Load summary arrays of event counts
# DEC_J_BSM_counts_arrays <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays.rds")
DEC_J_BSM_counts_arrays <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays.rds")


### 2.3/ Summary events per time stratum ####

## Can extract count for any time-frame, not the just the time-periods defined for the model !!!

## Make a loop to count events for all time-strata

# Extract the nb of time-periods
nb_strata <- length(DEC_J_fit$inputs$timeperiods)

# Extract the nb of maps
nb_maps <- length(DEC_J_clado_events_tables_with_unique_source)

# Initiate list of final arrays
DEC_J_BSM_counts_arrays_all_strata <- list()

# Loop per time-strata
for (i in 1:nb_strata)
{
  # i <- 2
  
  # Define time period
  if (i == 1)
  {
    start_time <- 0
  } else {
    start_time <- DEC_J_fit$inputs$timeperiods[i-1]
  }
  end_time <- DEC_J_fit$inputs$timeperiods[i]
  
  # Count all anagenetic and cladogenetic events for one time-stratum
  DEC_J_BSM_counts_list_per_stratum_i <- count_ana_clado_events(clado_events_tables = DEC_J_clado_events_tables_with_unique_source,
                                                                ana_events_tables = DEC_J_ana_events_tables_with_unique_source,
                                                                areanames = all_areas,
                                                                actual_names = all_areas,
                                                                timeperiod = c(start_time, end_time))
  
  str(DEC_J_BSM_counts_list_per_stratum_i, max.level = 1)
  
  # Convert in array
  DEC_J_BSM_counts_arrays_per_stratum_i <- convert_BSM_counts_lists_to_arrays(DEC_J_BSM_counts_list_per_stratum_i)
  
  ## Initialize final lists/arrays once
  if (i == 1)
  {
    # Unique subset events counts
    unique_sub_counts_all_strata <- list()
    # Unique vicariance events counts
    unique_vic_counts_all_strata <- list()
    # Unique range inheritance (sympatry) events counts
    unique_sym_counts_all_strata <- list()
    # Source/dest dispersal events count matrices
    dispersal_events_count_matrices_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$dispersal_events_count_matrices
    dispersal_events_count_matrices_all_strata <- array(data = NA,
                                                        dim = c(dim(dispersal_events_count_matrices_stratum_i), nb_strata),
                                                        dimnames = append(x = dimnames(dispersal_events_count_matrices_stratum_i), values = list(paste0("Stratum_", 1:nb_strata))))
    # Source/dest dispersal events summary matrices
    dispersal_events_summary_matrices_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$dispersal_events_summary_matrices
    dispersal_events_summary_matrices_all_strata <- array(data = NA,
                                                          dim = c(dim(dispersal_events_summary_matrices_stratum_i), nb_strata),
                                                          dimnames = append(x = dimnames(dispersal_events_summary_matrices_stratum_i), values = list(paste0("Stratum_", 1:nb_strata))))
    # Range extinction counts x states
    range_contraction_events_count_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$range_contraction_events_count
    range_contraction_events_count_all_strata <- array(data = NA,
                                                       dim = c(dim(range_contraction_events_count_stratum_i), nb_strata),
                                                       dimnames = append(x = dimnames(range_contraction_events_count_stratum_i), values = list(paste0("Stratum_", 1:nb_strata))))
    # Counts of events per maps
    all_events_count_per_maps_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$all_events_count_per_maps
    all_events_count_per_maps_all_strata <- array(data = NA,
                                                  dim = c(dim(all_events_count_per_maps_stratum_i), nb_strata),
                                                  dimnames = append(x = dimnames(all_events_count_per_maps_stratum_i), values = list(paste0("Stratum_", 1:nb_strata))))
    # Summary tables of all event types
    all_events_summary_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$all_events_summary_df
    all_events_summary_all_strata <- array(data = NA,
                                           dim = c(dim(all_events_summary_stratum_i), nb_strata),
                                           dimnames = append(x = dimnames(all_events_summary_stratum_i), values = list(paste0("Stratum_", 1:nb_strata))))
  }
  
  ## Save outputs in objects used to aggregate information across all strata
  
  # Unique subset events counts
  unique_sub_counts_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$unique_sub_counts
  old_names <- names(unique_sub_counts_all_strata)
  unique_sub_counts_all_strata <- append(x = DEC_J_BSM_counts_arrays_all_strata$unique_sub_counts_all_strata, values = list(unique_sub_counts_stratum_i))
  names(unique_sub_counts_all_strata) <- c(old_names, paste0("Stratum_", i))
  DEC_J_BSM_counts_arrays_all_strata$unique_sub_counts_all_strata <- unique_sub_counts_all_strata
  
  # Unique vicariance events counts
  unique_vic_counts_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$unique_vic_counts
  old_names <- names(DEC_J_BSM_counts_arrays_all_strata$unique_vic_counts_all_strata)
  DEC_J_BSM_counts_arrays_all_strata$unique_vic_counts_all_strata <- append(x = DEC_J_BSM_counts_arrays_all_strata$unique_vic_counts_all_strata, values = list(unique_vic_counts_stratum_i))
  names(DEC_J_BSM_counts_arrays_all_strata$unique_vic_counts_all_strata) <- c(old_names, paste0("Stratum_", i))
  
  # Unique range inheritance (sympatry) events counts
  unique_sym_counts_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$unique_sym_counts
  old_names <- names(DEC_J_BSM_counts_arrays_all_strata$unique_sym_counts_all_strata)
  DEC_J_BSM_counts_arrays_all_strata$unique_sym_counts_all_strata <- append(x = DEC_J_BSM_counts_arrays_all_strata$unique_sym_counts_all_strata, values = list(unique_sym_counts_stratum_i))
  names(DEC_J_BSM_counts_arrays_all_strata$unique_sym_counts_all_strata) <- c(old_names, paste0("Stratum_", i))
  
  # Source/dest dispersal events count matrices
  dispersal_events_count_matrices_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$dispersal_events_count_matrices
  dispersal_events_count_matrices_all_strata[ , , , ,i] <- dispersal_events_count_matrices_stratum_i
  DEC_J_BSM_counts_arrays_all_strata$dispersal_events_count_matrices_all_strata <- dispersal_events_count_matrices_all_strata
  
  # Source/dest dispersal events summary matrices
  dispersal_events_summary_matrices_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$dispersal_events_summary_matrices
  dispersal_events_summary_matrices_all_strata[ , , , ,i] <- dispersal_events_summary_matrices_stratum_i
  DEC_J_BSM_counts_arrays_all_strata$dispersal_events_summary_matrices_all_strata <- dispersal_events_summary_matrices_all_strata
  
  # Range extinction counts x states
  range_contraction_events_count_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$range_contraction_events_count
  range_contraction_events_count_all_strata[ , ,i] <- as.matrix(range_contraction_events_count_stratum_i)
  DEC_J_BSM_counts_arrays_all_strata$range_contraction_events_count_all_strata <- range_contraction_events_count_all_strata
  
  # Counts of events per maps
  all_events_count_per_maps_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$all_events_count_per_maps
  all_events_count_per_maps_all_strata[ , ,i] <- as.matrix(all_events_count_per_maps_stratum_i)
  DEC_J_BSM_counts_arrays_all_strata$all_events_count_per_maps_all_strata <- all_events_count_per_maps_all_strata
  
  # Summary tables of all event types
  all_events_summary_stratum_i <- DEC_J_BSM_counts_arrays_per_stratum_i$all_events_summary_df
  all_events_summary_all_strata[ , ,i] <- as.matrix(all_events_summary_stratum_i)
  DEC_J_BSM_counts_arrays_all_strata$all_events_summary_all_strata <- all_events_summary_all_strata
}

str(DEC_J_BSM_counts_arrays_all_strata)


# Save summary arrays of event counts across all time-strata
# saveRDS(object = DEC_J_BSM_counts_arrays_all_strata, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_strata.rds")
saveRDS(object = DEC_J_BSM_counts_arrays_all_strata, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_strata.rds")

# Load summary arrays of event counts across all time-strata
# DEC_J_BSM_counts_arrays_all_strata <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_strata.rds")
DEC_J_BSM_counts_arrays_all_strata <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_strata.rds")



##### 3/ Plot global summary of events counts #####

## List event types for plots

event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")
# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# Extract counts only for targeted type of events
all_events_count_per_maps <- t(DEC_J_BSM_counts_arrays$all_events_count_per_maps[event_types_list, ])

## Format data for ggplot
all_events_count_per_maps_ggplot <- all_events_count_per_maps %>%
  as.data.frame() %>% 
  mutate(Map = row.names(.)) %>%
  pivot_longer(cols = contains(" ("), names_to = "event_type", values_to = "counts")

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_ggplot$counts)

# Color scheme for events
event_types_col <- setNames(object = c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid"), nm = event_types_list)
event_types_labels <- str_to_sentence(names(event_types_col))

### 3.1/ Unique density curves per type of events plot ####

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_overall_single.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_overall_single.pdf", width = 10, height = 6)

density_curve_ggplot <- ggplot(data = all_events_count_per_maps_ggplot,
                               mapping = aes(x = counts, fill = event_type)) +
  # Plot density curve
  geom_density(alpha = 0.8) +
  
  # Set x limits to all be the same across event types
  xlim(0, max_counts) +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Counts of events") +
  ylab("Density") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", breaks = event_types_list, labels = event_types_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
density_curve_ggplot <- density_curve_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(density_curve_ggplot)

dev.off()

### 3.2/ Multifaceted density curves per type of events plot: 1 row = 1 type of events ####

## Generate ggplot
density_curve_ggplot_list <- list()
for ( i in seq_along(event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- event_types_list[i]
  
  # Extract event count data
  event_count_per_maps_i <- all_events_count_per_maps_ggplot[all_events_count_per_maps_ggplot$event_type == event_type_i, ]
  
  # Extract associated color
  event_col_i <- event_types_col[i]
  
  # Create density curve plot
  density_curve_ggplot <- ggplot(data = event_count_per_maps_i,
                                 mapping = aes(x = counts)) +
    # Plot density curve
    geom_density(fill = event_col_i, alpha = 0.8) +
    
    # Set x limits to all be the same across event types
    xlim(0, max_counts) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_ggplot <- density_curve_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_ggplot)
  
  # Store ggplot in a list
  density_curve_ggplot_list <- append(x = density_curve_ggplot_list, values = list(density_curve_ggplot))
}

## Multifaceted plot: 1 row = 1 type of events

# ggpubr::ggarrange(density_curve_ggplot_list,        # List of plots
#           # labels = c("A", "B", "C", "D"),   # Add labels
#           ncol = 1, nrow = length(density_curve_ggplot_list)) # Define window division

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_overall_multifaceted.pdf", width = 6, height = length(density_curve_ggplot_list) * 4)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_overall_multifaceted.pdf", width = 6, height = length(density_curve_ggplot_list) * 4)

gridExtra::grid.arrange(
  grobs = density_curve_ggplot_list, # List of ggplots
  widths = 4,  # Width of columns
  heights = rep(1, length(density_curve_ggplot_list)),
  nrow = length(density_curve_ggplot_list),
  ncol = 1)

dev.off()

### 3.3/ Plot mean counts of events ####

## List event types for plots

event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")
# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# Extract counts only for targeted type of events
all_events_count_per_maps <- t(DEC_J_BSM_counts_arrays$all_events_count_per_maps[event_types_list, ])

## Format mean data for ggplot
all_events_mean_count_per_maps_ggplot <- all_events_count_per_maps_ggplot %>%
  group_by(event_type) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  ungroup()

# Color scheme for events
event_types_col <- c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid")
event_types_legend_labels <- str_to_sentence(event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))

# Adjust order of event types
all_events_mean_count_per_maps_ggplot$event_type <- factor(all_events_mean_count_per_maps_ggplot$event_type, levels = event_types_list, labels = event_types_axis_labels)

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_overall.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_overall.pdf", width = 10, height = 6)

barplot_mean_count_per_events_ggplot <- ggplot(data = all_events_mean_count_per_maps_ggplot,
                                               mapping = aes(y = mean_counts, x = event_type, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Event types") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", breaks = event_types_axis_labels, labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_mean_count_per_events_ggplot <- barplot_mean_count_per_events_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_ggplot)

dev.off()


##### 4/ Plot global summary of events counts per time-strata #####

## List event types for plots

event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")
# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# Extract counts only for targeted type of events
all_events_count_per_maps_all_strata <- DEC_J_BSM_counts_arrays_all_strata$all_events_count_per_maps_all_strata[event_types_list, ,]

## Format data for ggplot
all_events_count_per_maps_all_strata_ggplot <- all_events_count_per_maps_all_strata %>%
  reshape2::melt()
names(all_events_count_per_maps_all_strata_ggplot) <- c("event_type", "map", "stratum", "counts")

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_all_strata_ggplot$counts)

# Color scheme for time-strata
time_strata_list <- unique(all_events_count_per_maps_all_strata_ggplot$stratum)
time_strata_col <- setNames(object = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"), nm = time_strata_list)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")

# Color scheme for events
event_types_col <- c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid")
event_types_legend_labels <- str_to_sentence(event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))


### 4.1/ Plot density curves of events, per time-strata ####

# Multifaceted: 
   # Rows = type of events
   # Fill of the overlaying histogram/density curve = geological epochs


## Generate ggplot
density_curve_all_strata_ggplot_list <- list()
for ( i in seq_along(event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- event_types_list[i]
  
  # Extract event count data
  event_count_all_strata_per_maps_i <- all_events_count_per_maps_all_strata_ggplot[all_events_count_per_maps_all_strata_ggplot$event_type == event_type_i, ]
  
  # Create density curve plot
  density_curve_all_strata_ggplot <- ggplot(data = event_count_all_strata_per_maps_i,
                                 mapping = aes(x = counts, fill = stratum)) +
    # Plot density curve
    geom_density(alpha = 0.8, adjust = 3) + # Use larger bandwidth to smooth density curve

    # Set x limits to all be the same across event types
    xlim(0, max_counts) +
    
    # Use logarithmic scale
    scale_x_log10() +
    
    # Adjust color scheme and legend
    scale_fill_manual("Geological epochs", breaks = time_strata_list, labels = time_strata_labels, values = time_strata_col) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_all_strata_ggplot <- density_curve_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_all_strata_ggplot)
  
  # Store ggplot in a list
  density_curve_all_strata_ggplot_list <- append(x = density_curve_all_strata_ggplot_list, values = list(density_curve_all_strata_ggplot))
}

## Multifaceted plot: 1 row = 1 type of events

# ggpubr::ggarrange(density_curve_ggplot_list,        # List of plots
#           # labels = c("A", "B", "C", "D"),   # Add labels
#           ncol = 1, nrow = length(density_curve_ggplot_list)) # Define window division

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_strata.pdf", width = 10, height = length(density_curve_all_strata_ggplot_list) * 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_strata.pdf", width = 10, height = length(density_curve_all_strata_ggplot_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_strata_ggplot_list, # List of ggplots
  widths = 4,  # Width of columns
  heights = rep(1, length(density_curve_all_strata_ggplot_list)),
  nrow = length(density_curve_all_strata_ggplot_list),
  ncol = 1)

dev.off()

### 4.2/ Plot stacked bars of raw mean event counts, per time-strata ####

# Plot the stacked bars of event types, per time-strata (raw counts)
   # Fill = type of events
   # X = Time-strata

## Format mean data for ggplot
all_events_mean_count_per_maps_all_strata_ggplot <- all_events_count_per_maps_all_strata_ggplot %>%
  group_by(stratum, event_type) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Adjust order of event types
all_events_mean_count_per_maps_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_strata_ggplot$event_type, levels = event_types_list, labels = event_types_axis_labels)

# Adjust order of time-strata
time_strata_labels <- c("PPHo", "Miocene", "Oligo.", "Eocene", "Paleo.", "Late Cr.", "Early Cr.")
all_events_mean_count_per_maps_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))


# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_strata.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_strata.pdf", width = 10, height = 6)

barplot_mean_count_per_events_all_strata_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_strata_ggplot,
                                                          mapping = aes(y = mean_counts, x = stratum, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Geological epochs") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", breaks = event_types_axis_labels, labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
                                   # color = rev(time_strata_col),
                                   # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_count_per_events_all_strata_ggplot <- barplot_mean_count_per_events_all_strata_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_all_strata_ggplot)

dev.off()

### 4.3/ Plot stacked bars of percentages of mean event counts, per time-strata ####

# Plot the stacked bars of event types, per time-strata (percentages)
   # Fill = type of events
   # X = Time-strata

## Format mean data for ggplot
all_events_mean_count_per_maps_all_strata_ggplot <- all_events_count_per_maps_all_strata_ggplot %>%
  group_by(stratum, event_type) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Adjust order of event types
all_events_mean_count_per_maps_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_strata_ggplot$event_type, levels = event_types_list, labels = event_types_axis_labels)

# Adjust order of time-strata
time_strata_labels <- c("PPHo", "Miocene", "Oligo.", "Eocene", "Paleo.", "Late Cr.", "Early Cr.")
all_events_mean_count_per_maps_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))

# Plot
# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_strata.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_strata.pdf", width = 10, height = 6)

barplot_mean_perc_per_events_all_strata_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_strata_ggplot,
                                                          mapping = aes(y = mean_perc, x = stratum, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Geological epochs") +
  ylab("% of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", breaks = event_types_axis_labels, labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_perc_per_events_all_strata_ggplot <- barplot_mean_perc_per_events_all_strata_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_perc_per_events_all_strata_ggplot)

dev.off()


##### 5/ Prepare date for counts of events per source/dest bioregions #####

### 5.1/ Function to convert unique events counts in source/dest matrices ####

convert_unique_events_to_source_dest_matrices <- function (unique_events_list, areas_list = NULL, event_type)
{
  # Input = List of dataframes with unique events to convert in 3D array
    # 1 dataframe per strata
  
  ## Check that the event type is a valid one
  if (!(event_type %in% c("subset speciation (s)", "vicariance (v)", "range inheritance (y)")))
  {
    stop(paste0('event_type must be either "subset speciation (s)", "vicariance (v)", or "range inheritance (y)"'))
  }
  
  ## If not provided, extract list of areas from data
  if (is.null(areas_list))
  {
    all_events_char <- unlist(lapply(X = unique_events_list, FUN = names))
    all_events_char_split <- unique(unlist(str_split(string = all_events_char, pattern = "")))
    # Remove non-coding characters
    areas_list <- all_events_char_split[!(all_events_char_split %in% c("-", ">", ","))]
  }
  
  ## Initiate final output array
  nb_bioregions <- length(areas_list)
  nb_event_types <- 1
  nb_maps <- nrow(unique_events_list[[1]])
  nb_strata <- length(unique_events_list)
  
  all_unique_events_source_dest_array <- array(data = 0,
                                               dim = c(nb_bioregions, nb_bioregions, nb_event_types, nb_maps, nb_strata),
                                               dimnames = list(areas_list, areas_list, event_type, paste0("Map_",row.names(unique_events_list[[1]])), names(unique_events_list)))
  
  ## Identify source/dest areas from unique event titles
  
  # Depend on the event type: subset speciation, vicariance, inheritance
  
  # For subset speciation: Ex: (AB) -> (AB),(A)
    # Source: spread weight among all the single areas of the parental multi-areas state
    # Destination: spread weight among all the single areas of the smallest children state (the one differing from the parental one)
  
  # For vicariance: Ex: (AB) -> (A),(B)
     # Source: spread weight among all the single areas of the parental multi-areas state
     # Destination: spread weight among all the single areas of the children states (both may be multi-areas in the case of wide-spread vicariance in DIVALIKE models)
  
  # For inheritance/sympatric speciation: Ex: (AB) -> (AB),(AB)
     # Source: spread weight among all the single areas of the parental state (may be multi-areas in the case of widespread sympatry in BAYAREALIKE models)
     # Destination: spread weight among all the single areas of the children states = parental multi-areas state (may be multi-areas in the case of widespread sympatry in BAYAREALIKE models)
  
  ## Loop per time-strata
  for (i in 1:length(unique_events_list))
  {
    # i <- 1
    
    # Extract data.frame for Strata i
    events_df_i <- unique_events_list[[i]]
    
    if (ncol(events_df_i) > 0)
    {
      # Extract parental and children states
      parental_states <- str_remove(string = names(events_df_i), pattern = "->.*")
      child1_states <- str_remove(string = str_remove(string = names(events_df_i), pattern = ".*->"), pattern = ",.*")
      # child2_states <- str_remove(string = names(events_df_i), pattern = ".*,")
      
      # Extract unique areas of parental and children states
      parental_areas <- str_split(string = parental_states, pattern = "")
      child1_areas <- str_split(string = child1_states, pattern = "")
      # child2_areas <- str_split(string = child2_states, pattern = "")
      
      # Loop per unique events
      for (j in 1:ncol(events_df_i))
      {
        # j <- 3
        
        # Extract parental areas of event j
        parental_areas_j <- parental_areas[[j]]
        child1_areas_j <- child1_areas[[j]]
        
        ## Generate combination of source/dest unique areas to use to split counts
        if (event_type %in% c("vicariance (v)", "range inheritance (y)"))
        {
          # In case of v and y, counts are split among combination of parental unique areas (as children areas are the same as parental areas)
          areas_comb_j <- expand.grid(parental_areas_j, parental_areas_j)
        } else {
          # In case of s, counts are split among combination of parental unique areas and smallest child unique areas
          areas_comb_j <- expand.grid(parental_areas_j, child1_areas_j)
        }
        # Compute splitting factors as the number of combination
        nb_comb_j <- nrow(areas_comb_j)
        
        ## Increment counts for each source/dest combination
        for (k in 1:nb_comb_j)
        {
          # k <- 1
          source_k <- areas_comb_j[k,1]
          dest_k <- areas_comb_j[k,2]
          counts_k <- events_df_i[,j] / nb_comb_j
          all_unique_events_source_dest_array[as.character(source_k), as.character(dest_k), event_type, , i] <- all_unique_events_source_dest_array[as.character(source_k), as.character(dest_k), event_type, , i] + counts_k
        }
        # print(all_unique_events_source_dest_array[ , as.character(dest_k) , event_type, , i])
      }
    }
    
    ## Print progress
    if (i %% 1 == 0)
    {
      cat(paste0(Sys.time(), " - Unique ",event_type," event counts recorded for Time-stratum n°", i, "/", length(unique_events_list),"\n"))
    }
  }
  
  ## Export output
  
  # Output: Array
    # Dim 1 = Source
    # Dim 2 = Dest
    # Dim 3 = Event_type
    # Dim 4 = Map
    # Dim 5 = Stratum
  
  return(all_unique_events_source_dest_array)
}



### 5.2/ Build array of events per source/dest areas ####


## 5.2.1/ Initiate final array

# Extract dimensions
areas_list <- dimnames(DEC_J_BSM_counts_arrays_all_strata$dispersal_events_count_matrices_all_strata)[[1]]
event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "range inheritance (y)", "jump-dispersal (j)")
maps_list <- dimnames(DEC_J_BSM_counts_arrays_all_strata$dispersal_events_count_matrices_all_strata)[[4]]
strata_list <- dimnames(DEC_J_BSM_counts_arrays_all_strata$dispersal_events_count_matrices_all_strata)[[5]]

DEC_J_BSM_all_unique_events_source_dest_array <- array(data = 0,
                                             dim = c(length(areas_list), length(areas_list), length(event_types_list), length(maps_list), length(strata_list)),
                                             dimnames = list(areas_list, areas_list, event_types_list, maps_list, strata_list))

## 5.2.2/ Get s, v, y from unique_events_df

## Get counts for subset speciation events
all_unique_sub_events_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_strata$unique_sub_counts_all_strata, areas_list = areas_list, event_type = "subset speciation (s)")

# all_unique_sub_events_source_dest_array[ , , , "Map_1" , "Stratum_1"]
# sum(all_unique_sub_events_source_dest_array[ , , , "Map_1" , "Stratum_1"])
# 
# DEC_J_BSM_counts_arrays_all_strata$unique_sub_counts_all_strata[[1]][1, ]
# sum(DEC_J_BSM_counts_arrays_all_strata$unique_sub_counts_all_strata[[1]][1, ])

## Counts need to be equal !

## Get counts for vicariance events
all_unique_vic_events_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_strata$unique_vic_counts_all_strata, areas_list = areas_list, event_type = "vicariance (v)")

## Should be symmetrical
# all_unique_vic_events_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Get counts for range inheritance events
all_unique_sym_events_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_strata$unique_sym_counts_all_strata, areas_list = areas_list, event_type = "range inheritance (y)")

## Should be only on the diagonal if no widespread sympatric speciation is allowed (in BAYAREALIKE models)
# all_unique_sym_events_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Fill the array with all event types

# Dim 1 = Source
# Dim 2 = Dest
# Dim 3 = Event_type
# Dim 4 = Map
# Dim 5 = Stratum

DEC_J_BSM_all_unique_events_source_dest_array[ , ,"subset speciation (s)", , ] <- all_unique_sub_events_source_dest_array
DEC_J_BSM_all_unique_events_source_dest_array[ , ,"vicariance (v)", , ] <- all_unique_vic_events_source_dest_array
DEC_J_BSM_all_unique_events_source_dest_array[ , ,"range inheritance (y)", , ] <- all_unique_sym_events_source_dest_array

## 5.2.3/ Get d and j from array

# Extract counts only for targeted type of events

dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")

# Extract counts only for targeted type of events
all_dispersal_events_count_matrices_per_maps_all_strata <- DEC_J_BSM_counts_arrays_all_strata$dispersal_events_count_matrices_all_strata[ , ,dispersal_event_types_list, , ]
dim(all_dispersal_events_count_matrices_per_maps_all_strata)
dimnames(all_dispersal_events_count_matrices_per_maps_all_strata)

## Fill the array with all event types
DEC_J_BSM_all_unique_events_source_dest_array[ , ,"range extension (d)", , ] <- all_dispersal_events_count_matrices_per_maps_all_strata[ , ,"range extension (d)", , ]
DEC_J_BSM_all_unique_events_source_dest_array[ , ,"jump-dispersal (j)", , ] <- all_dispersal_events_count_matrices_per_maps_all_strata[ , ,"jump-dispersal (j)", , ]

# DEC_J_BSM_all_unique_events_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Save the final array of counts of events per source/dest areas
# saveRDS(object = DEC_J_BSM_all_unique_events_source_dest_array, file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")
saveRDS(object = DEC_J_BSM_all_unique_events_source_dest_array, file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")


### 5.3/ Build ggplot dataframe of events per source/dest areas ####

## Load the array of counts of events per source/dest areas
# DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")
DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")

## List event types for plots

plot_event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")

# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# Extract counts only for targeted type of events
all_events_count_matrices_per_maps_all_strata <- DEC_J_BSM_all_unique_events_source_dest_array[ , ,plot_event_types_list, , ]

dim(all_events_count_matrices_per_maps_all_strata)
# Rows = Sources
# Cols = Dest

## Format data for ggplot
all_events_count_per_maps_all_bioregions_all_strata_ggplot <- all_events_count_matrices_per_maps_all_strata %>%
  reshape2::melt()
names(all_events_count_per_maps_all_bioregions_all_strata_ggplot) <- c("source", "dest", "event_type", "map", "stratum", "counts")

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
saveRDS(object = all_events_count_per_maps_all_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")


### 5.4/ Aggregate per ggplot dataframe of events per source areas ####

all_events_count_per_maps_all_source_bioregions_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  group_by(map, event_type, source) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_source_bioregions_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_ggplot.rds")
saveRDS(object = all_events_count_per_maps_all_source_bioregions_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_ggplot.rds")

### 5.5/ Aggregate per ggplot dataframe of events per destination areas ####

all_events_count_per_maps_all_dest_bioregions_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  group_by(map, event_type, dest) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()
  
# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_dest_bioregions_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_ggplot.rds")
saveRDS(object = all_events_count_per_maps_all_dest_bioregions_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_ggplot.rds")

### 5.6/ Aggregate per ggplot dataframe of events per source areas x time-strata ####

all_events_count_per_maps_all_source_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  group_by(map, event_type, source, stratum) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")
saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")

### 5.7/ Aggregate per ggplot dataframe of events per destination areas x time-strata ####

all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  group_by(map, event_type, dest, stratum) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")
saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot, file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")


##### 6/ Plot global summary of events counts per bioregion #####

### 6.1/ Plot density curves of events, per bioregion source ####

# Multifaceted: 
  # Rows = type of events
  # Fill of the overlaying histogram/density curve = bioregions

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_all_source_bioregions_ggplot$counts)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
all_events_count_per_maps_all_source_bioregions_ggplot$source <- factor(all_events_count_per_maps_all_source_bioregions_ggplot$source, levels = names(colors_list_for_areas))

## Generate ggplots per type of event
density_curve_all_source_bioregions_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  event_count_all_source_bioregions_per_maps_i <- all_events_count_per_maps_all_source_bioregions_ggplot[all_events_count_per_maps_all_source_bioregions_ggplot$event_type == event_type_i, ]
  
  # Create density curve plot
  density_curve_all_source_bioregions_ggplot <- ggplot(data = event_count_all_source_bioregions_per_maps_i,
                                            mapping = aes(x = counts, fill = source)) +
    # Plot density curve
    geom_density(alpha = 0.8, adjust = 3) + # Use larger bandwidth to smooth density curve
    
    # Set x limits to all be the same across event types
    # xlim(0, max_counts) +
    
    # Use logarithmic scale
    scale_x_log10(limits = c(1, max_counts)) +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregion\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_all_source_bioregions_ggplot <- density_curve_all_source_bioregions_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_all_source_bioregions_ggplot)
  
  # Store ggplot in a list
  density_curve_all_source_bioregions_ggplot_list <- append(x = density_curve_all_source_bioregions_ggplot_list, values = list(density_curve_all_source_bioregions_ggplot))
}

## Multifaceted plot: 1 row = 1 type of events

# ggpubr::ggarrange(density_curve_ggplot_list,        # List of plots
#           # labels = c("A", "B", "C", "D"),   # Add labels
#           ncol = 1, nrow = length(density_curve_ggplot_list)) # Define window division

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_source_bioregions.pdf", width = 10, height = length(density_curve_all_source_bioregions_ggplot_list) * 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_source_bioregions.pdf", width = 10, height = length(density_curve_all_source_bioregions_ggplot_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_source_bioregions_ggplot_list, # List of ggplots
  widths = 4,  # Width of columns
  heights = rep(1, length(density_curve_all_source_bioregions_ggplot_list)),
  nrow = length(density_curve_all_source_bioregions_ggplot_list),
  ncol = 1)

dev.off()


### 6.2/ Plot stacked bars of raw mean event counts, per bioregion source ####

# Plot the stacked bars of event types, per bioregion (raw counts)
   # Fill = type of events
   # X = bioregions

## Format mean data for ggplot
all_events_mean_count_per_maps_all_source_bioregions_ggplot <- all_events_count_per_maps_all_source_bioregions_ggplot %>%
  group_by(source, event_type) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Color scheme for events
event_types_col <- c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid")
event_types_legend_labels <- str_to_sentence(plot_event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))

# Adjust order of event types
all_events_mean_count_per_maps_all_source_bioregions_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_source_bioregions_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# Adjust order of Bioregions
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "East. PA", "West. PA")
all_events_mean_count_per_maps_all_source_bioregions_ggplot$source <- factor(all_events_mean_count_per_maps_all_source_bioregions_ggplot$source, levels = areas_list, labels = bioregion_names_reduced)


# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_source_bioregions.pdf", width = 10, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_source_bioregions.pdf", width = 10, height = 7)

barplot_mean_count_per_events_all_source_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_source_bioregions_ggplot,
                                                                     mapping = aes(y = mean_counts, x = source, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nper Source bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_count_per_events_all_source_bioregions_ggplot <- barplot_mean_count_per_events_all_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_all_source_bioregions_ggplot)

dev.off()


### 6.3/ Plot stacked bars of percentages of mean event counts, per bioregion source ####

# Plot the stacked bars of event types, per bioregion (percentages)
  # Fill = type of events
  # X = bioregions

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_source_bioregions.pdf", width = 10, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_source_bioregions.pdf", width = 10, height = 7)

barplot_mean_perc_per_events_all_source_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_source_bioregions_ggplot,
                                                                     mapping = aes(y = mean_perc, x = source, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of events\nper Source bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("% of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_perc_per_events_all_source_bioregions_ggplot <- barplot_mean_perc_per_events_all_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_perc_per_events_all_source_bioregions_ggplot)

dev.off()


### 6.4/ Plot density curves of events, per bioregion destination ####

# Multifaceted: 
  # Rows = type of events
  # Fill of the overlaying histogram/density curve = bioregions

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_all_dest_bioregions_ggplot$counts)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
all_events_count_per_maps_all_dest_bioregions_ggplot$dest <- factor(all_events_count_per_maps_all_dest_bioregions_ggplot$dest, levels = names(colors_list_for_areas))

# Multifaceted: 
# Rows = type of events
# Fill of the overlaying histogram/density curve = bioregions

## Generate ggplots per type of event
density_curve_all_dest_bioregions_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  event_count_all_dest_bioregions_per_maps_i <- all_events_count_per_maps_all_dest_bioregions_ggplot[all_events_count_per_maps_all_dest_bioregions_ggplot$event_type == event_type_i, ]
  
  # Create density curve plot
  density_curve_all_dest_bioregions_ggplot <- ggplot(data = event_count_all_dest_bioregions_per_maps_i,
                                                       mapping = aes(x = counts, fill = dest)) +
    # Plot density curve
    geom_density(alpha = 0.8, adjust = 3) + # Use larger bandwidth to smooth density curve
    
    # Set x limits to all be the same across event types
    # xlim(0, max_counts) +
    
    # Use logarithmic scale
    scale_x_log10(limits = c(1, max_counts)) +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregion\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_all_dest_bioregions_ggplot <- density_curve_all_dest_bioregions_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_all_dest_bioregions_ggplot)
  
  # Store ggplot in a list
  density_curve_all_dest_bioregions_ggplot_list <- append(x = density_curve_all_dest_bioregions_ggplot_list, values = list(density_curve_all_dest_bioregions_ggplot))
}

## Multifaceted plot: 1 row = 1 type of events

# ggpubr::ggarrange(density_curve_ggplot_list,        # List of plots
#           # labels = c("A", "B", "C", "D"),   # Add labels
#           ncol = 1, nrow = length(density_curve_ggplot_list)) # Define window division

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_dest_bioregions.pdf", width = 10, height = length(density_curve_all_dest_bioregions_ggplot_list) * 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_dest_bioregions.pdf", width = 10, height = length(density_curve_all_dest_bioregions_ggplot_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_dest_bioregions_ggplot_list, # List of ggplots
  widths = 4,  # Width of columns
  heights = rep(1, length(density_curve_all_dest_bioregions_ggplot_list)),
  nrow = length(density_curve_all_dest_bioregions_ggplot_list),
  ncol = 1)

dev.off()


### 6.5/ Plot stacked bars of raw mean event counts, per bioregion destination ####

# Plot the stacked bars of event types, per bioregion (raw and percentage)
  # Fill = type of events
  # X = bioregions

## Format mean data for ggplot
all_events_mean_count_per_maps_all_dest_bioregions_ggplot <- all_events_count_per_maps_all_dest_bioregions_ggplot %>%
  group_by(dest, event_type) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Color scheme for events
event_types_col <- c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid")
event_types_legend_labels <- str_to_sentence(plot_event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))

# Adjust order of event types
all_events_mean_count_per_maps_all_dest_bioregions_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_dest_bioregions_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# Adjust order of Bioregions
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "East. PA", "West. PA")
all_events_mean_count_per_maps_all_dest_bioregions_ggplot$dest <- factor(all_events_mean_count_per_maps_all_dest_bioregions_ggplot$dest, levels = areas_list, labels = bioregion_names_reduced)


# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_dest_bioregions.pdf", width = 10, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_dest_bioregions.pdf", width = 10, height = 7)

barplot_mean_count_per_events_all_dest_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_dest_bioregions_ggplot,
                                                                     mapping = aes(y = mean_counts, x = dest, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nper Destination bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_count_per_events_all_dest_bioregions_ggplot <- barplot_mean_count_per_events_all_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_all_dest_bioregions_ggplot)

dev.off()


### 6.6/ Plot stacked bars of percentages of mean event counts, per bioregion destination ####

# Plot the stacked bars of event types, per bioregion (raw and percentage)
  # Fill = type of events
  # X = bioregions

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_dest_bioregions.pdf", width = 10, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_dest_bioregions.pdf", width = 10, height = 7)

barplot_mean_perc_per_events_all_dest_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_dest_bioregions_ggplot,
                                                                    mapping = aes(y = mean_perc, x = dest, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of events\nper Destination bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("% of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_perc_per_events_all_dest_bioregions_ggplot <- barplot_mean_perc_per_events_all_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_perc_per_events_all_dest_bioregions_ggplot)

dev.off()


### 6.7/ Plot density curves of events, per bioregion Source & Destination ####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Multifaceted: 
  # Rows = type of events
  # Columns = Source and Destination bioregions
  # Fill of the overlaying histogram/density curve = bioregions


# Extract maximum count number across both Source & Dest bioregions
max_source_counts <- max(all_events_count_per_maps_all_source_bioregions_ggplot$counts)
max_dest_counts <- max(all_events_count_per_maps_all_dest_bioregions_ggplot$counts)
max_counts <- max(max_source_counts, max_dest_counts)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
all_events_count_per_maps_all_source_bioregions_ggplot$source <- factor(all_events_count_per_maps_all_source_bioregions_ggplot$source, levels = names(colors_list_for_areas))
all_events_count_per_maps_all_dest_bioregions_ggplot$dest <- factor(all_events_count_per_maps_all_dest_bioregions_ggplot$dest, levels = names(colors_list_for_areas))


## Generate ggplots per type of event for Source Bioregions (remove legend)
density_curve_all_source_bioregions_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  event_count_all_source_bioregions_per_maps_i <- all_events_count_per_maps_all_source_bioregions_ggplot[all_events_count_per_maps_all_source_bioregions_ggplot$event_type == event_type_i, ]
  
  # Create density curve plot
  density_curve_all_source_bioregions_ggplot <- ggplot(data = event_count_all_source_bioregions_per_maps_i,
                                                       mapping = aes(x = counts, fill = source)) +
    # Plot density curve
    geom_density(alpha = 0.8, adjust = 3, # Use larger bandwidth to smooth density curve
                 show.legend = F) + # Remove legend
    
    # Set x limits to all be the same across event types
    # xlim(0, max_counts) +
    
    # Use logarithmic scale
    scale_x_log10(limits = c(1, max_counts)) +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregion\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_all_source_bioregions_ggplot <- density_curve_all_source_bioregions_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_all_source_bioregions_ggplot)
  
  # Store ggplot in a list
  density_curve_all_source_bioregions_ggplot_list <- append(x = density_curve_all_source_bioregions_ggplot_list, values = list(density_curve_all_source_bioregions_ggplot))
}

## Generate ggplots per type of event for Destinatiopn bioregions
density_curve_all_dest_bioregions_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  event_count_all_dest_bioregions_per_maps_i <- all_events_count_per_maps_all_dest_bioregions_ggplot[all_events_count_per_maps_all_dest_bioregions_ggplot$event_type == event_type_i, ]
  
  # Create density curve plot
  density_curve_all_dest_bioregions_ggplot <- ggplot(data = event_count_all_dest_bioregions_per_maps_i,
                                                     mapping = aes(x = counts, fill = dest)) +
    # Plot density curve
    geom_density(alpha = 0.8, adjust = 3) + # Use larger bandwidth to smooth density curve
    
    # Set x limits to all be the same across event types
    # xlim(0, max_counts) +
    
    # Use logarithmic scale
    scale_x_log10(limits = c(1, max_counts)) +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregion\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Counts of events") +
    ylab("Density") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
  
  # Adjust aesthetics
  density_curve_all_dest_bioregions_ggplot <- density_curve_all_dest_bioregions_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(density_curve_all_dest_bioregions_ggplot)
  
  # Store ggplot in a list
  density_curve_all_dest_bioregions_ggplot_list <- append(x = density_curve_all_dest_bioregions_ggplot_list, values = list(density_curve_all_dest_bioregions_ggplot))
}


# Append together Source and Dest multifaceted plots of density curves
density_curve_all_source_dest_bioregions_ggplot_list <- append(x = density_curve_all_source_bioregions_ggplot_list, values = density_curve_all_dest_bioregions_ggplot_list)


# Export plot
# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_source_dest_bioregions.pdf", width = 15, height = length(density_curve_all_source_bioregions_ggplot_list) * 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_source_dest_bioregions.pdf", width = 15, height = length(density_curve_all_source_bioregions_ggplot_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_source_dest_bioregions_ggplot_list, # List of ggplots
  widths = c(6, 9),  # Width of columns
  heights = rep(1, length(density_curve_all_source_bioregions_ggplot_list)),
  nrow = length(density_curve_all_source_bioregions_ggplot_list),
  ncol = 2,
  layout_matrix = rbind(c(1, 5), # Position of ggplots in the layout
                        c(2, 6),
                        c(3, 7),
                        c(4, 8))
  )

dev.off()


### 6.8/ Plot stacked bars of raw mean event counts, per bioregion Source & Destination ####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Plot the stacked bars of event types, per bioregion (raw counts)
  # Fill = type of events
  # X = bioregions
  # Columns = Source and Destination bioregions

## Extract maximum cumulative count number across both Source & Dest bioregions

all_events_cum_counts_per_source_bioregions <- all_events_mean_count_per_maps_all_source_bioregions_ggplot %>% 
  group_by(source) %>%
  summarize(cum_counts = sum(mean_counts))
cum_source_counts <- max(all_events_cum_counts_per_source_bioregions$cum_counts)
all_events_cum_counts_per_dest_bioregions <- all_events_mean_count_per_maps_all_dest_bioregions_ggplot %>% 
  group_by(dest) %>%
  summarize(cum_counts = sum(mean_counts))
cum_dest_counts <- max(all_events_cum_counts_per_dest_bioregions$cum_counts)

max_cum_counts <- max(cum_source_counts, cum_dest_counts)


# Color scheme for events
event_types_col <- c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid")
event_types_legend_labels <- str_to_sentence(plot_event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))

# # Adjust order of event types
# all_events_mean_count_per_maps_all_source_bioregions_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_source_bioregions_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)
# all_events_mean_count_per_maps_all_dest_bioregions_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_dest_bioregions_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# # Adjust order of Bioregions
# bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "East. PA", "West. PA")
# all_events_mean_count_per_maps_all_source_bioregions_ggplot$source <- factor(all_events_mean_count_per_maps_all_source_bioregions_ggplot$source, levels = areas_list, labels = bioregion_names_reduced)


## Generate ggplot per type of event for Source Bioregions (remove legend)
barplot_mean_count_per_events_all_source_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_source_bioregions_ggplot,
                                                                     mapping = aes(y = mean_counts, x = source, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack",
           width = 0.8, col = "black",
           show.legend = F) + # Remove legend
  
  # Set Y-axis
  ylim(0, max_cum_counts) +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nper Source bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_count_per_events_all_source_bioregions_ggplot <- barplot_mean_count_per_events_all_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_all_source_bioregions_ggplot)


## Generate ggplot per type of event for Destination Bioregions
barplot_mean_count_per_events_all_dest_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_dest_bioregions_ggplot,
                                                                   mapping = aes(y = mean_counts, x = dest, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
  
  # Set Y-axis
  ylim(0, max_cum_counts) +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nper Destination bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_count_per_events_all_dest_bioregions_ggplot <- barplot_mean_count_per_events_all_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_count_per_events_all_dest_bioregions_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_mean_count_per_events_all_source_dest_bioregions_ggplot <- append(x = list(barplot_mean_count_per_events_all_source_bioregions_ggplot), values = list(barplot_mean_count_per_events_all_dest_bioregions_ggplot))


## Export multifaceted plot

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_source_dest_bioregions.pdf", width = 17, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_source_dest_bioregions.pdf", width = 17, height = 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_count_per_events_all_source_dest_bioregions_ggplot, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1,
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()


### 6.9/ Plot stacked bars of percentages of mean event counts, per bioregion Source & Destination ####

# Plot results in facing columns
# Remove legend from Source plots in the first columns

# Plot the stacked bars of event types, per bioregion (percentages)
  # Fill = type of events
  # X = bioregions
  # Columns = Source and Destination bioregions

## Generate ggplot per type of event for Source Bioregions (Remove legend)
barplot_mean_perc_per_events_all_source_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_source_bioregions_ggplot,
                                                                    mapping = aes(y = mean_perc, x = source, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack",
           width = 1.0, col = "black",
           show.legend = F) + # Remove legend
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of events\nper Source bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("% of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_perc_per_events_all_source_bioregions_ggplot <- barplot_mean_perc_per_events_all_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_perc_per_events_all_source_bioregions_ggplot)

## Generate ggplot per type of event for Destination Bioregions
barplot_mean_perc_per_events_all_dest_bioregions_ggplot <- ggplot(data = all_events_mean_count_per_maps_all_dest_bioregions_ggplot,
                                                                  mapping = aes(y = mean_perc, x = dest, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of events\nper Destination bioregion\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("% of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_legend_labels, values = event_types_col) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_mean_perc_per_events_all_dest_bioregions_ggplot <- barplot_mean_perc_per_events_all_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_perc_per_events_all_dest_bioregions_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_mean_perc_per_events_all_source_dest_bioregions_ggplot <- append(x = list(barplot_mean_perc_per_events_all_source_bioregions_ggplot), values = list(barplot_mean_perc_per_events_all_dest_bioregions_ggplot))


## Export multifaceted plot

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_source_dest_bioregions.pdf", width = 17, height = 7)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_source_dest_bioregions.pdf", width = 17, height = 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_perc_per_events_all_source_dest_bioregions_ggplot, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1,
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()


##### 7/ Plot global summary of events counts per time strata x bioregion #####

# Load ggplot df
# all_events_count_per_maps_all_source_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")
# all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_source_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")

### 7.1/ Plot the density curves of number of each type of events, per time strata x bioregion sources ####

# Multifaceted: 
   # Rows = type of events
   # Columns = Time-strata
   # Fill of the overlaying histogram/density curve = bioregions

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$counts)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$source <- factor(all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$source, levels = names(colors_list_for_areas))

# Set color scheme for time-strata
time_strata_list <- unique(all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum)
time_strata_col <- setNames(object = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"), nm = time_strata_list)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")

# Set levels of time-strata in defined order
all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum <- factor(all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))


## Generate ggplots per type of event x per time-strata
density_curve_all_source_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  for ( j in seq_along(time_strata_labels))
  {
  
    # j <- 2
    
    # Extract time-strata
    time_strata_j <- rev(time_strata_labels)[j]
    time_strata_full_name_j <- rev(time_strata_full_names)[j]
    
    # Extract event count data
    event_count_all_source_bioregions_all_strata_per_maps_ij <- all_events_count_per_maps_all_source_bioregions_all_strata_ggplot[(all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type == event_type_i) & (all_events_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum == time_strata_j), ]
    
    # Create density curve plot
    density_curve_all_source_bioregions_all_strata_ggplot <- ggplot(data = event_count_all_source_bioregions_all_strata_per_maps_ij,
                                                                    mapping = aes(x = counts, fill = source)) +
      # Plot density curve
      geom_density(alpha = 0.8, adjust = 2) + # Use larger bandwidth to smooth density curve
      
      # Set x limits to all be the same across event types
      # xlim(0, max_counts) +
      
      # Use logarithmic scale
      scale_x_log10(limits = c(1, max_counts)) +
      
      # Adjust color scheme and legend
      scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
      
      # Set plot title +
      ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregion\n", time_strata_full_name_j)) +
      
      # Set axes labels
      xlab("Counts of events") +
      ylab("Density") +
      
      # Adjust margins
      theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
    
    # Adjust aesthetics
    density_curve_all_source_bioregions_all_strata_ggplot <- density_curve_all_source_bioregions_all_strata_ggplot %>% 
      add_aesthetics_density_curve()
    
    # Plot
    # print(density_curve_all_source_bioregions_all_strata_ggplot)
    
    # Store ggplot in a list
    density_curve_all_source_bioregions_all_strata_ggplot_list <- append(x = density_curve_all_source_bioregions_all_strata_ggplot_list, values = list(density_curve_all_source_bioregions_all_strata_ggplot))
    
  }
}

## Multifaceted plot
  # Rows = type of events
  # Columns = Time-strata

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_source_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_source_bioregions_all_strata.pdf",
    width = length(time_strata_labels) * 10, height = length(plot_event_types_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_source_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = rep(1, length(time_strata_labels)),  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = length(time_strata_labels),
  layout_matrix = rbind(c(1:7), # Position of ggplots in the layout
                        c(8:14),
                        c(15:21),
                        c(22:28))
  )

dev.off()


### 7.2/ Plot the stacked bars of event types, per time strata x source bioregions (raw counts) ####
 
# Multifaceted: 
    # Rows = type of events
    # X = Time-strata
    # Fill = bioregions

## Format mean data for ggplot
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_source_bioregions_all_strata_ggplot %>%
  group_by(source, event_type, stratum) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  group_by(event_type, stratum) %>%
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Adjust order of event types
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# Set color scheme for time-strata
time_strata_list <- unique(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum)
time_strata_col <- setNames(object = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"), nm = time_strata_list)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")

# Set levels of time-strata in defined order
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum, levels = time_strata_list, labels = rev(time_strata_labels))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Adjust order of Bioregions
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$source <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$source, levels = areas_list, labels = bioregion_names)


## Generate ggplots per type of events
barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event count data
  barplot_mean_count_per_events_all_source_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_count_per_events_all_source_bioregions_all_strata_per_maps_i,
                                                                                  mapping = aes(y = mean_counts, x = stratum, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Geological epochs") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

  # Adjust aesthetics
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot <- barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()

  # Plot
  # print(barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list, values = list(barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot))

}

 
## Multifaceted plot
  # Rows = type of events
  # Columns = Time-strata

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_source_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_source_bioregions_all_strata.pdf",
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()



### 7.3/ Plot the stacked bars of event types, per time strata x source bioregions (percentages) ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-strata
  # Fill = bioregions


## Generate ggplots per type of events
barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event perc data
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_perc_per_events_all_source_bioregions_all_strata_per_maps_i,
                                                                                  mapping = aes(y = mean_perc, x = stratum, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot <- barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list, values = list(barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot))
  
}


## Multifaceted plot
  # Rows = type of events

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_source_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_source_bioregions_all_strata.pdf",
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.4/ Plot the density curves of number of each type of events, per time strata x dest bioregions  ####

# Multifaceted: 
  # Rows = type of events
  # Columns = Time-strata
  # Fill of the overlaying histogram/density curve = bioregions

# Load ggplot df
# all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")

# Extract maximum count number
max_counts <- max(all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$counts)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest <- factor(all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest, levels = names(colors_list_for_areas))

# Set color scheme for time-strata
time_strata_list <- unique(all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum)
time_strata_col <- setNames(object = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"), nm = time_strata_list)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")

# Set levels of time-strata in defined order
all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum <- factor(all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))


## Generate ggplots per type of event x per time-strata
density_curve_all_dest_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  for ( j in seq_along(time_strata_labels))
  {
    
    # j <- 2
    
    # Extract time-strata
    time_strata_j <- rev(time_strata_labels)[j]
    time_strata_full_name_j <- rev(time_strata_full_names)[j]
    
    # Extract event count data
    event_count_all_dest_bioregions_all_strata_per_maps_ij <- all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot[(all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type == event_type_i) & (all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum == time_strata_j), ]
    
    # Create density curve plot
    density_curve_all_dest_bioregions_all_strata_ggplot <- ggplot(data = event_count_all_dest_bioregions_all_strata_per_maps_ij,
                                                                    mapping = aes(x = counts, fill = dest)) +
      # Plot density curve
      geom_density(alpha = 0.8, adjust = 2) + # Use larger bandwidth to smooth density curve
      
      # Set x limits to all be the same across event types
      # xlim(0, max_counts) +
      
      # Use logarithmic scale
      scale_x_log10(limits = c(1, max_counts)) +
      
      # Adjust color scheme and legend
      scale_fill_manual("Bioregions", breaks = areas_list, labels = bioregion_names, values = colors_list_for_areas) +
      
      # Set plot title +
      ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregion\n", time_strata_full_name_j)) +
      
      # Set axes labels
      xlab("Counts of events") +
      ylab("Density") +
      
      # Adjust margins
      theme(plot.margin = margin(0.5, 1.0, 0.5, 0.5, "cm")) # trbl
    
    # Adjust aesthetics
    density_curve_all_dest_bioregions_all_strata_ggplot <- density_curve_all_dest_bioregions_all_strata_ggplot %>% 
      add_aesthetics_density_curve()
    
    # Plot
    # print(density_curve_all_dest_bioregions_all_strata_ggplot)
    
    # Store ggplot in a list
    density_curve_all_dest_bioregions_all_strata_ggplot_list <- append(x = density_curve_all_dest_bioregions_all_strata_ggplot_list, values = list(density_curve_all_dest_bioregions_all_strata_ggplot))
    
  }
}

## Multifaceted plot
# Rows = type of events
# Columns = Time-strata

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_density_curves_all_dest_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_density_curves_all_dest_bioregions_all_strata.pdf",
    width = length(time_strata_labels) * 10, height = length(plot_event_types_list) * 6)

gridExtra::grid.arrange(
  grobs = density_curve_all_dest_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = rep(1, length(time_strata_labels)),  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = length(time_strata_labels),
  layout_matrix = rbind(c(1:7), # Position of ggplots in the layout
                        c(8:14),
                        c(15:21),
                        c(22:28))
)

dev.off()


### 7.5/ Plot the stacked bars of event types, per time strata x dest bioregions (raw counts) ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-strata
  # Fill = bioregions

## Format mean data for ggplot
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot %>%
  group_by(dest, event_type, stratum) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  group_by(event_type, stratum) %>%
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

# Adjust order of event types
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# Set color scheme for time-strata
time_strata_list <- unique(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum)
time_strata_col <- setNames(object = c("#ffff99", "#f8eb88", "#fbb697", "#fac06a", "#f7923f", "#99cc00", "#0d731e"), nm = time_strata_list)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")

# Set levels of time-strata in defined order
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum, levels = time_strata_list, labels = rev(time_strata_labels))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Adjust order of Bioregions
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest, levels = areas_list, labels = bioregion_names)


## Generate ggplots per type of events
barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event count data
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_count_per_events_all_dest_bioregions_all_strata_per_maps_i,
                                                                                  mapping = aes(y = mean_counts, x = stratum, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot <- barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list, values = list(barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot))
  
}


## Multifaceted plot
  # Rows = type of events
  # Columns = Time-strata

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_dest_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_dest_bioregions_all_strata.pdf",
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.6/ Plot the stacked bars of event types, per time strata x dest bioregions (percentages) ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-strata
  # Fill = bioregions


## Generate ggplots per type of events
barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event perc data
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_per_maps_i,
                                                                                 mapping = aes(y = mean_perc, x = stratum, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot <- barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list, values = list(barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot))
  
}


## Multifaceted plot
  # Rows = type of events

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_dest_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_dest_bioregions_all_strata.pdf",
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.7/ Plot stacked bars of raw mean event counts, per time-strata x bioregion Sources & Destinations ####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Plot the stacked bars of event types, per bioregion (raw counts)
  # Rows = type of events
  # Columns = Source and Destination bioregions
  # X = Time-strata
  # Fill = bioregions
  
# Load ggplot df
# all_events_count_per_maps_all_source_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")
# all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_source_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot.rds")

## Format mean data for ggplot
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_source_bioregions_all_strata_ggplot %>%
  group_by(source, event_type, stratum) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  group_by(event_type, stratum) %>%
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_strata_ggplot %>%
  group_by(dest, event_type, stratum) %>% 
  summarize(mean_counts = mean(counts)) %>% 
  group_by(event_type, stratum) %>%
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%
  ungroup()

## Extract maximum cumulative count number per time-strata across both Source & Dest bioregions

all_events_maxcum_counts_per_source_bioregions_all_strata_ggplot <- all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot %>% 
  group_by(stratum, event_type) %>%
  summarize(cum_counts = sum(mean_counts)) %>%
  group_by(event_type) %>%
  summarize(max_counts_source = max(cum_counts))
all_events_maxcum_counts_per_dest_bioregions_all_strata_ggplot <- all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot %>% 
  group_by(stratum, event_type) %>%
  summarize(cum_counts = sum(mean_counts)) %>%
  group_by(event_type) %>%
  summarize(max_counts_dest = max(cum_counts))

maxcum_source_dest_counts <- left_join(x = all_events_maxcum_counts_per_source_bioregions_all_strata_ggplot, all_events_maxcum_counts_per_dest_bioregions_all_strata_ggplot)
maxcum_source_dest_counts
# Counts of events should be similar between source and destination
max_cum_counts <- all_events_maxcum_counts_per_source_bioregions_all_strata_ggplot$max_counts_source


# Set labels for event types
event_types_legend_labels <- str_to_sentence(plot_event_types_list)
# event_types_legend_labels[event_types_legend_labels == "Jump-dipsersal (j)"] <- "Jump-dispersal (j)"
event_types_axis_labels <- paste0("(", str_remove(string = event_types_legend_labels, pattern = ".* \\("))

# Adjust order of event types
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type, levels = plot_event_types_list, labels = event_types_axis_labels)

# Set labels for time-strata
time_strata_list <- unique(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum)
time_strata_labels <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")
time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")

# Set levels of time-strata in defined order
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$stratum, levels = rev(time_strata_list), labels = rev(time_strata_labels))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Adjust order of Bioregions
all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$source <- factor(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$source, levels = areas_list, labels = bioregion_names)
all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest <- factor(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$dest, levels = areas_list, labels = bioregion_names)


## Generate ggplots per type of events for Source Bioregions (remove legend)
barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event count data
  barplot_mean_count_per_events_all_source_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_count_per_events_all_source_bioregions_all_strata_per_maps_i,
                                                                                  mapping = aes(y = mean_counts, x = stratum, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack",
             width = 0.8, col = "black", 
             show.legend = F) +  # Remove legend
    
    # Set Y-axis
    ylim(0, max_cum_counts[i]) +

    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot <- barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list, values = list(barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event count data
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_count_per_events_all_dest_bioregions_all_strata_per_maps_i,
                                                                                mapping = aes(y = mean_counts, x = stratum, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 0.8, col = "black") +
    
    # Set Y-axis
    ylim(0, max_cum_counts[i]) +
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot <- barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list, values = list(barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_mean_count_per_events_all_source_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_count_per_events_all_source_bioregions_all_strata_ggplot_list, values = barplot_mean_count_per_events_all_dest_bioregions_all_strata_ggplot_list)


## Export multifaceted plot
  # Rows = type of events
  # Columns = Source and Destination bioregions

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_per_types_barplots_all_source_dest_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_count_per_types_barplots_all_source_dest_bioregions_all_strata.pdf",
    width = 17, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_count_per_events_all_source_dest_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = rep(1, length(plot_event_types_list)), # Height of rows
  nrow = length(plot_event_types_list),
  ncol = 2,
  layout_matrix = rbind(c(1, 5),
                        c(2, 6),
                        c(3, 7),
                        c(4, 8)) # Position of ggplots in the layout
)

dev.off()

### 7.8/ Plot stacked bars of percentages of mean event counts, per time-strata x bioregion Sources & Destinations ####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Plot the stacked bars of event types, per bioregion (raw counts)
  # Rows = type of events
  # Columns = Source and Destination bioregions
  # X = Time-strata
  # Fill = bioregions

## Generate ggplots per type of events for Source Bioregions (remove legend)
barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event perc data
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_source_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_perc_per_events_all_source_bioregions_all_strata_per_maps_i,
                                                                                  mapping = aes(y = mean_perc, x = stratum, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack",
             width = 1.0, col = "black", 
             show.legend = F) +  # Remove legend
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot <- barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list, values = list(barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- event_types_axis_labels[i]
  
  # Extract event perc data
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_per_maps_i <- all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot[(all_events_mean_count_per_maps_all_dest_bioregions_all_strata_ggplot$event_type == event_type_i), ]
  
  # Generate plot
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot <- ggplot(data = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_per_maps_i,
                                                                                mapping = aes(y = mean_perc, x = stratum, fill = dest)) +
    # Plot stacked bars
    geom_col(alpha = 1.0, position = "stack", width = 1.0, col = "black") +
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Geological epochs") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot <- barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot %>% 
    add_aesthetics_density_curve()
  
  # Plot
  # print(barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot)
  
  # Store ggplot in a list
  barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list, values = list(barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean percs per event types
barplot_mean_perc_per_events_all_source_dest_bioregions_all_strata_ggplot_list <- append(x = barplot_mean_perc_per_events_all_source_bioregions_all_strata_ggplot_list, values = barplot_mean_perc_per_events_all_dest_bioregions_all_strata_ggplot_list)


## Export multifaceted plot
  # Rows = type of events
  # Columns = Source and Destination bioregions

# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Events_perc_per_types_barplots_all_source_dest_bioregions_all_strata.pdf",
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Events_perc_per_types_barplots_all_source_dest_bioregions_all_strata.pdf",
    width = 18, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_mean_perc_per_events_all_source_dest_bioregions_all_strata_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = rep(1, length(plot_event_types_list)), # Height of rows
  nrow = length(plot_event_types_list),
  ncol = 2,
  layout_matrix = rbind(c(1, 5),
                        c(2, 6),
                        c(3, 7),
                        c(4, 8)) # Position of ggplots in the layout
)

dev.off()


##### 8/ Plot Emigration/Immigration counts between bioregions overall ####

### 8.1/ Aggregate count data across all time-strata ####

# Load imputed phylogeny with short names
# Ponerinae_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")
Ponerinae_MCC_phylogeny_1534t_short_names <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")

# Extract root age
# root_age <- max(phytools::nodeHeights(Ponerinae_phylogeny_1534t_short_names))
root_age <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t_short_names))

# Load ggplot df
# all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
names(all_events_count_per_maps_all_bioregions_all_strata_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Sum across event_types x strata
# Aggregate mean counts across maps
all_dispersal_events_mean_count_per_maps_overall_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(source, dest, map) %>% 
  summarize(counts = sum(counts)) %>% # Sum across event_types x strata
  group_by(source, dest) %>%
  summarize(mean_counts = mean(counts)) %>% # Aggregate mean counts across maps
  ungroup() %>% # Compute overall percentages
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts))
  
# Save ggplot df of mean counts of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_overall_ggplot, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")
saveRDS(object = all_dispersal_events_mean_count_per_maps_overall_ggplot, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")

## Convert back to array

# Pivot mean counts and percentages
all_dispersal_events_mean_count_per_maps_overall_melted <- all_dispersal_events_mean_count_per_maps_overall_ggplot %>% 
  pivot_longer(cols = c("mean_counts", "mean_perc"), names_to = "mean_stats", values_to = "values")

all_dispersal_events_mean_count_per_maps_overall_array <- reshape2::acast(data = all_dispersal_events_mean_count_per_maps_overall_melted,
                                                                          formula = source ~ dest ~ mean_stats)

dimnames(all_dispersal_events_mean_count_per_maps_overall_array)
all_dispersal_events_mean_count_per_maps_overall_array[, , "mean_counts"]
all_dispersal_events_mean_count_per_maps_overall_array[, , "mean_perc"]

# Save array of mean counts/perc of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_overall_array, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_array.rds")
saveRDS(object = all_dispersal_events_mean_count_per_maps_overall_array, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_array.rds")

# Rows = Sources
# Cols = Dest

### 8.2/ Compute node metadata, including richness per bioregions ####

# Load LTT data
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

bioregion_names_alpha_order <- c("Afrotropics", "Eastern Palearctic", "Indomalaya", "Neotropics", "Nearctic", "Australasia", "Western Palearctic")

nodes_metadata <- DEC_J_LTT_all_areas_mean_ggplot %>% 
  filter(time == 0) %>% 
  filter(areas != "total") %>% 
  rename(node_ID = areas) %>%
  mutate(bioregions = bioregion_names_alpha_order) %>%
  mutate(node_size = log(mean_counts)) %>%
  select(node_ID, bioregions, mean_counts, mean_percentages, node_size)

## Add label with richness
nodes_metadata$node_labels <- paste0(nodes_metadata$node_ID, "\n", round(nodes_metadata$mean_counts, 0))

nodes_metadata

log(nodes_metadata$mean_counts)

## Add spatial coordinates for each bioregion
nodes_metadata$bioregions
nodes_metadata$latitude <- c(0, 50, 10, -10, 40, -25, 47)
nodes_metadata$longitude <- c(25, 100, 100, -60, -100, 140, 20)
nodes_metadata$adjusted_latitude <- c(10, 35, 10, -10, 35, -15, 40)
nodes_metadata$adjusted_longitude <- c(15, 120, 100, -60, -100, 150, 15)

## Save node metadata
# saveRDS(object = nodes_metadata, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
saveRDS(object = nodes_metadata, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

### 8.3/ Plot as a network ####

# https://kateto.net/network-visualization

# Load node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

## 8.3.1/ Convert to igraph object ####

# Round counts and convert to edge size
all_dispersal_events_overall_df <- all_dispersal_events_mean_count_per_maps_overall_ggplot %>% 
  mutate(mean_counts = round(mean_counts, 0)) %>% 
  mutate(mean_perc = round(mean_perc, 1)) %>% 
  mutate(edge_width = log1p(mean_counts)) # Scale to log
  # mutate(edge_width = mean_counts / max(mean_counts)) # Standardize between 0 and 1
  
# Save edge metadata df for all dispersal events across all time-strata
# saveRDS(object = all_dispersal_events_overall_df, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
saveRDS(object = all_dispersal_events_overall_df, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")

## Convert to igraph
all_dispersal_events_overall_igraph <- igraph::graph_from_data_frame(d = all_dispersal_events_overall_df, vertices = nodes_metadata, directed = T)
str(all_dispersal_events_overall_igraph)

# Explore igraph
all_dispersal_events_overall_igraph
E(all_dispersal_events_overall_igraph) # Edge attributes
V(all_dispersal_events_overall_igraph) # Vertex/Node attributes
all_dispersal_events_overall_igraph[]  # Adjacency matrix of nodes x nodes

# Export adjacency matrix
as_adjacency_matrix(all_dispersal_events_overall_igraph, attr = "mean_counts")
as_adjacency_matrix(all_dispersal_events_overall_igraph, attr = "mean_perc")

## 8.3.2/ Adjust aesthetics ####

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[nodes_metadata$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[nodes_metadata$node_ID]

# Remove loops
all_dispersal_events_overall_igraph <- simplify(graph = all_dispersal_events_overall_igraph,
                                                remove.multiple = F, remove.loops = T) 

# Set node colors based on bioregions:
V(all_dispersal_events_overall_igraph)$color <- colors_list_for_areas
# Set node size based on current richness
V(all_dispersal_events_overall_igraph)$size <- V(all_dispersal_events_overall_igraph)$node_size * 5
# Set node font size based on current richness
V(all_dispersal_events_overall_igraph)$label.cex <- V(all_dispersal_events_overall_igraph)$node_size / 5

# Set edge width based on event counts
# E(all_dispersal_events_overall_igraph)$width <- E(all_dispersal_events_overall_igraph)$edge_width * 10
E(all_dispersal_events_overall_igraph)$width <- E(all_dispersal_events_overall_igraph)$edge_width * 2 + 1

# Set edge arrow size based on event counts
# E(all_dispersal_events_overall_igraph)$arrow.size <- E(all_dispersal_events_overall_igraph)$edge_width * 4 - 0.3
E(all_dispersal_events_overall_igraph)$arrow.size <- E(all_dispersal_events_overall_igraph)$edge_width / 4 + 0.1

# Set edge label as mean counts
# E(all_dispersal_events_overall_igraph)$edge.label <- E(all_dispersal_events_overall_igraph)$mean_counts

# Set edge color according to their source
edge_starts <- ends(all_dispersal_events_overall_igraph, es = E(all_dispersal_events_overall_igraph), names = F)[, 1]
edge_col <- V(all_dispersal_events_overall_igraph)$color[edge_starts]
E(all_dispersal_events_overall_igraph)$color <- edge_col

# Remove edges with low counts
table(E(all_dispersal_events_overall_igraph)$mean_counts)
cutoff_min_counts <- 5 
all_dispersal_events_overall_igraph <- delete_edges(graph = all_dispersal_events_overall_igraph,
                                                    edges = E(all_dispersal_events_overall_igraph)[mean_counts < cutoff_min_counts])

# Adjust layout to coordinates of vertex/nodes
layout_WGS84 <- cbind(V(all_dispersal_events_overall_igraph)$longitude, V(all_dispersal_events_overall_igraph)$latitude)
layout_pretty <- cbind(V(all_dispersal_events_overall_igraph)$adjusted_longitude, V(all_dispersal_events_overall_igraph)$adjusted_latitude)

## Save igraph for all dispersal events across all time-strata
# saveRDS(object = all_dispersal_events_overall_igraph, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.rds")
saveRDS(object = all_dispersal_events_overall_igraph, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.rds")

## Plot igraph

# pdf(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf",
pdf(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf",
        
    width = 6, height = 6)

plot.igraph(x = all_dispersal_events_overall_igraph,
            ## Node aesthetics
            # vertex.color = colors_list_for_areas, # Nodes color
            vertex.frame.color = "black", # Nodes border color
            vertex.shape = "circle",  # Node shape
            # vertex.size = 15, # Node size
            vertex.label	= nodes_metadata$node_labels, # Node labels
            vertex.label.color = "black",  # Node label color
            vertex.label.font = 2, # Node label font
            # vertex.label.cex = 1, # Node label cex
            ## Edge aesthetics
            # edge.color = "black",     # Edge color
            # edge.width = 1,         # Edge width
            # edge.arrow.size = 0.7,    # Terminal arrow size
            edge.label = E(all_dispersal_events_overall_igraph)$mean_counts,
            edge.label.cex = 0.7,
            edge.lty = 1,             # Edge line type
            edge.curved = 0.1,
            arrow.mode = "forward",    # Direction of arrows
            ## Other aesthetics
            # layout = layout_WGS84,
            layout = layout_pretty,
            main = paste0("Dispersal events between bioregions\n",
                          "All epochs (",round(root_age,0),"-0 My)")
            )

dev.off()


### 8.4/ With ggplot to overlay on continent maps ####

# Section 8 in https://kateto.net/network-visualization

## 8.4.1/ Build bioregion map ####

## Load worldwide adm0 map
Countries_NE_sf <- rnaturalearth::ne_load(file_name = "ne_10m_admin_0_countries", returnclass = 'sf', destdir = "./input_data/NaturalEarth_maps/Natural_Earth_quick_start/10m_cultural/")
plot(Countries_NE_sf["ISO_A3"])

## Load worldwide adm1 map
geoBoundaries_adm1_sf <- st_read(dsn = "./input_data/geoBoundaries/geoBoundariesCGAZ_ADM1/", layer = "geoBoundariesCGAZ_ADM1")
# plot(geoBoundaries_adm1_sf["shapeName"])
# plot(geoBoundaries_adm1_sf["shapeGroup"])
table(geoBoundaries_adm1_sf$shapeGroup)

## Intersect country and adm1 maps to keep territory without adm1
sf_use_s2(FALSE)
Bioregions_sf <- st_buffer(st_make_valid(st_intersection(Countries_NE_sf[, "ISO_A3"], geoBoundaries_adm1_sf)), 0)
sf_use_s2(TRUE)

# plot(Bioregions_sf[, "ISO_A3"])

Bioregions_sf <- Bioregions_sf %>% 
  rename(ISO_A3_for_GB = shapeGroup) %>% 
  rename(ISO_A3_for_NE = ISO_A3)

names(Bioregions_sf)

## Load metadata with bioregion assignment
Countries_NE_sf_metadata <- openxlsx::read.xlsx(xlsxFile = "./input_data/Biogeographic_data/Countries_NE_sf_metadata.xlsx")

Countries_NE_sf_metadata$ISO_A3[!(Countries_NE_sf_metadata$ISO_A3 %in% Bioregions_sf$ISO_A3_for_NE)]
Bioregions_sf$ISO_A3_for_NE[!(Bioregions_sf$ISO_A3_for_NE %in% Countries_NE_sf_metadata$ISO_A3)]

##  Assign bioregions using ISO_A3 from Geoboundaries
Bioregions_sf <- left_join(x = Bioregions_sf, y = Countries_NE_sf_metadata[, c("ISO_A3", "Bioregion_7_PaleA")], by = join_by(ISO_A3_for_GB == ISO_A3)) %>% 
  rename(Bioregion_for_GB = Bioregion_7_PaleA)

## Assign bioregions using ISO_A3 from NaturalEarth
Bioregions_sf <- left_join(x = Bioregions_sf, y = Countries_NE_sf_metadata[, c("ISO_A3", "Bioregion_7_PaleA")], by = join_by(ISO_A3_for_NE == ISO_A3)) %>% 
  rename(Bioregion_for_NE = Bioregion_7_PaleA)
View(Bioregions_sf)

# Fill missing Bioregions from Geoboundaries with match from NaturalEarth
Bioregions_sf$Bioregion <- Bioregions_sf$Bioregion_for_GB
Bioregions_sf$Bioregion[is.na(Bioregions_sf$Bioregion)] <- Bioregions_sf$Bioregion_for_NE[is.na(Bioregions_sf$Bioregion)]
View(Bioregions_sf)

Bioregions_sf[Bioregions_sf$ISO_A3 == "-99", ]

## Deal with the particular cases of countries overlapping bioregions

View(Bioregions_sf[is.na(Bioregions_sf$Bioregion), ] %>% arrange(ISO_A3_for_GB))

## Disputed regions in Himalaya between India and China
plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB %in% c("113", "112", "114", "116", "119", "121") | Bioregions_sf$shapeName %in% c("Him?chal Pradesh", "Jammu and Kashm?r", "Lad?kh"), "shapeID"])
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB %in% c("113", "112", "114", "116", "119", "121")] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB %in% c("113", "112", "114", "116", "119", "121")] <- "Eastern Palearctic"

## Code 120: Isla Brasilera = Disputed island between Brazil, Uruguay and Argentina
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB %in% c("120")] <- "Neotropics"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB %in% c("120")] <- "Neotropics"

## Code 125: Paracel Islands = Disputed archipelago between China and the Philippines in the South China Sea
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB %in% c("125")] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB %in% c("125")] <- "Indomalaya"

## Code CHN: China => Split between Eastern Palearctic and Indomalaya

China_adm1_in_Indomalaya <- c("Chongqing", "Chongqing Municipality", "Fujian", "Fujian Province", "Guangdong", "Guangzhou Province",
                              "Guangxi", "Guangxi Zhuang Autonomous Region", "Guizhou", "Guizhou Province", "Hainan", "Hainan Province",
                              "Hong Kong", "Hong Kong Special Administrative Region", "Hubei", "Hubei Province", "Hunan", "Hunan Province",
                              "Jiangxi", "Jiangxi Province", "Yunnan", "Yunnan Province", "Zhejiang", "Zhejiang Province")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "CHN"] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "CHN"] <- "Eastern Palearctic"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% China_adm1_in_Indomalaya] <- "Indomalaya"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% China_adm1_in_Indomalaya] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% China_adm1_in_Indomalaya] <- "Indomalaya"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "CHN", "Bioregion"])

## Code FRA: France => All Off-sea territories have their own code, so all FRA is in Western Palearctic
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "FRA"] <- "Western Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "FRA"] <- "Western Palearctic"


## Code IDN: Indonesia => Use adm1 to distinguish between IndoMalaya and Australasia according to AntWeb

Indonesia_adm1_in_Australasia <- c("Maluku", "North Maluku", "East Nusa Tenggara", "West Nusa Tenggara", "Papua", "West Papua")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "IDN"] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "IDN"] <- "Indomalaya"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% Indonesia_adm1_in_Australasia] <- "Australasia"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% Indonesia_adm1_in_Australasia] <- "Australasia"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% Indonesia_adm1_in_Australasia] <- "Australasia"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "IDN", "Bioregion"])

## Code IND: India => : Three Himalayan regions in Eastern Palearctic

India_adm1_in_Palearctic <- c("Him?chal Pradesh", "Himachal Pradesh", "Jammu and Kashm?r", "Jammu and Kashmir", "Ladakh", "Lad?kh")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "IND"] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "IND"] <- "Indomalaya"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% India_adm1_in_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% India_adm1_in_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% India_adm1_in_Palearctic] <- "Eastern Palearctic"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "IND", "Bioregion"])

## Code JPN: Japan => Okinawa Prefecture is in Indomalaya. Ogasawara Islands should be in Australasia, but since they are part of Tokyo Prefecture, we ignore them here.

Japan_adm1_in_IndoMalaya <- c("Okinawa Prefecture")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "JPN"] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "JPN"] <- "Eastern Palearctic"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% Japan_adm1_in_IndoMalaya] <- "Indomalaya"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% Japan_adm1_in_IndoMalaya] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% Japan_adm1_in_IndoMalaya] <- "Indomalaya"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "JPN", "Bioregion"])
# View(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "JPN", ])

## Code MEX: Mexico => Northern part is in Nearctic

Mexico_adm1_in_Nearctic <- c("Aguascalientes", "Baja California", "Baja California Sur", "Chihuahua", "Coahuila de Zaragoza", "Durango", "Guanajuato", "Hidalgo", "México", "MÃ©xico",
                             "Nuevo León", "Nuevo LeÃ³n", "Querétaro", "QuerÃ©taro de Arteaga", "Querétaro de Arteaga", "San Luis Potosí", "San Luis PotosÃ­", "San Luis PotosÃ", "Sonora", "Tamaulipas", "Tlaxcala", "Zacatecas")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "MEX"] <- "Neotropics"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "MEX"] <- "Neotropics"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% Mexico_adm1_in_Nearctic] <- "Nearctic"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% Mexico_adm1_in_Nearctic] <- "Nearctic"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% Mexico_adm1_in_Nearctic] <- "Nearctic"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "MEX", "Bioregion"])

## Code NLD: Netherlands => All Off-sea territories have their own code, so all NLD is in Western Palearctic
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "NLD"] <- "Western Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "NLD"] <- "Western Palearctic"

## Code PAK: Pakistan => Two Himalayan regions in Eastern Palearctic

Pakistan_adm1_in_Palearctic <- c("Azad Kashmir", "Gilgit-Baltistan")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "PAK"] <- "Indomalaya"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "PAK"] <- "Indomalaya"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% Pakistan_adm1_in_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% Pakistan_adm1_in_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% Pakistan_adm1_in_Palearctic] <- "Eastern Palearctic"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "PAK", "Bioregion"])

## Code RUS: Russia => Split between Eastern and Western Paleartic at Oural mountains

Russia_adm1_in_Eastern_Palearctic <- c("Altai Krai", "Altai Republic", "Amur Oblast", "Buryatia", "Chelyabinsk Oblast", "Chukotka Autonomous Okrug", "Irkutsk Oblast", "Jewish Autonomous Oblast", "Kamchatka Krai", "Kemerovo Oblast", "Khabarovsk Krai", "Khakassia", "Khanty-Mansiysk Autonomous Okrug ? Ugra", "Krasnoyarsk Krai", "Kurgan Oblast", "Magadan Oblast", "Novosibirsk Oblast", "Omsk Oblast", "Primorsky Krai", "Sakha Republic", "Sakhalin Oblast", "Sverdlovsk Oblast", "Tomsk Oblast", "Tuva", "Tyumen Oblast", "Yamalo-Nenets Autonomous Okrug", "Zabaykalsky Krai")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "RUS"] <- "Western Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "RUS"] <- "Western Palearctic"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% Russia_adm1_in_Eastern_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% Russia_adm1_in_Eastern_Palearctic] <- "Eastern Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% Russia_adm1_in_Eastern_Palearctic] <- "Eastern Palearctic"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "RUS", "Bioregion"])

## Code USA: United States of America => Hawaii is in Australasia

USA_adm1_in_Australasia <- c("Hawaii")

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "USA"] <- "Nearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "USA"] <- "Nearctic"

Bioregions_sf$Bioregion_for_GB[Bioregions_sf$shapeName %in% USA_adm1_in_Australasia] <- "Australasia"
Bioregions_sf$Bioregion_for_NE[Bioregions_sf$shapeName %in% USA_adm1_in_Australasia] <- "Australasia"
Bioregions_sf$Bioregion[Bioregions_sf$shapeName %in% USA_adm1_in_Australasia] <- "Australasia"

# plot(Bioregions_sf[Bioregions_sf$ISO_A3_for_GB == "USA", "Bioregion"])

## Code XKX: Kosovo in Western Palearctic
Bioregions_sf$Bioregion_for_GB[Bioregions_sf$ISO_A3_for_GB == "XKX"] <- "Western Palearctic"
Bioregions_sf$Bioregion[Bioregions_sf$ISO_A3_for_GB == "XKX"] <- "Western Palearctic"

# View(Bioregions_sf[is.na(Bioregions_sf$Bioregion), ] %>% arrange(ISO_A3_for_GB))
# plot(Bioregions_sf[, "Bioregion"])

## Save adm1 maps with bioregion assignments
Bioregions_sf_adm1_level <- Bioregions_sf
saveRDS(object = Bioregions_sf_adm1_level, file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level.rds")

## Union features at Country level

Bioregions_sf_adm1_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level.rds")

sf_use_s2(FALSE)
Bioregions_sf_Country_level <- Bioregions_sf_adm1_level %>%
  group_by(ISO_A3_for_GB, Bioregion) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Bioregions_sf_Country_level <- st_buffer(st_make_valid(Bioregions_sf_Country_level), 0)
sf_use_s2(TRUE)

plot(Bioregions_sf_Country_level[, "Bioregion"])

saveRDS(object = Bioregions_sf_Country_level, file = "./input_data/geoBoundaries/Bioregions_sf_Country_level.rds")

## Union features at Bioregion level

Bioregions_sf_adm1_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_adm1_level.rds")

sf_use_s2(FALSE)
Bioregions_sf_Bioregions_level <- Bioregions_sf_adm1_level %>%
  group_by(Bioregion) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Bioregions_sf_Bioregions_level <- st_buffer(st_make_valid(Bioregions_sf_Bioregions_level), 0)
sf_use_s2(TRUE)

plot(Bioregions_sf_Bioregions_level[, "Bioregion"])

saveRDS(object = Bioregions_sf_Bioregions_level, file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")


### 8.4.2/ Map network with bioregion map as background ####

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

## Load node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

# Load edge metadata
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")


## Load igraph for all dispersal events across all time-strata
# all_dispersal_events_overall_igraph <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.rds")
all_dispersal_events_overall_igraph <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.rds")

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[nodes_metadata$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
# names(colors_list_for_areas) <- nodes_metadata$bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
# colors_list_for_areas <- colors_list_for_areas[bioregion_names]

# Adjust order of bioregions
nodes_metadata$bioregions <- factor(nodes_metadata$bioregions, levels = bioregion_names, labels = bioregion_names)
Bioregions_sf_Bioregions_level$Bioregion <- factor(Bioregions_sf_Bioregions_level$Bioregion, levels = bioregion_names, labels = bioregion_names)

# Convert Long/Lat WGS84 to Mollweide
nodes_metadata_sf <- st_as_sf(x = nodes_metadata, coords = c("longitude", "latitude"), remove = F, crs = st_crs(Bioregions_sf_Bioregions_level))
nodes_metadata_sf <- st_transform(nodes_metadata_sf, crs = st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
Mollweide_coordinates <- st_as_text(nodes_metadata_sf$geometry)
longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
nodes_metadata$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
nodes_metadata$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))

## Inform edge metadata with coordinates and bioregion source labels
all_dispersal_events_overall_df <- left_join(x = all_dispersal_events_overall_df,
                                             y = nodes_metadata[, c("node_ID", "bioregions", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                             by = join_by(source == node_ID))
all_dispersal_events_overall_df <- all_dispersal_events_overall_df %>% 
  rename(source_labels = bioregions,
         source_latitude = latitude,
         source_longitude = longitude,
         source_latitude_Mollweide = latitude_Mollweide,
         source_longitude_Mollweide = longitude_Mollweide)
all_dispersal_events_overall_df <- left_join(x = all_dispersal_events_overall_df,
                                             y = nodes_metadata[, c("node_ID", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                             by = join_by(dest == node_ID))
all_dispersal_events_overall_df <- all_dispersal_events_overall_df %>% 
  rename(dest_latitude = latitude,
         dest_longitude = longitude,
         dest_latitude_Mollweide = latitude_Mollweide,
         dest_longitude_Mollweide = longitude_Mollweide)

# Adjust order of bioregions
all_dispersal_events_overall_df$source_labels <- factor(all_dispersal_events_overall_df$source_labels, levels = bioregion_names, labels = bioregion_names)


# Save edge metadata df
# saveRDS(object = all_dispersal_events_overall_df, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
saveRDS(object = all_dispersal_events_overall_df, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")


# May adjust curvature manually
# May adjust arrow size too

# Load edge metadata df
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")


# Set minimum nb of dispersal events to display
min_counts_threshold <- 5

## GGplot
# pdf(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf",
pdf(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf",
    width = 9, height = 4)

all_dispersal_events_overall_ggplot <- ggplot(data = nodes_metadata) +
  
  # Plot bioregion maps
  geom_sf(data = Bioregions_sf_Bioregions_level,
          mapping = aes(fill = Bioregion),
          colour = "black",
          alpha = 0.8) +
  
  # Adjust fill color scheme for bioregions
  scale_fill_manual("Bioregions", breaks = bioregion_names, labels = bioregion_names, values = colors_list_for_areas) +
  
  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of Wrodl map
           expand = FALSE) +
  
  # Plot nodes
  geom_point(mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide, size = mean_counts),
             alpha = 0.3, show.legend = F) +
  
  # Adjust size for points
  scale_size_continuous("Species richness", range = c(5.94, 30)) +
  
  # Plot vertices
  geom_curve(data = all_dispersal_events_overall_df[all_dispersal_events_overall_df$mean_counts >= min_counts_threshold, ],
             aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                 xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                 color = source_labels, linewidth = mean_counts),
             arrow = arrow(length = unit(0.03, "npc"), 
                           type = "closed"), # Describes arrow head (open or closed)
             angle = 90, # Anything other than 90 or 0 can look unusual
             alpha = 1.0
             ) + 
  
  # Adjust color scheme for arrows
  scale_color_manual("Source bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Adjust edge width legend
  scale_linewidth_continuous("Dispersal events",
                             breaks = c(5, 25, 75),
                             range = c(1, 6)) +
  
  # Add node labels
  ggnewscale::new_scale(new_aes = "size") +
  geom_text(mapping = aes(label = round(mean_counts, 0),
                          x = longitude_Mollweide, y = latitude_Mollweide,
                          size = mean_counts),
            # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
            color = "black",
            fontface = 2, show.legend = F) +
  
  # Adjust size for labels
  scale_size_continuous("Species richness", range = c(3, 10)) +
  
  # Add title
  ggtitle(label = paste0("Dispersal events between bioregions\n",
                         "All epochs (",round(root_age, 0),"-0 My)")) +
  
  # Adjust legend aesthetics
  guides(color = "none",
         size = "none",
         fill = guide_legend(order = 1),
         linewidth = guide_legend(order = 2,
                                  override.aes = list(arrow = 0.2))
         ) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        # legend.position = c(1.1, 0.9),
        legend.box.margin = margin(l = 15, b = 15),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())
 
print(all_dispersal_events_overall_ggplot)

dev.off()

# Map edges as geom_path between nodes or GreatCircle arc with geosphere::gcIntermediate

## Make an alternative version with labels just to ease conversion to PPT
# Add labels for node richness
# Add labels for edge values

# Can even become interactive with leaflet?
# To display richness by clicking?

### 8.5/ Plot as a Chord diagram ####
  
# Use bidirectional chord diagrams?

### 8.6/ Plot mean counts of Emigration events between (Source) bioregions ####

# Load melted df of mean counts of dispersal events between bioregions
# all_dispersal_events_mean_count_per_maps_overall_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")
all_dispersal_events_mean_count_per_maps_overall_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")

# Aggregate counts per source bioregions
all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot <- all_dispersal_events_mean_count_per_maps_overall_ggplot %>% 
  group_by(source) %>% 
  summarize(mean_counts = sum(mean_counts)) %>% 
  mutate(mean_percs = mean_counts / sum(mean_counts))

## Version with sd/95%HPD

# Load melted df of event counts per maps for all bioregions and strata
# all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
names(all_events_count_per_maps_all_bioregions_all_strata_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Aggregate counts across maps for sources: mean, perc, sd, HPD95%
all_dispersal_events_summary_counts_for_sources_overall_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(source, map) %>% 
  summarize(counts = sum(counts)) %>% # Sum across event_types x strata x destinations
  group_by(source) %>% # Aggregate counts across maps using summary statistics
  summarize(mean_counts = mean(counts),
            sd_counts = sd(counts),
            HPD2.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[1],
            HPD97.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[2]) %>% 
  ungroup() %>% # Compute overall percentages
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts))


# Adjust color scheme for bioregions
areas_list <- c("A", "U", "I", "R", "N", "E", "W")
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")

# Adjust order of bioregions
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "E-PA", "W-PA")
all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot$source <- factor(all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot$source, levels = levels(all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot$source), labels = bioregion_names_reduced)
all_dispersal_events_summary_counts_for_sources_overall_ggplot$source <- factor(all_dispersal_events_summary_counts_for_sources_overall_ggplot$source, levels = levels(all_dispersal_events_summary_counts_for_sources_overall_ggplot$source), labels = bioregion_names_reduced)

# Save melted df for summary of counts of dispersal events across source bioregions
saveRDS(object = all_dispersal_events_summary_counts_for_sources_overall_ggplot, file = "outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_summary_counts_for_sources_overall_ggplot.rds")


## GGplot without sd
# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)

barplot_mean_dispersal_events_count_per_source_bioregions_ggplot <- ggplot(data = all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot,
                                                                           mapping = aes(y = mean_counts, x = source, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Dispersal events per Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Source Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Source Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_mean_dispersal_events_count_per_source_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_dispersal_events_count_per_source_bioregions_ggplot)

dev.off()


### 8.7/ Plot mean counts of Immigration events between (Destination) bioregions ####

# Load melted df of mean counts of dispersal events between bioregions
# all_dispersal_events_mean_count_per_maps_overall_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")
all_dispersal_events_mean_count_per_maps_overall_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_overall_ggplot.rds")

# Aggregate counts per destination bioregions
all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot <- all_dispersal_events_mean_count_per_maps_overall_ggplot %>% 
  group_by(dest) %>% 
  summarize(mean_counts = sum(mean_counts)) %>% 
  mutate(mean_percs = mean_counts / sum(mean_counts))

## Version with sd/95%HPD

# Load melted df of event counts per maps for all bioregions and strata
# all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
names(all_events_count_per_maps_all_bioregions_all_strata_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Aggregate counts across maps for destinations: mean, perc, sd, HPD95%
all_dispersal_events_summary_counts_for_dest_overall_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(dest, map) %>% 
  summarize(counts = sum(counts)) %>% # Sum across event_types x strata x sources
  group_by(dest) %>% # Aggregate counts across maps using summary statistics
  summarize(mean_counts = mean(counts),
            sd_counts = sd(counts),
            HPD2.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[1],
            HPD97.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[2]) %>% 
  ungroup() %>% # Compute overall percentages
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts))


# Adjust color scheme for bioregions
areas_list <- c("A", "U", "I", "R", "N", "E", "W")
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")

# Adjust order of bioregions
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "E-PA", "W-PA")
all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot$dest <- factor(all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot$dest, levels = levels(all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot$dest), labels = bioregion_names_reduced)
all_dispersal_events_summary_counts_for_dest_overall_ggplot$dest <- factor(all_dispersal_events_summary_counts_for_dest_overall_ggplot$dest, levels = levels(all_dispersal_events_summary_counts_for_dest_overall_ggplot$dest), labels = bioregion_names_reduced)

# Save melted df for summary of counts of dispersal events across destination bioregions
saveRDS(object = all_dispersal_events_summary_counts_for_dest_overall_ggplot, file = "outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_summary_counts_for_dest_overall_ggplot.rds")

## GGplot without sd
# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_per_dest_bioregions_overall_barplot.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events_count_per_dest_bioregions_overall_barplot.pdf", width = 10, height = 6)

barplot_mean_dispersal_events_count_per_dest_bioregions_ggplot <- ggplot(data = all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot,
                                                                           mapping = aes(y = mean_counts, x = dest, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("Dispersal events per Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Destination Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Destination Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_mean_dispersal_events_count_per_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_mean_dispersal_events_count_per_dest_bioregions_ggplot)

dev.off()


### 8.8/ Plot mean counts of Emigration/Immigration events between (Source/Destination) bioregions ####

# Function to annotate using npc units
annotate_npc <- function(label, x, y, ...)
{
  ggplot2::annotation_custom(grid::textGrob(
    x = unit(x, "npc"), y = unit(y, "npc"), label = label, ...))
}

## GGplot without sd
# pdf(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_per_source_dest_bioregions_overall_barplot.pdf", width = 10, height = 6)
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events_count_per_source_dest_bioregions_overall_barplot.pdf", width = 10, height = 6)

barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- ggplot(data = all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot) +
  
  # Plot barplot for destination bioregions as positive values (Emigration)
  geom_col(data = all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot,
           aes(y = mean_counts, x = dest, fill = dest), alpha = 1.0, col = "black") +
  
  # Plot barplot for source bioregions as negative values (Immigration)
  geom_col(data = all_dispersal_events_mean_count_per_source_bioregions_overall_ggplot,
           aes(y = -1 * mean_counts, x = source, fill = source), alpha = 1.0, col = "black") +
  
  # Add Bioregion labels on the top of the bars
  geom_text(data = all_dispersal_events_mean_count_per_dest_bioregions_overall_ggplot,
            aes(y = mean_counts, x = dest), label = bioregion_names_reduced,
            alpha = 1.0, nudge_y = 15, size = 6, col = "black") +
  
  # # Add text to identify Source/Dest
  # annotate_npc(label = "Destination (Emigration)",
  #              x = 0.75, y = 0.92,
  #              gp = grid::gpar(fontsize = 19, fontface = "bold", col = "black")) +
  # annotate_npc(label = "Source (Immigration)",
  #              x = 0.75, y = 0.10,
  #              gp = grid::gpar(fontsize = 19, fontface = "bold", col = "black")) +
  
  # Add text to identify Source/Dest
  geom_label(label = "Destination (Emigration)",
             x = 5.6, y = 170,
             label.padding = unit(0.5, "lines"),
             size = 6, fontface = "bold", col = "black") +
  geom_label(label = "Source (Immigration)",
             x = 5.6, y = -170,
             label.padding = unit(0.5, "lines"),
             size = 6, fontface = "bold", col = "black") +
  
  # Add horizontal line for y = 0
  geom_hline(yintercept = 0, linewidth = 1.0) +
  
  # Adjust tick labels of Y-axis
  scale_y_continuous(breaks = c(-200, -100, 0, 100, 200), 
                     labels = c(200, 100, 0, 100, 200),
                     limits = c(-200, 200)) +
  
  # Set plot title +
  ggtitle(label = paste0("Dispersal events per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Remove x-axis
barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# # Remove vertical grid
# barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())

# Plot
print(barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot)

dev.off()


## GGplot with sd
pdf(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events_count_per_source_dest_bioregions_overall_barplot_with_sd.pdf", width = 10, height = 6)

barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- ggplot(data = all_dispersal_events_summary_counts_for_sources_overall_ggplot) +
  
  # Plot barplot for destination bioregions as positive values (Emigration)
  geom_col(data = all_dispersal_events_summary_counts_for_dest_overall_ggplot,
           aes(y = mean_counts, x = dest, fill = dest), alpha = 1.0, col = "black") +
  
  # Plot barplot for source bioregions as negative values (Immigration)
  geom_col(data = all_dispersal_events_summary_counts_for_sources_overall_ggplot,
           aes(y = -1 * mean_counts, x = source, fill = source), alpha = 1.0, col = "black") +
  
  # Add SD of destination as error bars
  geom_errorbar(data = all_dispersal_events_summary_counts_for_dest_overall_ggplot,
                aes(ymin = mean_counts, ymax = mean_counts + sd_counts, x = dest),
                width = 0.15,
                position = "identity") +
  
  # Add SD of sources as error bars
  geom_errorbar(data = all_dispersal_events_summary_counts_for_sources_overall_ggplot,
                aes(ymin = -1 *(mean_counts), ymax = -1 *(mean_counts + sd_counts), x = source),
                width = 0.15,
                position = "identity") +
  
  # Add Bioregion labels on the top of the bars (+ error bars)
  geom_text(data = all_dispersal_events_summary_counts_for_dest_overall_ggplot,
            aes(y = mean_counts + sd_counts, x = dest), label = bioregion_names_reduced,
            alpha = 1.0, nudge_y = 15, size = 6, col = "black") +
  
  # # Add text to identify Source/Dest
  # annotate_npc(label = "Destination (Emigration)",
  #              x = 0.75, y = 0.92,
  #              gp = grid::gpar(fontsize = 19, fontface = "bold", col = "black")) +
  # annotate_npc(label = "Source (Immigration)",
  #              x = 0.75, y = 0.10,
  #              gp = grid::gpar(fontsize = 19, fontface = "bold", col = "black")) +
  
  # Add text to identify Source/Dest
  geom_label(label = "Destination (Emigration)",
             x = 5.6, y = 170,
             label.padding = unit(0.5, "lines"),
             size = 6, fontface = "bold", col = "black") +
  geom_label(label = "Source (Immigration)",
             x = 5.6, y = -170,
             label.padding = unit(0.5, "lines"),
             size = 6, fontface = "bold", col = "black") +
  
  # Add horizontal line for y = 0
  geom_hline(yintercept = 0, linewidth = 1.0) +
  
  # Adjust tick labels of Y-axis
  scale_y_continuous(breaks = c(-200, -100, 0, 100, 200), 
                     labels = c(200, 100, 0, 100, 200),
                     limits = c(-200, 200)) +
  
  # Set plot title +
  ggtitle(label = paste0("Dispersal events per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Remove x-axis
barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# # Remove vertical grid
# barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot <- barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())

# Plot
print(barplot_mean_dispersal_events_count_per_source_dest_bioregions_ggplot)

dev.off()



##### 9/ Plot Emigration/Immigration counts between bioregions, per time-strata #####

### 9.1/ Aggregate counts for all dispersal events, per time-strata ####

# Load ggplot df
# all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")
all_events_count_per_maps_all_bioregions_all_strata_ggplot <- readRDS(file = "./outputs/Discrete_events_counts/Ponerinae_MCC_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_bioregions_all_strata_ggplot.rds")

names(all_events_count_per_maps_all_bioregions_all_strata_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Sum across event_types, per time-strata
# Aggregate mean counts across maps
all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot <- all_events_count_per_maps_all_bioregions_all_strata_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(source, dest, map, stratum) %>% 
  summarize(counts = sum(counts)) %>% # Sum across event_types
  group_by(source, dest, stratum) %>%
  summarize(mean_counts = mean(counts)) %>% # Aggregate mean counts across maps
  ungroup() %>% # Compute overall percentages
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts))

# Save ggplot df of mean counts of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot.rds")
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot.rds")

## Convert back to array

# Pivot mean counts and percentages
all_dispersal_events_mean_count_per_maps_all_time_strata_melted <- all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot %>% 
  pivot_longer(cols = c("mean_counts", "mean_perc"), names_to = "mean_stats", values_to = "values")

all_dispersal_events_mean_count_per_maps_all_time_strata_array <- reshape2::acast(data = all_dispersal_events_mean_count_per_maps_all_time_strata_melted,
                                                                          formula = source ~ dest ~ stratum ~ mean_stats)

dimnames(all_dispersal_events_mean_count_per_maps_all_time_strata_array)
all_dispersal_events_mean_count_per_maps_all_time_strata_array[, , , "mean_counts"]
all_dispersal_events_mean_count_per_maps_all_time_strata_array[, , , "mean_perc"]

# Save array of mean counts/perc of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_time_strata_array, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_array.rds")
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_time_strata_array, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_array.rds")

### 9.2/ Compute node metadata per time strata ####

# Load geological epoch time boundaries
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")

time_boundaries <- c(0, DEC_J_fit$inputs$timeperiods)

# Load LTT data
# DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot.rds")
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot.rds")

## Create table for bioregion coordinates
bioregion_names_alpha_order <- c("Afrotropics", "Eastern Palearctic", "Indomalaya", "Neotropics", "Nearctic", "Australasia", "Western Palearctic")
areas_list <- levels(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$areas)
areas_list <- areas_list[areas_list != "total"]
areas_alpha_order <- areas_list[order(areas_list)]
bioregion_names_alpha_order_df <- data.frame(areas = areas_alpha_order, bioregions = bioregion_names_alpha_order)

## Create table for bioregion coordinates (May be adjusted per time strata to account for continental drift)
bioregion_coordinates_df <- bioregion_names_alpha_order_df
bioregion_coordinates_df$latitude <- c(0, 50, 10, -10, 40, -25, 47)
bioregion_coordinates_df$longitude <- c(25, 100, 100, -60, -100, 140, 20)
bioregion_coordinates_df$adjusted_latitude <- c(10, 35, 10, -10, 35, -15, 40)
bioregion_coordinates_df$adjusted_longitude <- c(15, 120, 100, -60, -100, 150, 15)

## Compute mean/max species richness per geological epoch
nodes_metadata_per_time_strata <- list()
for (i in 1:(length(time_boundaries) - 1))
{
  # i <- 1
  
  # Extract time boundaries
  lower_boundary <- time_boundaries[i]
  upper_boundary <- time_boundaries[i+1]
  
  # Aggregate species richness within each time stratum
  nodes_metadata_i <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot %>% 
    filter((smoothed_time >= lower_boundary) & (smoothed_time <= upper_boundary)) %>% 
    filter(areas != "total") %>% 
    group_by(areas) %>%
    summarize(max_counts = max(mean_counts),  # May make more sense to plot the final richness of each stratum
              mean_counts = mean(mean_counts)) %>%  
    ungroup() %>%
    mutate(mean_percentages = 100 * mean_counts / sum(mean_counts),
           max_percentages = 100 * max_counts / sum(max_counts)) %>%
    mutate(bioregions = areas) %>%
    rename(node_ID = areas) %>%
    mutate(node_size = log(mean_counts),
           node_size_max = log(max_counts)) %>%
    mutate(stratum = paste0("Stratum_",i)) %>%
    dplyr::select(stratum, node_ID, mean_counts, mean_percentages, node_size, max_counts, max_percentages, node_size_max)
  
  # Add bioregion labels
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_names_alpha_order_df, by = join_by(node_ID == areas)) %>%
    dplyr::select(stratum, node_ID, bioregions, mean_counts, mean_percentages, node_size, max_counts, max_percentages, node_size_max)
  
  # Add label with richness
  nodes_metadata_i$node_labels <- paste0(nodes_metadata_i$node_ID, "\n", round(nodes_metadata_i$mean_counts, 0))
  nodes_metadata_i$node_labels_max <- paste0(nodes_metadata_i$node_ID, "\n", round(nodes_metadata_i$max_counts, 0))
  
  # Reorder in alphabetic order of node ID
  nodes_metadata_i <- nodes_metadata_i %>% 
   arrange(node_ID)
  
  ## Add spatial coordinates for each bioregion (May adjust to show continental drift)
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_coordinates_df,
                                by = join_by(node_ID == areas, bioregions == bioregions))

  # Store in final list
  nodes_metadata_per_time_strata <- append(x = nodes_metadata_per_time_strata, values = list(nodes_metadata_i))
}
names(nodes_metadata_per_time_strata) <- paste0("Stratum_",1:(length(time_boundaries) - 1))

## Update max counts for Stratum 1 (PPHo) because it is taken for time = 2My due to effect of sliding windows, and not time = 0My
nodes_metadata_per_time_strata$Stratum_1$max_counts <- nodes_metadata$mean_counts[match(x = nodes_metadata$bioregions, table = nodes_metadata_per_time_strata$Stratum_1$bioregions)]
nodes_metadata_per_time_strata$Stratum_1$node_labels_max <- paste0(nodes_metadata_per_time_strata$Stratum_1$node_ID, "\n", round(nodes_metadata_per_time_strata$Stratum_1$max_counts, 0))

## Save node metadata
# saveRDS(object = nodes_metadata_per_time_strata, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_time_strata.rds")
saveRDS(object = nodes_metadata_per_time_strata, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata.rds")


### 9.3/ Plot each time stratum as a network ####

# https://kateto.net/network-visualization

time_strata_full_names <- c("PPHo (5.33-0 My)", "Miocene (23.03-5.33 My)", "Oligocene (33.9-23.03 My)", "Eocene (56.0-33.9 My)", "Paleocene (66.0-56.0 My)", "Late Cretaceous (100.5-66.0 My)", "Early Cretaceous (145.0-100.5 My)")
all_dispersal_events_all_strata_df_list <- list()

# Load node metadata for all time strata
# nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_time_strata.rds")
nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata.rds")

# Load ggplot df of mean counts of dispersal events between bioregions per time strata
# all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot.rds")
all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot.rds")

# Adjust minimal number of events to be displayed for each time strata
min_counts_threshold_list <- c(5, 3, 1, 1, 1, 1, 1)

## Loop per time stratum
for (i in seq_along(nodes_metadata_per_time_strata))
{
  # i <- 1
  
  # Extract time stratum/Geological epoch names
  time_stratum_full_name <- time_strata_full_names[i]
  
  # Extract node metadata for the given time stratum
  nodes_metadata_stratum_i <- nodes_metadata_per_time_strata[[i]] %>%
    dplyr::select(-stratum)
  
  ## 9.3.1/ Convert to igraph object ####
  
  # Extract vertice metadata for the given time stratum
  # Round counts and convert to edge size
  all_dispersal_events_stratum_i_df <- all_dispersal_events_mean_count_per_maps_all_time_strata_ggplot %>% 
    filter(stratum == paste0("Stratum_",i)) %>% 
    mutate(mean_counts = round(mean_counts, 0)) %>% 
    mutate(mean_perc = round(100 * mean_counts / sum(mean_counts), 1)) %>% 
    mutate(edge_width = log1p(mean_counts)) # Adjust edge with to a log scale
    # mutate(edge_width = mean_counts / max(mean_counts)) # Standardize between 0 and 1
  
  ## Store egde metadata across all time-strata
  all_dispersal_events_all_strata_df_list <- append(x = all_dispersal_events_all_strata_df_list, values = list(all_dispersal_events_stratum_i_df))
  
  ## Convert to igraph
  all_dispersal_events_stratum_i_igraph <- igraph::graph_from_data_frame(d = all_dispersal_events_stratum_i_df, vertices = nodes_metadata_stratum_i, directed = T)
  # str(all_dispersal_events_stratum_i_igraph)
  
  # Explore igraph
  # all_dispersal_events_stratum_i_igraph
  E(all_dispersal_events_stratum_i_igraph) # Edge attributes
  V(all_dispersal_events_stratum_i_igraph) # Vertex/Node attributes
  all_dispersal_events_stratum_i_igraph[]  # Adjacency matrix of nodes x nodes
  
  # Export adjacency matrix
  as_adjacency_matrix(all_dispersal_events_stratum_i_igraph, attr = "mean_counts")
  as_adjacency_matrix(all_dispersal_events_stratum_i_igraph, attr = "mean_perc")
  
  ## 9.3.2/ Adjust aesthetics ####
  
  # Set color scheme for areas/bioregions (Use the BSM color scheme)
  colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
  # colors_list_for_areas <- colors_list_for_states[nodes_metadata_stratum_i$node_ID]
  colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
  colors_list_for_areas <- colors_list_for_areas[nodes_metadata_stratum_i$node_ID]
  
  # Remove loops
  all_dispersal_events_stratum_i_igraph <- simplify(graph = all_dispersal_events_stratum_i_igraph,
                                                  remove.multiple = F, remove.loops = T) 
  
  # Set node colors based on bioregions:
  V(all_dispersal_events_stratum_i_igraph)$color <- colors_list_for_areas
  # Set node size based on current richness
  # V(all_dispersal_events_stratum_i_igraph)$size <- V(all_dispersal_events_stratum_i_igraph)$node_size * 5 # Mean richness
  V(all_dispersal_events_stratum_i_igraph)$size <- V(all_dispersal_events_stratum_i_igraph)$node_size_max * 5 # Max richness
  V(all_dispersal_events_stratum_i_igraph)$size <- sapply(X = V(all_dispersal_events_stratum_i_igraph)$size, FUN = function (x) { max(x, 10) })
  # Set node font size based on current richness
  # V(all_dispersal_events_stratum_i_igraph)$label.cex <- V(all_dispersal_events_stratum_i_igraph)$node_size / 5 # Mean richness
  V(all_dispersal_events_stratum_i_igraph)$label.cex <- V(all_dispersal_events_stratum_i_igraph)$node_size_max / 5 # Max richness
  V(all_dispersal_events_stratum_i_igraph)$label.cex <- sapply(X = V(all_dispersal_events_stratum_i_igraph)$label.cex, FUN = function (x) { max(x, 0.5) })
    
  # Set edge width based on event counts
  E(all_dispersal_events_stratum_i_igraph)$width <- E(all_dispersal_events_stratum_i_igraph)$edge_width * 2 + 1
  # Set edge arrow size based on event counts
  E(all_dispersal_events_stratum_i_igraph)$arrow.size <- E(all_dispersal_events_stratum_i_igraph)$edge_width / 4 + 0.1
  
  # Set edge label as mean counts
  # E(all_dispersal_events_stratum_i_igraph)$edge.label <- E(all_dispersal_events_stratum_i_igraph)$mean_counts
  
  # Set edge color according to their source
  edge_starts <- ends(all_dispersal_events_stratum_i_igraph, es = E(all_dispersal_events_stratum_i_igraph), names = F)[, 1]
  edge_col <- V(all_dispersal_events_stratum_i_igraph)$color[edge_starts]
  E(all_dispersal_events_stratum_i_igraph)$color <- edge_col
  
  # Remove edges with low counts
  table(E(all_dispersal_events_stratum_i_igraph)$mean_counts)
  # cutoff_min_counts_i <- cutoff_min_counts
  cutoff_min_counts_i <- min_counts_threshold_list[i]
  all_dispersal_events_stratum_i_igraph <- delete_edges(graph = all_dispersal_events_stratum_i_igraph,
                                                      edges = E(all_dispersal_events_stratum_i_igraph)[mean_counts < cutoff_min_counts_i])
  
  # Adjust layout to coordinates of vertex/nodes
  layout_WGS84 <- cbind(V(all_dispersal_events_stratum_i_igraph)$longitude, V(all_dispersal_events_stratum_i_igraph)$latitude)
  layout_pretty <- cbind(V(all_dispersal_events_stratum_i_igraph)$adjusted_longitude, V(all_dispersal_events_stratum_i_igraph)$adjusted_latitude)
  
  ## Save igraph for all dispersal events across all time-strata
  # saveRDS(object = all_dispersal_events_stratum_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_igraph.rds"))
  saveRDS(object = all_dispersal_events_stratum_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_igraph.rds"))
  
  ## Plot igraph
  
  # pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_igraph.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_igraph.pdf"),
      width = 6, height = 6)
  
  plot.igraph(x = all_dispersal_events_stratum_i_igraph,
              ## Node aesthetics
              # vertex.color = colors_list_for_areas, # Nodes color
              vertex.frame.color = "black", # Nodes border color
              vertex.shape = "circle",  # Node shape
              # vertex.size = 15, # Node size
              # vertex.label	= nodes_metadata_stratum_i$node_labels, # Node labels based on mean richness
              vertex.label	= nodes_metadata_stratum_i$node_labels_max, # Node labels based on max/final richness
              vertex.label.color = "black",  # Node label color
              vertex.label.font = 2, # Node label font
              # vertex.label.cex = 1, # Node label cex
              ## Edge aesthetics
              # edge.color = "black",     # Edge color
              # edge.width = 1,         # Edge width
              # edge.arrow.size = 0.7,    # Terminal arrow size
              edge.label = E(all_dispersal_events_stratum_i_igraph)$mean_counts,
              edge.label.cex = 0.7,
              edge.lty = 1,             # Edge line type
              edge.curved = 0.1,
              arrow.mode = "forward",    # Direction of arrows
              ## Other aesthetics
              # layout = layout_WGS84,
              layout = layout_pretty,
              main = paste0("Dispersal events between bioregions\n",time_stratum_full_name)
  )
  
  dev.off()
  
}

## Save egde metadata across all time-strata
names(all_dispersal_events_all_strata_df_list) <- paste0("Stratum_",1:length(all_dispersal_events_all_strata_df_list))
# saveRDS(object = all_dispersal_events_all_strata_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))
saveRDS(object = all_dispersal_events_all_strata_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))


## 9.3.3/ Aggregate all igraph plots in a single pdf ####

# all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_time_strata/", pattern = "all_dispersal_events_stratum_", full.names = T)
all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata/", pattern = "all_dispersal_events_stratum_", full.names = T)

all_igraphs_path <- all_igraphs_path[str_detect(string = all_igraphs_path, pattern = "igraph.pdf")]
nb_igraphs <- length(all_igraphs_path)

# qpdf::pdf_combine(input = all_igraphs_path, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_igraph.pdf"))
qpdf::pdf_combine(input = all_igraphs_path, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_igraph.pdf"))

## 9.3.4/ Aggregate all igraph in forward timeline + overall ####

# overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"
overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"

# qpdf::pdf_combine(input = c(rev(all_igraphs_path), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_igraph_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_igraphs_path), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_igraph_forward.pdf"))


### 9.4/ Map networks over bioregion maps ####

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

## Load overall node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

## Load node metadata per time strata
# nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_time_strata.rds")
nodes_metadata_per_time_strata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_time_strata.rds")

## Load overall edge metadata
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")

# Load edge metadata
# all_dispersal_events_all_strata_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))
all_dispersal_events_all_strata_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds"))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[nodes_metadata_per_time_strata[[1]]$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[nodes_metadata_per_time_strata[[1]]$node_ID]
names(colors_list_for_areas) <- nodes_metadata_per_time_strata[[1]]$bioregions
bioregion_names <- nodes_metadata_per_time_strata[[1]]$bioregions

# Adjust order of bioregions
for (i in 1:length(nodes_metadata_per_time_strata))
{
  nodes_metadata_per_time_strata[[i]]$bioregions <- factor(nodes_metadata_per_time_strata[[i]]$bioregions, levels = bioregion_names, labels = bioregion_names)
}
Bioregions_sf_Bioregions_level$Bioregion <- factor(Bioregions_sf_Bioregions_level$Bioregion, levels = bioregion_names, labels = bioregion_names)

# Convert Long/Lat WGS84 to Mollweide
for (i in 1:length(nodes_metadata_per_time_strata))
{
  nodes_metadata_sf <- st_as_sf(x = nodes_metadata_per_time_strata[[i]], coords = c("longitude", "latitude"), remove = F, crs = st_crs(Bioregions_sf_Bioregions_level))
  nodes_metadata_sf <- st_transform(nodes_metadata_sf, crs = st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  Mollweide_coordinates <- st_as_text(nodes_metadata_sf$geometry)
  longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  nodes_metadata_per_time_strata[[i]]$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
  latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
  nodes_metadata_per_time_strata[[i]]$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))
}

## 9.4.1/ Inform edge metadata with coordinates and bioregion source labels ####

for (i in 1:length(all_dispersal_events_all_strata_df_list))
{
  all_dispersal_events_strata_i_df <- all_dispersal_events_all_strata_df_list[[i]]
  nodes_metadata_i <- nodes_metadata_per_time_strata[[i]]
  
  all_dispersal_events_strata_i_df <- left_join(x = all_dispersal_events_strata_i_df,
                                               y = nodes_metadata_i[, c("node_ID", "bioregions", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                               by = join_by(source == node_ID))
  all_dispersal_events_strata_i_df <- all_dispersal_events_strata_i_df %>% 
    rename(source_labels = bioregions,
           source_latitude = latitude,
           source_longitude = longitude,
           source_latitude_Mollweide = latitude_Mollweide,
           source_longitude_Mollweide = longitude_Mollweide)
  all_dispersal_events_strata_i_df <- left_join(x = all_dispersal_events_strata_i_df,
                                               y = nodes_metadata_i[, c("node_ID", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                               by = join_by(dest == node_ID))
  all_dispersal_events_strata_i_df <- all_dispersal_events_strata_i_df %>% 
    rename(dest_latitude = latitude,
           dest_longitude = longitude,
           dest_latitude_Mollweide = latitude_Mollweide,
           dest_longitude_Mollweide = longitude_Mollweide)
  
  # Adjust order of bioregions
  all_dispersal_events_strata_i_df$source_labels <- factor(all_dispersal_events_strata_i_df$source_labels, levels = bioregion_names, labels = bioregion_names)
  
  # Store updated node metadata df
  all_dispersal_events_all_strata_df_list[[i]] <- all_dispersal_events_strata_i_df
}

# Save edge metadata df
# saveRDS(object = all_dispersal_events_all_strata_df_list, file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds")
saveRDS(object = all_dispersal_events_all_strata_df_list, file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds")


# Load edge metadata df
# all_dispersal_events_all_strata_df_list <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds")
all_dispersal_events_all_strata_df_list <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_df_list.rds")


# May adjust curvature manually
# May adjust arrow size too

# Set minimum nb of dispersal events to display
# min_counts_threshold <- 1
min_counts_threshold_list <- c(5, 3, 1, 1, 1, 1, 1)

max_node_counts <- max(nodes_metadata$mean_counts)
max_edge_counts <- max(all_dispersal_events_overall_df$mean_counts)

# Set breaks for edge linewidth
breaks_linewidth_list <- list(Stratum_1 = c(6, 20, 45), Stratum_2 = c(5, 15, 30), Stratum_3 = c(1, 2, 3), Stratum_4 = c(1, 2, 5), Stratum_5 = c(1, 2), Stratum_6 = c(1, 2), Stratum_7 = c(1))

# Function to adjust range size of nods based on the range scale used for overall data
rescale_range_size <- function(x, data_min = 0, data_max = max_node_counts, range_min = 5, range_max = 30) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}
  
# Function to adjust range linewidth of edges based on the range scale used for overall data
rescale_range_linewidth <- function(x, data_min = 5, data_max = max_edge_counts, range_min = 1, range_max = 6) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}


## 9.4.2/ Plot ggplot for each time stratum ####

for (i in seq_along(all_dispersal_events_all_strata_df_list))
{
  # i <- 1
  # i <- 3
  
  # Extract node metadata
  nodes_metadata_strata_i <- nodes_metadata_per_time_strata[[i]]
  # Extract max node size to adjust range of points
  # max_node_size <- max(nodes_metadata_strata_i$mean_counts) # Mean richness
  max_node_size <- max(nodes_metadata_strata_i$max_counts) # Max richness
  max_range_size <- rescale_range_size(x = max_node_size)
  # Extract min node size to adjust range of points
  # min_node_size <- min(nodes_metadata_strata_i$mean_counts) # Mean richness
  min_node_size <- min(nodes_metadata_strata_i$max_counts) # Max richness
  min_range_size <- rescale_range_size(x = min_node_size)
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract edge metadata
  all_dispersal_events_strata_i_df <- all_dispersal_events_all_strata_df_list[[i]]
  # Extract max node size to adjust range of points
  max_edge_linewidth <- max(all_dispersal_events_strata_i_df$mean_counts)
  max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  # Extract min node size to adjust range of points
  min_edge_linewidth <- min(all_dispersal_events_strata_i_df$mean_counts)
  min_range_linewidth <- rescale_range_linewidth(x = min_edge_linewidth)
  # Extract manual breaks for edge linewidth per stratum
  breaks_linewidth_i <- breaks_linewidth_list[[i]]
  
  # Extract time stratum/Geological epoch names
  time_stratum_full_name <- time_strata_full_names[i]
  
  # Extract the minimum number of events to display an edge
  # min_counts_threshold_i <- min_counts_threshold
  min_counts_threshold_i <- min_counts_threshold_list[i]
  
  # Plot PDF
  #pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_ggplot.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata/all_dispersal_events_stratum_",i,"_ggplot.pdf"),
      width = 9, height = 4)
  
  all_dispersal_events_strata_i_ggplot <- ggplot(data = nodes_metadata_strata_i) +
    
    # Plot bioregion maps
    geom_sf(data = Bioregions_sf_Bioregions_level,
            mapping = aes(fill = Bioregion),
            colour = "black",
            alpha = 0.8) +
    
    # Adjust fill color scheme for bioregions
    scale_fill_manual("Bioregions", breaks = bioregion_names, labels = bioregion_names, values = colors_list_for_areas) +
    
    # Adjust CRS
    # coord_sf(default_crs = sf::st_crs(4326)) +
    coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
             clip = "off", # To allow plotting arrow outside of World map
             expand = FALSE) +
    
    # Plot nodes
    geom_point(mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide,
                             # size = mean_counts), # Mean richness
                             size = max_counts), # Max richness
               alpha = 0.3, show.legend = F) +
    
    # Adjust size legend
    scale_size_continuous("Species richness", 
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +
    
    # Plot vertices
    geom_curve(data = all_dispersal_events_strata_i_df[all_dispersal_events_strata_i_df$mean_counts >= min_counts_threshold_i, ],
               aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                   xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                   color = source_labels, linewidth = mean_counts),
               arrow = arrow(length = unit(0.03, "npc"), 
                             type = "closed"), # Describes arrow head (open or closed)
               angle = 90, # Anything other than 90 or 0 can look unusual
               alpha = 1.0
    ) + 
    
    # Adjust color scheme for arrows
    scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
    
    # Adjust edge width legend
    scale_linewidth_continuous("Dispersal events",
                               breaks = breaks_linewidth_i,
                               # range = c(5, 30),
                               range = c(min_range_linewidth, max_range_linewidth)) +
    
    # Add node labels
    ggnewscale::new_scale(new_aes = "size") +
    geom_text(data = nodes_metadata_strata_i,
              mapping = aes(# label = round(mean_counts, 0), # Mean richness
                            label = round(max_counts, 0), # Max richness
                            x = longitude_Mollweide, y = latitude_Mollweide,
                            # size = mean_counts), # Mean richness
                            size = max_counts), # Max richness
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +
    
    # Adjust size for labels
    scale_size_continuous("Species richness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
    
    # Add title
    ggtitle(label =  paste0("Dispersal events between bioregions\n",time_stratum_full_name)) +
    
    # Adjust legend aesthetics
    guides(color = "none",
           size = "none",
           fill = guide_legend(order = 1),
           linewidth = guide_legend(order = 2,
                                    override.aes = list(arrow(length = unit(0.2, "cm"))))
    ) +
    
    # Adjust aesthetics
    theme_classic() +
    
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(fill = NA, colour = NA),
          # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          legend.box.margin = margin(l = 15, b = 15),
          # axis.ticks = element_line(linewidth = 1.0),
          # axis.ticks.length = unit(10, "pt"),
          # axis.text = element_text(size = 21, color = "black", face = "bold"),
          # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
          # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title = element_blank())
  
  print(all_dispersal_events_strata_i_ggplot)
  
  dev.off()
  
}

## 9.4.3/ Aggregate all ggplot plots in a single pdf ####

# all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_time_strata/", pattern = "all_dispersal_events_stratum_", full.names = T)
all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_time_strata/", pattern = "all_dispersal_events_stratum_", full.names = T)

all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = "ggplot.pdf")]
nb_ggplots <- length(all_ggplots_path)

# qpdf::pdf_combine(input = all_ggplots_path, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot.pdf"))
qpdf::pdf_combine(input = all_ggplots_path, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot.pdf"))

## 9.4.3/ Aggregate all ggplot in forward timeline + overall ####

# overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"
overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"

# qpdf::pdf_combine(input = c(rev(all_ggplots_path), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_ggplots_path), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"))


### 9.5/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

# pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"),
pdf_pointer <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.pdf"),
                                      pages = NULL, density = 100)
magick::image_info(pdf_pointer)

image_resize_and_write_gif(image = pdf_pointer,
                           # path =  paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.gif"),
                           path =  paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_strata_ggplot_forward.gif"),
                           delay = 2, # Time between frames in seconds
                           width = 900, height = 400,
                           loop = FALSE,
                           progress = TRUE)


## Add arrow to legend (Does not work with guides. Need to add on PPT)
## Add bg color on sphere (does not clip out of the sphere... Looks better without grid and bg)

# Add Color contour in igraphs on areas that are adjacent at a given time: Section 4.3 and 4.4
# Adjust coordinates in each stratum to show continental drift?

# Create bioregion shapefile that can overlay the PALEOMAPS?

### Make an animation by plotting every geological epoch?
# Use gganimate for paths to make them fade in_out?



