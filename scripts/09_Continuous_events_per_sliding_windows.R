##### Script 09: Continuous events per sliding windows #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Aggregate event counts per sliding windows
   # Compute number/probabilities of each class of events (mean and sd), per sliding windows
   # Plot the number of each class of events, per sliding windows
      # Cum raw = lines/stacked bars
      # Raw = lines/stacked bars
      # Percentages = stacked bars
   # Plot the rates of event types, per sliding windows = divide per nb of lineages and width of sliding window (time)

# Plot Emigration/Immigration counts between bioregions through time (per sliding windows)
   # Raw = lines/stacked bars 
   # Percentages = stacked bars
# Plot rates of Emigration/Immigration between bioregions through time (per sliding windows) = divide per nb of lineages and width of sliding window (time)
   # Get mean and CI of rates from all BSM maps

## Type of events
  # Per type of transitions
     # Anagenetic:
        # Range extension (d)
        # Range extinction (e)
     # Cladogenetic
        # Inheritence (y)
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
  # Tables of number/probabilities of each class of events per sliding windows
  # Summary statistics (mean, sd, 95% HPD) of each class of events per sliding windows

## Event type plots
  # Stacked bars of cumulative mean number of event types along sliding windows (raw and percentage)
  # Stacked bars of mean number/perc of event types along sliding windows (raw and percentage)
  # Rates of event types along sliding windows = divide per nb of lineages and width of sliding window (time)
      # Overall vs. per bioregions (source/destination)

## Dispersal plots
  # Emigration/Immigration counts between bioregions
      # 1 network map per sliding window

###


####### May wish to focus on specific clades or bioregions to compute the same series of plots #######

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
library(BayesTwin) # To compute HPD intervals
library(igraph)
library(sf)
library(magick)   # For animated GIF


# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load modeling results from best model (DEC+J) #####

# Load time-stratified DEC+J model output
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Youngest_phylogeny_1534t/DEC_J_fit.rds")
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_Oldest_phylogeny_1534t/DEC_J_fit.rds")

# Load BSM outputs
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_output.rds")
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_output.rds")
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_output.rds")
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_output.rds")

# Extract records of events
DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
DEC_J_ana_events_tables <- DEC_J_BSM_output$RES_ana_events_tables

# DEC_J_clado_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
# DEC_J_ana_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")
DEC_J_clado_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
DEC_J_ana_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")
# DEC_J_clado_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
# DEC_J_ana_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")
# DEC_J_clado_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_clado_events_tables_with_unique_source.rds")
# DEC_J_ana_events_tables_with_unique_source <- readRDS("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_ana_events_tables_with_unique_source.rds")

# Extract areas list
returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
all_areas <- returned_mats$areanames

### 1.3/  Define the sliding windows ####

# Set sliding window characteristics
window_width <- 10 # Width = 10 My (Must be a multiple of window_steps)  # For barplots
window_width <- 5 # Width = 5 My (Must be a multiple of window_steps)    # For networks
window_width <- 1 # Width = 1 My (To use for cumulative counts!)         # For cumulative counts
window_steps <- 1 # Steps = 1 My (Must be a divider of window_width)

# Extract root age
root_age <- max(DEC_J_clado_events_tables[[1]]$time_bp)

# Generate start and end times of sliding windows
start_time_list <- seq(from = 0, to = root_age - window_width, by = window_steps)
end_time_list <- seq(from = window_width, to = root_age, by = window_steps)
mean_time_list <- (start_time_list + end_time_list) / 2

### 1.4/ Function for ggplot aesthetics ####

# Aesthetics for density curve plots
add_aesthetics_barplot <- function (ggplot)
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


##### 2/ Compute residence times per areas along sliding windows #####

source(file = "./functions/compute_residence_times_per_sliding_windows.R")

### 2.1/ Compute residence times x ranges per window time for all simmaps ####

# Load simmaps of BS maps
# DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps.rds")
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps.rds")
# DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_simmaps.rds")
# DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_simmaps.rds")

# Set time boundaries of sliding windows
time_boundaries_df <- data.frame(tipward_times = start_time_list,   # Tipward time boundaries
                                 rootward_times = end_time_list)  # Rootward time boundaries

residence_times_per_sliding_windows_in_ranges_array <- compute_residence_times_per_sliding_windows_on_MultiSimmap(multiSimmap = DEC_J_simmaps,
                                                                                                                  time_boundaries = time_boundaries_df,  # Data.frame of start/end time of sliding windows)
                                                                                                                  remove_empty_states = T, # Should states with no residence times be removed?
                                                                                                                  melted = F, # Should the final array be melted into a df?
                                                                                                                  verbose = T)  


### 2.2/ Aggregate residence times per unique areas instead of ranges ####

# Aggregate residence times per unique areas instead of ranges
# Time for ranges encompassing multiple areas is equally split across those areas

residence_times_per_sliding_windows_in_areas_array <- aggregate_residence_times_per_sliding_windows_in_unique_areas(residence_times_per_sliding_windows_in_ranges_array,
                                                                                                                    remove_empty_areas = T, # Should areas with no residence times be removed?
                                                                                                                    melted = F)
## Melt into data.frame
residence_times_per_sliding_windows_in_areas_melted_df <- reshape2::melt(residence_times_per_sliding_windows_in_areas_array)
names(residence_times_per_sliding_windows_in_areas_melted_df) <- c("area", "window_label", "map", "residence_time")

# Add window ID and mean (smoothed) time
time_boundaries_df$window_ID <- 1:nrow(time_boundaries_df)
time_boundaries_df$window_label <- paste0(time_boundaries_df$rootward_times, "_", time_boundaries_df$tipward_times) 
time_boundaries_df$smoothed_time <- (time_boundaries_df$tipward_times + time_boundaries_df$rootward_times) / 2

residence_times_per_sliding_windows_in_areas_melted_df <- left_join(x = residence_times_per_sliding_windows_in_areas_melted_df,
                                                                    y = time_boundaries_df,
                                                                    join_by(window_label == window_label))
  
## Save final array and melted data.frame
# saveRDS(object = residence_times_per_sliding_windows_in_areas_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_array_", window_width,"My.rds"))
# saveRDS(object = residence_times_per_sliding_windows_in_areas_melted_df, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
saveRDS(object = residence_times_per_sliding_windows_in_areas_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_array_", window_width,"My.rds"))
saveRDS(object = residence_times_per_sliding_windows_in_areas_melted_df, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# saveRDS(object = residence_times_per_sliding_windows_in_areas_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_array_", window_width,"My.rds"))
# saveRDS(object = residence_times_per_sliding_windows_in_areas_melted_df, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# saveRDS(object = residence_times_per_sliding_windows_in_areas_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_array_", window_width,"My.rds"))
# saveRDS(object = residence_times_per_sliding_windows_in_areas_melted_df, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))


##### 3/ Summarize events from stochastic maps per sliding time-window intervals to obtain continuous rates ####

## Can extract count for any time-frame, not the just the time-periods defined for the model !!!

### 3.1/ Function to reorganize event counts lists as arrays for clarity ####

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




### 3.3/  Make a loop to count events per sliding windows ####

# Extract the nb of sliding windows
nb_sliding_windows <- length(start_time_list)

# Extract the nb of maps
nb_maps <- length(DEC_J_clado_events_tables_with_unique_source)

# Initiate list of final arrays
DEC_J_BSM_counts_arrays_all_sliding_windows <- list()

# Loop per sliding windows
for (i in 1:nb_sliding_windows)
{
  # i <- 2
  
  # Define time period boundary
  start_time <- start_time_list[i]
  end_time <- end_time_list[i]
  
  # Count all anagenetic and cladogenetic events for one sliding window
  DEC_J_BSM_counts_list_per_sliding_window_i <- count_ana_clado_events(clado_events_tables = DEC_J_clado_events_tables_with_unique_source,
                                                                       ana_events_tables = DEC_J_ana_events_tables_with_unique_source,
                                                                       areanames = all_areas,
                                                                       actual_names = all_areas,
                                                                       timeperiod = c(start_time, end_time))
  
  str(DEC_J_BSM_counts_list_per_sliding_window_i, max.level = 1)
  
  # Convert in array
  DEC_J_BSM_counts_arrays_per_sliding_window_i <- convert_BSM_counts_lists_to_arrays(DEC_J_BSM_counts_list_per_sliding_window_i)
  
  ## Initialize final lists/arrays once
  if (i == 1)
  {
    # Unique subset events counts
    unique_sub_counts_all_sliding_windows <- list()
    # Unique vicariance events counts
    unique_vic_counts_all_sliding_windows <- list()
    # Unique range inheritance (sympatry) events counts
    unique_sym_counts_all_sliding_windows <- list()
    # Source/dest dispersal events count matrices
    dispersal_events_count_matrices_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$dispersal_events_count_matrices
    dispersal_events_count_matrices_all_sliding_windows <- array(data = NA,
                                                                 dim = c(dim(dispersal_events_count_matrices_sliding_window_i), nb_sliding_windows),
                                                                 dimnames = append(x = dimnames(dispersal_events_count_matrices_sliding_window_i), values = list(paste0("Stratum_", 1:nb_sliding_windows))))
    # Source/dest dispersal events summary matrices
    dispersal_events_summary_matrices_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$dispersal_events_summary_matrices
    dispersal_events_summary_matrices_all_sliding_windows <- array(data = NA,
                                                                   dim = c(dim(dispersal_events_summary_matrices_sliding_window_i), nb_sliding_windows),
                                                                   dimnames = append(x = dimnames(dispersal_events_summary_matrices_sliding_window_i), values = list(paste0("Stratum_", 1:nb_sliding_windows))))
    # Range extinction counts x states
    range_contraction_events_count_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$range_contraction_events_count
    range_contraction_events_count_all_sliding_windows <- array(data = NA,
                                                                dim = c(dim(range_contraction_events_count_sliding_window_i), nb_sliding_windows),
                                                                dimnames = append(x = dimnames(range_contraction_events_count_sliding_window_i), values = list(paste0("Stratum_", 1:nb_sliding_windows))))
    # Counts of events per maps
    all_events_count_per_maps_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$all_events_count_per_maps
    all_events_count_per_maps_all_sliding_windows <- array(data = NA,
                                                           dim = c(dim(all_events_count_per_maps_sliding_window_i), nb_sliding_windows),
                                                           dimnames = append(x = dimnames(all_events_count_per_maps_sliding_window_i), values = list(paste0("Stratum_", 1:nb_sliding_windows))))
    # Summary tables of all event types
    all_events_summary_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$all_events_summary_df
    all_events_summary_all_sliding_windows <- array(data = NA,
                                                    dim = c(dim(all_events_summary_sliding_window_i), nb_sliding_windows),
                                                    dimnames = append(x = dimnames(all_events_summary_sliding_window_i), values = list(paste0("Stratum_", 1:nb_sliding_windows))))
  }
  
  ## Save outputs in objects used to aggregate information across all sliding windows
  
  # Unique subset events counts
  unique_sub_counts_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$unique_sub_counts
  old_names <- names(unique_sub_counts_all_sliding_windows)
  unique_sub_counts_all_sliding_windows <- append(x = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sub_counts_all_sliding_windows, values = list(unique_sub_counts_sliding_window_i))
  names(unique_sub_counts_all_sliding_windows) <- c(old_names, paste0("Stratum_", i))
  DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sub_counts_all_sliding_windows <- unique_sub_counts_all_sliding_windows
  
  # Unique vicariance events counts
  unique_vic_counts_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$unique_vic_counts
  old_names <- names(DEC_J_BSM_counts_arrays_all_sliding_windows$unique_vic_counts_all_sliding_windows)
  DEC_J_BSM_counts_arrays_all_sliding_windows$unique_vic_counts_all_sliding_windows <- append(x = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_vic_counts_all_sliding_windows, values = list(unique_vic_counts_sliding_window_i))
  names(DEC_J_BSM_counts_arrays_all_sliding_windows$unique_vic_counts_all_sliding_windows) <- c(old_names, paste0("Stratum_", i))
  
  # Unique range inheritance (sympatry) events counts
  unique_sym_counts_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$unique_sym_counts
  old_names <- names(DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sym_counts_all_sliding_windows)
  DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sym_counts_all_sliding_windows <- append(x = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sym_counts_all_sliding_windows, values = list(unique_sym_counts_sliding_window_i))
  names(DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sym_counts_all_sliding_windows) <- c(old_names, paste0("Stratum_", i))
  
  # Source/dest dispersal events count matrices
  dispersal_events_count_matrices_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$dispersal_events_count_matrices
  dispersal_events_count_matrices_all_sliding_windows[ , , , ,i] <- dispersal_events_count_matrices_sliding_window_i
  DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows <- dispersal_events_count_matrices_all_sliding_windows
  
  # Source/dest dispersal events summary matrices
  dispersal_events_summary_matrices_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$dispersal_events_summary_matrices
  dispersal_events_summary_matrices_all_sliding_windows[ , , , ,i] <- dispersal_events_summary_matrices_sliding_window_i
  DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_summary_matrices_all_sliding_windows <- dispersal_events_summary_matrices_all_sliding_windows
  
  # Range extinction counts x states
  range_contraction_events_count_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$range_contraction_events_count
  range_contraction_events_count_all_sliding_windows[ , ,i] <- as.matrix(range_contraction_events_count_sliding_window_i)
  DEC_J_BSM_counts_arrays_all_sliding_windows$range_contraction_events_count_all_sliding_windows <- range_contraction_events_count_all_sliding_windows
  
  # Counts of events per maps
  all_events_count_per_maps_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$all_events_count_per_maps
  all_events_count_per_maps_all_sliding_windows[ , ,i] <- as.matrix(all_events_count_per_maps_sliding_window_i)
  DEC_J_BSM_counts_arrays_all_sliding_windows$all_events_count_per_maps_all_sliding_windows <- all_events_count_per_maps_all_sliding_windows
  
  # Summary tables of all event types
  all_events_summary_sliding_window_i <- DEC_J_BSM_counts_arrays_per_sliding_window_i$all_events_summary_df
  all_events_summary_all_sliding_windows[ , ,i] <- as.matrix(all_events_summary_sliding_window_i)
  DEC_J_BSM_counts_arrays_all_sliding_windows$all_events_summary_all_sliding_windows <- all_events_summary_all_sliding_windows
}

str(DEC_J_BSM_counts_arrays_all_sliding_windows, max.level = 2)

# Save summary arrays of event counts across all sliding windows
# saveRDS(object = DEC_J_BSM_counts_arrays_all_sliding_windows, file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
saveRDS(object = DEC_J_BSM_counts_arrays_all_sliding_windows, file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# saveRDS(object = DEC_J_BSM_counts_arrays_all_sliding_windows, file = paste0("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# saveRDS(object = DEC_J_BSM_counts_arrays_all_sliding_windows, file = paste0("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))


##### 4/ Plot overall event type counts (all bioregions aggregated) #####

# Load summary arrays of event counts across all sliding windows
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))

### 4.1/ Compute overall event counts and rates ####

## List event types for plots

event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")
# event_types_list_with_typo <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dipsersal (j)")

# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# dimnames(DEC_J_BSM_counts_arrays_all_sliding_windows$all_events_count_per_maps_all_sliding_windows)

## Extract counts only for targeted type of events
all_events_count_per_maps_all_sliding_windows <- DEC_J_BSM_counts_arrays_all_sliding_windows$all_events_count_per_maps_all_sliding_windows[event_types_list, ,]
# all_events_count_per_maps_all_sliding_windows <- DEC_J_BSM_counts_arrays_all_sliding_windows$all_events_count_per_maps_all_sliding_windows[event_types_list_with_typo, ,]

## Format data for ggplot
all_events_count_per_maps_all_sliding_windows_ggplot <- all_events_count_per_maps_all_sliding_windows %>%
  reshape2::melt()
names(all_events_count_per_maps_all_sliding_windows_ggplot) <- c("event_type", "map", "window", "counts")

# Match time with time-window
window_time_df <- data.frame(window = levels(all_events_count_per_maps_all_sliding_windows_ggplot$window),
                             time = mean_time_list)

all_events_count_per_maps_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_sliding_windows_ggplot, y = window_time_df)

# # Correct typo
# all_events_count_per_maps_all_sliding_windows_ggplot$event_type <- as.character(all_events_count_per_maps_all_sliding_windows_ggplot$event_type)
# all_events_count_per_maps_all_sliding_windows_ggplot$event_type[all_events_count_per_maps_all_sliding_windows_ggplot$event_type == "jump-dipsersal (j)"] <- "jump-dispersal (j)"
# all_events_count_per_maps_all_sliding_windows_ggplot$event_type <- as.factor(all_events_count_per_maps_all_sliding_windows_ggplot$event_type)

## Compute rates based on LTT
# Load LTT data
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Youngest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Oldest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

total_richness <- DEC_J_LTT_all_areas_mean_ggplot %>% 
  filter(areas == "total") %>% 
  select(time, mean_counts) %>% 
  rename(total_richness = mean_counts)
# Merge with LTT to get richness
all_events_count_per_maps_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_sliding_windows_ggplot, y = total_richness)
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_sliding_windows_ggplot <- all_events_count_per_maps_all_sliding_windows_ggplot %>% 
  mutate(LTT_rates = counts / total_richness / window_width)

## Compute rates based on residence times
# Load residence times data
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))

# Aggregated all bioregions
residence_times_per_sliding_windows_melted_df <- residence_times_per_sliding_windows_in_areas_melted_df %>% 
  group_by(smoothed_time, map) %>%
  summarize(total_time = sum(residence_time))
# Merge with residence times
all_events_count_per_maps_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_sliding_windows_ggplot,
                                                                  y = residence_times_per_sliding_windows_melted_df,
                                                                  join_by(map == map, time == smoothed_time))
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_sliding_windows_ggplot <- all_events_count_per_maps_all_sliding_windows_ggplot %>% 
  mutate(rates = counts / total_time)

## Aggregate across maps
all_events_mean_counts_overall_all_sliding_windows_ggplot <- all_events_count_per_maps_all_sliding_windows_ggplot %>%
  arrange(desc(time)) %>% 
  group_by(event_type, map) %>%
  mutate(cum_counts = cumsum(counts)) %>%
  group_by(event_type, time) %>% 
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_rates = mean(rates),
            mean_LTT_rates = mean(LTT_rates)) %>% 
  group_by(time) %>% 
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Color scheme for events
event_types_col <- setNames(object = c("limegreen", "darkgoldenrod1", "dodgerblue", "darkorchid"), nm = event_types_list)
event_types_labels <- str_to_sentence(names(event_types_col))

# Adjust order of event types
all_events_mean_counts_overall_all_sliding_windows_ggplot$event_type <- factor(all_events_mean_counts_overall_all_sliding_windows_ggplot$event_type, levels = event_types_list, labels = event_types_labels)
all_events_count_per_maps_all_sliding_windows_ggplot$event_type <- factor(all_events_count_per_maps_all_sliding_windows_ggplot$event_type, levels = event_types_list, labels = event_types_labels)

# Save overall (mean) event counts and rates

# saveRDS(all_events_count_per_maps_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_events_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(all_events_count_per_maps_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/all_events_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(all_events_count_per_maps_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/all_events_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(all_events_count_per_maps_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/all_events_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))

# saveRDS(all_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Event_types/all_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(all_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/all_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(all_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/all_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(all_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/all_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))


### 4.2/ Plot stacked bars of cumulative mean raw number of event types along sliding windows ####

# Plot the stacked bars of cumulative mean raw number of event types along sliding windows (raw counts)
  # Y = Cumulative counts
  # Fill = type of events
  # X = Time

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_cum_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_cum_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/Events_count_cum_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/Events_count_cum_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)

barplot_cum_count_per_events_all_sliding_windows_ggplot <- ggplot(data = all_events_mean_counts_overall_all_sliding_windows_ggplot,
                                                                   mapping = aes(y = mean_cum_counts, x = time, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           # col = NA
           ) +
  
  # Set plot title +
  ggtitle(label = paste0("Cumulative counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Cumulative counts of events") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_labels, values = unname(event_types_col)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_cum_count_per_events_all_sliding_windows_ggplot <- barplot_cum_count_per_events_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_cum_count_per_events_all_sliding_windows_ggplot)

dev.off()


### 4.3/ Plot stacked bars of mean raw number of event types along sliding windows ####

# Plot the stacked bars of mean raw number of event types along sliding windows (raw counts)
  # Y = Mean counts
  # Fill = type of events
  # X = Time

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_mean_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_mean_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/Events_count_mean_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/Events_count_mean_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)

barplot_raw_count_per_events_all_sliding_windows_ggplot <- ggplot(data = all_events_mean_counts_overall_all_sliding_windows_ggplot,
                                                                  mapping = aes(y = mean_counts, x = time, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           # col = NA
  ) +
  
  # Set plot title +
  ggtitle(label = paste0("Counts of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of events") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_labels, values = unname(event_types_col)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_raw_count_per_events_all_sliding_windows_ggplot <- barplot_raw_count_per_events_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_raw_count_per_events_all_sliding_windows_ggplot)

dev.off()


### 4.4/ Plot stacked bars of mean percentages of event types along sliding windows #### 

# Plot the stacked bars of mean raw number of event types along sliding windows (raw counts)
  # Y = Mean percentages
  # Fill = type of events
  # X = Time

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)

barplot_perc_count_per_events_all_sliding_windows_ggplot <- ggplot(data = all_events_mean_counts_overall_all_sliding_windows_ggplot,
                                                                  mapping = aes(y = mean_percs, x = time, fill = event_type)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           # col = NA
  ) +
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of events") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Event types", labels = event_types_labels, values = unname(event_types_col)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_perc_count_per_events_all_sliding_windows_ggplot <- barplot_perc_count_per_events_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_perc_count_per_events_all_sliding_windows_ggplot)

dev.off()


### 4.5/ Plot rates of event types along sliding windows ####

# Counts divided per nb of lineages and width of sliding window (time)
# Plot lines of rates of event types along sliding windows (raw counts)
  # Y = Mean rates
  # Group/Col = type of events
  # X = Time

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_rates_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/Events_count_perc_barplots_overall_all_sliding_windows_",window_width,"My.pdf"), width = 10, height = 6)

lines_rates_per_events_all_sliding_windows_ggplot <- ggplot(data = all_events_mean_counts_overall_all_sliding_windows_ggplot,
                                                                   mapping = aes(y = mean_rates, x = time, col = event_type)) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = all_events_count_per_maps_all_sliding_windows_ggplot,
  #             mapping = aes(y = rates, x = time, group = event_type, col = event_type, fill = event_type),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # Plot 1000 replicates
  geom_line(data = all_events_count_per_maps_all_sliding_windows_ggplot,
            mapping = aes(y = rates, x = time, group = interaction(event_type, map), col = event_type),
            alpha = 0.01,
            linewidth = 2.0
  ) +

  # Plot mean lines only
  geom_line(data = all_events_mean_counts_overall_all_sliding_windows_ggplot,
            mapping = aes(y = mean_rates, x = time, group = event_type, col = event_type),
            alpha = 1.0,
            linewidth = 2.0
            ) +

  # Set plot title +
  ggtitle(label = paste0("Rates of events\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Rates of events  [Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0,0.03)) +
  ylim(c(0,0.025)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     limits = c(100, 0) # Set limits
                     ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Event types", labels = event_types_labels, values = unname(event_types_col)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_events_all_sliding_windows_ggplot <- lines_rates_per_events_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(lines_rates_per_events_all_sliding_windows_ggplot)

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
  # Dim 5 = Stratum/Time-window
  
  return(all_unique_events_source_dest_array)
}



### 5.2/ Build array of events per source/dest areas ####

window_width <- 1  # For cumulative maps
window_width <- 5  # For network maps
window_width <- 10 # For barplots

# Load summary arrays of event counts across all sliding windows
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))

## 5.2.1/ Initiate final array

# Extract dimensions
areas_list <- dimnames(DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows)[[1]]
event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "range inheritance (y)", "jump-dispersal (j)")
maps_list <- dimnames(DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows)[[4]]
windows_list <- dimnames(DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows)[[5]]

areas_list
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- array(data = 0,
                                                                           dim = c(length(areas_list), length(areas_list), length(event_types_list), length(maps_list), length(windows_list)),
                                                                           dimnames = list(areas_list, areas_list, event_types_list, maps_list, windows_list))

## 5.2.2/ Get s, v, y from unique_events_df

## Get counts for subset speciation events
all_unique_sub_events_all_sliding_windows_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sub_counts_all_sliding_windows, areas_list = areas_list, event_type = "subset speciation (s)")

# all_unique_sub_events_all_sliding_windows_source_dest_array[ , , , "Map_1" , "Stratum_1"]
# sum(all_unique_sub_events_all_sliding_windows_source_dest_array[ , , , "Map_1" , "Stratum_1"])
# 
# DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sub_counts_all_sliding_windows[[1]][1, ]
# sum(DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sub_counts_all_sliding_windows[[1]][1, ])

## Counts need to be equal !

## Get counts for vicariance events
all_unique_vic_events_all_sliding_windows_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_vic_counts_all_sliding_windows, areas_list = areas_list, event_type = "vicariance (v)")

## Should be symmetrical
# all_unique_vic_events_all_sliding_windows_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Get counts for range inheritance events
all_unique_sym_events_all_sliding_windows_source_dest_array <- convert_unique_events_to_source_dest_matrices(unique_events_list = DEC_J_BSM_counts_arrays_all_sliding_windows$unique_sym_counts_all_sliding_windows, areas_list = areas_list, event_type = "range inheritance (y)")

## Should be only on the diagonal if no widespread sympatric speciation is allowed (in BAYAREALIKE models)
# all_unique_sym_events_all_sliding_windows_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Fill the array with all event types

# Dim 1 = Source
# Dim 2 = Dest
# Dim 3 = Event_type
# Dim 4 = Map
# Dim 5 = Stratum

DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"subset speciation (s)", , ] <- all_unique_sub_events_all_sliding_windows_source_dest_array
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"vicariance (v)", , ] <- all_unique_vic_events_all_sliding_windows_source_dest_array
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"range inheritance (y)", , ] <- all_unique_sym_events_all_sliding_windows_source_dest_array

## 5.2.3/ Get d and j from array

# Extract counts only for targeted type of events

# dispersal_event_types_list <- c("range extension (d)", "jump-dipsersal (j)")
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")

# Extract counts only for targeted type of events
all_dispersal_events_count_matrices_per_maps_all_sliding_windows <- DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows[ , ,dispersal_event_types_list, , ]
dim(all_dispersal_events_count_matrices_per_maps_all_sliding_windows)
dimnames(all_dispersal_events_count_matrices_per_maps_all_sliding_windows)

## Fill the array with all event types
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"range extension (d)", , ] <- all_dispersal_events_count_matrices_per_maps_all_sliding_windows[ , ,"range extension (d)", , ]
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"jump-dispersal (j)", , ] <- all_dispersal_events_count_matrices_per_maps_all_sliding_windows[ , ,"jump-dipsersal (j)", , ]
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,"jump-dispersal (j)", , ] <- all_dispersal_events_count_matrices_per_maps_all_sliding_windows[ , ,"jump-dispersal (j)", , ]

# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , , , "Map_1" , "Stratum_1"]

## Save the final array of counts of events per source/dest areas
# saveRDS(object = DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
saveRDS(object = DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# saveRDS(object = DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# saveRDS(object = DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))


### 5.3/ Build ggplot dataframe of events per source/dest areas ####

## Load the array of counts of events per source/dest areas
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))

## List event types for plots

plot_event_types_list <- c("range extension (d)", "subset speciation (s)", "vicariance (v)", "jump-dispersal (j)")

# Do not plot range inheritance (y) as this is the default event with no transition during speciation
# Anagenetic range-switching (a) are not allowed
# Range contraction (e) are not observed in this model (e = 0)

# Extract counts only for targeted type of events
all_events_count_matrices_per_maps_all_sliding_windows <- DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,plot_event_types_list, , ]

dim(all_events_count_matrices_per_maps_all_sliding_windows)
# Dim 1 = Rows = Sources
# Dim 2 = Cols = Dest
# Dim 3 = Event types
# Dim 4 = Maps
# Dim 5 = Sliding windows

## Format data for ggplot
all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- all_events_count_matrices_per_maps_all_sliding_windows %>%
  reshape2::melt()
names(all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot) <- c("source", "dest", "event_type", "map", "window", "counts")

# Match time with time-window
window_time_df <- data.frame(window = levels(all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot$window),
                             time = mean_time_list)

all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, y = window_time_df)

# Save ggplot df dataframe of events per areas transitions x time-windows
# saveRDS(object = all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))


### 5.4/ Aggregate per ggplot dataframe of events per source areas x time-windows ####

# all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot %>% 
  group_by(map, event_type, source, time) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()

## Compute rates based on LTT
# Load LTT data
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Youngest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Oldest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

total_richness_per_bioregions <- DEC_J_LTT_all_areas_mean_ggplot %>% 
  filter(areas != "total") %>% 
  select(areas, time, mean_counts) %>% 
  rename(area_richness = mean_counts)
# Merge with LTT to get richness
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, y = total_richness_per_bioregions, join_by(time == time, source == areas))
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  mutate(LTT_rates = counts / area_richness / window_width)

## Compute rates based on residence times
# Load residence times data
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))

# Merge with residence times
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                        y = residence_times_per_sliding_windows_in_areas_melted_df[, c("area", "smoothed_time", "map", "residence_time")],
                                                                                        join_by(source == area, time == smoothed_time, map == map))
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  mutate(rates = counts / residence_time)

# Replace NA (due to dividing by zero) with 0.000
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.na(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.nan(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Adjust order of Bioregions
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))


### 5.5/ Aggregate per ggplot dataframe of events per destination areas x time-windows ####

all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot %>% 
  group_by(map, event_type, dest, time) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup()

## Compute rates based on LTT
# Load LTT data
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Youngest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Oldest_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

total_richness_per_bioregions <- DEC_J_LTT_all_areas_mean_ggplot %>% 
  filter(areas != "total") %>% 
  select(areas, time, mean_counts) %>% 
  rename(area_richness = mean_counts)
# Merge with LTT to get richness
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot, y = total_richness_per_bioregions, join_by(time == time, dest == areas))
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  mutate(LTT_rates = counts / area_richness / window_width)

## Compute rates based on residence times
# Load residence times data
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))

# Merge with residence times
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- left_join(x = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                        y = residence_times_per_sliding_windows_in_areas_melted_df[ ,c("area", "smoothed_time", "map", "residence_time")],
                                                                                        join_by(dest == area, time == smoothed_time, map == map))
# Divide counts by richness and time span to get rates in events per lineage per My
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  mutate(rates = counts / residence_time)

# Replace NA (due to dividing by zero) with 0.000
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates[is.na(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates[is.nan(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates[(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Adjust order of Bioregions
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$dest <- factor(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot$dest, levels = areas_list, labels = bioregion_names)

# Save ggplot df
# saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))


##### 6/ Plot event type counts per source bioregions #####

### 6.1/ Aggregate counts and rates across maps ####

# Load ggplot df of event types counts per source bioregions
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

## Aggregate across maps ####

all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>%
  arrange(desc(time)) %>% 
  group_by(event_type, source, map) %>%
  mutate(cum_counts = cumsum(counts)) %>%
  group_by(event_type, source, time) %>% 
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_LTT_rates = mean(LTT_rates),
            mean_rates = mean(rates)) %>% 
  group_by(event_type, time) %>% 
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Adjust order of event types
all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$event_type <- factor(all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$event_type, levels = plot_event_types_list, labels = plot_event_types_list)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Adjust order of Bioregions
all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)

# Save ggplot df of of event types counts aggregated per source bioregions
# saveRDS(object = all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_",window_width,"My_ggplot.rds"))
saveRDS(object = all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_",window_width,"My_ggplot.rds"))

### 6.2/ Plot stacked bars of cumulative mean raw number of event types per source bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                  mapping = aes(y = mean_cum_counts, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Cumulative counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Cumulative counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Source bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}


## Multifaceted plot
# Rows = type of events
# Columns = Time-strata

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_cum_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_cum_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 6.3/ Plot stacked bars of mean raw number of event types per source bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                          mapping = aes(y = mean_counts, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Source bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_mean_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_mean_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 6.4/ Plot stacked bars of mean percentages of event types per source bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                           mapping = aes(y = mean_percs, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Source bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_perc_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_perc_per_types_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 6.5/ Plot rates of event types per source bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data across maps
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Extract mean event count data
  lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Set maximum of Y-axis
  # y_max <- max(lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i$mean_rates, na.rm = T) * 0.1
  
  # Generate plot
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = lines_rates_per_events_all_source_bioregions_all_sliding_windows_i) +
    
    # # Plot mean lines + 1000 replicates
    # geom_smooth(data = lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i ,
    #             mapping = aes(y = rates, x = time, group = source, col = source, fill = source),
    #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
    
    # # Plot 1000 replicates
    # geom_line(data = lines_rates_per_events_all_source_bioregions_all_sliding_windows_i,
    #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
    #           alpha = 0.01,
    #           linewidth = 2.0
    # ) +
    
    # Plot mean lines only
    geom_line(data = lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i ,
              mapping = aes(y = mean_rates, x = time, group = source, col = source),
              alpha = 1.0,
              linewidth = 2.0
    ) +
    
    # Set plot title +
    ggtitle(label = paste0("Rates of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Rates of events  [Events / lineage / My]") +
    
    # Set y-axis limits
    # ylim(c(0, 0.03)) +
    # ylim(c(0, y_max)) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse",
                       limits = c(100, 0) # Set limits
    ) + 
    
    # Adjust color scheme and legend
    scale_color_manual("Source bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Remove fill legend
    guides(fill = "none") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_rates_per_types_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_rates_per_types_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


##### 7/ Plot event type counts per destination bioregions #####

### 7.1/ Aggregate counts and rates across maps ####

# Load ggplot df of event types counts per destination bioregions
# all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

## Aggregate across maps ####
all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>%
  arrange(desc(time)) %>% 
  group_by(event_type, dest, map) %>%
  mutate(cum_counts = cumsum(counts)) %>%
  group_by(event_type, dest, time) %>% 
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_LTT_rates = mean(LTT_rates),
            mean_rates = mean(rates)) %>% 
  group_by(event_type, time) %>% 
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Adjust order of event types
all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$event_type <- factor(all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$event_type, levels = plot_event_types_list, labels = plot_event_types_list)

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Adjust order of Bioregions
all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest <- factor(all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest, levels = bioregion_names, labels = bioregion_names)

# Save ggplot df of of event types counts aggregated per destination bioregions
# saveRDS(object = all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot_",window_width,"My_ggplot.rds"))
saveRDS(object = all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot_",window_width,"My_ggplot.rds"))

### 7.2/ Plot stacked bars of cumulative mean raw number of event types per destination bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                          mapping = aes(y = mean_cum_counts, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Cumulative counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Cumulative counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Destination bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}


## Multifaceted plot
  # Rows = type of events
  # Columns = Time-strata

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_cum_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_cum_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.3/ Plot stacked bars of mean raw number of event types per destination bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                           mapping = aes(y = mean_counts, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Destination bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_mean_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_mean_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.4/ Plot stacked bars of mean percentages of event types per destination bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                           mapping = aes(y = mean_percs, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Destination bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_perc_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_perc_per_types_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


### 7.5/ Plot rates of event types per destination bioregions along sliding windows ####

# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

## Generate ggplots per type of events
lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data across maps
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Extract mean event count data
  lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Set maximum of Y-axis
  # y_max <- max(lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i$mean_rates, na.rm = T) * 0.1
  
  # Generate plot
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i) +
    
    # # Plot mean lines + 1000 replicates
    # geom_smooth(data = lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i ,
    #             mapping = aes(y = rates, x = time, group = dest, col = event_type, fill = dest),
    #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
    
    # # Plot 1000 replicates
    # geom_line(data = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i,
    #           mapping = aes(y = rates, x = time, group = interaction(dest, map), col = dest),
    #           alpha = 0.01,
    #           linewidth = 2.0
  # ) +
  
  # Plot mean lines only
  geom_line(data = lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i ,
            mapping = aes(y = mean_rates, x = time, group = dest, col = dest),
            alpha = 1.0,
            linewidth = 2.0
  ) +
    
    # Set plot title +
    ggtitle(label = paste0("Rates of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Rates of events  [Events / lineage / My]") +
    
    # Set y-axis limits
    # ylim(c(0, 0.03)) +
    # ylim(c(0, y_max)) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse",
                       limits = c(100, 0) # Set limits
    ) + 
    
    # Adjust color scheme and legend
    scale_color_manual("Destination bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Remove fill legend
    guides(fill = "none") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}


# Multifaceted: 
  # Rows = type of events
  # X = Time-windows
  # Fill/Groups = bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_rates_per_types_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_rates_per_types_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = 1,  # Width of columns
  heights = rep(1, length(plot_event_types_list)),
  nrow = length(plot_event_types_list),
  ncol = 1,
)

dev.off()


##### 8/ Plot event type counts per Source & Destination bioregions #####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Multifaceted: 
  # Rows = type of events
  # Columns = Source and Destination bioregions
  # Fill/Groups = bioregions

### 8.1/ Plot stacked bars of cumulative mean raw number of event types per source/destination bioregions along sliding windows ####

## Extract maximum cumulative count number per time-windows across both Source & Dest bioregions

all_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, event_type) %>%
  summarize(sum_cum_counts = sum(mean_cum_counts)) %>%
  group_by(event_type) %>%
  summarize(max_cum_counts_source = max(sum_cum_counts))
all_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, event_type) %>%
  summarize(sum_cum_counts = sum(mean_cum_counts)) %>%
  group_by(event_type) %>%
  summarize(max_cum_counts_dest = max(sum_cum_counts))

max_cum_source_dest_counts <- left_join(x = all_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot, all_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot)
max_cum_source_dest_counts
# Counts of events should be similar between source and destination
max_cum_source_dest_counts <- max_cum_source_dest_counts %>% 
  group_by(event_type) %>% 
  mutate(max_cum_counts = max(max_cum_counts_source, max_cum_counts_dest)) %>% 
  ungroup()

## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                          mapping = aes(y = mean_cum_counts, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3, 
             show.legend = F) + # Remove legend
    
    # Set Y-axis to all be the same across bioregion types
    ylim(0, max_cum_source_dest_counts$max_cum_counts[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +

    # Set plot title +
    ggtitle(label = paste0("Cumulative counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Cumulative counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                        mapping = aes(y = mean_cum_counts, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Set Y-axis to all be the same across bioregion types
    ylim(0, max_cum_source_dest_counts$max_cum_counts[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Cumulative counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Cumulative counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_cum_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_cum_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = barplot_count_cum_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list)

## Export multifaceted plot
  # Rows = type of events
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_cum_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_cum_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_cum_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
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


### 8.2/ Plot stacked bars of mean raw number of event types per source/destination bioregions along sliding windows ####

## Extract maximum mean count number per time-windows across both Source & Dest bioregions

all_events_max_mean_counts_per_source_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, event_type) %>%
  summarize(sum_mean_counts = sum(mean_counts)) %>%
  group_by(event_type) %>%
  summarize(max_mean_counts_source = max(sum_mean_counts))
all_events_max_mean_counts_per_dest_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, event_type) %>%
  summarize(sum_mean_counts = sum(mean_counts)) %>%
  group_by(event_type) %>%
  summarize(max_mean_counts_dest = max(sum_mean_counts))

max_mean_source_dest_counts <- left_join(x = all_events_max_mean_counts_per_source_bioregions_all_sliding_windows_ggplot, all_events_max_mean_counts_per_dest_bioregions_all_sliding_windows_ggplot)
max_mean_source_dest_counts
# Counts of events should be similar between source and destination
max_mean_source_dest_counts <- max_mean_source_dest_counts %>% 
  group_by(event_type) %>% 
  mutate(max_mean_counts = max(max_mean_counts_source, max_mean_counts_dest)) %>% 
  ungroup()

## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                           mapping = aes(y = mean_counts, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3,
             show.legend = F) + # Remove legend
    
    # Set Y-axis to all be the same across bioregion types
    ylim(0, max_mean_source_dest_counts$max_mean_counts[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                         mapping = aes(y = mean_counts, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Set Y-axis to all be the same across bioregion types
    ylim(0, max_mean_source_dest_counts$max_mean_counts[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Counts of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Counts of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_mean_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_mean_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = barplot_count_mean_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list)

## Export multifaceted plot
# Rows = type of events
# Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_mean_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_mean_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_mean_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
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


### 8.3/ Plot stacked bars of percentages of event types per source/destination bioregions along sliding windows ####

## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_i,
                                                                                           mapping = aes(y = mean_percs, x = time, fill = source)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3,
             show.legend = F) + # Remove legend
    
    # Set Y-axis to all be the same across bioregion types
    # ylim(0, 100) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 2
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Generate plot
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_i,
                                                                                         mapping = aes(y = mean_percs, x = time, fill = dest)) +
    # Plot barplot
    geom_col(alpha = 1.0, position = "stack", width = 1.0,
             col = "black", linewidth = 0.3) +
    
    # Set Y-axis to all be the same across bioregion types
    # ylim(0, 100) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse") +
    #                  limits = c(100, 0) # Set limits
    
    # Set plot title +
    ggtitle(label = paste0("Percentages of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("% of events") +
    
    # Adjust color scheme and legend
    scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_perc_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = barplot_count_perc_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = barplot_count_perc_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list)

## Export multifaceted plot
  # Rows = type of events
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_count_perc_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_count_perc_per_types_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = barplot_count_perc_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
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


### 8.4/ Plot rates of event types per source/destination bioregions along sliding windows ####

## Extract maximum rates per time-windows across both Source & Dest bioregions

all_events_max_rates_per_source_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(event_type) %>%
  summarize(max_rates_source = max(mean_rates, na.rm = T))
all_events_max_rates_per_dest_bioregions_all_sliding_windows_ggplot <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(event_type) %>%
  summarize(max_rates_dest = max(mean_rates, na.rm = T))


max_mean_source_dest_rates <- left_join(x = all_events_max_rates_per_source_bioregions_all_sliding_windows_ggplot, all_events_max_rates_per_dest_bioregions_all_sliding_windows_ggplot)
max_mean_source_dest_rates
# Rates of events does not need to be similar between source and destination
max_mean_source_dest_rates <- max_mean_source_dest_rates %>% 
  group_by(event_type) %>% 
  mutate(max_mean_rates = max(max_rates_source, max_rates_dest)) %>% 
  ungroup()

## Generate ggplots per type of events for Source Bioregions (Remove legend)
lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data across maps
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Extract mean event count data
  lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Set maximum of Y-axis
  # y_max <- max(lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i$mean_rates, na.rm = T) * 0.1
  
  # Generate plot
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = lines_rates_per_events_all_source_bioregions_all_sliding_windows_i) +
    
    # # Plot mean lines + 1000 replicates
    # geom_smooth(data = lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i ,
    #             mapping = aes(y = rates, x = time, group = source, col = event_type, fill = source),
    #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
    
    # # Plot 1000 replicates
    # geom_line(data = lines_rates_per_events_all_source_bioregions_all_sliding_windows_i,
    #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
    #           alpha = 0.01,
    #           linewidth = 2.0
  # ) +
  
  # Plot mean lines only
  geom_line(data = lines_mean_rates_per_events_all_source_bioregions_all_sliding_windows_i ,
            mapping = aes(y = mean_rates, x = time, group = source, col = source),
            alpha = 1.0,
            linewidth = 2.0,
            show.legend = F  # Remove legend
  ) +
    
    # Set plot title +
    ggtitle(label = paste0("Rates of ", event_type_i, " events\nper Source bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Rates of events  [Events / lineage / My]") +
    
    # Set Y-axis to all be the same across bioregion types
    # ylim(0, max_mean_source_dest_rates$max_mean_rates[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse",
                       limits = c(100, 0) # Set limits
    ) + 
    
    # Adjust color scheme and legend
    scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Remove fill legend
    guides(fill = "none") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list <- append(x = lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = list(lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot))
  
}

## Generate ggplots per type of events for Destination Bioregions
lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- list()
for ( i in seq_along(plot_event_types_list))
{
  # i <- 1
  
  # Extract event type
  event_type_i <- plot_event_types_list[i]
  
  # Extract event count data across maps
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Extract mean event count data
  lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i <- all_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    filter(event_type == event_type_i)
  
  # Set maximum of Y-axis
  # y_max <- max(lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i$mean_rates, na.rm = T) * 0.1
  
  # Generate plot
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i) +
    
    # # Plot mean lines + 1000 replicates
    # geom_smooth(data = lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i ,
    #             mapping = aes(y = rates, x = time, group = dest, col = event_type, fill = dest),
    #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
    
    # # Plot 1000 replicates
    # geom_line(data = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_i,
    #           mapping = aes(y = rates, x = time, group = interaction(dest, map), col = dest),
    #           alpha = 0.01,
    #           linewidth = 2.0
  # ) +
  
  # Plot mean lines only
  geom_line(data = lines_mean_rates_per_events_all_dest_bioregions_all_sliding_windows_i ,
            mapping = aes(y = mean_rates, x = time, group = dest, col = dest),
            alpha = 1.0,
            linewidth = 2.0
  ) +
    
    # Set plot title +
    ggtitle(label = paste0("Rates of ", event_type_i, " events\nper Destination bioregions\nacross 1000 BS maps")) +
    
    # Set axes labels
    xlab("Time  [My]") +
    ylab("Rates of events  [Events / lineage / My]") +
    
    # Set Y-axis to all be the same across bioregion types
    # ylim(0, max_mean_source_dest_rates$max_mean_rates[i]) +
    
    # Reverse time scale
    scale_x_continuous(transform = "reverse",
                       limits = c(100, 0) # Set limits
    ) + 
    
    # Adjust color scheme and legend
    scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
    
    # Remove fill legend
    guides(fill = "none") +
    
    # Adjust margins
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 
  
  # Adjust aesthetics
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot <- lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
    add_aesthetics_barplot()
  
  # Plot
  # print(lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot)
  
  # Store ggplot in a list
  lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list, values = list(lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot))
  
}

## Append together Source and Dest multifaceted plots of mean counts per event types
lines_rates_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = lines_rates_per_events_all_source_bioregions_all_sliding_windows_ggplot_list, values = lines_rates_per_events_all_dest_bioregions_all_sliding_windows_ggplot_list)

## Export multifaceted plot
  # Rows = type of events
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Events_rates_per_types_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/Events_rates_per_types_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = length(plot_event_types_list) * 7)

gridExtra::grid.arrange(
  grobs = lines_rates_per_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
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


##### 9/ Plot all dispersal event counts per source bioregions #####

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# dispersal_event_types_list_with_typo <- c("range extension (d)", "jump-dipsersal (j)")

### 9.1/ Aggregate dispersal event counts and rates per source bioregions ####

# Load ggplot df of of event types counts aggregated per source bioregions
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

names(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)

# Extract and sum counts of dispersal events
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(time, source, map) %>%
  summarize(counts = sum(counts), # Sum mean counts across dispersal event types
            area_richness = mean(area_richness),
            residence_time = mean(residence_time)) %>% 
  mutate(LTT_rates = counts / area_richness / window_width) %>% # Compute updated rates for all dispersal event types based on LTT data
  mutate(rates = counts / residence_time) %>% # Compute updated rates for all dispersal event types based on residence times
  arrange(desc(time)) %>% 
  group_by(source, map) %>% 
  mutate(cum_counts = cumsum(counts)) # Compute cumulative counts along time

# Replace NA (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates[is.na(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates[is.nan(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates[(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Aggregate mean counts of dispersal events across BSM maps
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(source, time) %>%  # Aggregate across maps
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_LTT_rates = mean(LTT_rates),
            mean_rates = mean(rates)) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Adjust order of Bioregions
# all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)

# Save ggplot df of (mean) counts per source bioregions across dispersal events
# saveRDS(object = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))


### 9.2/ Plot stacked bars of cumulative mean raw number of dispersal events per source bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_cum_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                  mapping = aes(y = mean_cum_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Cumulative counts of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Cumulative counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Source\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


### 9.3/ Plot stacked bars of mean raw number of dispersal events per source bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                  mapping = aes(y = mean_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Source\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


### 9.4/ Plot stacked bars of percentages of dispersal events per source bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_percs, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Source\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


### 9.5/ Plot rates of dispersal event types per source bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 2) %>%# Minimum 2 lineages to compute meaningful rates
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot,
  #             mapping = aes(y = rates, x = time, group = source, col = source, fill = source),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot,
  #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # Plot mean lines only
  geom_line(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
            mapping = aes(y = mean_rates, x = time, group = source, col = source),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("Rates of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Rates of dispersal events\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits for Oldest
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Source\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 10/ Plot all dispersal event counts per destination bioregions #####

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# dispersal_event_types_list_with_typo <- c("range extension (d)", "jump-dipsersal (j)")

### 10.1/ Aggregate dispersal event counts and rates per destination bioregions ####

# Load ggplot df of of event types counts aggregated per destination bioregions
# all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

names(all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot)

# Extract and sum counts of dispersal events 
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(time, dest, map) %>%
  summarize(counts = sum(counts), # Sum mean counts across dispersal event types
            area_richness = mean(area_richness),
            residence_time = mean(residence_time)) %>% 
  mutate(LTT_rates = counts / area_richness / window_width) %>% # Compute updated rates for all dispersal event types based on LTT data
  mutate(rates = counts / residence_time) %>% # Compute updated rates for all dispersal event types based on residence times
  arrange(desc(time)) %>% 
  group_by(dest, map) %>% 
  mutate(cum_counts = cumsum(counts)) # Compute cumulative counts along time
  
# Replace NA (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates[is.na(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates[is.nan(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates[(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates[(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Aggregate mean counts across dispersal events
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(dest, time) %>%  # Aggregate across maps
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_rates = mean(rates)) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Adjust order of Bioregions
all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest <- factor(all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest, levels = bioregion_names, labels = bioregion_names)
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest <- factor(all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$dest, levels = bioregion_names, labels = bioregion_names)

# Save ggplot df of (mean) counts per destination bioregions across dispersal events
# saveRDS(object = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

### 10.2/ Plot stacked bars of cumulative mean raw number of dispersal events per destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_cum_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                  mapping = aes(y = mean_cum_counts, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Cumulative counts of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Cumulative counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Destination\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)

dev.off()


### 10.3/ Plot stacked bars of mean raw number of dispersal events per destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_mean_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_counts, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Destination\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)

dev.off()

### 10.4/ Plot stacked bars of percentages of dispersal events per destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_perc_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"), 
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"), 
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"), 
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"), 
    width = 10, height = 6)

# Generate plot
barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_percs, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Destination\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)

dev.off()


### 10.5/ Plot rates of dispersal event types per destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 6) %>% # Minimum 6 lineages to compute meaningful rates
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_rates_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot,
  #             mapping = aes(y = rates, x = time, group = dest, col = dest, fill = dest),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot,
  #           mapping = aes(y = rates, x = time, group = interaction(dest, map), col = dest),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

# Plot mean lines only
geom_line(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
          mapping = aes(y = mean_rates, x = time, group = dest, col = dest),
          alpha = 1.0,
          linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("Rates of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Rates of dispersal events\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits for Oldest
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Destination\nbioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 11/ Plot all dispersal event counts per Source & Destination bioregions #####

# Plot results in facing columns with similar Y-axes
# Remove legend from Source plots in the first columns

# Multifaceted: 
  # Columns = Source and Destination bioregions
  # Fill/Groups = bioregions

### 11.1/ Plot stacked bars of cumulative mean raw number of dispersal event types per source/destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

## Extract maximum cumulative count number per time-windows across both Source & Dest bioregions

all_dispersal_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%
  summarize(sum_cum_counts = sum(mean_cum_counts)) # Sum across source bioregions
max_cum_counts_source <- max(all_dispersal_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot$sum_cum_counts)

all_dispersal_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%
  summarize(sum_cum_counts = sum(mean_cum_counts)) # Sum across destination bioregions
max_cum_counts_dest <- max(all_dispersal_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot$sum_cum_counts)

# Counts of dispersal events should be similar between source and destination
max_cum_source_dest_counts <- max(max_cum_counts_source, max_cum_counts_dest)



## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                  mapping = aes(y = mean_cum_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           show.legend = F) +  # Remove legend
  
  # Set Y-axis to all be the same across bioregion types
  ylim(0, max_cum_source_dest_counts) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Cumulative counts of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Cumulative counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

## Generate ggplots per type of events for Destination Bioregions
barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                mapping = aes(y = mean_cum_counts, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Set Y-axis to all be the same across bioregion types
  ylim(0, max_cum_source_dest_counts) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Cumulative counts of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Cumulative counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_cum_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = list(barplot_count_cum_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot), values = list(barplot_count_cum_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot))

## Export multifaceted plot
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_cum_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_cum_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = 6)

gridExtra::grid.arrange(
  grobs = barplot_count_cum_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1, # Height of rows
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()


### 11.2/ Plot stacked bars of mean raw number of dispersal events per Source & Destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

## Extract maximum cumulative count number per time-windows across both Source & Dest bioregions

all_dispersal_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%
  summarize(sum_counts = sum(mean_counts)) # Sum across source bioregions
max_counts_source <- max(all_dispersal_events_max_cum_counts_per_source_bioregions_all_sliding_windows_ggplot$sum_counts)

all_dispersal_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%
  summarize(sum_counts = sum(mean_counts)) # Sum across destination bioregions
max_counts_dest <- max(all_dispersal_events_max_cum_counts_per_dest_bioregions_all_sliding_windows_ggplot$sum_counts)

# Counts of dispersal events should be similar between source and destination
max_source_dest_counts <- max(max_counts_source, max_counts_dest)


## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           show.legend = F) +
  
  # Set Y-axis to all be the same across bioregion types
  ylim(0, max_source_dest_counts) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

## Generate ggplots per type of events for Destination Bioregions
barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                 mapping = aes(y = mean_counts, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Set Y-axis to all be the same across bioregion types
  ylim(0, max_source_dest_counts) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_mean_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = list(barplot_count_mean_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot), values = list(barplot_count_mean_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot))

## Export multifaceted plot
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_mean_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_mean_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
   width = 17, height = 6)

gridExtra::grid.arrange(
  grobs = barplot_count_mean_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1, # Height of rows
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()


### 11.3/ Plot stacked bars of percentages of dispersal events per Source & Destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

## Generate ggplots per type of events for Source Bioregions (Remove legend)
barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_percs, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3,
           show.legend = F) + # Remove legend
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

## Generate ggplots per type of events for Destination Bioregions
barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
                                                                                                 mapping = aes(y = mean_percs, x = time, fill = dest)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) +# Set limits for Oldest
                     limits = c(100, 0)) +# Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of dispersal events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
barplot_count_perc_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = list(barplot_count_perc_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot), values = list(barplot_count_perc_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot))

## Export multifaceted plot
  # Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_count_perc_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_count_perc_barplots_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
            width = 17, height = 6)

gridExtra::grid.arrange(
  grobs = barplot_count_perc_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1, # Height of rows
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()



### 11.4/ Plot rates of dispersal event types per destination bioregions along sliding windows ####

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

# Load ggplot df of mean counts per source bioregions across dispersal events
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Filter out records when cumulative mean counts of events is less than one to avoid accounting for outliers before colonization is agreed by the majority of BSM simulations
all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(mean_cum_counts >= 1) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup()

## Extract maximum rates per sliding windows across both Source & Dest bioregions

max_rates_source <- max(all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$mean_rates, na.rm = T)
max_rates_dest <- max(all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot$mean_rates, na.rm = T)

# Rates do NOT have to be similar between source and destination
max_source_dest_rates <- max(max_rates_source, max_rates_dest)


## Generate ggplots per type of events for Source Bioregions (Remove legend)
lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot,
  #             mapping = aes(y = rates, x = time, group = source, col = source, fill = source),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = all_dispersal_events_counts_all_source_bioregions_all_sliding_windows_ggplot,
  #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # Plot mean lines only
  geom_line(data = all_dispersal_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
            mapping = aes(y = mean_rates, x = time, group = source, col = source),
            alpha = 1.0,
            linewidth = 2.0,
            show.legend = F) + # Remove legend
  
  # Set plot title +
  ggtitle(label = paste0("Rates of dispersal events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Rates of dispersal events\n[Events / lineage / My]") +
  
  # Set Y-axis to all be the same across bioregion types
  # ylim(c(0, max_source_dest_rates)) +

  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits for Oldest
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot)

## Generate ggplots per type of events for Destination Bioregions
lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot,
  #             mapping = aes(y = rates, x = time, group = dest, col = dest, fill = dest),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = all_dispersal_events_counts_all_dest_bioregions_all_sliding_windows_ggplot,
  #           mapping = aes(y = rates, x = time, group = interaction(dest, map), col = dest),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # Plot mean lines only
  geom_line(data = all_dispersal_events_mean_counts_all_dest_bioregions_all_sliding_windows_ggplot,
            mapping = aes(y = mean_rates, x = time, group = dest, col = dest),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("Rates of dispersal events\nper Destination bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Rates of dispersal events\n[Events / lineage / My]") +
  
  # Set Y-axis to all be the same across bioregion types
  # ylim(c(0, max_source_dest_rates)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits for Oldest
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot <- lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
# print(lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot)


## Append together Source and Dest multifaceted plots of mean counts per event types
lines_rates_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list <- append(x = list(lines_rates_per_dispersal_events_all_source_bioregions_all_sliding_windows_ggplot), values = list(lines_rates_per_dispersal_events_all_dest_bioregions_all_sliding_windows_ggplot))

## Export multifaceted plot
# Columns = Source and Destination bioregions

# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/Dispersal_events_rates_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Dispersal_events/Dispersal_events_rates_all_source_dest_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 17, height = 6)

gridExtra::grid.arrange(
  grobs = lines_rates_per_dispersal_events_all_source_dest_bioregions_all_sliding_windows_ggplot_list, # List of ggplots
  widths = c(6.5, 9),  # Width of columns
  heights = 1, # Height of rows
  nrow = 1,
  ncol = 2,
  layout_matrix = rbind(c(1, 2)) # Position of ggplots in the layout
)

dev.off()


##### 12/ Map networks of mean counts of dispersal events along sliding windows ####

# Sliding series => use the 5 My window to have enough events per window

window_width <- 5

### 12.1/ Aggregate counts for all dispersal events, per sliding windows ####

# Load ggplot df
all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
names(all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Sum across event_types, per sliding windows
# Aggregate mean counts across maps
all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(source, dest, map, time) %>% 
  summarize(counts = sum(counts)) %>% # Sum counts across event_types
  arrange(desc(time)) %>% 
  group_by(source, dest, map) %>%
  mutate(cum_counts = cumsum(counts)) %>%  # Compute cumsum along time
  group_by(source, dest, time) %>%
  summarize(mean_counts = mean(counts),                # Aggregate mean counts across maps
            mean_cum_counts = mean(cum_counts)) %>%    # Aggregate mean sum counts across maps
  group_by(time) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%  # Compute overall percentages for a given time
  ungroup()


# Save ggplot df of mean counts of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))

## Convert back to array

# Pivot mean counts and percentages
all_dispersal_events_mean_count_per_maps_all_sliding_windows_melted <- all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot %>% 
  pivot_longer(cols = c("mean_counts", "mean_perc", "mean_cum_counts"), names_to = "mean_stats", values_to = "values")

all_dispersal_events_mean_count_per_maps_all_sliding_windows_array <- reshape2::acast(data = all_dispersal_events_mean_count_per_maps_all_sliding_windows_melted,
                                                                                  formula = source ~ dest ~ time ~ mean_stats)

dimnames(all_dispersal_events_mean_count_per_maps_all_sliding_windows_array)
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_counts"]
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_perc"]
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_cum_counts"]

# Save array of mean counts/perc of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_array.rds"))
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_array.rds"))

### 12.2/ Compute node metadata per sliding windows ####

# Extract root age
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/DEC_J_BSM_output.rds")
# DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
root_age <- max(DEC_J_clado_events_tables[[1]]$time_bp)

# Generate start and end times of sliding windows
start_time_list <- seq(from = 0, to = root_age - window_width, by = window_steps)
end_time_list <- seq(from = window_width, to = root_age, by = window_steps)
mean_time_list <- (start_time_list + end_time_list) / 2

# Load LTT data for species richness in time
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

## Create table for bioregion coordinates
bioregion_names_alpha_order <- c("Afrotropics", "Eastern Palearctic", "Indomalaya", "Neotropics", "Nearctic", "Australasia", "Western Palearctic")
areas_list <- levels(as.factor(DEC_J_LTT_all_areas_mean_ggplot$areas))
areas_list <- areas_list[areas_list != "total"]
areas_alpha_order <- areas_list[order(areas_list)]
bioregion_names_alpha_order_df <- data.frame(areas = areas_alpha_order, bioregions = bioregion_names_alpha_order)

## Create table for bioregion coordinates (May be adjusted per time strata to account for continental drift)
bioregion_coordinates_df <- bioregion_names_alpha_order_df
bioregion_coordinates_df$latitude <- c(0, 50, 10, -10, 40, -25, 47)
bioregion_coordinates_df$longitude <- c(25, 100, 100, -60, -100, 140, 20)
bioregion_coordinates_df$adjusted_latitude <- c(10, 35, 10, -10, 35, -15, 40)
bioregion_coordinates_df$adjusted_longitude <- c(15, 120, 100, -60, -100, 150, 15)

## Create node metadata with species richness at mean time of sliding windows
nodes_metadata_per_sliding_windows <- list()
for (i in 1:length(start_time_list))
{
  # i <- 1
  
  # Extract mean time boundaries of sliding windows
  mean_time <- mean_time_list[i]
  
  # Aggregate species richness within each time stratum
  nodes_metadata_i <- DEC_J_LTT_all_areas_mean_ggplot %>% 
    filter(time == mean_time) %>% 
    filter(areas != "total") %>% 
    mutate(mean_percentages = 100 * mean_counts / sum(mean_counts)) %>%
    mutate(bioregions = areas) %>%
    rename(node_ID = areas) %>%
    mutate(node_size = log(mean_counts)) %>%
    mutate(window = paste0("Window_",i)) %>%
    mutate(time = mean_time) %>%
    select(window, time, node_ID, mean_counts, mean_percentages, node_size)
  
  # Add bioregion labels
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_names_alpha_order_df, by = join_by(node_ID == areas)) %>%
    select(window, time, node_ID, bioregions, mean_counts, mean_percentages, node_size)
  
  # Add label with richness
  nodes_metadata_i$node_labels <- paste0(nodes_metadata_i$node_ID, "\n", round(nodes_metadata_i$mean_counts, 0))
  
  # Reorder in alphabetic order of node ID
  nodes_metadata_i <- nodes_metadata_i %>% 
    arrange(node_ID)
  
  ## Add spatial coordinates for each bioregion (May adjust to show continental drift)
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_coordinates_df,
                                by = join_by(node_ID == areas, bioregions == bioregions))
  
  # Store in final list
  nodes_metadata_per_sliding_windows <- append(x = nodes_metadata_per_sliding_windows, values = list(nodes_metadata_i))
}
names(nodes_metadata_per_sliding_windows) <- paste0("Window_",1:length(start_time_list))

## Save node metadata
# saveRDS(object = nodes_metadata_per_sliding_windows, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
saveRDS(object = nodes_metadata_per_sliding_windows, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))


### 12.3/ Plot each sliding windows as a network ####

# https://kateto.net/network-visualization

all_dispersal_events_all_sliding_windows_df_list <- list()

# Load node metadata for all sliding windows
# nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))

# Load ggplot df of mean counts of dispersal events between bioregions per sliding windows
# all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Adjust minimal number of events to be displayed for each sliding windows
min_counts_threshold_list <- rep(1, length(nodes_metadata_per_sliding_windows))

## Loop per sliding windows
for (i in seq_along(nodes_metadata_per_sliding_windows))
{
  # i <- 1
  
  # Extract node metadata for the given time stratum
  nodes_metadata_sliding_window_i <- nodes_metadata_per_sliding_windows[[i]] %>%
    select(node_ID, bioregions, time, mean_counts, mean_percentages, node_size, node_labels, latitude, longitude, adjusted_latitude, adjusted_longitude)
  
  # Extract mean time of the sliding window
  mean_time <- mean(nodes_metadata_sliding_window_i$time)
  start_time <- mean_time + window_width/2
  end_time <- mean_time - window_width/2
  
  ## 12.3.1/ Convert to igraph object ####
  
  # Extract vertice metadata for the given time stratum
  # Round counts and convert to edge size
  all_dispersal_events_sliding_window_i_df <- all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot %>% 
    filter(time == mean_time) %>% 
    mutate(mean_counts = round(mean_counts, 0)) %>% 
    mutate(mean_cum_counts = round(mean_cum_counts, 0)) %>% 
    mutate(mean_perc = round(100 * mean_counts / sum(mean_counts), 1)) %>% 
    mutate(edge_width = log1p(mean_counts)) # Adjust edge width to a log scale
  # mutate(edge_width = mean_counts / max(mean_counts)) # Standardize between 0 and 1
  
  ## Store egde metadata across all sliding window
  all_dispersal_events_all_sliding_windows_df_list <- append(x = all_dispersal_events_all_sliding_windows_df_list, values = list(all_dispersal_events_sliding_window_i_df))
  
  ## Convert to igraph
  all_dispersal_events_sliding_window_i_igraph <- igraph::graph_from_data_frame(d = all_dispersal_events_sliding_window_i_df, vertices = nodes_metadata_sliding_window_i, directed = T)
  # str(all_dispersal_events_sliding_window_i_igraph)
  
  # Explore igraph
  # all_dispersal_events_sliding_window_i_igraph
  E(all_dispersal_events_sliding_window_i_igraph) # Edge attributes
  V(all_dispersal_events_sliding_window_i_igraph) # Vertex/Node attributes
  all_dispersal_events_sliding_window_i_igraph[]  # Adjacency matrix of nodes x nodes
  
  # Export adjacency matrix
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_counts")
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_perc")
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_cum_counts")
  
  ## 12.3.2/ Adjust aesthetics ####
  
  # Set color scheme for areas/bioregions (Use the BSM color scheme)
  colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
  # colors_list_for_areas <- colors_list_for_states[nodes_metadata_sliding_window_i$node_ID]
  colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
  colors_list_for_areas <- colors_list_for_areas[nodes_metadata_sliding_window_i$node_ID]
  
  # Remove loops
  all_dispersal_events_sliding_window_i_igraph <- simplify(graph = all_dispersal_events_sliding_window_i_igraph,
                                                           remove.multiple = F, remove.loops = T) 
  
  # Set node colors based on bioregions:
  V(all_dispersal_events_sliding_window_i_igraph)$color <- colors_list_for_areas
  # Set node size based on current richness
  V(all_dispersal_events_sliding_window_i_igraph)$size <- V(all_dispersal_events_sliding_window_i_igraph)$node_size * 5
  V(all_dispersal_events_sliding_window_i_igraph)$size <- sapply(X = V(all_dispersal_events_sliding_window_i_igraph)$size, FUN = function (x) { max(x, 10) })
  # Set node font size based on current richness
  V(all_dispersal_events_sliding_window_i_igraph)$label.cex <- V(all_dispersal_events_sliding_window_i_igraph)$node_size / 5
  V(all_dispersal_events_sliding_window_i_igraph)$label.cex <- sapply(X = V(all_dispersal_events_sliding_window_i_igraph)$label.cex, FUN = function (x) { max(x, 0.5) })
  
  # Set edge width based on event counts
  E(all_dispersal_events_sliding_window_i_igraph)$width <- E(all_dispersal_events_sliding_window_i_igraph)$edge_width * 2 + 1
  # Set edge arrow size based on event counts
  E(all_dispersal_events_sliding_window_i_igraph)$arrow.size <- E(all_dispersal_events_sliding_window_i_igraph)$edge_width / 4 + 0.1
  
  # Set edge label as mean counts
  # E(all_dispersal_events_sliding_window_i_igraph)$edge.label <- E(all_dispersal_events_sliding_window_i_igraph)$mean_counts
  
  # Set edge color according to their source
  edge_starts <- ends(all_dispersal_events_sliding_window_i_igraph, es = E(all_dispersal_events_sliding_window_i_igraph), names = F)[, 1]
  edge_col <- V(all_dispersal_events_sliding_window_i_igraph)$color[edge_starts]
  E(all_dispersal_events_sliding_window_i_igraph)$color <- edge_col
  
  # Remove edges with low counts
  table(E(all_dispersal_events_sliding_window_i_igraph)$mean_counts)
  # cutoff_min_counts_i <- cutoff_min_counts
  cutoff_min_counts_i <- min_counts_threshold_list[i]
  all_dispersal_events_sliding_window_i_igraph <- delete_edges(graph = all_dispersal_events_sliding_window_i_igraph,
                                                        edges = E(all_dispersal_events_sliding_window_i_igraph)[mean_counts < cutoff_min_counts_i])
  
  # Adjust layout to coordinates of vertex/nodes
  layout_WGS84 <- cbind(V(all_dispersal_events_sliding_window_i_igraph)$longitude, V(all_dispersal_events_sliding_window_i_igraph)$latitude)
  layout_pretty <- cbind(V(all_dispersal_events_sliding_window_i_igraph)$adjusted_longitude, V(all_dispersal_events_sliding_window_i_igraph)$adjusted_latitude)
  
  ## Save igraph for all dispersal events across all sliding windows
  saveRDS(object = all_dispersal_events_sliding_window_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_window_",i,"_",window_width,"My_igraph.rds"))
  saveRDS(object = all_dispersal_events_sliding_window_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_window_",i,"_",window_width,"My_igraph.rds"))
  
  ## Plot igraph
  
  # pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_windows_",i,"_",window_width,"My_igraph.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_windows_",i,"_",window_width,"My_igraph.pdf"),
      width = 6, height = 6)
  
  plot.igraph(x = all_dispersal_events_sliding_window_i_igraph,
              ## Node aesthetics
              # vertex.color = colors_list_for_areas, # Nodes color
              vertex.frame.color = "black", # Nodes border color
              vertex.shape = "circle",  # Node shape
              # vertex.size = 15, # Node size
              vertex.label	= nodes_metadata_sliding_window_i$node_labels, # Node labels
              vertex.label.color = "black",  # Node label color
              vertex.label.font = 2, # Node label font
              # vertex.label.cex = 1, # Node label cex
              ## Edge aesthetics
              # edge.color = "black",     # Edge color
              # edge.width = 1,         # Edge width
              # edge.arrow.size = 0.7,    # Terminal arrow size
              edge.label = E(all_dispersal_events_sliding_window_i_igraph)$mean_counts,
              edge.label.cex = 0.7,
              edge.lty = 1,             # Edge line type
              edge.curved = 0.1,
              arrow.mode = "forward",    # Direction of arrows
              ## Other aesthetics
              # layout = layout_WGS84,
              layout = layout_pretty,
              main = paste0("Dispersal events between bioregions\n",start_time,"-",end_time," My")
  )
  
  dev.off()
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Igraph plotted for ",start_time,"-",end_time," My - Window n°", i, "/", length(nodes_metadata_per_sliding_windows),"\n"))
  }
  
}

## Save egde metadata across all sliding_window
names(all_dispersal_events_all_sliding_windows_df_list) <- paste0("Window_",1:length(all_dispersal_events_all_sliding_windows_df_list))
# saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))


## 12.3.3/ Aggregate all igraph plots in a single pdf ####

## With overlap (all windows)

# all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_sliding_windows_", full.names = T)
all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_sliding_windows_", full.names = T)
all_igraphs_path <- all_igraphs_path[str_detect(string = all_igraphs_path, pattern = paste0(window_width,"My"))]
all_igraphs_path <- all_igraphs_path[str_detect(string = all_igraphs_path, pattern = "igraph.pdf")]
nb_igraphs <- length(all_igraphs_path)

# Reorder paths in numerical order
indices <- as.character(1:nb_igraphs)
all_igraphs_path_prefix <- str_remove(string = all_igraphs_path, pattern = "_sliding_windows_.*")
all_igraphs_path_suffix <- paste0("_",window_width,"My_igraph.pdf")
all_igraphs_path_reordered <- paste0(all_igraphs_path_prefix, "_sliding_windows_", indices, all_igraphs_path_suffix)

# qpdf::pdf_combine(input = all_igraphs_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph.pdf"))
qpdf::pdf_combine(input = all_igraphs_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph.pdf"))

## Without overlap

window_coverage <- window_width / window_steps
indices_no_overlap <- seq(from = 1, to = nb_igraphs, by = (window_coverage + 1))

all_igraphs_path_reordered_no_overlap <- paste0(all_igraphs_path_prefix[1], "_sliding_windows_", indices_no_overlap, all_igraphs_path_suffix[1])

# qpdf::pdf_combine(input = all_igraphs_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_no_overlap.pdf"))
qpdf::pdf_combine(input = all_igraphs_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_no_overlap.pdf"))

## 12.3.4/ Aggregate all igraph in forward timeline + overall ####

## With overlap (all windows)

# overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"
overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"

# qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_forward.pdf"))

## Without overlap

# qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered_no_overlap), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_forward_no_overlap.pdf"))
qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered_no_overlap), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_igraph_forward_no_overlap.pdf"))


### 12.4/ Map networks over bioregion maps ####

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

## Load overall node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

## Load node metadata per sliding windows
# nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))

## Load overall edge metadata
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")

# Load edge metadata per sliding windows
# all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[nodes_metadata_per_time_strata[[1]]$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
# colors_list_for_areas <- colors_list_for_areas[nodes_metadata_per_sliding_windows[[1]]$node_ID]
# names(colors_list_for_areas) <- nodes_metadata_per_sliding_windows[[1]]$bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Adjust order of bioregions
for (i in 1:length(nodes_metadata_per_sliding_windows))
{
  nodes_metadata_per_sliding_windows[[i]]$bioregions <- factor(nodes_metadata_per_sliding_windows[[i]]$bioregions, levels = bioregion_names, labels = bioregion_names)
}
Bioregions_sf_Bioregions_level$Bioregion <- factor(Bioregions_sf_Bioregions_level$Bioregion, levels = bioregion_names, labels = bioregion_names)

# Convert Long/Lat WGS84 to Mollweide
for (i in 1:length(nodes_metadata_per_sliding_windows))
{
  nodes_metadata_sf <- st_as_sf(x = nodes_metadata_per_sliding_windows[[i]], coords = c("longitude", "latitude"), remove = F, crs = st_crs(Bioregions_sf_Bioregions_level))
  nodes_metadata_sf <- st_transform(nodes_metadata_sf, crs = st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  Mollweide_coordinates <- st_as_text(nodes_metadata_sf$geometry)
  longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  nodes_metadata_per_sliding_windows[[i]]$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
  latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
  nodes_metadata_per_sliding_windows[[i]]$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))
}

## 12.4.1/ Inform edge metadata with coordinates and bioregion source labels ####

for (i in 1:length(all_dispersal_events_all_sliding_windows_df_list))
{
  # i <- 1
  
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_all_sliding_windows_df_list[[i]]
  nodes_metadata_i <- nodes_metadata_per_sliding_windows[[i]]
  
  all_dispersal_events_sliding_windows_i_df <- left_join(x = all_dispersal_events_sliding_windows_i_df,
                                                y = nodes_metadata_i[, c("node_ID", "bioregions", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                                by = join_by(source == node_ID))
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_sliding_windows_i_df %>% 
    rename(source_labels = bioregions,
           source_latitude = latitude,
           source_longitude = longitude,
           source_latitude_Mollweide = latitude_Mollweide,
           source_longitude_Mollweide = longitude_Mollweide)
  all_dispersal_events_sliding_windows_i_df <- left_join(x = all_dispersal_events_sliding_windows_i_df,
                                                y = nodes_metadata_i[, c("node_ID", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                                by = join_by(dest == node_ID))
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_sliding_windows_i_df %>% 
    rename(dest_latitude = latitude,
           dest_longitude = longitude,
           dest_latitude_Mollweide = latitude_Mollweide,
           dest_longitude_Mollweide = longitude_Mollweide)
  
  # Adjust order of bioregions
  all_dispersal_events_sliding_windows_i_df$source_labels <- factor(all_dispersal_events_sliding_windows_i_df$source_labels, levels = bioregion_names, labels = bioregion_names)
  
  # Store updated node metadata df
  all_dispersal_events_all_sliding_windows_df_list[[i]] <- all_dispersal_events_sliding_windows_i_df
}

# Save edge metadata df
# saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))


# Load edge metadata df
# all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))


# May adjust curvature manually
# May adjust arrow size too

# Set minimum nb of dispersal events to display
# min_counts_threshold <- 1
min_counts_threshold_list <- rep(1, length(all_dispersal_events_all_sliding_windows_df_list))

max_node_counts <- max(nodes_metadata$mean_counts)
max_edge_counts <- max(all_dispersal_events_overall_df$mean_counts)

# Set breaks for edge linewidth
# breaks_linewidth_list <- list(Stratum_1 = c(6, 20, 45), Stratum_2 = c(5, 15, 30), Stratum_3 = c(1, 2, 3), Stratum_4 = c(1, 2, 5), Stratum_5 = c(1, 2), Stratum_6 = c(1, 2), Stratum_7 = c(1))
# breaks_linewidth_list <- rep(list(c(1, 1, 1)), length(all_dispersal_events_all_sliding_windows_df_list))

# Function to adjust range size of nodes based on the range scale used for overall data
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


## 12.4.2/ Plot ggplot for each sliding window ####

for (i in seq_along(all_dispersal_events_all_sliding_windows_df_list))
{
  # i <- 1
  # i <- 94
  
  # Extract the minimum number of events to display an edge
  # min_counts_threshold_i <- min_counts_threshold
  min_counts_threshold_i <- min_counts_threshold_list[i]
  
  # Extract node metadata
  nodes_metadata_sliding_windows_i <- nodes_metadata_per_sliding_windows[[i]]
  # Extract max node size to adjust range of points
  max_node_size <- max(nodes_metadata_sliding_windows_i$mean_counts)
  max_range_size <- rescale_range_size(x = max_node_size)
  # Extract min node size to adjust range of points
  min_node_size <- min(nodes_metadata_sliding_windows_i$mean_counts)
  min_range_size <- rescale_range_size(x = min_node_size)
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract mean time of the sliding window
  mean_time <- mean(nodes_metadata_sliding_windows_i$time)
  start_time <- mean_time + window_width/2
  end_time <- mean_time - window_width/2
  
  # Extract edge metadata
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_all_sliding_windows_df_list[[i]]
  # Extract max node size to adjust range of points
  max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_counts)
  max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  # Extract min node size to adjust range of points
  min_edge_linewidth <- min(all_dispersal_events_sliding_windows_i_df$mean_counts[all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i])
  min_range_linewidth <- rescale_range_linewidth(x = min_edge_linewidth)
  # Extract manual breaks for edge linewidth per stratum
  # breaks_linewidth_i <- breaks_linewidth_list[[i]]
  breaks_linewidth_i <- c(max_edge_linewidth, max_edge_linewidth, max_edge_linewidth) # Provide 3 valid entries that will be replaced in the legend
  
  # If min counts of edges = max counts of edges, need a fix by adding 0.001 to an edge count
  if (min_edge_linewidth == max_edge_linewidth)
  {
    min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_counts == min_edge_linewidth)
    edge_fix_index <- sample(x = min_egde_indices, size = 1)
    all_dispersal_events_sliding_windows_i_df$mean_counts[edge_fix_index] <- all_dispersal_events_sliding_windows_i_df$mean_counts[edge_fix_index] + 0.001
    max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_counts)
    max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  }
  
  # If number of edges to display is 1, create duplicates of that edge, and add 0.001 to the edge count
  min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i)
  if (length(min_egde_indices) == 1)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[min_egde_indices, ]
    fix_df$mean_counts <- fix_df$mean_counts + 0.001
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(fix_df$mean_counts, fix_df$mean_counts, fix_df$mean_counts)
  }
    # If number of egde to display is null, created two loop edges with different edge values and set linewidth range to zero
  if (length(min_egde_indices) == 0)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[c(1,9), ]
    fix_df$mean_counts <- c(1, 1.001)
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(1, 1, 1)
  }
  
  # Plot PDF
  # pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_windows_",i,"_",window_width,"My_ggplot.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_sliding_windows_",i,"_",window_width,"My_ggplot.pdf"),
      width = 9, height = 4)
  
  all_dispersal_events_sliding_windows_i_ggplot <- ggplot(data = nodes_metadata_sliding_windows_i) +
    
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
    geom_point(mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide, size = mean_counts),
               alpha = 0.3, show.legend = F) +
    
    # Adjust size legend
    scale_size_continuous("Species richness", 
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +
    
    # Plot vertices
    geom_curve(data = all_dispersal_events_sliding_windows_i_df[all_dispersal_events_sliding_windows_i_df$mean_counts >= min_counts_threshold_i, ],
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
                               labels = c("1", "15", "40"),
                               # range = c(5, 30),
                               range = c(min_range_linewidth, max_range_linewidth)) +
    
    # Add node labels
    ggnewscale::new_scale(new_aes = "size") +
    geom_text(data = nodes_metadata_sliding_windows_i,
              mapping = aes(label = round(mean_counts, 0),
                            x = longitude_Mollweide, y = latitude_Mollweide,
                            size = mean_counts),
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +
    
    # Adjust size for labels
    scale_size_continuous("Species richness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
    
    # Add title
    ggtitle(label =  paste0("Dispersal events between bioregions\n",start_time,"-",end_time," My")) +
    
    # Adjust legend aesthetics
    guides(color = "none",
           size = "none",
           fill = guide_legend(order = 1),
           linewidth = guide_legend(order = 2,
                                    override.aes = list(linewidth = c(0.77, 1.55, 2.95), # Values for 1, 15, 40
                                                        arrow(length = unit(0.2, "cm"))))
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
  
  print(all_dispersal_events_sliding_windows_i_ggplot)
  
  dev.off()
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Map created for ",start_time,"-",end_time," My - Window n°", i, "/", length(nodes_metadata_per_sliding_windows),"\n"))
  }
  
}

## 12.4.3/ Aggregate all ggplot plots in a single pdf ####

## With overlap (all windows)

# all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_sliding_windows_", full.names = T)
all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_sliding_windows_", full.names = T)

all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = paste0(window_width,"My"))]
all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = "ggplot.pdf")]
nb_ggplots <- length(all_ggplots_path)

# Reorder paths in numerical order
indices <- as.character(1:nb_ggplots)
all_ggplots_path_prefix <- str_remove(string = all_ggplots_path, pattern = "_sliding_windows_.*")
all_ggplots_path_suffix <- paste0("_",window_width,"My_ggplot.pdf")
all_ggplots_path_reordered <- paste0(all_ggplots_path_prefix, "_sliding_windows_", indices, all_ggplots_path_suffix)

# qpdf::pdf_combine(input = all_ggplots_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot.pdf"))
qpdf::pdf_combine(input = all_ggplots_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot.pdf"))

## Without overlap

window_coverage <- window_width / window_steps
indices_no_overlap <- seq(from = 1, to = nb_ggplots, by = (window_coverage + 1))

all_ggplots_path_reordered_no_overlap <- paste0(all_ggplots_path_prefix[1], "_sliding_windows_", indices_no_overlap, all_ggplots_path_suffix[1])

# qpdf::pdf_combine(input = all_ggplots_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_no_overlap.pdf"))
qpdf::pdf_combine(input = all_ggplots_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_no_overlap.pdf"))


## 12.4.4/ Aggregate all ggplot in forward timeline + overall ####

## With overlap (all windows)

# overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"
overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"

# qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"))

## Without overlap

# qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered_no_overlap), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward_no_overlap.pdf"))
qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered_no_overlap), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward_no_overlap.pdf"))


### 12.5/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

window_width <- 5

# pdf_pointer_mean_counts <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"),
pdf_pointer_mean_counts <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"),
                                                  pages = NULL, density = 100)
magick::image_info(pdf_pointer_mean_counts)

image_resize_and_write_gif(image = pdf_pointer_mean_counts,
                           # path =  paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.gif"),
                           path =  paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_",window_width,"My_ggplot_forward.gif"),
                           delay = 1/10, # Time between frames in seconds
                           width = 900, height = 400,
                           loop = FALSE,
                           progress = TRUE)


## Add arrow to legend (Does not work with guides. Need to add on PPT)
## Add bg color on sphere (does not clip out of the sphere... Looks better without grid and bg)

# Add Color contour in igraphs on areas that are adjacent at a given time: Section 4.3 and 4.4
# Adjust coordinates in each stratum to show continental drift?

# Create bioregion shapefile that can overlay the PALEOMAPS?

### Make an animation by plotting every sliding window?
# Use gganimate for paths to make them fade in_out?




##### 13/ Map networks of cumulative counts of dispersal events along sliding windows ####

# Cumulative series => Use non overlapping windows of 1My

window_width <- 1

### 13.1/ Aggregate counts for all dispersal events, per sliding windows ####

# Load ggplot df
# all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/all_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

names(all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot)

# Extract counts only for targeted type of events
dispersal_event_types_list <- c("range extension (d)", "jump-dispersal (j)")
# Only range extension (d) and jump dispersal (j) events

# Sum across event_types, per sliding windows
# Aggregate mean counts across maps
all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- all_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot %>% 
  filter(event_type %in% dispersal_event_types_list) %>% # Extract counts only for targeted type of events
  group_by(source, dest, map, time) %>% 
  summarize(counts = sum(counts)) %>% # Sum counts across event_types
  arrange(desc(time)) %>% 
  group_by(source, dest, map) %>%
  mutate(cum_counts = cumsum(counts)) %>%  # Compute cumsum along time
  group_by(source, dest, time) %>%
  summarize(mean_counts = mean(counts),                # Aggregate mean counts across maps
            mean_cum_counts = mean(cum_counts)) %>%    # Aggregate mean sum counts across maps
  group_by(time) %>% 
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts)) %>%  # Compute overall percentages for a given time
  ungroup()


# Save ggplot df of mean counts of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))

## Convert back to array

# Pivot mean counts and percentages
all_dispersal_events_mean_count_per_maps_all_sliding_windows_melted <- all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot %>% 
  pivot_longer(cols = c("mean_counts", "mean_perc", "mean_cum_counts"), names_to = "mean_stats", values_to = "values")

all_dispersal_events_mean_count_per_maps_all_sliding_windows_array <- reshape2::acast(data = all_dispersal_events_mean_count_per_maps_all_sliding_windows_melted,
                                                                                      formula = source ~ dest ~ time ~ mean_stats)

dimnames(all_dispersal_events_mean_count_per_maps_all_sliding_windows_array)
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_counts"]
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_perc"]
all_dispersal_events_mean_count_per_maps_all_sliding_windows_array[ , , ,"mean_cum_counts"]

# Save array of mean counts/perc of dispersal events between bioregions
# saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_array.rds"))
saveRDS(object = all_dispersal_events_mean_count_per_maps_all_sliding_windows_array, file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_array.rds"))


### 13.2/ Compute node metadata per sliding windows ####

# Extract root age
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/DEC_J_BSM_output.rds")
# DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
root_age <- max(DEC_J_clado_events_tables[[1]]$time_bp)

# Generate start and end times of sliding windows
start_time_list <- seq(from = 0, to = root_age - window_width, by = window_steps)
end_time_list <- seq(from = window_width, to = root_age, by = window_steps)
mean_time_list <- (start_time_list + end_time_list) / 2

# Load LTT data for species richness in time
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

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

## Create node metadata with species richness at mean time of sliding windows
nodes_metadata_per_sliding_windows <- list()
for (i in 1:length(start_time_list))
{
  # i <- 1
  
  # Extract mean time boundaries of sliding windows
  mean_time <- mean_time_list[i]
  
  # Aggregate species richness within each time stratum
  nodes_metadata_i <- DEC_J_LTT_all_areas_mean_ggplot %>% 
    filter(time == mean_time) %>% 
    filter(areas != "total") %>% 
    mutate(mean_percentages = 100 * mean_counts / sum(mean_counts)) %>%
    mutate(bioregions = areas) %>%
    rename(node_ID = areas) %>%
    mutate(node_size = log(mean_counts)) %>%
    mutate(window = paste0("Window_",i)) %>%
    mutate(time = mean_time) %>%
    select(window, time, node_ID, mean_counts, mean_percentages, node_size)
  
  # Add bioregion labels
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_names_alpha_order_df, by = join_by(node_ID == areas)) %>%
    select(window, time, node_ID, bioregions, mean_counts, mean_percentages, node_size)
  
  # Add label with richness
  nodes_metadata_i$node_labels <- paste0(nodes_metadata_i$node_ID, "\n", round(nodes_metadata_i$mean_counts, 0))
  
  # Reorder in alphabetic order of node ID
  nodes_metadata_i <- nodes_metadata_i %>% 
    arrange(node_ID)
  
  ## Add spatial coordinates for each bioregion (May adjust to show continental drift)
  nodes_metadata_i <- left_join(x = nodes_metadata_i, y = bioregion_coordinates_df,
                                by = join_by(node_ID == areas, bioregions == bioregions))
  
  # Store in final list
  nodes_metadata_per_sliding_windows <- append(x = nodes_metadata_per_sliding_windows, values = list(nodes_metadata_i))
}
names(nodes_metadata_per_sliding_windows) <- paste0("Window_",1:length(start_time_list))

## Save node metadata
# saveRDS(object = nodes_metadata_per_sliding_windows, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
saveRDS(object = nodes_metadata_per_sliding_windows, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))


### 13.3/ Plot each sliding windows as a network ####

# https://kateto.net/network-visualization

all_dispersal_events_all_sliding_windows_df_list <- list()

# Load node metadata for all sliding windows
# nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))

# Load ggplot df of mean counts of dispersal events between bioregions per sliding windows
# all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Dispersal_events/all_dispersal_events_mean_count_per_maps_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Adjust minimal number of events to be displayed for each sliding windows
min_counts_threshold_list <- rep(1, length(nodes_metadata_per_sliding_windows))

## Loop per sliding windows
for (i in seq_along(nodes_metadata_per_sliding_windows))
{
  # i <- 1
  
  # Extract node metadata for the given time stratum
  nodes_metadata_sliding_window_i <- nodes_metadata_per_sliding_windows[[i]] %>%
    select(node_ID, bioregions, time, mean_counts, mean_percentages, node_size, node_labels, latitude, longitude, adjusted_latitude, adjusted_longitude)
  
  # Extract mean time of the sliding window
  mean_time <- mean(nodes_metadata_sliding_window_i$time)
  start_time <- mean_time + window_width/2
  end_time <- mean_time - window_width/2
  
  ## 13.3.1/ Convert to igraph object ####
  
  # Extract vertice metadata for the given time stratum
  # Round counts and convert to edge size
  all_dispersal_events_sliding_window_i_df <- all_dispersal_events_mean_count_per_maps_all_sliding_windows_ggplot %>% 
    filter(time == mean_time) %>% 
    mutate(mean_counts = round(mean_counts, 0)) %>% 
    mutate(mean_cum_counts = round(mean_cum_counts, 0)) %>% 
    mutate(mean_perc = round(100 * mean_counts / sum(mean_counts), 1)) %>% 
    mutate(edge_width = log1p(mean_cum_counts)) # Adjust edge width to a log scale of cumulative counts
  # mutate(edge_width = mean_counts / max(mean_cum_counts)) # Standardize between 0 and 1
  
  ## Store egde metadata across all sliding window
  all_dispersal_events_all_sliding_windows_df_list <- append(x = all_dispersal_events_all_sliding_windows_df_list, values = list(all_dispersal_events_sliding_window_i_df))
  
  ## Convert to igraph
  all_dispersal_events_sliding_window_i_igraph <- igraph::graph_from_data_frame(d = all_dispersal_events_sliding_window_i_df, vertices = nodes_metadata_sliding_window_i, directed = T)
  # str(all_dispersal_events_sliding_window_i_igraph)
  
  # Explore igraph
  # all_dispersal_events_sliding_window_i_igraph
  E(all_dispersal_events_sliding_window_i_igraph) # Edge attributes
  V(all_dispersal_events_sliding_window_i_igraph) # Vertex/Node attributes
  all_dispersal_events_sliding_window_i_igraph[]  # Adjacency matrix of nodes x nodes
  
  # Export adjacency matrix
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_counts")
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_perc")
  as_adjacency_matrix(all_dispersal_events_sliding_window_i_igraph, attr = "mean_cum_counts")
  
  ## 13.3.2/ Adjust aesthetics ####
  
  # Set color scheme for areas/bioregions (Use the BSM color scheme)
  colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
  # colors_list_for_areas <- colors_list_for_states[nodes_metadata_sliding_window_i$node_ID]
  colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
  colors_list_for_areas <- colors_list_for_areas[nodes_metadata_sliding_window_i$node_ID]
  
  # Remove loops
  all_dispersal_events_sliding_window_i_igraph <- simplify(graph = all_dispersal_events_sliding_window_i_igraph,
                                                           remove.multiple = F, remove.loops = T) 
  
  # Set node colors based on bioregions:
  V(all_dispersal_events_sliding_window_i_igraph)$color <- colors_list_for_areas
  # Set node size based on current richness
  V(all_dispersal_events_sliding_window_i_igraph)$size <- V(all_dispersal_events_sliding_window_i_igraph)$node_size * 5
  V(all_dispersal_events_sliding_window_i_igraph)$size <- sapply(X = V(all_dispersal_events_sliding_window_i_igraph)$size, FUN = function (x) { max(x, 10) })
  # Set node font size based on current richness
  V(all_dispersal_events_sliding_window_i_igraph)$label.cex <- V(all_dispersal_events_sliding_window_i_igraph)$node_size / 5
  V(all_dispersal_events_sliding_window_i_igraph)$label.cex <- sapply(X = V(all_dispersal_events_sliding_window_i_igraph)$label.cex, FUN = function (x) { max(x, 0.5) })
  
  # Set edge width based on event cumulative counts
  E(all_dispersal_events_sliding_window_i_igraph)$width <- E(all_dispersal_events_sliding_window_i_igraph)$edge_width * 2 + 1
  # Set edge arrow size based on event cumulative counts
  E(all_dispersal_events_sliding_window_i_igraph)$arrow.size <- E(all_dispersal_events_sliding_window_i_igraph)$edge_width / 4 + 0.1
  
  # Set edge label as mean cumulative counts
  # E(all_dispersal_events_sliding_window_i_igraph)$edge.label <- E(all_dispersal_events_sliding_window_i_igraph)$mean_cum_counts
  
  # Set edge color according to their source
  edge_starts <- ends(all_dispersal_events_sliding_window_i_igraph, es = E(all_dispersal_events_sliding_window_i_igraph), names = F)[, 1]
  edge_col <- V(all_dispersal_events_sliding_window_i_igraph)$color[edge_starts]
  E(all_dispersal_events_sliding_window_i_igraph)$color <- edge_col
  
  # Remove edges with low cumulative counts
  table(E(all_dispersal_events_sliding_window_i_igraph)$mean_cum_counts)
  # cutoff_min_counts_i <- cutoff_min_counts
  cutoff_min_counts_i <- min_counts_threshold_list[i]
  all_dispersal_events_sliding_window_i_igraph <- delete_edges(graph = all_dispersal_events_sliding_window_i_igraph,
                                                               edges = E(all_dispersal_events_sliding_window_i_igraph)[mean_cum_counts < cutoff_min_counts_i])
  
  # Adjust layout to coordinates of vertex/nodes
  layout_WGS84 <- cbind(V(all_dispersal_events_sliding_window_i_igraph)$longitude, V(all_dispersal_events_sliding_window_i_igraph)$latitude)
  layout_pretty <- cbind(V(all_dispersal_events_sliding_window_i_igraph)$adjusted_longitude, V(all_dispersal_events_sliding_window_i_igraph)$adjusted_latitude)
  
  ## Save igraph for all dispersal events across all sliding windows
  # saveRDS(object = all_dispersal_events_sliding_window_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_window_",i,"_",window_width,"My_igraph.rds"))
  saveRDS(object = all_dispersal_events_sliding_window_i_igraph, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_window_",i,"_",window_width,"My_igraph.rds"))
  
  ## Plot igraph

  # pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_windows_",i,"_",window_width,"My_igraph.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_windows_",i,"_",window_width,"My_igraph.pdf"),
      width = 6, height = 6)
  
  plot.igraph(x = all_dispersal_events_sliding_window_i_igraph,
              ## Node aesthetics
              # vertex.color = colors_list_for_areas, # Nodes color
              vertex.frame.color = "black", # Nodes border color
              vertex.shape = "circle",  # Node shape
              # vertex.size = 15, # Node size
              vertex.label	= nodes_metadata_sliding_window_i$node_labels, # Node labels
              vertex.label.color = "black",  # Node label color
              vertex.label.font = 2, # Node label font
              # vertex.label.cex = 1, # Node label cex
              ## Edge aesthetics
              # edge.color = "black",     # Edge color
              # edge.width = 1,         # Edge width
              # edge.arrow.size = 0.7,    # Terminal arrow size
              edge.label = E(all_dispersal_events_sliding_window_i_igraph)$mean_cum_counts, # Use cumulative counts as labels
              edge.label.cex = 0.7,
              edge.lty = 1,             # Edge line type
              edge.curved = 0.1,
              arrow.mode = "forward",    # Direction of arrows
              ## Other aesthetics
              # layout = layout_WGS84,
              layout = layout_pretty,
              main = paste0("Dispersal events between bioregions\nTotal until ",mean_time + 0.5," My") # Use mean time to display time of the window
  )
  
  dev.off()
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Igraph of cumulative counts plotted for ",mean_time + 0.5," My - Window n°", i, "/", length(nodes_metadata_per_sliding_windows),"\n"))
  }
  
}

## Save egde metadata across all sliding window for cumulative counts
names(all_dispersal_events_all_sliding_windows_df_list) <- paste0("Window_",1:length(all_dispersal_events_all_sliding_windows_df_list))
# saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_df_list_",window_width,"My.rds"))
saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_df_list_",window_width,"My.rds"))


## 13.3.3/ Aggregate all igraph plots in a single pdf ####

## With overlap (all windows)

# all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_cum_counts_sliding_windows_", full.names = T)
all_igraphs_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_cum_counts_sliding_windows_", full.names = T)

all_igraphs_path <- all_igraphs_path[str_detect(string = all_igraphs_path, pattern = paste0(window_width,"My"))]
all_igraphs_path <- all_igraphs_path[str_detect(string = all_igraphs_path, pattern = "igraph.pdf")]
nb_igraphs <- length(all_igraphs_path)

# Reorder paths in numerical order
indices <- as.character(1:nb_igraphs)
all_igraphs_path_prefix <- str_remove(string = all_igraphs_path, pattern = "_sliding_windows_.*")
all_igraphs_path_suffix <- paste0("_",window_width,"My_igraph.pdf")
all_igraphs_path_reordered <- paste0(all_igraphs_path_prefix, "_sliding_windows_", indices, all_igraphs_path_suffix)

# qpdf::pdf_combine(input = all_igraphs_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph.pdf"))
qpdf::pdf_combine(input = all_igraphs_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph.pdf"))

## Without overlap

window_coverage <- window_width / window_steps
indices_no_overlap <- seq(from = 1, to = nb_igraphs, by = (window_coverage + 1))

all_igraphs_path_reordered_no_overlap <- paste0(all_igraphs_path_prefix[1], "_sliding_windows_", indices_no_overlap, all_igraphs_path_suffix[1])

# qpdf::pdf_combine(input = all_igraphs_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_no_overlap.pdf"))
qpdf::pdf_combine(input = all_igraphs_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_no_overlap.pdf"))

## 13.3.4/ Aggregate all igraph in forward timeline + overall ####

## With overlap (all windows)

# overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"
overall_igraph_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_igraph.pdf"

# qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_forward.pdf"))

## Without overlap

# qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered_no_overlap), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_forward_no_overlap.pdf"))
qpdf::pdf_combine(input = c(rev(all_igraphs_path_reordered_no_overlap), overall_igraph_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_igraph_forward_no_overlap.pdf"))


### 13.4/ Map networks over bioregion maps ####

## Load Bioregions map
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
# Remove Antarctica
Bioregions_sf_Bioregions_level <- Bioregions_sf_Bioregions_level[Bioregions_sf_Bioregions_level$Bioregion != "Antarctica", ]

## Load overall node metadata
# nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata.rds")
nodes_metadata <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata.rds")

## Load node metadata per sliding windows
# nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))
nodes_metadata_per_sliding_windows <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/nodes_metadata_per_sliding_windows_",window_width,"My.rds"))

## Load overall edge metadata
# all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_df.rds")
all_dispersal_events_overall_df <- readRDS(file = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_df.rds")

# Load edge metadata of cumulative counts per sliding windows
# all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_df_list_",window_width,"My.rds"))
all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_df_list_",window_width,"My.rds"))

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[nodes_metadata_per_sliding_windows[[1]]$node_ID]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
# colors_list_for_areas <- colors_list_for_areas[nodes_metadata_per_sliding_windows[[1]]$node_ID]
# names(colors_list_for_areas) <- nodes_metadata_per_sliding_windows[[1]]$bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Adjust order of bioregions
for (i in 1:length(nodes_metadata_per_sliding_windows))
{
  nodes_metadata_per_sliding_windows[[i]]$bioregions <- factor(nodes_metadata_per_sliding_windows[[i]]$bioregions, levels = bioregion_names, labels = bioregion_names)
}
Bioregions_sf_Bioregions_level$Bioregion <- factor(Bioregions_sf_Bioregions_level$Bioregion, levels = bioregion_names, labels = bioregion_names)

# Convert Long/Lat WGS84 to Mollweide
for (i in 1:length(nodes_metadata_per_sliding_windows))
{
  nodes_metadata_sf <- st_as_sf(x = nodes_metadata_per_sliding_windows[[i]], coords = c("longitude", "latitude"), remove = F, crs = st_crs(Bioregions_sf_Bioregions_level))
  nodes_metadata_sf <- st_transform(nodes_metadata_sf, crs = st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  Mollweide_coordinates <- st_as_text(nodes_metadata_sf$geometry)
  longitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  nodes_metadata_per_sliding_windows[[i]]$longitude_Mollweide <- as.numeric(str_remove(string = longitude_Mollweide, pattern = " .*"))
  latitude_Mollweide <- str_remove(string = Mollweide_coordinates, pattern = ".* \\(")
  latitude_Mollweide <- str_remove(string = latitude_Mollweide, pattern = ".* ")
  nodes_metadata_per_sliding_windows[[i]]$latitude_Mollweide <- as.numeric(str_remove(string = latitude_Mollweide, pattern = "\\)"))
}

## 13.4.1/ Inform edge metadata with coordinates and bioregion source labels ####

for (i in 1:length(all_dispersal_events_all_sliding_windows_df_list))
{
  # i <- 1
  
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_all_sliding_windows_df_list[[i]]
  nodes_metadata_i <- nodes_metadata_per_sliding_windows[[i]]
  
  all_dispersal_events_sliding_windows_i_df <- left_join(x = all_dispersal_events_sliding_windows_i_df,
                                                         y = nodes_metadata_i[, c("node_ID", "bioregions", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                                         by = join_by(source == node_ID))
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_sliding_windows_i_df %>% 
    rename(source_labels = bioregions,
           source_latitude = latitude,
           source_longitude = longitude,
           source_latitude_Mollweide = latitude_Mollweide,
           source_longitude_Mollweide = longitude_Mollweide)
  all_dispersal_events_sliding_windows_i_df <- left_join(x = all_dispersal_events_sliding_windows_i_df,
                                                         y = nodes_metadata_i[, c("node_ID", "latitude", "longitude", "latitude_Mollweide", "longitude_Mollweide")],
                                                         by = join_by(dest == node_ID))
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_sliding_windows_i_df %>% 
    rename(dest_latitude = latitude,
           dest_longitude = longitude,
           dest_latitude_Mollweide = latitude_Mollweide,
           dest_longitude_Mollweide = longitude_Mollweide)
  
  # Adjust order of bioregions
  all_dispersal_events_sliding_windows_i_df$source_labels <- factor(all_dispersal_events_sliding_windows_i_df$source_labels, levels = bioregion_names, labels = bioregion_names)
  
  # Store updated node metadata df
  all_dispersal_events_all_sliding_windows_df_list[[i]] <- all_dispersal_events_sliding_windows_i_df
}

# Save edge metadata df
# saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
saveRDS(object = all_dispersal_events_all_sliding_windows_df_list, file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))


# Load edge metadata df
# all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))
all_dispersal_events_all_sliding_windows_df_list <- readRDS(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_all_sliding_windows_df_list_",window_width,"My.rds"))

# May adjust curvature manually
# May adjust arrow size too

# Set minimum nb of dispersal events to display
# min_counts_threshold <- 1
min_counts_threshold_list <- rep(1, length(all_dispersal_events_all_sliding_windows_df_list))

max_node_counts <- max(nodes_metadata$mean_counts)
max_edge_counts <- max(all_dispersal_events_overall_df$mean_counts)

# Set breaks for edge linewidth
# breaks_linewidth_list <- list(Stratum_1 = c(6, 20, 45), Stratum_2 = c(5, 15, 30), Stratum_3 = c(1, 2, 3), Stratum_4 = c(1, 2, 5), Stratum_5 = c(1, 2), Stratum_6 = c(1, 2), Stratum_7 = c(1))
# breaks_linewidth_list <- rep(list(c(1, 1, 1)), length(all_dispersal_events_all_sliding_windows_df_list))

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


## 13.4.2/ Plot ggplot for each sliding window ####

for (i in seq_along(all_dispersal_events_all_sliding_windows_df_list))
{
  # i <- 1
  # i <- 94
  
  # Extract the minimum number of events to display an edge
  # min_counts_threshold_i <- min_counts_threshold
  min_counts_threshold_i <- min_counts_threshold_list[i]
  
  # Extract node metadata
  nodes_metadata_sliding_windows_i <- nodes_metadata_per_sliding_windows[[i]]
  # Extract max node size to adjust range of points
  max_node_size <- max(nodes_metadata_sliding_windows_i$mean_counts)
  max_range_size <- rescale_range_size(x = max_node_size)
  # Extract min node size to adjust range of points
  min_node_size <- min(nodes_metadata_sliding_windows_i$mean_counts)
  min_range_size <- rescale_range_size(x = min_node_size)
  # Extract range size for labels
  min_range_label_size <- rescale_range_size(x = min_node_size, range_min = 3, range_max = 10)
  max_range_label_size <- rescale_range_size(x = max_node_size, range_min = 3, range_max = 10)
  
  # Extract mean time of the sliding window
  mean_time <- mean(nodes_metadata_sliding_windows_i$time)
  start_time <- mean_time + window_width/2
  end_time <- mean_time - window_width/2
  
  # Extract edge metadata
  all_dispersal_events_sliding_windows_i_df <- all_dispersal_events_all_sliding_windows_df_list[[i]]
  # Extract max node size to adjust range of linewidth
  max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_cum_counts)
  max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  # Extract min node size to adjust range of linewidth
  min_edge_linewidth <- min(all_dispersal_events_sliding_windows_i_df$mean_cum_counts[all_dispersal_events_sliding_windows_i_df$mean_cum_counts >= min_counts_threshold_i])
  min_range_linewidth <- rescale_range_linewidth(x = min_edge_linewidth)
  # Extract manual breaks for edge linewidth per stratum
  # breaks_linewidth_i <- breaks_linewidth_list[[i]]
  breaks_linewidth_i <- c(max_edge_linewidth, max_edge_linewidth, max_edge_linewidth) # Provide 3 valid entries that will be replaced in the legend
  
  # If min counts of edges = max counts of edges, need a fix by adding 0.001 to an edge count
  if (min_edge_linewidth == max_edge_linewidth)
  {
    min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_cum_counts == min_edge_linewidth)
    edge_fix_index <- sample(x = min_egde_indices, size = 1)
    all_dispersal_events_sliding_windows_i_df$mean_cum_counts[edge_fix_index] <- all_dispersal_events_sliding_windows_i_df$mean_cum_counts[edge_fix_index] + 0.001
    max_edge_linewidth <- max(all_dispersal_events_sliding_windows_i_df$mean_cum_counts)
    max_range_linewidth <- rescale_range_linewidth(x = max_edge_linewidth)
  }
  
  # If number of edges to display is 1, create duplicates of that edge, and add 0.001 to the edge count
  min_egde_indices <- which(all_dispersal_events_sliding_windows_i_df$mean_cum_counts >= min_counts_threshold_i)
  if (length(min_egde_indices) == 1)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[min_egde_indices, ]
    fix_df$mean_cum_counts <- fix_df$mean_cum_counts + 0.001
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(fix_df$mean_cum_counts, fix_df$mean_cum_counts, fix_df$mean_cum_counts)
  }
  # If number of egde to display is null, created two loop edges with different edge values and set linewidth range to zero
  if (length(min_egde_indices) == 0)
  {
    fix_df <- all_dispersal_events_sliding_windows_i_df[c(1,9), ]
    fix_df$mean_cum_counts <- c(1, 1.001)
    all_dispersal_events_sliding_windows_i_df <- rbind(all_dispersal_events_sliding_windows_i_df, fix_df)
    breaks_linewidth_i <- c(1, 1, 1)
  }
  
  # Plot PDF
  # pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_windows_",i,"_",window_width,"My_ggplot.pdf"),
  pdf(file = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/all_dispersal_events_cum_counts_sliding_windows_",i,"_",window_width,"My_ggplot.pdf"),
      width = 9, height = 4)
  
  all_dispersal_events_sliding_windows_i_ggplot <- ggplot(data = nodes_metadata_sliding_windows_i) +
    
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
    geom_point(mapping = aes(x = longitude_Mollweide, y = latitude_Mollweide, size = mean_counts),
               alpha = 0.3, show.legend = F) +
    
    # Adjust size legend
    scale_size_continuous("Species richness", 
                          # range = c(5, 30),
                          range = c(min_range_size, max_range_size)) +
    
    # Plot vertices
    geom_curve(data = all_dispersal_events_sliding_windows_i_df[all_dispersal_events_sliding_windows_i_df$mean_cum_counts >= min_counts_threshold_i, ],
               aes(x = source_longitude_Mollweide, y = source_latitude_Mollweide,
                   xend = dest_longitude_Mollweide, yend = dest_latitude_Mollweide,
                   color = source_labels, linewidth = mean_cum_counts),
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
                               labels = c("1", "15", "40"),
                               # range = c(5, 30),
                               range = c(min_range_linewidth, max_range_linewidth)) +
    
    # Add node labels
    ggnewscale::new_scale(new_aes = "size") +
    geom_text(data = nodes_metadata_sliding_windows_i,
              mapping = aes(label = round(mean_counts, 0),
                            x = longitude_Mollweide, y = latitude_Mollweide,
                            size = mean_counts),
              # color = rgb(t(col2rgb("black")), alpha = 200, maxColorValue = 255),
              color = "black",
              fontface = 2, show.legend = F) +
    
    # Adjust size for labels
    scale_size_continuous("Species richness",
                          # range = c(3, 10),
                          range = c(min_range_label_size, max_range_label_size)) +
    
    # Add title
    ggtitle(label =  paste0("Dispersal events between bioregions\nTotal until ",mean_time + 0.5," My")) +
    
    # Adjust legend aesthetics
    guides(color = "none",
           size = "none",
           fill = guide_legend(order = 1),
           linewidth = guide_legend(order = 2,
                                    override.aes = list(linewidth = c(0.77, 1.55, 2.95), # Values for 1, 15, 40
                                                        arrow(length = unit(0.2, "cm"))))
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
  
  print(all_dispersal_events_sliding_windows_i_ggplot)
  
  dev.off()
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Map created for ",mean_time + 0.5," My - Window n°", i, "/", length(nodes_metadata_per_sliding_windows),"\n"))
  }
  
}

## 13.4.3/ Aggregate all ggplot plots in a single pdf ####

## With overlap (all windows)

# all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_cum_counts_sliding_windows_", full.names = T)
all_ggplots_path <- list.files(path = "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/plots_for_sliding_windows/", pattern = "all_dispersal_events_cum_counts_sliding_windows_", full.names = T)

all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = paste0(window_width,"My"))]
all_ggplots_path <- all_ggplots_path[str_detect(string = all_ggplots_path, pattern = "ggplot.pdf")]
nb_ggplots <- length(all_ggplots_path)

# Reorder paths in numerical order
indices <- as.character(1:nb_ggplots)
all_ggplots_path_prefix <- str_remove(string = all_ggplots_path, pattern = "_sliding_windows_.*")
all_ggplots_path_suffix <- paste0("_",window_width,"My_ggplot.pdf")
all_ggplots_path_reordered <- paste0(all_ggplots_path_prefix, "_sliding_windows_", indices, all_ggplots_path_suffix)

# qpdf::pdf_combine(input = all_ggplots_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot.pdf"))
qpdf::pdf_combine(input = all_ggplots_path_reordered, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot.pdf"))

## Without overlap

window_coverage <- window_width / window_steps
indices_no_overlap <- seq(from = 1, to = nb_ggplots, by = (window_coverage + 1))

all_ggplots_path_reordered_no_overlap <- paste0(all_ggplots_path_prefix[1], "_sliding_windows_", indices_no_overlap, all_ggplots_path_suffix[1])

# qpdf::pdf_combine(input = all_ggplots_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_no_overlap.pdf"))
qpdf::pdf_combine(input = all_ggplots_path_reordered_no_overlap, output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_no_overlap.pdf"))


## 13.4.4/ Aggregate all ggplot in forward timeline + overall ####

## With overlap (all windows)

# overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"
overall_ggplot_path <- "./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_overall_ggplot.pdf"

# qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"))
qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"))

## Without overlap

# qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered_no_overlap), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward_no_overlap.pdf"))
qpdf::pdf_combine(input = c(rev(all_ggplots_path_reordered_no_overlap), overall_ggplot_path), output = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward_no_overlap.pdf"))


### 13.5/ Convert to GIF ####

source("./functions/image_resize_and_write_gif.R")

window_width <- 1
fps <- 5
# fps <- 10

# pdf_pointer_cum_counts <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"),
pdf_pointer_cum_counts <- magick::image_read_pdf(path = paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward.pdf"),
                                                    pages = NULL, density = 100)
magick::image_info(pdf_pointer_cum_counts)

image_resize_and_write_gif(image = pdf_pointer_cum_counts,
                           # path =  paste0("./outputs/Network_events_counts/Ponerinae_rough_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward_",fps,"fps.gif"),
                           path =  paste0("./outputs/Network_events_counts/Ponerinae_MCC_phylogeny_1534t/all_dispersal_events_cum_counts_all_sliding_windows_",window_width,"My_ggplot_forward_",fps,"fps.gif"),
                           delay = 1/fps, # Time between frames in seconds
                           width = 900, height = 400,
                           loop = FALSE,
                           progress = TRUE)



## Add arrow to legend (Does not work with guides. Need to add on PPT)
## Add bg color on sphere (does not clip out of the sphere... Looks better without grid and bg)

# Add Color contour in igraphs on areas that are adjacent at a given time: Section 4.3 and 4.4
# Adjust coordinates in each stratum to show continental drift?

# Create bioregion shapefile that can overlay the PALEOMAPS?

### Make an animation by plotting every sliding window?
# Use gganimate for paths to make them fade in_out?



## All can be segregated per clades!
  # Across clades to compare clades dynamics
  # Within clades, to zoom on a particular subclade of interest
     # Try merging probability of each bioregions membership on the ContMap by using weighted color in RGB space (do that overall and subset specific subclades for plotting)




