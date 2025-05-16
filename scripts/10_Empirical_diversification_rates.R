##### Script 10: Empirical diversification rates #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compute empirical diversification rates based on number of cladogenetic/in situ speciation events along sliding windows
# Plot evolution of diversification rates in time
  # Overall
  # Per bioregions

## Quick temporary solution to see trends, but true BD and SSE models are needed

###

### Inputs

# Cladogenetic/in situ speciation events counts along sliding windows
  # Use counts of Inheritance (y), Vicariance (v), Subset speciation (s), Jump-dispersal (j) for Cladogenetic events
  # Use only Inheritance (y) for in situ speciation events
  # Per source bioregions to compute bioregion rates

###

### Outputs

## Diversification rate plots
  # Mean counts of cladogenetic/in situ speciation events in time
  # Percentages across bioregions in time 
  # Evolution of diversification rates in time
     # Overall
     # Per bioregions

###


####### May wish to focus on specific clades to compute the same series of plots #######

# All can be segregated per clades!


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(phytools)
library(ape)
library(MultinomialCI)    # For 95% CIs on BSM counts
library(BayesTwin) # To compute HPD intervals
library(magick)   # For animated GIF

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Set sliding window parameters ####

window_width <- 10 # Width = 10 My (Must be a multiple of window_steps)
# window_width <- 5 # Width = 5 My (Must be a multiple of window_steps)

window_steps <- 1 # Steps = 1 My (Must be a divider of window_width)

# Extract root age
# DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_ggplot.rds")
DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_ggplot.rds")
# DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Youngest_phylogeny_1534t/DEC_J_LTT_all_areas_ggplot.rds")
# DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_Oldest_phylogeny_1534t/DEC_J_LTT_all_areas_ggplot.rds")

root_age <- max(DEC_J_LTT_all_areas_ggplot$time)

# Generate start and end times of sliding windows
tipward_time_list <- seq(from = 0, to = root_age - window_width, by = window_steps)
rootward_time_list <- seq(from = window_width, to = root_age, by = window_steps)
mean_time_list <- (rootward_time_list + tipward_time_list) / 2

### 1.3/ Load event counts data along sliding windows ####

## Load the melted data.frame of counts of events per source areas
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_all_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

names(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)

### 1.4/ Load residence times per areas along sliding windows #####

## Load the melted data.frame of residence times per areas across BSM maps
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))

### 1.5/ Set bioregion color scheme ####

# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))
# DEC_J_BSM_counts_arrays_all_sliding_windows <- readRDS(file = paste0("./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_counts_arrays_all_sliding_windows_",window_width,"My.rds"))

areas_list <- dimnames(DEC_J_BSM_counts_arrays_all_sliding_windows$dispersal_events_count_matrices_all_sliding_windows)[[1]]

# Set color scheme for areas/bioregions (Use the BSM color scheme)
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas

bioregion_names_with_total <- c(bioregion_names, "Total")
colors_list_for_areas_with_total <- c(colors_list_for_areas, "grey80")
names(colors_list_for_areas_with_total) <- bioregion_names_with_total

### 1.6/ Function for ggplot aesthetics ####

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



##### 2/ Compute cladogenetic events counts and rates ####

### 2.1/ Compute counts/rates per source bioregions ####

# Extract counts only for targeted type of events
cladogenetic_event_types_list <- c("range inheritance (y)", "vicariance (v)", "subset speciation (s)", "jump-dispersal (j)")

names(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)

# # Adjust order and names of Bioregions
# all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)

# Extract and sum counts of cladogenetic events
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- all_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  filter(event_type %in% cladogenetic_event_types_list) %>% 
  group_by(time, source, map) %>%
  summarize(counts = sum(counts), # Sum mean counts across cladogenetic event types
            area_richness = mean(area_richness),
            residence_time = mean(residence_time)) %>% 
  mutate(LTT_rates = counts / area_richness / window_width) %>% # Compute updated rates for all cladogenetic event types based on LTT data
  mutate(rates = counts / residence_time) %>% # Compute updated rates for all cladogenetic event types based on residence times
  group_by(time, map) %>%
  mutate(percs = 100 * counts / sum(counts)) %>% # Compute % across source for a given map and time
  arrange(desc(time)) %>% 
  group_by(source, map) %>% 
  mutate(cum_counts = cumsum(counts)) %>% # Compute cumulative counts along time
  ungroup() %>% 
  select(time, source, map, counts, percs, cum_counts, LTT_rates, rates, area_richness, residence_time)

# Replace NA (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.na(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.nan(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Aggregate mean counts/percs/rates across maps
all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(source, time) %>%  # Aggregate across maps
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_LTT_rates = mean(LTT_rates),
            mean_rates = mean(rates),
            area_richness = mean(area_richness),
            residence_time = mean(residence_time)) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup() %>% 
  select(time, source, mean_counts, mean_percs, mean_cum_counts, mean_LTT_rates, mean_rates, area_richness, residence_time)

# Save ggplot df of (mean) counts per source bioregions across dispersal events
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

### 2.2/ Sum across all source bioregions ####

names(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)
names(all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot)

# For all maps
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot <- all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, map) %>%  # Summarize across source bioregions
  summarize(counts = sum(counts),
            cum_counts = sum(cum_counts),
            total_richness = sum(area_richness),
            total_time = sum(residence_time)) %>% 
  mutate(LTT_rates = counts / total_richness / window_width) %>% # Compute updated rates for all cladogenetic event types based on LTT data
  mutate(rates = counts / total_time) %>% # Compute updated rates for all cladogenetic event types based on residence times
  ungroup() %>% 
  select(time, map, counts, cum_counts, LTT_rates, rates, total_richness, total_time)

# Replace NA (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[is.na(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[is.na(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[is.nan(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates)] <- 0
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[is.nan(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[(all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Aggregated across maps
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot <- all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%  # Summarize across source bioregions
  summarize(mean_counts = sum(mean_counts),
            mean_cum_counts = sum(mean_cum_counts),
            total_richness = sum(area_richness),
            total_time = sum(residence_time)) %>% 
  mutate(mean_LTT_rates = mean_counts / total_richness / window_width) %>% # Compute updated rates for all cladogenetic event types on LTT data
  mutate(mean_rates = mean_counts / total_time) %>% # Compute updated rates for all cladogenetic event types on residence times
  ungroup() %>% 
  select(time, mean_counts, mean_cum_counts, mean_LTT_rates, mean_rates, total_richness, total_time)

# Replace NA (due to dividing by zero) with 0.000
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[is.na(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates)] <- 0
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[is.na(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[is.nan(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates)] <- 0
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[is.nan(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates >= Inf)] <- 0
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[(all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates >= Inf)] <- 0

# Save ggplot df of (mean) counts across dispersal events all source bioregions aggregated (overall)
# saveRDS(object = all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))


##### 3/ Plot stacked bars of mean raw number of cladogenetic events per source bioregions along sliding windows ####

# Load data.frame of mean counts per source bioregions across dispersal events
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/Cladogenetic_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/Cladogenetic_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/Cladogenetic_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/Cladogenetic_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_mean_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse", 
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of cladogenetic events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of cladogenetic events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_mean_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 4/ Plot stacked bars of percentages of cladogenetic events per source bioregions along sliding windows ####

# Load data.frame of mean counts per source bioregions across dispersal events
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/Cladogenetic_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/Cladogenetic_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/Cladogenetic_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/Cladogenetic_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_perc_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                   mapping = aes(y = mean_percs, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse", 
                     # limits = c(150, 0)) + # Set limits for Oldest
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of cladogenetic events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of cladogenetic events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_perc_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 5/ Plot rates of cladogenetic events per source bioregions along sliding windows ####

### 5.1/ Include overall rates in the data.frame ####

# Load data.frames of counts per maps across dispersal events
# all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Load data.frames of mean counts across dispersal events
# all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

## 5.1.1/ Counts per maps

all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot_to_bind <- all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot %>% 
  rename(area_richness = total_richness) %>% 
  rename(residence_time = total_time) %>% 
  mutate(source = "total",
         percs = 100) %>% 
  select(time, source, map, counts, percs, cum_counts, LTT_rates, rates, area_richness, residence_time)

all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total <- rbind(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, all_cladogenetic_events_count_per_maps_overall_all_sliding_windows_ggplot_to_bind)

# Save data.frame with all rates, including total
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))

## 5.1.2/ Mean counts aggregated across maps

all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot_to_bind <- all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot %>% 
  rename(area_richness = total_richness) %>% 
  rename(residence_time = total_time) %>% 
  mutate(source = "total",
         mean_percs = 100) %>% 
  select(time, source, mean_counts, mean_percs, mean_cum_counts, mean_LTT_rates, mean_rates, area_richness, residence_time)

all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total <- rbind(all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, all_cladogenetic_events_mean_counts_overall_all_sliding_windows_ggplot_to_bind)

# Save data.frame with all rates, including total
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))

### 5.2/ Plot rates of cladogenetic event types (= diversification rates) per source bioregions along sliding windows ####

# Extract max (mean) rates
max_rates <- max(all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total$rates)
max_mean_rates <- max(all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total$mean_rates, na.rm = T)

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/Cladogenetic_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/Cladogenetic_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/Cladogenetic_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/Cladogenetic_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
lines_rates_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total) +
  
  # Plot mean lines + 1000 replicates
  geom_smooth(data = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total,
              mapping = aes(y = rates, x = time, group = source, col = source, fill = source),
              method = "gam", se = TRUE, method.args = list(gamma = 3),
              linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = all_cladogenetic_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total,
  #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # # Plot mean lines only
  # geom_line(data = all_cladogenetic_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total,
  #         mapping = aes(y = mean_rates, x = time, group = source, col = source),
  #         alpha = 1.0,
  #         linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("Rates of cladogenetic events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  ylim(c(0, 0.06)) +
  # ylim(c(0, max_mean_rates)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names_with_total, values = unname(colors_list_for_areas_with_total)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(lines_rates_per_cladogenetic_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 6/ Compute in situ speciation counts and rates ####

### 6.1/ Compute counts/rates per source bioregions ####

## 6.1.1/ Extract counts of in situ speciation events

## Load the array of counts of events per source/dest areas
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))
# DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/Event_types/DEC_J_BSM_counts_all_unique_events_all_sliding_windows_",window_width,"My_source_dest_array.rds"))

# Extract counts only for targeted type of events
in_situ_speciation_event_types_list <- c("range inheritance (y)")

# Extract counts only for targeted type of events
in_situ_speciation_events_count_matrices_per_maps_all_sliding_windows <- DEC_J_BSM_all_unique_events_all_sliding_windows_source_dest_array[ , ,in_situ_speciation_event_types_list, , ]

dim(in_situ_speciation_events_count_matrices_per_maps_all_sliding_windows)
# Dim 1 = Rows = Sources
# Dim 2 = Cols = Dest
# Dim 3 = Maps
# Dim 4 = Sliding windows

## Format data for ggplot
in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_matrices_per_maps_all_sliding_windows %>%
  reshape2::melt()
names(in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot) <- c("source", "dest", "map", "window", "counts")

# Match time with time-window
window_time_df <- data.frame(window = levels(in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot$window),
                             time = mean_time_list)

in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot <- left_join(x = in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, y = window_time_df)

# Save ggplot df dataframe of in situ speciation events per areas transitions x time-windows
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

## 6.1.2/ Aggregate per source bioregions and merge with area richness and residence times

## Aggregate per source bioregions
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, source, map) %>%
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
  dplyr::select(areas, time, mean_counts) %>% 
  rename(area_richness = mean_counts)
# Merge with LTT to get richness
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- left_join(x = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, y = total_richness_per_bioregions, join_by(time == time, source == areas))
# Divide counts by richness and time span to get rates in events per lineage per My
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  mutate(LTT_rates = counts / area_richness / window_width)

## Compute rates based on residence times
# Load residence times data
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_rough_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_MCC_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Youngest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))
# residence_times_per_sliding_windows_in_areas_melted_df <- readRDS(file = paste0("./outputs/Continuous_events_counts/Ponerinae_Oldest_phylogeny_1534t/residence_times_per_sliding_windows_in_areas_melted_df_", window_width,"My.rds"))

# Merge with residence times
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- left_join(x = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                       y = residence_times_per_sliding_windows_in_areas_melted_df[, c("area", "smoothed_time", "map", "residence_time")],
                                                                                                       join_by(source == area, time == smoothed_time, map == map))
# Divide counts by richness and time span to get rates in events per lineage per My
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  mutate(rates = counts / residence_time)

names(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)

## 6.1.3/ Compute percentages, cumulative counts, and rates

# Compute percentages, cumulative counts, and rates
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, map) %>%
  mutate(percs = 100 * counts / sum(counts)) %>% # Compute % across source for a given map and time
  arrange(desc(time)) %>% 
  group_by(source, map) %>% 
  mutate(cum_counts = cumsum(counts)) %>% # Compute cumulative counts along time
  ungroup() %>% 
  dplyr::select(time, source, map, counts, percs, cum_counts, LTT_rates, rates, area_richness, residence_time)

# Replace NA (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.na(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.na(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[is.nan(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates)] <- 0
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[is.nan(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates[(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates[(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$rates >= Inf)] <- 0

## 6.1.4/ Aggregate mean counts/percs/rates across maps

in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(source, time) %>%  # Aggregate across maps
  summarize(mean_counts = mean(counts),
            mean_cum_counts = mean(cum_counts),
            mean_LTT_rates = mean(LTT_rates),
            mean_rates = mean(rates),
            area_richness = mean(area_richness),
            residence_time = mean(residence_time)) %>% 
  group_by(time) %>% # Compute % summing to 100% for a given time
  mutate(mean_percs = 100 * mean_counts / sum(mean_counts)) %>% 
  ungroup() %>% 
  dplyr::select(time, source, mean_counts, mean_percs, mean_cum_counts, mean_LTT_rates, mean_rates, area_richness, residence_time)

# Adjust order of Bioregions
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)
in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source <- factor(in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot$source, levels = areas_list, labels = bioregion_names)

# Save ggplot df of (mean) counts per source bioregions across dispersal events
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

### 6.2/ Sum across all source bioregions ####

names(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot)
names(in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot)

# For all maps
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time, map) %>%  # Summarize across source bioregions
  summarize(counts = sum(counts),
            cum_counts = sum(cum_counts),
            total_richness = sum(area_richness),
            total_time = sum(residence_time)) %>% 
  mutate(LTT_rates = counts / total_richness / window_width) %>% # Compute updated rates for all in situ speciation event types based on LTT data
  mutate(rates = counts / total_time) %>% # Compute updated rates for all in situ speciation event types based on residence times
  ungroup() %>% 
  dplyr::select(time, map, counts, cum_counts, LTT_rates, rates, total_richness, total_time)

# Replace NA (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[is.na(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates)] <- 0
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[is.na(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[is.nan(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates)] <- 0
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[is.nan(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates[(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$LTT_rates >= Inf)] <- 0
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates[(in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot$rates >= Inf)] <- 0

# Aggregated across maps
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot <- in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot %>% 
  group_by(time) %>%  # Summarize across source bioregions
  summarize(mean_counts = sum(mean_counts),
            mean_cum_counts = sum(mean_cum_counts),
            total_richness = sum(area_richness),
            total_time = sum(residence_time)) %>% 
  mutate(mean_LTT_rates = mean_counts / total_richness / window_width) %>% # Compute updated rates for all in situ speciation event types on LTT data
  mutate(mean_rates = mean_counts / total_time) %>% # Compute updated rates for all in situ speciation event types on residence times
  ungroup() %>% 
  dplyr::select(time, mean_counts, mean_cum_counts, mean_LTT_rates, mean_rates, total_richness, total_time)

# Replace NA (due to dividing by zero) with 0.000
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[is.na(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates)] <- 0
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[is.na(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates)] <- 0
# Replace NaN (due to dividing by zero) with 0.000
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[is.nan(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates)] <- 0
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[is.nan(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates)] <- 0
# Replace Inf (due to dividing by zero) with 0.000
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates[(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_LTT_rates >= Inf)] <- 0
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates[(in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot$mean_rates >= Inf)] <- 0

# Save ggplot df of (mean) counts across dispersal events all source bioregions aggregated (overall)
# saveRDS(object = in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
saveRDS(object = in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))


##### 7/ Plot stacked bars of mean raw number of in situ speciation events per source bioregions along sliding windows ####

# Load data.frame of mean counts per source bioregions across dispersal events
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/In_situ_speciation_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/In_situ_speciation_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/In_situ_speciation_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/In_situ_speciation_events_count_mean_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_mean_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                      mapping = aes(y = mean_counts, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0)) + # Set limits
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Counts of in situ speciation events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Counts of in situ speciation events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_mean_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_mean_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_mean_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 8/ Plot stacked bars of percentages of in situ speciation events per source bioregions along sliding windows ####

# Load data.frame of mean counts per source bioregions across dispersal events
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/In_situ_speciation_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/In_situ_speciation_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/In_situ_speciation_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/In_situ_speciation_events_count_perc_barplots_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
barplot_count_perc_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot,
                                                                                                      mapping = aes(y = mean_percs, x = time, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, position = "stack", width = 1.0,
           col = "black", linewidth = 0.3) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse", 
                     # limits = c(150, 0)) + # Set limits
                     limits = c(100, 0)) + # Set limits
  
  # Set plot title +
  ggtitle(label = paste0("Percentages of in situ speciation events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of in situ speciation events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
barplot_count_perc_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- barplot_count_perc_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(barplot_count_perc_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 9/ Plot rates of in situ speciation events per source bioregions along sliding windows ####

### 9.1/ Include overall rates in the data.frame ####

# Load data.frames of counts per maps across dispersal events
# in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

# Load data.frames of mean counts across dispersal events
# in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_overall_all_sliding_windows_",window_width,"My_ggplot.rds"))
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot <- readRDS(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot.rds"))

## 9.1.1/ Counts per maps

in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot_to_bind <- in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot %>% 
  rename(area_richness = total_richness) %>% 
  rename(residence_time = total_time) %>% 
  mutate(source = "total",
         percs = 100) %>% 
  dplyr::select(time, source, map, counts, percs, cum_counts, LTT_rates, rates, area_richness, residence_time)

in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total <- rbind(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot, in_situ_speciation_events_count_per_maps_overall_all_sliding_windows_ggplot_to_bind)

# Save data.frame with all rates, including total
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))

## 9.1.2/ Mean counts aggregated across maps

in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot_to_bind <- in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot %>% 
  rename(area_richness = total_richness) %>% 
  rename(residence_time = total_time) %>% 
  mutate(source = "total",
         mean_percs = 100) %>% 
  dplyr::select(time, source, mean_counts, mean_percs, mean_cum_counts, mean_LTT_rates, mean_rates, area_richness, residence_time)

in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total <- rbind(in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot, in_situ_speciation_events_mean_counts_overall_all_sliding_windows_ggplot_to_bind)

# Save data.frame with all rates, including total
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))
# saveRDS(object = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total, file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_",window_width,"My_ggplot_with_total.rds"))

### 9.2/ Plot rates of in situ speciation event types (= diversification rates) per source bioregions along sliding windows ####

# Extract max (mean) rates
max_rates <- max(in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total$rates)
max_mean_rates <- max(in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total$mean_rates, na.rm = T)

# # Remove rates of IndoMalaya before 85My as they are based on a single lineage
# in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered <- in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total %>% 
#   filter(!(source == "Indomalaya" & time > 80))
# in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered <- in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total %>% 
#   filter(!(source == "Indomalaya" & time > 80))

# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/In_situ_speciation_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/In_situ_speciation_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/In_situ_speciation_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
# pdf(file = paste0("./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/In_situ_speciation_events_rates_all_source_bioregions_all_sliding_windows_",window_width,"My.pdf"),
    width = 10, height = 6)

# Generate plot
lines_rates_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- ggplot(data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total
                                                                                                     # data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered
                                                                                                     ) +
  
  # Plot mean lines + 1000 replicates
  geom_smooth(data = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total,
              # data = in_situ_speciation_events_count_per_maps_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered,
              mapping = aes(y = rates, x = time, group = source, col = source, fill = source),
                            method = "gam", se = TRUE, method.args = list(gamma = 3),
                            linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(# data = in_situ_speciation_events_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total,
  #           data = in_situ_speciation_events_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered,
  #           mapping = aes(y = rates, x = time, group = interaction(source, map), col = source),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # # Plot mean lines only
  # geom_line(# data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total,
  #           data = in_situ_speciation_events_mean_counts_all_source_bioregions_all_sliding_windows_ggplot_with_total_filtered,
  #           mapping = aes(y = mean_rates, x = time, group = source, col = source),
  #                         alpha = 1.0,
  #                         linewidth = 2.0) +

  # Set plot title +
  ggtitle(label = paste0("Rates of in situ speciation events\nper Source bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.08)) +
  # ylim(c(0, max_mean_rates)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(150, 0) # Set limits for Oldest
                     limits = c(100, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names_with_total, values = unname(colors_list_for_areas_with_total)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
  # color = rev(time_strata_col),
  # angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_rates_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot <- lines_rates_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot %>% 
  add_aesthetics_barplot()

# Plot
print(lines_rates_per_in_situ_speciation_events_all_source_bioregions_all_sliding_windows_ggplot)

dev.off()


##### 10/ Overall counts of in situ speciation events per bioregions #####

### 10.1/ Compute counts of in situ speciation events per bioregions ####

## Load the array of counts of events per source/dest areas
# DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")
DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")
# DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_Youngest_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")
# DEC_J_BSM_all_unique_events_source_dest_array <- readRDS(file = "./outputs/BSM/Ponerinae_Oldest_phylogeny_1534t/DEC_J_BSM_counts_all_unique_events_source_dest_array.rds")

# Extract counts only for targeted type of events
in_situ_speciation_event_types_list <- c("range inheritance (y)")

# Extract counts only for targeted type of events
in_situ_speciation_events_count_matrices_per_maps_overall <- DEC_J_BSM_all_unique_events_source_dest_array[ , ,in_situ_speciation_event_types_list, , ]

dim(in_situ_speciation_events_count_matrices_per_maps_overall)
# Dim 1 = Rows = Sources
# Dim 2 = Cols = Dest
# Dim 3 = Maps
# Dim 4 = Time-strata

## Format data for ggplot
in_situ_speciation_events_count_per_maps_overall_ggplot <- in_situ_speciation_events_count_matrices_per_maps_overall %>%
  reshape2::melt()
names(in_situ_speciation_events_count_per_maps_overall_ggplot) <- c("source", "dest", "map", "stratum", "counts")

## Mean across maps
in_situ_speciation_events_mean_count_overall_ggplot <- in_situ_speciation_events_count_per_maps_overall_ggplot %>% 
  group_by(source, dest, stratum) %>% 
  summarize(mean_counts = mean(counts)) %>%
  ungroup()

## Sum across time-strata, and destination bioregions
in_situ_speciation_events_mean_count_per_source_bioregions_overall <- in_situ_speciation_events_mean_count_overall_ggplot %>% 
  group_by(source) %>%
  summarize(mean_counts = sum(mean_counts)) %>%
  ungroup()

## Version with sd and HPD95%

# Aggregate counts across maps for destinations: mean, perc, sd, HPD95%
in_situ_speciation_events_summary_counts_overall_ggplot <- in_situ_speciation_events_count_per_maps_overall_ggplot %>% 
  group_by(source, map) %>% 
  summarize(counts = sum(counts)) %>% # Sum across strata x dests
  group_by(source) %>% # Aggregate counts across maps using summary statistics
  summarize(mean_counts = mean(counts),
            sd_counts = sd(counts),
            HPD2.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[1],
            HPD97.5_counts = BayesTwin::HPD(sample = counts, cred_int = 0.95)[2]) %>% 
  ungroup() %>% # Compute overall percentages
  mutate(mean_perc = 100 * mean_counts / sum(mean_counts))


### 10.2/ Plot mean counts of in situ speciation events per bioregions ####

# Adjust color scheme for bioregions
areas_list <- c("A", "U", "I", "R", "N", "E", "W")
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names
colors_list_for_areas

# Adjust order of bioregions
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "E-PA", "W-PA")
in_situ_speciation_events_mean_count_per_source_bioregions_overall$source <- factor(in_situ_speciation_events_mean_count_per_source_bioregions_overall$source, levels = levels(in_situ_speciation_events_mean_count_per_source_bioregions_overall$source), labels = bioregion_names_reduced)
in_situ_speciation_events_summary_counts_overall_ggplot$source <- factor(in_situ_speciation_events_summary_counts_overall_ggplot$source, levels = levels(in_situ_speciation_events_summary_counts_overall_ggplot$source), labels = bioregion_names_reduced)


## GGplot without sd
# pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_rough_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)
pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)
# pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)
# pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot.pdf", width = 10, height = 6)

barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot <- ggplot(data = in_situ_speciation_events_mean_count_per_source_bioregions_overall,
                                                                           mapping = aes(y = mean_counts, x = source, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, col = "black") +
  
  # Set plot title +
  ggtitle(label = paste0("In situ speciation events per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot <- barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot)

dev.off()


## GGplot with sd
pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_MCC_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot_with_sd.pdf", width = 10, height = 6)
# pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_Youngest_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot_with_sd.pdf", width = 10, height = 6)
# pdf(file = "./outputs/Empirical_diversification_rates/Ponerinae_Oldest_phylogeny_1534t/In_situ_speciation_events_count_per_source_bioregions_overall_barplot_with_sd.pdf", width = 10, height = 6)

barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot <- ggplot(data = in_situ_speciation_events_summary_counts_overall_ggplot,
                                                                                    mapping = aes(y = mean_counts, x = source, fill = source)) +
  # Plot barplot
  geom_col(alpha = 1.0, col = "black") +
  
  # Add SD as error bars
  geom_errorbar(aes(ymin = mean_counts, ymax = mean_counts + sd_counts, x = source),
                width = 0.15,
                position = "identity") +
  
  # Set plot title +
  ggtitle(label = paste0("In situ speciation events per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Counts of events") +
  
  # Adjust color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot <- barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot %>% 
  add_aesthetics_density_curve()

# Plot
print(barplot_in_situ_speciation_events_mean_count_per_source_bioregions_ggplot)

dev.off()



# Make an animated GIF of evolution of diversification rates in time!
# Plot alongside evolution of LTT and dispersal source/dest rates



## All can be segregated per clades!
# Across clades to compare clades dynamics
# Within clades, to zoom on a particular subclade of interest
# Try merging probability of each bioregions membership on the ContMap by using weighted color in RGB space (do that overall and subset specific subclades for plotting)

