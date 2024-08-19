##### Script 07: Plot residence times and LTT #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compute residence time of each state across time-strata
# Compute lineages per bioregions across continous time (LTT)
# Summarize all information in array across time x bioregions x BSM maps

# Plot the evolution of residence times and lineages counts (LTT) across time 
   # Compute palodiversity across bioregions in discrete (time strata) and continuous time (LTT)
   # Raw residence times/diversity vs. relative diversity (%)

###

### Inputs

# Array of residence times per 
# Best fitted Biogeographic model

###

### Outputs

# Array of residence times x bioregions x time-strata x BSM maps 
# Summary stats of residence times x bioregions x time-strata across BSM maps 
# Stacked bar graph of residence times x bioregions x time-strata (Raw and %)

# Simmaps conversion of BSM maps
# Array of lineage counts x bioregions x time x BSM maps 
# Summary stats of lineage counts x bioregions x time across BSM maps 
# Cumulative plots of lineage counts x bioregions through time (Raw and %)

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

# Newer version of ggplot2 are not compatible with gganimate
# versions::install.versions(pkgs = "ggplot2", versions = "3.5.0") 
# old_version <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.5.0.tar.gz"
# utils::install.packages(old_version, repos = NULL, type = "source")
# remotes::install_version(package = "ggplot2", version = "3.5.0", repos = "http://cran.us.r-project.org")
# old_version <- "https://cran.r-project.org/package=ggplot2&version=3.5.0"
# utils::install.packages(pkgs = old_version, repos = NULL)
# utils::install.packages("./packages/ggplot2_3.5.0.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(phytools)
library(geiger)
library(ape)
library(BioGeoBEARS)
library(qpdf)
library(MultinomialCI)    # For 95% CIs on BSM counts
library(magick)           # For animated GIF
library(pdftools)         # To read PDF
library(BayesTwin)        # To compute HPD intervals
library(gganimate)

# devtools::install_github("nmatzke/BioGeoBEARS")


### 1.2/ Load BSM outputs ####

# Load time-stratified DEC+J model output
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")

# Load BSM outputs
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/DEC_J_BSM_output.rds")

# Extract records of events
DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
DEC_J_ana_events_tables <- DEC_J_BSM_output$RES_ana_events_tables


### 1.3/ Load custom functions ####

source("./functions/BSM_to_phytools_SM_custom.R") # Function to convert BSM output to phytools.simmap format
source("./functions/compute_residence_times_by_time_stratum.R") # Functions to compute residence time per state/range from BSMs in phytools.simmap format


##### 2/ Summarize residence times from simmaps #####

### 2.1/ Extract residence times from all simmaps ####

# Load simmaps of BS maps
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/DEC_J_simmaps.rds")

## Initiate output array

# Extract number of time-strata
nb_strata <- length(DEC_J_fit$inputs$timeperiods)

# Extract full range list
returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
# ranges_list <- returned_mats$ranges_list
reduced_ranges_list = returned_mats$ranges_list
max_range_size <- max(nchar(reduced_ranges_list))
ranges_list <- generate_list_ranges(areas_list = returned_mats$areanames, max_range_size = max_range_size, include_null_range = DEC_J_fit$inputs$include_null_range)

# Initiate final array for residence times
residence_times_all_BSM_maps_array <- array(data = NA,
                                            dim = c(2, length(ranges_list) + 1, nb_strata, length(DEC_J_clado_events_tables)),
                                            dimnames = list(c("raw", "perc"), c(ranges_list, "total"), paste0("Stratum_", 1:nb_strata), paste0("Map_", 1:length(DEC_J_clado_events_tables))))
dim(residence_times_all_BSM_maps_array) # Dim 1: stats type (raw/perc), Dim 2: Range/state, Dim 3: Time-strata, Dim 4: BSM maps


## Loop across all BSM maps 
for (i in 1:length(DEC_J_simmaps))
{
  # i <- 1
  
  # Extract simmap
  
  DEC_J_simmap_i <- DEC_J_simmaps[[i]]$simmap
  
  # Extract Residence times x states x Time-stratum
  residence_times_all_strata_i <- compute_residence_times_for_all_time_strata(simmap = DEC_J_simmap_i,
                                                                              model_fit = DEC_J_fit)
  residence_times_all_strata_i
  
  # Add to final array
  residence_times_all_BSM_maps_array[ , , ,i] <- residence_times_all_strata_i
  
  dim(residence_times_all_strata_i)
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Residence times computed for Time-stratifed Stochastic map n°", i, "/", length(DEC_J_clado_events_tables),"\n"))
  }
}

# Save residence times states x maps x Time-strata
saveRDS(object = residence_times_all_BSM_maps_array, file = "./outputs/BSM/residence_times_all_BSM_maps_array.rds")
residence_times_all_BSM_maps_array <- readRDS(file = "./outputs/BSM/residence_times_all_BSM_maps_array.rds")


### 2.2/ Compute summary statistics of residence time across all maps ####

# Load array of residence times states x maps x Time-strata
residence_times_all_BSM_maps_array <- readRDS(file = "./outputs/BSM/residence_times_all_BSM_maps_array.rds")

# (mean, sd, 95% CI, 95% HPD)
dim(residence_times_all_BSM_maps_array)
residence_times_all_BSM_maps_array

# Record summary statistics only for observed states across all BSMs
all_state_times <- apply(X = residence_times_all_BSM_maps_array["raw", , , ], MARGIN = 1, FUN = sum)
observed_states <- ranges_list[all_state_times[-length(all_state_times)] != 0]

# Dim 1: stats type (raw/perc), Dim 2: Range/state, Dim 3: Time-strata, Dim 4: BSM maps

residence_times_all_BSM_maps_summary_array <- array(data = NA,
                                                    dim = c(2, length(observed_states) + 1, nb_strata, 7),
                                                    dimnames = list(c("raw", "perc"), c(observed_states, "total"), paste0("Stratum_", 1:nb_strata), c("Mean", "Median", "sd", "2.5% CI", "97.5% CI", "2.5% HPD", "97.5% HPD")))
dim(residence_times_all_BSM_maps_summary_array)

residence_times_all_BSM_maps_summary_array[ , , , "Mean"] <- apply(X = residence_times_all_BSM_maps_array[,c(observed_states, "total"),,], MARGIN = c(1,2,3), FUN = mean)
residence_times_all_BSM_maps_summary_array[ , , , "Median"] <- apply(X = residence_times_all_BSM_maps_array[,c(observed_states, "total"),,], MARGIN = c(1,2,3), FUN = median)
residence_times_all_BSM_maps_summary_array[ , , , "sd"] <- apply(X = residence_times_all_BSM_maps_array[,c(observed_states, "total"),,], MARGIN = c(1,2,3), FUN = sd)
CI95 <- apply(X = residence_times_all_BSM_maps_array[,c(observed_states, "total"),,], MARGIN = c(1,2,3), FUN = quantile, probs = c(0.025, 0.975))
residence_times_all_BSM_maps_summary_array[ , , , c("2.5% CI", "97.5% CI")] <- aperm(CI95, c(2,3,4,1))
HPD95 <- apply(X = residence_times_all_BSM_maps_array[,c(observed_states, "total"),,], MARGIN = c(1,2,3), FUN = BayesTwin::HPD, cred_int = 0.95)
residence_times_all_BSM_maps_summary_array[ , , , c("2.5% HPD", "97.5% HPD")] <- aperm(HPD95, c(2,3,4,1))

aperm(residence_times_all_BSM_maps_summary_array["raw", , , ], c(1,3,2))

# Save summary statistics of residence times per states, per time strata across all maps
saveRDS(object = residence_times_all_BSM_maps_summary_array, file = "./outputs/BSM/residence_times_all_BSM_maps_summary_array.rds")


##### 3/ Aggregate residence times per unique areas #####

# Aggregate residence times per unique areas instead of ranges
# Time for ranges encompassing multiple areas is equally split across those areas

# Load array of residence times states x maps x Time-strata
residence_times_all_BSM_maps_array <- readRDS(file = "./outputs/BSM/residence_times_all_BSM_maps_array.rds")

# Extract ranges
all_ranges <- dimnames(residence_times_all_BSM_maps_array)[[2]]
all_ranges <- all_ranges[-length(all_ranges)] # Remove "total"

# Extract unique areas
all_areas <- all_ranges[nchar(all_ranges) == 1]
all_areas <- all_areas[!(all_areas == "_")] # Remove empty area

# Split ranges in unique areas
ranges_split_list <- str_split(string = all_ranges, pattern = "")

### 3.1/ Divide residence times of ranges equally across areas ####

residence_times_split_all_BSM_maps_array <- residence_times_all_BSM_maps_array
for (i in seq_along(all_ranges))
{
  # i <- 10
  
  # Extract focal range
  range_i <- all_ranges[i]
  
  # Extract number of areas composing the range
  nb_areas_i <- nchar(range_i)
  
  # Divide residence times accordingly
  residence_times_split_all_BSM_maps_array[,range_i,,] <- residence_times_split_all_BSM_maps_array[,range_i,,] / nb_areas_i
  
}

# residence_times_split_all_BSM_maps_array[,,,"Map_1"]

### 3.2/ Sum residence times per unique areas ####

## Initiate array

# Extract dimensions
nb_strata <- dim(residence_times_all_BSM_maps_array)[3]
nb_maps <- dim(residence_times_all_BSM_maps_array)[4]

# Initiate array for aggregated residence times per unique areas 
residence_times_per_areas_all_BSM_maps_array <- array(data = NA,
                                                      dim = c(2, length(all_areas) + 1, nb_strata, nb_maps),
                                                      dimnames = list(c("raw", "perc"), c(all_areas, "total"), paste0("Stratum_", 1:nb_strata), paste0("Map_", 1:nb_maps)))
dim(residence_times_per_areas_all_BSM_maps_array) # Dim 1: stats type (raw/perc), Dim 2: Areas, Dim 3: Time-strata, Dim 4: BSM maps

## Loop per areas to aggregate data
for (i in seq_along(all_areas))
{
  # i <- 1
  
  # Extract focal area
  areas_i <- all_areas[i]
  
  # Find matching ranges
  matching_ranges_i <- all_ranges[str_detect(string = ranges_split_list, pattern = areas_i)]
  
  # Sum counts across all matching ranges
  residence_times_per_areas_all_BSM_maps_array[,areas_i,,] <- apply(X = residence_times_split_all_BSM_maps_array[,matching_ranges_i,,], MARGIN = c(1,3,4), FUN = sum)
}

# Copy total data
residence_times_per_areas_all_BSM_maps_array[,"total",,] <- residence_times_all_BSM_maps_array[,"total",,]

# residence_times_per_areas_all_BSM_maps_array[,,,"Map_1"]

# Save summary statistics of residence times per unique areas, per time strata across all maps
saveRDS(object = residence_times_per_areas_all_BSM_maps_array, file = "./outputs/BSM/residence_times_per_areas_all_BSM_maps_array.rds")


### 3.3/ Compute summary statistics of residence time across all maps for unique areas ####

# Load array of residence times states x maps x Time-strata
residence_times_per_areas_all_BSM_maps_array <- readRDS(file = "./outputs/BSM/residence_times_per_areas_all_BSM_maps_array.rds")

# (mean, sd, 95% CI, 95% HPD)
residence_times_per_areas_all_BSM_maps_array
# Dim 1: stats type (raw/perc), Dim 2: Areas, Dim 3: Time-strata, Dim 4: BSM maps

# Extract dimensions
all_areas <- dimnames(residence_times_per_areas_all_BSM_maps_array)[[2]]
all_areas <- all_areas[-length(all_areas)] # Remove "total"
nb_strata <- dim(residence_times_per_areas_all_BSM_maps_array)[3]
nb_maps <- dim(residence_times_per_areas_all_BSM_maps_array)[4]

residence_times_per_areas_all_BSM_maps_summary_array <- array(data = NA,
                                                              dim = c(2, length(all_areas) + 1, nb_strata, 7),
                                                              dimnames = list(c("raw", "perc"), c(all_areas, "total"), paste0("Stratum_", 1:nb_strata), c("Mean", "Median", "sd", "2.5% CI", "97.5% CI", "2.5% HPD", "97.5% HPD")))
dim(residence_times_per_areas_all_BSM_maps_summary_array)

residence_times_per_areas_all_BSM_maps_summary_array[ , , , "Mean"] <- apply(X = residence_times_per_areas_all_BSM_maps_array[,c(all_areas, "total"),,], MARGIN = c(1,2,3), FUN = mean)
residence_times_per_areas_all_BSM_maps_summary_array[ , , , "Median"] <- apply(X = residence_times_per_areas_all_BSM_maps_array[,c(all_areas, "total"),,], MARGIN = c(1,2,3), FUN = median)
residence_times_per_areas_all_BSM_maps_summary_array[ , , , "sd"] <- apply(X = residence_times_per_areas_all_BSM_maps_array[,c(all_areas, "total"),,], MARGIN = c(1,2,3), FUN = sd)
CI95 <- apply(X = residence_times_per_areas_all_BSM_maps_array[,c(all_areas, "total"),,], MARGIN = c(1,2,3), FUN = quantile, probs = c(0.025, 0.975))
residence_times_per_areas_all_BSM_maps_summary_array[ , , , c("2.5% CI", "97.5% CI")] <- aperm(CI95, c(2,3,4,1))
HPD95 <- apply(X = residence_times_per_areas_all_BSM_maps_array[,c(all_areas, "total"),,], MARGIN = c(1,2,3), FUN = BayesTwin::HPD, cred_int = 0.95)
residence_times_per_areas_all_BSM_maps_summary_array[ , , , c("2.5% HPD", "97.5% HPD")] <- aperm(HPD95, c(2,3,4,1))

aperm(residence_times_per_areas_all_BSM_maps_summary_array["raw", , , ], c(1,3,2))
aperm(residence_times_per_areas_all_BSM_maps_summary_array["perc", , , ], c(1,3,2))

# Save summary statistics of residence times per states, per time strata across all maps
saveRDS(object = residence_times_per_areas_all_BSM_maps_summary_array, file = "./outputs/BSM/residence_times_per_areas_all_BSM_maps_summary_array.rds")


##### 4/ Plot residence times in areas per time-strata ####

# Load summary statistics of residence times per states, per time strata across all maps
residence_times_per_areas_all_BSM_maps_summary_array <- readRDS(file = "./outputs/BSM/residence_times_per_areas_all_BSM_maps_summary_array.rds")

### 4.1/ Raw residence times ####

## Extract data in ggplot format

aperm(residence_times_per_areas_all_BSM_maps_summary_array["raw", , , "Mean"], c(1,2))
residence_times_per_areas_raw_ggplot <- aperm(residence_times_per_areas_all_BSM_maps_summary_array["raw", , , "Mean"], c(1,2))

# Rename Time-strata with geological era names
colnames(residence_times_per_areas_raw_ggplot)  <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")

# Rename areas with bioregion names
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic", "Total")
row.names(residence_times_per_areas_raw_ggplot) <- bioregion_names

residence_times_per_areas_raw_ggplot <- reshape2::melt(residence_times_per_areas_raw_ggplot)
names(residence_times_per_areas_raw_ggplot) <- c("Bioregions", "Time_strata", "Lineage_counts")
residence_times_per_areas_raw_ggplot

# Save ggplot df for Raw residence times per bioregions, per geological era
saveRDS(residence_times_per_areas_raw_ggplot, file = "./outputs/LTT/residence_times_per_areas_raw_ggplot.rds")

## Set colors and legend

# Extract areas
all_areas <- dimnames(residence_times_per_areas_all_BSM_maps_summary_array)[[2]]
# Remove total
all_areas <- all_areas[all_areas != "total"]

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
names(colors_list_for_areas) <- bioregion_names[-length(bioregion_names)]

# Set levels of areas in defined order
residence_times_per_areas_raw_ggplot$Bioregions <- factor(residence_times_per_areas_raw_ggplot$Bioregions, levels = bioregion_names)

# Set levels of geological eras/time-strata in defined order (reverse order)
residence_times_per_areas_raw_ggplot$Time_strata <- factor(residence_times_per_areas_raw_ggplot$Time_strata, levels = rev(levels(residence_times_per_areas_raw_ggplot$Time_strata)))

# Remove the "Total" values
remove_total_indices <- residence_times_per_areas_raw_ggplot$Bioregions != "Total"

View(residence_times_per_areas_raw_ggplot)

## GGplot

pdf(file = "./outputs/LTT/Residence_times_per_areas_raw_ggplot.pdf", width = 12, height = 8)
ggplot_residence_times_per_areas_raw <- ggplot(data = residence_times_per_areas_raw_ggplot[remove_total_indices, ], mapping = aes(x = Time_strata, y = Lineage_counts, fill = Bioregions)) +
  
  ## Plot bars
  geom_bar(position = "stack", stat = "identity", width = 0.5) +
  
  # # Reverse time axis
  # scale_x_reverse() +
  
  # Set fill colors
  scale_fill_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "Residence times per bioregion") +
  xlab(label = "Geological periods") +
  ylab(label = "Residence time  [My]") +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, size = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.0),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 8)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_residence_times_per_areas_raw)

dev.off()


### 4.2/ Proportional residence times (in %) ####

## Extract data in ggplot format

aperm(residence_times_per_areas_all_BSM_maps_summary_array["perc", , , "Mean"], c(1,2))
residence_times_per_areas_perc_ggplot <- aperm(residence_times_per_areas_all_BSM_maps_summary_array["perc", , , "Mean"], c(1,2))

# Rename Time-strata with geological era names
colnames(residence_times_per_areas_perc_ggplot)  <- c("PPHo", "Miocene", "Oligocene", "Eocene", "Paleocene", "Late Cr.", "Early Cr.")

# Rename areas with bioregion names
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic", "Total")
row.names(residence_times_per_areas_perc_ggplot) <- bioregion_names

residence_times_per_areas_perc_ggplot <- reshape2::melt(residence_times_per_areas_perc_ggplot)
names(residence_times_per_areas_perc_ggplot) <- c("Bioregions", "Time_strata", "Lineage_counts")
residence_times_per_areas_perc_ggplot

# Save ggplot df for Raw residence times per bioregions, per geological era
saveRDS(residence_times_per_areas_perc_ggplot, file = "./outputs/LTT/residence_times_per_areas_perc_ggplot.rds")
residence_times_per_areas_perc_ggplot <- readRDS(file = "./outputs/LTT/residence_times_per_areas_perc_ggplot.rds")


## Set colors and legend

# Extract areas
all_areas <- dimnames(residence_times_per_areas_all_BSM_maps_summary_array)[[2]]
# Remove total
all_areas <- all_areas[all_areas != "total"]

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
names(colors_list_for_areas) <- bioregion_names[-length(bioregion_names)]

# Set levels of areas in defined order
residence_times_per_areas_perc_ggplot$Bioregions <- factor(residence_times_per_areas_perc_ggplot$Bioregions, levels = bioregion_names)

# Set levels of geological eras/time-strata in defined order (reverse order)
residence_times_per_areas_perc_ggplot$Time_strata <- factor(residence_times_per_areas_perc_ggplot$Time_strata, levels = rev(levels(residence_times_per_areas_perc_ggplot$Time_strata)))

# Remove the "Total" values
remove_total_indices <- residence_times_per_areas_perc_ggplot$Bioregions != "Total"

View(residence_times_per_areas_perc_ggplot)

## GGplot

pdf(file = "./outputs/LTT/Residence_times_per_areas_perc_ggplot.pdf", width = 12, height = 8)
ggplot_residence_times_per_areas_perc <- ggplot(data = residence_times_per_areas_perc_ggplot[remove_total_indices, ], mapping = aes(x = Time_strata, y = Lineage_counts, fill = Bioregions)) +
  
  ## Plot bars
  geom_bar(position = "stack", stat = "identity", width = 0.5) +
  
  # # Reverse time axis
  # scale_x_reverse() +
  
  # Set fill colors
  scale_fill_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "Residence times per bioregion") +
  xlab(label = "Geological periods") +
  ylab(label = "Residence time  [My]") +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, size = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 13)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 1.0),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 8)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_residence_times_per_areas_perc)

dev.off()


##### 5/ Extract and standardize data for lineages through time per bioregions across BS maps #####

### 5.1/ Convert simmap into LTT ####

# Load phytools.simmaps
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/DEC_J_simmaps.rds")

source(file = "./functions/ltt_simmap_custom.R")

# Convert simmaps into LTT plot data
DEC_J_LTT_all_maps <- ltt_multiSimmap_custom(DEC_J_simmaps)

# Save LTT plot data
saveRDS(object = DEC_J_LTT_all_maps, file = "./outputs/LTT/DEC_J_LTT_all_maps.rds")

### 5.2/ Extract LTT data across ranges ####

# Load LTT plot data
DEC_J_LTT_all_maps <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_maps.rds")

# Initiate df

DEC_J_LTT_all_ranges_ggplot <- data.frame()

# Loop per ltt to extract data

for (i in 1:length(DEC_J_LTT_all_maps)) 
{
  # i <- 2
  
  DEC_J_LTT_ggplot_i <- as.data.frame(DEC_J_LTT_all_maps[[i]]$ltt) 
  DEC_J_LTT_ggplot_i$branch_ID <- str_remove(string = row.names(DEC_J_LTT_ggplot_i), pattern = "X")
  DEC_J_LTT_ggplot_i$branch_ID[1] <- "root"
  DEC_J_LTT_ggplot_i$branch_ID[nrow(DEC_J_LTT_ggplot_i)] <- "tips"
  DEC_J_LTT_ggplot_i$time <- max(DEC_J_LTT_all_maps[[i]]$times) - DEC_J_LTT_all_maps[[i]]$times
  DEC_J_LTT_ggplot_i <- pivot_longer(data = DEC_J_LTT_ggplot_i, cols = 1:ncol(DEC_J_LTT_all_maps[[i]]$ltt), names_to = "states", values_to = "counts")
  DEC_J_LTT_ggplot_i$map <- i
  
  # Bind df
  DEC_J_LTT_all_ranges_ggplot <- rbind(DEC_J_LTT_all_ranges_ggplot, DEC_J_LTT_ggplot_i)
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - LTT data extracted from map n°", i, "/", length(DEC_J_LTT_all_maps),"\n"))
  }
}

dim(DEC_J_LTT_all_ranges_ggplot)
table(DEC_J_LTT_all_ranges_ggplot$states) # Nb of occurrences varies because not all states/ranges are present in all BSM

View(DEC_J_LTT_all_ranges_ggplot)

# Save data frame of LTT per ranges, per BS maps
saveRDS(object = DEC_J_LTT_all_ranges_ggplot, file = "./outputs/LTT/DEC_J_LTT_all_ranges_ggplot.rds")

# Remove raw LTT data to clean RAM
# rm(DEC_J_LTT_all_maps) ; gc()

### 5.3/ Standardize time steps across BS maps ####

# Each BS maps have different timing of transitions for anagenetic events, thus times in LTT are different
# Need to equalize them by sampling lineage counts to the same time steps
# For each map use the closest time steps that is younger than the targeted time step

# Subset for time every 0.1 My

time_step <- 0.1

root_age <- round(max(DEC_J_LTT_all_ranges_ggplot$time), 0)
target_times <- seq(from = 0, to = root_age, by = time_step)
maps_indices <- unique(DEC_J_LTT_all_ranges_ggplot$map)

DEC_J_LTT_all_ranges_time_std_ggplot <- data.frame()
for (i in maps_indices)
{
  # i <- 1
  
  # Extract map data
  LTT_df <- DEC_J_LTT_all_ranges_ggplot[DEC_J_LTT_all_ranges_ggplot$map == i, ]
  
  # Extract available times for this BS map
  available_times_i <- unique(LTT_df$time)
  
  # Compare available times to target times
  diff_times_df_i <- tidyr::expand_grid(target_times= target_times, available_times = available_times_i)
  diff_times_df_i$diff <- (diff_times_df_i$available_times - diff_times_df_i$target_times)
  
  # Find the closest available time that is >= to each target time
  selected_times_df_i <- diff_times_df_i %>%
    filter(!(diff < 0)) %>%
    group_by(target_times) %>%
    arrange(diff) %>%
    mutate(rank = row_number()) %>% 
    filter(rank == 1) %>%
    arrange(target_times)
  
  # Replace available times by target times in LTT data
  LTT_df_target_time <- selected_times_df_i %>% 
    left_join(y = LTT_df, by = c("available_times" = "time"), relationship = "many-to-many") %>%
    rename(time = target_times) %>%
    select(branch_ID, time, states, counts, map)
  
  # Check that all target times have one value per state
  # table(LTT_df_target_time$time)
  # as.numeric(names(table(table(LTT_df_target_time$time)))) == length(unique(LTT_df$states))
  
  # Append final df
  DEC_J_LTT_all_ranges_time_std_ggplot <- rbind(DEC_J_LTT_all_ranges_time_std_ggplot, LTT_df_target_time)
  
  # Print progress every 10 maps
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Time in LTT data standardized for map n°", i, "/", length(maps_indices),"\n"))
  }
  
}

# Arrange data in initial order
DEC_J_LTT_all_ranges_time_std_ggplot <- DEC_J_LTT_all_ranges_time_std_ggplot %>% 
  arrange(map, desc(time), states)

# Save data frame of LTT per ranges, per BS maps
saveRDS(object = DEC_J_LTT_all_ranges_time_std_ggplot, file = "./outputs/LTT/DEC_J_LTT_all_ranges_time_std_ggplot.rds")

# Remove raw LTT data with no time standardization to clean RAM
# rm(DEC_J_LTT_all_ranges_ggplot) ; gc()


### 5.4/ Aggregate ranges into unique areas ####

## Divide counts of lineages within ranges equally across areas

all_ranges <- unique(DEC_J_LTT_all_ranges_time_std_ggplot$states)

DEC_J_LTT_all_ranges_split_ggplot <- DEC_J_LTT_all_ranges_time_std_ggplot
for (i in seq_along(all_ranges))
{
  # i <- 10
  
  # Extract focal range
  range_i <- all_ranges[i]
  
  if (range_i != "total")
  {
    # Find matching entries
    matching_entries_i <- DEC_J_LTT_all_ranges_split_ggplot$states == range_i
    
    # Extract number of areas composing the range
    nb_areas_i <- nchar(range_i)
    
    # Divide nb of lineages accordingly
    DEC_J_LTT_all_ranges_split_ggplot$counts[matching_entries_i] <- DEC_J_LTT_all_ranges_split_ggplot$counts[matching_entries_i] / nb_areas_i
  }
  
  # Print progress every range
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Nb of lineages split across areas for range ",range_i,": n°", i, "/", length(all_ranges),"\n"))
  }
}

View(DEC_J_LTT_all_ranges_split_ggplot)

## Split ranges in unique areas
ranges_split_list <- str_split(string = all_ranges, pattern = "")

## Loop per areas to aggregate data
DEC_J_LTT_all_areas_ggplot <- data.frame()
for (i in seq_along(all_areas))
{
  # i <- 1
  
  # Extract focal area
  areas_i <- all_areas[i]
  
  # Find matching ranges
  matching_ranges_i <- all_ranges[str_detect(string = ranges_split_list, pattern = areas_i)]
  
  # Compute sum across all matching ranges for a given map x time
  LTT_sum_areas_i <- DEC_J_LTT_all_ranges_split_ggplot %>% 
    filter(states %in% matching_ranges_i) %>%
    group_by(map, time) %>%
    summarize(counts = sum(counts)) %>% 
    mutate(areas = areas_i) %>% 
    arrange(map, desc(time), areas) %>%
    ungroup()
    
  # View(LTT_sum_areas_i)
  
  # Append final data frame
  DEC_J_LTT_all_areas_ggplot <- rbind(DEC_J_LTT_all_areas_ggplot, LTT_sum_areas_i)
  
  # Print progress every range
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Nb of lineages aggregated across area ",areas_i,": n°", i, "/", length(all_areas),"\n"))
  }
  
}

# Copy total data
LTT_sum_total <- DEC_J_LTT_all_ranges_split_ggplot %>% 
  filter(states == "total") %>% 
  mutate(areas = states) %>%
  select(map, time, counts, areas) %>% 
  arrange(map, desc(time), areas)

DEC_J_LTT_all_areas_ggplot <- rbind(DEC_J_LTT_all_areas_ggplot, LTT_sum_total)

# Arrange data in initial order
DEC_J_LTT_all_areas_ggplot <- DEC_J_LTT_all_areas_ggplot %>% 
  arrange(map, desc(time), areas)

# Final nb of raws = nb of areas * nb of time steps * nb of maps
nrow(DEC_J_LTT_all_areas_ggplot) == (length(all_areas) + 1) * (length(target_times) - 1) * length(maps_indices)
View(DEC_J_LTT_all_areas_ggplot)

# Save data frame of LTT per areas, per BS maps
saveRDS(object = DEC_J_LTT_all_areas_ggplot, file = "./outputs/LTT/DEC_J_LTT_all_areas_ggplot.rds")


### 5.5/ Compute mean across all BS maps ####

# Compute mean across BS maps
DEC_J_LTT_all_areas_mean_ggplot <- DEC_J_LTT_all_areas_ggplot %>% 
  group_by(time, areas) %>%
  summarize(mean_counts = mean(counts)) %>% 
  ungroup

# Add percentages
DEC_J_LTT_all_areas_mean_ggplot <- DEC_J_LTT_all_areas_mean_ggplot %>% 
  group_by(time) %>%
  mutate(mean_percentages = 100 * mean_counts / sum(mean_counts) * 2) %>%  # Times 2 because total is already present in the sum
  ungroup

# Should have as many values as areas per time steps
# table(table(DEC_J_LTT_all_areas_mean_ggplot$time))

# Save data frame of LTT per areas, mean across all BS maps
saveRDS(object = DEC_J_LTT_all_areas_mean_ggplot, file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")


##### 6/ Plot lineages through time per bioregions #####

### 6.1/ Plot raw counts and CI using GGplot, with Total ####

# Load data frame of LTT per areas, per BS maps
DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
colors_list_for_areas <- c(colors_list_for_areas, "grey90") # Add grey for total
names(colors_list_for_areas)[length(colors_list_for_areas)] <- "total"
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic", "Total")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_ggplot$areas <- factor(DEC_J_LTT_all_areas_ggplot$areas, levels = c(all_areas, "total"))


View(DEC_J_LTT_all_areas_ggplot)

pdf(file = "./outputs/LTT/DEC_J_LTT_all_areas_raw_ggplot.pdf", width = 12, height = 8)

ggplot_DEC_J_LTT_all_areas_raw <- ggplot(data = DEC_J_LTT_all_areas_ggplot, mapping = aes(x = time, col = areas)) +
  
  ## Plot mean values
  # geom_step(data = DEC_J_LTT_all_areas_mean_ggplot, mapping = aes(y = mean_counts), linewidth = 1.5, alpha = 1) +
  # geom_smooth(mapping = aes(group = areas, y = counts), method = "loess", se = FALSE, linewidth = 1.5, alpha = 0.1) +
  # geom_smooth(mapping = aes(group = areas, y = counts), method = "loess", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  ## Plot all values
  geom_step(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 2, alpha = 0.1) +
  # geom_step(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 0.8, alpha = 0.1) +
  # geom_line(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 0.8, alpha = 0.1) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "red", linewidth = 1.5, linetype = 2) +
  
  # Set colors
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "LTT per bioregion") +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages") +
  
  # Adjust legend alpha
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_raw)

dev.off()

### 6.2/ Plot raw counts and CI using GGplot, without Total ####

# Load data frame of LTT per areas, per BS maps
DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_ggplot$areas <- factor(DEC_J_LTT_all_areas_ggplot$areas, levels = c(all_areas, "total"))

# Remove the "total" values
remove_total_indices <- DEC_J_LTT_all_areas_ggplot$areas != "total"


pdf(file = "./outputs/LTT/DEC_J_LTT_all_areas_raw_without_total_ggplot.pdf", width = 12, height = 8)

ggplot_DEC_J_LTT_all_areas_raw_without_total <- ggplot(data = DEC_J_LTT_all_areas_ggplot[remove_total_indices, ], mapping = aes(x = time, col = areas)) +
  
  ## Plot mean values
  # geom_step(data = DEC_J_LTT_all_areas_mean_ggplot, mapping = aes(y = mean_counts), linewidth = 1.5, alpha = 1) +
  # geom_smooth(mapping = aes(group = areas, y = counts), method = "loess", se = FALSE, linewidth = 1.5, alpha = 0.1) +
  # geom_smooth(mapping = aes(group = areas, y = counts), method = "loess", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  ## Plot all values
  geom_step(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 2, alpha = 0.1) +
  # geom_step(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 0.8, alpha = 0.1) +
  # geom_line(mapping = aes(group = interaction(map, areas), y = counts), linewidth = 0.8, alpha = 0.1) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "red", linewidth = 1.5, linetype = 2) +
  
  # Set colors
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "LTT per bioregion") +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages") +
  
  # Adjust legend alpha
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_raw_without_total)

dev.off()


### 6.3/ Plot cumulative percentages using GGplot ####

# Load data frame of LTT per areas, per BS maps
DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_mean_ggplot$areas <- factor(DEC_J_LTT_all_areas_mean_ggplot$areas, levels = c(all_areas, "total"))

# Remove the "total" values
remove_total_indices <- DEC_J_LTT_all_areas_mean_ggplot$areas != "total"


## GGplot

pdf(file = "./outputs/LTT/DEC_J_LTT_all_areas_cumperc_ggplot.pdf", width = 12, height = 8)

ggplot_DEC_J_LTT_all_areas_cumperc <- ggplot(data = DEC_J_LTT_all_areas_mean_ggplot[remove_total_indices, ], mapping = aes(x = time, y = mean_percentages, fill = areas)) +
  
  ## Plot bars
  geom_bar(position = "stack", stat = "identity", width = time_step) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "black", linewidth = 1.5, linetype = 2) +
  
  # Set fill colors
  scale_fill_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "LTT per bioregion") +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages  [%]") +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_cumperc)

dev.off()


### 6.4/ Plot smoothed cumulative percentages using GGplot ####

## Same than 6.3 but using a sliding window to compute smooth values

# Load data frame of LTT per areas, per BS maps
DEC_J_LTT_all_areas_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_mean_ggplot$areas <- factor(DEC_J_LTT_all_areas_mean_ggplot$areas, levels = c(all_areas, "total"))

## Set sliding window properties

time_step <- 0.1  # Steps = 0.1 My (Must be a divider of window_width)
window_width <- 4 # Width = 4 My (Must be a multiple of time_steps)

# Extract root age
root_age <- max(DEC_J_LTT_all_areas_mean_ggplot$time)

# Generate start and end times of sliding windows
start_time_list <- seq(from = 0, to = root_age - window_width, by = time_step)
end_time_list <- seq(from = window_width, to = root_age, by = time_step)

## Loop per sliding windows to compute average counts/percentages
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- data.frame()
for (i in seq_along(start_time_list))
{
  # i <- 1
  
  # Extract start/end times
  start_time_i <- start_time_list[i]
  end_time_i <- end_time_list[i]
  mean_time_i <- (end_time_i + start_time_i) / 2
  
  # Extract all mean values within the time frame
  DEC_J_LTT_all_areas_mean_ggplot_i <- DEC_J_LTT_all_areas_mean_ggplot %>% 
    filter((time >= start_time_i) & (time <= end_time_i)) %>% 
    group_by(areas) %>% 
    summarize(mean_counts = mean(mean_counts),
              mean_percentages = mean(mean_percentages)) %>%
    mutate(smoothed_time = mean_time_i)
  
  # Append to final data.frame
  DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- rbind(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot, DEC_J_LTT_all_areas_mean_ggplot_i)
}

View(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot)

# Save LTT smoothed data from sliding windows
saveRDS(object = DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot, file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot.rds")

# Remove the "total" values
remove_total_indices <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$areas != "total"


## GGplot

pdf(file = "./outputs/LTT/DEC_J_LTT_all_areas_cumperc_per_sliding_windows_ggplot.pdf", width = 12, height = 8)

ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows <- ggplot(data = DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot[remove_total_indices, ], mapping = aes(x = smoothed_time, y = mean_percentages, fill = areas)) +
  
  ## Plot bars
  geom_bar(position = "stack", stat = "identity", width = time_step) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "black", size = 1.5, linetype = 2) +
  
  # Set fill colors
  scale_fill_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = paste0("LTT per bioregion\nSliding window = ",window_width," My")) +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages  [%]") +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows)

dev.off()


##### 7/ Make animated plots #####

## gganimate

# transition_*() defines how the data should be spread out and how it relates to itself across time.
# view_*() defines how the positional scales should change along the animation.
# shadow_*() defines how data from other points in time should be presented in the given point in time.
# enter_*()/exit_*() defines how new data should appear and how old data should disappear during the course of the animation.
# ease_aes() defines how different aesthetics should be eased during transitions = type of trajectories between two frames

?ease_aes
?transition_time

### 7.1/ Animated LTT plot with total ####

## Use only the mean data, otherwise to long to generate!

# Load data frame of mean LTT per areas, per BS maps
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
colors_list_for_areas <- c(colors_list_for_areas, "grey90") # Add grey for total
names(colors_list_for_areas)[length(colors_list_for_areas)] <- "total"
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic", "Total")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_mean_ggplot$areas <- factor(DEC_J_LTT_all_areas_mean_ggplot$areas, levels = c(all_areas, "total"))


ggplot_DEC_J_LTT_all_areas_mean_raw <- ggplot(data = DEC_J_LTT_all_areas_mean_ggplot, mapping = aes(x = time, y = mean_counts, col = areas)) +
  
  ## Plot mean values
  # geom_step(mapping = aes(group = areas), linewidth = 2, alpha = 0.8) +
  geom_line(mapping = aes(group = areas), linewidth = 2, alpha = 0.8) +
  geom_point(mapping = aes(group = areas), size = 4) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "black", linewidth = 1.2, linetype = 2) +
  
  # Set colors
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "LTT per bioregion") +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages") +
  
  # Adjust legend alpha
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_mean_raw)

## Add animations

ggplot_DEC_J_LTT_all_areas_mean_raw_animated <- ggplot_DEC_J_LTT_all_areas_mean_raw +
  
  # Add transition between time
  # transition_time(-time) +
  transition_reveal(-time) +
  # transition_states(states = time, transition_length = 1, state_length = 0) +
  
  # Use quadratic transitions between frames
  ease_aes("quadratic-in-out") +
  
  # Adjust title label accordingly
  # labs(title = "LTT per bioregion\nTime: {-frame_time} My", x = "Time  [My]", y = "Lineages") +
  labs(title = "LTT per bioregion",
       subtitle = "Time: {round(-frame_along, 0)} My ago",
       x = "Time  [My]", y = "Lineages") +
  # labs(title = "LTT per bioregion\nTime: {-closest_frame} My", x = "Time  [My]", y = "Lineages") +
  
  # Aesthetics
  theme(plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 0, t = 5)),
        plot.subtitle = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines")) # trbl

print(ggplot_DEC_J_LTT_all_areas_mean_raw_animated)

## Save animation as a GIF
anim_save(filename = "./outputs/LTT/DEC_J_LTT_all_areas_mean_raw_ggplot.gif",
          animation = ggplot_DEC_J_LTT_all_areas_mean_raw_animated,
          height = 8, width = 12, units = "in", # In inches
          res = 300, # dpi = 300
          fps = 10, # 10 My/s Total = 10s
          renderer = gifski_renderer(loop = FALSE)) # Infinite loop?


### 7.2/ Animated LTT plot without total ####

## Use only the mean data, otherwise to long to generate!

# Load data frame of mean LTT per areas, per BS maps
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Use color scheme of BSM
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")

# Set colors for areas/bioregions
colors_list_for_areas <- colors_list_for_states[all_areas]
colors_list_for_areas <- c(colors_list_for_areas, "grey90") # Add grey for total
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")

# Set levels of areas in defined order
DEC_J_LTT_all_areas_mean_ggplot$areas <- factor(DEC_J_LTT_all_areas_mean_ggplot$areas, levels = c(all_areas, "total"))

# Remove the "total" values
remove_total_indices <- DEC_J_LTT_all_areas_mean_ggplot$areas != "total"


ggplot_DEC_J_LTT_all_areas_mean_raw_without_total <- ggplot(data = DEC_J_LTT_all_areas_mean_ggplot[remove_total_indices, ], mapping = aes(x = time, y = mean_counts, col = areas)) +
  
  ## Plot mean values
  # geom_step(mapping = aes(group = areas), linewidth = 2, alpha = 0.8) +
  geom_line(mapping = aes(group = areas), linewidth = 2, alpha = 0.8) +
  geom_point(mapping = aes(group = areas), size = 4) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "black", linewidth = 1.2, linetype = 2) +
  
  # Set colors
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = "LTT per bioregion") +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages") +
  
  # Adjust legend alpha
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_mean_raw_without_total)


## Add animations

ggplot_DEC_J_LTT_all_areas_mean_raw_without_total_animated <- ggplot_DEC_J_LTT_all_areas_mean_raw_without_total +
  
  # Add transition between time
  # transition_time(-time) +
  transition_reveal(-time) +
  # transition_states(states = time, transition_length = 1, state_length = 0) +
  
  # Use quadratic transitions between frames
  ease_aes("quadratic-in-out") +
  
  # Adjust title label accordingly
  # labs(title = "LTT per bioregion\nTime: {-frame_time} My", x = "Time  [My]", y = "Lineages") +
  labs(title = "LTT per bioregion",
       subtitle = "Time: {round(-frame_along, 0)} My ago",
       x = "Time  [My]", y = "Lineages") +
  # labs(title = "LTT per bioregion\nTime: {-closest_frame} My", x = "Time  [My]", y = "Lineages") +
  
  # Aesthetics
  theme(plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 0, t = 5)),
        plot.subtitle = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines")) # trbl

print(ggplot_DEC_J_LTT_all_areas_mean_raw_without_total_animated)

## Save animation as a GIF
anim_save(filename = "./outputs/LTT/DEC_J_LTT_all_areas_mean_raw_without_total_ggplot.gif",
          animation = ggplot_DEC_J_LTT_all_areas_mean_raw_without_total_animated,
          height = 8, width = 12, units = "in", # In inches
          res = 300, # dpi = 300
          fps = 10, # 10 My/s Total = 10s
          renderer = gifski_renderer(loop = FALSE)) # Infinite loop?

### 7.3/ Animated smoothed Cumulative %  ####

# Load LTT smoothed data from sliding windows
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot.rds")

## Add values for extreme times (close to present and root time)

# Create missing time sequence
root_time <- max(DEC_J_LTT_all_areas_mean_ggplot$time)
max_time <- max(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$smoothed_time)
min_time <- min(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$smoothed_time)
additional_time_past <- seq(from = max_time + 0.1, to = root_time, by = 0.1)
additional_time_present <- seq(from = 0, to = min_time - 0.1, by = 0.1)

# Create missing data for past time (110.9 to 112.9)
additional_time_past_to_add <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot[DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$smoothed_time == max_time, ] 
additional_time_past_df <- data.frame()
for (i in 1:length(additional_time_past))
{
  additional_time_past_df <- rbind(additional_time_past_df, additional_time_past_to_add)
}
additional_time_past_df$smoothed_time <- rep(x = additional_time_past, each = nrow(additional_time_past_to_add))

# Create missing data for present time (0 to 1.9)
additional_time_present_to_add <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot[DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$smoothed_time == min_time, ] 
additional_time_present_df <- data.frame()
for (i in 1:length(additional_time_present))
{
  additional_time_present_df <- rbind(additional_time_present_df, additional_time_present_to_add)
}
additional_time_present_df$smoothed_time <- rep(x = additional_time_present, each = nrow(additional_time_present_to_add))

# Append missing data for extreme times
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- rbind(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot, additional_time_past_df)
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- rbind(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot, additional_time_present_df)

# Arrange data in initial order
DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot %>% 
  arrange(desc(smoothed_time), areas)

# Remove the "total" values
remove_total_indices <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$areas != "total"

## Compute cumulative values for each bioregion in the order of the stack, to use for geom_point
DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot <- DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot
areas_lvl <- levels(DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot$areas)
areas_lvl <- areas_lvl[areas_lvl != "total"]
areas_lvl <- rev(areas_lvl)
current_counts <- 0 ; current_perc <- 0
DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$cumsum_counts <- NA
DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$cumsum_perc <- NA
for (i in seq_along(areas_lvl))
{
  # i <- 1
  
  # Extract area
  areas_i <- areas_lvl[i]
  areas_indices <- DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$areas == areas_i
  
  # Extract values for focal area
  focal_counts <- DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$mean_counts[areas_indices]
  focal_perc <- DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$mean_percentages[areas_indices]
  
  # Compute new cumulative sum
  current_counts <- current_counts + focal_counts
  current_perc <- current_perc + focal_perc
  
  # Store cumulative sum
  DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$cumsum_counts[areas_indices] <- current_counts
  DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot$cumsum_perc[areas_indices] <- current_perc
}


## GGplot
ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows <- ggplot(data = DEC_J_LTT_all_areas_mean_per_sliding_windows_ggplot[remove_total_indices, ], mapping = aes(x = smoothed_time, y = mean_percentages, fill = areas)) +
  
  ## Plot bars
  geom_bar(position = "stack", stat = "identity", width = time_step*12) + # Need to enlarge widht because only 1/10 time step are shown in the GIF
  
  ## Plot points
  # geom_line(mapping = aes(group = areas), linewidth = 2, alpha = 0.8) +
  geom_point(data = DEC_J_LTT_all_areas_mean_per_sliding_windows_cumsum_ggplot[remove_total_indices, ], mapping = aes(x = smoothed_time, y = cumsum_perc, group = areas, col = areas), size = 4) +
  
  # Reverse time axis
  scale_x_reverse() +
  
  # Add time periods
  geom_vline(xintercept = DEC_J_fit$inputs$timeperiods[-length(DEC_J_fit$inputs$timeperiods)], col = "black", size = 1.5, linetype = 2) +
  
  # Set fill colors
  scale_fill_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  # Set col colors
  scale_color_manual("Bioregions", labels = bioregion_names, values = colors_list_for_areas) +
  
  # Set labels
  ggtitle(label = paste0("LTT per bioregion")) +
  xlab(label = "Time  [My]") +
  ylab(label = "Lineages  [%]") +
  
  # Aesthetics
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(0.5, "line"),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows)

## Add animations

ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows_animated <- ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows +
  
  # Add transition between time
  transition_time(-smoothed_time) +
  # transition_reveal(-smoothed_time) +
  # transition_states(states = smoothed_time, transition_length = 1, state_length = 0) +
  
  # Keep previous values, 
  shadow_mark(exclude_layer = 2) +

  # Use quadratic transitions between frames
  ease_aes("quadratic-in-out") +
  
  # Adjust title label accordingly
  labs(title = "LTT per bioregion",
       subtitle = "Time: {round(-frame_time, 0)} My ago",
       x = "Time  [My]", y = "Lineages  [%]") +
  # labs(title = "LTT per bioregion",
  #      subtitle = "Time: {round(-frame_along, 0)} My ago",
  #      x = "Time  [My]", y = "Lineages  [%]") +
  # labs(title = "LTT per bioregion\nTime: {-closest_frame} My", x = "Time  [My]", y = "Lineages") +
  
  # Aesthetics
  theme(plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 0, t = 5)),
        plot.subtitle = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines")) # trbl

print(ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows_animated)

## Save animation as a GIF
anim_save(filename = "./outputs/LTT/DEC_J_LTT_all_areas_cumperc_per_sliding_windows_ggplot.gif",
          animation = ggplot_DEC_J_LTT_all_areas_mean_per_sliding_windows_animated,
          height = 8, width = 12, units = "in", # In inches
          res = 300, # dpi = 300
          fps = 10, # 10 My/s Total = 10s
          renderer = gifski_renderer(loop = FALSE)) # Infinite loop?
           

#### Improve color scheme for all plots !!! #####


