##### Script 12: Compare origin of lineages: in situ speciation vs. dispersal #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compare origin of lineages: in situ speciation vs. dispersal along time
# Plot evolution of proportion of origins in time, overall and per bioregions
# Map mean immigration age of taxa
# Test for latitudinal and longitudinal gradients of immigration age

###

### Inputs

# Simmaps of Biogeographic Stochastic Maps

###

### Outputs

## Plot evolution of proportion of origins in time
  # Current time = stacked barplot of Counts and Percentages per Bioregions
  # Along time = lines of % per bioregions, including total
## Map mean immigration age of taxa
## Plot latitudinal and longitudinal mean immigration age
  # Add null model values from permuted data

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
library(magick)    #  For animated GIF
library(phangorn)
library(BayesTwin)  # For 95% HPD
library(raster)

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load modeling results from best model (DEC+J) #####

# Load time-stratified DEC+J model output
# DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/Ponerinae_MCC_phylogeny_1534t/DEC_J_fit.rds")

# Extract areas list
returned_mats <- BioGeoBEARS::get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
areas_list <- returned_mats$areanames

# Load simmaps of BS maps
# DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps.rds")
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps.rds")

### 1.3/ Set time steps ####

window_steps <- 1 # Steps = 1 My

# Extract root age
root_age <- max(phytools::nodeHeights(DEC_J_simmaps[[1]]))

# Generate time steps
time_steps <- seq(from = 0, to = root_age, by = window_steps)

### 1.4/ Function for ggplot aesthetics ####

# Aesthetics for barplots
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


##### 2/ Extract edge states and origin from simmaps #####

source("./functions/compute_lineage_origins_per_time_steps.R")

### 2.1/ Extract edge states and origin from simmaps ####

areas_endemicity_all_time_all_maps_array <- find_egde_origins_on_time_steps_on_MultiSimmap(multiSimmap = DEC_J_simmaps,
                                                       time_steps = time_steps,
                                                       areas_list = areas_list,
                                                       melted = F,  # Should the final array be melted into a df?
                                                       verbose = T)

# Save final array
# saveRDS(areas_endemicity_all_time_all_maps_array, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/areas_endemicity_all_time_all_maps_array.rds")
saveRDS(areas_endemicity_all_time_all_maps_array, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/areas_endemicity_all_time_all_maps_array.rds")

### 2.2/ Melt results for ggplot ####

# areas_endemicity_all_time_all_maps_array <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/areas_endemicity_all_time_all_maps_array.rds")
areas_endemicity_all_time_all_maps_array <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/areas_endemicity_all_time_all_maps_array.rds")

areas_endemicity_all_time_all_maps_df <- reshape2::melt(areas_endemicity_all_time_all_maps_array)
names(areas_endemicity_all_time_all_maps_df) <- c("area", "time", "map", "endemicity")

# Save melted df of endemicity per bioregions x time x maps
# saveRDS(areas_endemicity_all_time_all_maps_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/areas_endemicity_all_time_all_maps_df.rds")
saveRDS(areas_endemicity_all_time_all_maps_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/areas_endemicity_all_time_all_maps_df.rds")

### 2.3/ Aggregate along maps ####

mean_areas_endemicity_all_time_df <- areas_endemicity_all_time_all_maps_df %>% 
  group_by(area, time) %>% 
  summarize(mean_endemicity = mean(endemicity, na.rm = T)) %>% 
  ungroup()

# Replace Nan by NA
mean_areas_endemicity_all_time_df$mean_endemicity[is.nan(mean_areas_endemicity_all_time_df$mean_endemicity)] <- NA

# Save melted df of mean endemicity per bioregions x time
# saveRDS(mean_areas_endemicity_all_time_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/mean_areas_endemicity_all_time_df.rds")
saveRDS(mean_areas_endemicity_all_time_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/mean_areas_endemicity_all_time_df.rds")


### 2.4/ Compute total counts and percentages in current time ####

# Load LTT data per bioregions along time
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

names(DEC_J_LTT_all_areas_mean_ggplot)

# Merge with LTT data
# Extract only current time
# Compute counts and percentage of regional taxa vs. exotic taxa
mean_areas_endemicity_current_time_df <- mean_areas_endemicity_all_time_df %>% 
  filter(time == 0) %>%  # Extract only current time
  left_join(x = ., y = DEC_J_LTT_all_areas_mean_ggplot, join_by(area == areas, time == time)) %>% # Merge with LTT data
  rename(area_richness = mean_counts) %>%
  mutate(regional_counts = round(mean_endemicity * area_richness, 1), # Compute counts and percentage of regional taxa vs. exotic taxa
         regional_perc = round(mean_endemicity * 100, 1),
         exotic_counts = round(area_richness - regional_counts, 1),
         exotic_perc = round(100 - regional_perc, 1)) %>%
  pivot_longer(cols = c(regional_counts, regional_perc, exotic_counts, exotic_perc, ), names_to = c("origin", ".value"), names_sep = "_", values_to = "value")
  

# Save melted df of mean counts and percentages of endemicity per bioregions in current time
# saveRDS(mean_areas_endemicity_current_time_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/mean_areas_endemicity_current_time_df.rds")
saveRDS(mean_areas_endemicity_current_time_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/mean_areas_endemicity_current_time_df.rds")


##### 3/ Plot counts of lineages per origin in current time ####

# Load melted df of mean counts and percentages of endemicity per bioregions in current time
# mean_areas_endemicity_current_time_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/mean_areas_endemicity_current_time_df.rds")
mean_areas_endemicity_current_time_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/mean_areas_endemicity_current_time_df.rds")

mean_areas_endemicity_current_time_df_no_total <- mean_areas_endemicity_current_time_df %>% 
  filter(area != "total")

# Add % as labels
mean_areas_endemicity_current_time_df_no_total$perc_labels <- paste0(round(mean_areas_endemicity_current_time_df_no_total$perc, 0), "%")
mean_areas_endemicity_current_time_df_no_total$perc_labels_filtered <- mean_areas_endemicity_current_time_df_no_total$perc_labels
mean_areas_endemicity_current_time_df_no_total$perc_labels_filtered[mean_areas_endemicity_current_time_df_no_total$origin == "Dispersal"] <- NA

# Adjust color scheme for bioregions
areas_list <- c("A", "U", "I", "R", "N", "E", "W")
colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[areas_list]

# Adjust color scheme for borders (interaction between origin and area)
colors_list_for_areas_in_situ <- colors_list_for_areas
colors_rgb_array <- col2rgb(colors_list_for_areas)
colors_list_for_areas_dispersal <- colors_list_for_areas
for (i in 1:ncol(colors_rgb_array))
{
  colors_list_for_areas_dispersal[i] <- rgb(red = colors_rgb_array[1,i], green = colors_rgb_array[2,i], blue = colors_rgb_array[3,i], alpha = 0.4*255, maxColorValue = 255)
}
colors_borders_area_origin <- c(colors_list_for_areas_dispersal, colors_list_for_areas_in_situ)

# Adjust order of bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic", "Neotropics", "Eastern Palearctic", "Western Palearctic")
bioregion_names_reduced <- c("Afrotr.", "Austr.", "IndoM", "Nearct.", "Neotr.", "E-PA", "W-PA")
mean_areas_endemicity_current_time_df_no_total$area <- factor(mean_areas_endemicity_current_time_df_no_total$area, levels = areas_list, labels = bioregion_names_reduced)

# Adjust order of origins
origin_names <- c("Dispersal", "In situ diversification")
mean_areas_endemicity_current_time_df_no_total$origin <- factor(mean_areas_endemicity_current_time_df_no_total$origin, levels = c("exotic", "regional"), labels = origin_names)


# GGplot
# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Barplot_lineage_origin_counts_per_bioregions_current_time.pdf", width = 10, height = 6)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Barplot_lineage_origin_counts_per_bioregions_current_time.pdf", width = 10, height = 6)

barplot_lineage_origin_counts_per_bioregions_current_time_ggplot <- ggplot(data = mean_areas_endemicity_current_time_df_no_total,
                                                                           mapping = aes(y = counts, x = area, fill = area)) +
  # Plot barplot
  geom_col(aes(alpha = origin,
               # col = interaction(area, origin)),
               col = area),
           position = "stack", width = 0.8, linewidth = 0.5) +
  
  # Add percentages as labels
  geom_text(aes(label = perc_labels_filtered),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = 2) +
  
  # Set plot title +
  ggtitle(label = paste0("Current lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("Number of lineages") +
  
  # Adjust fill color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # # Adjust border color scheme and legend
  # scale_color_manual("Origin", labels = origin_names, values = c("black", "grey80")) +
  scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  # scale_color_manual(values = unname(colors_borders_area_origin)) +
  
  # Adjust alpha scheme and legend
  scale_alpha_manual("Origin", labels = origin_names, values = c(0.4, 1)) +
  
  # Remove legend for interaction(area_origin)
  guides(color = "none",
         alpha = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_lineage_origin_counts_per_bioregions_current_time_ggplot <- barplot_lineage_origin_counts_per_bioregions_current_time_ggplot %>% 
  add_aesthetics_barplot()

# Remove vertical grid axes
barplot_lineage_origin_counts_per_bioregions_current_time_ggplot <- barplot_lineage_origin_counts_per_bioregions_current_time_ggplot + 
  theme(panel.grid.major.x = element_blank())

# Plot
print(barplot_lineage_origin_counts_per_bioregions_current_time_ggplot)

dev.off()


##### 4/ Plot % of lineages per origin in current time ####

# GGplot
# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Barplot_lineage_origin_percs_per_bioregions_current_time.pdf", width = 10, height = 6)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Barplot_lineage_origin_percs_per_bioregions_current_time.pdf", width = 10, height = 6)

barplot_lineage_origin_percs_per_bioregions_current_time_ggplot <- ggplot(data = mean_areas_endemicity_current_time_df_no_total,
                                                                           mapping = aes(y = perc, x = area, fill = area)) +
  # Plot barplot
  geom_col(aes(alpha = origin,
               # col = interaction(area, origin)),
               col = area),
           position = "stack", width = 0.8, linewidth = 0.5) +
  
  # Add percentages as labels
  geom_text(aes(label = perc_labels_filtered),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = 2) +
  
  # Set plot title +
  ggtitle(label = paste0("Current lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Bioregions") +
  ylab("% of lineages") +
  
  # Adjust fill color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  
  # Adjust border color scheme and legend
  # scale_color_manual("Origin", labels = origin_names, values = c("black", "grey80")) +
  scale_color_manual("Bioregions", labels = bioregion_names, values = unname(colors_list_for_areas)) +
  # scale_color_manual(values = unname(colors_borders_area_origin)) +
  
  # Adjust alpha scheme and legend
  scale_alpha_manual("Origin", labels = origin_names, values = c(0.4, 1)) +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # trbl

# Adjust aesthetics
barplot_lineage_origin_percs_per_bioregions_current_time_ggplot <- barplot_lineage_origin_percs_per_bioregions_current_time_ggplot %>% 
  add_aesthetics_barplot()

# Remove vertical grid axes
barplot_lineage_origin_percs_per_bioregions_current_time_ggplot <- barplot_lineage_origin_percs_per_bioregions_current_time_ggplot + 
  theme(panel.grid.major.x = element_blank())

# Plot
print(barplot_lineage_origin_percs_per_bioregions_current_time_ggplot)

dev.off()


##### 5/ Plot evolution of proportion of regional lineages along time #####

# Adjust Bioregions fill scale to add "Total"
colors_list_for_areas_with_total <- c(colors_list_for_areas, "grey80")
bioregion_names_with_total <- c(bioregion_names, "Total")

# Load LTT data per bioregions along time
# DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_rough_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/Ponerinae_MCC_phylogeny_1534t/DEC_J_LTT_all_areas_mean_ggplot.rds")

# Load melted df of endemicity per bioregions x time x maps
# areas_endemicity_all_time_all_maps_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/areas_endemicity_all_time_all_maps_df.rds")
areas_endemicity_all_time_all_maps_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/areas_endemicity_all_time_all_maps_df.rds")

# Load melted df of mean endemicity per bioregions x time
# mean_areas_endemicity_all_time_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/mean_areas_endemicity_all_time_df.rds")
mean_areas_endemicity_all_time_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/mean_areas_endemicity_all_time_df.rds")

# Remove values when there are less than 1 lineage in average in the bioregions (to avoid artifactual high percentage from anecdotal events)
min_richness_threshold <- 1

areas_endemicity_all_time_all_maps_df <- areas_endemicity_all_time_all_maps_df %>%  
  left_join(x = ., y = DEC_J_LTT_all_areas_mean_ggplot, join_by(area == areas, time == time)) %>% # Merge with LTT data
  rename(area_richness = mean_counts) %>%
  filter(area_richness >= min_richness_threshold) %>%
  dplyr::select(area, time, map, endemicity, area_richness)
mean_areas_endemicity_all_time_df <- mean_areas_endemicity_all_time_df %>%
  left_join(x = ., y = DEC_J_LTT_all_areas_mean_ggplot, join_by(area == areas, time == time)) %>% # Merge with LTT data
  rename(area_richness = mean_counts) %>%
  filter(area_richness >= min_richness_threshold) %>%
  dplyr::select(area, time, mean_endemicity, area_richness)

# Adjust min/max values to the plot limits to avoid artifact in polygon drawings
max_time <- 85
min_time <- 0

# Compute quantiles
areas_endemicity_all_time_all_maps_df_quantiles <- areas_endemicity_all_time_all_maps_df %>% 
  filter(time >= min_time) %>%
  filter(time <= max_time) %>%
  group_by(area, time) %>% 
  # Compute quantiles
  reframe(quant_percs = quantile(endemicity, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
  group_by(area, time) %>% 
  mutate(quantile = seq(from = 0, to = 1, by = 0.01),
         quantile_ID = c(1:50, 51, 50:1)) %>%
  filter(quantile_ID != 51) %>% # Remove median
  # Assign points ID (order for drawing the polygon)
  group_by(area, quantile_ID) %>%
  arrange(area, quantile_ID, quantile, time) %>%
  mutate(n_points = n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  # mutate(points_ID = c(1:361, 722:362)) %>%
  # Reorder by quantile_ID, poly_ID, points_ID to check conformity of polygons
  arrange(area, quantile_ID, poly_ID, points_ID) %>%
  ungroup()

# Remove most extreme quantiles to keep 95% CI (strictly a 94% CI)
areas_endemicity_all_time_all_maps_df_quantiles <- areas_endemicity_all_time_all_maps_df_quantiles %>% 
  filter(!(quantile %in% c(0, 0.01, 0.02, 0.98, 0.99, 1)))

# Order bioregions
areas_endemicity_all_time_all_maps_df$area <- factor(x = areas_endemicity_all_time_all_maps_df$area, levels = c(names(colors_list_for_areas), "total"), labels = c(bioregion_names, "Total"))
mean_areas_endemicity_all_time_df$area <- factor(x = mean_areas_endemicity_all_time_df$area, levels = c(names(colors_list_for_areas), "total"), labels = c(bioregion_names, "Total"))
areas_endemicity_all_time_all_maps_df_quantiles$area <- factor(x = areas_endemicity_all_time_all_maps_df_quantiles$area, levels = c(names(colors_list_for_areas), "total"), labels = c(bioregion_names, "Total"))

## 5.1/ GGplot: only mean lines ####

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps.pdf", width = 10, height = 6)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps.pdf", width = 10, height = 6)

lines_endemicity_per_bioregions_all_time_steps_ggplot <- ggplot(data = areas_endemicity_all_time_all_maps_df) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = areas_endemicity_all_time_all_maps_df,
  #             mapping = aes(y = endemicity*100, x = time, group = area, col = area, fill = area),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = areas_endemicity_all_time_all_maps_df,
  #           mapping = aes(y = endemicity*100, x = time, group = interaction(area, map), col = area),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # Plot mean lines only
  geom_line(data = mean_areas_endemicity_all_time_df,
            mapping = aes(y = mean_endemicity*100, x = time, group = area, col = area),
            alpha = 1.0,
            linewidth = 2.0,
            show.legend = T) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of regional lineages") +
  
  # Set Y-axis to all be the same across plots
  ylim(c(0, 100)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     limits = c(85, 0) # Set limits
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
lines_endemicity_per_bioregions_all_time_steps_ggplot <- lines_endemicity_per_bioregions_all_time_steps_ggplot %>% 
  add_aesthetics_barplot()

print(lines_endemicity_per_bioregions_all_time_steps_ggplot)

dev.off()


## 5.2/ GGplot: with all lines ####

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps.pdf", width = 10, height = 6)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps_fuzzy_lines.pdf", width = 10, height = 6)

lines_endemicity_per_bioregions_all_time_steps_ggplot <- ggplot(data = areas_endemicity_all_time_all_maps_df) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = areas_endemicity_all_time_all_maps_df,
  #             mapping = aes(y = endemicity*100, x = time, group = area, col = area, fill = area),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # Plot 1000 replicates
  geom_line(data = areas_endemicity_all_time_all_maps_df,
            mapping = aes(y = endemicity*100, x = time, group = interaction(area, map), col = area),
            alpha = 0.01,
            linewidth = 3.0) +

  # Plot mean lines only
  geom_line(data = mean_areas_endemicity_all_time_df,
          mapping = aes(y = mean_endemicity*100, x = time, group = area, col = area),
          alpha = 1.0,
          linewidth = 2.0,
          show.legend = T) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of regional lineages") +
  
  # Set Y-axis to all be the same across plots
  ylim(c(0, 100)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     limits = c(85, 0) # Set limits
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
lines_endemicity_per_bioregions_all_time_steps_ggplot <- lines_endemicity_per_bioregions_all_time_steps_ggplot %>% 
  add_aesthetics_barplot()

print(lines_endemicity_per_bioregions_all_time_steps_ggplot)

dev.off()


## 5.3/ GGplot: with quantiles ####

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps.pdf", width = 10, height = 6)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Lineage_origin_percs_per_bioregions_all_time_steps_fuzzy_quantiles.pdf", width = 10, height = 6)

lines_endemicity_per_bioregions_all_time_steps_ggplot <- ggplot(data = areas_endemicity_all_time_all_maps_df) +
  
  # Plot 50 quantile polygons of observed data
  geom_polygon(data = areas_endemicity_all_time_all_maps_df_quantiles,
               mapping = aes(y = quant_percs*100, x = time, group = interaction(area, quantile_ID), fill = area),
               alpha = 0.02,
               # alpha = 0.03,
               linewidth = 2.0) +
  
  # Plot mean lines only
  geom_line(data = mean_areas_endemicity_all_time_df,
            mapping = aes(y = mean_endemicity*100, x = time, group = area, col = area),
            alpha = 1.0,
            linewidth = 2.0,
            show.legend = T) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of regional lineages") +
  
  # Set Y-axis to all be the same across plots
  ylim(c(0, 100)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", labels = bioregion_names_with_total, values = unname(colors_list_for_areas_with_total)) +
  # Adjust fill color scheme and legend
  scale_fill_manual("Bioregions", labels = bioregion_names_with_total, values = unname(colors_list_for_areas_with_total)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust margins
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        axis.text.x = element_text(size = 14))
# color = rev(time_strata_col),
# angle = 45, margin = (margin(t = 15)))) 

# Adjust aesthetics
lines_endemicity_per_bioregions_all_time_steps_ggplot <- lines_endemicity_per_bioregions_all_time_steps_ggplot %>% 
  add_aesthetics_barplot()

print(lines_endemicity_per_bioregions_all_time_steps_ggplot)

dev.off()



##### 6/ Map mean lineage immigration age #####

## For each tip, retrieve the timing at which its ancestor emigrated from another bioregions on its tip-to-root path
# Use the BSM simulations to compute a mean value over the 1000 simulations
# Consider a shift only is the tip area is not found in the previous range
# For taxa with current multi-area ranges, consider the immigration age for all encompassed areas, and compute the mean.

### 6.1/ Function to retrieve immigration rates of tips from BSM simmaps ####

get_immigration_ages_from_Simmap <- function (simmap,
                                              verbose = F)
{
  # Extract number of tips/taxa
  nb_taxa <- length(simmap$tip.label)
  
  # Initiate final vector
  immigration_ages_per_taxa <- c()
  
  # Get node IDs of root to tip paths for all tips
  root_to_tip_node_paths <- ape::nodepath(phy = simmap)
  # Get edge IDs of root to tip paths for all tips
  root_to_tip_edge_paths <- lapply(X = root_to_tip_node_paths, FUN = function (x) { which(simmap$edge[, 2] %in% x) } )
  # Get edge IDs of tip  to root paths for all tips
  tip_to_root_edge_paths <- lapply(X = root_to_tip_edge_paths, FUN = rev)
  
  # Get root age
  root_age <- max(phytools::nodeHeights(tree = simmap))
  
  # Extract edge ages
  edges_data_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = simmap) + root_age, 5))
  names(edges_data_df) <- c("rootward_age", "tipward_age")
  edges_data_df$edge_ID <- 1:nrow(edges_data_df)
  edges_data_df <- edges_data_df[, c("edge_ID", "rootward_age", "tipward_age")]
  
  ## Loop per tip
  for (i in 1:nb_taxa)
  {
    # i <- 1
    # i <- 70
    
    # Extract edge ID of tip-to-root path
    tip_to_root_edge_path_i <- tip_to_root_edge_paths[[i]]

    # Extract tip edge
    tip_edge_ID <- tip_to_root_edge_path_i[1]
    
    # Identify current range/area(s)
    current_range_i <- names(simmap$maps[[tip_edge_ID]])
    current_areas_i <- unlist(str_split(string = current_range_i, pattern = ""))
    
    ## Loop per current area(s)
    immigration_ages_per_current_areas_i <- c()
    for (j in seq_along(current_areas_i))
    {
      # j <- 1
      
      # Extract the focal area
      focal_area_j <- current_areas_i[j]
      
      # Move backward on tip-to-root path until you find a range that does not encompass the focal area
      shift_detected <- F
      edge_explored <- 1
      nb_edges <- length(tip_to_root_edge_path_i)
      while (!shift_detected & (edge_explored <= nb_edges))
      {
        # Extract focal edge
        focal_edge_ID <- tip_to_root_edge_path_i[edge_explored]
        # Extract mapping of focal edge
        focal_edge_maps <- simmap$maps[[focal_edge_ID]]
        # Reverse mapping from tip to root
        focal_edge_maps <- rev(focal_edge_maps)
        
        # Check for a shift
        focal_edge_ranges <- names(focal_edge_maps)
        focal_edge_areas <- str_split(string = focal_edge_ranges, pattern = "")
        shift_check <- unlist(lapply(X = focal_edge_areas, FUN = function (x) { !(focal_area_j %in% x) } ))
        shift_detected <- any(shift_check)
        
        # Increment nb of edges explored
        edge_explored <- edge_explored + 1
      }
      
      ## Case with no shift detected
      if (!shift_detected)
      {
        # Record root age as the immigration age
        immigration_ages_per_current_areas_i <- c(immigration_ages_per_current_areas_i, root_age)
        
      } else { ## Case with shift detected
        
        # Find shifting state on the focal edge maps
        shift_state_ID <- which.max(shift_check)
        
        # Find shift timing (from tipward node) on the focal edge maps
        focal_edge_maps_cumtime_from_tips <- cumsum(focal_edge_maps)
        if (shift_state_ID == 1)
        {
          shift_timing <- 0
        } else {
          shift_timing <- focal_edge_maps_cumtime_from_tips[shift_state_ID - 1]
        }
        
        # Compute immigration age as the sum of the tipward node age and the shift timing
        immigration_age_j <- edges_data_df$tipward_age[focal_edge_ID] + shift_timing
        
        # Record immigration age for the focal area
        immigration_ages_per_current_areas_i <- c(immigration_ages_per_current_areas_i, immigration_age_j)
      }
    }
    
    names(immigration_ages_per_current_areas_i) <- current_areas_i
    
    # Get mean immigration from all current area(s)
    immigration_age_i <- mean(immigration_ages_per_current_areas_i)
    
    # Record immigration age for the tip
    immigration_ages_per_taxa <- c(immigration_ages_per_taxa, immigration_age_i)
    names(immigration_ages_per_taxa) <- simmap$tip.label[1:i]
    
    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Immigration age computed for taxa n°", i, "/", length(simmap$tip.label),"\n"))
    }
  }
  
  # Export output = Named vector of immigration age per taxa
  return(immigration_ages_per_taxa)
}

get_immigration_ages_from_MultiSimmap <- function (multiSimmap,
                                                   verbose = T)
{
  # Initiate final df
  immigration_ages_df <- data.frame()
  
  ## Loop per simmap
  for (i in seq_along(multiSimmap))
  {
    # i <- 1
    
    # Extract focal simmap (BSM simulation)
    simmap_i <- multiSimmap[[i]]
    
    # Get immigration ages per taxa
    immigration_ages_i <- get_immigration_ages_from_Simmap(simmap = simmap_i, verbose = F)
    
    
    # Store in final df
    immigration_ages_df <- rbind(immigration_ages_df, immigration_ages_i)
    names(immigration_ages_df) <- simmap_i$tip.label
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Immigration ages computed for simmap n°", i, "/", length(multiSimmap),"\n"))
    }
  }
  
  # Output = immigration age df: BSM simulations x taxa
  return(immigration_ages_df)
}

### 6.2/ Compute immigration ages along BSM simulations ####

## Compute ages per BSM simulations
taxa_immigration_ages_per_BSM_simulations_df <- get_immigration_ages_from_MultiSimmap(DEC_J_simmaps, verbose = T)

# Save immigration ages per BSM simulations
# saveRDS(object = taxa_immigration_ages_per_BSM_simulations_df, "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")
saveRDS(object = taxa_immigration_ages_per_BSM_simulations_df, "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")

## Compute summary statistics across all BSM simulations

taxa_immigration_ages_summary_stats_df <- data.frame(taxa = names(taxa_immigration_ages_per_BSM_simulations_df))
taxa_immigration_ages_summary_stats_df$mean <- apply(X = taxa_immigration_ages_per_BSM_simulations_df, MARGIN = 2, FUN = mean)  
taxa_immigration_ages_summary_stats_df$median <- apply(X = taxa_immigration_ages_per_BSM_simulations_df, MARGIN = 2, FUN = median)  
taxa_immigration_ages_summary_stats_df$sd <- apply(X = taxa_immigration_ages_per_BSM_simulations_df, MARGIN = 2, FUN = sd)

CI95 <- apply(X = taxa_immigration_ages_per_BSM_simulations_df, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
taxa_immigration_ages_summary_stats_df$CI_2.5 <- CI95[1, ]
taxa_immigration_ages_summary_stats_df$CI_97.5 <- CI95[2, ]

HPD95 <- apply(X = taxa_immigration_ages_per_BSM_simulations_df, MARGIN = 2, FUN = BayesTwin::HPD, cred_int = 0.95)
taxa_immigration_ages_summary_stats_df$HPD_2.5 <- HPD95[1, ]
taxa_immigration_ages_summary_stats_df$HPD_97.5 <- HPD95[2, ]

# Favor the median. Less affected by boundary effects
plot(taxa_immigration_ages_summary_stats_df$mean ~ taxa_immigration_ages_summary_stats_df$median)
plot(taxa_immigration_ages_summary_stats_df$median ~ taxa_immigration_ages_summary_stats_df$mean)
plot(taxa_immigration_ages_summary_stats_df$CI_2.5 ~ taxa_immigration_ages_summary_stats_df$HPD_2.5)
plot(taxa_immigration_ages_summary_stats_df$CI_97.5 ~ taxa_immigration_ages_summary_stats_df$HPD_97.5)

# Save immigration ages per taxa df
# saveRDS(object = taxa_immigration_ages_summary_stats_df, "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_summary_stats_df.rds")
saveRDS(object = taxa_immigration_ages_summary_stats_df, "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_summary_stats_df.rds")


### 6.3/ Compute mean immigration ages in space ####

## Replace range data with immigration age values in the stack

# Load alpha-hull range data
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

# Load terrestrial background
terrestrial_bg_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")

# Load immigration ages per taxa df
# taxa_immigration_ages_summary_stats_df <- readRDS("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_summary_stats_df.rds")
taxa_immigration_ages_summary_stats_df <- readRDS("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_summary_stats_df.rds")

# Reorder as in raster stack
taxa_immigration_ages_summary_stats_df_reordered <- taxa_immigration_ages_summary_stats_df[match(names(Ponerinae_alpha_hull_stack_WGS84), table = taxa_immigration_ages_summary_stats_df$taxa), ]
cbind(taxa_immigration_ages_summary_stats_df_reordered$taxa, names(Ponerinae_alpha_hull_stack_WGS84))

length(taxa_immigration_ages_summary_stats_df$taxa)
length(taxa_immigration_ages_summary_stats_df_reordered$taxa)
length(names(Ponerinae_alpha_hull_stack_WGS84))

# Assign median immigration ages to binary ranges 
Ponerinae_immigration_ages_stack_WGS84 <- Ponerinae_alpha_hull_stack_WGS84 * round(taxa_immigration_ages_summary_stats_df_reordered$median, 1)
plot(Ponerinae_immigration_ages_stack_WGS84[[1:9]])
rm(Ponerinae_alpha_hull_stack_WGS84) ; gc() # Remove alpha hull ranges to save RAM space
Ponerinae_immigration_ages_stack_WGS84 <- readAll(Ponerinae_immigration_ages_stack_WGS84)

## Save stack of median immigration ages
# saveRDS(Ponerinae_immigration_ages_stack_WGS84, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Ponerinae_immigration_ages_stack_WGS84.rds")
saveRDS(Ponerinae_immigration_ages_stack_WGS84, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_immigration_ages_stack_WGS84.rds")

## Compute mean immigration ages raster (do not account for zeros as they are absences)
Ponerinae_mean_immigration_ages_WGS84 <- raster::calc(x = Ponerinae_immigration_ages_stack_WGS84,
                                                   fun = function (x) { mean(x[x > 0], na.rm = T) })
Ponerinae_mean_immigration_ages_WGS84 <- readAll(Ponerinae_mean_immigration_ages_WGS84)

# Add terrestrial boundaries
temp <- terrestrial_bg_WGS84
temp[!is.na(Ponerinae_mean_immigration_ages_WGS84@data@values)] <- Ponerinae_mean_immigration_ages_WGS84@data@values[!is.na(Ponerinae_mean_immigration_ages_WGS84@data@values)]
Ponerinae_mean_immigration_ages_WGS84 <- temp

plot(Ponerinae_mean_immigration_ages_WGS84)
rm(Ponerinae_immigration_ages_stack_WGS84) ; gc() # Remove stack of median immigration ages to save RAM space

## Save raster of mean immigration ages
# saveRDS(Ponerinae_mean_immigration_ages_WGS84, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_immigration_ages_WGS84.rds")
saveRDS(Ponerinae_mean_immigration_ages_WGS84, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_immigration_ages_WGS84.rds")


### 6.4/ Filter out pixels with low richness ####

## Load raster of mean current net div rates
# Ponerinae_mean_net_div_rates_WGS84 <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_immigration_ages_WGS84.rds")
Ponerinae_mean_immigration_ages_WGS84 <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_immigration_ages_WGS84.rds")

## Load raw richness raster
Ponerinae_species_richness_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_species_richness_WGS84.rds")

plot(Ponerinae_species_richness_WGS84, col = pal_bl_red_brewer)
plot(Ponerinae_mean_immigration_ages_WGS84, col = pal_purple_green_brewer)

# Set minimum richness to use to compute/display mean rates
min_richness_threshold <- 3 

# source("./functions/contrasting_raster.R")
# Ponerinae_species_richness_WGS84_thresholded <- contrasting_raster(x = Ponerinae_species_richness_WGS84,
#                                                                    zmin = min_richness_threshold,
#                                                                    zmax = max(Ponerinae_species_richness_WGS84@data@values, na.rm = T))
# 
# pdf(paste0("./outputs/Species_richness_maps/Ponerinae_richness_maps_multiple_thresholds.pdf"),
#     width = 20, height = 10)
# 
# par(mfrow = c(2,2))
# plot(Ponerinae_species_richness_WGS84, col = pal_bl_red_brewer)
# # Repeat with different thresholds
# plot(Ponerinae_species_richness_WGS84_thresholded, col = pal_bl_red_brewer, main = paste0("Min theshold = ", min_richness_threshold))
# par(mfrow = c(1,1))
# 
# dev.off()

source("./functions/contrasting_raster.R")
Ponerinae_mean_immigration_ages_WGS84_thresholded <- Ponerinae_mean_immigration_ages_WGS84
Ponerinae_mean_immigration_ages_WGS84_thresholded@data@values[Ponerinae_species_richness_WGS84@data@values < min_richness_threshold] <- 0

plot(Ponerinae_mean_immigration_ages_WGS84, col = pal_purple_green_brewer)
plot(Ponerinae_mean_immigration_ages_WGS84_thresholded, col = pal_purple_green_brewer)

pdf(file = paste0("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_immigration_ages_maps_multiple_thresholds.pdf"),
    width = 20, height = 10)

par(mfrow = c(2,2))
plot(Ponerinae_mean_immigration_ages_WGS84, col = pal_purple_green_brewer)
# Repeat with different thresholds
for (i in 1:3)
{
  min_richness_threshold <- i 
  Ponerinae_mean_immigration_ages_WGS84_thresholded <- Ponerinae_mean_immigration_ages_WGS84
  Ponerinae_mean_immigration_ages_WGS84_thresholded@data@values[Ponerinae_species_richness_WGS84@data@values < min_richness_threshold] <- 0
  plot(Ponerinae_mean_immigration_ages_WGS84_thresholded, col = pal_purple_green_brewer, main = paste0("Min theshold = ", min_richness_threshold))
}
par(mfrow = c(1,1))

dev.off()



### 6.5/ Plot mean immigration ages geographic map ####

## Load bioregion sf
Bioregions_sf_Bioregions_level_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregions_sf_Bioregions_level_Mollweide.rds")

# Load color palette
pal_bl_red_brewer <- rev(RColorBrewer::brewer.pal(n = 11, "PRGn"))
pal_bl_red_brewer_fn <- colorRampPalette(colors = pal_bl_red_brewer)
pal_bl_red_brewer <- pal_bl_red_brewer_fn(n = 400)
blues <- round(seq(from = 1, to = 260, length.out = 99), 0)
reds <- round(seq(from = 261, to = 400, length.out = 100), 0)
pal_bl_red_brewer <- pal_bl_red_brewer[c(blues, reds)]
pal_bl_red_brewer <- c("grey90", pal_bl_red_brewer)

# Load color palette
pal_purple_green_brewer <- rev(RColorBrewer::brewer.pal(n = 11, "PRGn"))
pal_purple_green_brewer_fn <- colorRampPalette(colors = pal_purple_green_brewer)
pal_purple_green_brewer <- pal_purple_green_brewer_fn(n = 400)
pal_purple_green_brewer <- c("grey90", pal_purple_green_brewer)

# Contrast raster
source("./functions/contrasting_raster.R")
hist(Ponerinae_mean_immigration_ages_WGS84[])
table(Ponerinae_mean_immigration_ages_WGS84[])
hist(Ponerinae_mean_immigration_ages_WGS84_thresholded[])
table(Ponerinae_mean_immigration_ages_WGS84_thresholded[])
Ponerinae_mean_immigration_ages_WGS84_contrasted <- contrasting_raster(# x = Ponerinae_mean_immigration_ages_WGS84, # Use the non-thresholded full version
                                                                       x = Ponerinae_mean_immigration_ages_WGS84_thresholded, # Use the thresholded version 
                                                                       zmin = 0, zmax = 105)

# Convert raster of mean current net div rates to Mollweide
Ponerinae_mean_immigration_ages_Mollweide_contrasted <- raster::projectRaster(from = Ponerinae_mean_immigration_ages_WGS84_contrasted,
                                                                           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
                                                                           method = "ngb")
plot(Ponerinae_mean_immigration_ages_Mollweide_contrasted, col = pal_purple_green_brewer)

# Convert raster to SpatialPixelsDataFrame
Ponerinae_mean_immigration_ages_Mollweide_spdf <- as(Ponerinae_mean_immigration_ages_Mollweide_contrasted, "SpatialPixelsDataFrame")
Ponerinae_mean_immigration_ages_Mollweide_spdf <- as.data.frame(Ponerinae_mean_immigration_ages_Mollweide_spdf)
colnames(Ponerinae_mean_immigration_ages_Mollweide_spdf) <- c("value", "x", "y")

# GGplot
Ponerinae_mean_immigration_ages_ggplot <- ggplot() +
  
  # Plot species richness as raster background
  geom_tile(data = Ponerinae_mean_immigration_ages_Mollweide_spdf,
            aes(x = x, y = y, fill = value), alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Immigration age  [My]",
                       # colors = pal_bl_red_Mannion,
                       # colors = pal_bl_red_brewer) +
                       colors = pal_purple_green_brewer) +

  # Plot bioregion sf maps
  geom_sf(data = Bioregions_sf_Bioregions_level_Mollweide,
          fill = NA,
          colour = "black",
          alpha = 0.0) +
  
  # Adjust CRS
  # coord_sf(default_crs = sf::st_crs(4326)) +
  coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           clip = "off", # To allow plotting arrow outside of World map
           expand = FALSE) +
  
  
  # Add title
  ggtitle(label =  paste0("Mean immigration ages of current local taxa")) +
  
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
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.margin = margin(l = 5),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())


# ## No threshold
# # pdf(file = paste0("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_immigration_ages_map_no_threshold.pdf"),
# pdf(file = paste0("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_immigration_ages_map_no_threshold.pdf"),
#     width = 10, height = 5)
# 
# print(Ponerinae_mean_immigration_ages_ggplot)
# 
# dev.off()

## With threshold
pdf(file = paste0("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_immigration_ages_map_threshold_",min_richness_threshold,".pdf"),
    width = 10, height = 5)

print(Ponerinae_mean_immigration_ages_ggplot)

dev.off()


##### 7/ Test for a longitudinal gradient of immigration age #####

longitude_scale <- seq(from = -180, to = 180, by = 1)

### 7.1/ Compute the mean immigration age along longitudinal bands ####

# Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

# Load immigration ages per BSM simulations
# taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")
taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")

# Load longitudinal bands df
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df.rds")

## Loop per longitudinal bands
immigration_ages_across_BSM_all_longitudes_df <- data.frame()
for (i in 1:nrow(longitudinal_bands_df))
{
  # i <- 1
  
  # Extract list of taxa present in band
  longitude_binary_presence_i <- unlist(longitude_binary_presence_df[i,, drop = T])
  taxa_in_band_i <- names(longitude_binary_presence_df)[longitude_binary_presence_i]
  taxa_in_band_i_ID <- match(x = taxa_in_band_i, table = names(taxa_immigration_ages_per_BSM_simulations_df))
  taxa_in_band_i_ID <- taxa_in_band_i_ID[!is.na(taxa_in_band_i_ID)]
  
  # Extract immigration ages for focal taxa along BSM simulations
  taxa_immigration_ages_df_in_band_i <- taxa_immigration_ages_per_BSM_simulations_df[, taxa_in_band_i_ID, drop = F]
  
  # Compute mean across taxa
  immigration_ages_across_BSM_i <- apply(X = taxa_immigration_ages_df_in_band_i, MARGIN = 1, FUN = mean, na.rm = T)
  immigration_ages_across_BSM_i[is.nan(immigration_ages_across_BSM_i)] <- NA
  
  # Store result
  immigration_ages_across_BSM_all_longitudes_df <- rbind(immigration_ages_across_BSM_all_longitudes_df, immigration_ages_across_BSM_i)
}
names(immigration_ages_across_BSM_all_longitudes_df) <- paste0("BSM_", 1:ncol(immigration_ages_across_BSM_all_longitudes_df))
row.names(immigration_ages_across_BSM_all_longitudes_df) <- longitudinal_bands_df$longitude_dec

View(immigration_ages_across_BSM_all_longitudes_df)

# Save df of immigration ages along BSM for all longitudinal bands 
# saveRDS(object = immigration_ages_across_BSM_all_longitudes_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_df.rds")
saveRDS(object = immigration_ages_across_BSM_all_longitudes_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_df.rds")

# Melt output to use in ggplot
immigration_ages_across_BSM_all_longitudes_df$longitude <- row.names(immigration_ages_across_BSM_all_longitudes_df)
immigration_ages_across_BSM_all_longitudes_melted_df <- 
  reshape2::melt(immigration_ages_across_BSM_all_longitudes_df) 
names(immigration_ages_across_BSM_all_longitudes_melted_df) <- c("longitude_dec", "BSM_simul", "immigration_age")
immigration_ages_across_BSM_all_longitudes_melted_df <- immigration_ages_across_BSM_all_longitudes_melted_df %>% 
  mutate(longitude_dec = as.numeric(as.character(longitude_dec)))

# Save melted df of immigration ages along BSM for all longitudinal bands 
# saveRDS(object = immigration_ages_across_BSM_all_longitudes_melted_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_melted_df.rds")
saveRDS(object = immigration_ages_across_BSM_all_longitudes_melted_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_melted_df.rds")

## Compute median immigration ages across BSM
immigration_ages_all_longitudes_median_df <- immigration_ages_across_BSM_all_longitudes_melted_df %>% 
  group_by(longitude_dec) %>% 
  summarise(median_immigration_age = median(immigration_age)) %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  arrange(longitude_dec)

# Inform longitudinal bands df
longitudinal_bands_df <- left_join(longitudinal_bands_df, immigration_ages_all_longitudes_median_df)

# Save longitudinal bands df
# saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")


### 7.2/ Compute null distribution of median immigration ages ####

# Make a null model with species identity randomize (i.e., randomize immigration ages across tips)
# Compute null distribution of mean immigration age along longitudinal bands

# Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

# Load immigration ages per BSM simulations
# taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")
taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")

# Load longitudinal bands df
# longitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
longitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Set nb of permutations
nb_perm <- 1000
nb_taxa <- dim(taxa_immigration_ages_per_BSM_simulations_df)[2]

## Loop per permutations
set.seed(seed = 1234)
immigration_ages_median_all_longitudes_perm_df <- data.frame()
for (i in 1:nb_perm)
{
  # i <- 1
  
  # Randomize taxa identity
  taxa_perm_ID <- sample(x = 1:nb_taxa, size = nb_taxa, replace = F)
  
  # Randomize immigration ages in BSM simulations
  taxa_immigration_ages_per_BSM_simulations_perm_i <- taxa_immigration_ages_per_BSM_simulations_df[, taxa_perm_ID]
  names(taxa_immigration_ages_per_BSM_simulations_perm_i) <- names(taxa_immigration_ages_per_BSM_simulations_df)
  
  ## Loop per longitudinal bands
  immigration_ages_across_BSM_all_longitudes_df_i <- data.frame()
  for (j in 1:nrow(longitudinal_bands_df))
  {
    # j <- 1
    
    # Extract list of taxa present in band
    longitude_binary_presence_j <- unlist(longitude_binary_presence_df[j,, drop = T])
    taxa_in_band_j <- names(longitude_binary_presence_df)[longitude_binary_presence_j]
    taxa_in_band_j_ID <- match(x = taxa_in_band_j, table = names(taxa_immigration_ages_per_BSM_simulations_perm_i))
    taxa_in_band_j_ID <- taxa_in_band_j_ID[!is.na(taxa_in_band_j_ID)]
    
    # Extract immigration ages for focal taxa along BSM simulations
    taxa_immigration_ages_df_in_band_j <- taxa_immigration_ages_per_BSM_simulations_perm_i[, taxa_in_band_j_ID, drop = F]
    
    # Compute mean across taxa
    immigration_ages_across_BSM_j <- apply(X = taxa_immigration_ages_df_in_band_j, MARGIN = 1, FUN = mean, na.rm = T)
    immigration_ages_across_BSM_j[is.nan(immigration_ages_across_BSM_j)] <- NA
    
    # Store result
    immigration_ages_across_BSM_all_longitudes_df_i <- rbind(immigration_ages_across_BSM_all_longitudes_df_i, immigration_ages_across_BSM_j)
  }
  names(immigration_ages_across_BSM_all_longitudes_df_i) <- paste0("BSM_", 1:ncol(immigration_ages_across_BSM_all_longitudes_df_i))
  row.names(immigration_ages_across_BSM_all_longitudes_df_i) <- longitudinal_bands_df$longitude_dec
  
  # Melt output to be able to compute median ages
  immigration_ages_across_BSM_all_longitudes_df_i$longitude <- row.names(immigration_ages_across_BSM_all_longitudes_df_i)
  immigration_ages_across_BSM_all_longitudes_melted_df_i <- 
    reshape2::melt(immigration_ages_across_BSM_all_longitudes_df_i)
  names(immigration_ages_across_BSM_all_longitudes_melted_df_i) <- c("longitude_dec", "BSM_simul", "immigration_age")
  
  # Compute median across all BSM for perm i
  
  immigration_ages_median_all_longitudes_df_i <- immigration_ages_across_BSM_all_longitudes_melted_df_i %>% 
    group_by(longitude_dec) %>% 
    summarise(median_immigration_age = median(immigration_age)) %>% 
    mutate(longitude_dec = as.numeric(longitude_dec)) %>%
    arrange(longitude_dec)
  
  # Store results from perm i
  immigration_ages_median_all_longitudes_perm_df <- rbind(immigration_ages_median_all_longitudes_perm_df, immigration_ages_median_all_longitudes_df_i$median_immigration_age)
  names(immigration_ages_median_all_longitudes_perm_df) <- immigration_ages_median_all_longitudes_df_i$longitude_dec
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Median immigration ages along Longitudinal bands computed for permutation n°", i, "/", nb_perm,"\n"))
    
    # Save null data for median immigration ages along Longitudinal bands
    # saveRDS(object = immigration_ages_median_all_longitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_df.rds")
    saveRDS(object = immigration_ages_median_all_longitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_df.rds")
  }
}

# Save null data for median immigration ages along Longitudinal bands
# saveRDS(object = immigration_ages_median_all_longitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_df.rds")
saveRDS(object = immigration_ages_median_all_longitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_df.rds")

# Melt output to use in ggplot
immigration_ages_median_all_longitudes_perm_df$perm_ID <- row.names(immigration_ages_median_all_longitudes_perm_df)
immigration_ages_median_all_longitudes_perm_melted_df <- 
  reshape2::melt(immigration_ages_median_all_longitudes_perm_df)
names(immigration_ages_median_all_longitudes_perm_melted_df) <- c("perm_ID", "longitude_dec", "immigration_age")
immigration_ages_median_all_longitudes_perm_melted_df <- immigration_ages_median_all_longitudes_perm_melted_df %>% 
  mutate(longitude_dec = as.numeric(as.character(longitude_dec)))

# Save melted df of null data for immigration ages along BSM for all longitudinal bands 
# saveRDS(object = immigration_ages_median_all_longitudes_perm_melted_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_melted_df.rds")
saveRDS(object = immigration_ages_median_all_longitudes_perm_melted_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_melted_df.rds")

## Compute median immigration ages across Permutations
immigration_ages_median_all_longitudes_perm_median_df <- immigration_ages_median_all_longitudes_perm_melted_df %>% 
  group_by(longitude_dec) %>% 
  summarise(median_immigration_age = median(immigration_age)) %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  arrange(longitude_dec)

# Save median values of null data for immigration ages for all longitudinal bands 
# saveRDS(object = immigration_ages_median_all_longitudes_perm_median_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_median_df.rds")
saveRDS(object = immigration_ages_median_all_longitudes_perm_median_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_median_df.rds")


### 7.3/ Plot immigration ages along longitudinal bands ####

# Load longitudinal bands df with SR data and median immigration ages
# longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Load melted df of immigration ages along BSM for all longitudinal bands 
# immigration_ages_across_BSM_all_longitudes_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_melted_df.rds")
immigration_ages_across_BSM_all_longitudes_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_longitudes_melted_df.rds")

# Load melted df of null data for immigration ages along BSM for all longitudinal bands 
# immigration_ages_median_all_longitudes_perm_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_melted_df.rds")
immigration_ages_median_all_longitudes_perm_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_melted_df.rds")

# Load median values of null data for immigration ages for all longitudinal bands 
# immigration_ages_median_all_longitudes_perm_median_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_median_df.rds")
immigration_ages_median_all_longitudes_perm_median_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_longitudes_perm_median_df.rds")

## 7.3.1/ Detect NA bands to cut down polygons in sections

NA_bands <- immigration_ages_median_all_longitudes_perm_melted_df %>% 
  group_by(longitude_dec) %>%
  summarize(NA_band = any(is.na(immigration_age)))

# Use RLE to identify consecutive bands
RLE_output <- rle(NA_bands$NA_band)
# Identify the bands starting a new polygon section
poly_starts <- cumsum(c(1, RLE_output$lengths[-length(RLE_output$lengths)]))
# Label polygon sections
poly_ID <- 1
NA_bands$poly_ID <- NA
for (i in seq_along(poly_starts))
{
  poly_start_i <- poly_starts[i]
  NA_bands$poly_ID[poly_start_i:nrow(NA_bands)] <- poly_ID
  poly_ID <- poly_ID + 1
}

# Remove NA bands and update poly_ID
NA_bands <- NA_bands %>% 
  filter(!NA_band) %>% 
  group_by(poly_ID) %>%
  mutate(poly_ID = ceiling(poly_ID / 2 )) %>%
  dplyr::select(longitude_dec, poly_ID)

# Assign poly_ID to observed and null data before designing polygons
immigration_ages_median_all_longitudes_perm_melted_df <- immigration_ages_median_all_longitudes_perm_melted_df %>% 
  left_join(y = NA_bands)
immigration_ages_across_BSM_all_longitudes_melted_df <- immigration_ages_across_BSM_all_longitudes_melted_df %>% 
  left_join(y = NA_bands)

## 7.3.2/ Compute quantiles of immigration ages used to draw polygons

immigration_ages_across_BSM_all_longitudes_melted_df_quantiles <- immigration_ages_across_BSM_all_longitudes_melted_df %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  group_by(longitude_dec, poly_ID) %>% 
  # Compute quantiles
  reframe(quant_ages = quantile(immigration_age, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
  group_by(longitude_dec) %>%
  mutate(quantile = seq(from = 0, to = 1, by = 0.01),
         quantile_ID = c(1:50, 51, 50:1)) %>%
  filter(quantile_ID != 51) %>% # Remove median
  # Remove NA bands
  filter(!is.na(poly_ID)) %>%
  # Assign points ID (order for drawing the polygon)
  group_by(quantile_ID, poly_ID) %>%
  arrange(quantile_ID, quantile, longitude_dec) %>%
  mutate(n_points = n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  # mutate(points_ID = c(1:361, 722:362)) %>%
  # Reorder by quantile_ID, poly_ID, points_ID to check conformity of polygons
  arrange(quantile_ID, poly_ID, points_ID) %>%
  ungroup()

# Adjust min/max values to the plot limits to avoid artifact in polygon drawings
max_ages <- 100
min_ages <- 0
immigration_ages_across_BSM_all_longitudes_melted_df_quantiles$quant_ages[immigration_ages_across_BSM_all_longitudes_melted_df_quantiles$quant_ages > max_ages] <- max_ages
immigration_ages_across_BSM_all_longitudes_melted_df_quantiles$quant_ages[immigration_ages_across_BSM_all_longitudes_melted_df_quantiles$quant_ages < min_ages] <- min_ages

# Remove most extreme quantiles to keep 95% CI (strictly a 94% CI)
immigration_ages_across_BSM_all_longitudes_melted_df_quantiles <- immigration_ages_across_BSM_all_longitudes_melted_df_quantiles %>% 
  filter(!(quantile %in% c(0, 0.01, 0.02, 0.98, 0.99, 1)))

## 7.3.3/ Compute CI95% quantiles of null values for immigration ages used to draw polygon
immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles <- immigration_ages_median_all_longitudes_perm_melted_df %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  group_by(longitude_dec, poly_ID) %>% 
  # Compute CI95% quantiles
  reframe(CI95_ages = quantile(immigration_age, probs = c(0.025, 0.975), na.rm = T)) %>%
  group_by(longitude_dec) %>%
  mutate(quantile = c(0.025, 0.975)) %>%
  # Remove NA bands
  filter(!is.na(poly_ID)) %>%
  group_by(poly_ID) %>%
  # Assign points ID (order for drawing the polygon)
  arrange(quantile, longitude_dec) %>%
  mutate(n_points = n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  # mutate(points_ID = c(1:361, 722:362)) %>%
  # Reorder by poly_ID, points_ID to check conformity of polygons
  arrange(poly_ID, points_ID) %>%
  ungroup()

# Adjust min/max values to the plot limits to avoid artifact in polygon drawings
max_ages <- 100
min_ages <- 0
immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_ages[immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_ages > max_ages] <- max_ages
immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_ages[immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_ages < min_ages] <- min_ages


## 7.3.4/ Create df for legend

immigration_ages_per_longitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                                y = c(50, 50),
                                                                data_type = c("Observed", "Null"))

## 7.3.5/ Generate plot using lines

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_per_longitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_per_longitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

immigration_ages_per_longitudinal_bands_ggplot_fuzzy_lines <- ggplot(data = immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_ages, x = longitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = immigration_ages_median_all_longitudes_perm_median_df,
            mapping = aes(y = median_immigration_age, x = longitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 1000 line replicates of observed data
  geom_line(data = immigration_ages_across_BSM_all_longitudes_melted_df,
            mapping = aes(y = immigration_age, x = longitude_dec, group = BSM_simul),
            col = "red",
            alpha = 0.01,
            linewidth = 3.0) +
  
  # Plot median lines of observed data
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = median_immigration_age, x = longitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = immigration_ages_per_longitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            linewidth = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                    values = c("red", "grey90")) +
  
  # Set Y-axis limits and reverse it to fit a vertical geological scale
  # scale_y_continuous(limits = c(0, 80)) +
  scale_y_continuous(transform = "reverse",
                     # limits = c(0, 80)) +
                     # limits = c(min_ages, max_ages)) +
                     limits = c(max_ages, min_ages)) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Longitude", breaks = c(-120, -60, 0, 60, 120), labels = c("120°W", "60°W", "0°", "60°E", "120°E")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage immigration ages across longitudes")) +
  
  # Set axes labels
  xlab("Longitude") +
  ylab("Immigration age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        # legend.position.inside = c(0.8, 0.85),
        legend.position.inside = c(0.8, 0.20),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(immigration_ages_per_longitudinal_bands_ggplot_fuzzy_lines)

dev.off()


## 7.3.6/ Generate plot using quantiles

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_per_longitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_per_longitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

immigration_ages_per_longitudinal_bands_ggplot_fuzzy_quantiles <- ggplot(data = immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = immigration_ages_median_all_longitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_ages, x = longitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = immigration_ages_median_all_longitudes_perm_median_df,
            mapping = aes(y = median_immigration_age, x = longitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 50 quantile polygons of observed data
  geom_polygon(data = immigration_ages_across_BSM_all_longitudes_melted_df_quantiles,
               mapping = aes(y = quant_ages, x = longitude_dec, group = interaction(quantile_ID, poly_ID)),
               fill = "red",
               # alpha = 0.03,
               alpha = 0.02,
               linewidth = 2.0) +
  
  # Plot median lines of observed data
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = median_immigration_age, x = longitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = immigration_ages_per_longitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            linewidth = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Set Y-axis limits and reverse it to fit a vertical geological scale
  # scale_y_continuous(limits = c(0, 80)) +
  scale_y_continuous(transform = "reverse",
                     # limits = c(0, 80)) +
                     # limits = c(min_ages, max_ages)) +
                     limits = c(max_ages, min_ages)) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Longitude", breaks = c(-120, -60, 0, 60, 120), labels = c("120°W", "60°W", "0°", "60°E", "120°E")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage immigration ages across longitudes")) +
  
  # Set axes labels
  xlab("Longitude") +
  ylab("Immigration age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        # legend.position.inside = c(0.8, 0.85),
        legend.position.inside = c(0.8, 0.20),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(immigration_ages_per_longitudinal_bands_ggplot_fuzzy_quantiles)

dev.off()


##### 8/ Test for a latitudinal gradient of immigration age #####

latitude_scale <- seq(from = -60, to = 60, by = 1)

### 8.1/ Compute the mean immigration age along latitudinal bands ####

# Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

# Load immigration ages per BSM simulations
# taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")
taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")

# Load latitudinal bands df
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df.rds")

## Loop per latitudinal bands
immigration_ages_across_BSM_all_latitudes_df <- data.frame()
for (i in 1:nrow(latitudinal_bands_df))
{
  # i <- 1
  
  # Extract list of taxa present in band
  latitude_binary_presence_i <- unlist(latitude_binary_presence_df[i,, drop = T])
  taxa_in_band_i <- names(latitude_binary_presence_df)[latitude_binary_presence_i]
  taxa_in_band_i_ID <- match(x = taxa_in_band_i, table = names(taxa_immigration_ages_per_BSM_simulations_df))
  taxa_in_band_i_ID <- taxa_in_band_i_ID[!is.na(taxa_in_band_i_ID)]
  
  # Extract immigration ages for focal taxa along BSM simulations
  taxa_immigration_ages_df_in_band_i <- taxa_immigration_ages_per_BSM_simulations_df[, taxa_in_band_i_ID, drop = F]
  
  # Compute mean across taxa
  immigration_ages_across_BSM_i <- apply(X = taxa_immigration_ages_df_in_band_i, MARGIN = 1, FUN = mean, na.rm = T)
  immigration_ages_across_BSM_i[is.nan(immigration_ages_across_BSM_i)] <- NA
  
  # Store result
  immigration_ages_across_BSM_all_latitudes_df <- rbind(immigration_ages_across_BSM_all_latitudes_df, immigration_ages_across_BSM_i)
}
names(immigration_ages_across_BSM_all_latitudes_df) <- paste0("BSM_", 1:ncol(immigration_ages_across_BSM_all_latitudes_df))
row.names(immigration_ages_across_BSM_all_latitudes_df) <- latitudinal_bands_df$latitude_dec

View(immigration_ages_across_BSM_all_latitudes_df)

# Save df of immigration ages along BSM for all latitudinal bands 
# saveRDS(object = immigration_ages_across_BSM_all_latitudes_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_df.rds")
saveRDS(object = immigration_ages_across_BSM_all_latitudes_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_df.rds")

# Melt output to use in ggplot
immigration_ages_across_BSM_all_latitudes_df$latitude <- row.names(immigration_ages_across_BSM_all_latitudes_df)
immigration_ages_across_BSM_all_latitudes_melted_df <- 
  reshape2::melt(immigration_ages_across_BSM_all_latitudes_df)
names(immigration_ages_across_BSM_all_latitudes_melted_df) <- c("latitude_dec", "BSM_simul", "immigration_age")
immigration_ages_across_BSM_all_latitudes_melted_df <- immigration_ages_across_BSM_all_latitudes_melted_df %>% 
  mutate(latitude_dec = as.numeric(as.character(latitude_dec)))

# Save melted df of immigration ages along BSM for all latitudinal bands 
# saveRDS(object = immigration_ages_across_BSM_all_latitudes_melted_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_melted_df.rds")
saveRDS(object = immigration_ages_across_BSM_all_latitudes_melted_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_melted_df.rds")

## Compute median immigration ages across BSM
immigration_ages_all_latitudes_median_df <- immigration_ages_across_BSM_all_latitudes_melted_df %>% 
  group_by(latitude_dec) %>% 
  summarise(median_immigration_age = median(immigration_age)) %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  arrange(latitude_dec)

# Inform latitudinal bands df
latitudinal_bands_df <- left_join(latitudinal_bands_df, immigration_ages_all_latitudes_median_df)

# Save updated latitudinal bands df
# saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")


### 8.2/ Compute null distribution of median immigration ages ####

# Make a null model with species identity randomize (i.e., randomize immigration ages across tips)
# Compute null distribution of mean immigration age along latitudinal bands

# Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

# Load immigration ages per BSM simulations
# taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")
taxa_immigration_ages_per_BSM_simulations_df <- readRDS("./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/taxa_immigration_ages_per_BSM_simulations_df.rds")

# Load latitudinal bands df
# latitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
latitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Set nb of permutations
nb_perm <- 1000
nb_taxa <- dim(taxa_immigration_ages_per_BSM_simulations_df)[2]

## Loop per permutations
set.seed(seed = 1234)
immigration_ages_median_all_latitudes_perm_df <- data.frame()
for (i in 1:nb_perm)
{
  # i <- 1
  
  # Randomize taxa identity
  taxa_perm_ID <- sample(x = 1:nb_taxa, size = nb_taxa, replace = F)
  
  # Randomize immigration ages in BSM simulations
  taxa_immigration_ages_per_BSM_simulations_perm_i <- taxa_immigration_ages_per_BSM_simulations_df[, taxa_perm_ID]
  names(taxa_immigration_ages_per_BSM_simulations_perm_i) <- names(taxa_immigration_ages_per_BSM_simulations_df)
  
  ## Loop per latitudinal bands
  immigration_ages_across_BSM_all_latitudes_df_i <- data.frame()
  for (j in 1:nrow(latitudinal_bands_df))
  {
    # j <- 1
    
    # Extract list of taxa present in band
    latitude_binary_presence_j <- unlist(latitude_binary_presence_df[j,, drop = T])
    taxa_in_band_j <- names(latitude_binary_presence_df)[latitude_binary_presence_j]
    taxa_in_band_j_ID <- match(x = taxa_in_band_j, table = names(taxa_immigration_ages_per_BSM_simulations_perm_i))
    taxa_in_band_j_ID <- taxa_in_band_j_ID[!is.na(taxa_in_band_j_ID)]
    
    # Extract immigration ages for focal taxa along BSM simulations
    taxa_immigration_ages_df_in_band_j <- taxa_immigration_ages_per_BSM_simulations_perm_i[, taxa_in_band_j_ID, drop = F]
    
    # Compute mean across taxa
    immigration_ages_across_BSM_j <- apply(X = taxa_immigration_ages_df_in_band_j, MARGIN = 1, FUN = mean, na.rm = T)
    immigration_ages_across_BSM_j[is.nan(immigration_ages_across_BSM_j)] <- NA
    
    # Store result
    immigration_ages_across_BSM_all_latitudes_df_i <- rbind(immigration_ages_across_BSM_all_latitudes_df_i, immigration_ages_across_BSM_j)
  }
  names(immigration_ages_across_BSM_all_latitudes_df_i) <- paste0("BSM_", 1:ncol(immigration_ages_across_BSM_all_latitudes_df_i))
  row.names(immigration_ages_across_BSM_all_latitudes_df_i) <- latitudinal_bands_df$latitude_dec
  
  # Melt output to be able to compute mean ages
  immigration_ages_across_BSM_all_latitudes_df_i$latitude <- row.names(immigration_ages_across_BSM_all_latitudes_df_i)
  immigration_ages_across_BSM_all_latitudes_melted_df_i <- 
    reshape2::melt(immigration_ages_across_BSM_all_latitudes_df_i)
  names(immigration_ages_across_BSM_all_latitudes_melted_df_i) <- c("latitude_dec", "BSM_simul", "immigration_age")
  
  # Compute median across all BSM for perm i
  
  immigration_ages_median_all_latitudes_df_i <- immigration_ages_across_BSM_all_latitudes_melted_df_i %>% 
    group_by(latitude_dec) %>% 
    summarise(median_immigration_age = median(immigration_age)) %>% 
    mutate(latitude_dec = as.numeric(latitude_dec)) %>%
    arrange(latitude_dec)
  
  # Store results from perm i
  immigration_ages_median_all_latitudes_perm_df <- rbind(immigration_ages_median_all_latitudes_perm_df, immigration_ages_median_all_latitudes_df_i$median_immigration_age)
  names(immigration_ages_median_all_latitudes_perm_df) <- immigration_ages_median_all_latitudes_df_i$latitude_dec
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Median immigration ages along latitudinal bands computed for permutation n°", i, "/", nb_perm,"\n"))
    
    # Save null data for median immigration ages along latitudinal bands
    # saveRDS(object = immigration_ages_median_all_latitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_df.rds")
    saveRDS(object = immigration_ages_median_all_latitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_df.rds")
  }
}

# Save null data for median immigration ages along latitudinal bands
# saveRDS(object = immigration_ages_median_all_latitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_df.rds")
saveRDS(object = immigration_ages_median_all_latitudes_perm_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_df.rds")

# Melt output to use in ggplot
immigration_ages_median_all_latitudes_perm_df$perm_ID <- row.names(immigration_ages_median_all_latitudes_perm_df)
immigration_ages_median_all_latitudes_perm_melted_df <- 
  reshape2::melt(immigration_ages_median_all_latitudes_perm_df)
names(immigration_ages_median_all_latitudes_perm_melted_df) <- c("perm_ID", "latitude_dec", "immigration_age")
immigration_ages_median_all_latitudes_perm_melted_df <- immigration_ages_median_all_latitudes_perm_melted_df %>% 
  mutate(latitude_dec = as.numeric(as.character(latitude_dec)))

# Save melted df of null data for immigration ages along BSM for all latitudinal bands 
# saveRDS(object = immigration_ages_median_all_latitudes_perm_melted_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_melted_df.rds")
saveRDS(object = immigration_ages_median_all_latitudes_perm_melted_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_melted_df.rds")

## Compute median immigration ages across Permutations
immigration_ages_median_all_latitudes_perm_median_df <- immigration_ages_median_all_latitudes_perm_melted_df %>% 
  group_by(latitude_dec) %>% 
  summarise(median_immigration_age = median(immigration_age)) %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  arrange(latitude_dec)

# Save median values of null data for immigration ages for all latitudinal bands 
# saveRDS(object = immigration_ages_median_all_latitudes_perm_median_df, file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_median_df.rds")
saveRDS(object = immigration_ages_median_all_latitudes_perm_median_df, file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_median_df.rds")


### 8.3/ Plot immigration ages along latitudinal bands ####

# Load latitudinal bands df with SR data and median immigration ages
# latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Load melted df of immigration ages along BSM for all latitudinal bands 
# immigration_ages_across_BSM_all_latitudes_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_melted_df.rds")
immigration_ages_across_BSM_all_latitudes_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_across_BSM_all_latitudes_melted_df.rds")

# Load melted df of null data for immigration ages along BSM for all latitudinal bands 
# immigration_ages_median_all_latitudes_perm_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_melted_df.rds")
immigration_ages_median_all_latitudes_perm_melted_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_melted_df.rds")

# Load median values of null data for immigration ages for all latitudinal bands 
# immigration_ages_median_all_latitudes_perm_median_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_median_df.rds")
immigration_ages_median_all_latitudes_perm_median_df <- readRDS(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_median_all_latitudes_perm_median_df.rds")

## 8.3.1/ Detect NA bands to cut down polygons in sections

NA_bands <- immigration_ages_median_all_latitudes_perm_melted_df %>% 
  group_by(latitude_dec) %>%
  summarize(NA_band = any(is.na(immigration_age)))

# Use RLE to identify consecutive bands
RLE_output <- rle(NA_bands$NA_band)
# Identify the bands starting a new polygon section
poly_starts <- cumsum(c(1, RLE_output$lengths[-length(RLE_output$lengths)]))
# Label polygon sections
poly_ID <- 1
NA_bands$poly_ID <- NA
for (i in seq_along(poly_starts))
{
  poly_start_i <- poly_starts[i]
  NA_bands$poly_ID[poly_start_i:nrow(NA_bands)] <- poly_ID
  poly_ID <- poly_ID + 1
}

# Remove NA bands and update poly_ID
NA_bands <- NA_bands %>% 
  filter(!NA_band) %>% 
  group_by(poly_ID) %>%
  mutate(poly_ID = ceiling(poly_ID / 2 )) %>%
  dplyr::select(latitude_dec, poly_ID)

# Assign poly_ID to observed and null data before designing polygons
immigration_ages_median_all_latitudes_perm_melted_df <- immigration_ages_median_all_latitudes_perm_melted_df %>% 
  left_join(y = NA_bands)
immigration_ages_across_BSM_all_latitudes_melted_df <- immigration_ages_across_BSM_all_latitudes_melted_df %>% 
  left_join(y = NA_bands)

## 8.3.2/ Compute quantiles of immigration ages used to draw polygons

immigration_ages_across_BSM_all_latitudes_melted_df_quantiles <- immigration_ages_across_BSM_all_latitudes_melted_df %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  group_by(latitude_dec, poly_ID) %>% 
  # Compute quantiles
  reframe(quant_ages = quantile(immigration_age, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
  group_by(latitude_dec) %>%
  mutate(quantile = seq(from = 0, to = 1, by = 0.01),
         quantile_ID = c(1:50, 51, 50:1)) %>%
  filter(quantile_ID != 51) %>% # Remove median
  # Remove NA bands
  filter(!is.na(poly_ID)) %>%
  # Assign points ID (order for drawing the polygon)
  group_by(quantile_ID, poly_ID) %>%
  arrange(quantile_ID, quantile, latitude_dec) %>%
  mutate(n_points = n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  # mutate(points_ID = c(1:361, 722:362)) %>%
  # Reorder by quantile_ID, poly_ID, points_ID to check conformity of polygons
  arrange(quantile_ID, poly_ID, points_ID) %>%
  ungroup()

# Adjust min/max values to the plot limits to avoid artifact in polygon drawings
max_ages <- 100
min_ages <- 0
immigration_ages_across_BSM_all_latitudes_melted_df_quantiles$quant_ages[immigration_ages_across_BSM_all_latitudes_melted_df_quantiles$quant_ages > max_ages] <- max_ages
immigration_ages_across_BSM_all_latitudes_melted_df_quantiles$quant_ages[immigration_ages_across_BSM_all_latitudes_melted_df_quantiles$quant_ages < min_ages] <- min_ages

# Remove most extreme quantiles to keep 95% CI (strictly a 94% CI)
immigration_ages_across_BSM_all_latitudes_melted_df_quantiles <- immigration_ages_across_BSM_all_latitudes_melted_df_quantiles %>% 
  filter(!(quantile %in% c(0, 0.01, 0.02, 0.98, 0.99, 1)))

## 8.3.3/ Compute CI95% quantiles of null values for immigration ages used to draw polygon

immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles <- immigration_ages_median_all_latitudes_perm_melted_df %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  group_by(latitude_dec, poly_ID) %>% 
  # Compute CI95% quantiles
  reframe(CI95_ages = quantile(immigration_age, probs = c(0.025, 0.975), na.rm = T)) %>%
  group_by(latitude_dec) %>%
  mutate(quantile = c(0.025, 0.975)) %>%
  # Remove NA bands
  filter(!is.na(poly_ID)) %>%
  group_by(poly_ID) %>%
  # Assign points ID (order for drawing the polygon)
  arrange(quantile, latitude_dec) %>%
  mutate(n_points = n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  # mutate(points_ID = c(1:361, 722:362)) %>%
  # Reorder by poly_ID, points_ID to check conformity of polygons
  arrange(poly_ID, points_ID) %>%
  ungroup()

# Adjust min/max values to the plot limits to avoid artifact in polygon drawings
max_ages <- 100
min_ages <- 0
immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_ages[immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_ages > max_ages] <- max_ages
immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_ages[immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_ages < min_ages] <- min_ages


## 8.3.4/ Create df for legend

immigration_ages_per_latitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                                y = c(50, 50),
                                                                data_type = c("Observed", "Null"))

## 8.3.5/ Generate plot using lines

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_per_latitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_per_latitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

immigration_ages_per_latitudinal_bands_ggplot_fuzzy_lines <- ggplot(data = immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_ages, x = latitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = immigration_ages_median_all_latitudes_perm_median_df,
            mapping = aes(y = median_immigration_age, x = latitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 1000 line replicates of observed data
  geom_line(data = immigration_ages_across_BSM_all_latitudes_melted_df,
            mapping = aes(y = immigration_age, x = latitude_dec, group = BSM_simul),
            col = "red",
            alpha = 0.01,
            linewidth = 3.0) +
  
  # Plot median lines of observed data
  geom_line(data = latitudinal_bands_df,
            mapping = aes(y = median_immigration_age, x = latitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Flip coordinates to have vertical latitude
  coord_flip() +
  
  # Reverse age scale (Y-axis since the flip)
  scale_y_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(max_ages, min_ages) # Set limits
  ) + 

  # Adjust label on Latitude axis
  # scale_x_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  scale_x_continuous("Latitude", limits = c(-60, 80), breaks = c(-40, 0, 40, 80), labels = c("40°S", "0°", "40°N", "80°N")) +
  
  # Add fake data for legend
  geom_line(data = immigration_ages_per_latitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            size = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage immigration ages\nacross latitudes")) +
  
  # Set axes labels
  xlab("Latitude") +
  ylab("Immigration age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        # legend.position.inside = c(0.85, 0.50),
        legend.position.inside = c(0.85, 0.40),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 5, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(immigration_ages_per_latitudinal_bands_ggplot_fuzzy_lines)

dev.off()


## 8.3.6/ Generate plot using quantiles

# pdf(file = "./outputs/Lineage_origins/Ponerinae_rough_phylogeny_1534t/immigration_ages_per_latitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/Lineage_origins/Ponerinae_MCC_phylogeny_1534t/immigration_ages_per_latitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

immigration_ages_per_latitudinal_bands_ggplot_fuzzy_quantiles <- ggplot(data = immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = immigration_ages_median_all_latitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_ages, x = latitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = immigration_ages_median_all_latitudes_perm_median_df,
            mapping = aes(y = median_immigration_age, x = latitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 50 quantile polygons of observed data
  geom_polygon(data = immigration_ages_across_BSM_all_latitudes_melted_df_quantiles,
               mapping = aes(y = quant_ages, x = latitude_dec, group = interaction(quantile_ID, poly_ID)),
               fill = "red",
               alpha = 0.03,
               linewidth = 2.0) +
  
  # Plot median lines of observed data
  geom_line(data = latitudinal_bands_df,
            mapping = aes(y = median_immigration_age, x = latitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Flip coordinates to have vertical latitude
  coord_flip() +
  
  # Reverse age scale (Y-axis since the flip)
  scale_y_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(max_ages, min_ages) # Set limits
  ) + 
  
  # Adjust label on Latitude axis
  # scale_x_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  scale_x_continuous("Latitude", limits = c(-60, 80), breaks = c(-40, 0, 40, 80), labels = c("40°S", "0°", "40°N", "80°N")) +
  
  # Add fake data for legend
  geom_line(data = immigration_ages_per_latitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            size = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage immigration ages\nacross latitudes")) +
  
  # Set axes labels
  xlab("Latitude") +
  ylab("Immigration age  [My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        # legend.position.inside = c(0.85, 0.50),
        legend.position.inside = c(0.85, 0.40),
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 5, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(immigration_ages_per_latitudinal_bands_ggplot_fuzzy_quantiles)

dev.off()


