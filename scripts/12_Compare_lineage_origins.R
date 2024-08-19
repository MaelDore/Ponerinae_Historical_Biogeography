##### Script 12: Compare origin of lineages: in situ speciation vs. dispersal #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compare origin of lineages: in situ speciation vs. dispersal along time
# Plot evolution of proportion of origins in time, overall and per bioregions

###

### Inputs

# Simmaps of Biogeographic Stochastic Maps

###

### Outputs

### Plot evolution of proportion of origins in time
  # Current time = stacked barplot of Counts and Percentages per Bioregions
  # Along time = lines of % per bioregions, including total

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
library(phangorn)

 # devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load modeling results from best model (DEC+J) #####

# Load time-stratified DEC+J model output
DEC_J_fit <- readRDS(file = "./outputs/BioGeoBEARS_models/model_fits/DEC_J_fit.rds")

# Extract areas list
returned_mats <- BioGeoBEARS::get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = DEC_J_fit$inputs, include_null_range = DEC_J_fit$inputs$include_null_range)
areas_list <- returned_mats$areanames

# Load simmaps of BS maps
DEC_J_simmaps <- readRDS(file = "./outputs/BSM/DEC_J_simmaps.rds")

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
saveRDS(areas_endemicity_all_time_all_maps_array, file = "./outputs/Lineage_origins/areas_endemicity_all_time_all_maps_array.rds")

### 2.2/ Melt results for ggplot ####

areas_endemicity_all_time_all_maps_array <- readRDS(file = "./outputs/Lineage_origins/areas_endemicity_all_time_all_maps_array.rds")


areas_endemicity_all_time_all_maps_df <- reshape2::melt(areas_endemicity_all_time_all_maps_array)
names(areas_endemicity_all_time_all_maps_df) <- c("area", "time", "map", "endemicity")

# Save melted df of endemicity per bioregions x time x maps
saveRDS(areas_endemicity_all_time_all_maps_df, file = "./outputs/Lineage_origins/areas_endemicity_all_time_all_maps_df.rds")

### 2.3/ Aggregate along maps ####

mean_areas_endemicity_all_time_df <- areas_endemicity_all_time_all_maps_df %>% 
  group_by(area, time) %>% 
  summarize(mean_endemicity = mean(endemicity, na.rm = T)) %>% 
  ungroup()

# Replace Nan by NA
mean_areas_endemicity_all_time_df$mean_endemicity[is.nan(mean_areas_endemicity_all_time_df$mean_endemicity)] <- NA

# Save melted df of mean endemicity per bioregions x time
saveRDS(mean_areas_endemicity_all_time_df, file = "./outputs/Lineage_origins/mean_areas_endemicity_all_time_df.rds")


### 2.4/ Compute total counts and percentages in current time ####

# Load LTT data per bioregions along time
DEC_J_LTT_all_areas_mean_ggplot <- readRDS(file = "./outputs/LTT/DEC_J_LTT_all_areas_mean_ggplot.rds")

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
saveRDS(mean_areas_endemicity_current_time_df, file = "./outputs/Lineage_origins/mean_areas_endemicity_current_time_df.rds")


##### 3/ Plot counts of lineages per origin in current time ####

# Load melted df of mean counts and percentages of endemicity per bioregions in current time
mean_areas_endemicity_current_time_df <- readRDS(file = "./outputs/Lineage_origins/mean_areas_endemicity_current_time_df.rds")

mean_areas_endemicity_current_time_df_no_total <- mean_areas_endemicity_current_time_df %>% 
  filter(area != "total")

# Add % as labels
mean_areas_endemicity_current_time_df_no_total$perc_labels <- paste0(round(mean_areas_endemicity_current_time_df_no_total$perc, 0), "%")
mean_areas_endemicity_current_time_df_no_total$perc_labels_filtered <- mean_areas_endemicity_current_time_df_no_total$perc_labels
mean_areas_endemicity_current_time_df_no_total$perc_labels_filtered[mean_areas_endemicity_current_time_df_no_total$origin == "Dispersal"] <- NA

# Adjust color scheme for bioregions
areas_list <- c("A", "U", "I", "R", "N", "E", "W")
colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
colors_list_for_areas <- colors_list_for_states[areas_list]

# Adjust color scheme for borders (interaction between orgin and area)
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
pdf(file = "./outputs/Lineage_origins/Barplot_lineage_origin_counts_per_bioregions_current_time.pdf", width = 10, height = 6)

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
pdf(file = "./outputs/Lineage_origins/Barplot_lineage_origin_percs_per_bioregions_current_time.pdf", width = 10, height = 6)

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

### Add % as geom_label on the stack bars. See my scripts for CroCo outputs


##### 5/ Plot evolution of proportion of regional lineages along time #####

# Adjust Bioregions fill scale to add "Total"
colors_list_for_areas_with_total <- c(colors_list_for_areas, "grey80")
bioregion_names_with_total <- c(bioregion_names, "Total")

# Load melted df of endemicity per bioregions x time x maps
areas_endemicity_all_time_all_maps_df <- readRDS(file = "./outputs/Lineage_origins/areas_endemicity_all_time_all_maps_df.rds")
# Load melted df of mean endemicity per bioregions x time
mean_areas_endemicity_all_time_df <- readRDS(file = "./outputs/Lineage_origins/mean_areas_endemicity_all_time_df.rds")

# GGplot
pdf(file = "./outputs/Lineage_origins/Lineage_origin_percs_per_bioregions_all_time_steps.pdf", width = 10, height = 6)

lines_endemicity_per_bioregions_all_time_steps_ggplot <- ggplot(data = areas_endemicity_all_time_all_maps_df) +
  
  # # Plot mean lines + 1000 replicates
  # geom_smooth(data = areas_endemicity_all_time_all_maps_df,
  #             mapping = aes(y = endemicity, x = time, group = area, col = area, fill = area),
  #             method = "gam", se = TRUE, linewidth = 1.5, alpha = 0.1) +
  
  # # Plot 1000 replicates
  # geom_line(data = areas_endemicity_all_time_all_maps_df,
  #           mapping = aes(y = endemicity, x = time, group = interaction(area, map), col = area),
  #           alpha = 0.01,
  #           linewidth = 2.0) +

  # Plot mean lines only
  geom_line(data = mean_areas_endemicity_all_time_df,
            mapping = aes(y = mean_endemicity, x = time, group = area, col = area),
            alpha = 1.0,
            linewidth = 2.0,
            show.legend = T) +
  
  # Set plot title +
  ggtitle(label = paste0("Lineage origins per Bioregions\nacross 1000 BS maps")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("% of regional lineages") +
  
  # Set Y-axis to all be the same across bioregion types
  # ylim(c(0, max_source_dest_rates)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
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
lines_endemicity_per_bioregions_all_time_steps_ggplot <- lines_endemicity_per_bioregions_all_time_steps_ggplot %>% 
  add_aesthetics_barplot()

print(lines_endemicity_per_bioregions_all_time_steps_ggplot)

dev.off()



