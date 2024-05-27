##### Script 04: Fit models in BioGeoBEARS #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Fit different biogeographic models in BioGeoBEARS
# Compare model fits and select the best one
  # DEC, DEC + J, DIVALIKE, DIVALIKE + J
  # + X versions of the above (?)

###

### Inputs

# (Set of) Time-calibrated phylogeny/ies with imputed missing taxa
# Biogeographic ranges per species
# Hyperparameters for Biogeographic models:
  # Geological time boundaries
  # Time-stratified adjacency matrices (affects range space by canceling some combined areas from the range space)
  # Time-stratified dispersal multiplier matrices (for +X models) (modulates dispersal events (d and j))

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
library(ape)
library("BioGeoBEARS")

# devtools::install_github("nmatzke/BioGeoBEARS")

### 1.2/ Load biogeographic ranges per species ####

Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Temporary fixes to avoid having ranges across more than 3 bioregions
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Anochetus_graeffei", "Afrotropics"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Anochetus_graeffei", "Neotropics"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Anochetus_graeffei", "Eastern Palearctic"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Odontomachus_haematodus", "Indomalaya"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Odontomachus_haematodus", "Nearctic"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Odontomachus_haematodus", "Neotropics"] <- F
Taxa_bioregions_binary_table[Taxa_bioregions_binary_table$Current_name == "Parvaponera_darwinii", "Western Palearctic"] <- F

View(Taxa_bioregions_binary_table)

range_size <- apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 1, FUN = sum)
max(range_size)

Taxa_bioregions_binary_table$Current_name[range_size >= 3]

### 1.3/ Time-calibrated phylogeny(ies) ####

Ponerinae_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t.rds")

# Match names

Ponerinae_phylogeny_1534t_short_names <- Ponerinae_phylogeny_1534t
tips_short_names <- Ponerinae_phylogeny_1534t_short_names$tip.label

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

Ponerinae_phylogeny_1534t_short_names$tip.label <- tips_short_names

# Save imputed phylogeny with short names
saveRDS(object = Ponerinae_phylogeny_1534t_short_names, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")


##### 2/ Plot biogeographic tip data ####

# Extract binary matrix of presence/absence
binary_matrix <- as.matrix(Taxa_bioregions_binary_table[, -1])
# binary_matrix <- lapply(binary_matrix, factor) # Convert character strings into factors

dim(binary_matrix)

# Rename the columns of the matrix to correspond to our geographical areas
# colnames(binary_matrix) <- names(Taxa_bioregions_binary_table[, -1])
# Set the colors we’ll use for plotting
colors <- setNames(object = replicate(ncol(binary_matrix),
                                      setNames(c("white", "gray30"), 0:1),
                                      simplify = FALSE),
                   nm = colnames(binary_matrix))

pdf(file = "./outputs/Phylo_1534t_with_ranges.pdf", width = 20, height = 200)

# Set plotting parameters
par(mar = c(4.1,0.1,3.1,0.1))
# Graph presence/absence using plotTree.datamatrix
range_map <- phytools::plotTree.datamatrix(tree = Ponerinae_phylogeny_1534t_short_names,
                                           X = binary_matrix,
                                           fsize = 0.5, yexp = 1.1,
                                           header = TRUE, xexp = 1.25, colors = colors)
## Add a legend
legend(x = "topleft", legend = c("Absence", "Presence"),
       pch = 22, pt.bg = c("white","gray30"), pt.cex=  1.3,
       cex = 0.8, bty = "n")

dev.off()


# Kawahara et al., 2023

# Use seven time slices associated to geological periods: 
#   	0-5.33 My = Holocene + Pleistocene + Pliocene
#   	5.33-23.03 = Miocene 
#   	23.03-33.9 = Oligocene
#   	33.9-56.0 = Eocene
#   	56.0-66.0 = Paleocene
#   	66.0-100.5 = Late Cretaceous
#   	100.5-145.0 = Early Cretaceous (Root of Ponerinae = 124 My ?)
  
# Max areas per range = 3


