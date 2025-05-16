##### Script 22: Compare topologies across partition schemes  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Plot and compare topologies of ML phylogenies obtained with the different partition schemes
# 3 partition schemes: unpartitioned; by-locus; SWSC-EN

###

### Inputs

# ML phylogenies obtained with the 3 different partition schemes

###

### Outputs

# Plot showing discrepancies between the topologies

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(phytools)
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(Quartet)   # To compare trees using quartets
library(tidyverse)

### 1.2/ Load non-calibrated ML backbone phylogenies ####

Ponerinae_nopart_792t <- read.tree(file = "./input_data/Phylogenies/ponerinae-792t-nopart-iqtree-result-names-updated.tre")
Ponerinae_bylocus_792t <- read.tree(file = "./input_data/Phylogenies/ponerinae-792t-bylocus-iqtree-result-names-updated.tre")
Ponerinae_SWSC_792t <- read.tree(file = "./input_data/Phylogenies/ponerinae-792t-swsc-iqtree-result-names-updated.tre")

# Fix Neoponera bucki
Ponerinae_nopart_792t$tip.label[Ponerinae_nopart_792t$tip.label == "NewGenus_bucki_EX2455_DZUP549431"] <- "Neoponera_bucki_EX2455_DZUP549431"
Ponerinae_bylocus_792t$tip.label[Ponerinae_bylocus_792t$tip.label == "NewGenus_bucki_EX2455_DZUP549431"] <- "Neoponera_bucki_EX2455_DZUP549431"
Ponerinae_SWSC_792t$tip.label[Ponerinae_SWSC_792t$tip.label == "NewGenus_bucki_EX2455_DZUP549431"] <- "Neoponera_bucki_EX2455_DZUP549431"

# Root trees using outliers

Ponerinae_nopart_792t_rooted <- ape::root(phy = Ponerinae_nopart_792t, outgroup = "Amblyopone_australis_D0872_CASENT0106229", resolve.root = TRUE)
Ponerinae_bylocus_792t_rooted <- ape::root(phy = Ponerinae_bylocus_792t, outgroup = "Amblyopone_australis_D0872_CASENT0106229", resolve.root = TRUE)
Ponerinae_SWSC_792t_rooted <- ape::root(phy = Ponerinae_SWSC_792t, outgroup = "Amblyopone_australis_D0872_CASENT0106229", resolve.root = TRUE)

Ponerinae_792t_trees_list <- list(SWSC = Ponerinae_SWSC_792t_rooted,
                                  bylocus = Ponerinae_bylocus_792t_rooted,
                                  nopart = Ponerinae_nopart_792t_rooted)

plot(Ponerinae_nopart_792t)

plot(Ponerinae_nopart_792t_rooted)
plot(Ponerinae_bylocus_792t_rooted)
plot(Ponerinae_SWSC_792t_rooted)

##### 2/ Identify differences with ape::comparePhylo #####

?ape::comparePhylo

test_nopart <- ape::comparePhylo(x = Ponerinae_SWSC_792t_rooted, y = Ponerinae_nopart_792t_rooted, plot = FALSE, use.edge.length = TRUE)
test_nopart$messages

View(test_nopart$NODES)

shared_clades_ID_nopart <- str_remove(string = test_nopart$NODES[, 2], pattern = ".* \\(")
shared_clades_ID_nopart <- as.numeric(str_remove(string = shared_clades_ID_nopart, pattern = "\\)"))
length(shared_clades_ID_nopart)

nb_tips <- length(Ponerinae_nopart_792t_rooted$tip.label)
nb_nodes <- Ponerinae_nopart_792t_rooted$Nnode
all_nodes_ID <- (nb_tips + 1):(nb_tips + nb_nodes)
unique_clades_ID_nopart <- all_nodes_ID[!(all_nodes_ID %in% shared_clades_ID_nopart)]

# 784 shared splits among 790 ! 99.2% identity => Splits are based on internal edges of unrooted trees
# 786 shared clades among 791 ! 99.3% identity => Clades are based on internal nodes of rooted trees with resolved root

pdf(file = "./input_data/Phylogenies/Compare_toplogies_SWSC_vs_nopart.pdf", width = 100, height = 30)
ape::comparePhylo(x = Ponerinae_SWSC_792t, y = Ponerinae_nopart_792t, plot = T, commons = FALSE, use.edge.length = TRUE)
dev.off()

test_bylocus <- ape::comparePhylo(x = Ponerinae_SWSC_792t_rooted, y = Ponerinae_bylocus_792t_rooted, plot = FALSE, use.edge.length = TRUE)
test_bylocus

View(test_bylocus$NODES)

shared_clades_ID_bylocus <- str_remove(string = test_bylocus$NODES[, 2], pattern = ".* \\(")
shared_clades_ID_bylocus <- as.numeric(str_remove(string = shared_clades_ID_bylocus, pattern = "\\)"))
length(shared_clades_ID_bylocus)

nb_tips <- length(Ponerinae_bylocus_792t_rooted$tip.label)
nb_nodes <- Ponerinae_bylocus_792t_rooted$Nnode
all_nodes_ID <- (nb_tips + 1):(nb_tips + nb_nodes)
unique_clades_ID_bylocus <- all_nodes_ID[!(all_nodes_ID %in% shared_clades_ID_bylocus)]

extract.clade(phy = Ponerinae_bylocus_792t_rooted, node = 930)

# 788 shared splits among 790 ! 99.7% identity => Splits are based on internal edges of unrooted trees
# 790 shared clades among 791 ! 99.8% identity => Clades are based on internal nodes of rooted trees with resolved root

pdf(file = "./input_data/Phylogenies/Compare_toplogies_SWSC_vs_bylocus.pdf", width = 100, height = 30)
ape::comparePhylo(x = Ponerinae_SWSC_792t, y = Ponerinae_bylocus_792t, plot = T, commons = FALSE, use.edge.length = TRUE)
dev.off()


##### 3/ Plot phylograms with unique clades for No partitioning scheme (nopart) ####

### 3.1/ Prepare template for metadata df ####

Ponerinae_nopart_792t_treedata <- as.treedata(tree = Ponerinae_nopart_792t)

nb_tips <- length(Ponerinae_nopart_792t_treedata@phylo$tip.label)
nb_nodes <- length(Ponerinae_nopart_792t_treedata@phylo$node.label)

Ponerinae_nopart_792t_node_metadata_df <- as_tibble(data.frame(node = 1:(nb_tips + nb_nodes)))
Ponerinae_nopart_792t_treedata@data <- Ponerinae_nopart_792t_node_metadata_df

## Extract support values from node labels
UFBS_values <- str_remove(string = Ponerinae_nopart_792t_treedata@phylo$node.label, pattern = "/.*")
SH_aLRT_values <- str_remove(string = Ponerinae_nopart_792t_treedata@phylo$node.label, pattern = ".*/")

plot(UFBS_values, SH_aLRT_values)

# Add support values to metadata df
Ponerinae_nopart_792t_treedata@data$UFBS[(nb_tips+1):(nb_tips + nb_nodes)] <- UFBS_values
Ponerinae_nopart_792t_treedata@data$SH_aLRT[(nb_tips+1):(nb_tips + nb_nodes)] <- SH_aLRT_values
Ponerinae_nopart_792t_treedata@data$UFBS <- as.numeric(Ponerinae_nopart_792t_treedata@data$UFBS)
Ponerinae_nopart_792t_treedata@data$SH_aLRT <- as.numeric(Ponerinae_nopart_792t_treedata@data$SH_aLRT)

## Explore UFBS values

hist(Ponerinae_nopart_792t_treedata@data$UFBS)
table(Ponerinae_nopart_792t_treedata@data$UFBS)

## Use 95 and 100 as thresholds to define color scheme

Ponerinae_nopart_792t_treedata@data$UFBS_status <- NA
Ponerinae_nopart_792t_treedata@data$UFBS_status[Ponerinae_nopart_792t_treedata@data$UFBS > 0] <- "low"
Ponerinae_nopart_792t_treedata@data$UFBS_status[Ponerinae_nopart_792t_treedata@data$UFBS >= 95] <- "high"
Ponerinae_nopart_792t_treedata@data$UFBS_status[Ponerinae_nopart_792t_treedata@data$UFBS == 100] <- "max"

table(Ponerinae_nopart_792t_treedata@data$UFBS_status)

## Explore SH-aLRT values

hist(Ponerinae_nopart_792t_treedata@data$SH_aLRT)
table(Ponerinae_nopart_792t_treedata@data$SH_aLRT)

## Use 95 as threshold to define shape scheme

Ponerinae_nopart_792t_treedata@data$SH_aLRT_status <- NA
Ponerinae_nopart_792t_treedata@data$SH_aLRT_status[Ponerinae_nopart_792t_treedata@data$UFBS < 95] <- "failed"
Ponerinae_nopart_792t_treedata@data$SH_aLRT_status[Ponerinae_nopart_792t_treedata@data$UFBS >= 95] <- "success"

table(Ponerinae_nopart_792t_treedata@data$SH_aLRT_status)

table(Ponerinae_nopart_792t_treedata@data$UFBS_status, Ponerinae_nopart_792t_treedata@data$SH_aLRT_status)
# Redundant. Only shows UFBS values?

## Add info on clade status

Ponerinae_nopart_792t_treedata@data$clade_status <- NA

Ponerinae_nopart_792t_treedata@data$clade_status[1:nb_tips] <- "tip"
Ponerinae_nopart_792t_treedata@data$clade_status[shared_clades_ID_nopart[shared_clades_ID_nopart %in% (1:(nb_tips + nb_nodes))]] <- "shared"
Ponerinae_nopart_792t_treedata@data$clade_status[unique_clades_ID_nopart] <- "unique"

table(Ponerinae_nopart_792t_treedata@data$clade_status)

## Save treedata object
saveRDS(object = Ponerinae_nopart_792t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_nopart_792t_treedata.rds")


### 3.2/ Prepare genus-groups metadata #####

# 7 genus-groups from Schmidt & Shattuck, 2014
# Platythyrea
# Pachycondyla: Pachycondyla + Simopelta + Thaumatomyrmex + Belonopelta + Mayaponera + Raspone + Dinoponera + Neoponera
# Ponera: Ponera + Diacamma + Emeryopone + Austroponera + Pseudoponera + Parvaponera + Wadeura + Cryptoponera + Iroponera + Ectomomyrmex
# Harpegnathos
# Hypoponera 
# Plectroctena: Plectoponera + Centromyrmex + Psalidomyrmex + Loboponera + Boloponera
# Odontomachus: A lot of stuff... (ca. 20 genera)

## 3.3/ Initiate metadata df ####

# List Genus-groups from Schmidt & Shattuck, 2014
genus_groups_list <- c("Platythyrea", "Pachycondyla", "Ponera", "Harpegnathos", "Hypoponera", "Plectroctena", "Odontomachus")

# List associated taxa to use to retrieve MRCA
genus_groups_MRCA_taxa_list <- list(Platythyrea = c("Platythyrea_arnoldi_D2118_CASENT0842276", "Platythyrea_schultzei_D2414a_CASENT0815127"),
                                    Pachycondyla = c("Simopelta_minima_EX2354_ANTWEB1038199", "Belonopelta_attenuata_EX2501_ICN100255"),
                                    Ponera = c("Diacamma_magdalenae_EX2692_ZRC_ENT00007536", "Austroponera_castanea_EX2322_CASENT0618835"),
                                    Harpegnathos = c("Harpegnathos_saltator_D0887_CASENT0179535", "Harpegnathos_my01_EX2981_CASENT0389113"),
                                    Hypoponera = c("Hypoponera_inexorata_EX2317_CASENT0647598", "Hypoponera_comis_D2470_CASENT0250654"),
                                    Plectroctena =  c("Centromyrmex_fugator_D2431_CASENT0235867", "Boloponera_ikemkha_EX2977_CASENT0254323"),
                                    Odontomachus = c("Neoponera_bucki_EX2455_DZUP549431", "Anochetus_talpa_D2422_CASENT0817421"))

# Initiate metadata df
Ponerinae_nopart_792t_all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_node_ID <- NA
Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_edge_ID <- NA
Ponerinae_nopart_792t_all_genus_groups_metadata_df$crown_height <- NA
Ponerinae_nopart_792t_all_genus_groups_metadata_df$stem_height <- NA

# Extract edge heights (not age, as not time-calibrated. Unit = substitutions/site)
all_edge_heights <- phytools::nodeHeights(Ponerinae_nopart_792t_treedata@phylo)
root_height <- max(all_edge_heights[, 2])
all_edge_heights <- round(-1 * all_edge_heights + root_height, 5)

Ponerinae_nopart_792t_treedata@data$Genus_group <- NA
for (i in 1:nrow(Ponerinae_nopart_792t_all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- Ponerinae_nopart_792t_all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- genus_groups_MRCA_taxa_list[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_nopart_792t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant nodes and tips
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_nopart_792t_treedata@phylo, node = MRCA_node_ID_i)
  
  # Store information of genus-group membership in the tree metadata
  Ponerinae_nopart_792t_treedata@data$Genus_group[c(MRCA_node_ID_i, all_descendants_ID_i)] <- genus_group_i
  
  # Extract crown and stem heights
  MRCA_edge_ID_i <- which(Ponerinae_nopart_792t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_height_i <- all_edge_heights[MRCA_edge_ID_i, 2]
  stem_height_i <- all_edge_heights[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and heights
  Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  Ponerinae_nopart_792t_all_genus_groups_metadata_df$crown_height[i] <- round(crown_height_i, 4)
  Ponerinae_nopart_792t_all_genus_groups_metadata_df$stem_height[i] <- round(stem_height_i, 4)
  
}

View(Ponerinae_nopart_792t_all_genus_groups_metadata_df)
table(Ponerinae_nopart_792t_treedata@data$Genus_group)

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_nopart_792t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- genus_groups_list
Ponerinae_nopart_792t_all_genus_groups_metadata_df$color <- genus_groups_colors

## Save Genus-group metadata for 792t tree
saveRDS(Ponerinae_nopart_792t_all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/Ponerinae_nopart_792t_all_genus_groups_metadata_df.rds")

## Save treedata object 
saveRDS(Ponerinae_nopart_792t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_nopart_792t_treedata.rds")


### 3.4/ Rectangular plot ####

## Load treedata object for uncalibrated UCE 792t phylogeny with nopart
Ponerinae_nopart_792t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_nopart_792t_treedata.rds")

## Load Genus-group metadata for 792t tree
Ponerinae_nopart_792t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_nopart_792t_all_genus_groups_metadata_df.rds")

# Extract MRCA to flip position of Ponera and Pachycondyla
Ponera_MRCA_node_ID <- Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_nopart_792t_all_genus_groups_metadata_df$group_name == "Ponera"]
Pachycondyla_MRCA_node_ID <- Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_nopart_792t_all_genus_groups_metadata_df$group_name == "Pachycondyla"]


## Plot tree as background
plot_792t_nopart_rect <- ggtree(tr = Ponerinae_nopart_792t_treedata,
                               layout = "rectangular",
                               ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                               size = 1.2)
# aes(color = (pp1 > 0.5))

# Flip position of Ponera and Pachycondyla
plot_792t_nopart_rect <- ggtree::flip(tree_view = plot_792t_nopart_rect, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

## Add Genus-group cladelabels/colored rectangles

# Adjust colored rectangles
extend_rect_to <- 0.25

# Loop per genus-group
for (i in 1:nrow(Ponerinae_nopart_792t_all_genus_groups_metadata_df))
{
  plot_792t_nopart_rect <- plot_792t_nopart_rect +
    
    #  Add Background color on the whole clade
    geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_nopart_792t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
                 fill = Ponerinae_nopart_792t_all_genus_groups_metadata_df$color[i],
                 # type = "gradient", # gradient.direction = 'rt',
                 # align = "right",
                 extendto = extend_rect_to,
                 alpha = .5) # +
  
  # # Add a label to the median tip. Add a bar alongside all the clade tips
  # geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
  #                 # label = all_genus_groups_metadata_df$group_name[i], 
  #                 label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
  #                 align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
  #                 geom = 'text',     # Tip of label, with or without rectangle
  #                 fontsize = 7,      # Size of the text
  #                 fill = NA,         # Color of the label background
  #                 color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
  #                 angle = 270,
  #                 hjust = 'center',
  #                 offset = 5.0,
  #                 offset.text = 0.4,
  #                 barsize = 4)    # Width of the bar
}

## Add tree over genus-group background rectangles
plot_792t_nopart_rect <- plot_792t_nopart_rect +
  
  geom_tree(layout = "rectangular",
            size = 1.2) +
  # aes(color = (pp1 > 0.5))
  
  # ## Add nod support scores as side text
  # geom_text(aes(label = UFBS), hjust = -0.2, size = 3) +
  
  # ## Add nod support scores as label on the nod
  # geom_label(aes(label = UFBS, fill = UFBS_status),
  #            label.padding = unit(0.15, "lines"),
  #            hjust = 0.5, size = 3) +
  
  ## Add symbols as nod support
  geom_nodepoint(aes(fill = UFBS_status, shape = clade_status),
                 # shape = 22,
                 size = 5) +
  
  ### Add species tip labels
  geom_tiplab(aes(label = label),
              color = "black",
              align = F,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 4,
              offset = 0.00,    # To move to the exterior. Useful if you wish to put multiple layer of labels with different offsets.
              hjust = 0)  +   # To center, left-align, right-aligned text
  
  ### Set title
  labs(title = paste0("IQTREE ML tree - 792t")) +
  
  ### Set legend fill scheme
  scale_fill_manual("UFBS support",
                    values = c("orange", "lightblue", "limegreen"),
                    labels = c("< 95%", "95% - 99%", "100%"),
                    breaks = c("low", "high", "max")
                    
  ) +
  
  ### Set legend shape scheme
  scale_shape_manual("Clades",
                     values = c(22, 24),
                     labels = c("Shared", "Unique"),
                     breaks = c("shared", "unique")
  ) +
  
  
  ### Add scale in substitutions per sites
  geom_treescale(x = 0.01, y = 780,
                 width = 0.02, # In substitutions per sites
                 label = "Substitutions / sites",
                 color = "black",
                 offset = 3.0,
                 offset.label = -3.0,
                 linesize = 4.0,
                 fontsize = 20) +
  
  # Avoid labels to be cutted
  coord_cartesian(clip = "off") +
  
  # Set aesthetics
  theme(plot.margin = unit(c(20, 120, 5, 20), "mm"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        # axis.line.x = element_line(color = "black", size = 2, vjust = 3),
        legend.title = element_text(size = 48, vjust = 3, hjust = 0, face = "bold"),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 40, face = "plain"),
        legend.key.size = unit(5.0,"lines"),
        # legend.key.width = unit(2.0, 'lines'),
        # legend.key.height = unit(5.0, 'lines'),
        legend.position.inside = c(0.10, 0.92),
        # legend.position = unit(c(1.20, 0.50), "npc"),
        # legend.key = element_rect(fill = NA))
  )


## Plot final tree
pdf(file = "./outputs/Final_phylogenies/Ponerinae_IQTREE_phylogeny_792t_nopart_rect.pdf", height = 120, width = 40)

print(plot_792t_nopart_rect)

dev.off()


##### 4/ Plot phylograms with unique clades for No partitioning scheme (bylocus) ####

### 4.1/ Prepare template for metadata df ####

Ponerinae_bylocus_792t_treedata <- as.treedata(tree = Ponerinae_bylocus_792t)

nb_tips <- length(Ponerinae_bylocus_792t_treedata@phylo$tip.label)
nb_nodes <- length(Ponerinae_bylocus_792t_treedata@phylo$node.label)

Ponerinae_bylocus_792t_node_metadata_df <- as_tibble(data.frame(node = 1:(nb_tips + nb_nodes)))
Ponerinae_bylocus_792t_treedata@data <- Ponerinae_bylocus_792t_node_metadata_df

## Extract support values from node labels
UFBS_values <- str_remove(string = Ponerinae_bylocus_792t_treedata@phylo$node.label, pattern = "/.*")
SH_aLRT_values <- str_remove(string = Ponerinae_bylocus_792t_treedata@phylo$node.label, pattern = ".*/")

plot(UFBS_values, SH_aLRT_values)

# Add support values to metadata df
Ponerinae_bylocus_792t_treedata@data$UFBS[(nb_tips+1):(nb_tips + nb_nodes)] <- UFBS_values
Ponerinae_bylocus_792t_treedata@data$SH_aLRT[(nb_tips+1):(nb_tips + nb_nodes)] <- SH_aLRT_values
Ponerinae_bylocus_792t_treedata@data$UFBS <- as.numeric(Ponerinae_bylocus_792t_treedata@data$UFBS)
Ponerinae_bylocus_792t_treedata@data$SH_aLRT <- as.numeric(Ponerinae_bylocus_792t_treedata@data$SH_aLRT)

## Explore UFBS values

hist(Ponerinae_bylocus_792t_treedata@data$UFBS)
table(Ponerinae_bylocus_792t_treedata@data$UFBS)

## Use 95 and 100 as thresholds to define color scheme

Ponerinae_bylocus_792t_treedata@data$UFBS_status <- NA
Ponerinae_bylocus_792t_treedata@data$UFBS_status[Ponerinae_bylocus_792t_treedata@data$UFBS > 0] <- "low"
Ponerinae_bylocus_792t_treedata@data$UFBS_status[Ponerinae_bylocus_792t_treedata@data$UFBS >= 95] <- "high"
Ponerinae_bylocus_792t_treedata@data$UFBS_status[Ponerinae_bylocus_792t_treedata@data$UFBS == 100] <- "max"

table(Ponerinae_bylocus_792t_treedata@data$UFBS_status)

## Explore SH-aLRT values

hist(Ponerinae_bylocus_792t_treedata@data$SH_aLRT)
table(Ponerinae_bylocus_792t_treedata@data$SH_aLRT)

## Use 95 as threshold to define shape scheme

Ponerinae_bylocus_792t_treedata@data$SH_aLRT_status <- NA
Ponerinae_bylocus_792t_treedata@data$SH_aLRT_status[Ponerinae_bylocus_792t_treedata@data$UFBS < 95] <- "failed"
Ponerinae_bylocus_792t_treedata@data$SH_aLRT_status[Ponerinae_bylocus_792t_treedata@data$UFBS >= 95] <- "success"

table(Ponerinae_bylocus_792t_treedata@data$SH_aLRT_status)

table(Ponerinae_bylocus_792t_treedata@data$UFBS_status, Ponerinae_bylocus_792t_treedata@data$SH_aLRT_status)
# Redundant. Only shows UFBS values?

## Add info on clade status

Ponerinae_bylocus_792t_treedata@data$clade_status <- NA

Ponerinae_bylocus_792t_treedata@data$clade_status[1:nb_tips] <- "tip"
Ponerinae_bylocus_792t_treedata@data$clade_status[shared_clades_ID_bylocus[shared_clades_ID_bylocus %in% (1:(nb_tips + nb_nodes))]] <- "shared"
Ponerinae_bylocus_792t_treedata@data$clade_status[unique_clades_ID_bylocus] <- "unique"

table(Ponerinae_bylocus_792t_treedata@data$clade_status)

## Save treedata object
saveRDS(object = Ponerinae_bylocus_792t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_bylocus_792t_treedata.rds")


### 4.2/ Prepare genus-groups metadata #####

# 7 genus-groups from Schmidt & Shattuck, 2014
# Platythyrea
# Pachycondyla: Pachycondyla + Simopelta + Thaumatomyrmex + Belonopelta + Mayaponera + Raspone + Dinoponera + Neoponera
# Ponera: Ponera + Diacamma + Emeryopone + Austroponera + Pseudoponera + Parvaponera + Wadeura + Cryptoponera + Iroponera + Ectomomyrmex
# Harpegnathos
# Hypoponera 
# Plectroctena: Plectoponera + Centromyrmex + Psalidomyrmex + Loboponera + Boloponera
# Odontomachus: A lot of stuff... (ca. 20 genera)

## 4.3/ Initiate metadata df ####

# List Genus-groups from Schmidt & Shattuck, 2014
genus_groups_list <- c("Platythyrea", "Pachycondyla", "Ponera", "Harpegnathos", "Hypoponera", "Plectroctena", "Odontomachus")

# List associated taxa to use to retrieve MRCA
genus_groups_MRCA_taxa_list <- list(Platythyrea = c("Platythyrea_arnoldi_D2118_CASENT0842276", "Platythyrea_schultzei_D2414a_CASENT0815127"),
                                    Pachycondyla = c("Simopelta_minima_EX2354_ANTWEB1038199", "Belonopelta_attenuata_EX2501_ICN100255"),
                                    Ponera = c("Diacamma_magdalenae_EX2692_ZRC_ENT00007536", "Austroponera_castanea_EX2322_CASENT0618835"),
                                    Harpegnathos = c("Harpegnathos_saltator_D0887_CASENT0179535", "Harpegnathos_my01_EX2981_CASENT0389113"),
                                    Hypoponera = c("Hypoponera_inexorata_EX2317_CASENT0647598", "Hypoponera_comis_D2470_CASENT0250654"),
                                    Plectroctena =  c("Centromyrmex_fugator_D2431_CASENT0235867", "Boloponera_ikemkha_EX2977_CASENT0254323"),
                                    Odontomachus = c("Neoponera_bucki_EX2455_DZUP549431", "Anochetus_talpa_D2422_CASENT0817421"))

# Initiate metadata df
Ponerinae_bylocus_792t_all_genus_groups_metadata_df <- data.frame(group_name = genus_groups_list)
Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_node_ID <- NA
Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_edge_ID <- NA
Ponerinae_bylocus_792t_all_genus_groups_metadata_df$crown_height <- NA
Ponerinae_bylocus_792t_all_genus_groups_metadata_df$stem_height <- NA

# Extract edge heights (not age, as not time-calibrated. Unit = substitutions/site)
all_edge_heights <- phytools::nodeHeights(Ponerinae_bylocus_792t_treedata@phylo)
root_height <- max(all_edge_heights[, 2])
all_edge_heights <- round(-1 * all_edge_heights + root_height, 5)

Ponerinae_bylocus_792t_treedata@data$Genus_group <- NA
for (i in 1:nrow(Ponerinae_bylocus_792t_all_genus_groups_metadata_df))
{
  # i <- 1
  # i <- 2
  
  # Extract Genus-groups
  genus_group_i <- Ponerinae_bylocus_792t_all_genus_groups_metadata_df$group_name[i]
  
  # Extract decisive taxa used to identify MRCA
  MRCA_sp_i <- genus_groups_MRCA_taxa_list[[i]]
  
  # Get MRCA node
  MRCA_node_ID_i <- ape::getMRCA(phy = Ponerinae_bylocus_792t_treedata@phylo, tip = MRCA_sp_i)
  
  # Get all current descendant nodes and tips
  all_descendants_ID_i <- phytools::getDescendants(tree = Ponerinae_bylocus_792t_treedata@phylo, node = MRCA_node_ID_i)
  
  # Store information of genus-group membership in the tree metadata
  Ponerinae_bylocus_792t_treedata@data$Genus_group[c(MRCA_node_ID_i, all_descendants_ID_i)] <- genus_group_i
  
  # Extract crown and stem heights
  MRCA_edge_ID_i <- which(Ponerinae_bylocus_792t_treedata@phylo$edge[, 2] == MRCA_node_ID_i)
  crown_height_i <- all_edge_heights[MRCA_edge_ID_i, 2]
  stem_height_i <- all_edge_heights[MRCA_edge_ID_i, 1]
  
  # Record MRCA node/edge ID and heights
  Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_node_ID[i] <- MRCA_node_ID_i
  Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_edge_ID[i] <- MRCA_edge_ID_i
  Ponerinae_bylocus_792t_all_genus_groups_metadata_df$crown_height[i] <- round(crown_height_i, 4)
  Ponerinae_bylocus_792t_all_genus_groups_metadata_df$stem_height[i] <- round(stem_height_i, 4)
  
}

View(Ponerinae_bylocus_792t_all_genus_groups_metadata_df)
table(Ponerinae_bylocus_792t_treedata@data$Genus_group)

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(Ponerinae_bylocus_792t_all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]
names(genus_groups_colors) <- genus_groups_list
Ponerinae_bylocus_792t_all_genus_groups_metadata_df$color <- genus_groups_colors

## Save Genus-group metadata for 792t tree
saveRDS(Ponerinae_bylocus_792t_all_genus_groups_metadata_df, file = "./outputs/Grafting_missing_taxa/Ponerinae_bylocus_792t_all_genus_groups_metadata_df.rds")

## Save treedata object 
saveRDS(Ponerinae_bylocus_792t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_bylocus_792t_treedata.rds")


### 4.4/ Rectangular plot ####

## Load treedata object for uncalibrated UCE 792t phylogeny with bylocus
Ponerinae_bylocus_792t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_bylocus_792t_treedata.rds")

## Load Genus-group metadata for 792t tree
Ponerinae_bylocus_792t_all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_bylocus_792t_all_genus_groups_metadata_df.rds")

# Extract MRCA to flip position of Ponera and Pachycondyla
Ponera_MRCA_node_ID <- Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_bylocus_792t_all_genus_groups_metadata_df$group_name == "Ponera"]
Pachycondyla_MRCA_node_ID <- Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_node_ID[Ponerinae_bylocus_792t_all_genus_groups_metadata_df$group_name == "Pachycondyla"]


## Plot tree as background
plot_792t_bylocus_rect <- ggtree(tr = Ponerinae_bylocus_792t_treedata,
                                layout = "rectangular",
                                ladderize = TRUE, # Option FALSE to avoid ggtree rearranging label which makes it impossible to sort labels in a similar order across trees
                                size = 1.2)
# aes(color = (pp1 > 0.5))

# Flip position of Ponera and Pachycondyla
plot_792t_bylocus_rect <- ggtree::flip(tree_view = plot_792t_bylocus_rect, node1 = Ponera_MRCA_node_ID, node2 = Pachycondyla_MRCA_node_ID)

## Add Genus-group cladelabels/colored rectangles

# Adjust colored rectangles
extend_rect_to <- 0.25

# Loop per genus-group
for (i in 1:nrow(Ponerinae_bylocus_792t_all_genus_groups_metadata_df))
{
  plot_792t_bylocus_rect <- plot_792t_bylocus_rect +
    
    #  Add Background color on the whole clade
    geom_hilight(mapping = aes(subset = node %in% c(Ponerinae_bylocus_792t_all_genus_groups_metadata_df$MRCA_node_ID[i])),
                 fill = Ponerinae_bylocus_792t_all_genus_groups_metadata_df$color[i],
                 # type = "gradient", # gradient.direction = 'rt',
                 # align = "right",
                 extendto = extend_rect_to,
                 alpha = .5) # +
  
  # # Add a label to the median tip. Add a bar alongside all the clade tips
  # geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
  #                 # label = all_genus_groups_metadata_df$group_name[i], 
  #                 label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
  #                 align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
  #                 geom = 'text',     # Tip of label, with or without rectangle
  #                 fontsize = 7,      # Size of the text
  #                 fill = NA,         # Color of the label background
  #                 color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
  #                 angle = 270,
  #                 hjust = 'center',
  #                 offset = 5.0,
  #                 offset.text = 0.4,
  #                 barsize = 4)    # Width of the bar
}

## Add tree over genus-group background rectangles
plot_792t_bylocus_rect <- plot_792t_bylocus_rect +
  
  geom_tree(layout = "rectangular",
            size = 1.2) +
  # aes(color = (pp1 > 0.5))
  
  # ## Add nod support scores as side text
  # geom_text(aes(label = UFBS), hjust = -0.2, size = 3) +
  
  # ## Add nod support scores as label on the nod
  # geom_label(aes(label = UFBS, fill = UFBS_status),
  #            label.padding = unit(0.15, "lines"),
  #            hjust = 0.5, size = 3) +
  
  ## Add symbols as nod support
  geom_nodepoint(aes(fill = UFBS_status, shape = clade_status),
                 # shape = 22,
                 size = 5) +
  
  ### Add species tip labels
  geom_tiplab(aes(label = label),
              color = "black",
              align = F,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 4,
              offset = 0.00,    # To move to the exterior. Useful if you wish to put multiple layer of labels with different offsets.
              hjust = 0)  +   # To center, left-align, right-aligned text
  
  ### Set title
  labs(title = paste0("IQTREE ML tree - 792t")) +
  
  ### Set legend fill scheme
  scale_fill_manual("UFBS support",
                    values = c("orange", "lightblue", "limegreen"),
                    labels = c("< 95%", "95% - 99%", "100%"),
                    breaks = c("low", "high", "max")
                    
  ) +
  
  ### Set legend shape scheme
  scale_shape_manual("Clades",
                     values = c(22, 24),
                     labels = c("Shared", "Unique"),
                     breaks = c("shared", "unique")
  ) +
  
  
  ### Add scale in substitutions per sites
  geom_treescale(x = 0.01, y = 780,
                 width = 0.02, # In substitutions per sites
                 label = "Substitutions / sites",
                 color = "black",
                 offset = 3.0,
                 offset.label = -3.0,
                 linesize = 4.0,
                 fontsize = 20) +
  
  # Avoid labels to be cutted
  coord_cartesian(clip = "off") +
  
  # Set aesthetics
  theme(plot.margin = unit(c(20, 120, 5, 20), "mm"), # trbl
        plot.title = element_text(color = "black", size = 50, face = "bold", hjust = 0.6, vjust = 3),
        # axis.line.x = element_line(color = "black", size = 2, vjust = 3),
        legend.title = element_text(size = 48, vjust = 3, hjust = 0, face = "bold"),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 40, face = "plain"),
        legend.key.size = unit(5.0,"lines"),
        # legend.key.width = unit(2.0, 'lines'),
        # legend.key.height = unit(5.0, 'lines'),
        legend.position.inside = c(0.10, 0.92),
        # legend.position = unit(c(1.20, 0.50), "npc"),
        # legend.key = element_rect(fill = NA))
  )


## Plot final tree
pdf(file = "./outputs/Final_phylogenies/Ponerinae_IQTREE_phylogeny_792t_bylocus_rect.pdf", height = 120, width = 40)

print(plot_792t_bylocus_rect)

dev.off()


### To do in Postprod in Illustrator

# Add genus-group names as label
# Add cladelabel for Genera in alternating colors (beware of polyphyletic groups (Use numbers))
# Replace color rectangles with background fading rectangles
# Add nicer legend for nod symbols
# Add ant images with shadow effect




