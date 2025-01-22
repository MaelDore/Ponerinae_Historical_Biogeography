##### Script 03: Graft missing taxa on the phylogeny #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Graft randomly missing taxa on the phylogeny based on known taxonomic information
# Produce a set of phylogenies that encompasses the uncertainty to use for downstream analyses

###

### Inputs

# Dataframe of fossil calibration ages used to run time calibration
# (Set of) Time-calibrated backbone phylogeny/ies built with ML inference from UCE data
# Taxonomic information for grafting missing taxa

###

### Outputs

# One main imputed phylogeny to use for downstream analyses presented in main manuscript
# A set of imputed phylogenies used to test robustness of downstream analyses to phylogenetic uncertainty

###


# Clean environment
rm(list = ls())


##### 1/ Load stuff ####

# Install rgeos from Archive
# install.packages("./packages/rgeos_0.6-4.tar.gz", repos = NULL, type = "source")

# Install lwgeom from Source
# install.packages("lwgeom", type = "source")

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(openxlsx)  # Use Rccp. No need of Java
library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges
library(phytools)
library(ape)
library(TreeSim)
library(GGally)    # To plot nice correlograms
library(magick)
library(qpdf)

### 1.2/ Load backbone phylogeny ####

## 1.2.1/ Use temporary non-calibrated phylogeny ####

Ponerinae_uncalibrated_phylogeny_792t <- read.tree(file = "./input_data/Phylogenies/ponerinae-792t-spruce-75p-iqtree-swscmerge-mfp_v2_v2.tre")
Ponerinae_uncalibrated_phylogeny_792t$tip.label

# Prune outgroups
outgroups <- c("Amblyopone_australis_D0872_CASENT0106229", "Paraponera_clavata_EX1573_CASENT0633292", "Proceratium_google_MAMI0434_CASENT0035028")
Ponerinae_uncalibrated_phylogeny_789t <- tidytree::drop.tip(object = Ponerinae_uncalibrated_phylogeny_792t, tip = outgroups)

saveRDS(object = Ponerinae_uncalibrated_phylogeny_789t, file = "./input_data/Phylogenies/Ponerinae_uncalibrated_phylogeny_789t.rds")
Ponerinae_uncalibrated_phylogeny_789t <- readRDS(file = "./input_data/Phylogenies/Ponerinae_uncalibrated_phylogeny_789t.rds")

## 1.2.2/ Load Maximum Clade Credibility tree from TreePL "posteriors" ####

# Load as treedata to keep HPD ranges
Ponerinae_MCC_phylogeny_792t_treedata <- treeio::read.beast(file = "./input_data/Phylogenies/Maximum_Clade_Credibility_tree_from_TreePL_792t.tre")
Ponerinae_MCC_phylogeny_792t_treedata@phylo$tip.label
View(Ponerinae_MCC_phylogeny_792t_treedata@data) # View node metadata
dim(Ponerinae_MCC_phylogeny_792t_treedata@data)

# Remove data in useless slots
Ponerinae_MCC_phylogeny_792t_treedata@file <- ""
Ponerinae_MCC_phylogeny_792t_treedata@treetext <- ""

# Prune outgroups
outgroups <- c("Amblyopone_australis_D0872_CASENT0106229", "Paraponera_clavata_EX1573_CASENT0633292", "Proceratium_google_MAMI0434_CASENT0035028")
Ponerinae_MCC_phylogeny_789t_treedata <- drop.tip(object = Ponerinae_MCC_phylogeny_792t_treedata, tip = outgroups)
dim(Ponerinae_MCC_phylogeny_789t_treedata@data)

saveRDS(object = Ponerinae_MCC_phylogeny_789t_treedata, file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")
Ponerinae_MCC_phylogeny_789t_treedata <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")

# Export Nexus file
ape::write.tree(phy = Ponerinae_MCC_phylogeny_789t_treedata@phylo, file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t.tree")

## 1.2.3/ Load set of backbone phylogenies from TreePL to account for uncertainty in divergence dating ####

path_prefix <- "./BASH_scripts_for_treePL/treepl_configs/tree"
path_suffix <- ".dated.tre"
TreePL_posterior_dated_trees_path <- paste0(path_prefix, 1:1000, path_suffix)

Ponerinae_all_posteriors_phylogeny_792t <- list()
for (i in seq_along(TreePL_posterior_dated_trees_path))
{
  Ponerinae_posterior_phylogeny_792t_i <- read.tree(file = TreePL_posterior_dated_trees_path[i])
  Ponerinae_all_posteriors_phylogeny_792t[[i]] <- Ponerinae_posterior_phylogeny_792t_i
}
# Assign ape multiPhylo class
class(Ponerinae_all_posteriors_phylogeny_792t) <- "multiPhylo"
table(unlist(lapply(X = Ponerinae_all_posteriors_phylogeny_792t, FUN = function (x) { length(x$tip.label) } )))

# Prune outgroups
outgroups <- c("Amblyopone_australis_D0872_CASENT0106229", "Paraponera_clavata_EX1573_CASENT0633292", "Proceratium_google_MAMI0434_CASENT0035028")
Ponerinae_all_posteriors_phylogeny_789t <- ape::drop.tip(phy = Ponerinae_all_posteriors_phylogeny_792t, tip = outgroups)
table(unlist(lapply(X = Ponerinae_all_posteriors_phylogeny_789t, FUN = function (x) { length(x$tip.label) } )))

saveRDS(object = Ponerinae_all_posteriors_phylogeny_789t, file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.rds")
Ponerinae_all_posteriors_phylogeny_789t <- readRDS(file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.rds")

# Export in NEXUS format
write.tree(phy = Ponerinae_all_posteriors_phylogeny_789t, file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.tree")

### 1.3/ Load grafting information ####

Ponerinae_Macroevolution_taxa_database <- read_excel("input_data/Ponerinae_Macroevolution_taxa_database.xlsx", 
                                                     col_types = c("text", "text", "numeric", 
                                                                   "text", "text", "numeric", "text", 
                                                                   "text", "text", "numeric", "numeric", 
                                                                   "numeric", "numeric", "text", "text", 
                                                                   "text", "text", "numeric", "numeric", 
                                                                   "text", "text", "text", "text", "numeric", 
                                                                   "numeric", "text", "numeric", "text", 
                                                                   "text", "text", "text", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "numeric", "numeric", "numeric", 
                                                                   "text", "text", "numeric", "numeric", 
                                                                   "text", "text", "numeric", "text"))
### 1.4/ Load fossil ages for time-calibration ####

# Load time-calibrated phylogeny
Ponerinae_fossil_calibrations <- read_excel("input_data/Time_calibration/Ponerinae_fossil_calibrations.xlsx")[1:12, ]


##### 2/ Run coarse time-calibration on temporary non-calibrated backbone phylogeny ####

?ape::makeChronosCalib
?ape::chronos

# Have 3 clock-models: 
   # "clock" = single speciation rate
   # correlated = rates of each branch drawn from a unique distribution
   # discrete = multi-regime models with rate shifts. Need to fix a priori the number of rate categories
   # relaxed = ???

### 2.1/ Set calibration points ####

# Identify nodes of the least inclusive clade of each fossil
# Then create line in calibration_df using min/max ages

# Root calibration
root_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                               node = "root", 
                               age.min = 100, age.max = 126,
                               interactive = FALSE, soft.bounds = FALSE)

# Fossil Odontoponera_pseudotransversa in Odontoponera genus
Odontoponera_tips <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Odontoponera")]
Odontoponera_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = Odontoponera_tips)
Odontoponera_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                       node = Odontoponera_node, 
                                       age.min = 10, age.max = 126,
                                       interactive = FALSE, soft.bounds = FALSE)

# Fossil Neoponera_vejestoria in Neoponera_foetida group
Neoponera_foetida_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Neoponera_foetida")]
Neoponera_inversa_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Neoponera_inversa")]
Neoponera_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = c(Neoponera_foetida_tip, Neoponera_inversa_tip))
Neoponera_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                       node = Neoponera_node, 
                                       age.min = 15, age.max = 126,
                                       interactive = FALSE, soft.bounds = FALSE)

# Fossil Platythyrea_dlusskyi in Platythyrea genus
Platythyrea_tips <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Platythyrea")]
Platythyrea_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = Platythyrea_tips)
Platythyrea_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                       node = Platythyrea_node, 
                                       age.min = 52, age.max = 126,
                                       interactive = FALSE, soft.bounds = FALSE)

# Fossil Ponera mayri in Ponera genus
Ponera_tips <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Ponera")]
Ponera_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = Ponera_tips)
Ponera_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                 node = Ponera_node, 
                                 age.min = 34, age.max = 126,
                                 interactive = FALSE, soft.bounds = FALSE)

# Fossil Hypoponera atavia in Hypoponera genus
Hypoponera_tips <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Hypoponera")]
Hypoponera_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = Hypoponera_tips)
Hypoponera_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                      node = Hypoponera_node, 
                                      age.min = 34, age.max = 126,
                                      interactive = FALSE, soft.bounds = FALSE)

# Fossil Odontomachus_spinifer in Odontomachus_haematodus group
Odontomachus_haematodus_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Odontomachus_haematodus")]
Odontomachus_troglodytes_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Odontomachus_troglodytes")]
Odontomachus_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = c(Odontomachus_haematodus_tip, Odontomachus_troglodytes_tip))
Odontomachus_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                    node = Odontomachus_node, 
                                    age.min = 15, age.max = 126,
                                    interactive = FALSE, soft.bounds = FALSE)

# Fossil Anochetus_ambiguus in Anochetus_emarginatus group
Anochetus_emarginatus_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Anochetus_emarginatus")]
Anochetus_striatulus_tip <- Ponerinae_uncalibrated_phylogeny_789t$tip.label[str_detect(string = Ponerinae_uncalibrated_phylogeny_789t$tip.label, pattern = "Anochetus_striatulus")]
Anochetus_node <- getMRCA(phy = Ponerinae_uncalibrated_phylogeny_789t, tip = c(Anochetus_emarginatus_tip, Anochetus_striatulus_tip))
Anochetus_calib <- makeChronosCalib(phy = Ponerinae_uncalibrated_phylogeny_789t, 
                                       node = Anochetus_node, 
                                       age.min = 15, age.max = 126,
                                       interactive = FALSE, soft.bounds = FALSE)

# Bind all calib to create calib_df
calib_df <- rbind(root_calib, Odontoponera_calib, Neoponera_calib, Platythyrea_calib, Ponera_calib, Hypoponera_calib, Odontomachus_calib, Anochetus_calib)

# Save df with calibration information from fossils
saveRDS(object = calib_df, file = "./input_data/Time_calibration/calib_df.rds")

### 2.2/ Run penalized-likelihood clock models ####

# Normally, should use different for the smoothing parameter lambda
# Should also use different number of discrete regimes
# Ideally, use rjMCMC to jump form different hyperparametrization

# Run correlated clock model with rates of each branch drawn from a unique distribution
correlated_model_output <- ape::chronos(phy = Ponerinae_uncalibrated_phylogeny_789t,
                                        lambda = 1,
                                        model = "correlated", 
                                        calibration = calib_df,
                                        control = chronos.control())
saveRDS(object = correlated_model_output, file = "./input_data/Time_calibration/correlated_model_output.rds")

# Run clock model with a unique rate
clock_model_output <- ape::chronos(phy = Ponerinae_uncalibrated_phylogeny_789t,
                                        lambda = 1,
                                        model = "clock", 
                                        calibration = calib_df,
                                        control = chronos.control())
saveRDS(object = clock_model_output, file = "./input_data/Time_calibration/clock_model_output.rds")

# Run discrete clock model with 3 regimes
discrete_model_output <- ape::chronos(phy = Ponerinae_uncalibrated_phylogeny_789t,
                                   lambda = 1,
                                   model = "discrete", 
                                   calibration = calib_df,
                                   control = chronos.control(nb.rate.cat = 3))
saveRDS(object = discrete_model_output, file = "./input_data/Time_calibration/discrete_model_output.rds")

# Run relaxed clock model
relaxed_model_output <- ape::chronos(phy = Ponerinae_uncalibrated_phylogeny_789t,
                                   lambda = 1,
                                   model = "relaxed", 
                                   calibration = calib_df,
                                   control = chronos.control())
saveRDS(object = relaxed_model_output, file = "./input_data/Time_calibration/relaxed_model_output.rds")

### 2.3/ Compare model outputs ####

# Plot histograms of rates
hist(attr(x = correlated_model_output, which = "rates")) # Should be unimodal
hist(attr(x = clock_model_output, which = "rates")) # Should be unique
hist(attr(x = discrete_model_output, which = "rates")) # Should be discrete in N categories
hist(attr(x = relaxed_model_output, which = "rates")) # Could be multimodal

# Compare rates (make a correlogram)
all_rates_df <- data.frame(correlated = attr(x = correlated_model_output, which = "rates"),
                           clock = attr(x = clock_model_output, which = "rates"),
                           # discrete = attr(discrete_model_output, which = "rates"),
                           relaxed = attr(x = relaxed_model_output, which = "rates"))
GGally::ggpairs(all_rates_df)

# Compare node ages (make a correlogram)
?ape::node.height
all_node_ages_df <- data.frame(correlated = phytools::nodeHeights(tree = correlated_model_output)[,1],
                               clock = phytools::nodeHeights(tree = clock_model_output)[,1],
                               discrete = phytools::nodeHeights(tree = discrete_model_output)[,1],
                               relaxed = phytools::nodeHeights(tree = relaxed_model_output)[,1])
GGally::ggpairs(all_node_ages_df)

# Plot phylogenies
par(mfrow = c(2,2))
plot(correlated_model_output, main = "Correlated model", show.tip.label = F)
plot(clock_model_output, main = "Clock model", show.tip.label = F)
plot(discrete_model_output, main = "Discrete model - N = 3", show.tip.label = F)
plot(relaxed_model_output, main = "Relaxed model", show.tip.label = F)


### 2.4/ Compare model fits ####

# Compare PHIIC !
correlated_PHIIC <- attr(x = correlated_model_output, which = "PHIIC")$PHIIC
clock_PHIIC <- attr(x = clock_model_output, which = "PHIIC")$PHIIC
discrete_PHIIC <- attr(x = discrete_model_output, which = "PHIIC")$PHIIC
relaxed_PHIIC <- attr(x = relaxed_model_output, which = "PHIIC")$PHIIC

models_comparison <- data.frame(model = c("correlated", "clock", "discrete", "relaxed"),
                                PHIIC = c(correlated_PHIIC, clock_PHIIC, discrete_PHIIC, relaxed_PHIIC))
# Compute delta PHIIC
models_comparison$delta_PHIIC <- models_comparison$PHIIC - min(models_comparison$PHIIC)

# Compute Akaike's weights (in %) from PHIIC
models_comparison$Akaike_weights <- round(phytools::aic.w(models_comparison$PHIIC), 3) * 100

# Display result
models_comparison

# Extract best fitted time-calibrated tree

Ponerinae_phylogeny_789t_calibrated <- correlated_model_output

# Save time-calibrated tree
saveRDS(object = Ponerinae_phylogeny_789t_calibrated, "./input_data/Phylogenies/Ponerinae_phylogeny_789t_calibrated.rds")


##### 3/ Build list of missing taxa #####

### 3.1/ Clean phylogeny tip names ####

## 3.1.1/ For rough calibrated phylogeny

Ponerinae_phylogeny_789t_calibrated <- readRDS("./input_data/Phylogenies/Ponerinae_phylogeny_789t_calibrated.rds")

Ponerinae_phylogeny_789t_calibrated$tip.label

Ponerinae_phylogeny_789t_calibrated_short_names <- Ponerinae_phylogeny_789t_calibrated
tips_short_names <- Ponerinae_phylogeny_789t_calibrated_short_names$tip.label

# Remove extraction and specimen codes
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_EX.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_CASENT.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_MAMI.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BBX\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BEB\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_D\\d{4}.*")

Ponerinae_phylogeny_789t_calibrated_short_names$tip.label <- tips_short_names

# Save roughly calibrated phylo with short names
saveRDS(object = Ponerinae_phylogeny_789t_calibrated_short_names, file = "./input_data/Phylogenies/Ponerinae_phylogeny_789t_calibrated_short_names.rds")

## 3.1.2/ For MCC phylogeny

Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label

Ponerinae_MCC_phylogeny_789t_treedata_short_names <- Ponerinae_MCC_phylogeny_789t_treedata
tips_short_names <- Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label

# Remove extraction and specimen codes
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_EX.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_CASENT.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_MAMI.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BBX\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BEB\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_D\\d{4}.*")

Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label <- tips_short_names

# Save MCC phylo with short names
saveRDS(object = Ponerinae_MCC_phylogeny_789t_treedata_short_names, file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata_short_names.rds")

## 3.1.3/ For TreePL posterior phylogenies

Ponerinae_all_posteriors_phylogeny_789t[[1]]$tip.label

Ponerinae_all_posteriors_phylogeny_789t_short_names <- Ponerinae_all_posteriors_phylogeny_789t
tips_short_names <- Ponerinae_all_posteriors_phylogeny_789t_short_names[[1]]$tip.label

# Remove extraction and specimen codes
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_EX.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_CASENT.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_MAMI.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BBX\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_BEB\\d{3}.*")
tips_short_names <- str_remove_all(string = tips_short_names, pattern = "_D\\d{4}.*")

Ponerinae_all_posteriors_phylogeny_789t_short_names <- list()
for (i in seq_along(Ponerinae_all_posteriors_phylogeny_789t))
{
  Phylo_i <- Ponerinae_all_posteriors_phylogeny_789t[[i]]
  Phylo_i$tip.label <- tips_short_names
  Ponerinae_all_posteriors_phylogeny_789t_short_names[[i]] <- Phylo_i
}
class(Ponerinae_all_posteriors_phylogeny_789t) <- "multiPhylo" 

Ponerinae_all_posteriors_phylogeny_789t[[i]]$tip.label
Ponerinae_all_posteriors_phylogeny_789t_short_names[[i]]$tip.label

# Save TreePL posterior phylogenies with short names
saveRDS(object = Ponerinae_all_posteriors_phylogeny_789t_short_names, file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t_short_names.rds")

View(cbind(Ponerinae_phylogeny_789t_calibrated$tip.label, Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label, Ponerinae_all_posteriors_phylogeny_789t[[1]]$tip.label))
View(cbind(Ponerinae_phylogeny_789t_calibrated_short_names$tip.label, Ponerinae_MCC_phylogeny_789t_treedata_short_names@phylo$tip.label, Ponerinae_all_posteriors_phylogeny_789t_short_names[[1]]$tip.label))


### 3.2/ Clean grafting MRCA tips information ####

# Extract grafting information
Missing_taxa_df <- Ponerinae_Macroevolution_taxa_database %>% 
  filter(In_phylogeny == F) %>%
  select(Current_name, Conservative_clade_Terminals_with_MRCA, Exclusive_subclades_Terminals_with_MRCA) %>% 
  rename(Taxa_name = Current_name,
         MRCA_tips = Conservative_clade_Terminals_with_MRCA, 
         Exclusive_subclades = Exclusive_subclades_Terminals_with_MRCA)

# Add Taxa name to the list of missing taxa
Missing_taxa_list <- list(Taxa_name = Missing_taxa_df$Taxa_name)

# Extract MRCA tips
MRCA_tips_df <- as.data.frame(str_split(string = Missing_taxa_df$MRCA_tips, pattern = ",", simplify = T))
# Remove space
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = " "))
# Remove extraction and specimen codes
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_EX.*"))
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_CASENT.*"))
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_MAMI.*"))
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_BBX\\d{3}.*"))
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_BEB\\d{3}.*"))
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_D\\d{4}.*"))

# Remove CAPITALS
MRCA_tips_df <- as.data.frame(apply(X = MRCA_tips_df, MARGIN = 2, FUN = str_to_title))

# Replace "" by NA
MRCA_tips_df[replace_na(data = (MRCA_tips_df[, 1] == ""), replace = F), 1] <- NA
MRCA_tips_df[replace_na(data = (MRCA_tips_df[, 2] == ""), replace = F), 2] <- NA

# Match names with the short name tip labels
MRCA_updated_tips_df <- MRCA_tips_df
for(i in 1:nrow(MRCA_tips_df))
{
  # Tip n°1
  focal_tip <- MRCA_tips_df[i, 1]
  if (is.na(focal_tip))
  {
    MRCA_updated_tips_df[i, 1] <- NA
  } else {
    tip_ID <- which(Ponerinae_phylogeny_789t_calibrated_short_names$tip.label == focal_tip)
    if (length(tip_ID) == 1)
    {
      updated_tip <- Ponerinae_phylogeny_789t_calibrated$tip.label[tip_ID]
      MRCA_updated_tips_df[i, 1] <- updated_tip
    }
    if (length(tip_ID) == 0)
    {
      MRCA_updated_tips_df[i, 1] <- "no match"
    }
    if (length(tip_ID) > 1)
    {
      MRCA_updated_tips_df[i, 1] <- "multiple matches"
    }
  }
  
  # Tip n°2
  focal_tip <- MRCA_tips_df[i, 2]
  if (is.na(focal_tip))
  {
    MRCA_updated_tips_df[i, 2] <- NA
  } else {
    tip_ID <- which(Ponerinae_phylogeny_789t_calibrated_short_names$tip.label == focal_tip)
    if (length(tip_ID) == 1)
    {
      updated_tip <- Ponerinae_phylogeny_789t_calibrated$tip.label[tip_ID]
      MRCA_updated_tips_df[i, 2] <- updated_tip
    }
    if (length(tip_ID) == 0)
    {
      MRCA_updated_tips_df[i, 2] <- "no match"
    }
    if (length(tip_ID) > 1)
    {
      MRCA_updated_tips_df[i, 2] <- "multiple matches"
    }
  }
}

table(MRCA_updated_tips_df[, 1] == "multiple matches")
table(MRCA_updated_tips_df[, 1] == "no match")

table(MRCA_updated_tips_df[, 2] == "multiple matches")
table(MRCA_updated_tips_df[, 2] == "no match")
# No mismatch

View(MRCA_tips_df)
View(MRCA_updated_tips_df)

# Ponerinae_phylogeny_789t_calibrated$tip.label[order(Ponerinae_phylogeny_789t_calibrated$tip.label)]


### 3.3/ Retrieve Genus MRCA for missing taxa with no grafting information ####

names(MRCA_updated_tips_df) <- c("Tip_1", "Tip_2")

# Flag missing taxa with no grafting information
MRCA_updated_tips_df$Default_grafting <- is.na(MRCA_updated_tips_df$Tip_1)
table(MRCA_updated_tips_df$Default_grafting)

# Provide tips to retrieve MRCA of the Genus by default
All_genera_tips <- str_split(string = Ponerinae_phylogeny_789t_calibrated$tip.label, pattern = "_", simplify = T)[,1]
nb_tips <- length(Ponerinae_phylogeny_789t_calibrated$tip.label)

for (i in 1:nrow(MRCA_updated_tips_df))
{
  # i <- 4
  
  Need_default_grafting_info <- MRCA_updated_tips_df$Default_grafting[i]
  
  # Only for taxa with no grafting information
  if (Need_default_grafting_info)
  {
    # Extract taxa name and Genus
    Focal_taxon <- Missing_taxa_list$Taxa_name[i]
    Focal_genus <- str_split(string = Focal_taxon, pattern = "_", simplify = T)[1,1]
    # Find all Genus tips in the phylogeny
    All_tips_in_genus <- Ponerinae_phylogeny_789t_calibrated$tip.label[All_genera_tips == Focal_genus]
    
    # Case where the MRCA is a tip because there is only one terminal of this Genus in the phylogeny
    if (length(All_tips_in_genus) == 1)
    {
      # Provide only this tip name
      Genus_MRCA_tips <- c(All_tips_in_genus, NA)
    } else {
    # Case where the MRCA is an internal node, thus there are multiple tips
      
      # Extract MRCA
      Genus_MRCA_node <- ape::getMRCA(phy = Ponerinae_phylogeny_789t_calibrated, tip = All_tips_in_genus)
      # Find next two descending nodes (i.e., partitions)
      nodes_RL <- Ponerinae_phylogeny_789t_calibrated$edge[(Ponerinae_phylogeny_789t_calibrated$edge[,1] == Genus_MRCA_node),2]
      # Identify descending tips from each partition
      descendant_tips_R <- phytools::getDescendants(tree = Ponerinae_phylogeny_789t_calibrated, node = nodes_RL[1])
      descendant_tips_R <- descendant_tips_R[descendant_tips_R <= nb_tips]
      descendant_tips_L <- phytools::getDescendants(tree = Ponerinae_phylogeny_789t_calibrated, node = nodes_RL[2])
      descendant_tips_L <- descendant_tips_L[descendant_tips_L <= nb_tips]
      # Select one tip from each partition
      Genus_MRCA_tips <- c(descendant_tips_R[1], descendant_tips_L[1]) 
      Genus_MRCA_tips <- Ponerinae_phylogeny_789t_calibrated$tip.label[Genus_MRCA_tips]
    }
    
    # Fill MRCA pairs information
    MRCA_updated_tips_df[i, c("Tip_1", "Tip_2")] <- Genus_MRCA_tips
  }
}

View(MRCA_updated_tips_df)

# Add MRCA_tips to the list of missing taxa
MRCA_tips_list <- split(x = MRCA_updated_tips_df[, c("Tip_1", "Tip_2")], f = seq(nrow(MRCA_updated_tips_df))) 
MRCA_tips_list <- lapply(X = MRCA_tips_list, FUN = unlist)
# Remove NA in the case of single tips
MRCA_tips_list <- lapply(X = MRCA_tips_list, FUN = function (x) {x[!is.na(x)]})
names(MRCA_tips_list) <- Missing_taxa_list$Taxa_name
Missing_taxa_list <- append(x = Missing_taxa_list, values = list(MRCA_tips = MRCA_tips_list))

# Save Missing_taxa_list
saveRDS(object = Missing_taxa_list, file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")


### 3.4/ Clean excluding subclades MRCA tips information ####

## Currently work only for a single exclusive subclade per missing taxa

# Extract MRCA tips
Excluding_subclades_tips_df <- as.data.frame(str_split(string = Missing_taxa_df$Exclusive_subclades, pattern = ",", simplify = T))
# Remove space
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = " "))
# Remove extraction and specimen codes
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_EX.*"))
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_CASENT.*"))
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_MAMI.*"))
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_BBX\\d{3}.*"))
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_BEB\\d{3}.*"))
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_remove_all, pattern = "_D\\d{4}.*"))

# Remove CAPITALS
Excluding_subclades_tips_df <- as.data.frame(apply(X = Excluding_subclades_tips_df, MARGIN = 2, FUN = str_to_title))

# Replace "" by NA
Excluding_subclades_tips_df[replace_na(data = (Excluding_subclades_tips_df[, 1] == ""), replace = F), 1] <- NA
Excluding_subclades_tips_df[replace_na(data = (Excluding_subclades_tips_df[, 2] == ""), replace = F), 2] <- NA

# Match names with the short name tip labels
Excluding_subclades_updated_tips_df <- Excluding_subclades_tips_df
for(i in 1:nrow(Excluding_subclades_tips_df))
{
  # Tip n°1
  focal_tip <- Excluding_subclades_tips_df[i, 1]
  if (is.na(focal_tip))
  {
    Excluding_subclades_updated_tips_df[i, 1] <- NA
  } else {
    tip_ID <- which(Ponerinae_phylogeny_789t_calibrated_short_names$tip.label == focal_tip)
    if (length(tip_ID) == 1)
    {
      updated_tip <- Ponerinae_phylogeny_789t_calibrated$tip.label[tip_ID]
      Excluding_subclades_updated_tips_df[i, 1] <- updated_tip
    }
    if (length(tip_ID) == 0)
    {
      Excluding_subclades_updated_tips_df[i, 1] <- "no match"
    }
    if (length(tip_ID) > 1)
    {
      Excluding_subclades_updated_tips_df[i, 1] <- "multiple matches"
    }
  }
  
  # Tip n°2
  focal_tip <- Excluding_subclades_tips_df[i, 2]
  if (is.na(focal_tip))
  {
    Excluding_subclades_updated_tips_df[i, 2] <- NA
  } else {
    tip_ID <- which(Ponerinae_phylogeny_789t_calibrated_short_names$tip.label == focal_tip)
    if (length(tip_ID) == 1)
    {
      updated_tip <- Ponerinae_phylogeny_789t_calibrated$tip.label[tip_ID]
      Excluding_subclades_updated_tips_df[i, 2] <- updated_tip
    }
    if (length(tip_ID) == 0)
    {
      Excluding_subclades_updated_tips_df[i, 2] <- "no match"
    }
    if (length(tip_ID) > 1)
    {
      Excluding_subclades_updated_tips_df[i, 2] <- "multiple matches"
    }
  }
}

table(Excluding_subclades_updated_tips_df[, 1] == "multiple matches")
table(Excluding_subclades_updated_tips_df[, 1] == "no match")

table(Excluding_subclades_updated_tips_df[, 2] == "multiple matches")
table(Excluding_subclades_updated_tips_df[, 2] == "no match")
# No mismatch

View(Excluding_subclades_tips_df)
View(Excluding_subclades_updated_tips_df)

# Ponerinae_phylogeny_789t_calibrated$tip.label[order(Ponerinae_phylogeny_789t_calibrated$tip.label)]

Excluding_subclades_updated_tips_df

names(Excluding_subclades_updated_tips_df) <- paste0("Tip_", 1:ncol(Excluding_subclades_updated_tips_df))

# Add Excluding_subclades_tips to the list of missing taxa
Excluding_subclades_tips_list <- split(x = Excluding_subclades_updated_tips_df, f = seq(nrow(Excluding_subclades_updated_tips_df))) 
Excluding_subclades_tips_list <- lapply(X = Excluding_subclades_tips_list, FUN = unlist)
# Remove NA in the case of single tips
Excluding_subclades_tips_list <- lapply(X = Excluding_subclades_tips_list, FUN = function (x) {x[!is.na(x)]})
names(Excluding_subclades_tips_list) <- Missing_taxa_list$Taxa_name
# Regroup per pairs associated to a subclade j
for (i in 1:length(Excluding_subclades_tips_list))
{
  # i <- 1
  # i <- 638
  Excluding_subclades_tips_list_i <- Excluding_subclades_tips_list[[i]]
  
  # If no exclusive subclades, fill with NULL
  if (length(Excluding_subclades_tips_list_i) == 0)
  {
    Excluding_subclades_tips_list[[i]] <- list()
  } else {
  # If exclusive subclades, regroup tips per pairs/subclade
    Excluding_subclades_tips_paired_list_i <- list()
    j <- 1
    while(length(Excluding_subclades_tips_list_i) > 0)
    {
      Excluding_subclades_tips_paired_list_j <- Excluding_subclades_tips_list_i[1:2]
      old_names_i <- names(Excluding_subclades_tips_paired_list_i)
      Excluding_subclades_tips_paired_list_i <- append(x = Excluding_subclades_tips_paired_list_i, values = list(Excluding_subclades_tips_paired_list_j))
      names(Excluding_subclades_tips_paired_list_i) <- c(old_names_i, paste0("Subclade_",j))
      
      # Remove pair j
      Excluding_subclades_tips_list_i <- Excluding_subclades_tips_list_i[-c(1:2)]
      
      # Increment number of exclusive subclades
      j <- j + 1
    }
    Excluding_subclades_tips_list[[i]] <- Excluding_subclades_tips_paired_list_i
  }
}

Excluding_subclades_tips_list

# Add exclusive subclades list to the main missing taxa list object
Missing_taxa_list <- append(x = Missing_taxa_list, values = list(Excluding_subclades = Excluding_subclades_tips_list))

Missing_taxa_list$Excluding_subclades

# Save Missing_taxa_list
saveRDS(object = Missing_taxa_list, file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")


##### 4/ Functions to graft missing taxa #####

### 4.1/ Function to draw randomly a branch within a clade ####

# phy <- sim.bdtree(b = 1, d = 0, n = 10)
# plot(phy)
# nodelabels()
# edgelabels()
# 
# MRCA_node <- 18
# MRCA_node <- 5
# MRCA_node <- 12
# exclusive_subclades_MRCA_nodes <- c(13, 16)

drawing_random_branch <- function (phy, 
                                   MRCA_node, # To provide the MRCA node of the total clades where to draw the selected branch from
                                   which_node = "tipward", # Chose to provide the rootward node or the tipward node of the selected branch
                                   include_stem = T, # Include the stem branch among the set of possible branches
                                   exclusive_subclades_MRCA_nodes = NULL, # To provide MRCA nodes of exclusive subclades where not to draw the selected branch
                                   seed = 1234) # Set the seed
{
  # Node ID: N_tips (N = N_tips), root (N = 1), N_internal (N = N_tips - 2) 
  # Total = N_tips - 1
  
  # N_edges = N_tips - 2 because the root node does not have an edge
  # By default, edges are not ordered as internal/terminal, but rather cladewise

  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }
  
  # Extract descendant nodes from MRCA
  descendant_nodes <- phytools::getDescendants(tree = phy, node = MRCA_node)
  
  # If MRCA is the root, prevent from drawing the stem (does not exist)
  root_node <- length(phy$tip.label) + 1
  if (MRCA_node == root_node)
  {
    include_stem <- F
  }
  
  # Include node of the stem if needed (= MRCA node)
  if (include_stem)
  {
    descendant_nodes <- unique(c(descendant_nodes, MRCA_node))
  }
  # Extract branches ID
  branches_ID <- match(x = descendant_nodes, table = phy$edge[, 2])
  # Extract branch length
  branches_length <- setNames(object = phy$edge.length[branches_ID], nm = branches_ID)

  # Remove branch from exclusive subclades from the draw
  if (!is.null(exclusive_subclades_MRCA_nodes))
  {
    # Initiate vector of forbidden branches
    excluded_branches_ID <- c()
    # Loop per exclusive subclade
    for (i in 1:length(exclusive_subclades_MRCA_nodes))
    {
      # Extract exclusive subclade MRCA node
      exclusive_subclade_MRCA_node_i <- exclusive_subclades_MRCA_nodes[i]
      # Extract descendant nodes from exclusive subclade MRCA
      exclusive_subclade_descendant_nodes_i <- phytools::getDescendants(tree = phy, node = exclusive_subclade_MRCA_node_i)
      # Extract forbidden branches ID
      excluded_branches_ID_i <- match(x = exclusive_subclade_descendant_nodes_i, table = phy$edge[, 2])
      # Store in vector of forbidden branches
      excluded_branches_ID <- c(excluded_branches_ID, excluded_branches_ID_i)
    }
    
    # Remove excluded branches from the draw
    branches_ID <- setdiff(branches_ID, excluded_branches_ID)
    branches_length <- branches_length[!(names(branches_length) %in% excluded_branches_ID)]
  }
  
  # Draw randomly the branch
  if (length(branches_ID) == 1)
  {
    selected_branch_ID <- branches_ID
  } else {
    selected_branch_ID <- sample(x = branches_ID, size = 1, prob = branches_length)
  }
  
  # Provide the tipward or rootward node
  if (which_node == "tipward")
  {
    selected_node <- phy$edge[selected_branch_ID, 2]
  } else {
    selected_node <- phy$edge[selected_branch_ID, 1]
  }
  
  # Return selected node ID
  return(selected_node)
}
  

### 4.2/ Main function to graft missing taxa ####


# Arguments for MRCA, include stem or not, exclusive subclades 
# Argument for the type of grafting
  # Purely random = "uniform"
     # branch selection probability based on branch length (excluding exclusive subclades), position on the branch drawn follows a uniform prior
  # Accounting for subclade diversification dynamics = fit birth-death models per subclades
     # Per taxa: "BD_taxa"
         # Fit BD with one missing taxa
         # branch selection probability based on branch length (excluding exclusive subclades), but position on the branch drawn from exponential distribution according to subclade speciation rate (not a uniform prior!)
     # Per subclades/MRCA: "BD_clades"
         # Account for the proper number of missing taxa during the BD fit
         # Estimate age of missing taxa within a clade with TreeSim::corsim, then select randomly the branch to graft (Prevent branching to occur within the exclusive(s) subclade(s))
         # Need to graft taxa per MRCA (if nested taxa, need to graft the nested clades first)

grafting_missing_taxa <- function (phy, # Backbone phylogeny
                                   missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                   include_stem = T, # Include or not the stem of the total clade when grafting?
                                   use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                   method = "uniform", # Either "uniform", "BD_taxa", or "BD_clades"
                                   verbose = T, # To print progress every 100 taxa
                                   seed = 1234) # Set the seed
{
  # Check that tips for locating MRCA are all found on the phylogeny
  all_MRCA_tips <- unique(unlist(missing_taxa_list$MRCA_tips))
  if(!all(all_MRCA_tips %in% phy$tip.label))
  {
    non_matching_tips <- all_MRCA_tips[!(all_MRCA_tips %in% phy$tip.label)]
    stop(paste0("Some tip labels provided to locate MRCA nodes of grafting clades are missing from the phylogeny\n",
                "Missing tips: ", paste(non_matching_tips, collapse = ", ")))
  }
  
  # Check that tips for locating excluding_subclades are all found on the phylogeny
  all_MRCA_tips <- unique(unlist(missing_taxa_list$Excluding_subclades))
  if(!all(all_MRCA_tips %in% phy$tip.label))
  {
    non_matching_tips <- all_MRCA_tips[!(all_MRCA_tips %in% phy$tip.label)]
    stop(paste0("Some tip labels provided to locate nodes of excluding subclades are missing from the phylogeny\n",
                "Missing tips: ", paste(non_matching_tips, collapse = ", ")))
  }
  
  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }
  
  # Randomize the order of taxa in the missing_taxa_list so it has no influence on the pattern of grafting
  random_order <- sample(x = 1:length(missing_taxa_list$Taxa_name), size = length(missing_taxa_list$Taxa_name), replace = F)
  missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

  # For "uniform" and BD_clades => Loop per taxa
  
  if (method %in% c("uniform", "BD_taxa"))
  {
    for (i in 1:length(missing_taxa_list$Taxa_name))
    {
      # i <- 1
      
      # Extract taxon name
      focal_taxon_name <- missing_taxa_list$Taxa_name[i]
      
      # Extract MRCA node
      if (length(missing_taxa_list$MRCA_tips[[i]]) == 2) # If using pairs
      {
        MRCA_node <- ape::getMRCA(phy = phy, tip = missing_taxa_list$MRCA_tips[[i]])
      }
      if (length(missing_taxa_list$MRCA_tips[[i]]) == 1) # If using a single tip
      {
        MRCA_node <- which(phy$tip.label == missing_taxa_list$MRCA_tips[[i]])
      }
      
      # If using exclusive subclades, retrieve the MRCA nodes of forbidden subclades
      if (use_exclusive_subclades)
      {
        # Extract exclusive subclades pairs for taxon i
        exclusive_subclades_pairs <- missing_taxa_list$Exclusive_subclades[[i]]
        
        # Initiate vector of excluding subclades MRCA nodes
        exclusive_subclades_MRCA_nodes <- c()
        
        # Loop per pairs of tips defining excluding subclades
        for (j in 1:length(exclusive_subclades_pairs))
        {
          # j <- 1
          
          # Extract pair j
          exclusive_subclades_pair_j <- exclusive_subclades_pairs[[j]]
          # Retrieve MRCA node of subclade j
          exclusive_subclades_MRCA_node_j <- ape::getMRCA(phy = phy, tip = exclusive_subclades_pair_j)
          # Store MRCA node of subclade j
          exclusive_subclades_MRCA_nodes <- c(exclusive_subclades_MRCA_nodes, exclusive_subclades_MRCA_node_j)
        }
        
      } else {
        exclusive_subclades_MRCA_nodes <- NULL
      }
      
      # Draw node of random branch used to graft the taxa
      selected_node_ID <- drawing_random_branch(phy = phy, MRCA_node = MRCA_node, include_stem = include_stem, exclusive_subclades_MRCA_nodes = exclusive_subclades_MRCA_nodes, seed = NULL)
      # Extract branch length
      selected_branch_ID <- match(x = selected_node_ID, table = phy$edge[, 2])
      selected_branch_length <- phy$edge.length[selected_branch_ID]
      
      # Draw position of the new tip on the selected branch
      if (method == "uniform")
      {
        # For "uniform" method, draw a random position following a uniform PDF
        position <- runif(n = 1, min = 0, max = selected_branch_length)
      } else {
        # For "BD_taxa" method, draw a random position following an exponential law in accordance to the fitted speciation rate in the clade
        
        # Identify number of descending taxa in the stem clade
        nb_desc_tips <- length(geiger::tips(phy = phy, node = MRCA_node))
        
        ## Run BD model according to the number of descending taxa
        # Run BD model on a subclade of at least 3 tips encomapssing the stem group
        
        if (nb_desc_tips == 1)
        {
          ## Case with a single tip
          # Need to move backward twice to ensure having a clade of at least 3 tips
          
          # Retrieve MRCA node of the extended substitute clade
          BD_clade_node <- phy$edge[(phy$edge[,2] == MRCA_node), 1]
          BD_clade_node <- phy$edge[(phy$edge[,2] == BD_clade_node), 1]
        }
        if (nb_desc_tips == 2)
        {
          ## Case with two tips
          # Need to move backward once to ensure having a clade of at least 3 tips
          
          # Retrieve MRCA node of the extended substitute clade
          BD_clade_node <- phy$edge[(phy$edge[,2] == MRCA_node), 1]
        }
        if (nb_desc_tips > 2)
        {
          ## Case with more than two tips
          # Can use the MRCA node
          BD_clade_node <- MRCA_node
        }
        
        # Extract the crown group or the extended substitute tree
        BD_clade_tree <- ape::extract.clade(phy = phy, node = BD_clade_node)
        
        # plot(BD_clade_tree)
        
        # Reidentify number of terminals in the BD clade
        nb_desc_tips_BD_clade <- length(BD_clade_tree$tip.label)
        
        # Estimate sampling fraction accounting for one missing taxa (the one we are grafting)
        sampling_fraction <- nb_desc_tips_BD_clade / (nb_desc_tips_BD_clade + 1) 
        # Fit BD model in the BD clade, accounting for one missing taxa
        BD_clade_tree <- force.ultrametric(tree = BD_clade_tree, method = "nnls", message = F)
        lambda_BD <- phytools::fit.bd(tree = BD_clade_tree, rho = sampling_fraction)$b
        
        # Draw a random position following an exponential law in accordance to the fitted speciation rate in the clade
        # plot(x = seq(from = 0, to = 5, by = 0.01), y = dexp(x = seq(from = 0, to = 5, by = 0.01), rate = lambda_BD), type = "l", xlab = "position", ylab = "probability density")
        # Position cannot excess the branch length
        position <- Inf
        while (position > selected_branch_length)
        {
          position <- rexp(n = 1, rate = lambda_BD)
        }
      }
      
      # Extract root age
      root_age <- max(phytools::nodeHeights(tree = phy)[, 2])
      
      # Compute branch length needed to remain ultrametric
      new_branch_length <- position + (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
      
      # Bind new taxa in the tree
      phy <- phytools::bind.tip(tree = phy, tip.label = focal_taxon_name, edge.length = new_branch_length, where = selected_node_ID, position = position)
      
      # Check phylo
      # plot(phy)
      
      # Print progress every 100 tips
      if (verbose & (i %% 100 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted n°",i,"/",length(missing_taxa_list$Taxa_name),"\n"))
      }
    }
  }
  
  # For "BD_clades" => Loop per stem clades used for grafting, starting with the nested ones

  if (method == "BD_clades")
  {
    ## Get info per stem clades
    
    # Identify all stem clades used for grafting
    MRCA_tips_merged <- unlist(lapply(missing_taxa_list$MRCA_tips, FUN = paste0, collapse = ","))
    # Group missing taxa per stem clades
    MRCA_tips_groups <- as.factor(MRCA_tips_merged)
    as.numeric(MRCA_tips_groups)
    
    # Get info per stem clades
    groups_ID <- 1:nlevels(MRCA_tips_groups)
    group_indices <- match(x = groups_ID, table = as.numeric(MRCA_tips_groups))
    
    stem_clades_list <- list(ID = groups_ID,
                             MRCA_tips = missing_taxa_list$MRCA_tips[group_indices])
    names(stem_clades_list$MRCA_tips) <- groups_ID
    
    # Identify MRCA nodes & descending nodes
    stem_clades_list$MRCA_nodes <- list()
    stem_clades_list$descendant_nodes <- list()
    for (i in 1:length(groups_ID))
    {
      if (length(stem_clades_list$MRCA_tips[[i]]) > 1) 
      {
        # Case with pairs of tips
        stem_clades_list$MRCA_nodes[[i]] <- ape::getMRCA(phy = phy, tip = stem_clades_list$MRCA_tips[[i]])
        stem_clades_list$descendant_nodes[[i]] <- phytools::getDescendants(tree = phy, node = stem_clades_list$MRCA_nodes[[i]])
      } else { 
        # Case with a single tip => MRCA node = tip node ; no descendants
        stem_clades_list$MRCA_nodes[[i]] <- which(phy$tip.label == stem_clades_list$MRCA_tips[[i]])
        stem_clades_list$descendant_nodes[[i]] <- c()
      }
      
    }
    names(stem_clades_list$MRCA_nodes) <- stem_clades_list$ID
    names(stem_clades_list$descendant_nodes) <- stem_clades_list$ID
    
    ## Order stem groups such as the nested ones are selected fist

    # Initiate vector of stem clades order
    stem_clades_order <- c()
    stem_clades_list_temp <- stem_clades_list
    while (length(stem_clades_order) < length(groups_ID))
    {
      # Initiate stopping criterion as check for nestedness of the stem clade
      # Selected stem clade must have no other stem clade nested within
      no_nested_subclades_within <- F
      j <- 0 # Initiate counter
      while (!no_nested_subclades_within)
      {
        j <- j + 1
        target_stem_clade_descendant_nodes <- unlist(stem_clades_list_temp$descendant_nodes[j])
        target_MCRA_other_stem_clades <- unlist(stem_clades_list_temp$MRCA_nodes)[(j+1):length(stem_clades_list_temp$MRCA_nodes)]
        no_nested_subclades_within <- !any(target_MCRA_other_stem_clades %in% target_stem_clade_descendant_nodes)
      }
      # Extract selected stem clade ID
      selected_clade_ID <- stem_clades_list_temp$ID[j]
      # Assign stem clade in the order
      stem_clades_order <- c(stem_clades_order, selected_clade_ID)
      # Remove it from selection
      stem_clades_list_temp <- lapply(X = stem_clades_list_temp, FUN = function (x) { x[-j] })

      # Stop when all stem clades are ordered
    }
    
    ## Loop per stem clades used for grafting, starting with the nested ones
    
    # Initiate number of taxa/stem clades grafted
    nb_taxa_grafted <- 0
    nb_stem_clades_grafted <- 0
    
    for (i in stem_clades_order)
    {
      # i <- stem_clades_order[1]
      
      ## Subset missing_taxa_list for the given stem clade
      missing_taxa_ID_stem_clade_i <- which(as.numeric(MRCA_tips_groups) == i)
      missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[missing_taxa_ID_stem_clade_i] })
      
      # Extract MRCA tips
      MRCA_tips_i <- stem_clades_list$MRCA_tips[[i]]
      
      ## Draw divergence times according to the number of tips used to define the crown clade, and within the crown clade
      if (length(MRCA_tips_i) > 1)
      {
        ## Case with pairs of tips used to define the crown clade
        
        # Extract MRCA node using tips (cannot use previous extract as the tree as been updated with the missing taxa from the previous stem group)
        MRCA_node_i <- ape::getMRCA(phy = phy, tip = MRCA_tips_i)
        
        # Extract the MRCA crown clade tree
        crown_clade_tree_i <- ape::extract.clade(phy = phy, node = MRCA_node_i)
        
        # plot(crown_clade_tree_i)
        
        # Extract divergence times in the crown clade tree
        crown_age_i <- max(phytools::nodeHeights(crown_clade_tree_i)[,2])
        divergence_times_i <- crown_age_i - phytools::nodeHeights(crown_clade_tree_i)[,1]
        
        # Set upper bound to drawn divergence age as the crown age
        upper_age_i <- crown_age_i
        
        # Add stem branch by using the stem age as the upper limit to draw divergence times instead of the MRCA node age
        if (include_stem)
        {
          full_tree_age <- max(phytools::nodeHeights(phy)[,2])
          stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
          stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
          upper_age_i <- stem_age_i # Set upper bound to drawn divergence age as the stem age, thus include the possibility of branching on the stem branch
        }
        
        # Compute sampling fraction
        ntips_present <- length(crown_clade_tree_i$tip.label)
        ntips_missing <- length(missing_taxa_list_i$Taxa_name)
        rho <- ntips_present / (ntips_present + ntips_missing)
        
        ## Fit BD on crown clade, or substitute extended crown clade according to number of tips in the crown clade
        if (ntips_present > 2)
        {
          ## Case with more than 2 lineages in the crown clade: fit a BD on the crown clade
          
          # Estimate BD rates
          crown_clade_tree_i <- force.ultrametric(tree = crown_clade_tree_i, method = "nnls", message = F)
          crown_clade_BD_i <- phytools::fit.bd(tree = crown_clade_tree_i, rho = rho)
          
          ## Draw divergence times
          
          # Draw divergence ages of missing taxa
          # ?TreeSim::corsim
          new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = crown_clade_BD_i$b, mu = crown_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
          divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)
          
          # divergence_times_missing_taxa_i
          
        } else {
          
          ## Case with only 2 lineages in the crown clade: fit a BD on a substitute clade of at least 3 tips, encompassing the crown clade
          
          # Retrieve MRCA node of the extended substitute clade
          substitute_clade_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i), 1]
          
          # Extract the substitute tree
          substitute_clade_tree_i <- ape::extract.clade(phy = phy, node = substitute_clade_node_i)
          
          # plot(crown_clade_tree_i)
          # plot(substitute_clade_tree_i)
          
          # Reidentify number of terminals in the BD clade
          nb_desc_tips_substitute_clade <- length(substitute_clade_tree_i$tip.label)
          # Adjust sampling fraction
          rho_BD_clade <- nb_desc_tips_substitute_clade / (nb_desc_tips_substitute_clade + ntips_missing)
          
          # Estimate BD rates from the substitute tree including at least 3 tips
          substitute_clade_tree_i <- force.ultrametric(tree = substitute_clade_tree_i, method = "nnls", message = F)
          substitute_clade_BD_i <- phytools::fit.bd(tree = substitute_clade_tree_i, rho = rho_BD_clade)
          
          ## Draw divergence times
          
          # Draw divergence ages of missing taxa
          # ?TreeSim::corsim
          new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
          divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)
          
          # divergence_times_missing_taxa_i

        }
        
      } else {
        
        ## Case with a single tip => MRCA node = tip node ; no descendants
        
        # MRCA node = unique tip node
        MRCA_node_i <- which(phy$tip.label == MRCA_tips_i)
        
        # Stem clade tree is a single branch
        # Divergence times = the stem age
        stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
        stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
        divergence_times_i <- stem_age_i
  
        # Set upper bound to drawn divergence age as the stem age
        upper_age_i <- stem_age_i
        
        # Compute sampling fraction
        ntips_present <- 1
        ntips_missing <- length(missing_taxa_list_i$Taxa_name)
        rho <- ntips_present / (ntips_present + ntips_missing)
        
        ## Fit a BD on a substitute clade of at least 3 tips, encompassing the target grafting tip
        
        # Retrieve MRCA node of the extended substitute clade
        substitute_clade2_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i), 1]
        substitute_clade3_node_i <- phy$edge[(phy$edge[,2] == substitute_clade2_node_i), 1]
        
        # Extract the substitute tree
        substitute_clade_tree_i <- ape::extract.clade(phy = phy, node = substitute_clade3_node_i)
        
        # plot(substitute_clade_tree_i)
        
        # Reidentify number of terminals in the BD clade
        nb_desc_tips_substitute_clade <- length(substitute_clade_tree_i$tip.label)
        # Adjust sampling fraction
        rho_BD_clade <- nb_desc_tips_substitute_clade / (nb_desc_tips_substitute_clade + ntips_missing)
        
        # Estimate BD rates from the substitute tree including at least 3 tips
        substitute_clade_tree_i <- force.ultrametric(tree = substitute_clade_tree_i, method = "nnls", message = F)
        substitute_clade_BD_i <- phytools::fit.bd(tree = substitute_clade_tree_i, rho = rho_BD_clade)
        
        ## Draw divergence times
        
        # Draw divergence ages of missing taxa
        # ?TreeSim::corsim
        new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
        divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)
        
        # divergence_times_missing_taxa_i
        
      }
      
      ## Graft missing taxa using the divergence times
      
      # Graft new taxa from oldest to youngest to allow the possibility to graft missing taxa on missing taxa branches!
      for (j in 1:length(missing_taxa_list_i$Taxa_name))
      {
        divergence_age <- divergence_times_missing_taxa_i[j]
        tip_label <- missing_taxa_list_i$Taxa_name[j]
        MRCA_tips <- missing_taxa_list_i$MRCA_tips[[j]]
        excluding_subclades_tips_list <- missing_taxa_list_i$Excluding_subclades[[j]]
        
        phy <- grafting_from_divergence_age(phy = phy, # Backbone phylogeny
                                         divergence_age = divergence_age, # Divergence age
                                         tip_label = tip_label, # Missing taxa label
                                         MRCA_tips = MRCA_tips, # To provide the pair of tips defining the MRCA node of the total clade within which to graft the missing taxa
                                         excluding_subclades_tips_list = excluding_subclades_tips_list, # List of pairs defining subclades to exclude from grafting options
                                         seed = NULL) # Set seed outside of the loop
      }
        
      # Update number of missing taxa / stem clades grafted
      nb_taxa_grafted <- nb_taxa_grafted + length(missing_taxa_list_i$Taxa_name)
      nb_stem_clades_grafted <- nb_stem_clades_grafted + 1
      
      # Print progress every stem groups
      if (verbose & (nb_stem_clades_grafted %% 10 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted within clades n°",nb_stem_clades_grafted,"/",length(groups_ID)," - Taxa n°",nb_taxa_grafted,"/",length(missing_taxa_list$Taxa_name),"\n"))
      }
    }
  }
  
  
  # Force ultrametry to avoid issue with bad rounding
  phy <- force.ultrametric(tree = phy, method = "nnls", message = F)
  
  # Return new phylogeny
  return(phy)
}

### 4.3/ Grafting missing taxa from divergence age #####

grafting_from_divergence_age <- function (phy, # Backbone phylogeny
                                       divergence_age, # Divergence age
                                       tip_label, # Missing taxa label
                                       MRCA_tips, # To provide the pair of tips defining the MRCA node of the total clade within which to graft the missing taxa. Can also provide a single tip used to graft the missing taxa on.
                                       excluding_subclades_tips_list, # List of pairs defining subclades to exclude from grafting options
                                       seed = 1234) 
  
{
  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }
  
  # Extract root age
  root_age <- max(phytools::nodeHeights(tree = phy)[, 2])
  
  # Subset tipward nodes of branches crossing the target divergence age
  branch_times_df <- as.data.frame(root_age - phytools::nodeHeights(tree = phy))
  branch_times_df$rootward_age_match <- branch_times_df[,1] >= divergence_age
  branch_times_df$tipward_age_match <- branch_times_df[,2] <= divergence_age
  branch_times_df$compatible <- branch_times_df$rootward_age_match & branch_times_df$tipward_age_match
  compatible_branches_ID <- which(branch_times_df$compatible)
  compatible_nodes_ID <- phy$edge[compatible_branches_ID, 2]
  
  # Case with two tips used to define MRCA of grafting clade
  if (length(MRCA_tips) == 2)
  {
    # Get MRCA node
    MRCA_node <- ape::getMRCA(phy = phy, tip = MRCA_tips)
    # Extract descendant nodes from MRCA
    descendant_nodes <- phytools::getDescendants(tree = phy, node = MRCA_node)
    # Include node of the stem (= MRCA node)
    stem_clade_nodes <- unique(c(descendant_nodes, MRCA_node))
  } else {
    # Case with a single grafting tip
    MRCA_node <- which(phy$tip.label == MRCA_tips)
    # Stem clade nodes = stem node of the unique grafting tip
    stem_clade_nodes <- MRCA_node
  }

  # Subset nodes of clade to use for grafting as the intersect between the nodes within the stem clades, and the nodes of branching crossing the divergence age
  valid_nodes_ID <- intersect(compatible_nodes_ID, stem_clade_nodes)
  
  # Exclude nodes from select options based on list of pairs defining excluding subclades
  if (length(excluding_subclades_tips_list) > 0)
  {
    # Initiate vector of forbidden nodes
    excluded_nodes_ID <- c()
    # Loop per exclusive subclade
    for (j in 1:length(excluding_subclades_tips_list))
    {
      # j <- 1
      
      # Extract exclusive subclade tips
      excluding_subclade_tips_j <- excluding_subclades_tips_list[[j]]
      # Retrieve exclusive subclade MRCA node
      exclusive_subclade_MRCA_node_j <- ape::getMRCA(phy = phy, tip = excluding_subclade_tips_j)
      # Extract descendant nodes from exclusive subclade MRCA
      excluded_nodes_ID_i <- phytools::getDescendants(tree = phy, node = exclusive_subclade_MRCA_node_j)
      # Store in vector of forbidden branches
      excluded_nodes_ID <- c(excluded_nodes_ID, excluded_nodes_ID_i)
    }
    
    # Remove excluded nodes from the draw
    valid_nodes_ID <- setdiff(valid_nodes_ID, excluded_nodes_ID)
  }
  
  # Select randomly the branch/node among the valid options
  selection_index <- sample(x = 1:length(valid_nodes_ID), size = 1)
  selected_node_ID <- valid_nodes_ID[selection_index]
  
  # Compute position needed for the tree to remain ultrametric such as position + selected node age = divergence age
  position <- divergence_age - (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
  
  # Bind new taxa in the tree
  phy <- phytools::bind.tip(tree = phy, tip.label = tip_label, edge.length = divergence_age, where = selected_node_ID, position = position)
  
  # Force ultrametry to avoid issue with bad rounding
  phy <- force.ultrametric(tree = phy, method = "nnls", message = F)
  
  # Return phylogeny with the new grafted taxa
  return(phy)
}


##### 5/ Run one random grafting to create main imputed phylogeny #####

### 5.1/ On reduced toyset phylogeny ####

phy <- sim.bdtree(b = 1, d = 0, n = 10)
plot(phy)

test_phy <- phy

plot(test_phy)
nodelabels()
edgelabels()

MRCA_node <- 18
MRCA_node <- 5
MRCA_node <- 12
exclusive_subclades_MRCA_nodes <- c(13, 16)

missing_taxa_list <- list(Taxa_name = c("A", "B", "C"),
                          MRCA_tips = list(c("s3", "s10"), c("s3"), c("s3", "s8")),
                          Exclusive_subclades = list(list(c("s3", "s8"), c("s9", "s10")), NULL, list(c("s7", "s3"))))



grafted_phy <- grafting_missing_taxa(phy = test_phy, # Backbone phylogeny
                                     missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                     include_stem = T, # Include or not the stem of the total clade when grafting?
                                     use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                     method = "uniform", # Either "uniform", "BD_taxa", or "BD_clades"
                                     verbose = T,
                                     seed = 1234) # Set the seed

plot(grafted_phy)


### 5.2/ On Ectomyrmex-Cryptopone phylogeny ####

# Case study on missing Ectomomyrmex in the Ectomomyrmex-Cryptopone clade

Ecto_Crypto_MRCA_node <- getMRCA(phy = Ponerinae_phylogeny_789t_calibrated, tip = c("Ectomomyrmex_lobocarenus_EX2957_CASENT0235335", "Cryptopone_my07_EX2746_CASENT0700813"))
Ecto_Crypto_tree <- ape::extract.clade(phy = Ponerinae_phylogeny_789t_calibrated, node = Ecto_Crypto_MRCA_node)
Ecto_MRCA_node <- getMRCA(phy = Ponerinae_phylogeny_789t_calibrated, tip = c("Ectomomyrmex_lobocarenus_EX2957_CASENT0235335", "Ectomomyrmex_ph02_EX2959_CASENT0266733"))
Ecto_tree <- ape::extract.clade(phy = Ponerinae_phylogeny_789t_calibrated, node = Ecto_MRCA_node)

Ecto_missing_taxa_ID <- which(str_detect(string = Missing_taxa_list$Taxa_name, pattern = "Ectomomyrmex"))
Ecto_missing_taxa_names <- Missing_taxa_list$Taxa_name[str_detect(string = Missing_taxa_list$Taxa_name, pattern = "Ectomomyrmex")]

Ecto_missing_taxa_list <- list()
Ecto_missing_taxa_list$Taxa_name <- Missing_taxa_list$Taxa_name[Ecto_missing_taxa_ID]
Ecto_missing_taxa_list$MRCA_tips <- Missing_taxa_list$MRCA_tips[Ecto_missing_taxa_names]
Ecto_missing_taxa_list$Excluding_subclades <- Missing_taxa_list$Excluding_subclades[Ecto_missing_taxa_names]

## Adjust MRCA data to test fro different cases
Ecto_missing_taxa_list_test <- Ecto_missing_taxa_list

# Test with stem subclade of only 3 tips
Ecto_missing_taxa_list_test$MRCA_tips$Ectomomyrmex_insulanus[1] <- "Ectomomyrmex_ruficornis_EX2694_JTLC000008604"

# Test with stem subclade of only 2 tips
Ecto_missing_taxa_list_test$MRCA_tips$Ectomomyrmex_insulanus[1] <- "Ectomomyrmex_simillimus_EX3022_CASENT0650200"

# Test with stem subclade of only 1 tip
Ecto_missing_taxa_list_test$MRCA_tips$Ectomomyrmex_insulanus <- Ecto_missing_taxa_list$MRCA_tips$Ectomomyrmex_insulanus[1]

missing_taxa_list = Ecto_missing_taxa_list_test

## Run grafting
phy_test_BD_clades <- grafting_missing_taxa(phy = Ecto_Crypto_tree, # Backbone phylogeny
                                            missing_taxa_list = Ecto_missing_taxa_list_test, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                            include_stem = T, # Include or not the stem of the total clade when grafting?
                                            use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                            method = "BD_clades", # Either "uniform", "BD_taxa", or "BD_clades"
                                            verbose = T, # To print progress every 100 taxa or every 10 stem clades
                                            seed = 1234) # Set the seed

## Plot results
par(mfrow = c(2,2))
plot(phy)

# Adjust tip color scheme
tip_col <- c("black", "red")[as.numeric(phy_test_unif$tip.label %in% Ecto_missing_taxa_list_test$Taxa_name) + 1]
names(tip_col) <- phy_test_unif$tip.label
plot(phy_test_unif, tip.color = tip_col)

# Adjust tip color scheme
tip_col <- c("black", "red")[as.numeric(phy_test_BD_taxa$tip.label %in% Ecto_missing_taxa_list_test$Taxa_name) + 1]
names(tip_col) <- phy_test_BD_taxa$tip.label
plot(phy_test_BD_taxa, tip.color = tip_col)

# Adjust tip color scheme
tip_col <- c("black", "red")[as.numeric(phy_test_BD_clades$tip.label %in% Ecto_missing_taxa_list_test$Taxa_name) + 1]
names(tip_col) <- phy_test_BD_clades$tip.label
plot(phy_test_BD_clades, tip.color = tip_col)


### 5.3/ On roughly calibrated Ponerinae backbone phylogeny ####

# Load backbone phylogeny, calibrated with rough ages
Ponerinae_phylogeny_789t_calibrated <- readRDS(file = "./input_data/Phylogenies/Ponerinae_phylogeny_789t_calibrated.rds")

# Load Missing_taxa_list
Missing_taxa_list <- readRDS(file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")

Ponerinae_phylogeny_1534t <- grafting_missing_taxa(phy = Ponerinae_phylogeny_789t_calibrated, # Backbone phylogeny
                                                   missing_taxa_list = Missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                                   include_stem = T, # Include or not the stem of the total clade when grafting?
                                                   use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                                   method = "BD_clades", # Either "uniform", "BD_taxa", or "BD_clades"
                                                   verbose = T, # To print progress every 100 taxa or every 10 stem clades
                                                   seed = 1234) # Set the seed
# May get warnings from fit.bd when BD parameter optimization fails at first attempts

length(Ponerinae_phylogeny_1534t$tip.label)

# Save imputed phylogeny
saveRDS(object = Ponerinae_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t.rds")

## Plot result

# Detect grafted tips on the phylogeny
missing_taxa_ID <- which(Ponerinae_phylogeny_1534t$tip.label %in% Missing_taxa_list$Taxa_name)
# Adjust tip color scheme
tip_col <- c("black", "red")[as.numeric(Ponerinae_phylogeny_1534t$tip.label %in% Missing_taxa_list$Taxa_name) + 1]
names(tip_col) <- Ponerinae_phylogeny_1534t$tip.label

pdf(file = "./outputs/Grafting_missing_taxa/grafted_phy_BD_clades.pdf", width = 20, height = 200)
par(mar = c(5, 5, 5, 5))
plot(Ponerinae_phylogeny_1534t, tip.color = tip_col)
dev.off()

### 5.4/ On MCC phylogeny calibrated with TreePL ####

# Load MCC phylogeny calibrated with TreePL
Ponerinae_MCC_phylogeny_789t_treedata <- readRDS(file = "./input_data/Phylogenies/Ponerinae_MCC_phylogeny_789t_treedata.rds")

# Load Missing_taxa_list
Missing_taxa_list <- readRDS(file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")

# Graft missing taxa
Ponerinae_MCC_phylogeny_1534t <- grafting_missing_taxa(phy = Ponerinae_MCC_phylogeny_789t_treedata@phylo, # Backbone phylogeny
                                                       missing_taxa_list = Missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                                       include_stem = T, # Include or not the stem of the total clade when grafting?
                                                       use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                                       method = "BD_clades", # Either "uniform", "BD_taxa", or "BD_clades"
                                                       verbose = T, # To print progress every 100 taxa or every 10 stem clades
                                                       seed = 1234) # Set the seed
# May get warnings from fit.bd when BD parameter optimization fails at first attempts

Ponerinae_MCC_phylogeny_1534t
length(Ponerinae_MCC_phylogeny_1534t$tip.label)

# Save imputed phylogeny
saveRDS(object = Ponerinae_MCC_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_MCC_1534t.rds")
# Save NEXUS format
write.tree(phy = Ponerinae_MCC_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_MCC_1534t.tree")

## Manual tweaking of volatile taxa to avoid weird placements

# None required?

## Initiate new treedata object for the grafted tree
Ponerinae_MCC_phylogeny_1534t_treedata <- Ponerinae_MCC_phylogeny_789t_treedata
Ponerinae_MCC_phylogeny_1534t_treedata@phylo <- Ponerinae_MCC_phylogeny_1534t # Replace phylogeny
Ponerinae_MCC_phylogeny_1534t_treedata@data$backbone_node <- Ponerinae_MCC_phylogeny_1534t_treedata@data$node # Save backbone node ID
Ponerinae_MCC_phylogeny_1534t_treedata@data$node <- NA # Initiate new node ID

## Add metadata for backbone nodes

# Need to find correspondence between previous nodes and current nodes using matching descendant tips and age

# Get node ages for backbone tree
root_age_789t <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_789t_treedata@phylo))
Ponerinae_MCC_phylogeny_789t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_MCC_phylogeny_789t_treedata@phylo) + root_age_789t, 5))
names(Ponerinae_MCC_phylogeny_789t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_MCC_phylogeny_789t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_MCC_phylogeny_789t_edges_ages_df)
Ponerinae_MCC_phylogeny_789t_edges_ages_df <- Ponerinae_MCC_phylogeny_789t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_MCC_phylogeny_789t_edges_ages_df <- cbind(Ponerinae_MCC_phylogeny_789t_edges_ages_df, Ponerinae_MCC_phylogeny_789t_treedata@phylo$edge)
names(Ponerinae_MCC_phylogeny_789t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_MCC_phylogeny_789t_nodes_ages_df <- Ponerinae_MCC_phylogeny_789t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_MCC_phylogeny_789t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label),
                                                         node_age = 0)
Ponerinae_MCC_phylogeny_789t_nodes_ages_df <- rbind(Ponerinae_MCC_phylogeny_789t_tips_ages_df, Ponerinae_MCC_phylogeny_789t_nodes_ages_df)

# Get node ages for grafted tree
root_age_1534t <- max(phytools::nodeHeights(Ponerinae_MCC_phylogeny_1534t))
Ponerinae_MCC_phylogeny_1534t_edges_ages_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = Ponerinae_MCC_phylogeny_1534t) + root_age_1534t, 5))
names(Ponerinae_MCC_phylogeny_1534t_edges_ages_df) <- c("rootward_age", "tipward_age")
Ponerinae_MCC_phylogeny_1534t_edges_ages_df$edge_ID <- 1:nrow(Ponerinae_MCC_phylogeny_1534t_edges_ages_df)
Ponerinae_MCC_phylogeny_1534t_edges_ages_df <- Ponerinae_MCC_phylogeny_1534t_edges_ages_df[, c("edge_ID", "rootward_age", "tipward_age")]
Ponerinae_MCC_phylogeny_1534t_edges_ages_df <- cbind(Ponerinae_MCC_phylogeny_1534t_edges_ages_df, Ponerinae_MCC_phylogeny_1534t$edge)
names(Ponerinae_MCC_phylogeny_1534t_edges_ages_df) <- c("edge_ID", "rootward_age", "tipward_age", "rootward_node_ID", "tipward_node_ID")
Ponerinae_MCC_phylogeny_1534t_nodes_ages_df <- Ponerinae_MCC_phylogeny_1534t_edges_ages_df %>% 
  select(rootward_node_ID, rootward_age) %>% 
  arrange(rootward_node_ID) %>% 
  distinct(rootward_node_ID, rootward_age) %>%
  rename(node_ID = rootward_node_ID,
         node_age = rootward_age)
Ponerinae_MCC_phylogeny_1534t_tips_ages_df <- data.frame(node_ID = 1:length(Ponerinae_MCC_phylogeny_1534t$tip.label),
                                                         node_age = 0)
Ponerinae_MCC_phylogeny_1534t_nodes_ages_df <- rbind(Ponerinae_MCC_phylogeny_1534t_tips_ages_df, Ponerinae_MCC_phylogeny_1534t_nodes_ages_df)

# Extract descendants for all nodes in backbone tree
Ponerinae_MCC_phylogeny_789t_descendants_list <- list()
nb_nodes_789t <- dim(Ponerinae_MCC_phylogeny_789t_treedata@data)[1]
for (i in 1:nb_nodes_789t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_MCC_phylogeny_789t_treedata@phylo, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label[all_descendants_i]
  
  # Store in list
  Ponerinae_MCC_phylogeny_789t_descendants_list[[i]] <- all_descendants_i
}

# Extract descendants for all nodes in grafted tree
Ponerinae_MCC_phylogeny_1534t_descendants_list <- list()
nb_nodes_1534t <- Ponerinae_MCC_phylogeny_1534t$Nnode + length(Ponerinae_MCC_phylogeny_1534t$tip.label)
for (i in 1:nb_nodes_1534t)
{
  # Extract nodes and tips
  all_descendants_i <- getDescendants(Ponerinae_MCC_phylogeny_1534t, node = i)
  # Keep only tips
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_MCC_phylogeny_1534t$tip.label)] 
  # Convert to labels
  all_descendants_i <- Ponerinae_MCC_phylogeny_1534t$tip.label[all_descendants_i]
  # Keep only tips from the backbone
  all_descendants_i <- all_descendants_i[all_descendants_i %in% Ponerinae_MCC_phylogeny_789t_treedata@phylo$tip.label]
  # Reorder alphabetically
  all_descendants_i <- all_descendants_i[order(all_descendants_i)]
  
  # Store in list
  Ponerinae_MCC_phylogeny_1534t_descendants_list[[i]] <- all_descendants_i
}

# Find corresponding node in the grafted tree based on list of descending tips and node age
for (i in 1:nb_nodes_789t)
{
  # Extract descending tips from backbone node
  all_descendants_backbone_i <- Ponerinae_MCC_phylogeny_789t_descendants_list[[i]]
  # Reorder alphabetically
  all_descendants_backbone_i <- all_descendants_backbone_i[order(all_descendants_backbone_i)]
  # Find the matches among list of descendants from grafted tree
  descendant_match_test_i <- unlist(lapply(X = Ponerinae_MCC_phylogeny_1534t_descendants_list, FUN = function (x) { identical(x, all_descendants_backbone_i) }))
  descendant_match_nodes_ID <- which(descendant_match_test_i)

  # Extract node age in the backbone tree
  focal_node_age <- Ponerinae_MCC_phylogeny_789t_nodes_ages_df$node_age[i]
  # Extract node ages of the matching nodes in the grafted tree
  descendant_match_nodes_ages <- Ponerinae_MCC_phylogeny_1534t_nodes_ages_df$node_age[descendant_match_nodes_ID]
  # Find the matching node age
  age_matching_node_ID <- descendant_match_nodes_ID[abs(focal_node_age - descendant_match_nodes_ages) < 0.001]
  
  # Store matching ID in treedata object
  Ponerinae_MCC_phylogeny_1534t_treedata@data$node[i] <- age_matching_node_ID
}

# Create template for all nodes (height = age, length of parental edge, posterior (tips = NA, non-tips = 1))
Ponerinae_MCC_phylogeny_1534t_nodes_data_df <- Ponerinae_MCC_phylogeny_1534t_nodes_ages_df %>% 
  rename(node = node_ID,
         height = node_age)
Ponerinae_MCC_phylogeny_1534t_nodes_data_df$posterior <- NA
Ponerinae_MCC_phylogeny_1534t_nodes_data_df$posterior[Ponerinae_MCC_phylogeny_1534t_nodes_data_df$node > length(Ponerinae_MCC_phylogeny_1534t$tip.label)] <- 1
Ponerinae_MCC_phylogeny_1534t_nodes_data_df$length <- Ponerinae_MCC_phylogeny_1534t$edge.length[match(x = Ponerinae_MCC_phylogeny_1534t_nodes_data_df$node, table = Ponerinae_MCC_phylogeny_1534t$edge[,2])]
  
# Join backbone data with template data
Ponerinae_MCC_phylogeny_1534t_treedata@data <- left_join(Ponerinae_MCC_phylogeny_1534t_nodes_data_df, Ponerinae_MCC_phylogeny_1534t_treedata@data, by = "node") %>% 
  rename(height = height.x) %>% # Keep the height/age as measured in the grafted tree
  select(-height.y) %>% # Remove height as measured in the backbone tree (should be equal modulo numerical precision)
  rename(posterior = posterior.x) %>% # Keep info on non-tips from the grafted tree
  select(-posterior.y) %>% # Remove info on non-tips from the backbone tree as some nodes are missing
  rename(length = length.x) %>% # Keep parental branch length from the grafted tree as the current length
  rename(length_backbone = length.y) %>% # Record parental branch length from the backbone tree (may be different due to grafting)
  rename(length_median_backbone = length_median,
         length_0.95_HPD_backbone = length_0.95_HPD,
         length_range_backbone = length_range,
         height_median_backbone = height_median,
         height_0.95_HPD_backbone = height_0.95_HPD,
         height_range_backbone = height_range) %>% 
  mutate(new_node = is.na(backbone_node)) %>%
  select(node, backbone_node, posterior, height, height_median_backbone, height_0.95_HPD_backbone, height_range_backbone, length, length_backbone, length_median_backbone, length_0.95_HPD_backbone, length_range_backbone, new_node) %>% 
  as_tibble()

# Record nodes with new parental edge
all_missing_taxa_ID <- which(Ponerinae_MCC_phylogeny_1534t_treedata@data$new_node)
all_missing_taxa_ID <- all_missing_taxa_ID[all_missing_taxa_ID <= length(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label)]
Ponerinae_MCC_phylogeny_1534t_treedata@data$new_parental_edge <- FALSE

for (i in 1:nrow(Ponerinae_MCC_phylogeny_1534t_treedata@data))
{
  # i <- 1
  
  # Extract node ID
  node_ID_i <- Ponerinae_MCC_phylogeny_1534t_treedata@data$node[i]
  
  # Extract all current descendants
  all_descendants_i <- phytools::getDescendants(tree = Ponerinae_MCC_phylogeny_1534t_treedata@phylo, node = node_ID_i)
  all_descendants_i <- all_descendants_i[all_descendants_i <= length(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label)]
  
  # Check if all current descendants are missing taxa
  missing_i <- all(all_descendants_i %in% all_missing_taxa_ID)
  
  # Record status
  Ponerinae_MCC_phylogeny_1534t_treedata@data$new_parental_edge[i] <- missing_i
}
table(Ponerinae_MCC_phylogeny_1534t_treedata@data$new_parental_edge)

View(Ponerinae_MCC_phylogeny_1534t_treedata@data)
View(Ponerinae_MCC_phylogeny_789t_treedata@data)

## Save grafted MCC phylogeny with node metadata
saveRDS(object = Ponerinae_MCC_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")

## Plot result

# Detect grafted tips on the phylogeny
missing_taxa_ID <- which(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label %in% Missing_taxa_list$Taxa_name)
# Adjust tip color scheme
tip_col <- c("black", "red")[as.numeric(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label %in% Missing_taxa_list$Taxa_name) + 1]
names(tip_col) <- Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label

pdf(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_BD_clades.pdf", width = 20, height = 200)
par(mar = c(5, 5, 5, 5))
plot(Ponerinae_MCC_phylogeny_1534t_treedata@phylo, tip.color = tip_col)
dev.off()

## Plot result with node/tip labels

# Use ggtree
pdf(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_node_labels.pdf", height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_MCC_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  coord_cartesian(clip = 'off') + 
  
  geom_tree(mapping = aes(colour = new_parental_edge),
            # color = branches_color,
            linewidth = 0.5) +
  
  geom_tiplab(# mapping = aes(colour = missing_node),
    color = tip_col,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Add ID labels to nodes
  geom_label(aes(label = node),
             fill = 'steelblue', size = 1) +
  
  # # Add ID labels to edges
  # geom_label(aes(x = branch, label = parental_edge_ID),
  #            fill = 'lightgreen', size = 1) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.002,0.12,0.002,0), "npc"), # trbl
        legend.position = c(1.09, 0.5),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )

print(Ponerinae_phylogeny_plot)

dev.off()

## Plot result with 95% HPD interval on backbone nodes 

# Use ggtree
pdf(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_95HPD.pdf", height = 100, width = 30)

Ponerinae_phylogeny_plot <- ggtree(Ponerinae_MCC_phylogeny_1534t_treedata,
                                   layout = "rectangular") +
  
  coord_cartesian(clip = 'off') + 
  
  geom_tree(mapping = aes(colour = new_parental_edge),
            # color = branches_color,
            linewidth = 0.5) +
  
  geom_tiplab(# mapping = aes(colour = missing_node),
    color = tip_col,
    align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
    linetype = 0, 
    size = 2, 
    offset = 0.3,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
    hjust = 0) +    # To center, left-align, right-aligned text
  
  # Adjust tip/edge color scheme
  scale_colour_manual("Taxa source", breaks = c(F, T), values = c("black", "red"), labels = c("UCE data", "Grafted")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 8))) +
  
  # Add HPD 95% intervals on backbone nodes
  geom_range(range = "height_0.95_HPD_backbone", center = "height",
             color = "steelblue", size = 2, alpha = 0.5) +
  
  # Adjust legend for taxa source
  theme(plot.margin = unit(c(0.002,0.12,0.002,0), "npc"), # trbl
        legend.position = c(1.09, 0.5),
        legend.title = element_text(size  = 25, margin = margin(b = 10)), 
        legend.text = element_text(size = 20, margin = margin(b = 5, t = 5, l = 20)),
        # legend.key.size = unit(1.8, "lines"),
        # legend.spacing.y = unit(1.0, "lines")
  )

print(Ponerinae_phylogeny_plot)

dev.off()


##### 6/ Run 1000 random grafting to create robustness set of imputed phylogenies ####

### 6.1/ On the same backbone phylogeny ####

N_phylo <- 2

imputed_phylogenies_1000 <- list()

# Set seed for reproducibility outside of the loop!
set.seed(seed = 1234)

# Loop per imputed phylogeny to generate
for (i in 1:N_phylo)
{
  # Graft randomly missing taxa on the phylogeny
  imputed_phylo_i <- grafting_missing_taxa(phy = Ponerinae_MCC_phylogeny_1534t_treedata@phylo, # Backbone phylogeny
                                           missing_taxa_list = Missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                           include_stem = T, # Include or not the stem of the total clade when grafting?
                                           use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                           method = "BD_clades", # Either "uniform", "BD_taxa", or "BD_clades"
                                           verbose = T, # To print progress every 100 taxa or every 10 stem clades
                                           seed = NULL) # DO NOT set the seed inside the function, but outside of the loop, otherwise all phylogenies will be identical!

  # Store output
  old_names <- names(imputed_phylogenies_1000)
  imputed_phylogenies_1000 <- append(x = imputed_phylogenies_1000, values = list(imputed_phylo_i))
  names(c(imputed_phylogenies_1000, paste0("Phylo_", i)))
  
  # Print progress every 1 phylogenies
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - Missing taxa grafted on phylogeny n°",i,"/",N_phylo,"\n"))
  }
}

plot(imputed_phylogenies_1000[[1]])
plot(imputed_phylogenies_1000[[2]])

# Descendants nodes should (not) be similar if (not) using the same seed
test1 <- phytools::getDescendants(tree = imputed_phylogenies_1000[[1]], node = 2000)
test1 <- test1[test1 <= length(imputed_phylogenies_1000[[1]]$tip.label)]
  
test2 <- phytools::getDescendants(tree = imputed_phylogenies_1000[[2]], node = 2000)
test2 <- test2[test2 <= length(imputed_phylogenies_1000[[2]]$tip.label)]

### 6.2/ On the 1000 posterior backbone phylogenies ####

# Load the 1000 posterior backbone phylogenies from TreePL
Ponerinae_all_posteriors_phylogeny_789t <- readRDS(file = "./input_data/Phylogenies/Ponerinae_all_posteriors_phylogeny_789t.rds")

N_phylo <- length(Ponerinae_all_posteriors_phylogeny_789t)

Ponerinae_all_posteriors_phylogeny_1534t <- list()

# Set seed for reproducibility outside of the loop!
set.seed(seed = 1234)

# Loop per imputed phylogeny to generate
for (i in 1:N_phylo)
# for (i in 1:2)
{
  # Graft randomly missing taxa on the phylogeny
  imputed_phylo_i <- grafting_missing_taxa(phy = Ponerinae_all_posteriors_phylogeny_789t[[i]], # Backbone phylogeny n°i
                                           missing_taxa_list = Missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                           include_stem = T, # Include or not the stem of the total clade when grafting?
                                           use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                           method = "BD_clades", # Either "uniform", "BD_taxa", or "BD_clades"
                                           verbose = T, # To print progress every 100 taxa or every 10 stem clades
                                           seed = NULL) # DO NOT set the seed inside the function, but outside of the loop, otherwise all phylogenies will be identical (if using the same backbone)!
  
  # Store output
  old_names <- names(Ponerinae_all_posteriors_phylogeny_1534t)
  Ponerinae_all_posteriors_phylogeny_1534t <- append(x = Ponerinae_all_posteriors_phylogeny_1534t, values = list(imputed_phylo_i))
  names(c(Ponerinae_all_posteriors_phylogeny_1534t, paste0("Phylo_", i)))
  
  # Print progress every 1 phylogenies
  if (i %% 1 == 0)
  {
    # Save imputed posterior phylogenies
    saveRDS(object = Ponerinae_all_posteriors_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.rds")
    
    cat(paste0(Sys.time(), " - Missing taxa grafted on phylogeny n°",i,"/",N_phylo,"\n"))
  }
}

# Assign ape multiPhylo class
class(Ponerinae_all_posteriors_phylogeny_1534t) <- "multiPhylo"
table(unlist(lapply(X = Ponerinae_all_posteriors_phylogeny_1534t, FUN = function (x) { length(x$tip.label) } )))

# Save imputed posterior phylogenies
saveRDS(object = Ponerinae_all_posteriors_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.rds")
# Save NEXUS format
write.tree(phy = Ponerinae_all_posteriors_phylogeny_1534t, file = "./outputs/Grafting_missing_taxa/Ponerinae_all_posteriors_phylogeny_1534t.tree")

plot(Ponerinae_all_posteriors_phylogeny_1534t[[1]])
plot(Ponerinae_all_posteriors_phylogeny_1534t[[200]])

# Descendants nodes should (not) be similar if (not) using the same seed
test1 <- phytools::getDescendants(tree = Ponerinae_all_posteriors_phylogeny_1534t[[1]], node = 2000)
test1 <- test1[test1 <= length(Ponerinae_all_posteriors_phylogeny_1534t[[1]]$tip.label)]
test1

test2 <- phytools::getDescendants(tree = Ponerinae_all_posteriors_phylogeny_1534t[[200]], node = 2000)
test2 <- test2[test2 <= length(Ponerinae_all_posteriors_phylogeny_1534t[[200]]$tip.label)]
test2

### 6.3/ Plot all grafted phylogenies (base plot version) ####

for (i in 1:N_phylo)
# for (i in 1:2)
{
  # Detect grafted tips on the phylogeny
  missing_taxa_ID <- which(Ponerinae_all_posteriors_phylogeny_1534t[[i]]$tip.label %in% Missing_taxa_list$Taxa_name)
  # Adjust tip color scheme
  tip_col <- c("black", "red")[as.numeric(Ponerinae_all_posteriors_phylogeny_1534t[[i]]$tip.label %in% Missing_taxa_list$Taxa_name) + 1]
  names(tip_col) <- Ponerinae_all_posteriors_phylogeny_1534t[[i]]$tip.label
  
  pdf(file = paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t/Ponerinae_phylogeny_1534t_posterior_",i,".pdf"), width = 20, height = 200)
  par(mar = c(5, 5, 5, 5))
  plot(Ponerinae_all_posteriors_phylogeny_1534t[[i]], tip.color = tip_col)
  
  dev.off()
}

## Aggregate all posterior phylo in a unique PDF

Ponerinae_all_posteriors_phylogeny_1534t_path <- list.files(path = "./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t/", pattern = "Ponerinae_phylogeny_1534t_posterior_", full.names = T)
Ponerinae_all_posteriors_phylogeny_1534t_path <- Ponerinae_all_posteriors_phylogeny_1534t_path[str_detect(string = Ponerinae_all_posteriors_phylogeny_1534t_path, pattern = ".pdf")]
nb_phylos <- length(Ponerinae_all_posteriors_phylogeny_1534t_path)

# Reorder in numerical rather than alphabetic order
prefix <- str_remove(string = Ponerinae_all_posteriors_phylogeny_1534t_path[1], pattern = "_\\d+.pdf")
Ponerinae_all_posteriors_phylogeny_1534t_path <- paste0(prefix,"_",1:nb_phylos,".pdf")

qpdf::pdf_combine(input = Ponerinae_all_posteriors_phylogeny_1534t_path, output = paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t_",nb_phylos,"_phylos.pdf"))

## Do it by smaller batch to avoid saturating RAM

batch_size <- 100
start_i <- seq(from = 1, to = nb_phylos, by = batch_size)

for (i in seq_along(start_i))
{
  subset_indices_i <- start_i[i]:(start_i[i] + batch_size - 1)
  
  all_grafted_phylos_path_i <- Ponerinae_all_posteriors_phylogeny_1534t_path[subset_indices_i]
  nb_maps_i <- length(all_grafted_phylos_path_i)
  subset_label <- paste0(subset_indices_i[1],"_",subset_indices_i[nb_maps_i])
  
  qpdf::pdf_combine(input = all_grafted_phylos_path_i, output = paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t/All_posterior_phylo_1534t_batch_",subset_label,"_phylos.pdf"))
  
  # Print progress
  cat(paste0(Sys.time(), " - Grafted phylos aggregated from n°", subset_indices_i[1], " to ", subset_indices_i[nb_maps_i],"\n"))
}

## Aggregate batches in a single PDF

all_grafted_phylos_batches_path <- list.files(path = "./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t/", pattern = "All_posterior_phylo_1534t_batch_", full.names = T)

nb_batches <- length(all_grafted_phylos_batches_path)

# # Reorder in numerical rather than alphabetic order
# prefix <- str_remove(string = all_grafted_phylos_batches_path[1], pattern = "_\\d+.pdf")
# all_grafted_phylos_batches_path <- paste0(prefix,"_",1:nb_maps,".pdf")

qpdf::pdf_combine(input = all_grafted_phylos_batches_path, output = paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t_",nb_phylos,"_phylos.pdf"))


## Make a GIF

# Limit number of phylos to reduce GIF size and avoid crashing because of memory allocation limits
nb_phylos_in_GIF <- nb_phylos
nb_phylos_in_GIF <- 100

source("./functions/image_resize_and_write_gif.R")

# Set fps
# fps <- 1
fps <- 5
# fps <- 10

# Load pdf as image frames
pdf_pointer_posterior_phylos <- magick::image_read_pdf(path = paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t_",nb_phylos_in_GIF,"_phylos.pdf"),
                                                  pages = NULL, density = 75)
magick::image_info(pdf_pointer_posterior_phylos)

image_resize_and_write_gif(image = pdf_pointer_posterior_phylos,
                           path =  paste0("./outputs/Grafting_missing_taxa/All_posterior_phylo_1534t_",nb_phylos_in_GIF,"_phylos.gif"),
                           delay = 1/fps, # Time between frames in seconds
                           width = 300, height = 3000,
                           loop = FALSE,
                           progress = TRUE)



