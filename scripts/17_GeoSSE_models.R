##### Script 17: GeoSSE models  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Test for State-dependent Diversification dynamics between Old World and New World

###

### Inputs

# (Set of) Time-calibrated phylogeny(ies)
# Binary table of tip bioregion membership (including polymorphic states)

###

### Outputs

# Binary table for OW vs. NW tips data (including widespread state encompassing OW and NW)
# Fits of Geographic SSE models: CID-1, CID-2, CID-4, GeoSSE, GeoHiSSE 
  # Table of AICc model comparison
  # Parameter estimates plots
# Estimate residence times from ARE marginal likelihoods
  # Plot Residence times ~ states
# Schematic summary of model parameter estimates: lambda - mu per state, vicariance, range extension, range extirpation
  # Plot Igraph networks for each model + AICc

###

### Bonus in other scripts (?)

# For a more smooth model with many small cladogenetic changes, see CLaDS


# Clean environment
rm(list = ls())

library(ape)
library(phytools)
library(tidyverse)
library(diversitree)  # For time-dependent BD models
library(hisse) # For models with hidden-states
library(igraph)
library(xlsx)



##### 1/ Load input files ####

### 1.1/ Load the phylogeny ####

# Load phylogeny
# Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE.rds")

### 1.2/ Load binary tip data for bioregion membership  ####

# Load binary biogeographic tip data
Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

### 1.3/ Home-made function to convert turnover rates and extinction fractions in speciation and extinctions rates ####

convert_SSE_rates <- function (turnover, epsilon)
{
  # Convert speciation rates
  if (length(turnover) == 1) 
  {
    speciation <- turnover / (1 + epsilon)
  } else {
    speciation <- unlist(purrr::map2(.x = turnover, .y = epsilon, .f = function (d, eps) { d / (1 + eps) }))
  }
  names(speciation) <- str_replace(string = names(speciation), pattern = "turnover", replacement = "lambda_")
  
  # Convert extinction rates
  if (length(epsilon) == 1) 
  {
    extinction <- (epsilon * turnover) / (1 + epsilon)
  } else {
    extinction <- unlist(purrr::map2(.x = turnover, .y = epsilon, .f = function (d, eps) { (eps * d) / (1 + eps) }))
  }
  names(extinction) <- str_replace(string = names(extinction), pattern = "turnover", replacement = "mu_")
  
  SE_rates <- list(speciation = speciation, extinction = extinction)
  return(SE_rates)
}

convert_rates_from_GeoHiSSE_output <- function(GeoSSE_output, keep_initial_pars = F)
{
  k <- nrow(GeoSSE_output$trans.matrix) # Extract nb of states
  states_names <- colnames(GeoSSE_output$trans.matrix) # Extract state names
  states_names_no_parenthesis <- str_remove_all(string = states_names, pattern = "\\(|\\)")
  # nb_regimes <- sum(colSums(apply(X = matrix(LETTERS[1:10]), MARGIN = 1, FUN = str_detect, string = states_names)) > 0) # Extract number of hidden states/regimes
  nb_regimes <- k/3 # Extract number of hidden states/regimes
  
  pars <- GeoSSE_output$solution # Extract parameters
  # tau_pars <- pars[grep("tau",names(pars))][1:k] # Extract turnover parameters
  # ef_pars <- pars[grep("ef",names(pars))][1:(2*nb_regimes)] # Extract extinction fraction parameters
  
  all_tau_pars <- pars[grep("tau",names(pars))] # Extract turnover parameters
  tau_pars <- all_tau_pars[match(x = states_names_no_parenthesis, table = str_remove(string = names(all_tau_pars), pattern = "tau"))]
  all_ef_pars <- pars[grep("ef",names(pars))] # Extract extinction fraction parameters
  ef_pars <- all_ef_pars[match(x = states_names_no_parenthesis, table = str_remove(string = names(all_ef_pars), pattern = "ef"))]
  # Replace NA of 01 states with 0
  ef_pars[is.na(ef_pars)] <- 0
    
  # Compute rates
  lambda <- tau_pars / (1 + ef_pars) # Compute speciation rates
  mu <- tau_pars - lambda # Compute extinction rates
  
  # Build output matrix
  if (!keep_initial_pars)
  {
    matrix(data = c(lambda, mu),
           nrow = k, ncol = 2,
           dimnames = list(states_names, c("lambda","mu")))
  } else {
    matrix(data = c(lambda, mu, tau_pars, ef_pars),
           nrow = k, ncol = 4,
           dimnames = list(states_names, c("lambda","mu","turnover","extinct fraction")))
  }
}

### 1.4/ Home-made function to extract Q matrix of transition rates ####

extract_Q_matrix_from_GeoHiSSE_output <- function(GeoSSE_output)
{
  Q_mat <- GeoSSE_output$trans.matrix # Extract the design Q matrix
  states_names <- colnames(Q_mat) # Extract state names
  
  pars <- GeoSSE_output$solution # Extract parameters
  q_pars <- pars[grep("d",names(pars))] # Extract transition parameters
  pars_indices <- GeoSSE_output$index.par[grep("d",names(pars))] # Extract transition parameter indices
  pars_indices <- pars_indices - min(pars_indices) + 1 # Adjust so they match indices in the design Q matrix
  
  # Fill Q matrix
  Q_mat_pars <- Q_mat
  for (i in 1:nrow(Q_mat))
  {
    for (j in 1:ncol(Q_mat))
    {
      q_index <- Q_mat[i,j]
      if ((is.na(q_index)) | (q_index != 0))
      {
        Q_mat_pars[i,j] <- unique(q_pars[q_index == pars_indices])
      }
    }
  }
  
  return(Q_mat_pars)
}


##### 2/ Convert tip data to GeoSSE data ####

Taxa_bioregions_binary_table$New_World <- (Taxa_bioregions_binary_table$Neotropics + Taxa_bioregions_binary_table$Nearctic) > 0
Taxa_bioregions_binary_table$Old_World <- (Taxa_bioregions_binary_table$Afrotropics + Taxa_bioregions_binary_table$Australasia + Taxa_bioregions_binary_table$Indomalaya + Taxa_bioregions_binary_table$`Eastern Palearctic` + Taxa_bioregions_binary_table$`Western Palearctic`) > 0
Taxa_bioregions_binary_table$Both <- (Taxa_bioregions_binary_table$New_World + Taxa_bioregions_binary_table$Old_World) == 2

table(Taxa_bioregions_binary_table$New_World, Taxa_bioregions_binary_table$Old_World, Taxa_bioregions_binary_table$Both)

# Values for the areas need to be 0, 1, or 2, where 
  # 0 is the widespread area ’NW + OW’, but parameters labeled as 01  # 1 in ClaSSE # 01 in GeoHiSSE
  # 1 is endemic area ’NW’, but parameters labeled as 0               # 2 in ClaSSE # 0 in GeoHiSSE
  # 2 is endemic area ’OW’, but parameters labeled as 1               # 3 in ClaSSE # 1 in GeoHiSSE

OW_NW_GeoSSE_tip_data <- as.numeric(Taxa_bioregions_binary_table$New_World)
OW_NW_GeoSSE_tip_data[Taxa_bioregions_binary_table$Old_World] <- 2
OW_NW_GeoSSE_tip_data[Taxa_bioregions_binary_table$Both] <- 0

table(OW_NW_GeoSSE_tip_data)

# Convert to df
OW_NW_GeoSSE_tip_data_df <- data.frame(taxa = Taxa_bioregions_binary_table$Current_name,
                                       area = OW_NW_GeoSSE_tip_data)

# Save GeoSSE tip data for OW vs. NW
saveRDS(object = OW_NW_GeoSSE_tip_data_df, file = "./outputs/GeoSSE/OW_NW_GeoSSE_tip_data_df.rds")


##### 3/ Run models #####

?diversitree::make.geosse # Alternative for GeoSSE. Handle time-dependent dynamics, but does not handle Hidden states

?hisse::GeoHiSSE # Handle hidden states, but not time-dependent dynamics (only constant state-dependent rates)
?hisse::TransMatMakerGeoHiSSE # To set the design Q matrix of transition parameters, fixing transition rates nullity and equality (between states, including hidden states)

# Load custom version of the GeoHiSSE function that allows keeping track of progress
source("./functions/GeoHiSSE_with_print.R")

# To show correspondence of parameters with a diversitree::classe model
hisse::TranslateParsMakerGeoHiSSE(k = 0,
                                  add.extinction = T, # Make a distinction between local extirpation in area 
                                  add.jumps = T)
hisse::TranslateParsMakerGeoHiSSE(k = 1,
                                  add.extinction = T, # Make a distinction between local extirpation in area 
                                  add.jumps = T)
# ClaSSE model have transition parameters for every possible combination of parental and pair of descendant branches


## What are extirpation parameters in a GeoSSE model?
  # Anagenetic transitions from the large area 01 to one of the endemic areas (0 and 1)
  # Extirpation = q12 and q13 (= x1A, x0A) = parameter labels designate area where extirpation occurred
  # Must be different from local extinction in endemic area
  # Local extinction = mu2, mu3 (= x*0A, x*1A) = parameter labels designate area where extinction occurred
  # In practice, the Transition Q design matrix indicates that extirpation are simply not allowed: q12 and q13 set to 0, thus mu2, mu3 <=> x1A, x0A

## What are jump dispersals in a GeoSSE model?
  # Anagenetic transitions from one small area (0 or 1) to another endemic area (0 or 1) without going through the intermediate step of the wide-area 01
  # Jumps = q23 and q32 (= jd0A, jd1A) = parameter labels designate source area NOT destination area
  # These are NOT "jump-dispersal" sensu BioGeoBEARS. They are anagenetic events of range-shift (a), not cladogenetic events (j)

## Range extensions
  # Anagenetic transitions from one of the endemic areas (0 and 1) to large area 01
  # Range extension = q21 and q31 (= d0A, d1A) = parameter labels designate source area

## Cladogenetic events
  # No "jump dispersal at speciation" = Cladogenetic transitions !
     # Transition from 2 to 2 and 3 = lambda223 set to zero
     # Transition from 3 to 2 and 3 = lambda323 set to zero
  # No "wide-spread" sympatry (as in DEC)
     # Transition from 1 to 1 and 1 = lambda111 need to be set to zero
  # Narrow vicariance is allowed
     # Transition from 1 to 2 and 3 = lambda123 = s01A in GeoSSE => parameter labels designate the source area
  # No distinction between endemic speciation (speciation without transitions) and sympatric subset speciation (speciation with transition/extirpation for one descendant)
    # lambda112 AND lambda222 = s0A in GeoSSE => parameter labels designate the destination area
      # Transition from 1 to 1 and 2 = lambda112 = sympatric subset speciation in endemic area 2
      # Transition from 2 to 2 and 2 = lambda222 = endemic speciation in endemic area 2
    # lambda113 AND lambda333 = s1A in GeoSSE => parameter labels designate the destination area
      # Transition from 1 to 1 and 3 = lambda113 = sympatric subset speciation in endemic area 3
      # Transition from 3 to 3 and 3 = lambda333 = endemic speciation in endemic area 3

# Use distinction between extirpation (source = 01) and local extinction (source = 0 or 1) as otherwise their is no anagenetic events allowing to leave the wide-spread area state (01)
# However, since it is impossible to distinguish between endemic speciation (source = 0 or 1) and sympatric subset speciation (source = 01), it remains impossible to distinguish net diversification rates in endemic areas (0 or 1) vs. wide-spread area (01)
# Use range-shift ("jumps") as cladogenetic jump-dispersal are not allowed, thus this is the only way to transition from area 0 to 1 without going through an intermediate wide-spread state (01)

## Can GeoHiSSE be used to run GeoSSE and their equivalent null models? CID-1, CID-2 and CID-4?
  # Set s0A = s1A => speciation rates independent from source/dest area 
  # Set x0A = x1A => extirpation rates independent from extirpation area 
  # Can apply this for hidden states B in case of CID-4 with effects of hidden states (A/B), but not observed ranges (0/1)


### 3.1/ Create design Q matrices for models ####

# Define equality and nullity of anagenetic transitions

## 3.1.1/ Create CID-1/GeoSSE design matrix (1 observed trait)

CID1_Qmat <- GeoSSE_Qmat <- hisse::TransMatMakerGeoHiSSE(hidden.traits = 0,
                                                         include.jumps = TRUE,
                                                         separate.extirpation = TRUE)
GeoSSE_Qmat

# # Version with no range-shift = jumps: jd0A (d00A_11A) and jd1A (d11A_00A) set to zero (as in a DEC model)
# CID1_no_jump_Qmat <- GeoSSE_no_jump_Qmat <- hisse::TransMatMakerGeoHiSSE(hidden.traits = 0,
#                                                                          include.jumps = FALSE,
#                                                                          separate.extirpation = TRUE)


# Only 1 trait with three states (00), (11), and (01)
# ARD model:
  # 1 forward rate and 1 backward rate for jumps/range-shifts => 3/ jd0A (d00A_11A)  and 1/ jd1A (d11A_00A) => 2 rates
  # 1 forward rate for range extension per source area => 5/ d0A (d00A_01A) and 6/ d1A (d11A_01A) => 2 rates
  # If separate.extirpation = F: 
     # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
     # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
     # than endemic speciation, which is an issue.
  # If separate.extirpation = T: 
     # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) and 4/ x0A (d01A_11A) => 2 rates
# Total = 6 transition parameters

## 3.1.2/ Create GeoHiSSE design matrix (1 binary observed trait x 1 binary hidden trait)

## Version with Uncorrelated design (the simplest null model)

GeoHiSSE_Qmat <- TransMatMakerGeoHiSSE(hidden.traits = 1,
                                       make.null = TRUE, # Transition between observed ranges are independent from hidden traits. Transitions between hidden traits are independent from observed areas.
                                       include.jumps = TRUE,
                                       separate.extirpation = TRUE)
GeoHiSSE_Qmat

# 2 traits with 3x2 crossed states: (00A), (11A), (01A), (00B), (11B) and (01B)
  # Observed states = (00), (11) and (01)
  # Hidden states = diversification regimes = (A) and (B)
  # No direct transitions across regimes x states: No (00A) <-> (11B) ; No (11A) <-> (00B)
# Uncorrelated design:
  # Hidden trait DOES NOT affect transition rates across observed states
    # q00A <-> q11A are equal to q00B <-> q11B
  # Observed trait DOES NOT affect transition rates across hidden states
    # q00A <-> q00B are equal to q11A <-> q11B
# ARD model for observed trait: 
  # 1 forward rate and 1 backward rate for jumps/range-shifts (independent from hidden states) => 3/ jd0A (d00A_11A) = jd0B (d00A_11A) and 1/ jd1A (d11A_00A) = jd1B (d11B_00B) => 2 rates
  # 1 forward rate for range extension per source area => 5/ d0A (d00A_01A) = d0B (d00B_01B) and 6/ d1A (d11A_01A) = d1B (d11B_01B) => 2 rates
  # If separate.extirpation = F: 
     # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
     # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
     # than endemic speciation, which is an issue.
  # If separate.extirpation = T: 
     # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) = x1B (d01B_00B) and 4/ x0A (d01A_11A) = x0B (d01B_11B) => 2 rates
# ER model for hidden trait: 1 unique rate (independent from observed state)
  # 7/ q0AB (d00A_00B) = q0BA (d00B_00A) = q1AB (d11A_11B) = q1BA (d11B_11A) = q01AB (d01A_01B) = q01BA (d01B_01A)
# Total = 7 transition parameters

# Idea = keep it the simplest: both trait transitions are independent from states in the other trait


## 3.1.3/ Create CID-2 design matrix (1 binary observed trait x 1 two-states hidden trait)

# Use CID-2 as null model against GeoSSE as their are truely two diversification regimes in GeoSSE
# (lambda123 is s01 is vicariance, but not diversification (no extinction under wide-spread state, only extirpation))

## Version with Uncorrelated design (the simplest null model)

CID2_Qmat <- TransMatMakerGeoHiSSE(hidden.traits = 1, 
                                   make.null = TRUE, # Transition between observed ranges are independent from hidden traits. Transitions between hidden traits are independent from observed areas.
                                   include.jumps = TRUE,
                                   separate.extirpation = TRUE)
CID2_Qmat

# 2 traits with 3x2 crossed states
  # Observed states = (00), (11) and (01)
  # Hidden states = diversification regimes = (A), (B)
  # No direct transitions across regimes x states: No (00A) <-> (11B) ; No (11A) <-> (01B) ; etc
# Uncorrelated design:
  # Hidden trait DOES NOT affect transition rates across observed states
    # q00A <-> q11A are equal to q00C <-> q11C
  # Observed trait DOES NOT affect transition rates across hidden states
    # q00A <-> q00C are equal to q11A <-> q11C
# ARD model for observed trait: 
  # 1 forward rate and 1 backward rate for jumps/range-shifts (independent from hidden states) => 3/ jd0A (d00A_11A) = jd0B (d00A_11A) and 1/ jd1A (d11A_00A) = jd1B (d11B_00B) => 2 rates
  # 1 forward rate for range extension per source area => 5/ d0A (d00A_01A) = d0B (d00B_01B) and 6/ d1A (d11A_01A) = d1B (d11B_01B) => 2 rates
  # If separate.extirpation = F: 
    # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
    # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
    # than endemic speciation, which is an issue.
  # If separate.extirpation = T: 
    # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) = x1B (d01B_00B) and 4/ x0A (d01A_11A) = x0B (d01B_11B) => 2 rates
# ER model for hidden trait: 1 unique rate (independent from observed state)
  # 7/ q0AB (d00A_00B) = q0BA (d00B_00A) = q1AB (d11A_11B) = q1BA (d11B_11A) = q01AB (d01A_01B) = q01BA (d01B_01A)
# Total = 7 transition parameters

# Idea = keep it the simplest: both trait transitions are independent from states in the other trait


# ## 3.1.3/ Create CID-3 design matrix (1 binary observed trait x 1 three-states hidden trait)
# 
# ## Should not be used as the number of true diversification regimes in GeoSSE is 2, not 3!
# 
# ## Version with Uncorrelated design (the simplest null model)
# 
# CID3_Qmat <- TransMatMakerGeoHiSSE(hidden.traits = 2, 
#                                    make.null = TRUE, # Transition between observed ranges are independent from hidden traits. Transitions between hidden traits are independent from observed areas.
#                                    include.jumps = TRUE,
#                                    separate.extirpation = TRUE)
# CID3_Qmat
# 
# # 2 traits with 3x3 crossed states
#   # Observed states = (00), (11) and (01)
#   # Hidden states = diversification regimes = (A), (B), (C)
#   # No direct transitions across regimes x states: No (00A) <-> (11B) ; No (11A) <-> (01B) ; etc
# # Uncorrelated design:
#   # Hidden trait DOES NOT affect transition rates across observed states
#     # q00A <-> q11A are equal to q00C <-> q11C
#   # Observed trait DOES NOT affect transition rates across hidden states
#     # q00A <-> q00C are equal to q11A <-> q11C
# # ARD model for observed trait: 
#   # 1 forward rate and 1 backward rate for jumps/range-shifts (independent from hidden states) => 3/ jd0A (d00A_11A) = jd0B (d00A_11A) = jd0C (d00C_11C) and 1/ jd1A (d11A_00A) = jd1B (d11B_00B) = jd1C (d11C_00C) => 2 rates
#   # 1 forward rate for range extension per source area => 5/ d0A (d00A_01A) = d0B (d00B_01B) = d0C (d00C_01C) and 6/ d1A (d11A_01A) = d1B (d11B_01B) = d1C (d11C_01C) => 2 rates
#   # If separate.extirpation = F: 
#     # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
#     # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
#     # than endemic speciation, which is an issue.
#   # If separate.extirpation = T: 
#     # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) = x1B (d01B_00B) = x1C (d01C_00C) and 4/ x0A (d01A_11A) = x0B (d01B_11B) = x0C (d01C_11C) => 2 rates
# # ER model for hidden trait: 1 unique rate (independent from observed state)
#   # 7/ q0AB (d00A_00B) = q0BA (d00B_00A) = q1AB (d11A_11B) = q1BA (d11B_11A) = q01AB (d01A_01B) = q01BA (d01B_01A) = q0AC (d00A_00C) = q0CA (d00C_00A) = ...
# # Total = 7 transition parameters
# 
# # Idea = keep it the simplest: both trait transitions are independent from states in the other trait


## 3.1.4/ Create CID-4 design matrix (1 binary observed trait x 1 four-states hidden trait)

# Use CID-4 as null model against GeoHiSSE as their are truely four diversification regimes in GeoHiSSE
# (lambda123 is s01 is vicariance, but not diversification (no extinction under wide-spread state, only extirpation))

## Version with Uncorrelated design (the simplest null model)

CID4_Qmat <- TransMatMakerGeoHiSSE(hidden.traits = 3, 
                                   make.null = TRUE, # Transition between observed ranges are independent from hidden traits. Transitions between hidden traits are independent from observed areas.
                                   include.jumps = TRUE,
                                   separate.extirpation = TRUE)
CID4_Qmat

# 2 traits with 3x4 crossed states
  # Observed states = (00), (11) and (01)
  # Hidden states = diversification regimes = (A), (B), (C), (D)
  # No direct transitions across regimes x states: No (00A) <-> (11B) ; No (11A) <-> (01B) ; etc
# Uncorrelated design:
  # Hidden trait DOES NOT affect transition rates across observed states
    # q00A <-> q11A are equal to q00C <-> q11C
  # Observed trait DOES NOT affect transition rates across hidden states
    # q00A <-> q00C are equal to q11A <-> q11C
# ARD model for observed trait: 
  # 1 forward rate and 1 backward rate for jumps/range-shifts (independent from hidden states) => 3/ jd0A (d00A_11A) = jd0B (d00A_11A) = ... = jd0D (d00D_11D) = and 1/ jd1A (d11A_00A) = jd1B (d11B_00B) = ... = jd1D (d11D_00D) => 2 rates
  # 1 forward rate for range extension per source area =>  5/ d0A (d00A_01A) = d0B (d00B_01B) = ... = d0D (d00D_01D) and 6/ d1A (d11A_01A) = d1B (d11B_01B) = ... = d1F (d11D_01D) => 2 rates
  # If separate.extirpation = F: 
    # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
    # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
    # than endemic speciation, which is an issue.
  # If separate.extirpation = T: 
    # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) = x1B (d01B_00B) = ... = x1D (d01D_00D) and 4/ x0A (d01A_11A) = x0B (d01B_11B) = ... = x0F (d01D_11D) => 2 rates
# ER model for hidden trait: 1 unique rate (independent from observed state)
  # 7/ q0AB (d00A_00B) = q0BA (d00B_00A) = q1AB (d11A_11B) = q1BA (d11B_11A) = q01AB (d01A_01B) = q01BA (d01B_01A) = ... = q0AD (d00A_00D) = q0DA (d00D_00A) = ...
# Total = 7 transition parameters

# Idea = keep it the simplest: both trait transitions are independent from states in the other trait


# ## 3.1.4/ Create CID-6 design matrix (1 binary observed trait x 1 six-states hidden trait)
# 
# ## Should not be used as the number of true diversification regimes in GeoHiSSE is 4, not 6!
# 
# ## Version with Uncorrelated design (the simplest null model)
# 
# CID6_Qmat <- TransMatMakerGeoHiSSE(hidden.traits = 5, 
#                                    make.null = TRUE, # Transition between observed ranges are independent from hidden traits. Transitions between hidden traits are independent from observed areas.
#                                    include.jumps = TRUE,
#                                    separate.extirpation = TRUE)
# CID6_Qmat
# 
# # 2 traits with 3x6 crossed states
#   # Observed states = (00), (11) and (01)
#   # Hidden states = diversification regimes = (A), (B), (C), (D), (E), (F)
#   # No direct transitions across regimes x states: No (00A) <-> (11B) ; No (11A) <-> (01B) ; etc
# # Uncorrelated design:
#   # Hidden trait DOES NOT affect transition rates across observed states
#     # q00A <-> q11A are equal to q00C <-> q11C
#   # Observed trait DOES NOT affect transition rates across hidden states
#     # q00A <-> q00C are equal to q11A <-> q11C
# # ARD model for observed trait: 
#   # 1 forward rate and 1 backward rate for jumps/range-shifts (independent from hidden states) => 3/ jd0A (d00A_11A) = jd0B (d00A_11A) = ... = jd0F (d00F_11F) = and 1/ jd1A (d11A_00A) = jd1B (d11B_00B) = ... = jd1F (d11F_00F) => 2 rates
#   # 1 forward rate for range extension per source area =>  5/ d0A (d00A_01A) = d0B (d00B_01B) = ... = d0F (d00F_01F) and 6/ d1A (d11A_01A) = d1B (d11B_01B) = ... = d1F (d11F_01F) => 2 rates
#   # If separate.extirpation = F: 
#     # No backward rates as range extirpation along branches is forbidden (only extinction is allowed).
#     # But then, the only way to escape the wide-range state is through cladogenetic subset speciation which is controlled by the same parameter
#     # than endemic speciation, which is an issue.
#   # If separate.extirpation = T: 
#     # 1 backward rate for range contraction/extirpation per extirpation area => 2/ x1A (d01A_00A) = x1B (d01B_00B) = ... = x1F (d01F_00F) and 4/ x0A (d01A_11A) = x0B (d01B_11B) = ... = x0F (d01F_11F) => 2 rates
# # ER model for hidden trait: 1 unique rate (independent from observed state)
#   # 7/ q0AB (d00A_00B) = q0BA (d00B_00A) = q1AB (d11A_11B) = q1BA (d11B_11A) = q01AB (d01A_01B) = q01BA (d01B_01A) = ... = q0AF (d00A_00F) = q0FA (d00F_00A) = ...
# # Total = 7 transition parameters
# 
# # Idea = keep it the simplest: both trait transitions are independent from states in the other trait


### 3.2/ Set diversification rate design for models ####

## Create CID-1 diversification rate design: 

# # Version with non-independent vicariance
# CID1_turnover_rates <- c(1,1,1) # s0A = s1A = s01A
# CID1_eps_rates <- c(1,1) # x0A = x1A # No extirpation from 01 as this would be non-allowed transitions

# Version with independent vicariance
CID1_turnover_rates <- c(1,1,2) # s0A = s1A. Vicariance = s01A is independent
CID1_eps_rates <- c(1,1) # x0A = x1A # No extirpation from 01 as this would be non-allowed transitions

# We want the observed states to DO NOT affect turnover and extinction fraction
# Thus, only a single turnover rate and a single extinction fraction
# However, we do not want the vicariance to depend on turnover rates in unique areas as it does not have an extinction rates, but extirpation rates

## Create GeoSSE diversification rate design:
GeoSSE_turnover_rates <- c(1,2,3) #  s0A, s1A, s01A
GeoSSE_eps_rates <- c(1,2) # x0A and x1A # No extirpation from 01 as this would be non-allowed transitions

# We want the observed states to affect both turnover and extinction fraction
# If the extinction fraction is fixed to 0, we model speciation alone

## Create CID-2 diversification rate design:

# # Version with non-independent vicariances
# CID2_turnover_rates <- c(1,1,1,2,2,2) # s0A = s1A = s01A and s0B = s1B = s01B
# CID2_eps_rates <- c(1,1,2,2) # x0A = x1A and x0B = x1B # No extirpation from 01 as this would be non-allowed transitions

# Version with independent vicariances
CID2_turnover_rates <- c(1,1,3,2,2,4) # s0A = s1A and s0B = s1B. Vicariances are independent from turnover in unique areas and different across regimes: s01A =/= s01B
CID2_eps_rates <- c(1,1,2,2) # x0A = x1A and x0B = x1B # No extirpation from 01 as this would be non-allowed transitions

# We want the hidden states to affect both turnover and extinction fraction
# If the extinction fraction is fixed to 0, we model speciation alone
# However, we do not want the vicariance to depend on turnover rates in unique areas as it does not have an extinction rates, but extirpation rates

# ## Create CID-3 diversification rate design:
# CID3_turnover_rates <- c(1,1,1,2,2,2,3,3,3) # s0A = s1A = s01A, s0B = s1B = s01B and s0C = s1C = s01C
# CID3_eps_rates <- c(1,1,2,2,3,3) # x0A = x1A, x0B = x1B, and x0C = x1C # No extirpation from 01 as this would be non-allowed transitions
# 
# # We want the hidden states to affect both turnover and extinction fraction
# # If the extinction fraction is fixed to 0, we model speciation alone

## Create GeoHiSSE diversification rate design:
GeoHiSSE_turnover_rates <- c(1,2,3,4,5,6) # s0A, s1A, s01A and s0B, s1B, s01B
GeoHiSSE_eps_rates <- c(1,2,3,4) # x0A, x1A and x0B, x1B # No extirpation from 01 as this would be non-allowed transitions 

# We want the observed states to affect both turnover and extinction fraction
# We want the hidden states to affect both turnover and extinction fraction
# If the extinction fraction is fixed to 0, we model speciation alone

## Create CID-4 diversification rate design:

# # Version with non-independent vicariances
# CID4_turnover_rates <- c(1,1,1,2,2,2,3,3,3,4,4,4) # s0A = s1A = s01A, s0B = s1B = s01B, ..., s0D = s1D = s01D
# CID4_eps_rates <- c(1,1,2,2,3,3,4,4) # x0A = x1A, x0B = x1B, ..., x0D = x1D # No extirpation from 01 as this would be non-allowed transitions

# Version with independent vicariances
CID4_turnover_rates <- c(1,1,5,2,2,6,3,3,7,4,4,8) # s0A = s1A, s0B = s1B, ..., s0D = s1D. Vicariances are independent from turnover in unique areas and different across regimes: s01A =/= s01B ... =/= s01D
CID4_eps_rates <- c(1,1,2,2,3,3,4,4) # x0A = x1A, x0B = x1B, ..., x0D = x1D # No extirpation from 01 as this would be non-allowed transitions


# We want the observed states to DO NOT affect turnover and extinction fraction
# We want the hidden states to affect both turnover and extinction fraction
# If the extinction fraction is fixed to 0, we the model speciation alone

# ## Create CID-6 diversification rate design:
# CID6_turnover_rates <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6) # s0A = s1A = s01A, s0B = s1B = s01B, ..., s0F = s1F = s01F
# CID6_eps_rates <- c(1,1,2,2,3,3,4,4,5,5,6,6) # x0A = x1A, x0B = x1B, ..., x0F = x1F # No extirpation from 01 as this would be non-allowed transitions
# 
# # We want the observed states to DO NOT affect turnover and extinction fraction
# # We want the hidden states to affect both turnover and extinction fraction
# # If the extinction fraction is fixed to 0, we the model speciation alone

### 3.3/ Fit models with ML optimization in hisse ####

## 3.3.1/ Fit CID-1 model using GeoHiSSE ####

CID1_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
                                data = OW_NW_GeoSSE_tip_data_df,
                                turnover = CID1_turnover_rates, # One single rate
                                eps = CID1_eps_rates, # One single rate
                                hidden.states = FALSE, # No hidden states
                                # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
                                trans.rate = CID1_Qmat, # 6 transitions rates between observed states:
                                  # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
                                  # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
                                  # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
                                sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
CID1_MLE

# # Test version with no range-shift (jumps)
# CID1_no_jump_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
#                             data = OW_NW_GeoSSE_tip_data_df,
#                             turnover = CID1_turnover_rates, # One single rate
#                             eps = CID1_eps_rates, # One single rate
#                             hidden.states = FALSE, # No hidden states
#                             # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
#                             trans.rate = CID1_no_jump_Qmat, # 4 transitions rates between observed states:
#                                # 2 extirpations: 1/ x1A (d01A_00A), 2/ x0A (d01A_11A)
#                                # 2 range extensions: 3/ d0A (d00A_01A), 4/ d1A (d11A_01A)
#                             sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
# CID1_no_jump_MLE
# CID1_no_jump_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID1_no_jump_MLE)
# CID1_no_jump_MLE_div_rates
# 
# CID1_no_jump_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID1_no_jump_MLE)
# CID1_no_jump_MLE_Q_matrix
#
# # Does not affect diversification rates, but lead to even more extreme transition rates. Favor the version with range-shift allowed.

# Save model outputs
# saveRDS(object = CID1_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID1_MLE.rds")
saveRDS(object = CID1_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID1_MLE.rds")

# Compute speciation and extinction rates per states (in event / My / lineage)
CID1_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID1_MLE)
CID1_MLE_div_rates

# Compare to rates estimated with BD model 
# Likelihood are different because the BD model does not account for trait data evolution
# But rate estimates should be similar/compatible
phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)

# Compute mean expected time (in My) before speciation/extinction event on a branch
1/CID1_MLE_div_rates
# Extremely short expected time in Wide-spread state to compensate for the absence of jump-dispersal

# Extract transition rates = Q matrix
CID1_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID1_MLE)
CID1_MLE_Q_matrix


## 3.3.2/ Fit GeoSSE model using GeoHiSSE ####

GeoSSE_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
                                  data = OW_NW_GeoSSE_tip_data_df,
                                  turnover = GeoSSE_turnover_rates, # One turnover rate per observed state = 3 rates
                                  eps = GeoSSE_eps_rates, # One extinct fraction per endemic area = 2 rates
                                  hidden.states = FALSE, # No hidden states
                                  # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
                                  trans.rate = GeoSSE_Qmat, # 6 transitions rates between states:
                                    # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
                                    # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
                                    # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
                                  sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
GeoSSE_MLE

# Save model outputs
# saveRDS(object = GeoSSE_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_MLE.rds")
saveRDS(object = GeoSSE_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_MLE.rds")

# Compute speciation and extinction rates per states (in event / My / lineage)
GeoSSE_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(GeoSSE_MLE)
GeoSSE_MLE_div_rates

# Compare to rates estimated with BD model 
# Likelihood are different because the BD model does not account for trait data evolution
# But rate estimates should be similar/compatible
phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)

# Compute mean expected time (in My) before speciation/extinction event on a branch
1/GeoSSE_MLE_div_rates

# Extract transition rates = Q matrix
GeoSSE_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(GeoSSE_MLE)
GeoSSE_MLE_Q_matrix


## 3.3.3/ Fit CID-2 model using GeoHiSSE ####

CID2_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
                                data = OW_NW_GeoSSE_tip_data_df,
                                turnover = CID2_turnover_rates, # 2 equal turnover rates across 2 hidden states + 2 vicariance rates = 4 rates
                                eps = CID2_eps_rates, # 2 equal extinction fraction across 2 hidden states
                                hidden.states = TRUE, # 2 hidden states
                                # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
                                trans.rate = CID2_Qmat, # 7 transitions rates between states:
                                # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
                                # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
                                # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
                                # 1 regime transition: 7/ regime transition rates between A and B
                                sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
CID2_MLE

# Save model outputs
# saveRDS(object = CID2_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID2_MLE.rds")
saveRDS(object = CID2_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID2_MLE.rds")

# Compute speciation and extinction rates per states (in event / My / lineage)
CID2_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID2_MLE)
CID2_MLE_div_rates

# Compare to rates estimated with BD model 
# Likelihood are different because the BD model does not account for trait data evolution
# But rate estimates should be similar/compatible
phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)

# Compute mean expected time (in My) before speciation/extinction event on a branch
1/CID2_MLE_div_rates

# Extract transition rates = Q matrix
CID2_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID2_MLE)
CID2_MLE_Q_matrix


# ## 3.3.3/ Fit CID-3 model using GeoHiSSE ####
# 
# CID3_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
#                                 data = OW_NW_GeoSSE_tip_data_df,
#                                 turnover = CID3_turnover_rates, # 3 turnover rates across 3 hidden states
#                                 eps = CID3_eps_rates, # 3 extinction fraction across 3 hidden states
#                                 hidden.states = TRUE, # 3 hidden states
#                                 # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
#                                 trans.rate = CID3_Qmat, # 7 transitions rates between states:
#                                   # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
#                                   # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
#                                   # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
#                                   # 1 regime transition: 7/ regime transition rates between A and B
#                                 sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
# CID3_MLE
# 
# # Save model outputs
# # saveRDS(object = CID3_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID3_MLE.rds")
# saveRDS(object = CID3_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID3_MLE.rds")
# 
# # Compute speciation and extinction rates per states (in event / My / lineage)
# CID3_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID3_MLE)
# CID3_MLE_div_rates
# 
# # Compare to rates estimated with BD model 
# # Likelihood are different because the BD model does not account for trait data evolution
# # But rate estimates should be similar/compatible
# phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)
# 
# # Compute mean expected time (in My) before speciation/extinction event on a branch
# 1/CID3_MLE_div_rates
# 
# # Extract transition rates = Q matrix
# CID3_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID3_MLE)
# CID3_MLE_Q_matrix


## 3.3.4/ Fit GeoHiSSE model using GeoHiSSE ####

GeoHiSSE_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
                                    data = OW_NW_GeoSSE_tip_data_df,
                                    turnover = GeoHiSSE_turnover_rates, # 6 turnover rates: 3 observed states x 2 hidden states
                                    eps = GeoHiSSE_eps_rates, # 4 extinct fractions: 2 endemic areas x 2 hidden states
                                    hidden.states = TRUE, # 2 hidden states
                                    # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
                                    trans.rate = GeoHiSSE_Qmat, # 7 transitions rates between states:
                                      # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
                                      # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
                                      # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
                                      # 1 regime transition: 7/ regime transition rates between A and B
                                    sann = TRUE, # Simulated annealing step can be shut down to save time (but lower chance of convergence)
                                    print_level_subplex = 1,
                                    verbose_genSA = TRUE) 
GeoHiSSE_MLE

# Save model outputs
# saveRDS(object = GeoHiSSE_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoHiSSE_MLE.rds")
saveRDS(object = GeoHiSSE_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoHiSSE_MLE.rds")

# Compute speciation and extinction rates per states (in event / My / lineage)
GeoHiSSE_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(GeoHiSSE_MLE)
GeoHiSSE_MLE_div_rates

# Compare to rates estimated with BD model 
# Likelihood are different because the BD model does not account for trait data evolution
# But rate estimates should be similar/compatible
phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)

# Compute mean expected time (in My) before speciation/extinction event on a branch
1/GeoHiSSE_MLE_div_rates

# Extract transition rates = Q matrix
GeoHiSSE_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(GeoHiSSE_MLE)
GeoHiSSE_MLE_Q_matrix

## 3.3.5/ Fit CID-4 model using GeoHiSSE ####

CID4_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
                                data = OW_NW_GeoSSE_tip_data_df,
                                turnover = CID4_turnover_rates, # 2 equal turnover rates across 4 hidden states + 4 vicariance rates = 8 rates
                                eps = CID4_eps_rates, # 2 equal extinct fractions across 4 hidden states
                                hidden.states = TRUE, # 4 hidden states
                                # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
                                trans.rate = CID4_Qmat, # 7 transitions rates between states:
                                # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
                                # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
                                # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
                                # 1 regime transition: 7/ regime transition rates between A and B
                                sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
CID4_MLE

# Save model outputs
# saveRDS(object = CID4_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID4_MLE.rds")
saveRDS(object = CID4_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID4_MLE.rds")

# Compute speciation and extinction rates per states (in event / My / lineage)
CID4_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID4_MLE)
CID4_MLE_div_rates

# Compare to rates estimated with BD model 
# Likelihood are different because the BD model does not account for trait data evolution
# But rate estimates should be similar/compatible
phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)

# Compute mean expected time (in My) before speciation/extinction event on a branch
1/CID4_MLE_div_rates

# Extract transition rates = Q matrix
CID4_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID4_MLE)
CID4_MLE_Q_matrix

# ## 3.3.5/ Fit CID-6 model using GeoHiSSE ####
# 
# CID6_MLE <- GeoHiSSE_with_print(phy = Ponerinae_phylogeny_1534t_treedata@phylo, 
#                                 data = OW_NW_GeoSSE_tip_data_df,
#                                 turnover = CID6_turnover_rates, # 6 turnover rates across hidden states
#                                 eps = CID6_eps_rates, # 6 extinct fractions across hidden states
#                                 hidden.states = TRUE, # 6 hidden states
#                                 # starting.vals = c(0.1, 0.2, 0.01), # Taus, Efs, Dispersal = Turnovers, Extinct fractions, Transitions
#                                 trans.rate = CID6_Qmat, # 7 transitions rates between states:
#                                   # 2 range-shifts: 1/ jd1A (d11A_00A), 3/ jd0A (d00A_11A)
#                                   # 2 extirpations: 2/ x1A (d01A_00A), 4/ x0A (d01A_11A)
#                                   # 2 range extensions: 5/ d0A (d00A_01A), 6/ d1A (d11A_01A)
#                                   # 1 regime transition: 7/ regime transition rates between A and B
#                                 sann = TRUE) # Simulated annealing step can be shut down to save time (but lower chance of convergence)
# CID6_MLE
# 
# # Save model outputs
# # saveRDS(object = CID6_MLE, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID6_MLE.rds")
# saveRDS(object = CID6_MLE, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID6_MLE.rds")
# 
# # Compute speciation and extinction rates per states (in event / My / lineage)
# CID6_MLE_div_rates <- convert_rates_from_GeoHiSSE_output(CID6_MLE)
# CID6_MLE_div_rates
# 
# # Compare to rates estimated with BD model 
# # Likelihood are different because the BD model does not account for trait data evolution
# # But rate estimates should be similar/compatible
# phytools::fit.bd(Ponerinae_phylogeny_1534t_treedata@phylo)
# 
# # Compute mean expected time (in My) before speciation/extinction event on a branch
# 1/CID6_MLE_div_rates
# 
# # Extract transition rates = Q matrix
# CID6_MLE_Q_matrix <- extract_Q_matrix_from_GeoHiSSE_output(CID6_MLE)
# CID6_MLE_Q_matrix


##### 4/ Compare ML models with AICc and Akaike's weights ####

# Load model outputs
# CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID1_MLE.rds")
# GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID3_MLE.rds")
# GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoHiSSE_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID6_MLE.rds")

CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID1_MLE.rds")
GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_MLE.rds")
CID2_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID2_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID3_MLE.rds")
GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoHiSSE_MLE.rds")
CID4_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID4_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID6_MLE.rds")

### 4.1/ Home-made functions to extract model evaluations ####

?MuMIn::AICc() # Does not work because need a method to extract nobs

# Home-made function to compute AIC from logLk and number of parameters (k)
compute_AIC <- function (logLk, k)
{
  AIC <- - 2 * (logLk - k)
}

# Home-made function to compute AICc from AIC, number of observations, and number of parameters (k)
compute_AICc <- function (AIC, nobs, k)
{
  AICc <- AIC + (2 * k) * (k - 1) / (nobs - k - 1)
}

### 4.2/ Compare models ####

# Extract number of observations for all models
nobs <- length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)

# Combine models' fit n a list
# GeoSSE_models_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE)
GeoSSE_models_list <- list(CID1_MLE, GeoSSE_MLE, CID2_MLE, GeoHiSSE_MLE, CID4_MLE)
# GeoSSE_models_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE, CID6_MLE)

attributes(GeoSSE_models_list[[1]])
max(GeoSSE_models_list[[1]]$index.par) - 1

# Extract ln-likelihood, number of parameters, and AIC from models' list
GeoSSE_models_comparison_df <- data.frame(# model = c("CID-1", "GeoSSE", "CID-2", "GeoHiSSE"),
                                          model = c("CID-1", "GeoSSE", "CID-2", "GeoHiSSE", "CID-4"),
                                          # model = c("CID-1", "GeoSSE", "CID-3", "GeoHiSSE", "CID-6"),
                                          logLk = sapply(X = GeoSSE_models_list, FUN = function(x) x$loglik),
                                          k = sapply(X = GeoSSE_models_list, FUN = function(x) max(x$index.par) - 1))

# Inform type of models, nb hidden states and diversification regimes
GeoSSE_models_comparison_df$Model_type <- c("Dependent", "Independent")[(str_detect(string = GeoSSE_models_comparison_df$model, pattern = "CID") + 1)]
GeoSSE_models_comparison_df$Hidden_states <- unlist(lapply(X = GeoSSE_models_list, FUN = function (x) { length(unique(str_extract(string = colnames(x$trans.matrix), pattern = "[A-Z]")))} ))
GeoSSE_models_comparison_df$Div_regimes <- GeoSSE_models_comparison_df$Hidden_states * c(1,2)[(GeoSSE_models_comparison_df$Model_type == "Dependent") + 1]
GeoSSE_models_comparison_df$Hidden_states[GeoSSE_models_comparison_df$Hidden_states == 1] <- 0

# Compute AIC from logLk and nobs
GeoSSE_models_comparison_df$AIC <- unlist(purrr::map2(.x = GeoSSE_models_comparison_df$logLk, .y = GeoSSE_models_comparison_df$k, .f = compute_AIC))

# Compute AICc from AIC, k and nobs
GeoSSE_models_comparison_df$AICc <- unlist(purrr::map2(.x = GeoSSE_models_comparison_df$AIC, nobs = nobs, .y = GeoSSE_models_comparison_df$k, .f = compute_AICc))

# Compute delta AICc
best_AICc <- min(GeoSSE_models_comparison_df$AICc)
GeoSSE_models_comparison_df$delta_AICc <- GeoSSE_models_comparison_df$AICc - best_AICc

# Compute Akaike's weights from AIC
GeoSSE_models_comparison_df$Akaike_weights <- round(phytools::aic.w(GeoSSE_models_comparison_df$AICc), 3) * 100

# Add ranks
GeoSSE_models_comparison_df$rank <- rank(GeoSSE_models_comparison_df$AICc)

# Reorder columns
GeoSSE_models_comparison_df <- GeoSSE_models_comparison_df %>% 
  dplyr::select(model, Model_type, Hidden_states, Div_regimes, k, logLk, AICc, delta_AICc, Akaike_weights, rank)

# Display result
GeoSSE_models_comparison_df

# Save GeoSSE model comparison summary table
# saveRDS(object = GeoSSE_models_comparison_df, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.rds")
saveRDS(object = GeoSSE_models_comparison_df, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.rds")

# Export model comparison df
# xlsx::write.xlsx(x = GeoSSE_models_comparison_df, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.xlsx")
xlsx::write.xlsx(x = GeoSSE_models_comparison_df, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.xlsx")


##### 5/ Compute residence time in each state ####

# Quick version:
  # Obtain marginal likelihoods of nodes/tips
  # Compute edge proba as mean of rootward/tipward nodes/tips
  # Get residence times as the weighted sum of branch length
  # Aggregate across regimes and observed states
  # Limit: does not account for transition rates along branches

# Ideal version: (time-consuming? overkill?)
  # Obtain marginal likelihoods of nodes/tips
  # Run BSM
  # Generate simmap
  # Get residence times from simmap
  # Aggregate residence times across BSM
  # Pros: Does account for transition rates along branches

?hisse::MarginReconGeoSSE() # Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm

# Load GeoSSE model comparison summary table
# GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.rds")
GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.rds")

# Load model outputs
# CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID1_MLE.rds")
# GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID3_MLE.rds")
# GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoHiSSE_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID6_MLE.rds")

CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID1_MLE.rds")
GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_MLE.rds")
CID2_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID2_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID3_MLE.rds")
GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoHiSSE_MLE.rds")
CID4_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID4_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID6_MLE.rds")

# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE)
# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE)
model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID2_MLE, GeoHiSSE_MLE, CID4_MLE)
# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE, CID6_MLE)

# model_list <- c("CID-1", "GeoSSE", "CID-3")
# model_list <- c("CID-1", "GeoSSE", "CID-3", "GeoHiSSE")
model_list <- c("CID-1", "GeoSSE", "CID-2", "GeoHiSSE", "CID-4")
# model_list <- c("CID-1", "GeoSSE", "CID-3", "GeoHiSSE", "CID-6")

## Loop per model
residence_times_all_models_df <- data.frame()
for (i in seq_along(model_list))
{
  # i <- 1
  
  model_i <- model_list[i]
  model_MLE_i <- model_fit_list[[i]]
  
  ### 5.1/ Get marginal likelihood of ARE ####

  # Extract number of shifts (= number of hidden states - 1)
  pars <- model_MLE_i$solution[model_MLE_i$solution != 0]
  par_regimes <- str_sub(names(pars), start = nchar(names(pars)))
  nb_regimes <- length(unique(par_regimes))

  # SSE_model_marginal_ARE_i <- hisse::MarginReconGeoSSE(phy = model_MLE_i$phy,
  #                                                      data = model_MLE_i$data,
  #                                                      pars = model_MLE_i$solution,
  #                                                      hidden.states = nb_regimes,
  #                                                      # 1 for CID-1 and GeoSSE ; 2 for GeoHiSSE ; 2 for CID-2 ; 3 for CID-3 ; 4 for CID-4 ; 6 for CID-6
  #                                                      f = model_MLE_i$f,
  #                                                      AIC = model_MLE_i$AIC,
  #                                                      root.type = model_MLE_i$root.type) # "madfitz" by default
  # 
  # dim(model_MLE_i$data)
  # 
  # View(SSE_model_marginal_ARE_i$node.mat)
  # dim(SSE_model_marginal_ARE_i$tip.mat)
  # dim(SSE_model_marginal_ARE_i$rates.mat)
  # 
  # # Save marginal probabilities
  # # saveRDS(object = SSE_model_marginal_ARE_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Marginal_ARE_",model_i,".rds"))
  # saveRDS(object = SSE_model_marginal_ARE_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Marginal_ARE_",model_i,".rds"))
  
  ## 5.2/ Approximate edge probabilities as mean of rootward/tipward nodes ####

  # Load marginal probabilities of states
  # SSE_model_marginal_ARE_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Marginal_ARE_",model_i,".rds"))
  SSE_model_marginal_ARE_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Marginal_ARE_",model_i,".rds"))

  # Extract marginal probability for all nodes (internal + tips)
  node_states_df_i <- rbind(SSE_model_marginal_ARE_i$tip.mat, SSE_model_marginal_ARE_i$node.mat)
  node_states_df_i <- node_states_df_i[, -1]

  # Restrict states to the used ones
  used_regimes <- LETTERS[1:nb_regimes]
  used_states_indices <- str_detect(string = colnames(node_states_df_i), pattern = paste(used_regimes, collapse = "|"))
  node_states_df_i <- node_states_df_i[, used_states_indices]

  # From node probabilities to edge probabilities as the mean of probabilities of the two tips (rootward and backward) of each edge
  edge_states_df_i <- node_states_df_i[-1, ]
  for (i in 1:nrow(edge_states_df_i))
  {
    rootward_node_prob <- node_states_df_i[SSE_model_marginal_ARE_i$phy$edge[i, 1], ]
    tipward_node_prob <- node_states_df_i[SSE_model_marginal_ARE_i$phy$edge[i, 2], ]
    edge_prob <- apply(X = rbind(rootward_node_prob, tipward_node_prob), MARGIN = 2, FUN = mean)
    edge_states_df_i[i, ] <- edge_prob
  }

  ## 5.3/ Compute residence times per state as weighted sum of branch length ####

  edge_lengths <- SSE_model_marginal_ARE_i$phy$edge.length
  branch_lengths_per_state <- apply(X = edge_states_df_i, MARGIN = 2, FUN = function(x) { y <- as.numeric(t(edge_lengths) %*% x) })
  branch_lengths_perc_per_state <- branch_lengths_per_state / sum(branch_lengths_per_state) * 100

  residence_times_df_model_i <- data.frame(state = names(branch_lengths_per_state),
                                           residence_time = branch_lengths_per_state,
                                           residence_time_perc = branch_lengths_perc_per_state)

  ## 5.4/ Aggregate residence times across regimes and observed states ####

  # Remove parenthesis
  residence_times_df_model_i$state <- str_remove(string = residence_times_df_model_i$state, pattern = "\\(")
  residence_times_df_model_i$state <- str_remove(string = residence_times_df_model_i$state, pattern = "\\)")
  # Extract range (observed states)
  residence_times_df_model_i$range <- str_sub(string = residence_times_df_model_i$state, end = nchar(residence_times_df_model_i$state) - 1)
  # Extract regime (hidden states)
  residence_times_df_model_i$regime <- str_sub(string = residence_times_df_model_i$state, start = nchar(residence_times_df_model_i$state))

  # Rename ranges
  residence_times_df_model_i$range[residence_times_df_model_i$range == "00"] <- "NW"
  residence_times_df_model_i$range[residence_times_df_model_i$range == "11"] <- "OW"
  residence_times_df_model_i$range[residence_times_df_model_i$range == "01"] <- "WS"
  # Update states
  residence_times_df_model_i$state <- paste0(residence_times_df_model_i$range, "-", residence_times_df_model_i$regime)


  # Compute percentages per range
  residence_times_df_model_i <- residence_times_df_model_i %>%
    group_by(range) %>%
    mutate(residence_time_perc_per_range = residence_time / sum(residence_time) * 100)

  # Aggregate per ranges
  residence_times_ranges_df_model_i <- residence_times_df_model_i %>%
    group_by(range) %>%
    summarize(residence_time = sum(residence_time),
              residence_time_perc = sum(residence_time_perc),
              residence_time_perc_per_range = sum(residence_time_perc_per_range)) %>%
    mutate(state = range,
           regime = NA)

  # Aggregate per regimes
  residence_times_regimes_df_model_i <- residence_times_df_model_i %>%
    group_by(regime) %>%
    summarize(residence_time = sum(residence_time),
              residence_time_perc = sum(residence_time_perc),
              residence_time_perc_per_range = sum(residence_time_perc_per_range)) %>%
    mutate(state = regime,
           range = NA)

  # Merge all residence times
  residence_times_df_model_i <- rbind(residence_times_df_model_i, residence_times_ranges_df_model_i, residence_times_regimes_df_model_i)
  residence_times_df_model_i$model <- model_i

  ## Merge residence times for all models
  residence_times_all_models_df <- rbind(residence_times_all_models_df, residence_times_df_model_i)

  ## Save residence times summary table
  # saveRDS(residence_times_all_models_df, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/residence_times_all_models_df.rds")
  saveRDS(residence_times_all_models_df, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/residence_times_all_models_df.rds")
  
  ## Print progress
  cat(paste0(Sys.time(), " - Residence times computed for ", model_i," - n° ", i, "/", length(model_list),"\n"))
  
}

## Save residence times summary table
# saveRDS(residence_times_all_models_df, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/residence_times_all_models_df.rds")
saveRDS(residence_times_all_models_df, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/residence_times_all_models_df.rds")


### 5.5/ Plot residence times ####

# Load residence times summary table
# residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/residence_times_all_models_df.rds")
residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/residence_times_all_models_df.rds")

# residence_times_all_models_df <- residence_times_df_model_i

# Remove aggregated times
residence_times_all_models_df_for_plot <- residence_times_all_models_df %>% 
  filter(!is.na(range)) %>%
  filter(!is.na(regime))

## Set color scheme for states

# Color = Ranges x Regimes
  # Purple shades for Old World
  # Beige shades for New World
  # Blue shades for Both
  # Lighter = A
  # Darker = F

colors_for_OW_NW <- c("peachpuff2", "mediumpurple2")

# Set colors for New World states
col_NW_fn <- colorRampPalette(colors = c("plum", "mediumpurple2", "purple3"))
col_NW_scale <- col_NW_fn(n = 6)
names(col_NW_scale) <- paste0("NW-", LETTERS[1:6])

# Set colors for Old World states
col_OW_fn <- colorRampPalette(colors = c("papayawhip", "peachpuff2", "sandybrown"))
col_OW_scale <- col_OW_fn(n = 6)
names(col_OW_scale) <- paste0("OW-", LETTERS[1:6])

# Set colors for Wide-spread states
col_WS_fn <- colorRampPalette(colors = c("lightskyblue", "royalblue", "royalblue4"))
col_WS_scale <- col_WS_fn(n = 6)
names(col_WS_scale) <- paste0("WS-", LETTERS[1:6])

color_scheme_states <- c(col_NW_scale, col_OW_scale, col_WS_scale)

# Set color scheme for regimes
regimes_col <- rev(RColorBrewer::brewer.pal(name = "Spectral", n = 6))

# Adjust model order
residence_times_all_models_df_for_plot$model <- factor(x = residence_times_all_models_df_for_plot$model,
                                                       levels = model_list, labels = model_list)

# pdf(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Residence_times_all_models_ggplot.pdf", width = 12, height = 10)
pdf(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Residence_times_all_models_ggplot.pdf", width = 12, height = 10)

ggplot_residence_times_per_areas_perc <- ggplot(data = residence_times_all_models_df_for_plot) +
  
  # Plot bars
  geom_bar(mapping = aes(x = range, y = residence_time_perc, fill = state),
           position = "stack", stat = "identity", width = 0.8) +
  
  # Add percentages as text
  geom_text(mapping = aes(y = residence_time_perc, x = range, group = regime,
                          label = paste0(round(residence_time_perc_per_range, 1), " %")),
            position = position_stack(vjust = 0.5), stat = "identity",
            col = "black", size = 4.5, fontface = "bold",
            alpha = 1.0) +
  
  # Set fill colors
  scale_fill_manual("States", values = color_scheme_states) +
  
  # Make a faceted plot
  facet_wrap(facets = ~ model,
             # scales = "fixed",
             scales = "free_x",
             nrow = 2, ncol = 3) +
  
  # Add AICc and Akaike's weights
  geom_text(data = GeoSSE_models_comparison_df,
            mapping = aes(label = paste0("AICc = ",round(AICc, 0),"\n",
                                         "Ak.w = ", round(Akaike_weights, 1), "%")),
            y = 54, x = 3.5, size = 5,
            hjust = 1, fontface = "plain") +
  
  # Set labels
  ggtitle(label = "Residence times per ranges") +
  xlab(label = "Ranges") +
  ylab(label = "Residence time  [%]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        plot.tag = element_text(size = 16, hjust = 1, face = "bold", color = "black"),
        plot.tag.position = c(0.98, 0.77),
        legend.title = element_text(size  = 20, margin = margin(b = 8)), 
        legend.text = element_text(size = 15),
        legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = unit(1.8, "line"),
        legend.spacing.y = unit(1.0, "line"),
        strip.text.x = element_text(size = 14, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(ggplot_residence_times_per_areas_perc)

dev.off()



##### 6/ Explore parameter estimates and uncertainties #####

# Explore the likelihood surface using point sampling
# to find the threshold values where different parameter estimates leads to a Delta lnLk > 2
# Such values are considered equivalent to a significant test, thus represents the boundaries of the confidence interval of parameters
# Output define the boundaries of the "supporting region" of the likelihood surface where Delta LnLk < 2
# The smaller is this region (i.e., the confidence intervals), the better are the parameter estimates

?hisse::SupportRegionGeoSSE()

# Load model outputs
# CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID1_MLE.rds")
# GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID3_MLE.rds")
# GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoHiSSE_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/CID6_MLE.rds")

CID1_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID1_MLE.rds")
GeoSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_MLE.rds")
CID2_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID2_MLE.rds")
# CID3_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID3_MLE.rds")
GeoHiSSE_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoHiSSE_MLE.rds")
CID4_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID4_MLE.rds")
# CID6_MLE <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/CID6_MLE.rds")

# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE)
# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE)
model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID2_MLE, GeoHiSSE_MLE, CID4_MLE)
# model_fit_list <- list(CID1_MLE, GeoSSE_MLE, CID3_MLE, GeoHiSSE_MLE, CID6_MLE)

# model_list <- c("CID-1", "GeoSSE", "CID-3")
# model_list <- c("CID-1", "GeoSSE", "CID-3", "GeoHiSSE")
model_list <- c("CID-1", "GeoSSE", "CID-2", "GeoHiSSE", "CID-4")
# model_list <- c("CID-1", "GeoSSE", "CID-3", "GeoHiSSE", "CID-6")


## Loop per model
GeoSSE_pars_all_models_df <- data.frame()
for (i in seq_along(model_list))
# for (i in 5:5)
{
  # i <- 1
  
  model_i <- model_list[i]
  model_MLE_i <- model_fit_list[[i]]
  
  # ### 6.1/ Extract parameters CI ####
  # 
  # ## 6.1.1/ Sample alternative parameter sets from the support region ####
  # 
  # model_MLE_CI_samples_i <- hisse::SupportRegionGeoSSE(geohisse.obj = model_MLE_i,
  #                                                      scale.int = 0.01, # Interval around MLE to sample from
  #                                                      # If most points fall within the desired.delta, need to increase the sampling interval
  #                                                      # As it would indicate the likelihood surface is very flat and CI are larger than the sampled region
  #                                                      n.points = 1000, # Number of parameter sets sampled around the optimal point in the likelihood surface
  #                                                      min.number.points = 10, # Minimum number of points meeting the criteria to find.
  #                                                      # If this number is not reached, likely because the likelihood surface is steep, will sample more points
  #                                                      # Alternative = reduce the scale.int to sample in a smaller region
  #                                                      desired.delta = 2) # Use Delta loglk = 2 as a threshold
  # 
  # str(model_MLE_CI_samples_i)
  # 
  # # Check the number of points selected is not too high (indicating the sampling region is smaller than the support region)
  # perc_samples_in_support_region <- round(nrow(model_MLE_CI_samples_i$points.within.region) / (nrow(model_MLE_CI_samples_i$all.points) - 1) * 100, 1)
  # cat(paste0("Proportion of samples in the support region for ",model_i," = ", perc_samples_in_support_region, " %\n"))
  # 
  # # Save samples from the support region
  # # saveRDS(object = model_MLE_CI_samples_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/MLE_CI_samples_", model_i, ".rds"))
  # saveRDS(object = model_MLE_CI_samples_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/MLE_CI_samples_", model_i, ".rds"))
  
  ## 6.1.2/ Compute lambda/mu ####

  # Read samples from the support region
  # model_MLE_CI_samples_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/MLE_CI_samples_", model_i, ".rds"))
  model_MLE_CI_samples_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/MLE_CI_samples_", model_i, ".rds"))

  # Extract turnover pars
  turnover_indices <- str_detect(colnames(model_MLE_CI_samples_i$points.within.region), pattern = "tau")
  turnover_samples <- model_MLE_CI_samples_i$points.within.region[, turnover_indices]

  # Extract extinct fraction pars
  ef_indices <- str_detect(colnames(model_MLE_CI_samples_i$points.within.region), pattern = "ef")
  ef_samples <- model_MLE_CI_samples_i$points.within.region[, ef_indices]

  # Match extinct fraction pars with turnover pars
  turnover_names <- str_remove(string = colnames(turnover_samples), pattern = "tau")
  ef_names <- str_remove(string = colnames(ef_samples), pattern = "ef")
  ef_samples_updated <- ef_samples[, match(x = turnover_names, table = ef_names)]
  colnames(ef_samples_updated) <- paste0("ef", turnover_names)
  ef_samples_updated[is.na(ef_samples_updated)] <- 0

  # Compute lambda and mu in samples
  SSE_rates_sampled_i <- data.frame()
  for (j in 1:nrow(model_MLE_CI_samples_i$points.within.region))
  {
    # j <- 1

    SSE_rates_j <- convert_SSE_rates(turnover = turnover_samples[j,], epsilon = ef_samples_updated[j,])
    SSE_rates_lambda_j <- SSE_rates_j$speciation
    names(SSE_rates_lambda_j) <- paste0("lambda", turnover_names)
    SSE_rates_mu_j <- SSE_rates_j$extinction
    names(SSE_rates_mu_j) <- paste0("mu", turnover_names)
    SSE_rates_j <- c(SSE_rates_lambda_j, SSE_rates_mu_j)

    # Merge in df for all samples
    SSE_rates_sampled_i <- rbind(SSE_rates_sampled_i, SSE_rates_j)
    names(SSE_rates_sampled_i) <- names(SSE_rates_j)
  }

  # Merge with sample parameters
  model_MLE_CI_samples_i$points.within.region <- cbind(model_MLE_CI_samples_i$points.within.region, SSE_rates_sampled_i)

  # Compute lambda and mu for MLE
  MLE_div_rates_i <- convert_rates_from_GeoHiSSE_output(model_MLE_i)
  paste0(colnames(MLE_div_rates_i), row.names(MLE_div_rates_i))
  states_names <- str_remove(string = row.names(MLE_div_rates_i), pattern = "\\(")
  states_names <- str_remove(string = states_names, pattern = "\\)")
  div_par_names <- as.vector(outer(X = states_names, Y = colnames(MLE_div_rates_i), FUN = function (x, y) { paste0(y,x) } ))
  if (!model_MLE_i$hidden.states) { div_par_names <- paste0(div_par_names, "A") }
  div_pars <- c(MLE_div_rates_i[,1], MLE_div_rates_i[,2])
  names(div_pars) <- div_par_names
  # Add lambda and mu to MLE pars
  MLE_i <- model_MLE_i$solution[model_MLE_i$solution != 0]
  MLE_i <- c(MLE_i, div_pars)

  # Add MLE to the samples before extracted range (needed to prevent MLE to be outside CI...)
  MLE_to_add_to_samples_i <- MLE_i[match(names(model_MLE_CI_samples_i$points.within.region), table = names(MLE_i))]
  names(MLE_to_add_to_samples_i) <- names(model_MLE_CI_samples_i$points.within.region)
  pars_samples <- rbind(MLE_to_add_to_samples_i, model_MLE_CI_samples_i$points.within.region)

  ## 6.1.3/ Extract range of sampled parameters in the support area as confidence intervals ####

  # Get parameter ranges from samples
  CI_range_i <- apply(X = pars_samples[,-1], MARGIN = 2, FUN = range, na.rm = T)
  row.names(CI_range_i) <- c("lower_CI", "upper_CI")
  # Subset only free parameters
  CI_range_i <- CI_range_i[, (CI_range_i[2,] != 0)]

  # Extract only matching pars from MLE
  MLE_i <- MLE_i[match(x = colnames(CI_range_i), table = names(MLE_i))]

  # Merge with estimates
  pars_df_i <- rbind(MLE_i, CI_range_i)

  ## Rename parameters according to model to account for equality

  # Transitions parameters
  # 2 range-shifts: 1/ jd1 (d11_00), 3/ jd0 (d00_11)
  # 2 extirpations: 2/ x1 (d01_00), 4/ x0 (d01_11)
  # 2 range extensions: 5/ d0 (d00_01), 6/ d1 (d11_01)
  # 1 regime transition: 7/ regime transition rates between A:F

  # Diversification parameters: 3 Observed states x 6 hidden states
  # 18 Turnover rates
  # 18 Extinction fraction rates
  # 18 Speciation rates
  # 18 Extinction rates

  # Accounting for equality in rates
  # For CID1: tau
  # For GeoSSE: tau0, tau1, tau01
  # For CID2: tauA, tauB
  ## For CID3: tauA, tauB, tauC
  # For GeoHiSSE: tau0A, tau0B, tau1A, tau1B, tau01A, tau01B
  # For CID4: tauA, tauB, tauC, tauD
  ## For CID6: tauA, tauB, tauC, tauD, tauE, tauF

  if (model_i == "CID-1")
  {
    pars_df_i <- pars_df_i %>%
      as.data.frame() %>%
      mutate(# lbd = lambda00A,        # Unique speciation rate for endemic areas
             lbd0 = lambda00A,       # Speciation in endemic area 0
             lbd1 = lambda11A,       # Speciation in endemic area 1
             # lbd01 = lambda01A,      # Speciation in wide-spread area 01 = Vicariance (as wide-spread speciation is not allowed)
             v01 = lambda01A,      # Speciation in wide-spread area 01 = Vicariance (as wide-spread speciation is not allowed)
             # mu = mu00A,             # Unique extinction rate for endemic areas
             mu0 = mu00A,            # Extinction in endemic area 0
             mu1 = mu11A,            # Extinction in endemic area 1
             # tau = tau00A,           # Unique turnover rate for all areas
             tau0 = tau00A,          # Turnover in endemic area 0
             tau1 = tau11A,          # Turnover in endemic area 1
             # tau01 = tau01A,         # Turnover in wide-spread area 01 = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
             # ef = ef00A,             # Unique effective fraction for endemic areas
             ef0 = ef00A,            # Effective fraction in endemic area 0
             ef1 = ef11A,            # Effective fraction in endemic area 1
             rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
             rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
             d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
             d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
             x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
             x01_1 = d01A_11A)       # Range extirpation from wide-spread area 01 to endemic area 1
    pars_df_i <- pars_df_i %>%
      # dplyr::select(lbd, lbd0, lbd1, lbd01, mu, mu0, mu1, tau, tau0, tau1, tau01, ef, ef0, ef1, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1)
      # dplyr::select(lbd0, lbd1, lbd01, mu0, mu1, tau0, tau1, tau01, ef0, ef1, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1)
      dplyr::select(lbd0, lbd1, v01, mu0, mu1, tau0, tau1, ef0, ef1, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1)

  }

  if (model_i == "GeoSSE")
  {
    pars_df_i <- pars_df_i %>%
      as.data.frame() %>%
      mutate(lbd0 = lambda00A,       # Speciation in endemic area 0
             lbd1 = lambda11A,       # Speciation in endemic area 1
             # lbd01 = lambda01A,      # Speciation in wide-spread area 01 = Vicariance (as wide-spread speciation is not allowed)
             v01 = lambda01A,      # Speciation in wide-spread area 01 = Vicariance (as wide-spread speciation is not allowed)
             mu0 = mu00A,            # Extinction in endemic area 0
             mu1 = mu11A,            # Extinction in endemic area 1
             tau0 = tau00A,          # Turnover in endemic area 0
             tau1 = tau11A,          # Turnover in endemic area 1
             # tau01 = tau01A,         # Turnover in wide-spread area 01 = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
             ef0 = ef00A,            # Extinct fraction in endemic area 0
             ef1 = ef11A,            # Extinct fraction in endemic area 1
             rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
             rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
             d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
             d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
             x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
             x01_1 = d01A_11A)       # Range extirpation from wide-spread area 01 to endemic area 1
    pars_df_i <- pars_df_i %>%
      # dplyr::select(lbd0, lbd1, lbd01, mu0, mu1, tau0, tau1, tau01, ef0, ef1, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1)
      dplyr::select(lbd0, lbd1, v01, mu0, mu1, tau0, tau1, ef0, ef1, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1)
  }

  if (model_i == "CID-2")
  {
    pars_df_i <- pars_df_i %>%
      as.data.frame() %>%
      mutate(# lbdA = lambda00A,       # Unique speciation for endemic areas in regime A
        # lbdB = lambda00B,       # Unique speciation for endemic areas in regime B
        lbd0A = lambda00A,      # Speciation for endemic area 0 in regime A
        lbd0B = lambda00B,      # Speciation for endemic area 0 in regime B
        lbd1A = lambda11A,      # Speciation for endemic area 1 in regime A
        lbd1B = lambda11B,      # Speciation for endemic area 1 in regime B
        # lbd01A = lambda01A,     # Speciation in wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
        # lbd01B = lambda01B,     # Speciation in wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
        v01A = lambda01A,     # Speciation in wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
        v01B = lambda01B,     # Speciation in wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
        # muA = mu00A,            # Unique extinction for endemic areas in regime A
        # muB = mu00B,            # Unique extinction for endemic areas in regime B
        mu0A = mu00A,           # Extinction for endemic area 0 in regime A
        mu0B = mu00B,           # Extinction for endemic area 0 in regime B
        mu1A = mu11A,           # Extinction for endemic area 1 in regime A
        mu1B = mu11B,           # Extinction for endemic area 1 in regime B
        # tauA = tau00A,          # Unique turnover rate for all areas in regime A
        # tauB = tau00B,          # Unique turnover rate for all areas in regime B
        tau0A = tau00A,         # Turnover rate for endemic area 0 in regime A
        tau0B = tau00B,         # Turnover rate for endemic area 0 in regime B
        tau1A = tau11A,         # Turnover rate for endemic area 1 in regime A
        tau1B = tau11B,         # Turnover rate for endemic area 1 in regime B
        # tau01A = tau01A,        # Turnover rate for wide-spread area 01 in regime A = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # tau01B = tau01B,        # Turnover rate for wide-spread area 01 in regime B = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # efA = ef00A,            # Unique effective fraction for endemic areas in regime A
        # efB = ef00B,            # Unique effective fraction for endemic areas in regime B
        ef0A = ef00A,           # Extinct fraction for endemic area 0 in regime A
        ef0B = ef00B,           # Extinct fraction for endemic area 0 in regime B
        ef1A = ef11A,           # Extinct fraction for endemic area 1 in regime A
        ef1B = ef11B,           # Extinct fraction for endemic area 1 in regime B
        rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
        rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
        d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
        d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
        x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
        x01_1 = d01A_11A,       # Range extirpation from wide-spread area 01 to endemic area 1
        regime_shift = d00A_00B)  # Unique regime shift rate
    pars_df_i <- pars_df_i %>%
      # dplyr::select(lbdA, lbdB, lbd01A, lbd01B, muA, muB, tauA, tauB, efA, efB, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
      # dplyr::select(lbd0A, lbd0B, lbd1A, lbd1B, lbd01A, lbd01B, mu0A, mu0B, mu1A, mu1B, tau0A, tau0B, tau1A, tau1B, tau01A, tau01B, ef0A, ef0B, ef1A, ef1B, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
      dplyr::select(lbd0A, lbd0B, lbd1A, lbd1B, v01A, v01B, mu0A, mu0B, mu1A, mu1B, tau0A, tau0B, tau1A, tau1B, ef0A, ef0B, ef1A, ef1B, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  }

  # if (model_i == "CID-3")
  # {
  #   pars_df_i <- pars_df_i %>%
  #     as.data.frame() %>%
  #     mutate(# lbdA = lambda00A,       # Unique speciation for endemic areas in regime A
  #            # lbdB = lambda00B,       # Unique speciation for endemic areas in regime B
  #            # lbdC = lambda00C,       # Unique speciation for endemic areas in regime C
  #            lbd0A = lambda00A,      # Speciation for endemic area 0 in regime A
  #            lbd0B = lambda00B,      # Speciation for endemic area 0 in regime B
  #            lbd0C = lambda00C,      # Speciation for endemic area 0 in regime C
  #            lbd1A = lambda11A,      # Speciation for endemic area 1 in regime A
  #            lbd1B = lambda11B,      # Speciation for endemic area 1 in regime B
  #            lbd1C = lambda11C,      # Speciation for endemic area 1 in regime C
  #            lbd01A = lambda01A,     # Speciation in wide-spread area 01 in regime A
  #            lbd01B = lambda01B,     # Speciation in wide-spread area 01 in regime B
  #            lbd01C = lambda01C,     # Speciation in wide-spread area 01 in regime C
  #            # muA = mu00A,            # Unique extinction for endemic areas in regime A
  #            # muB = mu00B,            # Unique extinction for endemic areas in regime B
  #            # muC = mu00C,            # Unique extinction for endemic areas in regime C
  #            mu0A = mu00A,           # Extinction for endemic area 0 in regime A
  #            mu0B = mu00B,           # Extinction for endemic area 0 in regime B
  #            mu0C = mu00C,           # Extinction for endemic area 0 in regime C
  #            mu1A = mu11A,           # Extinction for endemic area 1 in regime A
  #            mu1B = mu11B,           # Extinction for endemic area 1 in regime B
  #            mu1C = mu11C,           # Extinction for endemic area 1 in regime C
  #            # tauA = tau00A,          # Unique turnover rate for all areas in regime A
  #            # tauB = tau00B,          # Unique turnover rate for all areas in regime B
  #            # tauC = tau00C,          # Unique turnover rate for all areas in regime C
  #            tau0A = tau00A,         # Turnover rate for endemic area 0 in regime A
  #            tau0B = tau00B,         # Turnover rate for endemic area 0 in regime B
  #            tau0C = tau00C,         # Turnover rate for endemic area 0 in regime C
  #            tau1A = tau11A,         # Turnover rate for endemic area 1 in regime A
  #            tau1B = tau11B,         # Turnover rate for endemic area 1 in regime B
  #            tau1C = tau11C,         # Turnover rate for endemic area 1 in regime C
  #            tau01A = tau01A,        # Turnover rate for wide-spread area 01 in regime A
  #            tau01B = tau01B,        # Turnover rate for wide-spread area 01 in regime B
  #            tau01C = tau01C,        # Turnover rate for wide-spread area 01 in regime C
  #            # efA = ef00A,            # Unique effective fraction for endemic areas in regime A
  #            # efB = ef00B,            # Unique effective fraction for endemic areas in regime B
  #            # efC = ef00C,            # Unique effective fraction for endemic areas in regime C
  #            ef0A = ef00A,           # Extinct fraction for endemic area 0 in regime A
  #            ef0B = ef00B,           # Extinct fraction for endemic area 0 in regime B
  #            ef0C = ef00C,           # Extinct fraction for endemic area 0 in regime C
  #            ef1A = ef11A,           # Extinct fraction for endemic area 1 in regime A
  #            ef1B = ef11B,           # Extinct fraction for endemic area 1 in regime B
  #            ef1C = ef11C,           # Extinct fraction for endemic area 1 in regime C
  #            rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
  #            rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
  #            d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
  #            d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
  #            x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
  #            x01_1 = d01A_11A,       # Range extirpation from wide-spread area 01 to endemic area 1
  #            regime_shift = d00A_00B)  # Unique regime shift rate
  #   pars_df_i <- pars_df_i %>%
  #     # dplyr::select(lbdA, lbdB, lbdC, lbd01A, lbd01B, lbd01C, muA, muB, muC, tauA, tauB, tauC, efA, efB, efC, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  #     dplyr::select(lbd0A, lbd0B, lbd0C, lbd1A, lbd1B, lbd1C, lbd01A, lbd01B, lbd01C, mu0A, mu0B, mu0C, mu1A, mu1B, mu1C, tau0A, tau0B, tau0C, tau1A, tau1B, tau1C, tau01A, tau01B, tau01C, ef0A, ef0B, ef0C, ef1A, ef1B, ef1C, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  # }

  if (model_i == "GeoHiSSE")
  {
    pars_df_i <- pars_df_i %>%
      as.data.frame() %>%
      mutate(lbd0A = lambda00A,      # Speciation for endemic area 0 in regime A
             lbd0B = lambda00B,      # Speciation for endemic area 0 in regime B
             lbd1A = lambda11A,      # Speciation for endemic area 1 in regime A
             lbd1B = lambda11B,      # Speciation for endemic area 1 in regime B
             # lbd01A = lambda01A,     # Speciation for wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
             # lbd01B = lambda01B,     # Speciation for wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
             v01A = lambda01A,     # Speciation for wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
             v01B = lambda01B,     # Speciation for wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
             mu0A = mu00A,           # Extinction for endemic area 0 in regime A
             mu0B = mu00B,           # Extinction for endemic area 0 in regime B
             mu1A = mu11A,           # Extinction for endemic area 1 in regime A
             mu1B = mu11B,           # Extinction for endemic area 1 in regime B
             tau0A = tau00A,         # Turnover for endemic area 0 in regime A
             tau0B = tau00B,         # Turnover for endemic area 0 in regime B
             tau1A = tau11A,         # Turnover for endemic area 1 in regime A
             tau1B = tau11B,         # Turnover for endemic area 1 in regime B
             # tau01A = tau01A,        # Turnover for wide-spread area 01 in regime A = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
             # tau01B = tau01B,        # Turnover for wide-spread area 01 in regime B = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
             ef0A = ef00A,           # Extinct fraction for endemic area 0 in regime A
             ef0B = ef00B,           # Extinct fraction for endemic area 0 in regime B
             ef1A = ef11A,           # Extinct fraction for endemic area 1 in regime A
             ef1B = ef11B,           # Extinct fraction for endemic area 1 in regime B
             rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
             rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
             d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
             d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
             x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
             x01_1 = d01A_11A,       # Range extirpation from wide-spread area 01 to endemic area 1
             regime_shift = d00A_00B)  # Unique regime shift rate
    pars_df_i <- pars_df_i %>%
      # dplyr::select(lbd0A, lbd0B, lbd1A, lbd1B, lbd01A, lbd01B, mu0A, mu0B, mu1A, mu1B, tau0A, tau0B, tau1A, tau1B, tau01A, tau01B, ef0A, ef0B, ef1A, ef1B, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
      dplyr::select(lbd0A, lbd0B, lbd1A, lbd1B, v01A, v01B, mu0A, mu0B, mu1A, mu1B, tau0A, tau0B, tau1A, tau1B, ef0A, ef0B, ef1A, ef1B, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  }

  if (model_i == "CID-4")
  {
    pars_df_i <- pars_df_i %>%
      as.data.frame() %>%
      mutate(# lbdA = lambda00A,       # Unique speciation for endemic areas in regime A
        # lbdB = lambda00B,       # Unique speciation for endemic areas in regime B
        # lbdC = lambda00C,       # Unique speciation for endemic areas in regime C
        # lbdD = lambda00D,       # Unique speciation for endemic areas in regime D
        lbd0A = lambda00A,      # Speciation for endemic area 0 in regime A
        lbd0B = lambda00B,      # Speciation for endemic area 0 in regime B
        lbd0C = lambda00C,      # Speciation for endemic area 0 in regime C
        lbd0D = lambda00D,      # Speciation for endemic area 0 in regime D
        lbd1A = lambda11A,      # Speciation for endemic area 1 in regime A
        lbd1B = lambda11B,      # Speciation for endemic area 1 in regime B
        lbd1C = lambda11C,      # Speciation for endemic area 1 in regime C
        lbd1D = lambda11D,      # Speciation for endemic area 1 in regime D
        # lbd01A = lambda01A,     # Speciation in wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
        # lbd01B = lambda01B,     # Speciation in wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
        # lbd01C = lambda01C,     # Speciation in wide-spread area 01 in regime C = Vicariance (as wide-spread speciation is not allowed)
        # lbd01D = lambda01D,     # Speciation in wide-spread area 01 in regime D = Vicariance (as wide-spread speciation is not allowed)
        v01A = lambda01A,     # Speciation in wide-spread area 01 in regime A = Vicariance (as wide-spread speciation is not allowed)
        v01B = lambda01B,     # Speciation in wide-spread area 01 in regime B = Vicariance (as wide-spread speciation is not allowed)
        v01C = lambda01C,     # Speciation in wide-spread area 01 in regime C = Vicariance (as wide-spread speciation is not allowed)
        v01D = lambda01D,     # Speciation in wide-spread area 01 in regime D = Vicariance (as wide-spread speciation is not allowed)
        # muA = mu00A,            # Unique extinction for endemic areas in regime A
        # muB = mu00B,            # Unique extinction for endemic areas in regime B
        # muC = mu00C,            # Unique extinction for endemic areas in regime C
        # muD = mu00D,            # Unique extinction for endemic areas in regime D
        mu0A = mu00A,           # Extinction for endemic area 0 in regime A
        mu0B = mu00B,           # Extinction for endemic area 0 in regime B
        mu0C = mu00C,           # Extinction for endemic area 0 in regime C
        mu0D = mu00D,           # Extinction for endemic area 0 in regime D
        mu1A = mu11A,           # Extinction for endemic area 1 in regime A
        mu1B = mu11B,           # Extinction for endemic area 1 in regime B
        mu1C = mu11C,           # Extinction for endemic area 1 in regime C
        mu1D = mu11D,           # Extinction for endemic area 1 in regime D
        # tauA = tau00A,          # Unique turnover rate for all areas in regime A
        # tauB = tau00B,          # Unique turnover rate for all areas in regime B
        # tauC = tau00C,          # Unique turnover rate for all areas in regime C
        # tauD = tau00D,          # Unique turnover rate for all areas in regime D
        tau0A = tau00A,         # Turnover for endemic area 0 in regime A
        tau0B = tau00B,         # Turnover for endemic area 0 in regime B
        tau0C = tau00C,         # Turnover for endemic area 0 in regime C
        tau0D = tau00D,         # Turnover for endemic area 0 in regime D
        tau1A = tau11A,         # Turnover for endemic area 1 in regime A
        tau1B = tau11B,         # Turnover for endemic area 1 in regime B
        tau1C = tau11C,         # Turnover for endemic area 1 in regime C
        tau1D = tau11D,         # Turnover for endemic area 1 in regime D
        # tau01A = tau01A,        # Turnover for wide-spread area 01 in regime A = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # tau01B = tau01B,        # Turnover for wide-spread area 01 in regime B = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # tau01C = tau01C,        # Turnover for wide-spread area 01 in regime C = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # tau01D = tau01D,        # Turnover for wide-spread area 01 in regime D = Vicariance (as extinction from wide-spread state AND wide-spread speciation is not allowed)
        # efA = ef00A,            # Unique extinct fraction for endemic areas in regime A
        # efB = ef00B,            # Unique extinct fraction for endemic areas in regime B
        # efC = ef00C,            # Unique extinct fraction for endemic areas in regime C
        # efD = ef00D,            # Unique extinct fraction for endemic areas in regime D
        ef0A = ef00A,           # Extinct fraction for endemic area 0 in regime A
        ef0B = ef00B,           # Extinct fraction for endemic area 0 in regime B
        ef0C = ef00C,           # Extinct fraction for endemic area 0 in regime C
        ef0D = ef00D,           # Extinct fraction for endemic area 0 in regime D
        ef1A = ef11A,           # Extinct fraction for endemic area 1 in regime A
        ef1B = ef11B,           # Extinct fraction for endemic area 1 in regime B
        ef1C = ef11C,           # Extinct fraction for endemic area 1 in regime C
        ef1D = ef11D,           # Extinct fraction for endemic area 1 in regime D
        rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
        rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
        d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
        d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
        x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
        x01_1 = d01A_11A,       # Range extirpation from wide-spread area 01 to endemic area 1
        regime_shift = d00A_00B)  # Unique regime shift rate
    pars_df_i <- pars_df_i %>%
      # dplyr::select(lbdA, lbdB, lbdC, lbdD, lbd01A, lbd01B, lbd01C, lbd01D, muA, muB, muC, muD, tauA, tauB, tauC, tauD, efA, efB, efC, efD, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
      # dplyr::select(lbd0A, lbd0B, lbd0C, lbd0D, lbd1A, lbd1B, lbd1C, lbd1D, lbd01A, lbd01B, lbd01C, lbd01D, mu0A, mu0B, mu0C, mu0D, mu1A, mu1B, mu1C, mu1D, tau0A, tau0B, tau0C, tau0D, tau1A, tau1B, tau1C, tau1D, tau01A, tau01B, tau01C, tau01D, ef0A, ef0B, ef0C, ef0D, ef1A, ef1B, ef1C, ef1D, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
      dplyr::select(lbd0A, lbd0B, lbd0C, lbd0D, lbd1A, lbd1B, lbd1C, lbd1D, v01A, v01B, v01C, v01D, mu0A, mu0B, mu0C, mu0D, mu1A, mu1B, mu1C, mu1D, tau0A, tau0B, tau0C, tau0D, tau1A, tau1B, tau1C, tau1D, ef0A, ef0B, ef0C, ef0D, ef1A, ef1B, ef1C, ef1D, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  }

  # if (model_i == "CID-6")
  # {
  #   pars_df_i <- pars_df_i %>%
  #     as.data.frame() %>%
  #     mutate(# lbdA = lambda00A,       # Unique speciation for endemic areas in regime A
  #            # lbdB = lambda00B,       # Unique speciation for endemic areas in regime B
  #            # lbdC = lambda00C,       # Unique speciation for endemic areas in regime C
  #            # lbdD = lambda00D,       # Unique speciation for endemic areas in regime D
  #            # lbdE = lambda00E,       # Unique speciation for endemic areas in regime E
  #            # lbdF = lambda00F,       # Unique speciation for endemic areas in regime F
  #            lbd0A = lambda00A,      # Speciation for endemic area 0 in regime A
  #            lbd0B = lambda00B,      # Speciation for endemic area 0 in regime B
  #            lbd0C = lambda00C,      # Speciation for endemic area 0 in regime C
  #            lbd0D = lambda00D,      # Speciation for endemic area 0 in regime D
  #            lbd0E = lambda00E,      # Speciation for endemic area 0 in regime E
  #            lbd0F = lambda00F,      # Speciation for endemic area 0 in regime F
  #            lbd1A = lambda11A,      # Speciation for endemic area 1 in regime A
  #            lbd1B = lambda11B,      # Speciation for endemic area 1 in regime B
  #            lbd1C = lambda11C,      # Speciation for endemic area 1 in regime C
  #            lbd1D = lambda11D,      # Speciation for endemic area 1 in regime D
  #            lbd1E = lambda11E,      # Speciation for endemic area 1 in regime E
  #            lbd1F = lambda11F,      # Speciation for endemic area 1 in regime F
  #            lbd01A = lambda01A,     # Speciation in wide-spread area 01 in regime A
  #            lbd01B = lambda01B,     # Speciation in wide-spread area 01 in regime B
  #            lbd01C = lambda01C,     # Speciation in wide-spread area 01 in regime C
  #            lbd01D = lambda01D,     # Speciation in wide-spread area 01 in regime D
  #            lbd01E = lambda01E,     # Speciation in wide-spread area 01 in regime E
  #            lbd01F = lambda01F,     # Speciation in wide-spread area 01 in regime F
  #            # muA = mu00A,            # Unique extinction for endemic areas in regime A
  #            # muB = mu00B,            # Unique extinction for endemic areas in regime B
  #            # muC = mu00C,            # Unique extinction for endemic areas in regime C
  #            # muD = mu00D,            # Unique extinction for endemic areas in regime D
  #            # muE = mu00E,            # Unique extinction for endemic areas in regime E
  #            # muF = mu00F,            # Unique extinction for endemic areas in regime F
  #            mu0A = mu00A,           # Extinction for endemic area 0 in regime A
  #            mu0B = mu00B,           # Extinction for endemic area 0 in regime B
  #            mu0C = mu00C,           # Extinction for endemic area 0 in regime C
  #            mu0D = mu00D,           # Extinction for endemic area 0 in regime D
  #            mu0E = mu00E,           # Extinction for endemic area 0 in regime E
  #            mu0F = mu00F,           # Extinction for endemic area 0 in regime F
  #            mu1A = mu11A,           # Extinction for endemic area 1 in regime A
  #            mu1B = mu11B,           # Extinction for endemic area 1 in regime B
  #            mu1C = mu11C,           # Extinction for endemic area 1 in regime C
  #            mu1D = mu11D,           # Extinction for endemic area 1 in regime D
  #            mu1E = mu11E,           # Extinction for endemic area 1 in regime E
  #            mu1F = mu11F,           # Extinction for endemic area 1 in regime F
  #            # tauA = tau00A,          # Unique turnover rate for all areas in regime A
  #            # tauB = tau00B,          # Unique turnover rate for all areas in regime B
  #            # tauC = tau00C,          # Unique turnover rate for all areas in regime C
  #            # tauD = tau00D,          # Unique turnover rate for all areas in regime D
  #            # tauE = tau00E,          # Unique turnover rate for all areas in regime E
  #            # tauF = tau00F,          # Unique turnover rate for all areas in regime F
  #            tau0A = tau00A,         # Turnover for endemic area 0 in regime A
  #            tau0B = tau00B,         # Turnover for endemic area 0 in regime B
  #            tau0C = tau00C,         # Turnover for endemic area 0 in regime C
  #            tau0D = tau00D,         # Turnover for endemic area 0 in regime D
  #            tau0E = tau00E,         # Turnover for endemic area 0 in regime E
  #            tau0F = tau00F,         # Turnover for endemic area 0 in regime F
  #            tau1A = tau11A,         # Turnover for endemic area 1 in regime A
  #            tau1B = tau11B,         # Turnover for endemic area 1 in regime B
  #            tau1C = tau11C,         # Turnover for endemic area 1 in regime C
  #            tau1D = tau11D,         # Turnover for endemic area 1 in regime D
  #            tau1E = tau11E,         # Turnover for endemic area 1 in regime E
  #            tau1F = tau11F,         # Turnover for endemic area 1 in regime F
  #            tau01A = tau01A,        # Turnover for wide-spread area 01 in regime A
  #            tau01B = tau01B,        # Turnover for wide-spread area 01 in regime B
  #            tau01C = tau01C,        # Turnover for wide-spread area 01 in regime C
  #            tau01D = tau01D,        # Turnover for wide-spread area 01 in regime D
  #            tau01E = tau01E,        # Turnover for wide-spread area 01 in regime E
  #            tau01F = tau01F,        # Turnover for wide-spread area 01 in regime F
  #            # efA = ef00A,            # Unique extinct fraction for endemic areas in regime A
  #            # efB = ef00B,            # Unique extinct fraction for endemic areas in regime B
  #            # efC = ef00C,            # Unique extinct fraction for endemic areas in regime C
  #            # efD = ef00D,            # Unique extinct fraction for endemic areas in regime D
  #            # efE = ef00E,            # Unique extinct fraction for endemic areas in regime E
  #            # efF = ef00F,            # Unique extinct fraction for endemic areas in regime F
  #            ef0A = ef00A,           # Extinct fraction for endemic area 0 in regime A
  #            ef0B = ef00B,           # Extinct fraction for endemic area 0 in regime B
  #            ef0C = ef00C,           # Extinct fraction for endemic area 0 in regime C
  #            ef0D = ef00D,           # Extinct fraction for endemic area 0 in regime D
  #            ef0E = ef00E,           # Extinct fraction for endemic area 0 in regime E
  #            ef0F = ef00F,           # Extinct fraction for endemic area 0 in regime F
  #            ef1A = ef11A,           # Extinct fraction for endemic area 1 in regime A
  #            ef1B = ef11B,           # Extinct fraction for endemic area 1 in regime B
  #            ef1C = ef11C,           # Extinct fraction for endemic area 1 in regime C
  #            ef1D = ef11D,           # Extinct fraction for endemic area 1 in regime D
  #            ef1E = ef11E,           # Extinct fraction for endemic area 1 in regime E
  #            ef1F = ef11F,           # Extinct fraction for endemic area 1 in regime F
  #            rs0_1 = d00A_11A,       # Range-shift (jump) from endemic area 0 to endemic area 1
  #            rs1_0 = d11A_00A,       # Range-shift (jump) from endemic area 1 to endemic area 0
  #            d0_01 = d00A_01A,       # Range extension from endemic area 0 to wide-spread area 01
  #            d1_01 = d11A_01A,       # Range extension from endemic area 1 to wide-spread area 01
  #            x01_0 = d01A_00A,       # Range extirpation from wide-spread area 01 to endemic area 0
  #            x01_1 = d01A_11A,       # Range extirpation from wide-spread area 01 to endemic area 1
  #            regime_shift = d00A_00B)  # Unique regime shift rate
  #   pars_df_i <- pars_df_i %>%
  #     # dplyr::select(lbdA, lbdB, lbdC, lbdD, lbdE, lbdF, lbd01A, lbd01B, lbd01C, lbd01D, lbd01E, lbd01F, muA, muB, muC, muD, muE, muF, tauA, tauB, tauC, tauD, tauE, tauF, efA, efB, efC, efD, efE, efF, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  #     dplyr::select(lbd0A, lbd0B, lbd0C, lbd0D, lbd0E, lbd0F, lbd1A, lbd1B, lbd1C, lbd1D, lbd1E, lbd1F, lbd01A, lbd01B, lbd01C, lbd01D, lbd01E, lbd01F, mu0A, mu0B, mu0C, mu0D, mu0E, mu0F, mu1A, mu1B, mu1C, mu1D, mu1E, mu1F, tau0A, tau0B, tau0C, tau0D, tau0E, tau0F, tau1A, tau1B, tau1C, tau1D, tau1E, tau1F, tau01A, tau01B, tau01C, tau01D, tau01E, tau01F, ef0A, ef0B, ef0C, ef0D, ef0E, ef0F, ef1A, ef1B, ef1C, ef1D, ef1E, ef1F, rs0_1, rs1_0, d0_01, d1_01, x01_0, x01_1, regime_shift)
  # }

  # Melt dataframe
  pars_df_i <- reshape2::melt(as.matrix(pars_df_i))
  names(pars_df_i) <- c("stats", "pars", "value")
  pars_df_i$model <- model_i
  pars_df_i <- pars_df_i[, c("model", "stats", "pars", "value")]

  # Merge parameter df
  GeoSSE_pars_all_models_df <- rbind(GeoSSE_pars_all_models_df, pars_df_i)
  
  ## Print progress
  cat(paste0(Sys.time(), " - Parameter CI identified for ", model_i," - n° ", i, "/", length(model_list),"\n"))
}

## Add parameter types
GeoSSE_pars_all_models_df$parameter_type <- NA
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "lbd")] <- "Speciation"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "v")] <- "Vicariance"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "mu")] <- "Extinction"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "tau")] <- "Turnover"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "ef")] <- "Extinct fraction"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "rs")] <- "Range-shift"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "^d")] <- "Range extension"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "^x")] <- "Range extirpation"
GeoSSE_pars_all_models_df$parameter_type[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "regime_shift")] <- "Regime shift"

## Add states
GeoSSE_pars_all_models_df$state <- "All"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "0")] <- "NW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "1")] <- "OW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "01")] <- "WS" # WS = Wide-Spread
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "0_1")] <- "NW to OW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "1_0")] <- "OW to NW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "0_01")] <- "To OW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "1_01")] <- "To NW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "01_0")] <- "From OW"
GeoSSE_pars_all_models_df$state[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "01_1")] <- "From NW"

table(GeoSSE_pars_all_models_df$state)

## Add regimes

GeoSSE_pars_all_models_df$regime <- ""
GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "A")] <- "A"
GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "B")] <- "B"
GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "C")] <- "C"
GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "D")] <- "D"
# GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "E")] <- "E"
# GeoSSE_pars_all_models_df$regime[str_detect(string = GeoSSE_pars_all_models_df$pars, pattern = "F")] <- "F"

table(GeoSSE_pars_all_models_df$regime)

## Save GeoSSE models parameter summary table
# saveRDS(GeoSSE_pars_all_models_df, file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")
saveRDS(GeoSSE_pars_all_models_df, file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")


### 6.2/ Plot parameter estimates ####

# Question of interest: do the values of the CI of the state-dependent diversification rates overlap?

## GGplot SE to compare rate estimates

# Load GeoSSE models parameter summary table
# GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")
GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")

# Load GeoSSE model comparison summary table
# GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.rds")
GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.rds")

## One plot per model

# Remove turnover and extinct fraction parameters
GeoSSE_pars_all_models_df_for_plot <- GeoSSE_pars_all_models_df %>% 
  filter(parameter_type %in% c("Speciation", "Extinction", "Vicariance", "Range-shift", "Range extension", "Range extirpation", "Regime shift"))

# Set color scheme for parameter types
colors_list_for_parameter_types <- c("limegreen", "coral1", "tan4", "dodgerblue", "purple", "grey20", "grey")
parameter_types_names <- c("Speciation", "Extinction", "Vicariance", "Range-shift", "Range extension", "Range extirpation", "Regime shift")
names(colors_list_for_parameter_types) <- parameter_types_names
parameter_types_labels <- c("lambda", "mu", "v", "rs", "d", "x", "q")

# Order parameter_type so they are properly ordered
GeoSSE_pars_all_models_df_for_plot$parameter_type <- factor(GeoSSE_pars_all_models_df_for_plot$parameter_type,
                                                            levels = parameter_types_names,
                                                            labels = parameter_types_names)

## Loop per model
for (i in seq_along(model_list))
{
  # i <- 1
  
  # Extract model
  model_i <- model_list[i]
  
  # Extract parameter data
  GeoSSE_pars_df_model_i <- GeoSSE_pars_all_models_df_for_plot %>% 
    filter(model == model_i) %>%
    # filter(!(state == "All" & regime == "")) %>%
    pivot_wider(names_from = stats, values_from = "value")
  
  # Extract evaluation results
  GeoSSE_models_comparison_df_i <- GeoSSE_models_comparison_df %>% 
    filter(model == model_i)
  
  # # Define axis limits (ignore outliers with values > 1)
  # ylim <- max(GeoSSE_pars_df_model_i$upper_CI[GeoSSE_pars_df_model_i$upper_CI <= 1]) * 1.1
  
  nb_par <- length(unique(as.character(GeoSSE_pars_df_model_i$pars)))
  
  ## GGplot for all parameters
  # pdf(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_all_pars_plot_",model_i,".pdf"), height = 6, width = nb_par*1.5)
  pdf(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_all_pars_plot_",model_i,".pdf"), height = 6, width = nb_par*1.5)
  
  GeoSSE_model_pars_plot_i <- ggplot(data = GeoSSE_pars_df_model_i) +
    
    # Add horizontal line: par = 0
    geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
    
    # Add error bars
    geom_errorbar(mapping = aes(x = interaction(state, regime), y = MLE_i, ymin = lower_CI, ymax = upper_CI, col = parameter_type),
                  width = 0.2, linewidth = 1.2,
                  alpha = 1.0, show.legend = F) +
    
    # Add mean points
    geom_point(mapping = aes(x = interaction(state, regime), y = MLE_i, fill = parameter_type),
               shape = 21, size = 5, col = "black", alpha = 1.0, show.legend = F) +
    
    # Add AICc and Akaike's weights
    # annotate(geom = "text", x = 1, y = max(GeoSSE_pars_df_model_i$upper_CI) * 0.95,
    #          label = paste0("AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
    #                         "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights * 100, 1)),
    #          hjust = 0, fontface = "bold", size = 5) +
    
    # Add AICc and Akaike's weights as tag
    labs(tag = paste0("AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
                      "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights, 1), " %")) +
    
    # # Convert Y-scale to log
    # scale_y_continuous(transform = "log1p",
    #                    # breaks = c(7.389, 20.085, 54.598, 148.41),
    #                    # labels = c("7", "20", "55", "150")
    #                    ) + 
    
    # Set Y-axis limits
    # scale_y_continuous(limits =  c(0, ylim)) + 
    
    
    # Make a faceted plot
    facet_wrap(facets = ~ parameter_type,
               # scales = "fixed",
               scales = "free",
               nrow = 1, ncol = 7
    ) +
    
    # Adjust color scheme and legend
    scale_color_manual("Parameter types", labels = parameter_types_labels, values = unname(colors_list_for_parameter_types)) +
    
    # Adjust fill scheme and legend
    scale_fill_manual("Parameter types", labels = parameter_types_labels, values = unname(colors_list_for_parameter_types)) +
    
    # # Adjust x-axis labels
    # scale_x_discrete(labels = c("lambda" = bquote(lambda), "mu" = bquote(mu), "v" = "v", "jd" = "jd", "d" = "d", "x" = "x", "q" = "q")) +
    
    # Flip axes
    # coord_flip() + 
    
    # Set plot title +
    ggtitle(label = paste0("Parameter estimates for ",model_i," model")) +
    
    # Set axes labels
    xlab("Parameters") +
    ylab("Rates  [Events / lineage / My]") +
    
    # Adjust aesthetics
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
          plot.tag = element_text(size = 16, hjust = 0, face = "bold", color = "black"),
          plot.tag.position = c(0.065, 0.97),
          legend.title = element_text(size  = 20, margin = margin(b = 8)), 
          legend.text = element_text(size = 15),
          legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
          legend.key.size = unit(1.8, "line"),
          legend.spacing.y = unit(1.0, "line"),
          strip.text.x = element_text(size = 14, colour = "black", face = "bold"),
          axis.title = element_text(size = 20, color = "black"),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 12)),
          axis.line = element_line(linewidth = 1.5),
          axis.ticks.length = unit(8, "pt"),
          axis.text = element_text(size = 14, color = "black"),
          axis.text.x = element_text(margin = margin(t = 5)),
          axis.text.y = element_text(margin = margin(r = 5)))
  
  # Plot
  print(GeoSSE_model_pars_plot_i)
  
  dev.off()
  
  
  ## Extract only diversification parameters
  
  # Extract parameter data
  GeoSSE_div_pars_df_model_i <- GeoSSE_pars_df_model_i %>% 
    filter(parameter_type %in% c("Speciation", "Extinction"))

  nb_par <- length(unique(as.character(GeoSSE_div_pars_df_model_i$pars)))
  
  ## GGplot for diversification parameters
  # pdf(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_div_pars_plot_",model_i,".pdf"), height = 6, width = nb_par*1.5)
  pdf(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_div_pars_plot_",model_i,".pdf"), height = 6, width = nb_par*1.5)
  
  GeoSSE_model_div_pars_plot_i <- ggplot(data = GeoSSE_div_pars_df_model_i) +
    
    # Add horizontal line: par = 0
    geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed", col = "grey") +
    
    # Add error bars
    geom_errorbar(mapping = aes(x = interaction(state, regime), y = MLE_i, ymin = lower_CI, ymax = upper_CI, col = parameter_type),
                  width = 0.2, linewidth = 1.2,
                  alpha = 1.0, show.legend = F) +
    
    # Add mean points
    geom_point(mapping = aes(x = interaction(state, regime), y = MLE_i, fill = parameter_type),
               shape = 21, size = 5, col = "black", alpha = 1.0, show.legend = F) +
    
    # Add AICc and Akaike's weights
    # annotate(geom = "text", x = 1, y = max(GeoSSE_pars_df_model_i$upper_CI) * 0.95,
    #          label = paste0("AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
    #                         "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights * 100, 1)),
    #          hjust = 0, fontface = "bold", size = 5) +
    
    # Add AICc and Akaike's weights as tag
    labs(tag = paste0("AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
                      "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights, 1), " %")) +
    
    # Convert Y-scale to log
    scale_y_continuous(transform = "log1p",
                       # breaks = c(7.389, 20.085, 54.598, 148.41),
                       # labels = c("7", "20", "55", "150")
                       ) +
    
    # Set Y-axis limits
    # scale_y_continuous(limits =  c(0, ylim)) + 
    
    
    # Make a facetted plot
    facet_wrap(facets = ~ parameter_type,
               # scales = "fixed",
               scales = "free_x",
               nrow = 1, ncol = 2) +
    
    # Adjust color scheme and legend
    scale_color_manual("Parameter types", labels = parameter_types_labels, values = unname(colors_list_for_parameter_types)) +
    
    # Adjust fill scheme and legend
    scale_fill_manual("Parameter types", labels = parameter_types_labels, values = unname(colors_list_for_parameter_types)) +
    
    # # Adjust x-axis labels
    # scale_x_discrete(labels = c("lambda" = bquote(lambda), "mu" = bquote(mu), "jd" = "jd", "d" = "d", "x" = "x", "q" = "q")) +
    
    # Flip axes
    # coord_flip() + 
    
    # Set plot title +
    ggtitle(label = paste0("Parameter estimates for ",model_i," model")) +
    
    # Set axes labels
    xlab("Parameters") +
    ylab("Rates  [Events / lineage / My]") +
    
    # Adjust aesthetics
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
          panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.5),
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
          plot.tag = element_text(size = 16, hjust = 1, face = "bold", color = "black"),
          plot.tag.position = c(0.98, 0.77),
          legend.title = element_text(size  = 20, margin = margin(b = 8)), 
          legend.text = element_text(size = 15),
          legend.key = element_rect(colour = NA, fill = NA, linewidth = 5),
          legend.key.size = unit(1.8, "line"),
          legend.spacing.y = unit(1.0, "line"),
          strip.text.x = element_text(size = 14, colour = "black", face = "bold"),
          axis.title = element_text(size = 20, color = "black"),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 12)),
          axis.line = element_line(linewidth = 1.5),
          axis.ticks.length = unit(8, "pt"),
          axis.text = element_text(size = 14, color = "black"),
          axis.text.x = element_text(margin = margin(t = 5)),
          axis.text.y = element_text(margin = margin(r = 5)))
  
  # Plot
  print(GeoSSE_model_div_pars_plot_i)
  
  dev.off()
  
  ## Print progress
  cat(paste0(Sys.time(), " - Parameter estimates plotted for ", model_i," - n° ", i, "/", length(model_list),"\n"))
  
}

### 6.3/ Aggregate parameter plots of all models in unique PDFs ####

## 6.3.1/ All parameters plots ####

# all_pars_plots_path <- list.files(path = "./outputs/GeoSSE/", pattern = "Model_all_pars_plot_", full.names = T)
# nb_ggplots <- length(all_pars_plots_path) ; nb_ggplots

# all_pars_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_all_pars_plot_", model_list, ".pdf")
all_pars_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_all_pars_plot_", model_list, ".pdf")

# qpdf::pdf_combine(input = all_pars_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_all_pars_plot_all_models.pdf"))
qpdf::pdf_combine(input = all_pars_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_all_pars_plot_all_models.pdf"))

## 6.3.2/ Diversification parameters plots ####

# all_div_pars_plots_path <- list.files(path = "./outputs/GeoSSE/", pattern = "Model_div_pars_plot_", full.names = T)
# nb_ggplots <- length(all_div_pars_plots_path) ; nb_ggplots

# all_div_pars_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_div_pars_plot_", model_list, ".pdf")
all_div_pars_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_div_pars_plot_", model_list, ".pdf")

# qpdf::pdf_combine(input = all_div_pars_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/Model_div_pars_plot_all_models.pdf"))
qpdf::pdf_combine(input = all_div_pars_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/Model_div_pars_plot_all_models.pdf"))


##### 7/ Plot relational igraphs of model parameters across states #####

# Schematic summary of model parameter estimates: lambda - mu per state, vicariance, range extension, range extirpation
# Use residence times as node size
# Plot Igraph networks for each model + AICc/Akaike's weight

# Load residence times summary table
# residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/residence_times_all_models_df.rds")
residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/residence_times_all_models_df.rds")

# Load GeoSSE models parameter summary table
# GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")
GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")

# Load GeoSSE model comparison summary table
# GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.rds")
GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.rds")


## Set color scheme for nodes and edges

# Edge = rates
# Nodes = residence times
# Color = Ranges x Regimes
  # Purple shades for Old World
  # Beige shades for New World
  # Blue shades for Both
  # Lighter = A
  # Darker = F

colors_for_OW_NW <- c("peachpuff2", "mediumpurple2")

# Set colors for New World states
col_NW_fn <- colorRampPalette(colors = c("plum", "mediumpurple2", "purple3"))
col_NW_scale <- col_NW_fn(n = 6)
names(col_NW_scale) <- paste0("NW-", LETTERS[1:6])

# Set colors for Old World states
col_OW_fn <- colorRampPalette(colors = c("papayawhip", "peachpuff2", "sandybrown"))
col_OW_scale <- col_OW_fn(n = 6)
names(col_OW_scale) <- paste0("OW-", LETTERS[1:6])

# Set colors for Wide-spread states
col_WS_fn <- colorRampPalette(colors = c("lightskyblue", "royalblue", "royalblue4"))
col_WS_scale <- col_WS_fn(n = 6)
names(col_WS_scale) <- paste0("WS-", LETTERS[1:6])

color_scheme_states <- c(col_NW_scale, col_OW_scale, col_WS_scale)

# Function to adjust node size of nodes based on the range scale used for overall data
rescale_node_size <- function(x, data_min = 0, data_max, range_min = 5, range_max = 30) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}

# Function to adjust edge size of nodes based on the range scale used for overall data
rescale_edge_size <- function(x, data_min = 0, data_max, range_min = 1, range_max = 6) 
{
  y <- x - data_min
  y_0_1 <- y / data_max
  y_0_max <- y_0_1 * (range_max - range_min)
  y_min_max <- y_0_max + range_min
  return(y_min_max)
}


# Keep only non-aggregated states
residence_times_all_models_df_for_igraph <- residence_times_all_models_df %>% 
  filter(!is.na(range)) %>%
  filter(!is.na(regime))

# Extract maximum residence time to scale node sizes
max_time <- max(residence_times_all_models_df_for_igraph$residence_time)

# Extract maximum log_rates to scale edge width
GeoSSE_pars_all_models_df_for_igraph <- GeoSSE_pars_all_models_df %>%
  filter(stats == "MLE_i") %>%
  filter(parameter_type %in% c("Speciation", "Extinction", "Vicariance", "Range-shift", "Range extension", "Range extirpation", "Regime shift")) %>%
  filter() %>%
  mutate(log_rates = log1p(value))
max_rate <- max(GeoSSE_pars_all_models_df_for_igraph$log_rates)


## One plot per model

## Loop per model
for (i in seq_along(model_list))
{
  # i <- 1
  
  model_i <- model_list[i]
  model_MLE_i <- model_fit_list[[i]]
 
  # Extract model evaluation
  GeoSSE_models_comparison_df_i <- GeoSSE_models_comparison_df %>% 
    filter(model == model_i)
  
  ### 7.1/ Prepare node data ####
  
  # Extract residence time data
  node_metadata_df_i <- residence_times_all_models_df_for_igraph %>% 
    filter(model == model_i)
  
  # Compute node size ranging from 5 to 30
  node_metadata_df_i$node_size <- rescale_node_size(x = node_metadata_df_i$residence_time, data_max = max_time)
  
  # Create node labels with % of residence time
  node_metadata_df_i$node_label <- paste0(node_metadata_df_i$state, "\n",
                                          round(node_metadata_df_i$residence_time_perc, 1), " %")
  
  # Define manually the node positions
  # node_layout_df <- data.frame(state = names(color_scheme_states),
  #                              x = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75),
  #                              y = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 1, 1, 1, 1, 1, 1))
  
  node_layout_df <- data.frame(state = names(color_scheme_states),
                               x = c(0, 0.2, 0, 0.2, 0, 0.2, 0.8, 1, 0.8, 1, 0.8, 1, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75),
                               y = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 1, 0.95, 1, 0.95, 1, 0.95))
  
  node_metadata_df_i <- left_join(node_metadata_df_i, node_layout_df)
  
  node_metadata_df_i
  
  
  ### 7.2/ Prepare edge data ####
  
  # Extract rates data
  edge_metadata_df_i <- GeoSSE_pars_all_models_df_for_igraph %>% 
    filter(model == model_i) %>% 
    rename(range = state)
  # Add regime "A" as default
  edge_metadata_df_i$regime[(edge_metadata_df_i$regime == "") & (edge_metadata_df_i$parameter_type != c("Regime shift"))] <- "A"
  # Update ranges and states
  # edge_metadata_df_i$range[edge_metadata_df_i$range == "Both"] <- "WS"
  edge_metadata_df_i$state <- paste0(edge_metadata_df_i$range, "-", edge_metadata_df_i$regime)
  
  # Compute net div. rates
  net_div_rates_df_i <- edge_metadata_df_i %>% 
    select(-pars, -log_rates) %>%
    filter(parameter_type %in% c("Speciation", "Extinction")) %>%
    pivot_wider(names_from = parameter_type, values_from = value) %>%
    mutate(Extinction = replace_na(Extinction, replace = 0)) %>%
    mutate(net_div = Speciation - Extinction) %>%
    pivot_longer(cols = c(net_div), names_to = "parameter_type", values_to = "value") %>%
    select(-Speciation, -Extinction)
  
  # Aggregate net div. rates with transition rates
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    select(-pars, -log_rates) %>%
    rbind(net_div_rates_df_i) %>%
    filter(!parameter_type %in% c("Speciation", "Extinction"))
  
  # Compute log_rates as log(1 + rate)
  edge_metadata_df_i$log_rates = log1p(edge_metadata_df_i$value)
  # Rescale to edge_width
  edge_metadata_df_i$edge_width <- rescale_edge_size(x = edge_metadata_df_i$log_rates, data_max = max_rate)
  
  ## Assign source/dest states for diversification rates
  edge_metadata_df_i$source[edge_metadata_df_i$parameter_type == "net_div"] <- edge_metadata_df_i$state[edge_metadata_df_i$parameter_type == "net_div"]
  edge_metadata_df_i$dest[edge_metadata_df_i$parameter_type == "net_div"] <- edge_metadata_df_i$state[edge_metadata_df_i$parameter_type == "net_div"]
  
  ## Assign source/dest states for range-shift rates
  range_shift_df_i <- edge_metadata_df_i %>% 
    filter(parameter_type == "Range-shift")
  range_shift_sources_dests <- range_shift_df_i$range
  # Extract source from range label
  range_shift_df_i$source <- str_remove(string = range_shift_sources_dests, pattern = " to .*") 
  # Extract dest from range label
  range_shift_df_i$dest <- str_remove(string = range_shift_sources_dests, pattern = ".* to ") 
  # Extract all regimes
  all_regimes <- unique(node_metadata_df_i$regime) 
  # Replicate edge data across all regimes
  range_shift_sources_across_regimes <- as.vector(outer(X = range_shift_df_i$source, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_shift_dests_across_regimes <- as.vector(outer(X = range_shift_df_i$dest, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_shift_df_across_regimes_i <- data.frame(source = range_shift_df_i$source, dest = range_shift_df_i$dest,
                                                source_across_regimes = range_shift_sources_across_regimes, dest_across_regimes = range_shift_dests_across_regimes,
                                                regime = rep(all_regimes, each = length(range_shift_sources_dests)))
  # Merge with initial range shift data
  range_shift_df_i <- range_shift_df_i %>% 
    select(-regime)
  range_shift_df_i <- left_join(range_shift_df_across_regimes_i, range_shift_df_i) %>% 
    select(-source, -dest) %>%
    rename(source = source_across_regimes,
           dest = dest_across_regimes) %>%
    select(model, stats, value, parameter_type, range, regime, state, log_rates, edge_width, source, dest)
  # Merge with initial edge data
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    filter(!parameter_type == "Range-shift") %>%
    rbind(range_shift_df_i)
  
  ## Assign source/dest states for range extension rates
  range_extension_df_i <- edge_metadata_df_i %>% 
    filter(parameter_type == "Range extension")
  range_extension_sources_dests <- range_extension_df_i$range
  # Extract source from range label
  range_extension_df_i$source <- str_remove(string = range_extension_sources_dests, pattern = "To ")
  range_extension_df_i$source <- c("OW", "NW")[-1 * (match(x = range_extension_df_i$source, table = c("OW", "NW")) - 1) + 2]
  # Assign destination as the Wide-spread range (WS)
  range_extension_df_i$dest <- "WS"
  # Extract all regimes
  all_regimes <- unique(node_metadata_df_i$regime)
  # Replicate edge data across all regimes
  range_extension_sources_across_regimes <- as.vector(outer(X = range_extension_df_i$source, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_extension_dests_across_regimes <- as.vector(outer(X = range_extension_df_i$dest, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_extension_df_across_regimes_i <- data.frame(source = range_extension_df_i$source, dest = range_extension_df_i$dest,
                                                    source_across_regimes = range_extension_sources_across_regimes, dest_across_regimes = range_extension_dests_across_regimes,
                                                    regime = rep(all_regimes, each = length(range_extension_sources_dests)))
  # Merge with initial range extension data
  range_extension_df_i <- range_extension_df_i %>% 
    select(-regime)
  range_extension_df_i <- left_join(range_extension_df_across_regimes_i, range_extension_df_i) %>% 
    select(-source, -dest) %>%
    rename(source = source_across_regimes,
           dest = dest_across_regimes) %>%
    select(model, stats, value, parameter_type, range, regime, state, log_rates, edge_width, source, dest)
  # Merge with initial edge data
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    filter(!parameter_type == "Range extension") %>%
    rbind(range_extension_df_i)
  
  ## Assign source/dest states for range extirpation rates
  range_extirpation_df_i <- edge_metadata_df_i %>% 
    filter(parameter_type == "Range extirpation")
  range_extirpation_sources_dests <- range_extirpation_df_i$range
  # Assign source as the Wide-spread range (WS)
  range_extirpation_df_i$source <- "WS"
  # Extract dest from range label
  range_extirpation_df_i$dest <- str_remove(string = range_extirpation_sources_dests, pattern = "From ")
  range_extirpation_df_i$dest <- c("OW", "NW")[-1 * (match(x = range_extirpation_df_i$dest, table = c("OW", "NW")) - 1) + 2]
  # Extract all regimes
  all_regimes <- unique(node_metadata_df_i$regime)
  # Replicate edge data across all regimes
  range_extirpation_sources_across_regimes <- as.vector(outer(X = range_extirpation_df_i$source, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_extirpation_dests_across_regimes <- as.vector(outer(X = range_extirpation_df_i$dest, Y = all_regimes, FUN = function (x, y) { paste0(x, "-", y) } ))
  range_extirpation_df_across_regimes_i <- data.frame(source = range_extirpation_df_i$source, dest = range_extirpation_df_i$dest,
                                                      source_across_regimes = range_extirpation_sources_across_regimes, dest_across_regimes = range_extirpation_dests_across_regimes,
                                                      regime = rep(all_regimes, each = length(range_extirpation_sources_dests)))
  # Merge with initial range extirpation data
  range_extirpation_df_i <- range_extirpation_df_i %>% 
    select(-regime)
  range_extirpation_df_i <- left_join(range_extirpation_df_across_regimes_i, range_extirpation_df_i) %>% 
    select(-source, -dest) %>%
    rename(source = source_across_regimes,
           dest = dest_across_regimes) %>%
    select(model, stats, value, parameter_type, range, regime, state, log_rates, edge_width, source, dest)
  # Merge with initial edge data
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    filter(!parameter_type == "Range extirpation") %>%
    rbind(range_extirpation_df_i)
  
  ## Assign source/dest states for regime shift rates
  if (sum(edge_metadata_df_i$parameter_type == "Regime shift") > 0) # Only if their are multiple regimes!
  {
    regime_shift_df_i <- edge_metadata_df_i %>% 
      filter(parameter_type == "Regime shift") %>%
      mutate(key = 1)
    # Find all combination of regime shifts among states
    all_pairwise_transitions <- combn(x = node_metadata_df_i$state, m = 2)
    all_regime_pairs <- combn(x = node_metadata_df_i$regime, m = 2)
    all_range_pairs <- combn(x = node_metadata_df_i$range, m = 2)
    different_regimes_indices <- all_regime_pairs[1, ] != all_regime_pairs[2, ]
    same_ranges_indices <- all_range_pairs[1, ] == all_range_pairs[2, ]
    all_pairwise_transitions <- all_pairwise_transitions[, different_regimes_indices & same_ranges_indices ]
    # Create regime shift df with all valid transitions
    regime_shift_df_across_regimes_i <- data.frame(source_across_regimes = c(all_pairwise_transitions[1, ], all_pairwise_transitions[2, ]),
                                                   dest_across_regimes = c(all_pairwise_transitions[2, ], rev(all_pairwise_transitions[1, ])),
                                                   key = 1)
    # Merge with initial range extirpation data
    regime_shift_df_i <- left_join(regime_shift_df_across_regimes_i, regime_shift_df_i) %>% 
      select(-source, -dest) %>%
      rename(source = source_across_regimes,
             dest = dest_across_regimes) %>%
      select(model, stats, value, parameter_type, range, regime, state, log_rates, edge_width, source, dest)
    # Merge with initial edge data
    edge_metadata_df_i <- edge_metadata_df_i %>% 
      filter(!parameter_type == "Regime shift") %>%
      rbind(regime_shift_df_i)
  }
  
  ## Adjust net_div in Wide-spread state (WS) to be summed with range extirpation and remove net_div WS edges
  # Extract vicariance data
  vicariance_df_i <- edge_metadata_df_i %>% 
    # filter(parameter_type == "net_div" & range == "WS") %>%
    filter(parameter_type == "Vicariance" & range == "WS") %>%
    rename(value_to_add = value) %>%
    select(value_to_add, regime) 
  # Extract range extirpation data
  range_extirpation_df_i <- edge_metadata_df_i %>% 
    filter(parameter_type == "Range extirpation")
  # Merge by summing values
  range_extirpation_df_i <- range_extirpation_df_i %>% 
    left_join(vicariance_df_i) %>%
    mutate(value = value + value_to_add) %>%
    select(-value_to_add)
  # Update log_rates and edge_width
  range_extirpation_df_i$log_rates <- log1p(range_extirpation_df_i$value)
  range_extirpation_df_i$edge_width <- rescale_edge_size(x = range_extirpation_df_i$log_rates, data_max = max_rate)
  # Merge with initial edge data
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    # filter(!(parameter_type == "net_div" & range == "WS")) %>%
    filter(!(parameter_type == "Vicariance" & range == "WS")) %>%
    filter(!(parameter_type == "Range extirpation")) %>%
    rbind(range_extirpation_df_i)
  
  # Select and order variables
  edge_metadata_df_i <- edge_metadata_df_i %>% 
    rename(rate = value) %>%
    select(source, dest, parameter_type, rate, log_rates, edge_width, model) %>%
    mutate(source_range = str_remove(string = source, pattern = "-.*")) %>%
    mutate(source_regime = str_remove(string = source, pattern = ".*-")) %>%
    mutate(dest_range = str_remove(string = dest, pattern = "-.*")) %>%
    mutate(dest_regime = str_remove(string = dest, pattern = ".*-"))
  
  # Create edge labels
  edge_metadata_df_i$edge_label <- format(round(edge_metadata_df_i$rate, 3), nsmall = 3)
  
  edge_metadata_df_i
  
  # Merge and save metadata used to generate igraph
  igraph_metadata_df_i <- list(node_metadata_df = node_metadata_df_i,
                               edge_metadata_df = edge_metadata_df_i)
  # saveRDS(object = igraph_metadata_df_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  saveRDS(object = igraph_metadata_df_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  
  
  ### 7.3/ Plot relational igraphs of model parameters ####
  
  # Load metadata used to generate igraph
  # igraph_metadata_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  igraph_metadata_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  
  ## 7.3.1/ Convert to igraph object ####
  
  ?igraph::graph_from_data_frame
  
  ## Convert to igraph
  SSE_model_igraph_i <- igraph::graph_from_data_frame(d = igraph_metadata_df_i$edge_metadata_df, vertices = igraph_metadata_df_i$node_metadata_df, directed = T)
  # str(SSE_model_igraph_i)
  
  # Explore igraph
  # SSE_model_igraph_i
  E(SSE_model_igraph_i) # Edge attributes
  V(SSE_model_igraph_i) # Vertex/Node attributes
  SSE_model_igraph_i[]  # Adjacency matrix of nodes x nodes
  
  # Export adjacency matrix
  as_adjacency_matrix(SSE_model_igraph_i, attr = "rate")
  as_adjacency_matrix(SSE_model_igraph_i, attr = "log_rates")
  as_adjacency_matrix(SSE_model_igraph_i, attr = "edge_width")
  
  ## 7.3.2/ Adjust aesthetics ####
  
  # # Remove loops
  # SSE_model_igraph_i <- simplify(graph = SSE_model_igraph_i,
  #                                remove.multiple = F, remove.loops = F) 
  
  # Set node colors based on states:
  V(SSE_model_igraph_i)$color <- color_scheme_states[names(V(SSE_model_igraph_i))]
  # Set node size based on scaled residence times
  V(SSE_model_igraph_i)$size <- V(SSE_model_igraph_i)$node_size * 3
  # V(SSE_model_igraph_i)$size <- sapply(X = V(SSE_model_igraph_i)$size, FUN = function (x) { max(x, 10) })
  # Set node font size based on scaled residence times
  V(SSE_model_igraph_i)$label.cex <- V(SSE_model_igraph_i)$node_size / 20
  # V(SSE_model_igraph_i)$label.cex <- sapply(X = V(SSE_model_igraph_i)$label.cex, FUN = function (x) { max(x, 0.5) })
    
  # Set edge width based on scaled log rates
  E(SSE_model_igraph_i)$width <- E(SSE_model_igraph_i)$edge_width * 3
  # Set edge arrow size based on event counts
  E(SSE_model_igraph_i)$arrow.size <- E(SSE_model_igraph_i)$edge_width / 4
  
  # Set edge label as rates
  # E(SSE_model_igraph_i)$edge.label <- E(SSE_model_igraph_i)$edge_label
  
  # Set edge color according to their source
  edge_starts <- ends(SSE_model_igraph_i, es = E(SSE_model_igraph_i), names = F)[, 1]
  edge_col <- V(SSE_model_igraph_i)$color[edge_starts]
  E(SSE_model_igraph_i)$color <- edge_col
  
  # # Remove edges with low rates
  # table(E(SSE_model_igraph_i)$rate)
  # cutoff_min_rates_i <- 0.001
  # SSE_model_igraph_i <- delete_edges(graph = SSE_model_igraph_i,
  #                                    edges = E(SSE_model_igraph_i)[rate < cutoff_min_counts_i])
  
  # Adjust layout to coordinates of vertex/nodes
  layout_manual <- cbind(V(SSE_model_igraph_i)$x, V(SSE_model_igraph_i)$y)
  
  ## Save igraph object
  # saveRDS(object = SSE_model_igraph_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_object_all_states_", model_i, ".rds"))
  saveRDS(object = SSE_model_igraph_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_object_all_states_", model_i, ".rds"))
  
  ## 7.3.3/ Plot igraph ####
  
  # pdf(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_all_states_",model_i,".pdf"), width = 6, height = 6)
  pdf(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_all_states_",model_i,".pdf"), width = 6, height = 6)
  
  par(mar = c(0.5, 2.0, 7.0, 2.0)) # bltr
  plot.igraph(x = SSE_model_igraph_i,
              ## Node aesthetics
              # vertex.color = colors_list_for_areas, # Nodes color
              vertex.frame.color = "black", # Nodes border color
              vertex.shape = "circle",  # Node shape
              # vertex.size = 15, # Node size
              vertex.label	= igraph_metadata_df_i$node_metadata_df$node_label, # Node labels
              vertex.label.color = "black",  # Node label color
              vertex.label.font = 2, # Node label font
              # vertex.label.cex = 1, # Node label cex
              ## Edge aesthetics
              # edge.color = "black",     # Edge color
              # edge.width = 1,           # Edge width
              # edge.arrow.size = 0.7,    # Terminal arrow size
              edge.label = E(SSE_model_igraph_i)$edge_label,
              edge.label.cex = 1.0,
              edge.lty = 1,             # Edge line type
              edge.curved = 0.1,
              loop.size = 2.0,
              arrow.mode = "forward",    # Direction of arrows
              ## Other aesthetics
              layout = layout_manual,
              # layout = layout_nicely, # Default layout
              margin = c(0,0,-0.4,0.5), # To adjust plot margins # bltr
              main = paste0(model_i," model\n\n",
                            "AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
                            "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights, 1), " %")
  )
  
  dev.off()
  
  ## Print progress
  cat(paste0(Sys.time(), " - Igraph relational scheme plotted for ", model_i," - n° ", i, "/", length(model_list),"\n"))
  
}


### 7.4/ Merge Igraph plots of all models in a single PDF ####

# all_igraph_network_all_states_plots_path <- list.files(path = "./outputs/GeoSSE/", pattern = "igraph_network_all_states_", full.names = T)
# nb_ggplots <- length(all_igraph_network_all_states_plots_path) ; nb_ggplots

# all_igraph_network_all_states_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_all_states_", model_list, ".pdf")
all_igraph_network_all_states_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_all_states_", model_list, ".pdf")

# qpdf::pdf_combine(input = all_igraph_network_all_states_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_all_states_all_models.pdf"))
qpdf::pdf_combine(input = all_igraph_network_all_states_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_all_states_all_models.pdf"))


##### 8/ Plot relational igraphs of model parameters across ranges #####

# Aggregate regimes (hidden states) within each range (observed state) as a pie proportional to residence times
# Edge rates should be aggregated across regimes using a weighted mean of summed regime residence times
 # Ex: rate from 0 to 1 = weighted mean of 0A to 1A and 0B to 1B with weights being residence times in 0A+1A and 0B+1B.

# Load residence times summary table
# residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/residence_times_all_models_df.rds")
residence_times_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/residence_times_all_models_df.rds")

# Load GeoSSE models parameter summary table
# GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")
GeoSSE_pars_all_models_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_pars_all_models_df.rds")

# Load GeoSSE model comparison summary table
# GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/GeoSSE_models_comparison_df.rds")
GeoSSE_models_comparison_df <- readRDS(file = "./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/GeoSSE_models_comparison_df.rds")

## Set color scheme for nodes and edges

color_scheme_ranges <- c("mediumpurple2", "peachpuff2", "lightskyblue")
names(color_scheme_ranges) <- c("NW", "OW", "WS")

## Loop per model
for (i in seq_along(model_list))
{
  # i <- 1
  
  model_i <- model_list[i]
  model_MLE_i <- model_fit_list[[i]]
  
  # Extract model evaluation
  GeoSSE_models_comparison_df_i <- GeoSSE_models_comparison_df %>% 
    filter(model == model_i)
  
  # Load metadata used to generate igraph with all states as ranges x regimes
  # igraph_metadata_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  igraph_metadata_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_metadata_df_all_states_", model_i, ".rds"))
  
  # Initiate metadata for aggregated states
  igraph_metadata_aggregated_states_df_i <- igraph_metadata_df_i
  
  ### 8.1/ Aggregate node metadata across regimes ####
  
  node_metadata_aggregated_states_df_i <- igraph_metadata_df_i$node_metadata_df %>% 
    group_by(range) %>%
    summarize(residence_time = sum(residence_time)) %>%
    mutate(residence_time_perc = residence_time / sum(residence_time) * 100) %>%
    mutate(model = model_i)
  
  # Compute node size ranging from 5 to 30
  node_metadata_aggregated_states_df_i$node_size <- rescale_node_size(x = node_metadata_aggregated_states_df_i$residence_time, data_max = max_time)
  
  # Create node labels with % of residence time
  node_metadata_aggregated_states_df_i$node_label <- paste0(node_metadata_aggregated_states_df_i$range, "\n",
                                                            round(node_metadata_aggregated_states_df_i$residence_time_perc, 1), " %")
  
  # Define manually the node positions
  
  node_layout_df <- data.frame(range = names(color_scheme_ranges),
                               x = c(0, 1, 0.5),
                               y = c(0, 0.5, 1))
  
  node_metadata_aggregated_states_df_i <- left_join(node_metadata_aggregated_states_df_i, node_layout_df)
  
  # Store updated node metadata
  igraph_metadata_aggregated_states_df_i$node_metadata_df <- node_metadata_aggregated_states_df_i 
  
  ### 8.2/ Aggregate edge metadata across regimes ####
  
  # Edge rates should be aggregated across regimes using a weighted mean of summed regime residence times
  # Ex: rate from 0 to 1 = weighted mean of 0A to 1A and 0B to 1B with weights being residence times in 0A+1A and 0B+1B.
  
  # Join residence times of source and dest states
  edge_metadata_aggregated_states_df_i <- igraph_metadata_df_i$edge_metadata_df %>%
    left_join(igraph_metadata_df_i$node_metadata_df[, c("state", "residence_time")], by = join_by(source == state)) %>%
    rename(source_residence_time = residence_time) %>%
    left_join(igraph_metadata_df_i$node_metadata_df[, c("state", "residence_time")], by = join_by(dest == state)) %>%
    rename(dest_residence_time = residence_time) %>%
    mutate(weights = source_residence_time + dest_residence_time)
  
  # Compute mean weighted rates per ranges using residence time as weights
  edge_metadata_aggregated_states_df_i <- edge_metadata_aggregated_states_df_i %>% 
    group_by(parameter_type, source_range, dest_range) %>%
    summarize(rate = weighted.mean(rate, w = weights)) %>%
    filter(parameter_type != "Regime shift")
  
  # Update source/dest and log_rates
  edge_metadata_aggregated_states_df_i <- edge_metadata_aggregated_states_df_i %>% 
    mutate(source = source_range,
           dest = dest_range,
           log_rates = log1p(rate),
           model = model_i)
  
  # Compute edge size ranging from 1 to 6
  edge_metadata_aggregated_states_df_i$edge_width <- rescale_edge_size(x = edge_metadata_aggregated_states_df_i$log_rates, data_max = max_rate)
  
  # Create edge labels
  edge_metadata_aggregated_states_df_i$edge_label <- format(round(edge_metadata_aggregated_states_df_i$rate, 3), nsmall = 3)
  
  # Reorder columns
  edge_metadata_aggregated_states_df_i <- edge_metadata_aggregated_states_df_i %>% 
    dplyr::select(source, dest, parameter_type, rate, log_rates, edge_width, model, source_range, dest_range, edge_label)
  
  # Store updated node metadata
  igraph_metadata_aggregated_states_df_i$edge_metadata_df <- edge_metadata_aggregated_states_df_i
  
  # Save metadata for igraph with aggregated states
  # saveRDS(object = igraph_metadata_aggregated_states_df_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_metadata_df_aggregated_states_", model_i, ".rds"))
  saveRDS(object = igraph_metadata_aggregated_states_df_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_metadata_df_aggregated_states_", model_i, ".rds"))
  
  ### 8.3/ Plot relational igraphs of model parameters across aggregated states ####
  
  # Load metadata used to generate igraph
  # igraph_metadata_aggregated_states_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_metadata_df_aggregated_states_", model_i, ".rds"))
  igraph_metadata_aggregated_states_df_i <- readRDS(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_metadata_df_aggregated_states_", model_i, ".rds"))
  
  ## 8.3.1/ Convert to igraph object ####
  
  ?igraph::graph_from_data_frame
  
  ## Convert to igraph
  SSE_model_igraph_i <- igraph::graph_from_data_frame(d = igraph_metadata_aggregated_states_df_i$edge_metadata_df, vertices = igraph_metadata_aggregated_states_df_i$node_metadata_df, directed = T)
  # str(SSE_model_igraph_i)
  
  # Explore igraph
  # SSE_model_igraph_i
  E(SSE_model_igraph_i) # Edge attributes
  V(SSE_model_igraph_i) # Vertex/Node attributes
  SSE_model_igraph_i[]  # Adjacency matrix of nodes x nodes
  
  # Export adjacency matrix
  as_adjacency_matrix(SSE_model_igraph_i, attr = "rate")
  as_adjacency_matrix(SSE_model_igraph_i, attr = "log_rates")
  as_adjacency_matrix(SSE_model_igraph_i, attr = "edge_width")
  
  ## 8.3.2/ Adjust aesthetics ####
  
  # # Remove loops
  # SSE_model_igraph_i <- simplify(graph = SSE_model_igraph_i,
  #                                remove.multiple = F, remove.loops = F) 
  
  # Set node colors based on states:
  V(SSE_model_igraph_i)$color <- color_scheme_ranges[names(V(SSE_model_igraph_i))]
  # Set node size based on scaled residence times
  V(SSE_model_igraph_i)$size <- V(SSE_model_igraph_i)$node_size * 3
  # V(SSE_model_igraph_i)$size <- sapply(X = V(SSE_model_igraph_i)$size, FUN = function (x) { max(x, 10) })
  # Set node font size based on scaled residence times
  V(SSE_model_igraph_i)$label.cex <- V(SSE_model_igraph_i)$node_size / 20
  V(SSE_model_igraph_i)$label.cex <- sapply(X = V(SSE_model_igraph_i)$label.cex, FUN = function (x) { max(x, 0.5) })
  
  # Extract pie fractions as residence times per range
  pie_fractions_list <- split(x = igraph_metadata_df_i$node_metadata_df$residence_time_perc_per_range, f = igraph_metadata_df_i$node_metadata_df$range)
  V(SSE_model_igraph_i)$pie <- pie_fractions_list
  # Extract pie colors from states color scheme
  states_colors <- color_scheme_states[igraph_metadata_df_i$node_metadata_df$state]
  pie_colors_list <- split(x = states_colors, f = igraph_metadata_df_i$node_metadata_df$range)
  V(SSE_model_igraph_i)$pie.color <- pie_colors_list
  
  # Set edge width based on scaled log rates
  E(SSE_model_igraph_i)$width <- E(SSE_model_igraph_i)$edge_width * 3
  # Set edge arrow size based on event counts
  E(SSE_model_igraph_i)$arrow.size <- E(SSE_model_igraph_i)$edge_width / 4
  
  # Set edge label as rates
  # E(SSE_model_igraph_i)$edge.label <- E(SSE_model_igraph_i)$edge_label
  
  # Set edge color according to their source
  edge_starts <- ends(SSE_model_igraph_i, es = E(SSE_model_igraph_i), names = F)[, 1]
  edge_col <- V(SSE_model_igraph_i)$color[edge_starts]
  E(SSE_model_igraph_i)$color <- edge_col
  
  # # Remove edges with low rates
  # table(E(SSE_model_igraph_i)$rate)
  # cutoff_min_rates_i <- 0.001
  # SSE_model_igraph_i <- delete_edges(graph = SSE_model_igraph_i,
  #                                    edges = E(SSE_model_igraph_i)[rate < cutoff_min_counts_i])
  
  # Adjust layout to coordinates of vertex/nodes
  layout_manual <- cbind(V(SSE_model_igraph_i)$x, V(SSE_model_igraph_i)$y)
  
  ## Save igraph object
  # saveRDS(object = SSE_model_igraph_i, file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_object_aggregated_states_", model_i, ".rds"))
  saveRDS(object = SSE_model_igraph_i, file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_object_aggregated_states_", model_i, ".rds"))
  
  
  ## 8.3.3/ Plot igraph ####
  
  # pdf(file = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_aggregated_states_",model_i,".pdf"), width = 6, height = 6)
  pdf(file = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_aggregated_states_",model_i,".pdf"), width = 6, height = 6)
  
  par(mar = c(0.5, 2.0, 7.0, 2.0)) # bltr
  plot.igraph(x = SSE_model_igraph_i,
              ## Node aesthetics
              # vertex.color = colors_list_for_areas, # Nodes color
              vertex.frame.color = "black", # Nodes border color
              # vertex.shape = "circle",  # Node shape
              vertex.shape = "pie", # Node shape as pies
              # vertex.size = 15, # Node size
              # vertex.pie = list(), # Set size of pie fractions
              # vertex.pie.color = list(), # Set colors of pie fractions
              vertex.label	= igraph_metadata_aggregated_states_df_i$node_metadata_df$node_label, # Node labels
              vertex.label.color = "black",  # Node label color
              vertex.label.font = 2, # Node label font
              # vertex.label.cex = 1, # Node label cex
              ## Edge aesthetics
              # edge.color = "black",     # Edge color
              # edge.width = 1,           # Edge width
              # edge.arrow.size = 0.7,    # Terminal arrow size
              edge.label = E(SSE_model_igraph_i)$edge_label,
              edge.label.cex = 1.0,
              edge.lty = 1,             # Edge line type
              edge.curved = 0.1,
              loop.size = 2.0,
              arrow.mode = "forward",    # Direction of arrows
              ## Other aesthetics
              layout = layout_manual,
              # layout = layout_nicely, # Default layout
              margin = c(0,0,-0.4,0.5), # To adjust plot margins # bltr
              main = paste0(model_i," aggregated model\n\n",
                            "AICc = ",round(GeoSSE_models_comparison_df_i$AICc, 1),"\n",
                            "Akaike's weight = ", round(GeoSSE_models_comparison_df_i$Akaike_weights, 1), " %")
  )
  
  dev.off()
  
  ## Print progress
  cat(paste0(Sys.time(), " - Igraph aggregated relational scheme plotted for ", model_i," - n° ", i, "/", length(model_list),"\n"))
}

### 8.4/ Merge Igraph plots of all models in a single PDF ####

# all_igraph_network_aggregated_states_plots_path <- list.files(path = "./outputs/GeoSSE/", pattern = "igraph_network_aggregated_states_", full.names = T)
# nb_ggplots <- length(all_igraph_network_aggregated_states_plots_path) ; nb_ggplots

# all_igraph_network_aggregated_states_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_aggregated_states_", model_list, ".pdf")
all_igraph_network_aggregated_states_plots_path <- paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_aggregated_states_", model_list, ".pdf")

# qpdf::pdf_combine(input = all_igraph_network_aggregated_states_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_rough_phylogeny_1534t/igraph_network_aggregated_states_all_models.pdf"))
qpdf::pdf_combine(input = all_igraph_network_aggregated_states_plots_path, output = paste0("./outputs/GeoSSE/Ponerinae_MCC_phylogeny_1534t/igraph_network_aggregated_states_all_models.pdf"))






