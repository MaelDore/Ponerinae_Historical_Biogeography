##### Script 15: BAMM diversification models  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Estimate diversification rates and rate shift number and location along branches
# Test for differences in diversification rates across Bioregions
  # STRAPP Test: traitDependentBAMM function
  # All bioregions
  # Old World vs. New world
# Map the net dispersal rates based on tip rates in space
# Investigate a relationship between regime shift and dispersal events
  # Plot phylogeny scaled by posterior probability of rate shift to observe relationship between regime shift and bioregion membership/dispersal event
  # Test relationship between PP of rate shifts vs. PP of dispersal events
# Plot and test for latitudinal and longitudinal gradients of net diversification rates


###

### Inputs

# (Set of) Time-calibrated phylogeny(ies)
# Probability of bioregion membership per branch
# Stacked-SDM of taxa ranges

###

### Outputs

# BAMM models
# Rates of diversification along branches
# Localisation and probabilities of rate shift
# Shift-scaled phylogeny to observe relationship between regime shift and dispersal events
# LTT per bioregions (compare BAMM rates with empirical rates)
# STRAPP test for differences in diversification rates across Bioregions
# Test relationship between PP of rate shifts vs. PP of dispersal events
# Map of net dispersal rates based on tip rates
# Plot and test for latitudinal and longitudinal gradients of net diversification rates
   # Add null model values from block-permuted data


###

### Bonus in other scripts (?)

# For a more smooth model with many small cladogenetic changes, see CLaDS
# Compare in situ radiations
# GeoSSE, HiSSE, BiSSE for Old World vs. New world


# Clean environment
rm(list = ls())

library(ape)
library(phytools)
library(BAMMtools)
library(coda)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggtree)
library(parallel)
library(raster)
library(sf)
library(vcd)  # For association plots


##### 1/ Load input files ####

### 1.1/ Inform path to BAMM folder ####

BAMM_path <- "./software/bamm-2.5.0/"

### 1.2/ Load phylogeny and check if it is valid

# Ponerinae_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")
Ponerinae_MCC_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")

# Check if ultrametic
# is.ultrametric(Ponerinae_phylogeny_1534t)
is.ultrametric(Ponerinae_MCC_phylogeny_1534t)
# Check if fully resolved
# is.binary(Ponerinae_phylogeny_1534t)
is.binary(Ponerinae_MCC_phylogeny_1534t)
# Check that all branch have positive length
# min(Ponerinae_phylogeny_1534t$edge.length) > 0 ; min(Ponerinae_phylogeny_1534t$edge.length)
min(Ponerinae_MCC_phylogeny_1534t$edge.length) > 0 ; min(Ponerinae_MCC_phylogeny_1534t$edge.length)

### 1.3/ Export phylogeny as .tree file for the analyses
# write.tree(phy = Ponerinae_phylogeny_1534t, file = paste0(BAMM_path, "Ponerinae_phylogeny.tree"))
write.tree(phy = Ponerinae_MCC_phylogeny_1534t, file = paste0(BAMM_path, "Ponerinae_MCC_phylogeny.tree"))

# phy_path <- paste0(BAMM_path, "Ponerinae_phylogeny.tree")
phy_path <- paste0(BAMM_path, "Ponerinae_MCC_phylogeny.tree")

##### 2/ Configure the control file for the run #####

### 2.1/ Load control file template for diversification analyses ####

div_config_file_template <- readLines(con = paste0(BAMM_path, "template_diversification.txt"))
my_div_control_file <- div_config_file_template

# For phenotypic trait evolution analyses
# trait_config_file_template <- readLines(con = paste0(BAMM_path, "template_trait.txt"))

### 2.2/ Set general settings and data input ####

# Set path to the phylogenetic tree file
phy_path_line <- which(str_detect(string = my_div_control_file, pattern = "treefile ="))
my_div_control_file[phy_path_line] <- paste0("treefile = ", phy_path)

# Set path to the file storing all information on this run
run_info_path <- "run_info.txt"
run_info_path_line <- which(str_detect(string = my_div_control_file, pattern = "runInfoFilename = "))
my_div_control_file[run_info_path_line] <- paste0("runInfoFilename = ", run_info_path)

# Should the run sample from prior only?
  # To check that the prior settings are fine
  # To compare with the posterior, to check the influence of the prior on the posterior
sampleFromPriorOnly <- 0 # 0 = No. 1 = Yes
sampleFromPriorOnly_line <- which(str_detect(string = my_div_control_file, pattern = "sampleFromPriorOnly = "))
my_div_control_file[sampleFromPriorOnly_line] <- paste0("sampleFromPriorOnly = ", sampleFromPriorOnly)

# Should the MCMC simulation be performed?
  # If runMCMC = 0, the program will only check whether the data file can be read and the initial likelihood computed
runMCMC <- 1  # 0 = No. 1 = Yes                            
runMCMC_line <- which(str_detect(string = my_div_control_file, pattern = "runMCMC = "))[1]
my_div_control_file[runMCMC_line] <- paste0("runMCMC = ", runMCMC)

# Should the prior distribution of the number of shift events, given the hyperprior on the Poisson rate parameter, be simulated?
  # This was necessary to compute Bayes factors
  # But is now disabled as the exact (analytical) prior is now implemented in BAMMtools.
simulatePriorShifts <- 0  # 0 = No. 1 = Yes 
simulatePriorShifts_line <- which(str_detect(string = my_div_control_file, pattern = "simulatePriorShifts = "))[1]
my_div_control_file[simulatePriorShifts_line] <- paste0("simulatePriorShifts = ", simulatePriorShifts)

# Whether to load a previous event data file 
loadEventData <- 0  # 0 = No. 1 = Yes                        
loadEventData_line <- which(str_detect(string = my_div_control_file, pattern = "loadEventData = "))[1]
my_div_control_file[loadEventData_line] <- paste0("loadEventData = ", loadEventData)

# Provides file name of the event data file to load, used only if loadEventData = 1
eventDataInfile_path <- "event_data_in.txt"
eventDataInfile_path_line <- which(str_detect(string = my_div_control_file, pattern = "eventDataInfile = "))[1]
my_div_control_file[eventDataInfile_path_line] <- paste0("eventDataInfile = ", eventDataInfile_path)

# Whether to initialize (but not run) the MCMC. 
  # If initializeModel = 0, the program will only ensure that the data files (e.g., treefile) can be read
initializeModel <- 1  # 0 = No. 1 = Yes 
initializeModel_line <- which(str_detect(string = my_div_control_file, pattern = "initializeModel = "))[1]
my_div_control_file[initializeModel_line] <- paste0("initializeModel = ", initializeModel)

# Whether to use a "global" sampling probability to assign the proportion of terminal represented in the tree.
  # If False (0), expects a file path for species-specific sampling probabilities (see sampleProbsFilename)
useGlobalSamplingProbability = 1  # 0 = No. 1 = Yes         
useGlobalSamplingProbability_line <- which(str_detect(string = my_div_control_file, pattern = "useGlobalSamplingProbability = "))[1]
my_div_control_file[useGlobalSamplingProbability_line] <- paste0("useGlobalSamplingProbability = ", useGlobalSamplingProbability)

# Provides the global sampling fraction
  # If useGlobalSamplingFraction = 0, this is ignored and BAMM looks for a file path to species-specific sampling fractions
globalSamplingFraction <- 1.0            
globalSamplingFraction_line <- which(str_detect(string = my_div_control_file, pattern = "globalSamplingFraction = "))[1]
my_div_control_file[globalSamplingFraction_line] <- paste0("globalSamplingFraction = ", globalSamplingFraction)

# Provides path to the file containing species-specific sampling fractions
sampleProbsFilename_path <- "taxa_sampling_probs.txt"
sampleProbsFilename_line <- which(str_detect(string = my_div_control_file, pattern = "sampleProbsFilename = "))[1]
my_div_control_file[sampleProbsFilename_line] <- paste0("sampleProbsFilename = ", sampleProbsFilename_path)

# Set the seed for the random number generator.
  # If not specified (or is -1), a seed is obtained from the system clock
seed <- 12345
seed_line <- which(str_detect(string = my_div_control_file, pattern = "seed = "))[1]
my_div_control_file[seed_line] <- paste0("seed = ", seed)

# Should the output files be overwritten?
  # If True (1), the program will overwrite any output files in the current directory (if present)
# overwrite <- 0
overwrite <- 1
overwrite_line <- which(str_detect(string = my_div_control_file, pattern = "overwrite = "))[1]
my_div_control_file[overwrite_line] <- paste0("overwrite = ", overwrite)

# Set limits to valid configurations
  # If 1, rejects proposals that cause a branch and both of its direct descendants to have at least one event. 
  # Such an event configuration may cause the parameters of the parent event to change to unrealistic values.
  # If 0, no such proposals are immediately rejected. The default value is 0.
# validateEventConfiguration <- 0
validateEventConfiguration <- 1
validateEventConfiguration_line <- which(str_detect(string = my_div_control_file, pattern = "validateEventConfiguration = "))[1]
my_div_control_file[validateEventConfiguration_line] <- paste0("validateEventConfiguration = ", validateEventConfiguration)


### 2.3/ Set (hyper)prior settings ####

## Can use this help function to automatically tune prior adapted to your data by scaling the prior distributions based on the age (root depth) of your tree
# In practice, setBAMMpriors first estimates the rate of speciation for your full tree under a pure birth model of diversification. 
# Then assume, arbitrarily, that a reasonable prior distribution for the initial lambda0/mu0 rate parameters is an exponential distribution with a mean five times greater than this pure birth value.
# Rationale = having a weakly informative prior that is still in the order of magnitude of the true rate
# For the shift parameter (alpha), the sd of the normal prior is set such as mean +/- 2s gives an alpha parameter that results in
# either a 90% decline in the evolutionary rate or a 190% increase in rate on the interval of time from the root to the tips of the tree.
setBAMMpriors(read.tree(paste0(BAMM_path, "my_phy.tree")))
default_tuned_priors <- readLines(con = "./myPriors.txt")
file.remove("./myPriors.txt")

# Set the expected number of shifts used to set the exponential hyperprior for nb of rate shifts (from which the Λ is drawn)
# Suggested values: 
  #  expectedNumberOfShifts = 1.0 for small trees (< 500 tips)
  #	 expectedNumberOfShifts = 10 or even 50 for large trees (> 5000 tips) 
# Good practice = set several runs with a range of values to be sure it does not affect results
expectedNumberOfShifts <- 1.0
expectedNumberOfShifts_default_line <- which(str_detect(string = default_tuned_priors, pattern = "expectedNumberOfShifts = "))[1]
expectedNumberOfShifts <- as.numeric(str_remove(string = default_tuned_priors[expectedNumberOfShifts_default_line], pattern = "expectedNumberOfShifts = "))
expectedNumberOfShifts_line <- which(str_detect(string = my_div_control_file, pattern = "expectedNumberOfShifts = "))[1]
my_div_control_file[expectedNumberOfShifts_line] <- paste0("expectedNumberOfShifts = ", expectedNumberOfShifts)

# Set the rate parameter of the exponential prior(s) of initial lambda parameters (lambda0) of speciation rate regimes
  # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
lambdaInitPrior <- 1.0
lambdaInitPrior <- 6.256328
lambdaInitPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaInitPrior = "))[1]
lambdaInitPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaInitPrior_default_line], pattern = "lambdaInitPrior = "))
lambdaInitPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaInitPrior = "))[1]
my_div_control_file[lambdaInitPrior_line] <- paste0("lambdaInitPrior = ", lambdaInitPrior)

# Set the standard deviation of the normal distribution prior(s) of rate variation parameters (alpha) of speciation rate regimes
  # alpha in lambda(t) = lamba0 x exp(alpha*t)
  # Mean of this prior(s) are fixed to zero such as a constant rate diversification process is the most probable a priori
lambdaShiftPrior <- 0.05
lambdaShiftPrior <- 0.01190669
lambdaShiftPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaShiftPrior = "))[1]
lambdaShiftPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaShiftPrior_default_line], pattern = "lambdaShiftPrior = "))
lambdaShiftPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaShiftPrior = "))[1]
my_div_control_file[lambdaShiftPrior_line] <- paste0("lambdaShiftPrior = ", lambdaShiftPrior)

# Set the rate parameter of the exponential prior(s) of initial lambda parameters (mu0) of extinction rate regimes
  # mu0 in mu(t) = mu0 x exp(alpha*t)
  # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
muInitPrior <- 1.0
muInitPrior <- 6.256328
muInitPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "muInitPrior = "))[1]
muInitPrior <- as.numeric(str_remove(string = default_tuned_priors[muInitPrior_default_line], pattern = "muInitPrior = "))
muInitPrior_line <- which(str_detect(string = my_div_control_file, pattern = "muInitPrior = "))[1]
my_div_control_file[muInitPrior_line] <- paste0("muInitPrior = ", muInitPrior)

# No prior for rate variation parameters (alpha) of extinction rate regimes as they are assumed to follow constant rates

# Set the prior (probability) of the time mode (of speciation?) being time-variable (vs. time-constant)
  # By default, allows all regimes to be time-variable, as their rate can still be estimated as constant with alpha = 0
lambdaIsTimeVariablePrior <- 1
lambdaIsTimeVariablePrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaIsTimeVariablePrior = "))[1]
my_div_control_file[lambdaIsTimeVariablePrior_line] <- paste0("lambdaIsTimeVariablePrior = ", lambdaIsTimeVariablePrior)


### 2.4/ Set the MCMC simulation settings, MCMC logs and output options ####

# Set the number of generations to perform MCMC simulation
numberOfGenerations = format(10000000, scientific = F) # 10^7
# numberOfGenerations = format(100000, scientific = F) # For the test run
# numberOfGenerations = 1000 # 10^3
numberOfGenerations_line <- which(str_detect(string = my_div_control_file, pattern = "numberOfGenerations = "))[1]
my_div_control_file[numberOfGenerations_line] <- paste0("numberOfGenerations = ", numberOfGenerations)

# Set the path to the MCMC output file
  # Includes only summary information about MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
mcmcOutfile_path <- "mcmc_log.txt"
mcmcOutfile_path_line <- which(str_detect(string = my_div_control_file, pattern = "mcmcOutfile = "))[1]
my_div_control_file[mcmcOutfile_path_line] <- paste0("mcmcOutfile = ", mcmcOutfile_path)

# Set the frequency in which to write the MCMC output to the log file
  # Aim for 500-5000 posterior samples ideally
  # Will need to remove some to account for the burn-in
mcmcWriteFreq <- round(as.numeric(numberOfGenerations) / 2000)
# mcmcWriteFreq <- 10 
mcmcWriteFreq_line <- which(str_detect(string = my_div_control_file, pattern = "mcmcWriteFreq = "))[1]
my_div_control_file[mcmcWriteFreq_line] <- paste0("mcmcWriteFreq = ", mcmcWriteFreq)

# Set the path to the main output file
  # The raw event data. ALL of the results are contained in this file,
  # All branch-specific speciation rates, shift positions, marginal distributions, etc. can be reconstructed from this output.
  # See R package BAMMtools for working with this output
eventDataOutfile_path <- "event_data.txt"
eventDataOutfile_path_line <- which(str_detect(string = my_div_control_file, pattern = "eventDataOutfile = "))[1]
my_div_control_file[eventDataOutfile_path_line] <- paste0("eventDataOutfile = ", eventDataOutfile_path)

# Set frequency in which to write the event data to the output file = the sampling frequency of posterior samples
  # Aim for 500-5000 posterior samples ideally
  # Will need to remove some to account for the burn-in
eventDataWriteFreq <- as.numeric(numberOfGenerations) / 2000 # Sample every 5*10^3 generations
# eventDataWriteFreq <- 100 # Sample every 100 generations
eventDataWriteFreq_line <- which(str_detect(string = my_div_control_file, pattern = "eventDataWriteFreq = "))[1]
my_div_control_file[eventDataWriteFreq_line] <- paste0("eventDataWriteFreq = ", eventDataWriteFreq)

# Set frequency in which to print MCMC status to the screen
# printFreq <- 1000 # Print status every 10^3 generations for short test runs
printFreq <- 10000 # Print status every 10^4 generations for long runs
# printFreq <- 100 # Print status every 100 generations
printFreq_line <- which(str_detect(string = my_div_control_file, pattern = "printFreq = "))[1]
my_div_control_file[printFreq_line] <- paste0("printFreq = ", printFreq)

# Whether acceptance/proposal history should be saved.
  # If 1, outputs whether each proposal was accepted. The number identifying the proposal matches the one in the code.
  # The default value is 0 (i.e., do not output this information).
outputAcceptanceInfo <- 0 
outputAcceptanceInfo_line <- which(str_detect(string = my_div_control_file, pattern = "outputAcceptanceInfo = "))[1]
my_div_control_file[outputAcceptanceInfo_line] <- paste0("outputAcceptanceInfo = ", outputAcceptanceInfo)

# Set the path to the acceptance/proposal history file
  # The path of the file to which to write whether each proposal was accepted. 
  # outputAcceptedInfo must be set to 1 for this information to be written. 
  # The default value is acceptance_info.txt.
acceptanceInfoFileName_path <- "acceptance_info.txt"
acceptanceInfoFileName_line <- which(str_detect(string = my_div_control_file, pattern = "acceptanceInfoFileName = "))[1]
my_div_control_file[acceptanceInfoFileName_line] <- paste0("acceptanceInfoFileName = ", acceptanceInfoFileName_path)

# Set frequency in which to update the acceptance rate calculation
  # Acceptance rate = how often new proposal are accepted as the next step in the chain
  # Give information on how efficient is the movement of the chain in the parameter space
  # The acceptance rate is output to both the MCMC data file and print to the screen
acceptanceResetFreq <- 1000 # Reset every 10^3 generations
acceptanceResetFreq_line <- which(str_detect(string = my_div_control_file, pattern = "acceptanceResetFreq = "))[1]
my_div_control_file[acceptanceResetFreq_line] <- paste0("acceptanceResetFreq = ", acceptanceResetFreq)

# Set prefix to add to all output files (separated with "_")
  # Need to be commented out to do not add any prefix
# run_prefix_for_output_files <- "BAMM_Ponerinae_test_run"
# run_prefix_for_output_files <- "BAMM_Ponerinae"
run_prefix_for_output_files <- "BAMM_Ponerinae_MCC"
outName_line <- which(str_detect(string = my_div_control_file, pattern = "outName = "))[1]
my_div_control_file[outName_line] <- paste0("outName = ", run_prefix_for_output_files)
# my_div_control_file[outName_line] <- paste0("# outName = ", run_prefix_for_output_files) # Commented version to add any prefix


### 2.5/ Set the scaling operators = temperatures, to propose new values for sampled parameters ####

# The highest scaling operators = temperatures = the bigger changes can be implemented
  # Advantages = allows to exact suboptimum
  # Cons = may be unstable / harder to reach convergence

# Set scale parameter used for updating the initial speciation rate (lambda0) for each regime/process
  # Updated as lambda0_new ~ lambda0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
updateLambdaInitScale <- 2.0
updateLambdaInitScale_line <- which(str_detect(string = my_div_control_file, pattern = "updateLambdaInitScale = "))[1]
my_div_control_file[updateLambdaInitScale_line] <- paste0("updateLambdaInitScale = ", updateLambdaInitScale)

# Set window size parameter used for updating the rate variation parameter (alpha) of speciation rates for each regime/process
  # Updated as alpha_new ~ alpha_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
updateLambdaShiftScale <- 0.1
updateLambdaShiftScale_line <- which(str_detect(string = my_div_control_file, pattern = "updateLambdaShiftScale = "))[1]
my_div_control_file[updateLambdaShiftScale_line] <- paste0("updateLambdaShiftScale = ", updateLambdaShiftScale)

# Set scale parameter used for updating the initial extinction rate (mu0) for each regime/process
   # Updated as mu0_new ~ mu0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
updateMuInitScale <- 2.0
updateMuInitScale_line <- which(str_detect(string = my_div_control_file, pattern = "updateMuInitScale = "))[1]
my_div_control_file[updateMuInitScale_line] <- paste0("updateMuInitScale = ", updateMuInitScale)

# Set window size parameter used for updating LOCAL moves of the position of shifts on the tree
  # Updated as position_new ~ position_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
  # Unit = fraction of root_to_tip length. May lead to jump of the shift position to a new branch
  # Ex: For a tree of 100My, with parameter set to 0.05, the proposal window for local position change is +/- 5My around the previous value
updateEventLocationScale <- 0.05
updateEventLocationScale_line <- which(str_detect(string = my_div_control_file, pattern = "updateEventLocationScale = "))[1]
my_div_control_file[updateEventLocationScale_line] <- paste0("updateEventLocationScale = ", updateEventLocationScale)

# Set scale parameter used for updating the LAMBDA rate parameter of the Poisson process controlling the number of shifts in the submodels
  # Updated as LAMBDA_new ~ LAMBDA_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
updateEventRateScale <- 4.0
updateEventRateScale_line <- which(str_detect(string = my_div_control_file, pattern = "updateEventRateScale = "))[1]
my_div_control_file[updateEventRateScale_line] <- paste0("updateEventRateScale = ", updateEventRateScale)


### 2.6/ Set the relative frequencies of operator uses (frequency of parameter updates) at each generation ####

# Set the relative frequency of MCMC moves that change the number of events in the submodel (shift from Mk to Mk+1 or Mk-1)
updateRateEventNumber <- 0.1 # 1/10 generations
updateRateEventNumber_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateEventNumber = "))[1]
my_div_control_file[updateRateEventNumber_line] <- paste0("updateRateEventNumber = ", updateRateEventNumber)

# Set the relative frequency of MCMC moves that change the location of an event on the tree (update position parameters)
updateRateEventPosition <- 1 # Every generation !
updateRateEventPosition_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateEventPosition = "))[1]
my_div_control_file[updateRateEventPosition_line] <- paste0("updateRateEventPosition = ", updateRateEventPosition)

# Set the relative frequency of MCMC moves that change the rate at which events occur (update LAMBDA parameter)
updateRateEventRate <- 1 # Every generation !
updateRateEventRate_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateEventRate = "))[1]
my_div_control_file[updateRateEventRate_line] <- paste0("updateRateEventRate = ", updateRateEventRate)

# Set the relative frequency of MCMC moves that change the initial speciation rates (lambda0) associated with a regime
   # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
updateRateLambda0 <- 1
updateRateLambda0_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateLambda0 = "))[1]
my_div_control_file[updateRateLambda0_line] <- paste0("updateRateLambda0 = ", updateRateLambda0)

# Set the relative frequency of MCMC moves that change the exponential shift parameter (alpha) of the speciation rate associated with a regime
   # alpha in lambda(t) = lamba0 x exp(alpha*t)
updateRateLambdaShift <- 1
updateRateLambda0_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateLambda0 = "))[1]
my_div_control_file[updateRateLambda0_line] <- paste0("updateRateLambda0 = ", updateRateLambda0)

# Set the relative frequency of MCMC moves that change the (initial) extinction rate associated with a regime
   # mu0 in mu(t) = mu0 x exp(alpha*t)
   # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
updateRateMu0 <- 1
updateRateMu0_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateMu0 = "))[1]
my_div_control_file[updateRateMu0_line] <- paste0("updateRateMu0 = ", updateRateMu0)

# Set the relative frequency of MCMC moves that flip the time mode (time-constant <=> time-variable)
updateRateLambdaTimeMode <- 0 # By default, only use time-variable mode
updateRateLambdaTimeMode_line <- which(str_detect(string = my_div_control_file, pattern = "updateRateLambdaTimeMode = "))[1]
my_div_control_file[updateRateLambdaTimeMode_line] <- paste0("updateRateLambdaTimeMode = ", updateRateLambdaTimeMode)

# Set the ratio of local to global moves used to propose new location of events on the tree (update position parameters)
localGlobalMoveRatio <- 10.0 # Ten times more local changes than global changes
localGlobalMoveRatio_line <- which(str_detect(string = my_div_control_file, pattern = "localGlobalMoveRatio = "))[1]
my_div_control_file[localGlobalMoveRatio_line] <- paste0("localGlobalMoveRatio = ", localGlobalMoveRatio)


### 2.7/ Set the initial parameter values to start the MCMC chain(s) ####

# The MCMC chain start with a model with no shift (M0 submodel). So the initial parameter values are for this unique regime
# (But probably also for the new regime if added?)

# Run a BD model to obtain credible starting value for the root process
# BD_fit <- phytools::fit.bd(tree = Ponerinae_phylogeny_1534t)
BD_fit <- phytools::fit.bd(tree = Ponerinae_MCC_phylogeny_1534t)

# Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
  # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
lambdaInit0 <- 0.032
lambdaInit0 <- 0.1568
lambdaInit0 <- BD_fit$b
lambdaInit0_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaInit0 = "))[1]
my_div_control_file[lambdaInit0_line] <- paste0("lambdaInit0 = ", lambdaInit0)

# Set the initial shift parameter (alpha) for the root process (regime 0)
  # alpha in lambda(t) = lamba0 x exp(alpha*t)
  # Initial value set to 0 such as the process is a constant rate
lambdaShift0 <- 0
lambdaShift0_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaShift0 = "))[1]
my_div_control_file[lambdaShift0_line] <- paste0("lambdaShift0 = ", lambdaShift0)

# Set the intial extinction rate (mu0) for the first regime starting at the root of the tree (regime 0)
  # mu0 in mu(t) = mu0 x exp(alpha*t)
  # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
muInit0 <- 0.005
muInit0 <- max(0.005, BD_fit$d)
muInit0_line <- which(str_detect(string = my_div_control_file, pattern = "muInit0 = "))[1]
my_div_control_file[muInit0_line] <- paste0("muInit0 = ", muInit0)

# Set the initial number of non-root processes/shifts = M0 submodel
initialNumberEvents <- 0
initialNumberEvents_line <- which(str_detect(string = my_div_control_file, pattern = "initialNumberEvents = "))[1]
my_div_control_file[initialNumberEvents_line] <- paste0("initialNumberEvents = ", muInit0)


### 2.8/ Set the MCMC chain behavior ####

# Set the number of Markov chains to run
  # Each chain will have a different temperature to favor different exploration behavior of the parameter space
numberOfChains <- 4
numberOfChains_line <- which(str_detect(string = my_div_control_file, pattern = "numberOfChains = "))[1]
my_div_control_file[numberOfChains_line] <- paste0("numberOfChains = ", numberOfChains)

# Set the temperature increment parameter that control the difference of temperatures between the chains.
  # This value should be > 0
  # The temperature for the i-th chain is computed as 1 / [1 + deltaT * (i - 1)]
  # Chain 1 is the coldest. Highest chain is the hottest
  # For large trees, over 100 tips, deltaT = 0.05 or deltaT = 0.01 works better
deltaT <- 0.05
deltaT_line <- which(str_detect(string = my_div_control_file, pattern = "deltaT = "))[1]
my_div_control_file[deltaT_line] <- paste0("deltaT = ", deltaT)

# Set the frequency of generations at which to propose a chain swap
  # The coupled-MCMC algorithm will check chain state and eventually swap for the one having reach the highest likelihood
swapPeriod <- 1000 # Check swapping every 10^3 generations
swapPeriod_line <- which(str_detect(string = my_div_control_file, pattern = "swapPeriod = "))[1]
my_div_control_file[swapPeriod_line] <- paste0("swapPeriod = ", swapPeriod)

# Set the path to the file where to store information about each chain swap proposal
  # The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
  # where [generation] is the generation in which the swap proposal was made,
  # [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
  # whether the swap was made. The cold chain has a rank of 1.
chainSwapFileName_path <- "chain_swap_log.txt"
chainSwapFileName_line <- which(str_detect(string = my_div_control_file, pattern = "chainSwapFileName = "))[1]
my_div_control_file[chainSwapFileName_line] <- paste0("chainSwapFileName = ", chainSwapFileName_path)


### 2.9/ Set other parameters ####

# Set the minimum size of a clade to allow a shift to occur
  # Constrain location of possible rate-change events to occur only on branches with at least this many descendant tips.
  # The default value of 1 allows shifts to occur on all branches.
# minCladeSizeForShift <- 1
minCladeSizeForShift <- 3
minCladeSizeForShift_line <- which(str_detect(string = my_div_control_file, pattern = "minCladeSizeForShift = "))[1]
my_div_control_file[minCladeSizeForShift_line] <- paste0("minCladeSizeForShift = ", minCladeSizeForShift)

# Set the "grain" at which time-continuous calculations are discretized
  # The continuous-time change in diversification rates are approximated by breaking each branch into constant-rate diversification segments
  # with each segment given a length determined by the segLength parameter.
  # segLength is in fraction of the root-to-tip distance of the tree.
  # Ex: For an ultrametric tree of 100My, a segLength of 0.02 lead to a step size of 2My
  # If the value is greater than a given branch length BAMM will not break the branch into segments but use the mean rate across the entire branch.
# segLength = 0.02
segLength = 0.01
segLength_line <- which(str_detect(string = my_div_control_file, pattern = "segLength = "))[1]
my_div_control_file[segLength_line] <- paste0("segLength = ", segLength)


### 2.10/ Export the updated custom control file ####
# writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt"))
# writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt"))
writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt"))


##### 3/ Run BAMM with calls to command lines from within r using system() #####

### 3.1/ Run BAMM ####

?system

# Version
system(paste0(BAMM_path,"bamm --version"))

# Test run
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt"))

# Full run
# system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt"))
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt"))

### Outputs
# The run_info.txt file, containing a summary of your parameters/settings
# An mcmc_log.txt containing raw MCMC information useful in diagnosing convergence
# An event_data.txt file containing all of evolutionary rate parameters and their topological mappings
# A chain_swap.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).

### 3.2/ Clean outputs (move them to a dedicated folder) ####

# BAMM_output_folder_path <- "./outputs/BAMM/BAMM_Ponerinae_test_run/"
# BAMM_output_folder_path <- "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/"
BAMM_output_folder_path <- "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/"

# Detect output files
output_files_path <- list.files(path = "./", pattern = "BAMM_")
# Move output files to dedicated folder
file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "my_phy.tree")) # Phylo file
# file.rename(from = paste0(BAMM_path, "my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "my_div_control_file.txt")) # Control file
file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt")) # Control file


### 3.3/ Run multiple runs with different prior for the LAMBDA controlling the expected number of shifts ####

expectedNumberOfShifts_range <- c(1, 5, 10, 20)

## Load control files

# my_div_control_file <- readLines(con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt")) # Test run
# my_div_control_file <- readLines(con = paste0(BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt")) # Full run

# my_div_control_file <- readLines(con = paste0(BAMM_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt")) # Full MCC run
my_div_control_file <- readLines(con = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt")) # Full MCC run


for (i in seq_along(expectedNumberOfShifts_range))
{
  # i <- 2
  
  # Extract expected nb of shifts
  expectedNumberOfShifts_i <- expectedNumberOfShifts_range[i]
  
  # Load control file
  my_div_control_file_i <- my_div_control_file
  
  # Update control file for expectedNumberOfShifts
  expectedNumberOfShifts_line <- which(str_detect(string = my_div_control_file_i, pattern = "expectedNumberOfShifts = "))[1]
  my_div_control_file_i[expectedNumberOfShifts_line] <- paste0("expectedNumberOfShifts = ", expectedNumberOfShifts_i)
  
  # Update control file for output prefix
  # run_prefix_for_output_files <- paste0("BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i) # Test run
  # run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i) # Full run
  run_prefix_for_output_files <- paste0("BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i) # Full MCC run
  outName_line <- which(str_detect(string = my_div_control_file_i, pattern = "outName = "))[1]
  my_div_control_file_i[outName_line] <- paste0("outName = ", run_prefix_for_output_files)
  
  # Write control file
  # writeLines(text = my_div_control_file_i, con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Test run
  # writeLines(text = my_div_control_file_i, con = paste0(BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full run
  writeLines(text = my_div_control_file_i, con = paste0(BAMM_path, "BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full MCC run
  
  # Lauch BAMM run
  # system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Test run
  # system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full run
  system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full MCC run
  
  # Set output folder
  # BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_test_runs/BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i, "/") # Test run
  # BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i, "/") # Full run
  BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_Ponerinae_MCC_full_runs/BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i, "/") # Full MCC run
  
  dir.create(path = BAMM_output_folder_path, recursive = T, showWarnings = F)
  
  # Detect output files
  output_files_path <- list.files(path = "./", pattern = run_prefix_for_output_files)
  # Move output files to dedicated folder
  file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
  # file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "Ponerinae_phylogeny.tree")) # Phylo file
  file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "Ponerinae_MCC_phylogeny.tree")) # Phylo file
  # file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Control file for Test run
  # file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Control file for Full run
  file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Control file for Full MCC run
  
  cat(paste0(Sys.time(), " - BAMM run finished for expectedNumberOfShifts = ",expectedNumberOfShifts_i," - n°", i, "/", length(expectedNumberOfShifts_range),"\n"))
}


##### 4/ Check chain convergence #####

## Compare MCMC chains of runs with different expectedNumberofShifts

### 4.1/ Load the MCMC log files ####

MCMC_log_list <- list()
for (i in seq_along(expectedNumberOfShifts_range))
{
  # i <- 1
  
  # Extract expected nb of shifts
  expectedNumberOfShifts_i <- expectedNumberOfShifts_range[i]
  
  # Set path to output folder
  # BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_test_runs/BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i, "/") # Test run
  # BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i, "/") # Full run
  BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_Ponerinae_MCC_full_runs/BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i, "/") # Full MCC run
  
  # Extract prefix of output files
  # run_prefix_for_output_files <- paste0("BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i) # Test run
  run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i) # Full run
  run_prefix_for_output_files <- paste0("BAMM_Ponerinae_MCC_nbshifts",expectedNumberOfShifts_i) # Full MCC run
  
  # Load the MCMC log file
  MCMC_log_i <- read.csv(paste0(BAMM_output_folder_path, run_prefix_for_output_files, "_mcmc_log.txt"), header = T)
  
  # Add expectedNumberOfShifts
  MCMC_log_i$ExpNbShifts <- expectedNumberOfShifts_i
  
  # Store output
  MCMC_log_list[[i]] <- MCMC_log_i
  
}
  

### 4.2/ Plot traces of the sampled generations ####

# One trace following the main chain as a mixed of chains using swapping

## One plot per BAMM model

par(mfrow = c(2,2))

for (i in seq_along(expectedNumberOfShifts_range))
{
  # i <- 1
  
  # Extract expected nb of shifts
  expectedNumberOfShifts_i <- expectedNumberOfShifts_range[i]
  
  # Extract MCMC log
  MCMC_log_i <- MCMC_log_list[[i]]
  
  # Plot
  plot(MCMC_log_i$logLik ~ MCMC_log_i$generation,
       type = "o",
       main = paste0("expectedNumberOfShifts = ",expectedNumberOfShifts_i))
}


## A single ggplot with all BAMM models

MCMC_log_all_ExpNbShifts <- data.frame()
for (i in seq_along(expectedNumberOfShifts_range))
{
  # i <- 1
  
  # Extract MCMC log
  MCMC_log_i <- MCMC_log_list[[i]]
  
  # Merge in ggplot df
  MCMC_log_all_ExpNbShifts <- rbind(MCMC_log_all_ExpNbShifts, MCMC_log_i)
}
MCMC_log_all_ExpNbShifts$ExpNbShifts <- as.factor(MCMC_log_all_ExpNbShifts$ExpNbShifts)

# GGplot
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MCMC_traces_per_expected_nb_of_shifts.pdf", width = 10, height = 6)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MCMC_traces_per_expected_nb_of_shifts.pdf", width = 10, height = 6)

MCMC_log_ggplot <- ggplot(data = MCMC_log_all_ExpNbShifts) +
  
  geom_line(mapping = aes(y = logLik, x = generation, col = ExpNbShifts),
            linewidth = 2.0, alpha = 0.5) +
  
  scale_color_manual("ExpNbShifts", values = RColorBrewer::brewer.pal(n = 4, name = "Spectral")) +
  
  labs(x = "Generations") +
  
  ggtitle("MCMC Traces per Expected nb of Shifts") +
  
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
  
print(MCMC_log_ggplot)

dev.off()

### 4.4/ Remove burn-in ####

# Choose burn-in according to position of the plateau in the MCMC traces

burn_in <- 0.1
# burn_in <- 0.25

# # Or select burn-in to sample a define number of posterior samples
# 
# # Targeted number of samples
# nb_samples <- 1000
# # Nb of samples in the posterior before removing burn-in
# nb_samples_full <- max(table(MCMC_log_all_ExpNbShifts$ExpNbShifts))
# # Burn-in proportion
# burn_in <- 1 - (nb_samples/nb_samples_full)

# Better to subset the required nb of posterio samples randomly after the removing the burn-in
# Will improve computation of ESS

# Remove burn-in from MCMC logs
post_burn_MCMC_log_list <- MCMC_log_list
for (i in seq_along(post_burn_MCMC_log_list))
{
  # i <- 1
  
  # Extract MCMC log
  MCMC_log_i <- MCMC_log_list[[i]]
  
  # Set threshold
  burn_in_threshold <- ceiling(burn_in * MCMC_log_i$generation[nrow(MCMC_log_i)])
  
  # Remove burn-in
  post_burn_MCMC_log_i <- MCMC_log_i[MCMC_log_i$generation >= burn_in_threshold, ]
  
  # Store updated MCMC log
  post_burn_MCMC_log_list[[i]] <- post_burn_MCMC_log_i
}


### 4.5/ Explore effective sample sizes ####

ESS_df <- NULL
for (i in seq_along(post_burn_MCMC_log_list))
{
  # i <- 1
  
  # Extract post burn-in MCMC log
  post_burn_MCMC_log_i <- post_burn_MCMC_log_list[[i]]
  
  # Store SSE info
  ESS_df_i <- data.frame(ExpNbShifts = post_burn_MCMC_log_i$ExpNbShifts[1],
                         SSE_N_shifts = coda::effectiveSize(post_burn_MCMC_log_i$N_shifts),
                         SSE_eventRate = coda::effectiveSize(post_burn_MCMC_log_i$eventRate),
                         SSE_logLik = coda::effectiveSize(post_burn_MCMC_log_i$logLik))
  
  # Merge
  ESS_df <- rbind(ESS_df, ESS_df_i)
  row.names(ESS_df) <- 1:nrow(ESS_df)
}

ESS_df

# Proper ESS > 200


### 4.6/ Compare prior and posterior distribution of LAMBDA (parameter controlling the nb of rates) ####

# Prior distribution is a Poisson, so we use a bar graph for a discrete distribution

# Compare prior and posterior of the different runs with different expectedNumberofShifts

par(mfrow = c(2,2), mar = c(6.1, 4.1, 4.1, 2.1))

for (i in seq_along(expectedNumberOfShifts_range))
{
  # i <- 1
  
  # Extract expected nb of shifts
  expectedNumberOfShifts_i <- expectedNumberOfShifts_range[i]
  
  # Extract MCMC log
  MCMC_log_i <- MCMC_log_list[[i]]
  
  # Plot
  plotPrior(mcmc = MCMC_log_i,
            expectedNumberOfShifts = expectedNumberOfShifts_i,
            burnin = burn_in,
            main = paste0("expectedNumberOfShifts = ",expectedNumberOfShifts_i))
}

# Merge all BAMM runs to explore overall posterior probabilities
post_burn_MCMC_log_all <- bind_rows(post_burn_MCMC_log_list)
table(post_burn_MCMC_log_all$N_shifts) / nrow(post_burn_MCMC_log_all) * 100
summary(post_burn_MCMC_log_all$N_shifts) # Posterior distribution of the Nb of shifts
summary(post_burn_MCMC_log_all$eventRate) # Posterior distribution of LAMBDA = parameter controlling the expected nb of shifts


## Choose the expectedNumberofShifts the closest to the overall median
selected_expectedNumberofShifts <- 10


##### 5/ Load and format selected model output files ####

?BAMMtools::getEventData

### 5.1/ Load BAMM output data ####

# Choose the expectedNumberofShifts the closest to the overall median
selected_expectedNumberofShifts <- 10

# Set path to output folder
# BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_test_runs/BAMM_Ponerinae_test_run_nbshifts",selected_expectedNumberofShifts, "/") # Test run
# BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",selected_expectedNumberofShifts, "/") # Full run
BAMM_output_folder_path <- paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_Ponerinae_MCC_full_runs/BAMM_Ponerinae_MCC_nbshifts",selected_expectedNumberofShifts, "/") # Full MCC run

# Extract prefix of output files
# run_prefix_for_output_files <- paste0("BAMM_Ponerinae_test_run_nbshifts",selected_expectedNumberofShifts) # Test run
# run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",selected_expectedNumberofShifts) # Full run
run_prefix_for_output_files <- paste0("BAMM_Ponerinae_MCC_nbshifts",selected_expectedNumberofShifts) # Full MCC run

# Load the phylogeny
# phy <- read.tree(paste0(BAMM_output_folder_path, "Ponerinae_phylogeny.tree")) # Full run
phy <- read.tree(paste0(BAMM_output_folder_path, "Ponerinae_MCC_phylogeny.tree")) # Full MCC run

# Create the bammdata summarizing BAMM outputs
BAMM_data_output <- getEventData(phy = phy, 
                            eventdata = paste0(BAMM_output_folder_path, run_prefix_for_output_files, "_event_data.txt"),
                            burnin = burn_in,
                            type = "diversification")

### 5.2/ Select the subset of posterior samples ####

# Set the targeted number of posterior samples
nb_samples <- 1000

# Get a subset of a selected number of posterior samples
set.seed(seed = 1234)
sample_indices <- sample(x = 1:length(BAMM_data_output$eventData), size = nb_samples)
BAMM_posterior_samples_data <- subsetEventData(BAMM_data_output, index = sample_indices)

# Save the BAMM posterior samples object
# saveRDS(BAMM_posterior_samples_data, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
saveRDS(BAMM_posterior_samples_data, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

### 5.3/ Explore the subset of posterior samples ####

# Basically an extended ape phylogeny object storing BAMM output data
str(BAMM_posterior_samples_data, max.level = 1)
BAMM_posterior_samples_data$begin # Absolute time since root of edge/branch start
BAMM_posterior_samples_data$end # Absolute time since root of edge/branch end
BAMM_posterior_samples_data$numberEvents # Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts
BAMM_posterior_samples_data$eventData # List of dataframes recording shift events in each posterior configuration
BAMM_posterior_samples_data$eventData[[1]] # Dataframe recording shift events and macroevolutionary regimes in the focal posterior configuration. 1st line = Background root regime
BAMM_posterior_samples_data$eventData[[1]]$node # Tipward node ID of the branch where the shift occured
BAMM_posterior_samples_data$eventData[[1]]$time # Absolute time since root when the shift occurred on the branch
BAMM_posterior_samples_data$eventData[[1]]$lam1 # Initial rate of speciation (lambda0) of the regime
BAMM_posterior_samples_data$eventData[[1]]$lam2 # Speciation rate change parameter (alpha/z) of the regime. If = 0, rate is constant. If < 0 = decay. If > 0 = growth
BAMM_posterior_samples_data$eventData[[1]]$mu1 # Initial rate of extinction (mu0) of the regime. Should be the constant extinction rate of the regime.
BAMM_posterior_samples_data$eventData[[1]]$mu2 # Extinction rate change parameter (alpha/z) of the regime. Should be fixed to 0 as extinction rates are constant in BAMM.
BAMM_posterior_samples_data$eventVectors # List of integer vectors of regime membership per branches in each posterior configuration
BAMM_posterior_samples_data$eventVectors[[1]] # Integer vectors of regime membership per branches
BAMM_posterior_samples_data$eventBranchSegs[[1]] # Same but with matrix including tipward node ID (NOT the edge ID) and begin/end ages of the branches
BAMM_posterior_samples_data$tipStates[[1]] # Integer vectors of regime membership per tips
BAMM_posterior_samples_data$tipLambda[[1]] # Integer vectors of final speciation rates at tips = current speciation rates
BAMM_posterior_samples_data$meanTipLambda # Mean current tip speciation rates across all posterior configurations. Better to use the median and 95% HPD
BAMM_posterior_samples_data$tipMu[[1]] # Integer vectors of final extinction rates at tips = current extinction rates
BAMM_posterior_samples_data$tipLambda[[1]] - BAMM_posterior_samples_data$tipMu[[1]] # Integer vectors of final diversification  rates at tips = current diversification rates
BAMM_posterior_samples_data$meanTipMu # Mean current tip extinction rates across all posterior configurations. Better to use the median and 95% HPD

# Explore number of regimes
summary(BAMM_posterior_samples_data$numberEvents)
BayesTwin::HPD(sample = BAMM_posterior_samples_data$numberEvents, cred_int = 0.95)
sd(BAMM_posterior_samples_data$numberEvents)

##### 6/ Plot regime shifts #####

# Marginal probabilities. If high for a sequence of branches, good practice could be to look for in the posterior samples for the joint probability of a shift happening on 0/1/2/… of these branches in a single sample.
# Goal: avoid misinterpreting those independent marginal probabilities has evidence for a sequence of shifts as what it actually supports is the likely presence of one shift among those sequence of branches, whit a bit of uncertainty on the exact location among those successive branches

# Caution: showing all marginal PP of shift location may be biases since shift location may not be independent
# The presence of one shift may be associated with the absence of another. Be careful with the interpretation of marginal PP.
# Try to investigate joint probabilities of some events.

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

### 6.1/ Plot the set of credible shift configurations ####

# To extract the set of credible shift configuration found in the posterior samples
# Extract configurations until cum prob reach 95%
# Probability of shift as proportion in the posterior samples
# Use a threshold to select only shifts with a minimal marginal odd ratio that depict substantially higher PP than prior probability = core-shifts.
# If not ignored, many configurations will contain non-core shifts that are not more expected than according to the prior, thus meaningless to model our data

?credibleShiftSet
?plot.credibleshiftset

cset <- credibleShiftSet(BAMM_posterior_samples_data,
                         expectedNumberOfShifts = selected_expectedNumberofShifts, # Initial parameter used in the config file to set the exponential prior of the lambda parameter for number of shifts
                         threshold = 5) # Significance threshold of odd-ratios between prior and posterior probabilities of rate shift on branches
# Display summary with PP of each configuration and associated number of core shifts (with odd-ratio > threshold)
summary(cset)
# 800 configurations needed to reach 95% in cumulative probability!
# MAP with only PP of only 1.8% 

# Plot the set of credible shift configurations
plot.credibleshiftset(cset, lwd = 2.5)

### 6.2/ Scale phylogeny according to the Marginal Probability of shift ####

?marginalShiftProbsTree

branch_marg_posterior_probs <- marginalShiftProbsTree(BAMM_posterior_samples_data)
branch_marg_posterior_probs$edge.length

par(mfrow = c(1,1))
plot.phylo(branch_marg_posterior_probs)

# Highlight credible location of shifts
# Not independent! Successive high values do not mean that shifts are likely on all branches
# Rather they indicate uncertainty in the location of one single shift

## Save phylogeny scaled according to the marginal posterior probability of shift 
saveRDS(object = branch_marg_posterior_probs, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/branch_marg_posterior_probs.rds")

### 6.3/ Compute marginal posterior odd-ratios ####

?getBranchShiftPriors
?marginalOddsRatioBranches

## Compare Marginal posterior probability of shift with Prior probabilities

# Because the expected number of rate shifts under the prior among all branches is uniform, 
# we expect to observe more shifts on long branches than short branches, just by chance alone.
#	Thus, a good practice to find support for a rate shift is not only to look at the marginal probability of rate shift,
# BUT to account for a substantial change between the prior probability of shift on branch
# (which depend on the prior and the branch length) and the posterior probability observed in the sample
# Odds-ratio = an estimate of the “density” of shifts on a particular branch, independent of the length of the branch

# Get prior probabilities of shift per branches
branch_prior_probs <- getBranchShiftPriors(BAMM_posterior_samples_data, expectedNumberOfShifts = selected_expectedNumberofShifts)
branch_prior_probs$edge.length
# Should be directly proportional to branch length

# Branch length are scaled to prior probabilities of shift per branches
sum(BAMM_posterior_samples_data$edge.length)
sum(branch_prior_probs$edge.length)

# Compute marginal odds ratio to account for prior probabilities
# MO =  posterior proba of rate shift / prior proba of rate shift. 
# The highest, the more substantial is the support for a shift occurring on that branch.
branch_marginal_odd_ratios <- marginalOddsRatioBranches(ephy = BAMM_posterior_samples_data, expectedNumberOfShifts = selected_expectedNumberofShifts)

# Branch length are scaled to odd-ratios of rate shift per branch
branch_marginal_odd_ratios$edge.length
summary(branch_marginal_odd_ratios$edge.length)
sum(branch_marginal_odd_ratios$edge.length)
table(branch_marginal_odd_ratios$edge.length > 5) # OR > 5 = typical threshold for significance = core shift
table(branch_marginal_odd_ratios$edge.length > 50) # Threshold needed to obtain equivalent number of core shifts vs. posterior expected number of shift (LAMBDA)

## Save phylogeny scaled according to the Odd-Ratio of marginal posterior probability / prior probability of shift (proportional to branch length)
saveRDS(object = branch_marginal_odd_ratios, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/branch_marginal_odd_ratios.rds")

## Compare prior, posterior, and Odd-ratios
plot(branch_marg_posterior_probs) # Posterior probabilities
plot(branch_prior_probs) # Prior probabilities
plot(branch_marginal_odd_ratios) # Odd-ratios = Posterior probabilities / Prior probabilities

# Save prior, posterior, and odd-ratio in a single df
edge_rate_shift_probs_df <- data.frame(edge_ID = 1:nrow(branch_marg_posterior_probs$edge),
                                       edge_length = phy$edge.length,
                                       marg_posterior_probs = branch_marg_posterior_probs$edge.length,
                                       prior_probs = branch_prior_probs$edge.length,
                                       marginal_odd_ratios = branch_marginal_odd_ratios$edge.length)
# saveRDS(edge_rate_shift_probs_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/edge_rate_shift_probs_df.rds")
saveRDS(edge_rate_shift_probs_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/edge_rate_shift_probs_df.rds")

### 6.4/ Map bioregion membership on phylogeny scaled by shift probability ####

# Load treedata with bioregion membership information
# Ponerinae_phylogeny_1534t_treedata_for_ARE <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata_for_ARE.rds")
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE.rds")

# Assign a unique bioregion per edge
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
# Bioregion_binary_df <- Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Afrotropics", "mean_state_Australasia", "mean_state_Eastern Palearctic", "mean_state_Indomalaya", "mean_state_Nearctic", "mean_state_Neotropics", "mean_state_Western Palearctic")]
Bioregion_binary_df <- Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo[, c("mean_state_Afrotropics", "mean_state_Australasia", "mean_state_Eastern Palearctic", "mean_state_Indomalaya", "mean_state_Nearctic", "mean_state_Neotropics", "mean_state_Western Palearctic")]
Bioregion_binary_df <- bioregion_names[unlist(apply(X = Bioregion_binary_df, MARGIN = 1, FUN = which.max))]
# Bioregion_binary_df <- as.data.frame(cbind(Ponerinae_phylogeny_1534t_treedata@extraInfo$node, Bioregion_binary_df))
Bioregion_binary_df <- as.data.frame(cbind(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo$node, Bioregion_binary_df))
names(Bioregion_binary_df) <- c("node", "Bioregion")
Bioregion_binary_df$node <- as.numeric(Bioregion_binary_df$node)
Bioregion_binary_df$Bioregion <- factor(x = Bioregion_binary_df$Bioregion, levels = bioregion_names, labels = bioregion_names)

# Ponerinae_phylogeny_1534t_treedata_for_ARE@extraInfo <- left_join(Ponerinae_phylogeny_1534t_treedata_for_ARE@extraInfo, y = Bioregion_binary_df)
Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo <- left_join(Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE@extraInfo, y = Bioregion_binary_df)

# Save updated treedata with unique Bioregion membership (the most likely one)
# saveRDS(object = Ponerinae_phylogeny_1534t_treedata_for_ARE, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata_for_ARE.rds")
saveRDS(object = Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE, file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE.rds")


# Load Genus-groups metadata
# all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/all_genus_groups_metadata_df.rds")
all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_all_genus_groups_metadata_df.rds")

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]

# Set color scheme for areas/bioregions (Use the BSM color scheme)
# colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
areas_list <- c("A", "U", "E", "I", "R", "N", "W")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
colors_list_for_areas <- colors_list_for_areas[areas_list]
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Load Missing_taxa_list
Missing_taxa_list <- readRDS(file = "./outputs/Grafting_missing_taxa/Missing_taxa_list.rds")

# Set color scheme for tips
# nb_tips <- length(Ponerinae_phylogeny_1534t_treedata@phylo$tip.label)
# tips_color <- c("black", "red")[Ponerinae_phylogeny_1534t_treedata@extraInfo$missing_node + 1]
# tips_color <- tips_color[1:nb_tips]
nb_tips <- length(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label)
tips_color <- c("black", "red")[as.numeric(Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label %in% Missing_taxa_list$Taxa_name) + 1]
names(tips_color) <- Ponerinae_MCC_phylogeny_1534t_treedata@phylo$tip.label

# Initiate new treedata
# Ponerinae_phylogeny_scaled_PP_shifts <- Ponerinae_phylogeny_1534t_treedata 
Ponerinae_phylogeny_scaled_PP_shifts <- Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE 

# Update branch length for marginal posterior probability of rate shift
Ponerinae_phylogeny_scaled_PP_shifts@phylo$edge.length <- branch_marg_posterior_probs$edge.length

# Plot PDF to aggregate all bioregions
# pdf(file = paste0("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_phylogeny_scaled_PP_shifts.pdf"), height = 16, width = 17)
pdf(file = paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_phylogeny_scaled_PP_shifts.pdf"), height = 16, width = 17)


# Initiate plot
Ponerinae_phylogeny_scaled_PP_shifts_plot <- ggtree(Ponerinae_phylogeny_scaled_PP_shifts,
                                                    color = NA,
                                                    # color = "grey90",
                                                    layout = "rectangular") +
  
  coord_cartesian(clip = "off") +
  
  # Add scale
  geom_treescale(linesize = 1.5, offset = 12.0, offset.label = -18.0, width = 0.5,
                 label = "PP shift", fontsize = 8) +
  
  # # Add tip label
  # geom_tiplab(# mapping = aes(colour = missing_node),
  #   color = tips_color,
  #   align = F,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
  #   linetype = 0, 
  #   size = 2, 
  #   offset = 0.05,   # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
  #   hjust = 0) +    # To center, left-align, right-aligned text
  
  # Add title
  ggtitle(label = paste0("Posterior probability of diversification rate shift")) 

## Loop per bioregions
for (i in seq_along(bioregion_names))
{
  # i <- 1
  
  # Extract bioregion
  bioregion_i <- bioregion_names[i]
  
  # Find index in metadata df
  bioregion_ID_i <- which(names(Ponerinae_phylogeny_scaled_PP_shifts@extraInfo) == paste0("mean_state_",bioregion_i))
  
  Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
    
    geom_tree(color = colors_list_for_areas[i],
              alpha = Ponerinae_phylogeny_scaled_PP_shifts@extraInfo[, bioregion_ID_i, drop = T],
              linewidth = 1.5)
}

# Add Genus-group cladelabels
for (i in 1:nrow(all_genus_groups_metadata_df))
{
  Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
    
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = all_genus_groups_metadata_df$MRCA_node_ID[i],
                    # label = all_genus_groups_metadata_df$group_name[i], 
                    label = paste0('bold(',all_genus_groups_metadata_df$group_name[i],')'), parse = T,
                    align = T,         # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',     # Tip of label, with or without rectangle
                    fontsize = 7,      # Size of the text
                    fill = NA,         # Color of the label background
                    color = c(genus_groups_colors[i], genus_groups_colors[i]),   # Color of the bar and the text
                    angle = 270,
                    hjust = 'center',
                    # offset = c(-2.6, -2.6, -1.6, -1.1, -0.9, -1.1, 0.0)[i],       # For rough phylogeny
                    # offset.text = c(0.10, 0.03, 0.05, 0.03, 0.05, 0.10, 0.05)[i], # For rough phylogeny
                    offset = c(-1.38, -0.97, -0.12, -0.8, -0.05, -0.9, 0.0)[i],         # For MCC phylogeny
                    offset.text = c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)[i],   # For MCC phylogeny
                    barsize = 4)    # Width of the bar
}
  
# Add transparent tree with Bioregion colors (needed for the legend)
Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
  geom_tree(data = Ponerinae_phylogeny_scaled_PP_shifts,
            mapping = aes(col = Bioregion),
            alpha = 0.0, linewidth = 5.0,
            layout = "rectangular") +

  # Add Bioregion legend
  scale_color_manual("Bioregion", breaks = bioregion_names, labels = bioregion_names, values = colors_list_for_areas) +
  guides(color = guide_legend(override.aes = list(alpha = 1.0)))

# Adapt margins
Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
  theme(plot.margin = unit(c(0, 0, 20, 0), "mm"), # trbl?
        # legend.position = c(0.85, 0.2),
        legend.position = c(0.15, 0.75),
        legend.title = element_text(size = 24, face = "bold", margin = margin(b = 25)),
        legend.text = element_text(size = 18, margin = margin(l = 15)),
        legend.key.spacing.y = unit(1.5, 'line'),
        plot.title = element_text(face = "bold", colour = "black", size = 30, hjust = 0.5,
                                  # vjust = -10,
                                  vjust = 0,
                                  angle = 0),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(color = "transparent", fill = NA),
        panel.background = element_rect(color = NA, fill = "transparent"))

print(Ponerinae_phylogeny_scaled_PP_shifts_plot)

dev.off()



##### 7/ Identify the most likely configuration ####

### 7.1/ Maximum A posteriori Probability (MAP) configuration ###

## Get the most frequent configuration = maximum a posteriori probability (MAP) shift configuration = the single configuration of shift location showing up the most in the posterior sample

# Idea: find the posterior configuration that is the most frequent in the posterior samples

MAP <- getBestShiftConfiguration(BAMM_posterior_samples_data,
                                 expectedNumberOfShifts = selected_expectedNumberofShifts,
                                 threshold = 5) # Odd-ratio threshold used to select credible shifts
plot.bammdata(MAP, lwd = 1.25)
addBAMMshifts(MAP, cex=2)

# Only found in 18 posterior samples across 1000!
str(cset$indices, max.level = 1)

### 7.2/ Maximum Shift Credibility configuration (MSC) ####

## Get the most likely configuration = maximum shift credibility configuration (MSC)

# Idea: find the posterior configuration that has the highest probability according to the marginal probability of the branches and number of shifts

# Useful for big phylogenies when there are no clear MAP because of the many possible configurations

# Identify the MSC configs
MSC_detection <- maximumShiftCredibility(BAMM_posterior_samples_data)
# Extract the MSC config
MSC_tree <- subsetEventData(BAMM_posterior_samples_data, index = MSC_detection$sampleindex)
plot.bammdata(MSC_tree, lwd = 1.25)
addBAMMshifts(MSC_tree, par.reset = FALSE, cex = 2)

# Save the MSC tree
# saveRDS(object = MSC_tree, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_tree.rds")
saveRDS(object = MSC_tree, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_tree.rds")

### 7.3/ Extract metadata on rate shifts in the MSC ####

root_age <- max(phytools::nodeHeights(phy))

# Extract event table
MSC_shifts_df <- BAMM_posterior_samples_data$eventData[[MSC_detection$sampleindex]][, -6]
names(MSC_shifts_df) <- c("tipward_node", "start_time_since_root", "lambda0", "alpha", "mu", "regime_ID")

# Compute start age
MSC_shifts_df$start_age <- root_age - MSC_shifts_df$start_time_since_root

# Get associated branch ID
MSC_shifts_df$edge_ID <- match(x = MSC_shifts_df$tipward_node, table = phy$edge[, 2])
# Extract rootward node ID
MSC_shifts_df$rootward_node_ID <- phy$edge[MSC_shifts_df$edge_ID, 1]

# Get the associated probabilities/odd-ratios
MSC_shifts_df$branch_prior_prob <- branch_prior_probs$edge.length[MSC_shifts_df$edge_ID]
MSC_shifts_df$branch_marg_posterior_prob <- branch_marg_posterior_probs$edge.length[MSC_shifts_df$edge_ID]
MSC_shifts_df$branch_odd_ratio <- branch_marginal_odd_ratios$edge.length[MSC_shifts_df$edge_ID]

# Identify dynamics
MSC_shifts_df$trend <- c("decrease", "increase")[as.numeric(MSC_shifts_df$alpha > 0)+1]
MSC_shifts_df$trend_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$alpha > 0)+1]

## Identify previous regime before shift

# Extract edge regime membership
MSC_eventBranchSegs_df <- as.data.frame(BAMM_posterior_samples_data$eventBranchSegs[[MSC_detection$sampleindex]])
names(MSC_eventBranchSegs_df) <- c("tipward_node", "start_time_since_root", "end_time_since_root", "regime_ID")
MSC_eventBranchSegs_df <- MSC_eventBranchSegs_df[MSC_eventBranchSegs_df$tipward_node %in% MSC_shifts_df$tipward_node, ]

# Detect start/ends of regimes
MSC_eventBranchSegs_df <- MSC_eventBranchSegs_df %>% 
  arrange(regime_ID, start_time_since_root) %>%
  group_by(regime_ID) %>%
  mutate(detect_type = row_number()) %>%
  mutate(detect_type = (detect_type > 1) + 1) %>%
  mutate(regime_status = c("start", "end")[detect_type]) %>%
  dplyr::select(-detect_type) %>%
  arrange(tipward_node, regime_ID)

# Find regime with end mathcing with the start of the focal regime
MSC_regime_start_df <- MSC_eventBranchSegs_df %>% 
  filter(regime_status == "start")
MSC_regime_end_df <- MSC_eventBranchSegs_df %>% 
  rename(previous_regime_ID = regime_ID)
MSC_regime_start_end_df <- MSC_regime_start_df %>%
  left_join(y = MSC_regime_end_df[, c("tipward_node", "end_time_since_root", "previous_regime_ID")],
            by = join_by(tipward_node == tipward_node, start_time_since_root == end_time_since_root))

# Join previous regime ID
MSC_shifts_df <- MSC_shifts_df %>% 
  left_join(y = MSC_regime_start_end_df[, c("regime_ID", "previous_regime_ID")])

# Compute rate before/at shift (final rate of previous regime along the focal branch)
previous_regimes_df <- MSC_shifts_df %>% 
  dplyr::select(regime_ID, previous_regime_ID, start_time_since_root) %>%
  rename(regime_shift_time = start_time_since_root,
         current_regime_ID = regime_ID) %>%
  # Extract previous regime parameters
  left_join(y = MSC_shifts_df[, c("regime_ID", "start_time_since_root", "lambda0", "alpha", "mu")],
            join_by(previous_regime_ID == regime_ID))  %>%
  rename(start_time_since_root_previous_regime = start_time_since_root,
         lambda0_previous_regime = lambda0,
         previous_alpha = alpha,
         mu_previous_regime = mu) %>%
  # Compute time length of previous regime until the regime shift
  mutate(time_previous_regime = regime_shift_time - start_time_since_root_previous_regime) %>%
  # Compute rate before/at shift
  mutate(previous_rate = lambda0_previous_regime * exp(previous_alpha * time_previous_regime))

# Join previous regime rate at regime shift
MSC_shifts_df <- MSC_shifts_df %>% 
  left_join(y = previous_regimes_df[, c("current_regime_ID", "previous_rate", "previous_alpha")],
            join_by(regime_ID == current_regime_ID))

# Compute and identify rate shift type by comparing last rate of previous regime and new rate of new regime
MSC_shifts_df$rate_shift <- MSC_shifts_df$lambda0 - MSC_shifts_df$previous_rate
MSC_shifts_df$rate_shift_type <- c("decrease", "increase")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]
MSC_shifts_df$rate_shift_type_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]

# Compute and identify rate trend shift by comparing previous regime trend (alpha) and new regime trend
MSC_shifts_df$alpha_shift <- MSC_shifts_df$alpha - MSC_shifts_df$previous_alpha
MSC_shifts_df$alpha_shift_type <- c("decrease", "increase")[as.numeric(MSC_shifts_df$alpha_shift > 0)+1]
MSC_shifts_df$alpha_shift_type_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$alpha_shift > 0)+1]

MSC_shifts_df

# ## Identify shift types by comparing previous branch mean rate to current branch mean rate
# # Extract mean branch rates
# div_rates_tree <- getMeanBranchLengthTree(BAMM_posterior_samples_data,
#                                           rate = "ndr") # Net Diversification Rates
# MSC_shifts_df$mean_branch_rate <- div_rates_tree$phy$edge.length[MSC_shifts_df$edge_ID]
# MSC_shifts_df$rootward_node_ID <- phy$edge[MSC_shifts_df$edge_ID, 1]
# MSC_shifts_df$previous_branch_ID <- match(x = MSC_shifts_df$rootward_node_ID, table = phy$edge[, 2])
# MSC_shifts_df$previous_branch_mean_rate <- div_rates_tree$phy$edge.length[MSC_shifts_df$previous_branch_ID]
# MSC_shifts_df$rate_shift <- MSC_shifts_df$mean_branch_rate - MSC_shifts_df$previous_branch_mean_rate
# MSC_shifts_df$rate_shift_type <- c("decrease", "increase")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]
# MSC_shifts_df$rate_shift_type_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]
  
# Identify clades/paraphyletic groups associated with regimes

# Plot phylogeny with rate shift location
# pdf("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_tree_with_labels_rect.pdf", width = 20, height = 200)
pdf("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_tree_with_labels_rect.pdf", width = 20, height = 200)

plot.phylo(phy)
#nodelabels(node = MSC_shift_nodes, pch = 21, col = "black", bg = "red", cex = 1.5)
edgelabels(edge = MSC_shifts_df$edge_ID, pch = 21, col = "black", bg = MSC_shifts_df$rate_shift_type_col[-1], cex = MSC_shifts_df$branch_marg_posterior_prob * 3)

dev.off()

View(MSC_shifts_df)

# # Provide name to the regimes for rough tree
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1535] <- "Root_process"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2668] <- "Ponera_group"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2615] <- "Hypoponera_fast_subgroup1"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2491] <- "Hypoponera_fast_subgroup2"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2426] <- "Hypoponera_fast_subgroup3"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2352] <- "Plectroctena_Genus"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1537] <- "PlOHH"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1763] <- "Myopias_Leptogenys"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1767] <- "Leptogenys_subgroup"
# MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1705] <- "Odontomachus_subgroup"

# Provide name to the regimes for MCC tree
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1535] <- "Platythyrea"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1536] <- "Ponerini_background"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1540] <- "Trap_jaw_ants"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1749] <- "Relict_group"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1762] <- "Leptogenys_subgroup"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2393] <- "Hypoponera_Afrotropics"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2819] <- "Diacamma"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2126] <- "Myopias"


# Save metadata on rate shifts in the MSC tree
# saveRDS(object = MSC_shifts_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_shifts_df.rds")
saveRDS(object = MSC_shifts_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_shifts_df.rds")

# Export MSC shift table
MSC_shifts_df_extract <- MSC_shifts_df %>% 
  dplyr::select(regime_ID, regime_group_name, start_age, lambda0, alpha, mu, branch_prior_prob, branch_marg_posterior_prob, branch_odd_ratio)
write.xlsx(x = MSC_shifts_df_extract, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_shifts_df.xlsx")

##### 8/ Plot branch diversification rates on the phylogeny ####

?plot.bammdata
?BAMMtools::dtRates # To compute discretized rates on small segments along branches

### 8.1/ Extract overall rates ####

# Get overall speciation/extinction rates
all_branch_div_rates <- getCladeRates(BAMM_posterior_samples_data)

# Compute overall rate as the mean of mean branch rates
mean(all_branch_div_rates$lambda - all_branch_div_rates$mu)
median(all_branch_div_rates$lambda - all_branch_div_rates$mu)
quantile(all_branch_div_rates$lambda - all_branch_div_rates$mu, c(0.05, 0.95)) # 95% CI interval of the mean branch rate across the whole phylogeny

### 8.2/ Scale branch length to diversification rates ####

?getMeanBranchLengthTree

div_rates_tree <- getMeanBranchLengthTree(BAMM_posterior_samples_data,
                                          rate = "ndr") # Net Diversification Rates

# Explore mean branch rates
div_rates_tree$phy$edge.length

# Explore overall rates
div_rates_tree$mean
div_rates_tree$median

# Plot mean rates
par(mfrow = c(1,1))
plot(div_rates_tree$phy)

# Save phylogeny scaled with mean diversification rates
# saveRDS(div_rates_tree, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/div_rates_tree.rds")
saveRDS(div_rates_tree, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/div_rates_tree.rds")

### 8.3/ Plot mean sliding rates + credible shift locations ####

## Good practice would be to use median rather than mean
# Median can be computed from BAMM_output$eventData or using BAMMtools::dtRates to compute discretize rates

# Use the MSC configuration for plotting the shifts
# Load metadata on rate shifts in the MSC tree
# MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_shifts_df.rds")
MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_shifts_df.rds")

## Plot histogram of different color palette settings for selection

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Histograms_branch_rates.pdf", height = 12, width = 6)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Histograms_branch_rates.pdf", height = 12, width = 6)

par(mfrow = c(4,1), mar = c(3.0, 1.5, 1.0, 0.5), xpd = FALSE) # bltr

# Fully linear
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       tau = 0.001, breaksmethod = 'linear',
                                       show = F)
ratesHistogram(branch_div_rates_plot, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Linear', cex.main = 1, line = -1)

# Linear with a max bound
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       tau = 0.001, breaksmethod = 'linear',
                                       color.interval = c(NA, 0.3),
                                       show = F)
ratesHistogram(branch_div_rates_plot, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Linear with max interval', cex.main = 1, line = -1)

# Quantile method to discretize rates in categories of equal frequencies
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       tau = 0.001, breaksmethod = 'quantile',
                                       show = F)
ratesHistogram(branch_div_rates_plot, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Quantile method', cex.main = 1, line = -1)

# Jenks method to discretize rates in categories min/max variance within/between groups
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       tau = 0.001, breaksmethod = 'jenks',
                                       show = F)
ratesHistogram(branch_div_rates_plot, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Jenks method', cex.main = 1, line = -1)

dev.off()


## Plot branch rates on the phylogeny: phylogram

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_diversity_dynamics_rect.pdf", height = 10, width = 6)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_diversity_dynamics_rect.pdf", height = 10, width = 6)

par(mfrow = c(1,1), xpd = TRUE)

# tau = length of segments used to discretize the continuous rates
# tau = 0.001 => each segment is 1/1000 of the tree height (root age)

# Quantile method to discretize rates in categories of equal frequencies
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       method = "phylogram",
                                       tau = 0.001,
                                       # breaksmethod = 'quantile',
                                       breaksmethod = 'jenks',
                                       lwd = 2)

# Use the configuration from the Maximum Credibility Shift (MSC)
addBAMMshifts(BAMM_posterior_samples_data, index = MSC_detection$sampleindex,
              method = "phylogram",
              col = "black", 
              # bg = MSC_shifts_df$trend_col[-1], # Color shift by trend fo the new regime
              # bg = MSC_shifts_df$rate_shift_type_col[-1], # Color shift by discrete rate shift on the edge
              bg = "grey20",
              par.reset = FALSE,
              cex = MSC_shifts_df$branch_marg_posterior_prob[-1] * 3)
title(main = 'Ponerinae diversity dynamics', cex = 2, line = 0)
text(x = 0.8, y = 400, labels = "Net Div. Rates", font = 2)
addBAMMlegend(branch_div_rates_plot, location = c(-1.0, 0, 140, 350), font = 2)

dev.off()


## Plot branch rates on the phylogeny: polar

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_diversity_dynamics_polar.pdf", height = 7, width = 6)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_diversity_dynamics_polar.pdf", height = 7, width = 6)

par(mfrow = c(1,1), xpd = TRUE)

# tau = length of segments used to discretize the continuous rates
# tau = 0.001 => each segment is 1/1000 of the tree height (root age)

# Quantile method to discretize rates in categories of equal frequencies
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       method = "polar",
                                       tau = 0.001, 
                                       # breaksmethod = 'quantile',
                                       breaksmethod = 'jenks',
                                       lwd = 2)

# Need to identify a configuration with the credible shifts!
addBAMMshifts(BAMM_posterior_samples_data, index = MSC_detection$sampleindex,
              method = "polar",
              col = "black",               
              # bg = MSC_shifts_df$trend_col[-1], # Color shift by trend fo the new regime
              # bg = MSC_shifts_df$rate_shift_type_col[-1], # Color shift by discrete rate shift on the edge
              bg = "grey20",
              par.reset = FALSE, cex = 2)
title(main = 'Ponerinae diversity dynamics', cex = 2, line = 0)
# text(x = 0.1, y = 0.9, labels = "Net Div. Rates", font = 2)
addBAMMlegend(branch_div_rates_plot, location = "topleft", longFrac = 0.2, font = 2)

dev.off()


### 8.4/ Plot uncertainty in diversification rates ####

# Need to define the length of segments used to discretize time
# Compute continuous rates on branches from the regime parameters
# Extract summary stats on rates per branches for each time slice (mean, median, CI, HPD)
# Map those summary stats on phylo

?BAMMtools::dtRates # To compute discretized rates on small segments along branches

# Extract it from the plot.bammdata / branch_div_rates_plot?

# Median can be computed from BAMM_output_primates$eventData or using BAMMtools::dtRates to compute discretize rates



##### 9/ Plot diversification rates through time per clades #####

?branching.times # To extract divergence/branching times
?getRateThroughTimeMatrix  # To extract rates though time data
?plotRateThroughTime

# Load regime shifts metadata for MSC tree 
# MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/MSC_shifts_df.rds")
MSC_shifts_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/MSC_shifts_df.rds")

View(MSC_shifts_df)

# Load phylogeny with branch length scaled to mean net div rates
div_rates_tree <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/div_rates_tree.rds")

## Groups to focus on for rough phylogeny: 
 # Pachycondyla/Platythyrea: excluded: 1537, 2668
 # Ponera-group: included: 2668
 # Myopias/Leptogenys: included: 1763
 # Other taxa: included: 1537 + excluded: 1763 
 # Overall

## Groups to focus on for MCC phylogeny: 
# A: Platythyrea: excluded: 1536 (with root); included: 3025 (without root)
# Ponerini: included: 1536
# B: Background Ponerini: included 1536, excluded: 2819, 2393, 1540, 1749, 2126, 1762
# C: Diacamma: included: 2819 
# D: Afrotropical_Hypoponera: included: 2393
# E: Anochetus_Odontomachus (trap jaw ants): included: 1540
# F: Relict group: included: 1749
# G: Myopias: included: 2126
# H: Leptogenys subclade: included: 1762

max(div_rates_tree$phy$edge.length)

### 9.1/ Plot of rates through time with fuzzy intervals: median + all trajectories ####

# ## For rough phylogeny
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_LLT_per_clades_original_plot_fuzzy.pdf", height = 10, width = 12)
# 
# par(mfrow = c(2, 2))
# 
# root_age <- max(phytools::nodeHeights(phy))
# 
# # A/ All branches
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = seq(from = 0, to = 1, by = 0.01),
#                     intervalCol = "black", avgCol = "black",
#                     start.time = root_age,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 2)
# text(x = 100, y = 0.4, label = "All Ponerinae", col = "black", font = 4, cex = 2.0, pos = 4)
# 
# # B/ Only Pachycondyla/Platythyrea
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = seq(from = 0, to = 1, by = 0.01),
#                     intervalCol = "red", avgCol = "red", 
#                     start.time = root_age, 
#                     node = c(1537, 2668), nodetype = "exclude",
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 2)
# text(x = 100, y = 0.4, label = "Pachycondyla/Platythyrea", col = "red", font = 4, cex = 2.0, pos = 4)
# 
# # C/ Only Ponera
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = seq(from = 0, to = 1, by = 0.01), 
#                     intervalCol = "orange", avgCol = "orange",
#                     start.time = root_age, 
#                     node = 2668,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 1.5)
# text(x = 100, y = 0.4, label = "Ponera group", col = "orange", font = 4, cex = 2.0, pos = 4)
# 
# # D/ Only Myopias/Leptogenys
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = seq(from = 0, to = 1, by = 0.01), 
#                     intervalCol = "dodgerblue3", avgCol = "dodgerblue3",
#                     start.time = root_age, 
#                     node = 1763,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 1.5)
# text(x = 100, y = 0.4, label = "Myopias/Leptogenys", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)
# 
# dev.off()


## For MCC phylogeny
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_LLT_per_clades_original_plot_fuzzy.pdf", height = 18, width = 12)

par(mfrow = c(4, 2))

root_age <- max(phytools::nodeHeights(phy))

# A/ Only Platythyrea
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01),
                    intervalCol = "dodgerblue3", avgCol = "dodgerblue3", 
                    start.time = root_age, 
                    # node = c(1536), nodetype = "exclude", # With the root
                    node = c(3025), nodetype = "include", # Without the root
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "A/ Platythyrea", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)

# B/ All Ponerini (not just background)
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01),
                    intervalCol = "grey70", avgCol = "grey70",
                    start.time = root_age,
                    node = c(1536), nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "All Ponerini", col = "grey70", font = 4, cex = 2.0, pos = 4)

# C/ Diacamma
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "brown4", avgCol = "brown4",
                    start.time = root_age, 
                    node = 2819, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "C/ Diacamma", col = "brown4", font = 4, cex = 2.0, pos = 4)

# D/ Afrotropical Hypoponera
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "darkblue", avgCol = "darkblue",
                    start.time = root_age, 
                    node = 2393, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "D/ Afrotropical Hypoponera", col = "darkblue", font = 4, cex = 2.0, pos = 4)

# E/ Only trap jaw ants: Anochetus_Odontomachus: included: 1540
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "orange", avgCol = "orange",
                    start.time = root_age, 
                    node = 1540, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "E/ Trap jaw ants", col = "orange", font = 4, cex = 2.0, pos = 4)

# F/ Relict group: included: 1749
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "tan", avgCol = "tan",
                    start.time = root_age, 
                    node = 1749, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "F/ Relict group", col = "tan", font = 4, cex = 2.0, pos = 4)

# G/ Myopias
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "red", avgCol = "red",
                    start.time = root_age, 
                    node = 2126, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "G/ Myopias", col = "red", font = 4, cex = 2.0, pos = 4)

# H/ Leptogenys subclade: included: 1762
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "limegreen", avgCol = "limegreen",
                    start.time = root_age, 
                    node = 1762, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "H/ Leptogenys subclade", col = "limegreen", font = 4, cex = 2.0, pos = 4)

dev.off()


### 9.2/ Plot of rates through time with strict intervals: mean + 95% CI polygon ####

# Can be partitioned per clades (and per bioregions?)

# ## For rough phylogeny
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_LLT_per_clades_original_plot_CI95.pdf", height = 10, width = 12)
# 
# par(mfrow = c(2, 2))
# 
# root_age <- max(phytools::nodeHeights(phy))
# 
# # A/ All branches
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = c(0.025, 0.975),
#                     opacity = 0.2,
#                     intervalCol = "black", avgCol = "black",
#                     start.time = root_age,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 2)
# text(x = 100, y = 0.4, label = "All Ponerinae", col = "black", font = 4, cex = 2.0, pos = 4)
# 
# # B/ Only Pachycondyla/Platythyrea
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = c(0.025, 0.975),
#                     opacity = 0.2,
#                     intervalCol = "red", avgCol = "red", 
#                     start.time = root_age, 
#                     node = c(1537, 2668), nodetype = "exclude",
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 2)
# text(x = 100, y = 0.4, label = "Pachycondyla/Platythyrea", col = "red", font = 4, cex = 2.0, pos = 4)
# 
# # C/ Only Ponera
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = c(0.025, 0.975),
#                     opacity = 0.2, 
#                     intervalCol = "orange", avgCol = "orange",
#                     start.time = root_age, 
#                     node = 2668,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 1.5)
# text(x = 100, y = 0.4, label = "Ponera group", col = "orange", font = 4, cex = 2.0, pos = 4)
# 
# # D/ Only Myopias/Leptogenys
# plotRateThroughTime(BAMM_posterior_samples_data,
#                     intervals = c(0.025, 0.975),
#                     opacity = 0.2, 
#                     intervalCol = "dodgerblue3", avgCol = "dodgerblue3",
#                     start.time = root_age, 
#                     node = 1763,
#                     ylim = c(0, 0.5),
#                     xlim = c(root_age, 0),
#                     cex.axis = 1.5)
# text(x = 100, y = 0.4, label = "Myopias/Leptogenys", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)
# 
# dev.off()


## For MCC phylogeny
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_LLT_per_clades_original_plot_CI95.pdf", height = 18, width = 12)

par(mfrow = c(4, 2))

root_age <- max(phytools::nodeHeights(phy))

# A/ Only Platythyrea
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "dodgerblue3", avgCol = "dodgerblue3", 
                    start.time = root_age, 
                    # node = c(1536), nodetype = "exclude", # With the root
                    node = c(3025), nodetype = "include", # Without the root
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "A/ Platythyrea", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)

# B/ All Ponerini (not just background)
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "grey70", avgCol = "grey70",
                    start.time = root_age,
                    node = c(1536), nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "All Ponerini", col = "grey70", font = 4, cex = 2.0, pos = 4)

# C/ Diacamma
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "brown4", avgCol = "brown4",
                    start.time = root_age, 
                    node = 2819, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "C/ Diacamma", col = "brown4", font = 4, cex = 2.0, pos = 4)

# D/ Afrotropical Hypoponera
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "darkblue", avgCol = "darkblue",
                    start.time = root_age, 
                    node = 2393, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "D/ Afrotropical Hypoponera", col = "darkblue", font = 4, cex = 2.0, pos = 4)

# E/ Only trap jaw ants: Anochetus_Odontomachus: included: 1540
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "orange", avgCol = "orange",
                    start.time = root_age, 
                    node = 1540, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "E/ Trap jaw ants", col = "orange", font = 4, cex = 2.0, pos = 4)

# F/ Relict group: included: 1749
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "tan", avgCol = "tan",
                    start.time = root_age, 
                    node = 1749, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "F/ Relict group", col = "tan", font = 4, cex = 2.0, pos = 4)

# G/ Myopias
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "red", avgCol = "red",
                    start.time = root_age, 
                    node = 2126, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "G/ Myopias", col = "red", font = 4, cex = 2.0, pos = 4)

# H/ Leptogenys subclade: included: 1762
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "limegreen", avgCol = "limegreen",
                    start.time = root_age, 
                    node = 1762, nodetype = "include",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "H/ Leptogenys subclade", col = "limegreen", font = 4, cex = 2.0, pos = 4)

dev.off()


### 9.3/ Extract rates through time data ####

## Can set a time window and the size of the sliding window
# Use start.time, end.time, and nslices

find_closest_value <- function(x, y)
{ 
  ID <- which(abs(y - x) == min(abs(y - x)))
  target <- y[ID]
  return(target)
}

# ## Extract rates for rough phylogeny groups
# 
# # A/ Extract for the whole tree
# all_Ponerinae_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100)
# 
# all_Ponerinae_rates_through_time_matrix$net_div <- all_Ponerinae_rates_through_time_matrix$lambda - all_Ponerinae_rates_through_time_matrix$mu
# 
# str(all_Ponerinae_rates_through_time_matrix, max.level = 1)
# 
# dim(all_Ponerinae_rates_through_time_matrix$lambda) # Speciation rates per posterior sample (raw) x time window (columns)
# dim(all_Ponerinae_rates_through_time_matrix$mu) # Extinction rates per posterior sample (raw) x time window (columns)
# all_Ponerinae_rates_through_time_matrix$times # Mean time of the time windows. Names = age. Values = time since root.
# 
# all_Ponerinae_times <- round(as.numeric(names(all_Ponerinae_rates_through_time_matrix$times)), 1)
# 
# # B/ Extract for Only Pachycondyla/Platythyrea
# Pachycondyla_Platythyrea_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
#                                                                                node = c(1537, 2668), nodetype = "exclude")
# 
# Pachycondyla_Platythyrea_rates_through_time_matrix$net_div <- Pachycondyla_Platythyrea_rates_through_time_matrix$lambda - Pachycondyla_Platythyrea_rates_through_time_matrix$mu
# 
# 
# # C/ Extract for Only Ponera
# Ponera_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
#                                                                                node = 2668, nodetype = "include")
# 
# Ponera_rates_through_time_matrix$net_div <- Ponera_rates_through_time_matrix$lambda - Ponera_rates_through_time_matrix$mu
# 
# # Adjust to overall time scale
# Ponera_times <- round(as.numeric(names(Ponera_rates_through_time_matrix$times)), 1)
# Ponera_times_updated <- sapply(X = Ponera_times, FUN = find_closest_value, y = all_Ponerinae_times)
# 
# # Update matrix data
# Ponera_rates_through_time_matrix_updated <- Ponera_rates_through_time_matrix
# Ponera_update_indices <- match(all_Ponerinae_times, Ponera_times_updated)
# Ponera_rates_through_time_matrix_updated$lambda <- Ponera_rates_through_time_matrix$lambda[, Ponera_update_indices]
# Ponera_rates_through_time_matrix_updated$mu <- Ponera_rates_through_time_matrix$mu[, Ponera_update_indices]
# Ponera_rates_through_time_matrix_updated$times <- all_Ponerinae_times
# Ponera_rates_through_time_matrix_updated$net_div <- Ponera_rates_through_time_matrix$net_div[, Ponera_update_indices]
# 
# # D/ Extract for Only Myopias/Leptogenys
# Myopias_Leptogenys_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
#                                                              node = 1763, nodetype = "include")
# 
# Myopias_Leptogenys_rates_through_time_matrix$net_div <- Myopias_Leptogenys_rates_through_time_matrix$lambda - Myopias_Leptogenys_rates_through_time_matrix$mu
# 
# # Adjust to overall time scale
# Myopias_Leptogenys_times <- round(as.numeric(names(Myopias_Leptogenys_rates_through_time_matrix$times)), 1)
# Myopias_Leptogenys_times_updated <- sapply(X = Myopias_Leptogenys_times, FUN = find_closest_value, y = all_Ponerinae_times)
# 
# # Update matrix data
# Myopias_Leptogenys_rates_through_time_matrix_updated <- Myopias_Leptogenys_rates_through_time_matrix
# Myopias_Leptogenys_update_indices <- match(all_Ponerinae_times, Myopias_Leptogenys_times_updated)
# Myopias_Leptogenys_rates_through_time_matrix_updated$lambda <- Myopias_Leptogenys_rates_through_time_matrix$lambda[, Myopias_Leptogenys_update_indices]
# Myopias_Leptogenys_rates_through_time_matrix_updated$mu <- Myopias_Leptogenys_rates_through_time_matrix$mu[, Myopias_Leptogenys_update_indices]
# Myopias_Leptogenys_rates_through_time_matrix_updated$times <- all_Ponerinae_times
# Myopias_Leptogenys_rates_through_time_matrix_updated$net_div <- Myopias_Leptogenys_rates_through_time_matrix$net_div[, Myopias_Leptogenys_update_indices]


## Extract rates for MCC phylogeny groups

# 0/ Extract for the whole tree
all_Ponerinae_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100)

all_Ponerinae_rates_through_time_matrix$net_div <- all_Ponerinae_rates_through_time_matrix$lambda - all_Ponerinae_rates_through_time_matrix$mu

str(all_Ponerinae_rates_through_time_matrix, max.level = 1)

dim(all_Ponerinae_rates_through_time_matrix$lambda) # Speciation rates per posterior sample (raw) x time window (columns)
dim(all_Ponerinae_rates_through_time_matrix$mu) # Extinction rates per posterior sample (raw) x time window (columns)
all_Ponerinae_rates_through_time_matrix$times # Mean time of the time windows. Names = age. Values = time since root.

all_Ponerinae_times <- round(as.numeric(names(all_Ponerinae_rates_through_time_matrix$times)), 1)

# A/ Extract for Platythyrea
Platythyrea_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                               node = c(3025), nodetype = "include")

Platythyrea_rates_through_time_matrix$net_div <- Platythyrea_rates_through_time_matrix$lambda - Platythyrea_rates_through_time_matrix$mu

# Adjust to overall time scale
Platythyrea_times <- round(as.numeric(names(Platythyrea_rates_through_time_matrix$times)), 1)
Platythyrea_times_updated <- sapply(X = Platythyrea_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Platythyrea_rates_through_time_matrix_updated <- Platythyrea_rates_through_time_matrix
Platythyrea_update_indices <- match(all_Ponerinae_times, Platythyrea_times_updated)
Platythyrea_rates_through_time_matrix_updated$lambda <- Platythyrea_rates_through_time_matrix$lambda[, Platythyrea_update_indices]
Platythyrea_rates_through_time_matrix_updated$mu <- Platythyrea_rates_through_time_matrix$mu[, Platythyrea_update_indices]
Platythyrea_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Platythyrea_rates_through_time_matrix_updated$net_div <- Platythyrea_rates_through_time_matrix$net_div[, Platythyrea_update_indices]

# B/ Extract for all Ponerini
all_Ponerini_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                               node = c(1536), nodetype = "include")

all_Ponerini_rates_through_time_matrix$net_div <- all_Ponerini_rates_through_time_matrix$lambda - all_Ponerini_rates_through_time_matrix$mu

# Adjust to overall time scale
all_Ponerini_times <- round(as.numeric(names(all_Ponerini_rates_through_time_matrix$times)), 1)
all_Ponerini_times_updated <- sapply(X = all_Ponerini_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
all_Ponerini_rates_through_time_matrix_updated <- all_Ponerini_rates_through_time_matrix
all_Ponerini_update_indices <- match(all_Ponerinae_times, all_Ponerini_times_updated)
all_Ponerini_rates_through_time_matrix_updated$lambda <- all_Ponerini_rates_through_time_matrix$lambda[, all_Ponerini_update_indices]
all_Ponerini_rates_through_time_matrix_updated$mu <- all_Ponerini_rates_through_time_matrix$mu[, all_Ponerini_update_indices]
all_Ponerini_rates_through_time_matrix_updated$times <- all_Ponerinae_times
all_Ponerini_rates_through_time_matrix_updated$net_div <- all_Ponerini_rates_through_time_matrix$net_div[, all_Ponerini_update_indices]

# C/ Extract for Only Diacamma
Diacamma_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                               node = 2819, nodetype = "include")

Diacamma_rates_through_time_matrix$net_div <- Diacamma_rates_through_time_matrix$lambda - Diacamma_rates_through_time_matrix$mu

# Adjust to overall time scale
Diacamma_times <- round(as.numeric(names(Diacamma_rates_through_time_matrix$times)), 1)
Diacamma_times_updated <- sapply(X = Diacamma_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Diacamma_rates_through_time_matrix_updated <- Diacamma_rates_through_time_matrix
Diacamma_update_indices <- match(all_Ponerinae_times, Diacamma_times_updated)
Diacamma_rates_through_time_matrix_updated$lambda <- Diacamma_rates_through_time_matrix$lambda[, Diacamma_update_indices]
Diacamma_rates_through_time_matrix_updated$mu <- Diacamma_rates_through_time_matrix$mu[, Diacamma_update_indices]
Diacamma_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Diacamma_rates_through_time_matrix_updated$net_div <- Diacamma_rates_through_time_matrix$net_div[, Diacamma_update_indices]

# D/ Extract for Afrotropical_Hypoponera
Afrotropical_Hypoponera_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                             node = 2393, nodetype = "include")

Afrotropical_Hypoponera_rates_through_time_matrix$net_div <- Afrotropical_Hypoponera_rates_through_time_matrix$lambda - Afrotropical_Hypoponera_rates_through_time_matrix$mu

# Adjust to overall time scale
Afrotropical_Hypoponera_times <- round(as.numeric(names(Afrotropical_Hypoponera_rates_through_time_matrix$times)), 1)
Afrotropical_Hypoponera_times_updated <- sapply(X = Afrotropical_Hypoponera_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Afrotropical_Hypoponera_rates_through_time_matrix_updated <- Afrotropical_Hypoponera_rates_through_time_matrix
Afrotropical_Hypoponera_update_indices <- match(all_Ponerinae_times, Afrotropical_Hypoponera_times_updated)
Afrotropical_Hypoponera_rates_through_time_matrix_updated$lambda <- Afrotropical_Hypoponera_rates_through_time_matrix$lambda[, Afrotropical_Hypoponera_update_indices]
Afrotropical_Hypoponera_rates_through_time_matrix_updated$mu <- Afrotropical_Hypoponera_rates_through_time_matrix$mu[, Afrotropical_Hypoponera_update_indices]
Afrotropical_Hypoponera_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Afrotropical_Hypoponera_rates_through_time_matrix_updated$net_div <- Afrotropical_Hypoponera_rates_through_time_matrix$net_div[, Afrotropical_Hypoponera_update_indices]

# E/ Extract for Only Anochetus_Odontomachus (trap jaw ants)
Anochetus_Odontomachus_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                             node = 1540, nodetype = "include")

Anochetus_Odontomachus_rates_through_time_matrix$net_div <- Anochetus_Odontomachus_rates_through_time_matrix$lambda - Anochetus_Odontomachus_rates_through_time_matrix$mu

# Adjust to overall time scale
Anochetus_Odontomachus_times <- round(as.numeric(names(Anochetus_Odontomachus_rates_through_time_matrix$times)), 1)
Anochetus_Odontomachus_times_updated <- sapply(X = Anochetus_Odontomachus_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Anochetus_Odontomachus_rates_through_time_matrix_updated <- Anochetus_Odontomachus_rates_through_time_matrix
Anochetus_Odontomachus_update_indices <- match(all_Ponerinae_times, Anochetus_Odontomachus_times_updated)
Anochetus_Odontomachus_rates_through_time_matrix_updated$lambda <- Anochetus_Odontomachus_rates_through_time_matrix$lambda[, Anochetus_Odontomachus_update_indices]
Anochetus_Odontomachus_rates_through_time_matrix_updated$mu <- Anochetus_Odontomachus_rates_through_time_matrix$mu[, Anochetus_Odontomachus_update_indices]
Anochetus_Odontomachus_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Anochetus_Odontomachus_rates_through_time_matrix_updated$net_div <- Anochetus_Odontomachus_rates_through_time_matrix$net_div[, Anochetus_Odontomachus_update_indices]

# F/ Extract for Relict group
Relict_group_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                             node = 1749, nodetype = "include")

Relict_group_rates_through_time_matrix$net_div <- Relict_group_rates_through_time_matrix$lambda - Relict_group_rates_through_time_matrix$mu

# Adjust to overall time scale
Relict_group_times <- round(as.numeric(names(Relict_group_rates_through_time_matrix$times)), 1)
Relict_group_times_updated <- sapply(X = Relict_group_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Relict_group_rates_through_time_matrix_updated <- Relict_group_rates_through_time_matrix
Relict_group_update_indices <- match(all_Ponerinae_times, Relict_group_times_updated)
Relict_group_rates_through_time_matrix_updated$lambda <- Relict_group_rates_through_time_matrix$lambda[, Relict_group_update_indices]
Relict_group_rates_through_time_matrix_updated$mu <- Relict_group_rates_through_time_matrix$mu[, Relict_group_update_indices]
Relict_group_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Relict_group_rates_through_time_matrix_updated$net_div <- Relict_group_rates_through_time_matrix$net_div[, Relict_group_update_indices]

# G/ Extract for Only Myopias
Myopias_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                         node = 2126, nodetype = "include")

Myopias_rates_through_time_matrix$net_div <- Myopias_rates_through_time_matrix$lambda - Myopias_rates_through_time_matrix$mu

# Adjust to overall time scale
Myopias_times <- round(as.numeric(names(Myopias_rates_through_time_matrix$times)), 1)
Myopias_times_updated <- sapply(X = Myopias_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Myopias_rates_through_time_matrix_updated <- Myopias_rates_through_time_matrix
Myopias_update_indices <- match(all_Ponerinae_times, Myopias_times_updated)
Myopias_rates_through_time_matrix_updated$lambda <- Myopias_rates_through_time_matrix$lambda[, Myopias_update_indices]
Myopias_rates_through_time_matrix_updated$mu <- Myopias_rates_through_time_matrix$mu[, Myopias_update_indices]
Myopias_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Myopias_rates_through_time_matrix_updated$net_div <- Myopias_rates_through_time_matrix$net_div[, Myopias_update_indices]

# H: Leptogenys subclade
Leptogenys_subclade_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                             node = 1762, nodetype = "include")

Leptogenys_subclade_rates_through_time_matrix$net_div <- Leptogenys_subclade_rates_through_time_matrix$lambda - Leptogenys_subclade_rates_through_time_matrix$mu

# Adjust to overall time scale
Leptogenys_subclade_times <- round(as.numeric(names(Leptogenys_subclade_rates_through_time_matrix$times)), 1)
Leptogenys_subclade_times_updated <- sapply(X = Leptogenys_subclade_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Leptogenys_subclade_rates_through_time_matrix_updated <- Leptogenys_subclade_rates_through_time_matrix
Leptogenys_subclade_update_indices <- match(all_Ponerinae_times, Leptogenys_subclade_times_updated)
Leptogenys_subclade_rates_through_time_matrix_updated$lambda <- Leptogenys_subclade_rates_through_time_matrix$lambda[, Leptogenys_subclade_update_indices]
Leptogenys_subclade_rates_through_time_matrix_updated$mu <- Leptogenys_subclade_rates_through_time_matrix$mu[, Leptogenys_subclade_update_indices]
Leptogenys_subclade_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Leptogenys_subclade_rates_through_time_matrix_updated$net_div <- Leptogenys_subclade_rates_through_time_matrix$net_div[, Leptogenys_subclade_update_indices]


View(Leptogenys_subclade_rates_through_time_matrix_updated)


# ## Aggregate LTT data in an array for rough phylogeny
# LTT_per_clades_array <- array(data = NA,
#                               dim = c(dim(all_Ponerinae_rates_through_time_matrix$net_div), 4),
#                               dimnames = list(sample = paste0("sample_",1:dim(all_Ponerinae_rates_through_time_matrix$net_div)[1]),
#                                               time = names(all_Ponerinae_rates_through_time_matrix$times),
#                                               group = c("All_Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys")))
# 
# LTT_per_clades_array[,,1] <- all_Ponerinae_rates_through_time_matrix$net_div                           
# LTT_per_clades_array[,,2] <- Pachycondyla_Platythyrea_rates_through_time_matrix$net_div
# LTT_per_clades_array[,,3] <- Ponera_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
# LTT_per_clades_array[,,4] <- Myopias_Leptogenys_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale



## Aggregate LTT data in an array for MCC phylogeny
LTT_per_clades_array <- array(data = NA,
                              dim = c(dim(all_Ponerinae_rates_through_time_matrix$net_div), 10, 2),
                              dimnames = list(sample = paste0("sample_",1:dim(all_Ponerinae_rates_through_time_matrix$net_div)[1]),
                                              time = names(all_Ponerinae_rates_through_time_matrix$times),
                                              group = c("All_Ponerinae", "Platythyrea", "All_Ponerini", "Background_Ponerini", "Diacamma", "Afrotropical_Hypoponera", "Anochetus_Odontomachus", "Relict_group", "Myopias", "Leptogenys_subclade"),
                                              data_type = c("net_div_rates", "richness")))

LTT_per_clades_array[,,1,"net_div_rates"] <- all_Ponerinae_rates_through_time_matrix$net_div                           
LTT_per_clades_array[,,2,"net_div_rates"] <- Platythyrea_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,3,"net_div_rates"] <- all_Ponerini_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,5,"net_div_rates"] <- Diacamma_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,6,"net_div_rates"] <- Afrotropical_Hypoponera_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,7,"net_div_rates"] <- Anochetus_Odontomachus_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,8,"net_div_rates"] <- Relict_group_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,9,"net_div_rates"] <- Myopias_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,10,"net_div_rates"] <- Leptogenys_subclade_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale

### 9.4/ Compute Background Ponerini rates as the remaining rates once other focal subclades have been subtracted from the mean rates ####

## 9.4.1/ Get LTT data for each group ###

# Function to detect the ID of the next higher values in vector Y compared to value X
# Used to match sliding window times with LTT ages
find_next_higher_value_ID <- function(x, y)
{ 
  # Detect higher values
  higher_values_ID <- which((y - x) > 0)
  higher_values <- y[higher_values_ID]
  # Find the closest one
  target_ID <- higher_values_ID[which.min(higher_values)]
  return(target_ID)
}

root_age <- max(phytools::nodeHeights(phy))
all_edges_ages_df <- round(root_age - phytools::nodeHeights(phy), 5)

# 0/ For all Ponerinae
LTT_all_Ponerinae <- phytools::ltt(tree = phy) # Get LTT data
# Convert to df
LTT_all_Ponerinae_df <- data.frame(richness = LTT_all_Ponerinae$ltt, time = LTT_all_Ponerinae$times, age = root_age - LTT_all_Ponerinae$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_all_Ponerinae_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_all_Ponerinae_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_all_Ponerinae_df_updated$richness <- LTT_all_Ponerinae_df$richness[matching_ID]
# Include stem edge
LTT_all_Ponerinae_df_updated$richness[LTT_all_Ponerinae_df_updated$time > root_age] <- 1

# A/ For Platythyrea
phy_Platythyrea <- ape::extract.clade(phy = phy,
                                       node = 3025) # Without the root
root_age_Platythyrea <- max(phytools::nodeHeights(phy_Platythyrea))
LTT_Platythyrea <- phytools::ltt(tree = phy_Platythyrea) # Get LTT data
# Convert to df
LTT_Platythyrea_df <- data.frame(richness = LTT_Platythyrea$ltt, time = LTT_Platythyrea$times, age = root_age_Platythyrea - LTT_Platythyrea$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Platythyrea_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Platythyrea_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Platythyrea_df_updated$richness <- LTT_Platythyrea_df$richness[matching_ID]
# Include stem edge
LTT_Platythyrea_df_updated$richness[LTT_Platythyrea_df_updated$time > root_age_Platythyrea & LTT_Platythyrea_df_updated$time < root_age] <- 1

# For all Ponerini
phy_all_Ponerini <- ape::extract.clade(phy = phy,
                                       node = 1536)
root_age_all_Ponerini <- max(phytools::nodeHeights(phy_all_Ponerini))
LTT_all_Ponerini <- phytools::ltt(tree = phy_all_Ponerini) # Get LTT data
# Convert to df
LTT_all_Ponerini_df <- data.frame(richness = LTT_all_Ponerini$ltt, time = LTT_all_Ponerini$times, age = root_age_all_Ponerini - LTT_all_Ponerini$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_all_Ponerini_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_all_Ponerini_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_all_Ponerini_df_updated$richness <- LTT_all_Ponerini_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Poneriini_background"]
LTT_all_Ponerini_df_updated$richness[(LTT_all_Ponerini_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_all_Ponerini_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# C/ For Diacamma
phy_Diacamma <- ape::extract.clade(phy = phy,
                                       node = 2819)
root_age_Diacamma <- max(phytools::nodeHeights(phy_Diacamma))
LTT_Diacamma <- phytools::ltt(tree = phy_Diacamma) # Get LTT data
# Convert to df
LTT_Diacamma_df <- data.frame(richness = LTT_Diacamma$ltt, time = LTT_Diacamma$times, age = root_age_Diacamma - LTT_Diacamma$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Diacamma_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Diacamma_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Diacamma_df_updated$richness <- LTT_Diacamma_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Diacamma"]
LTT_Diacamma_df_updated$richness[(LTT_Diacamma_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Diacamma_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# D/ For Afrotropical_Hypoponera
phy_Afrotropical_Hypoponera <- ape::extract.clade(phy = phy,
                                   node = 2393)
root_age_Afrotropical_Hypoponera <- max(phytools::nodeHeights(phy_Afrotropical_Hypoponera))
LTT_Afrotropical_Hypoponera <- phytools::ltt(tree = phy_Afrotropical_Hypoponera) # Get LTT data
# Convert to df
LTT_Afrotropical_Hypoponera_df <- data.frame(richness = LTT_Afrotropical_Hypoponera$ltt, time = LTT_Afrotropical_Hypoponera$times, age = root_age_Afrotropical_Hypoponera - LTT_Afrotropical_Hypoponera$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Afrotropical_Hypoponera_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Afrotropical_Hypoponera_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Afrotropical_Hypoponera_df_updated$richness <- LTT_Afrotropical_Hypoponera_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Hypoponera_Afrotropics"]
LTT_Afrotropical_Hypoponera_df_updated$richness[(LTT_Afrotropical_Hypoponera_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Afrotropical_Hypoponera_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# E/ For Anochetus_Odontomachus (trap jaw ants)
phy_Anochetus_Odontomachus <- ape::extract.clade(phy = phy,
                                                  node = 1540)
root_age_Anochetus_Odontomachus <- max(phytools::nodeHeights(phy_Anochetus_Odontomachus))
LTT_Anochetus_Odontomachus <- phytools::ltt(tree = phy_Anochetus_Odontomachus) # Get LTT data
# Convert to df
LTT_Anochetus_Odontomachus_df <- data.frame(richness = LTT_Anochetus_Odontomachus$ltt, time = LTT_Anochetus_Odontomachus$times, age = root_age_Anochetus_Odontomachus - LTT_Anochetus_Odontomachus$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Anochetus_Odontomachus_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Anochetus_Odontomachus_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Anochetus_Odontomachus_df_updated$richness <- LTT_Anochetus_Odontomachus_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Trap_jaw_ants"]
LTT_Anochetus_Odontomachus_df_updated$richness[(LTT_Anochetus_Odontomachus_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Anochetus_Odontomachus_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# F/ For Relict_group
phy_Relict_group <- ape::extract.clade(phy = phy,
                                       node = 1749)
root_age_Relict_group <- max(phytools::nodeHeights(phy_Relict_group))
LTT_Relict_group <- phytools::ltt(tree = phy_Relict_group) # Get LTT data
# Convert to df
LTT_Relict_group_df <- data.frame(richness = LTT_Relict_group$ltt, time = LTT_Relict_group$times, age = root_age_Relict_group - LTT_Relict_group$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Relict_group_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Relict_group_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Relict_group_df_updated$richness <- LTT_Relict_group_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Relict_group"]
LTT_Relict_group_df_updated$richness[(LTT_Relict_group_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Relict_group_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# G/ For Myopias
phy_Myopias <- ape::extract.clade(phy = phy,
                                  node = 2126)
root_age_Myopias <- max(phytools::nodeHeights(phy_Myopias))
LTT_Myopias <- phytools::ltt(tree = phy_Myopias) # Get LTT data
# Convert to df
LTT_Myopias_df <- data.frame(richness = LTT_Myopias$ltt, time = LTT_Myopias$times, age = root_age_Myopias - LTT_Myopias$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Myopias_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Myopias_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Myopias_df_updated$richness <- LTT_Myopias_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Myopias"]
LTT_Myopias_df_updated$richness[(LTT_Myopias_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Myopias_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# H/ For Leptogenys_subclade
phy_Leptogenys_subclade <- ape::extract.clade(phy = phy,
                                  node = 1762)
root_age_Leptogenys_subclade <- max(phytools::nodeHeights(phy_Leptogenys_subclade))
LTT_Leptogenys_subclade <- phytools::ltt(tree = phy_Leptogenys_subclade) # Get LTT data
# Convert to df
LTT_Leptogenys_subclade_df <- data.frame(richness = LTT_Leptogenys_subclade$ltt, time = LTT_Leptogenys_subclade$times, age = root_age_Leptogenys_subclade - LTT_Leptogenys_subclade$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_Leptogenys_subclade_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_Leptogenys_subclade_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_Leptogenys_subclade_df_updated$richness <- LTT_Leptogenys_subclade_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Leptogenys_subgroup"]
LTT_Leptogenys_subclade_df_updated$richness[(LTT_Leptogenys_subclade_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_Leptogenys_subclade_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

# D/ For background Ponerini

## Need the other trees to drop tips from the subclades and compute LTT
tips_to_drop <- c(phy_Diacamma$tip.label, phy_Afrotropical_Hypoponera$tip.label, phy_Anochetus_Odontomachus$tip.label,
                  phy_Relict_group$tip.label, phy_Myopias$tip.label, phy_Leptogenys_subclade$tip.label)
phy_background_Ponerini <- ape::extract.clade(phy = phy,
                                              node = 1536)
phy_background_Ponerini <- ape::drop.tip(phy = phy_background_Ponerini, tip = tips_to_drop)
root_age_background_Ponerini <- max(phytools::nodeHeights(phy_background_Ponerini))
LTT_background_Ponerini <- phytools::ltt(tree = phy_background_Ponerini) # Get LTT data
# Convert to df
LTT_background_Ponerini_df <- data.frame(richness = LTT_background_Ponerini$ltt, time = LTT_background_Ponerini$times, age = root_age_background_Ponerini - LTT_background_Ponerini$times)
# Assign updated richness accroding to standard time scale by matching times
LTT_background_Ponerini_df_updated <- data.frame(richness = NA, time = all_Ponerinae_times)
matching_ID <- sapply(X = all_Ponerinae_times, FUN = find_next_higher_value_ID, y = LTT_background_Ponerini_df$age)
matching_ID <- unlist(lapply(X = matching_ID, FUN = function (x) { if (purrr::is_empty(x)) { y <- NA } else { y <- x } } ))
LTT_background_Ponerini_df_updated$richness <- LTT_background_Ponerini_df$richness[matching_ID]
# Include stem edge
stem_edge_ID <- MSC_shifts_df$edge_ID[MSC_shifts_df$regime_group_name == "Poneriini_background"]
LTT_background_Ponerini_df_updated$richness[(LTT_background_Ponerini_df_updated$time < all_edges_ages_df[stem_edge_ID,1]) & (LTT_background_Ponerini_df_updated$time > all_edges_ages_df[stem_edge_ID,2])] <- 1

## Inform the LTT array
LTT_per_clades_array[,,"All_Ponerinae","richness"] <- matrix(data = LTT_all_Ponerinae_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Platythyrea","richness"] <- matrix(data = LTT_Platythyrea_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"All_Ponerini","richness"] <- matrix(data = LTT_all_Ponerini_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Background_Ponerini","richness"] <- matrix(data = LTT_background_Ponerini_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Diacamma","richness"] <- matrix(data = LTT_Diacamma_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Afrotropical_Hypoponera","richness"] <- matrix(data = LTT_Afrotropical_Hypoponera_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Anochetus_Odontomachus","richness"] <- matrix(data = LTT_Anochetus_Odontomachus_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Relict_group", "richness"] <- matrix(data = LTT_Relict_group_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Myopias","richness"] <- matrix(data = LTT_Myopias_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)
LTT_per_clades_array[,,"Leptogenys_subclade","richness"] <- matrix(data = LTT_Leptogenys_subclade_df_updated$richness, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2], byrow = T)

## Inform MSC shift metadata df with current richness and crown ages
all_regime_groups_phy_list <- list(Platythyrea = phy_Platythyrea, Poneriini_background = phy_background_Ponerini, Trap_jaw_ants = phy_Anochetus_Odontomachus, Relict_group = phy_Relict_group,
                            Leptogenys_subgroup = phy_Leptogenys_subclade, Hypoponera_Afrotropics = phy_Afrotropical_Hypoponera, Diacamma = phy_Diacamma, Myopias = phy_Myopias)

MSC_shifts_df$current_richness <- lapply(X = all_regime_groups_phy_list, FUN = function (x) { length(x$tip.label) })
MSC_shifts_df$crown_age <- lapply(X = all_regime_groups_phy_list, FUN = function (x) { max(phytools::nodeHeights(x)) })

## 9.4.2/ Compute background Ponerini rates

compute_background_rate <- function(rate_overall, rate_subclade, N_overall, N_subclade)
{
  rate_background <- ((rate_overall * N_overall) - (rate_subclade * N_subclade)) / (N_overall - N_subclade)
  return(rate_background)
}

Subclades_names <- c("Diacamma", "Afrotropical_Hypoponera", "Anochetus_Odontomachus", "Relict_group", "Myopias", "Leptogenys_subclade")
Subclades_rates <- matrix(data = NA, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2])

# Compute weighted mean rates of all other subclades
for (i in 1:dim(LTT_per_clades_array)[1])
{
  for (j in 1:dim(LTT_per_clades_array)[2])
  {
    rates_data <- LTT_per_clades_array[i,j,Subclades_names,"net_div_rates"]
    rates_data_no_NA <- replace_na(data = rates_data, replace = 0)
    richness_data <- LTT_per_clades_array[i,j,Subclades_names,"richness"]
    richness_data_no_NA <- replace_na(data = richness_data, replace = 0)
    Subclades_rates[i,j] <- weighted.mean(x = rates_data_no_NA, w = richness_data_no_NA)
  }
}
row.names(Subclades_rates) <- dimnames(LTT_per_clades_array)$sample
colnames(Subclades_rates) <- dimnames(LTT_per_clades_array)$time

# Compute richness of all other subclades
Subclades_richness <- apply(X = LTT_per_clades_array[,,Subclades_names,"richness"], MARGIN = c(1,2), FUN = sum, na.rm = T)

Background_Ponerini_rates <- matrix(data = NA, nrow = dim(LTT_per_clades_array)[1], ncol = dim(LTT_per_clades_array)[2])
# Compute background rates for Ponerini
for (i in 1:dim(LTT_per_clades_array)[1])
{
  for (j in 1:dim(LTT_per_clades_array)[2])
  {
    rate_overall_data <- LTT_per_clades_array[i,j,"All_Ponerini","net_div_rates"]
    rate_subclades_data <- Subclades_rates[i,j]
    rate_subclades_data_no_NA <- replace_na(data = rate_subclades_data, replace = 0)
    N_overall <- LTT_per_clades_array[i,j,"All_Ponerini","richness"]
    N_subclades_data <- Subclades_richness[i,j]
    N_subclades_data_no_NA <- replace_na(data = N_subclades_data, replace = 0)
    Background_Ponerini_rates[i,j] <- compute_background_rate(rate_overall = rate_overall_data,
                                                    rate_subclade = rate_subclades_data_no_NA,
                                                    N_overall = N_overall,
                                                    N_subclade = N_subclades_data_no_NA)
  }
}
row.names(Background_Ponerini_rates) <- dimnames(LTT_per_clades_array)$sample
colnames(Background_Ponerini_rates) <- dimnames(LTT_per_clades_array)$time

# Inform LTT array
LTT_per_clades_array[,,"Background_Ponerini","net_div_rates"] <- Background_Ponerini_rates

## Save LTT array of richness/RTT per clades
# saveRDS(LTT_per_clades_array, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/LTT_per_clades_array.rds")
saveRDS(LTT_per_clades_array, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/LTT_per_clades_array.rds")


## To include stem rates, need to use the exclude option to exclude the whole sister-clade of the focal clade...


### 9.5/ Melt array of RTT per clades in df ####

# Aggregate LTT data in a unique melted data.frame
LTT_per_clades_melted_df <- reshape2::melt(LTT_per_clades_array)
LTT_per_clades_melted_df$time <- round(as.numeric(LTT_per_clades_melted_df$time), 1)

# Remove richness data
LTT_per_clades_melted_df <- LTT_per_clades_melted_df %>% 
  filter(data_type != "richness") %>%
  dplyr::select(-data_type)
  
# Rename rates data
LTT_per_clades_melted_df <- LTT_per_clades_melted_df %>% 
  rename(net_div_rates = value)

# Order groups
# LTT_per_clades_melted_df$group <- factor(LTT_per_clades_melted_df$group, levels = c("All_Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys"), labels = c("All Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys"))
# LTT_per_clades_melted_df$group <- factor(LTT_per_clades_melted_df$group, levels = c("All_Ponerinae", "Platythyrea", "Ponerini", "Anochetus_Odontomachus", "Myopias", "Diacamma"), labels = c("All Ponerinae", "Platythyrea", "All Ponerini", "Anochetus/Odontomachus", "Myopias", "Diacamma"))
LTT_per_clades_melted_df$group <- factor(LTT_per_clades_melted_df$group, levels = c("All_Ponerinae", "Platythyrea", "All_Ponerini", "Background_Ponerini", "Diacamma", "Afrotropical_Hypoponera", "Anochetus_Odontomachus", "Relict_group", "Myopias", "Leptogenys_subclade"),
                                         labels = c("All Ponerinae", "A/ Platythyrea", "All Ponerini", "B/ Background Ponerini", "C/ Diacamma", "D/ Afrotropical Hypoponera", "E/ Trap-jaw ants", "F/ Relict group", "G/ Myopias", "H/ Leptogenys subclade"))

## Save LTT data per clades
# saveRDS(LTT_per_clades_melted_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/LTT_per_clades_melted_df.rds")
saveRDS(LTT_per_clades_melted_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/LTT_per_clades_melted_df.rds")


### 9.6/ Plot RTT per clades: ggplot version ####

## Load RTT data per clades
# LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/LTT_per_clades_melted_df.rds")
LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/LTT_per_clades_melted_df.rds")

# Remove All Ponerinae (redundant with All Ponerini)
# Remove All Ponerini (redundant with Background Ponerini)
LTT_per_clades_melted_df <- LTT_per_clades_melted_df %>% 
  filter(group != "All Ponerinae") %>% 
  filter(group != "All Ponerini")
  # filter(group != "Diacamma")

## Compute medians
LTT_per_clades_melted_df_median <- LTT_per_clades_melted_df %>% 
  group_by(time, group) %>% 
  summarize(median_rates = median(net_div_rates)) %>%
  ungroup()

## Compute quantiles

LTT_per_clades_melted_df_quantiles <- LTT_per_clades_melted_df %>% 
  filter(time <= 85)

LTT_per_clades_melted_df_quantiles <- LTT_per_clades_melted_df_quantiles %>% 
  group_by(time, group) %>% 
  # Compute quantiles
  reframe(quant_rates = quantile(net_div_rates, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
  group_by(time, group) %>%
  mutate(quantile = seq(from = 0, to = 1, by = 0.01),
         quantile_ID = c(1:50, 51, 50:1)) %>%
  filter(quantile_ID != 51) %>% # Remove median
  # Assign points ID (order for drawing the polygon)
  group_by(group, quantile_ID) %>%
  arrange(group, quantile_ID, quantile) %>%
  # mutate(points_ID = c(1:100, 200:101)) %>%
  mutate(n_points = dplyr::n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  dplyr::select(-n_points) %>%
  # Reorder by points ID
  arrange(group, quantile_ID, points_ID) %>%
  filter(!is.na(quant_rates)) %>%
  mutate(points_ID = row_number()) %>%
  ungroup()

## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_LLT_per_clades_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_LLT_per_clades_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

LTT_per_clades_ggplot_lines <- ggplot(data = LTT_per_clades_melted_df) +
  
  # Plot 1000 line replicates
  geom_line(data = LTT_per_clades_melted_df,
            mapping = aes(y = net_div_rates, x = time, group = interaction(group, sample), col = group),
            alpha = 0.01,
            linewidth = 5.0) +

  # Plot mean lines
  geom_line(data = LTT_per_clades_melted_df_median,
           mapping = aes(y = median_rates, x = time, group = group, col = group),
           alpha = 1.0,
           linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  # scale_color_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) + # For rough phylogeny
  # scale_color_manual("Groups", values = c("dodgerblue3", "grey70", "orange", "red")) +
  scale_color_manual("Regimes", values = c("dodgerblue3", "grey70", "brown4", "darkblue", "orange", "tan", "red", "limegreen")) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(LTT_per_clades_ggplot_lines)

dev.off()


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_LLT_per_clades_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_LLT_per_clades_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

LTT_per_clades_ggplot_quantiles <- ggplot(data = LTT_per_clades_melted_df_quantiles) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = LTT_per_clades_melted_df_quantiles,
               mapping = aes(y = quant_rates, x = time, group = interaction(group, quantile_ID), fill = group),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = LTT_per_clades_melted_df_median,
            mapping = aes(y = median_rates, x = time, group = group, col = group),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  # scale_fill_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) +
  # scale_fill_manual("Groups", values = c("dodgerblue3", "grey70", "orange", "red")) +
  scale_fill_manual("Regimes", values = c("dodgerblue3", "grey70", "brown4", "darkblue", "orange", "tan", "red", "limegreen")) +
  
  # Adjust color scheme and legend
  # scale_color_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) + # For rough phylogeny
  # scale_color_manual("Groups", values = c("dodgerblue3", "grey70", "orange", "red")) +
  scale_color_manual("Regimes", values = c("dodgerblue3", "grey70", "brown4", "darkblue", "orange", "tan", "red", "limegreen")) +

  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(LTT_per_clades_ggplot_quantiles)

dev.off()



##### 10/ Plot diversification rates through time per bioregions #####

# Need to get data for each posterior sample to account for uncertainty in diversification rates
# Need to extract rates per branch to be able to weight rates per bioregion
# Easy version = use mean rates and membership weight per branches
# Tricky version = get rates and membership per time segment as in densityMap
# Then, use membership probabilities as weights to compute weighted mean of rate per bioregions

# May also remove rates with a minimum number of lineages (sum of proba/weights)

### 10.1/ Extract rates per branches across all posterior samples ####

BAMM_div_rates_trees <- list()
for (i in seq_along(BAMM_posterior_samples_data$eventData))
{
  # i <- 1
  
  # Subset configuration i
  BAMM_tree_i <- subsetEventData(BAMM_posterior_samples_data, index = i)
  
  # Obtain tree scaled by net diversification rates
  BAMM_div_rates_tree_i <- getMeanBranchLengthTree(BAMM_tree_i,
                                                   rate = "ndr") # Net Diversification Rates
  
  # plot(BAMM_div_rates_tree_i$phy)
  
  # Store scaled tree
  BAMM_div_rates_trees[[i]] <- BAMM_div_rates_tree_i$phy
  
  ## Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - BAMM tree scaled by div rates created for n°", i, "/", length(BAMM_posterior_samples_data$eventData),"\n"))
  }
}

# Save BAMM trees scaled by div rates
# saveRDS(BAMM_div_rates_trees, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_div_rates_trees.rds")
saveRDS(BAMM_div_rates_trees, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_div_rates_trees.rds")

### 10.2/ Compute weighted mean rates along phylogenies ####

# Easy version = use mean rates and membership weight per branches
# Tricky version = get rates and membership per time segment as in densityMap

## Load BAMM trees scaled by div rates
# BAMM_div_rates_trees <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_div_rates_trees.rds")
BAMM_div_rates_trees <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_div_rates_trees.rds")

## Load treedata with bioregion membership information
# Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata_for_ARE.rds")
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata_for_ARE.rds")

# Mean state = PP of bioregion membership
Ponerinae_phylogeny_1534t_treedata@extraInfo 

## Aggregate bioregion membership for OW vs. NW
Ponerinae_phylogeny_1534t_treedata@extraInfo$mean_state_OW <- apply(X = Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Afrotropics", "mean_state_Eastern Palearctic", "mean_state_Indomalaya",
                                                                                                               "mean_state_Western Palearctic", "mean_state_Australasia")],
                                                                    MARGIN = 1, FUN = sum, na.rm = T)
Ponerinae_phylogeny_1534t_treedata@extraInfo$mean_state_NW <- apply(X = Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Neotropics", "mean_state_Nearctic")],
                                                                    MARGIN = 1, FUN = sum, na.rm = T)

## Aggregate bioregion membership for Tropics vs. Temperate
Ponerinae_phylogeny_1534t_treedata@extraInfo$mean_state_Tropics <- apply(X = Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Afrotropics", "mean_state_Australasia",
                                                                                                                              "mean_state_Indomalaya", "mean_state_Neotropics")],
                                                                         MARGIN = 1, FUN = sum, na.rm = T)
Ponerinae_phylogeny_1534t_treedata@extraInfo$mean_state_Temperate <- apply(X = Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Eastern Palearctic", "mean_state_Western Palearctic", "mean_state_Nearctic")],
                                                                           MARGIN = 1, FUN = sum, na.rm = T)


## Function to compute weighted mean values along phylogenies 
compute_weighted_mean_values_along_phylo <- function(treedata_obj, values, weight_groups, time_scale, keep_weight_sums = T, melt = F)
{
  # Initiate final array
  if(keep_weight_sums)
  { 
    value_types <- c("weighted_mean", "weight_sum")
  } else {
    value_types <-"weighted_mean"
  }
  weighted_mean_per_groups_along_phylo_array <- array(data = NA,
                                                      dim = c(length(time_scale), length(weight_groups), length(value_types)),
                                                      dimnames = list(time = time_scale, group = weight_groups, value_type = value_types)
                                                      )
  
  # Get node ages per branch (no root edge)
  all_node_heights_per_edge <- phytools::nodeHeights(treedata_obj@phylo)
  root_age <- max(phytools::nodeHeights(treedata_obj@phylo)[,2])
  all_node_ages_per_edge <- as.data.frame(round(root_age - all_node_heights_per_edge, 2))
  names(all_node_ages_per_edge) <- c("rootward_node_age", "tipward_node_age")
  
  # Detect root node ID
  root_node_ID <- length(treedata_obj@phylo$tip.label) + 1
  
  # Loop per time points
  for (i in seq_along(time_scale))
  {
    # i <- 1
    # i <- 5
    
    # Extract time point
    time_i <- time_scale[i]
    
    # Initiate node age df 
    all_node_ages_per_edge_i <- all_node_ages_per_edge
    
    # Identify tipward nodes of branch present at time i
    if (time_i == round(root_age,1))
    {
      selected_tipward_nodes <- root_node_ID
    } else {
      all_node_ages_per_edge_i$rootward_test <- all_node_ages_per_edge_i$rootward_node_age > time_i
      all_node_ages_per_edge_i$tipward_test <- all_node_ages_per_edge_i$tipward_node_age <= time_i
      all_node_ages_per_edge_i$present <- all_node_ages_per_edge_i$rootward_test & all_node_ages_per_edge_i$tipward_test
      
      selected_tipward_nodes <- treedata_obj@phylo$edge[all_node_ages_per_edge_i$present,2]
    }
    
    # Extract values data
    values_data_i <- treedata_obj@extraInfo[selected_tipward_nodes, values, drop = T]
    
    # Extract weights data
    weights_data_i <- treedata_obj@extraInfo[selected_tipward_nodes, weight_groups]
    
    # Compute weighted mean per weight groups
    mean_weighted_values_i <- apply(X = weights_data_i, MARGIN = 2, FUN = weighted.mean, x = values_data_i)
    mean_weighted_values_i[is.nan(mean_weighted_values_i)] <- NA
    
    # Store data
    weighted_mean_per_groups_along_phylo_array[as.character(time_i),,"weighted_mean"] <- mean_weighted_values_i
    
    # Compute and store sum of weights
    if (keep_weight_sums)
    {
      weight_sum_i <- apply(X = weights_data_i, MARGIN = 2, FUN = sum)
      weighted_mean_per_groups_along_phylo_array[as.character(time_i),,"weight_sum"] <- weight_sum_i
    }
  }

  output <- weighted_mean_per_groups_along_phylo_array
    
  # Melt output array if needed
  if (melt)
  {
    output <- reshape2::melt(output)
  }
  
  # Export output
  return(output)
}


# Extract time scale from LTT data
time_scale <- unique(LTT_per_clades_melted_df$time)
time_scale <- time_scale[order(time_scale, decreasing = T)]

# Inform bioregion variable names
Bioregions_group_var <- c("mean_state_Neotropics", "mean_state_Afrotropics", "mean_state_Eastern Palearctic",
                          "mean_state_Indomalaya", "mean_state_Western Palearctic", "mean_state_Australasia", 
                          "mean_state_Nearctic")
OW_NW_group_var <- c("mean_state_OW", "mean_state_NW")
Trop_Temp_group_var <- c("mean_state_Tropics", "mean_state_Temperate")

# Inform type of data needed
value_types <- c("weighted_mean", "weight_sum")

# Initiate output array
weighted_mean_per_bioregions_all_BAMM_phylo_array <- array(data = NA,
                                                    dim = c(length(BAMM_div_rates_trees), length(time_scale), length(c(Bioregions_group_var, OW_NW_group_var, Trop_Temp_group_var)), length(value_types)),
                                                    dimnames = list(sample = paste0("Sample_",1:length(BAMM_div_rates_trees)), time = time_scale, bioregion = c(Bioregions_group_var, OW_NW_group_var, Trop_Temp_group_var), value_type = value_types))

## Loop per Posterior sample 

for (i in seq_along(BAMM_div_rates_trees))
{
  # i <- 1
  
  # Extract Posterior sample i
  BAMM_div_rates_tree_i <- BAMM_div_rates_trees[[i]]
  
  # Get diversification rates per branches
  BAMM_div_rates_df_i <- data.frame(tipward_node_ID = BAMM_div_rates_tree_i$edge[,2],
                                  div_rates = BAMM_div_rates_tree_i$edge.length) %>% 
                         arrange(tipward_node_ID)
  
  # Assign diversification rates per branches to treedata
  treedata_obj_i <- Ponerinae_phylogeny_1534t_treedata
  treedata_obj_i <- treedata_obj_i %>% 
    left_join(y = BAMM_div_rates_df_i, by = join_by(node == tipward_node_ID))
  
  # Compute weighted rates data per bioregions along the phylogeny
  weighted_rates_per_bioregions_array <- compute_weighted_mean_values_along_phylo(treedata_obj = treedata_obj_i,
                                                                                  values = "div_rates",
                                                                                  weight_groups = c("mean_state_Neotropics", "mean_state_Afrotropics", "mean_state_Eastern Palearctic",
                                                                                                    "mean_state_Indomalaya", "mean_state_Western Palearctic", "mean_state_Australasia", 
                                                                                                    "mean_state_Nearctic"),
                                                                                  time_scale = time_scale,
                                                                                  keep_weight_sums = T,
                                                                                  melt = F)
  # Compute weighted rates data per OW vs. NW along the phylogeny
  weighted_rates_per_OW_vs_NW_array <- compute_weighted_mean_values_along_phylo(treedata_obj = treedata_obj_i,
                                                                                  values = "div_rates",
                                                                                  weight_groups = c("mean_state_OW", "mean_state_NW"),
                                                                                  time_scale = time_scale,
                                                                                  keep_weight_sums = T,
                                                                                  melt = F)
  # Compute weighted rates data per Tropics vs. Temperate along the phylogeny
  weighted_rates_per_Trop_vs_Temp_array <- compute_weighted_mean_values_along_phylo(treedata_obj = treedata_obj_i,
                                                                                values = "div_rates",
                                                                                weight_groups = c("mean_state_Tropics", "mean_state_Temperate"),
                                                                                time_scale = time_scale,
                                                                                keep_weight_sums = T,
                                                                                melt = F)
  
  # Store outputs
  weighted_mean_per_bioregions_all_BAMM_phylo_array[i,,Bioregions_group_var,] <- weighted_rates_per_bioregions_array
  weighted_mean_per_bioregions_all_BAMM_phylo_array[i,,OW_NW_group_var,] <- weighted_rates_per_OW_vs_NW_array
  weighted_mean_per_bioregions_all_BAMM_phylo_array[i,,Trop_Temp_group_var,] <- weighted_rates_per_Trop_vs_Temp_array
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Weighted rates per bioregions computed for posterior sample n°", i, "/", length(BAMM_div_rates_trees),"\n"))
  } 
}

# Rename bioregions
dimnames(weighted_mean_per_bioregions_all_BAMM_phylo_array)$bioregion <- c("Neotropics", "Afrotropics", "Eastern Palearctic",
                                                                           "Indomalaya", "Western Palearctic", "Australasia", 
                                                                           "Nearctic", "OW", "NW", "Tropics", "Temperate")
weighted_mean_per_bioregions_all_BAMM_phylo_array["Sample_150",,,"weighted_mean"]

weighted_mean_per_bioregions_all_BAMM_phylo_array[,"0",,"weighted_mean"]
weighted_mean_per_bioregions_all_BAMM_phylo_array[,"0",,"weight_sum"]

# Save array of weighted mean rates along phylogenies
# saveRDS(weighted_mean_per_bioregions_all_BAMM_phylo_array, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_array.rds")
saveRDS(weighted_mean_per_bioregions_all_BAMM_phylo_array, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_array.rds")


### 10.3/ Melt and filter rates for minimum branch counts ####

weighted_mean_per_bioregions_all_BAMM_phylo_melted <- reshape2::melt(weighted_mean_per_bioregions_all_BAMM_phylo_array, value.name = "div_rate")
weighted_mean_per_bioregions_all_BAMM_phylo_melted <- pivot_wider(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted, names_from = "value_type", values_from = "div_rate")

# Set threshold of minimum number of lineages to compute rates
min_richness_threshold <- 1

weighted_mean_per_bioregions_all_BAMM_phylo_melted$weighted_mean[weighted_mean_per_bioregions_all_BAMM_phylo_melted$weight_sum < min_richness_threshold] <- NA

# Order bioregions
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic",  "Neotropics", "Eastern Palearctic", "Western Palearctic", "OW", "NW", "Tropics", "Temperate")
weighted_mean_per_bioregions_all_BAMM_phylo_melted$bioregion <- factor(weighted_mean_per_bioregions_all_BAMM_phylo_melted$bioregion, levels = bioregion_names, labels = bioregion_names)

# Save melted dataframe
# saveRDS(object = weighted_mean_per_bioregions_all_BAMM_phylo_melted, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")
saveRDS(object = weighted_mean_per_bioregions_all_BAMM_phylo_melted, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")


### 10.4/ Compute median and quantiles rates across posterior samples ####

## Compute medians
weighted_mean_per_bioregions_median <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  group_by(time, bioregion) %>% 
  summarize(median_rates = median(weighted_mean, na.rm = T)) %>%
  ungroup()

## Compute quantiles

# Filter time to avoid failed polygons when plotting
weighted_mean_per_bioregions_quantiles <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  filter(time <= 85)

weighted_mean_per_bioregions_quantiles <- weighted_mean_per_bioregions_quantiles %>% 
  group_by(time, bioregion) %>% 
  # Compute quantiles
  reframe(quant_rates = quantile(weighted_mean, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
  group_by(time, bioregion) %>%
  mutate(quantile = seq(from = 0, to = 1, by = 0.01),
         quantile_ID = c(1:50, 51, 50:1)) %>%
  filter(quantile_ID != 51) %>% # Remove median
  # Assign points ID (order for drawing the polygon)
  group_by(bioregion, quantile_ID) %>%
  arrange(bioregion, quantile_ID, quantile) %>%
  # mutate(points_ID = c(1:100, 200:101)) %>%
  mutate(n_points = dplyr::n()) %>%
  mutate(points_ID = c(1:(first(n_points)/2), first(n_points):((first(n_points)/2) + 1))) %>%
  dplyr::select(-n_points) %>%
  # Reorder by points ID
  arrange(bioregion, quantile_ID, points_ID) %>%
  filter(!is.na(quant_rates)) %>%
  mutate(points_ID = row_number()) %>%
  ungroup()

View(weighted_mean_per_bioregions_quantiles)

## Save medians and quantiles
# saveRDS(object = weighted_mean_per_bioregions_median, "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_median.rds")
# saveRDS(object = weighted_mean_per_bioregions_quantiles, "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_quantiles.rds")
saveRDS(object = weighted_mean_per_bioregions_median, "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_median.rds")
saveRDS(object = weighted_mean_per_bioregions_quantiles, "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_quantiles.rds")

  
### 10.5/ Plot rates through time per bioregions ####

## Set bioregion color scheme
# colors_list_for_states <- readRDS(file = "./outputs/BSM/colors_list_for_states.rds")
# areas_list <- c("A", "U", "E", "I", "R", "N", "W")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Indomalaya", "Nearctic",  "Neotropics", "Eastern Palearctic", "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Filter df to extract bioregion rates (and not aggregated rates)
weighted_mean_per_bioregions_all_BAMM_phylo_melted_all_bioregions <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  filter(bioregion %in% bioregion_names)
weighted_mean_per_bioregions_median_all_bioregions <- weighted_mean_per_bioregions_median %>% 
  filter(bioregion %in% bioregion_names)
weighted_mean_per_bioregions_quantiles_all_bioregions <- weighted_mean_per_bioregions_quantiles %>% 
  filter(bioregion %in% bioregion_names)

## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_per_bioregions_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_per_bioregions_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

ratesT_per_bioregions_ggplot_lines <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_all_bioregions) +
  
  # Plot 1000 line replicates
  geom_line(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_all_bioregions,
            mapping = aes(y = weighted_mean, x = time, group = interaction(bioregion, sample), col = bioregion),
            alpha = 0.01,
            linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_all_bioregions,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = bioregion_names, values = colors_list_for_areas) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_per_bioregions_ggplot_lines)

dev.off()


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_per_bioregions_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_per_bioregions_fuzzy_quantiles.pdf", height = 8, width = 12)

ratesT_per_bioregions_ggplot_quantiles <- ggplot(data = weighted_mean_per_bioregions_quantiles_all_bioregions) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = weighted_mean_per_bioregions_quantiles_all_bioregions,
               mapping = aes(y = quant_rates, x = time, group = interaction(bioregion, quantile_ID), fill = bioregion),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_all_bioregions,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = bioregion_names, values = colors_list_for_areas) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = bioregion_names, values = colors_list_for_areas) +

  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_per_bioregions_ggplot_quantiles)

dev.off()


### 10.6/ Plot rates through time for OW vs. NW ####

## Set OW vs. NW color scheme
colors_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_factors <- c("OW", "NW")
OW_NW_names <- c("Old World", "New World")
names(colors_for_OW_NW) <- OW_NW_names

# Filter df to extract bioregion rates (and not aggregated rates)
weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  filter(bioregion %in% OW_NW_factors)
weighted_mean_per_bioregions_median_OW_vs_NW <- weighted_mean_per_bioregions_median %>% 
  filter(bioregion %in% OW_NW_factors)
weighted_mean_per_bioregions_quantiles_OW_vs_NW <- weighted_mean_per_bioregions_quantiles %>% 
  filter(bioregion %in% OW_NW_factors)

## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

ratesT_OW_vs_NW_ggplot_lines <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW) +
  
  # Plot 1000 line replicates
  geom_line(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW,
            mapping = aes(y = weighted_mean, x = time, group = interaction(bioregion, sample), col = bioregion),
            alpha = 0.01,
            linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_OW_vs_NW,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_OW_vs_NW_ggplot_lines)

dev.off()


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_fuzzy_quantiles.pdf", height = 8, width = 12)

ratesT_OW_vs_NW_ggplot_quantiles <- ggplot(data = weighted_mean_per_bioregions_quantiles_OW_vs_NW) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = weighted_mean_per_bioregions_quantiles_OW_vs_NW,
               mapping = aes(y = quant_rates, x = time, group = interaction(bioregion, quantile_ID), fill = bioregion),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_OW_vs_NW,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_OW_vs_NW_ggplot_quantiles)

dev.off()


### 10.7/ Plot rates through time for Tropics vs. Temperate ####

## Set Tropics vs. Temperate color scheme
colors_for_Trop_Temp <- c("limegreen", "steelblue1")
Trop_Temp_names <- c("Tropics", "Temperate")
names(colors_for_Trop_Temp) <- Trop_Temp_names

# Filter df to extract bioregion rates (and not aggregated rates)
weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  filter(bioregion %in% Trop_Temp_names)
weighted_mean_per_bioregions_median_Trop_vs_Temp <- weighted_mean_per_bioregions_median %>% 
  filter(bioregion %in% Trop_Temp_names)
weighted_mean_per_bioregions_quantiles_Trop_vs_Temp <- weighted_mean_per_bioregions_quantiles %>% 
  filter(bioregion %in% Trop_Temp_names)

## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

ratesT_Trop_vs_Temp_ggplot_lines <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp) +
  
  # Plot 1000 line replicates
  geom_line(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp,
            mapping = aes(y = weighted_mean, x = time, group = interaction(bioregion, sample), col = bioregion),
            alpha = 0.01,
            linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_Trop_vs_Temp,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = Trop_Temp_names, labels = Trop_Temp_names, values = colors_for_Trop_Temp) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_Trop_vs_Temp_ggplot_lines)

dev.off()


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_fuzzy_quantiles.pdf", height = 8, width = 12)

ratesT_Trop_vs_Temp_ggplot_quantiles <- ggplot(data = weighted_mean_per_bioregions_quantiles_Trop_vs_Temp) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = weighted_mean_per_bioregions_quantiles_Trop_vs_Temp,
               mapping = aes(y = quant_rates, x = time, group = interaction(bioregion, quantile_ID), fill = bioregion),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_Trop_vs_Temp,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = Trop_Temp_names, labels = Trop_Temp_names, values = colors_for_Trop_Temp) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = Trop_Temp_names, labels = Trop_Temp_names, values = colors_for_Trop_Temp) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(ratesT_Trop_vs_Temp_ggplot_quantiles)

dev.off()



##### 11/ STRAPP test for differences in diversification rates across Bioregions #####

?traitDependentBAMM

# STRAPP = STructured Rate Permutations on Phylogenies
# Alternative to SSE models (BiSSE, HiSSE, QuSSE, GeoSSE) that account for the number of independent shifts in diversifciation regimes
# SSE models may be significative with only a single shift (but if this shift coincidence with a trait shift, why not... 
# The issue is more about false positive Type I error leading to significant support for trait that are not synapomorphy of the focla clade, but of a wider clade including the focal clade)
# See also the Canonical Phylogenetic Ordination (CPO) in Title et al., 2024 for an alternative dedicated to detecting single diversification shift associated with trait shift

#	Calculates a correlation coefficient and test for significance by permutating diversification rates on the phylogeny 
# while preserving the underlying covariation in macroevolutionary rates and species traits.

# Work with discrete and continuous traits
  # 1/ Run BAMM (or any other diversification models that involve discrete shifts) to estimate mean branch rates and rate shift locations
  # 2/ Compute a statistic of relationship between tip rates ~ tip trait values
  # 3/ Use a structured permutation to build a null model of tip rates by permutating tip rates by taxa-blocks based on macroevolutionary regimes
  # 4/ Compare observed statistic to its null distribution on the permuted tip values

### 11.1/ Assign monomorphic Biogeographic states to tips ####

# Load biogeographic data
Taxa_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Taxa_bioregions_binary_table_for_analyses_7_regions_PaleA.rds")

# Load updated range table
range_size <- apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 1, FUN = sum)
table(range_size) # Proportion of multiarea ranges is limited

## Assign tips to the bioregion with the largest area intersecting the taxa range

# Load alpha-hull range data
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

# Load bioregion sf
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")

## Transform binary rasters in binary sf ranges

# Initiate final sf object
Ponerinae_taxa_ranges_sf <- data.frame()
# Loop per layer
for (i in 1:raster::nlayers(Ponerinae_alpha_hull_stack_WGS84))
{
  # i <- 1
  
  # Extract raster
  raster_i <- Ponerinae_alpha_hull_stack_WGS84[[i]]
  # Convert to terra
  terra_raster_i <- as(raster_i, "SpatRaster")
  # Convert in Polygons
  terra_poly_i <- terra::as.polygons(terra_raster_i)
  # Convert to sf
  sf_poly_i <- st_as_sf(terra_poly_i)
  if (nrow(sf_poly_i) > 1) 
  {
    # Remove absence range
    names(sf_poly_i)[1] <- "Taxa"
    sf_poly_i <- sf_poly_i %>% 
      filter(Taxa == 1)
    sf_poly_i$Taxa <- names(raster_i)
    
    # Store
    Ponerinae_taxa_ranges_sf <- rbind(Ponerinae_taxa_ranges_sf, sf_poly_i)
  } 
  
  # Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Range polygonized for taxa n°", i, "/", raster::nlayers(Ponerinae_alpha_hull_stack_WGS84),"\n"))
  }
}
# Save range sf file (missing some taxa with no range in the alpha-hull)
saveRDS(Ponerinae_taxa_ranges_sf, file = "./outputs/Species_richness_maps/Ponerinae_taxa_ranges_sf.rds")

## Rank areas of intersection with bioregions

sf_use_s2(FALSE)
Ponerinae_taxa_ranges_per_bioregions_sf <- st_intersection(Ponerinae_taxa_ranges_sf, Bioregions_sf_Bioregions_level)
Ponerinae_taxa_ranges_per_bioregions_sf$area <- st_area(Ponerinae_taxa_ranges_per_bioregions_sf)
sf_use_s2(TRUE)

Ponerinae_taxa_ranges_per_bioregions_sf <- Ponerinae_taxa_ranges_per_bioregions_sf %>% 
  group_by(Taxa) %>% 
  arrange(Taxa, area) %>%
  mutate(Bioregion_ranks = row_number()) %>%
  ungroup()

# Save range sf file per bioregions (missing some taxa with no range in the alpha-hull)
saveRDS(Ponerinae_taxa_ranges_per_bioregions_sf, file = "./outputs/Species_richness_maps/Ponerinae_taxa_ranges_per_bioregions_sf.rds")

## Extract the top ranked bioregions

Taxa_bioregions_per_area <- Ponerinae_taxa_ranges_per_bioregions_sf %>% 
  st_drop_geometry() %>%
  filter(Bioregion_ranks == 1) %>%
  dplyr::select(Taxa, Bioregion) %>%
  rename(Bioregion_from_area = Bioregion)
  
# Save bioregion membership based on area
saveRDS(Taxa_bioregions_per_area, file = "./outputs/Species_richness_maps/Taxa_bioregions_per_area.rds")

## Get Bioregion tips data for monomorphic taxa

range_size <- apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 1, FUN = sum)

Taxa_bioregions_first <- names(Taxa_bioregions_binary_table[, -1])[apply(X = Taxa_bioregions_binary_table[, -1], MARGIN = 1, FUN = which.max)]

Taxa_bioregions_tips_data <- data.frame(Taxa = Taxa_bioregions_binary_table$Current_name, Bioregion = NA)
Taxa_bioregions_tips_data$Bioregion[range_size == 1] <- Taxa_bioregions_first[range_size == 1]

## Update tip data for multi-areas cases, based on larger area membership

Taxa_bioregions_tips_data <- left_join(Taxa_bioregions_tips_data, Taxa_bioregions_per_area)
Taxa_bioregions_tips_data$Bioregion[is.na(Taxa_bioregions_tips_data$Bioregion)] <- Taxa_bioregions_tips_data$Bioregion_from_area[is.na(Taxa_bioregions_tips_data$Bioregion)]

Taxa_bioregions_tips_data <- Taxa_bioregions_tips_data %>% 
  dplyr::select(Taxa, Bioregion)
Taxa_bioregions_tips_data$Bioregion <- as.factor(Taxa_bioregions_tips_data$Bioregion)

## Aggregate bioregions in higher classification

NW_bioregions <- c("Neotropics", "Nearctic")

Taxa_bioregions_tips_data$OW_vs_NW <- "Old_World"
Taxa_bioregions_tips_data$OW_vs_NW[Taxa_bioregions_tips_data$Bioregion %in% NW_bioregions] <- "New_World"
Taxa_bioregions_tips_data$OW_vs_NW <- as.factor(Taxa_bioregions_tips_data$OW_vs_NW)

table(Taxa_bioregions_tips_data$OW_vs_NW)

Tropical_bioregions <- c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")

Taxa_bioregions_tips_data$Trop_vs_Temp <- "Temperate"
Taxa_bioregions_tips_data$Trop_vs_Temp[Taxa_bioregions_tips_data$Bioregion %in% Tropical_bioregions] <- "Tropical"
Taxa_bioregions_tips_data$Trop_vs_Temp <- as.factor(Taxa_bioregions_tips_data$Trop_vs_Temp)

table(Taxa_bioregions_tips_data$Trop_vs_Temp)

# Save Bioregion tip data
saveRDS(object = Taxa_bioregions_tips_data, file = "./outputs/Species_richness_maps/Taxa_bioregions_tips_data.rds")


### 11.2/ Test across all bioregions ####

# Load Bioregion tip data
Taxa_bioregions_tips_data <- readRDS(file = "./outputs/Species_richness_maps/Taxa_bioregions_tips_data.rds")

# Assign tip names to tip data
Ponerinae_bioregion_tip_data_all_bioregions <- Taxa_bioregions_tips_data$Bioregion
names(Ponerinae_bioregion_tip_data_all_bioregions) <- Taxa_bioregions_tips_data$Taxa

# Reorder as in phylogeny
Ponerinae_bioregion_tip_data_all_bioregions <- Ponerinae_bioregion_tip_data_all_bioregions[match(phy$tip.label, names(Ponerinae_bioregion_tip_data_all_bioregions))]
identical(names(Ponerinae_bioregion_tip_data_all_bioregions), phy$tip.label)

set.seed(seed = 1234)
# Run STRAPP test
STRAPP_test_all_bioregions <- traitDependentBAMM(ephy = BAMM_posterior_samples_data,
                                                 traits = Ponerinae_bioregion_tip_data_all_bioregions,
                                                 reps = 1000, rate = "net diversification",
                                                 return.full = T,
                                                 method = "kruskal", # For categorical multinomial data (G > 2)
                                                 logrates = F, # Do not use log for Kruskal-Wallis as it is a rank test, so there is no need, and it prevents removing the negative rates
                                                 two.tailed = T)

STRAPP_test_all_bioregions$estimate # Mean tip rates per categories
STRAPP_test_all_bioregions$p.value
STRAPP_test_all_bioregions$obs.corr # Observed statistic for each posterior sample
STRAPP_test_all_bioregions$gen # Generation ID of the selected posterior sample
STRAPP_test_all_bioregions$null # Null statistic for each posterior sample

# Save STRAPP test output
# saveRDS(STRAPP_test_all_bioregions, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_all_bioregions.rds")
saveRDS(STRAPP_test_all_bioregions, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_test_all_bioregions.rds")

## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_all_bioregions$test_stats <- STRAPP_test_all_bioregions$obs.corr - STRAPP_test_all_bioregions$null
summary(STRAPP_test_all_bioregions$test_stats)
table(STRAPP_test_all_bioregions$test_stats > 0)
STRAPP_test_all_bioregions$stat_median <- median(STRAPP_test_all_bioregions$test_stats)
STRAPP_test_all_bioregions$stat_Q5 <- quantile(STRAPP_test_all_bioregions$test_stats, p = 0.05)

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_all_bioregions_hist.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_test_all_bioregions_hist.pdf", height = 6, width = 8)

par(mar = c(5.1, 5.1, 4.1, 2.1)) # bltr

hist_data <- hist(x = STRAPP_test_all_bioregions$test_stats,
                  breaks = 30, freq = TRUE, col = "gray", 
                  main = "STRAPP test for\nCurrent rates difference across Bioregions", 
                  # xlab = "Δ Kruskal-Wallis's Khi²",
                  xlab = expression(paste(Delta,"  Kruskal-Wallis's Khi²")),
                  cex.main = 1.3, cex.axis = 1.3, cex.lab = 1.4,
                  lwd = 2)

# Identify axis ranges
x_min_max <- range(hist_data$breaks)
x_range <- x_min_max[2] - x_min_max[1]
y_min_max <- range(hist_data$counts)
y_range <- y_min_max[2] - y_min_max[1]

# Draw control lines
abline(v = 0, lty = 2, lwd = 2, col = "black")
abline(v = as.numeric(STRAPP_test_all_bioregions$stat_Q5), lty = 2, lwd = 2, col = "red")

# Draw blank rectangle for the legend
rect(xleft = x_min_max[1], x_min_max[1] + x_range * 0.45,
     ybottom = y_min_max[1] + y_range * 0.50, ytop = y_min_max[2] * 1.05,
     col = "white", border = NA)

# Reorder estimates
rate_estimates_per_bioregions <- STRAPP_test_all_bioregions$estimate[bioregion_names]

# Insert estimates of net div rates
for (i in seq_along(STRAPP_test_all_bioregions$estimate))
{
  text(labels = paste0(names(STRAPP_test_all_bioregions$estimate)[i], " rate = ", round(STRAPP_test_all_bioregions$estimate[[i]], 3)),
       cex = 1.2, adj = 0,
       x = x_min_max[1] + x_range * 0.00, y = y_min_max[1] + y_range * (1.05 - (i * 0.07)))
}

# Insert quantiles legend
legend(legend = c(paste0("CI 5% = ", round(as.numeric(STRAPP_test_all_bioregions$stat_Q5), 1)), 
                  paste0("p = ", STRAPP_test_all_bioregions$p.value)),
       x = "topleft", inset = c(0.02, 0.55), lty = 2 , lwd = 2,
       col = c("red", NA), cex = 1.2, bty ="n")

dev.off()

### 11.3/ Test Old World vs. New World ####

# Assign tip names to tip data
Ponerinae_bioregion_tip_data_OW_vs_NW <- Taxa_bioregions_tips_data$OW_vs_NW
names(Ponerinae_bioregion_tip_data_OW_vs_NW) <- Taxa_bioregions_tips_data$Taxa

# Reorder as in phylogeny
Ponerinae_bioregion_tip_data_OW_vs_NW <- Ponerinae_bioregion_tip_data_OW_vs_NW[match(phy$tip.label, names(Ponerinae_bioregion_tip_data_OW_vs_NW))]
identical(names(Ponerinae_bioregion_tip_data_OW_vs_NW), phy$tip.label)

set.seed(seed = 1234)
# Run STRAPP test
STRAPP_test_OW_vs_NW <- traitDependentBAMM(ephy = BAMM_posterior_samples_data,
                                                 traits = Ponerinae_bioregion_tip_data_OW_vs_NW,
                                                 reps = 1000, rate = "net diversification",
                                                 return.full = T,
                                                 method = "mann-whitney", # For categorical binomial data (G = 2)
                                                 logrates = F, # Do not use log for Mann-Whitney as it is a rank test, so there is no need, and it prevents removing the negative rates
                                                 two.tailed = T)

STRAPP_test_OW_vs_NW$estimate # Mean tip rates per categories
STRAPP_test_OW_vs_NW$p.value
STRAPP_test_OW_vs_NW$obs.corr # Observed statistic for each posterior sample
STRAPP_test_OW_vs_NW$gen # Generation ID of the selected posterior sample
STRAPP_test_OW_vs_NW$null # Null statistic for each posterior sample

# Save STRAPP test output
# saveRDS(STRAPP_test_OW_vs_NW, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_OW_vs_NW.rds")
saveRDS(STRAPP_test_OW_vs_NW, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_test_OW_vs_NW.rds")

## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_OW_vs_NW$test_stats <- STRAPP_test_OW_vs_NW$null - STRAPP_test_OW_vs_NW$obs.corr
summary(STRAPP_test_OW_vs_NW$test_stats)
table(STRAPP_test_OW_vs_NW$test_stats > 0)
STRAPP_test_OW_vs_NW$stat_median <- median(STRAPP_test_OW_vs_NW$test_stats)
STRAPP_test_OW_vs_NW$stat_Q5 <- quantile(STRAPP_test_OW_vs_NW$test_stats, p = 0.05)

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_OW_vs_NW_hist.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_test_OW_vs_NW_hist.pdf", height = 6, width = 8)

par(mar = c(5.1, 5.1, 4.1, 2.1)) # bltr

hist_data <- hist(x = STRAPP_test_OW_vs_NW$test_stats,
                  breaks = 30, freq = TRUE, col = "gray", 
                  main = "STRAPP test\nDifference in current net div. rates\nOW vs. NW", 
                  # xlab = "Δ Mann-Whitney's U",
                  xlab = expression(paste(Delta,"  Mann-Whitney's U")),
                  cex.main = 1.3, cex.axis = 1.3, cex.lab = 1.4,
                  lwd = 2)

# Identify axis ranges
x_min_max <- range(hist_data$breaks)
x_range <- x_min_max[2] - x_min_max[1]
y_min_max <- range(hist_data$counts)
y_range <- y_min_max[2] - y_min_max[1]

# Draw control lines
abline(v = 0, lty = 2, lwd = 2, col = "black")
abline(v = as.numeric(STRAPP_test_OW_vs_NW$stat_Q5), lty = 2, lwd = 2, col = "red")

# Draw blank rectangle for the legend
rect(xleft = x_min_max[1] + x_range * 0.85, x_min_max[2] ,
     ybottom = y_min_max[1] + y_range * 0.68, ytop = y_min_max[2] * 1.05,
     col = "white", border = NA)

# Insert estimates of net div rates
text(labels = c(paste0("OW rate = ", round(STRAPP_test_OW_vs_NW$estimate[[1]], 3)), 
                paste0("\n\nNW rate = ", round(STRAPP_test_OW_vs_NW$estimate[[2]], 3))),
     cex = 1.2,
     x = x_min_max[1] + x_range * 0.89, y = y_min_max[1] + y_range * 1.00)

# Insert quantiles legend
legend(legend = c(paste0("CI 5% = ", round(as.numeric(STRAPP_test_OW_vs_NW$stat_Q5), 0)), 
                  paste0("p = ", STRAPP_test_OW_vs_NW$p.value)),
       x = "topright", inset = c(0.02, 0.11), lty = 2 , lwd = 2,
       col = c("red", NA), cex = 1.2, bty ="n")

dev.off()

### 11.4/ Test Tropics vs. Temperate bioregions ####

# Assign tip names to tip data
Ponerinae_bioregion_tip_data_Trop_vs_Temp <- Taxa_bioregions_tips_data$Trop_vs_Temp
names(Ponerinae_bioregion_tip_data_Trop_vs_Temp) <- Taxa_bioregions_tips_data$Taxa

# Reorder as in phylogeny
Ponerinae_bioregion_tip_data_Trop_vs_Temp <- Ponerinae_bioregion_tip_data_Trop_vs_Temp[match(phy$tip.label, names(Ponerinae_bioregion_tip_data_Trop_vs_Temp))]
identical(names(Ponerinae_bioregion_tip_data_Trop_vs_Temp), phy$tip.label)

set.seed(seed = 1234)
# Run STRAPP test
STRAPP_test_Trop_vs_Temp <- traitDependentBAMM(ephy = BAMM_posterior_samples_data,
                                           traits = Ponerinae_bioregion_tip_data_Trop_vs_Temp,
                                           reps = 1000, rate = "net diversification",
                                           return.full = T,
                                           method = "mann-whitney", # For categorical multinomial data (G > 2)
                                           logrates = F, # Do not use log for Mann-Whitney as it is a rank test, so there is no need, and it prevents removing the negative rates 
                                           two.tailed = T)

STRAPP_test_Trop_vs_Temp$estimate # Mean tip rates per categories
STRAPP_test_Trop_vs_Temp$p.value
STRAPP_test_Trop_vs_Temp$obs.corr # Observed statistic for each posterior sample
STRAPP_test_Trop_vs_Temp$gen # Generation ID of the selected posterior sample
STRAPP_test_Trop_vs_Temp$null # Null statistic for each posterior sample

# Save STRAPP test output
# saveRDS(STRAPP_test_Trop_vs_Temp, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_Trop_vs_Temps.rds")
saveRDS(STRAPP_test_Trop_vs_Temp, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_Trop_vs_Temps.rds")


## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_Trop_vs_Temp$test_stats <- STRAPP_test_Trop_vs_Temp$null - STRAPP_test_Trop_vs_Temp$obs.corr
summary(STRAPP_test_Trop_vs_Temp$test_stats)
table(STRAPP_test_Trop_vs_Temp$test_stats > 0)
STRAPP_test_Trop_vs_Temp$stat_median <- median(STRAPP_test_Trop_vs_Temp$test_stats)
STRAPP_test_Trop_vs_Temp$stat_Q5 <- quantile(STRAPP_test_Trop_vs_Temp$test_stats, p = 0.05)

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_test_Trop_vs_Temp_hist.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_test_Trop_vs_Temp_hist.pdf", height = 6, width = 8)

par(mar = c(5.1, 5.1, 4.1, 2.1)) # bltr

hist_data <- hist(x = STRAPP_test_Trop_vs_Temp$test_stats,
                  breaks = 30, freq = TRUE, col = "gray", 
                  main = "STRAPP test\nDifference in current net div. rates\nTropics vs. Temperate", 
                  # xlab = "Δ Mann-Whitney's U",
                  xlab = expression(paste(Delta,"  Mann-Whitney's U")),
                  cex.main = 1.3, cex.axis = 1.3, cex.lab = 1.4,
                  lwd = 2)

# Identify axis ranges
x_min_max <- range(hist_data$breaks)
x_range <- x_min_max[2] - x_min_max[1]
y_min_max <- range(hist_data$counts)
y_range <- y_min_max[2] - y_min_max[1]

# Draw control lines
abline(v = 0, lty = 2, lwd = 2, col = "black")
abline(v = as.numeric(STRAPP_test_Trop_vs_Temp$stat_Q5), lty = 2, lwd = 2, col = "red")

# Draw blank rectangle for the legend
rect(xleft = x_min_max[1], x_min_max[1] + x_range * 0.33,
     ybottom = y_min_max[1] + y_range * 0.68, ytop = y_min_max[2] * 1.05,
     col = "white", border = NA)

# Insert estimates of net div rates
text(labels = c(paste0("Trop. rate = ", round(STRAPP_test_Trop_vs_Temp$estimate[[1]], 3)), 
                paste0("\n\nTemp. rate = ", round(STRAPP_test_Trop_vs_Temp$estimate[[2]], 3))),
     cex = 1.2, adj = 0,
     x = x_min_max[1] + x_range * 0.01, y = y_min_max[1] + y_range * 0.96)

# Insert quantiles legend
legend(legend = c(paste0("CI 5% = ", round(as.numeric(STRAPP_test_Trop_vs_Temp$stat_Q5), 0)), 
                  paste0("p = ", STRAPP_test_Trop_vs_Temp$p.value)),
       x = "topleft", inset = c(0.02, 0.17), lty = 2 , lwd = 2,
       col = c("red", NA), cex = 1.2, bty ="n")

dev.off()


##### 12/ STRAPP tests along evolutionary time #####

# The STRAPP test can in theory be run for any point in time
# Cut posterior tree shift configurations to a specific age
# Run STRAPP test using blocks found in the cut trees

## Estimates for mean rates per categories may vary from what is shown on the plot
# as "unique areas" are used for categories in STRAPP tests
# while weighted area membership are used for plotting rates through time

### 12.1/ Compute density maps of aggregated bioregions

# Load simmaps with only unique areas
# DEC_J_simmaps_unique_areas <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")
DEC_J_simmaps_unique_areas <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_simmaps_unique_areas.rds")

# Load a custom function to produce densityMap
source(file = "./functions/densityMap_custom.R")

# Load labels of areas
all_areas <- c("N", "A", "E", "I", "W", "U", "R")
all_areas_labels <- c("Neotropics", "Afrotropics", "Eastern Palearctic", "Indomalaya", "Western Palearctic", "Australasia", "Nearctic")

# DEC_J_simmaps_unique_areas[[1]]$maps

## 12.1.1/ Compute density maps for Old World vs. New World ####

# Define groups of bioregions
NW_bioregions <- c("Neotropics", "Nearctic")
OW_bioregions <- c("Afrotropics", "Australasia", "Indomalaya", "Eastern Palearctic", "Western Palearctic")

# Extract the focal areas for NW
focal_areas <- all_areas[all_areas_labels %in% NW_bioregions]
focal_area_labels <- all_areas_labels[all_areas_labels %in% NW_bioregions]

# Define other areas by contrast
other_areas <- all_areas[!(all_areas %in% focal_areas)]

## Aggregate states in simmaps
DEC_J_simmaps_binary_OW_vs_NW <- DEC_J_simmaps_unique_areas
DEC_J_simmaps_binary_OW_vs_NW <- lapply(X = DEC_J_simmaps_binary_OW_vs_NW, FUN = phytools::mergeMappedStates, old.states = other_areas, new.state = "0")
DEC_J_simmaps_binary_OW_vs_NW <- lapply(X = DEC_J_simmaps_binary_OW_vs_NW, FUN = phytools::mergeMappedStates, old.states = focal_areas, new.state = "1")

class(DEC_J_simmaps_binary_OW_vs_NW) <- c("list", "multiSimmap", "multiPhylo")

# Check that remaining states are all binary
# unique(unlist(lapply(X = DEC_J_simmaps_binary_OW_vs_NW, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
# DEC_J_simmaps_binary_OW_vs_NW[[1]]$maps

## Estimate the posterior probabilities of states along all branches (from the set of simulated maps) ####

# Use a custom version to have better control
source(file = "./functions/densityMap_custom.R")

# # Create a "densityMap" object
# DEC_J_densityMap_OW_vs_NW <- phytools::densityMap(trees = DEC_J_simmaps_binary_OW_vs_NW,
#                                          plot = TRUE)

# Create a "densityMap" object
DEC_J_density_map_OW_vs_NW <- densityMap_custom(trees = DEC_J_simmaps_binary_OW_vs_NW,
                                                tol = 1e-5, verbose = T,
                                                col_scale = NULL,
                                                plot = FALSE)

## Update color gradient and states

colors_for_OW_NW <- c("peachpuff2", "mediumpurple2")

# Set colors for areas/bioregions
col_fn <- colorRampPalette(colors = colors_for_OW_NW)
col_scale <- col_fn(n = 1001)

# Update color gradient and states
DEC_J_density_map_OW_vs_NW <- setMap(DEC_J_density_map_OW_vs_NW, col_scale)
DEC_J_density_map_OW_vs_NW$states <- c("Old World", "New World")

## Plot result
plot(DEC_J_density_map_OW_vs_NW)

## Save the resulting Density map with updated color gradient
# saveRDS(object = DEC_J_density_map_OW_vs_NW, file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))
saveRDS(object = DEC_J_density_map_OW_vs_NW, file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))

## Load the resulting Density map with updated color gradient
# DEC_J_density_map_OW_vs_NW <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))
DEC_J_density_map_OW_vs_NW <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))


## 12.1.2/ Compute density maps for Tropics vs. Temperate ####

# Define groups of bioregions
Trop_bioregions <- c("Afrotropics", "Australasia", "Indomalaya", "Neotropics")
Temp_bioregions <- c("Nearctic", "Eastern Palearctic", "Western Palearctic")

# Extract the focal areas for Tropics
focal_areas <- all_areas[all_areas_labels %in% Trop_bioregions]
focal_area_labels <- all_areas_labels[all_areas_labels %in% Trop_bioregions]

# Define other areas by contrast
other_areas <- all_areas[!(all_areas %in% focal_areas)]

## Aggregate states in simmaps
DEC_J_simmaps_binary_Trop_vs_Temp <- DEC_J_simmaps_unique_areas
DEC_J_simmaps_binary_Trop_vs_Temp <- lapply(X = DEC_J_simmaps_binary_Trop_vs_Temp, FUN = phytools::mergeMappedStates, old.states = other_areas, new.state = "0")
DEC_J_simmaps_binary_Trop_vs_Temp <- lapply(X = DEC_J_simmaps_binary_Trop_vs_Temp, FUN = phytools::mergeMappedStates, old.states = focal_areas, new.state = "1")

class(DEC_J_simmaps_binary_Trop_vs_Temp) <- c("list", "multiSimmap", "multiPhylo")

# Check that remaining states are all binary
# unique(unlist(lapply(X = DEC_J_simmaps_binary_Trop_vs_Temp, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
# DEC_J_simmaps_binary_Trop_vs_Temp[[1]]$maps

## Estimate the posterior probabilities of states along all branches (from the set of simulated maps) ####

# Use a custom version to have better control
source(file = "./functions/densityMap_custom.R")

# # Create a "densityMap" object
# DEC_J_densityMap_Trop_vs_Temp <- phytools::densityMap(trees = DEC_J_simmaps_binary_Trop_vs_Temp,
#                                                       plot = TRUE)

# Create a "densityMap" object
DEC_J_density_map_Trop_vs_Temp <- densityMap_custom(trees = DEC_J_simmaps_binary_Trop_vs_Temp,
                                                    tol = 1e-5, verbose = T,
                                                    col_scale = NULL,
                                                    plot = FALSE)

## Update color gradient and states

colors_for_Temp_vs_Trop <- c("steelblue1", "limegreen")

# Set colors for areas/bioregions
col_fn <- colorRampPalette(colors = colors_for_Temp_vs_Trop)
col_scale <- col_fn(n = 1001)

# Update color gradient and states
DEC_J_density_map_Trop_vs_Temp <- setMap(DEC_J_density_map_Trop_vs_Temp, col_scale)
DEC_J_density_map_Trop_vs_Temp$states <- c("Temperate", "Tropics")

## Plot result
plot(DEC_J_density_map_Trop_vs_Temp)

## Save the resulting Density map with updated color gradient
# saveRDS(object = DEC_J_density_map_Trop_vs_Temp, file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_Trop_vs_Temp.rds"))
saveRDS(object = DEC_J_density_map_Trop_vs_Temp, file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_Trop_vs_Temp.rds"))


### 12.2/ Functions to run STRAPP test along evolutionary time ####

## 12.2.1/ Function to update tipStates and tipRates according to time ####

update_tipStates_and_tipRates_for_focal_time <- function (BAMM_object, time, update_rates = T, verbose = T)
{
  # Add "phylo" class to the BAMM_object
  class(BAMM_object) <- c("bammdata", "phylo")
  
  ## Identify edges present at focal time
  
  # Edge, rootward_node, tipward_node, length (once cut)
  
  # Get node ages per branch (no root edge)
  all_edges_df <- phytools::nodeHeights(BAMM_object)
  root_age <- max(phytools::nodeHeights(BAMM_object)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)
  
  all_edges_ID_df <- BAMM_object$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df) %>% 
    select("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")
  
  # Detect root node ID
  root_node_ID <- length(BAMM_object$tip.label) + 1
  
  # Identify edges present at time i
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= time
  all_edges_df$time_test <- all_edges_df$rootward_test & all_edges_df$tipward_test
  all_edges_df$length <- all_edges_df$rootward_node_age - time
  # all_edges_df$length[!all_edges_df$time_test] <- NA
  
  # Initiate regime ID
  all_edges_df$regime_ID <- NA
  
  # Extract only edge that are present
  # present_edges_df <- all_edges_df[all_edges_df$time_test, ]
  
  ## Loop per Posterior sample
  for (i in seq_along(BAMM_object$eventData))
  {
    # i <- 1
    
    # Extract eventData records
    eventData_i <- BAMM_object$eventData[[i]] # Tipward node ID
    
    # Compute regime age and updated length
    eventData_i$age <- root_age - eventData_i$time
    eventData_i$updated_length <- eventData_i$age - time
      
    ## Identify edge ID per regimes
    ## Loop per regime
    for (j in 1:nrow(eventData_i))
    {
      # j <- 4
      
      tipward_node_ID_j <- eventData_i$node[j]
      
      # Get descendant tipward nodes of regime j
      regime_nodes_j <- getDescendants(tree = BAMM_object, node = tipward_node_ID_j)
      
      # Assign regime ID
      all_edges_df$regime_ID[all_edges_df$tipward_node_ID %in% regime_nodes_j] <- j
      
      # Deal with special case of the edge where the process starts
      # Should the edge where the process starts be included in the regime at the focal time?
      if (j != 1) # No need for the root process
      {
        # Identify the starting edge
        starting_edge_j <- as.numeric(all_edges_df$edge_ID[all_edges_df$tipward_node_ID == tipward_node_ID_j])
        
        # Get relative position of the regime shift
        relative_position_shift_j <- all_edges_df$rootward_node_age[starting_edge_j] - eventData_i$age[j]
        # Assign starting edge to process only if the regime shift happen before the time cut
        if (relative_position_shift_j < all_edges_df$length[starting_edge_j])
        {
          all_edges_df$regime_ID[starting_edge_j] <- j
        }
      }
    }
    
    # Update tipStates by providing only regimes for tips that are present at the focal time
    tipStates_i <- all_edges_df$regime_ID[all_edges_df$time_test]
    names(tipStates_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
    BAMM_object$tipStates[[i]] <- tipStates_i
    
    ## If needed, also update tipRates
    if (update_rates)
    {
      eventData_i$tip_speciation_rates <- NA
      eventData_i$tip_extinction_rates <- NA
      
      ## Loop per regime
      for (j in 1:nrow(eventData_i))
      {
        # Compute new tip speciation rates based on regime parameters 
        lambda_0_j <- eventData_i$lam1[j]
        alpha_j <- eventData_i$lam2[j]
        time_j <- eventData_i$updated_length[j]
        
        if (alpha_j <= 0) # If alpha <= 0 (decrease): lambda_t = lambda_0 * exp(alpha*t)
        {
          eventData_i$tip_speciation_rates[j] <- lambda_0_j * exp(alpha_j*time_j)
        } else { # If alpha > 0 (increase): lambda_t = lambda_0 * (2 - exp(-alpha*t))
          eventData_i$tip_speciation_rates[j] <- lambda_0_j * (2 - exp(-alpha_j*time_j))
        }
        
        # Compute new tip extinction rates based on regime parameters 
          # All extinction rates are constant within regime in the current BAMM settings
        eventData_i$tip_extinction_rates[j] <- eventData_i$mu1[j]
        
        if (time_j < 0)
        {
          eventData_i$tip_speciation_rates[j] <- NA
          eventData_i$tip_extinction_rates[j] <- NA
        }
      }
      
      # Assign rates to edge according to regime ID
      all_edges_df$tipLambda <- NA
      all_edges_df$tipLambda <- eventData_i$tip_speciation_rates[match(x = all_edges_df$regime_ID, table = eventData_i$index)]
      all_edges_df$tipMu <- NA
      all_edges_df$tipMu <- eventData_i$tip_extinction_rates[match(x = all_edges_df$regime_ID, table = eventData_i$index)]
      
      # Update tipLambda and tipMu by providing only regimes for tips that are present at the focal time
      tipLambda_i <- all_edges_df$tipLambda[all_edges_df$time_test]
      names(tipLambda_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
      BAMM_object$tipLambda[[i]] <- tipLambda_i
      
      tipMu_i <- all_edges_df$tipMu[all_edges_df$time_test]
      names(tipMu_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
      BAMM_object$tipMu[[i]] <- tipMu_i
    }
    
    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Tip states/rates updated for BAMM posterior sample n°", i, "/", length(BAMM_object$eventData),"\n"))
    }
  }
  
  # Update tip labels
  BAMM_object$tip.label <- all_edges_df$edge_ID[all_edges_df$time_test]
  # Inform focal time
  BAMM_object$focal_time <- time
  
  # Export updated BAMM_object
  return(BAMM_object)
}


## 12.2.2/ Function to update states PP according to time ####

get_most_likely_binary_states_for_focal_time <- function (density_map, time)
{
  
  ## Identify edges present at focal time
  
  # Edge, rootward_node, tipward_node, length (once cut)
  
  # Get node ages per branch (no root edge)
  all_edges_df <- phytools::nodeHeights(density_map$tree)
  root_age <- max(phytools::nodeHeights(density_map$tree)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)
  
  all_edges_ID_df <- density_map$tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df) %>% 
    select("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")
  
  # Detect root node ID
  root_node_ID <- length(density_map$tree$tip.label) + 1
  
  # Identify edges present at time i
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= time
  all_edges_df$time_test <- all_edges_df$rootward_test & all_edges_df$tipward_test
  all_edges_df$length <- all_edges_df$rootward_node_age - time
  # all_edges_df$length[!all_edges_df$time_test] <- NA
  
  # If no edge present, send warning
  if (sum(all_edges_df$time_test) == 0)
  {
    warning(paste0("No branch are present for time = ", time, ". Return an empty table.\n"))
    
    # Return a NULL object
    trait_data <- NULL
    return(trait_data)
    
  } else {
    
    # Extract only edge that are present
    present_edges_df <- all_edges_df[all_edges_df$time_test, ]
    present_edges_df$current_state_binary <- NA
    present_edges_df$current_state <- NA
    
    # Loop per edge
    for (i in 1:nrow(present_edges_df))
    {
      # i <- 1
      
      # Extract edge ID
      edge_ID_i <- as.numeric(present_edges_df$edge_ID[i])
      # Extract time to cut
      edge_length_i <- present_edges_df$length[i]
      # Extract edge map
      edge_map_i <- density_map$tree$maps[[edge_ID_i]]
      
      # Identify segment matching the time cut
      cut_position_i <- which.min(cumsum(edge_map_i) < edge_length_i)
      # Identify current state
      current_state_i <- names(edge_map_i)[cut_position_i]
      present_edges_df$current_state_binary[i] <- current_state_i
      
      # Convert into state names selecting the most likely state
      if (as.numeric(current_state_i) < 500)
      {
        present_edges_df$current_state[i] <- density_map$states[1]
      } else {
        present_edges_df$current_state[i] <- density_map$states[2]
      }
    }
    
    # table(present_edges_df$current_state)
    
    # Format output = named vector of most likely states
    trait_data <- present_edges_df$current_state
    names(trait_data) <- present_edges_df$edge_ID
    
    # Export trait data
    return(trait_data)
  }
  
  
}


### 12.3/ STRAPP test along evolutionary time for OW_vs_NW #### 

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

## Load the binary density map for OW vs. NW
# DEC_J_density_map_OW_vs_NW <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))
DEC_J_density_map_OW_vs_NW <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_OW_vs_NW.rds"))


# Load the custom STRAPP test function
source("./functions/run_STRAPP_test.R")

# Extract time scale from LTT data
# LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/LTT_per_clades_melted_df.rds")
LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/LTT_per_clades_melted_df.rds")
time_scale <- unique(LTT_per_clades_melted_df$time)
time_scale <- time_scale[order(time_scale, decreasing = T)]

## May need to unload 'raster' package to avoid conflict
detach("package:raster", unload = TRUE)

## Set seed for reproductibility of permutation tests
set.seed(seed = 12345)

## Loop per time scale
STRAPP_data_OW_vs_NW_all_time <- list()
STRAPP_tests_OW_vs_NW_all_time <- list()
for (i in seq_along(time_scale))
# for (i in 10:length(time_scale))
{
  # i <- 1
  
  # Extract focal time
  time_i <- time_scale[i]
  
  ## 12.3.1/ Get STRAPP data ####
  
  # Get tip regimes and tip rates for focal time i
  BAMM_data_updated_i <- update_tipStates_and_tipRates_for_focal_time(BAMM_object = BAMM_posterior_samples_data, time = time_i, update_rates = T, verbose = T)
  
  # Extract only $tipStates, $tipLambda and $tipMu in the loop per time to save place
  BAMM_data_updated_i <- BAMM_data_updated_i[c("tipStates", "tipLambda", "tipMu", "type", "tip.label")]
  
  # Get tip states for focal time i
  bioregion_data_OW_vs_NW_i <- get_most_likely_binary_states_for_focal_time(density_map = DEC_J_density_map_OW_vs_NW, time = time_i)
  
  ## Store input data for STRAPP test
  STRAPP_data_OW_vs_NW_all_time[[i]] <- list(BAMM_data = BAMM_data_updated_i,
                                             Bioregion_data_OW_vs_NW = bioregion_data_OW_vs_NW_i)
  
  # Save input data for STRAPP test
  # saveRDS(STRAPP_data_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_data_OW_vs_NW_all_time.rds")
  saveRDS(STRAPP_data_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_data_OW_vs_NW_all_time.rds")

  ## 12.3.2/ Run STRAPP test ####
  
  if ((is.null(bioregion_data_OW_vs_NW_i) | (any(!c("Old World", "New World") %in% bioregion_data_OW_vs_NW_i))))
  { # If not edge data, provide NA
    
    STRAPP_test_OW_vs_NW_i <- list()
    STRAPP_test_OW_vs_NW_i$estimate <- NA
    STRAPP_test_OW_vs_NW_i$p.value <- NA
    STRAPP_test_OW_vs_NW_i$obs.corr <- NA
    STRAPP_test_OW_vs_NW_i$gen <- NA
    STRAPP_test_OW_vs_NW_i$null <- NA
    STRAPP_test_OW_vs_NW_i$test_stats <- NA
    STRAPP_test_OW_vs_NW_i$stat_median <- NA
    STRAPP_test_OW_vs_NW_i$stat_Q5 <- NA
    
  } else {
    
    # Run STRAPP test if edge data is present
    STRAPP_test_OW_vs_NW_i <- run_STRAPP_test(BAMM_data = BAMM_data_updated_i,
                                              trait_data = bioregion_data_OW_vs_NW_i,
                                              reps = 1000, rate = "net diversification",
                                              return.full = T,
                                              method = "mann-whitney", # For categorical binomial data (G = 2)
                                              logrates = F, # Do not use log for Mann-Whitney as it is a rank test, so there is no need, and it prevents removing the negative rates
                                              two.tailed = T,
                                              replace = FALSE,
                                              nthreads = 1)
    
    # # Explore output
    # STRAPP_test_OW_vs_NW_i$estimate # Mean tip rates per categories
    # STRAPP_test_OW_vs_NW_i$p.value
    # STRAPP_test_OW_vs_NW_i$obs.corr # Observed statistic for each posterior sample
    # STRAPP_test_OW_vs_NW_i$gen # Generation ID of the selected posterior sample
    # STRAPP_test_OW_vs_NW_i$null # Null statistic for each posterior sample
    
    # Compute info for histogram 
    STRAPP_test_OW_vs_NW_i$test_stats <- STRAPP_test_OW_vs_NW_i$null - STRAPP_test_OW_vs_NW_i$obs.corr
    # summary(STRAPP_test_OW_vs_NW_i$test_stats)
    # table(STRAPP_test_OW_vs_NW_i$test_stats > 0)
    STRAPP_test_OW_vs_NW_i$stat_median <- median(STRAPP_test_OW_vs_NW_i$test_stats)
    STRAPP_test_OW_vs_NW_i$stat_Q5 <- quantile(STRAPP_test_OW_vs_NW_i$test_stats, p = 0.05)
  }
  
  # Store focal time information
  STRAPP_test_OW_vs_NW_i$focal_time <- time_i
  
  ## Store test output
  STRAPP_tests_OW_vs_NW_all_time[[i]] <- STRAPP_test_OW_vs_NW_i
  
  # Save STRAPP test outputs
  # saveRDS(STRAPP_tests_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")
  saveRDS(STRAPP_tests_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")
  
  ## Print progress
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - STRAPP test ran for time = ",time_i," My - n°", i, "/", length(time_scale),"\n"))
  }
}


## 12.3.3/ Plot evolution of p-value in time ####

# Load STRAPP test outputs
# STRAPP_tests_OW_vs_NW_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")
STRAPP_tests_OW_vs_NW_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")

# Extract p-value and time
STRAPP_tests_OW_vs_NW_df <- data.frame(time = unlist(lapply(STRAPP_tests_OW_vs_NW_all_time, FUN = function (x) { x$focal_time }) ),
                                       p_value = unlist(lapply(STRAPP_tests_OW_vs_NW_all_time, FUN = function (x) { x$p.value } )))
STRAPP_tests_OW_vs_NW_df

# Save STRAPP tests through evolutionary time df for OW vs NW
# saveRDS(STRAPP_tests_OW_vs_NW_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df.rds")
saveRDS(STRAPP_tests_OW_vs_NW_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df.rds")

# Extract root age
root_age <- time_scale[1]
# Set X-axis limit
max_age <- 85

# Prepare data for significance area
significance_area_poly_df <- data.frame(p_value = c(0.00, 0.00, 0.05, 0.05),
                                        # time = c(root_age, 0, 0, root_age), 
                                        time = c(max_age, 0, 0, max_age),
                                        poly_ID = rep(1, 4))

# Prepare data for significance gradient
  # One polygon per Y/p-value
p_value_y <- 0.05
p_value_y_data <- c()
nb_poly <- 0
max_gradient <- 0.10
y_increment <- 0.001
while (p_value_y < max_gradient)
{
  p_value_y_data_i <- c(p_value_y, p_value_y, p_value_y + y_increment, p_value_y + y_increment)
  p_value_y_data <- c(p_value_y_data, p_value_y_data_i)
  
  p_value_y <- p_value_y + y_increment
  nb_poly <- nb_poly + 1
}

significance_area_gradient_df <- data.frame(p_value = p_value_y_data,
                                            # time = rep(x = c(root_age, 0, 0, root_age), times = nb_poly), 
                                            time = rep(x = c(max_age, 0, 0, max_age), times = nb_poly), 
                                            poly_ID = rep(1:nb_poly, each = 4))

# Identify time window of significance (p-value <= 0.10)
signif_test <- STRAPP_tests_OW_vs_NW_df$p_value <= 0.10
RLE_output <- rle(signif_test)
ends <- cumsum(RLE_output$lengths)
starts <- c(1, ends + 1)[-(length(RLE_output$lengths) + 1)]
start_signif <- starts[replace_na(data = RLE_output$values, replace = F)]
end_signif <- ends[replace_na(data = RLE_output$values, replace = F)]
time_signif <- STRAPP_tests_OW_vs_NW_df$time[unique(c(start_signif, end_signif))]


## GGplot
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_pvalTT.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_pvalTT.pdf", height = 6, width = 8)

STRAPP_tests_OW_vs_NW_plot <- ggplot(data = STRAPP_tests_OW_vs_NW_df,
                                     mapping = aes(y = p_value, x = time)) +
  
  # Plot significance area above p-value < 0.05
  geom_polygon(data = significance_area_poly_df,
               mapping = aes(y = p_value, x = time, group = poly_ID),
               fill = "limegreen", col = NA,
               alpha = 1.00,
               linewidth = 1.0) +
  
  # Plot significance gradient (0.05 < p-value < 0.10)
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = p_value, x = time,
                             group = poly_ID, alpha = poly_ID),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = F) +
  
  # Plot vertical lines of significance
  geom_vline(xintercept = time_signif,
             col = "black", lty = 2,
             linewidth = 1.0) +
  
  # # Add significance legend
  # annotate(geom = "text", x = 21.5, y = 0.5,
  #          label = "Significant\nperiod",
  #          size = 5.0, hjust = 0.5,
  #          fontface = "bold",
  #          color = "limegreen") +
  
  # Plot mean lines
  geom_line(col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add special times to investigate differences
  geom_vline(xintercept = 23, linewidth = 1.0, linetype = "dashed") + # T = 23 My (Oligocene-Miocene) => Peak of differences
  geom_vline(xintercept = 10, linewidth = 1.0, linetype = "dashed") + # T = 10 My (Mid Miocene) => Inflexion/Crossing point
  geom_vline(xintercept = 0, linewidth = 1.0, linetype = "dashed") + # T = 0 My (Present) => Classic BAMM test
  
  # Set plot title +
  ggtitle(label = paste0("STRAPP test\nDifference in net div. rates through time\nOld World vs. New World")) +

  # Set axes labels
  xlab("Time  [My]") +
  ylab("P-value") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     expand = c(0,0),
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Reverse p-value scale
  scale_y_continuous(transform = "reverse",
                     limits = c(1, 0) # Set limits
  ) +
  
  # Allow drawing outside plot inner margins
  coord_cartesian(clip = "off") +
  
  # Adjust alpha scale
  scale_alpha_continuous(range = c(1,0)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(STRAPP_tests_OW_vs_NW_plot)

dev.off()


## 12.3.4/ Add significance period to the RatesTT plot ####

# Function to use annotate in ggplot with npc units
annotate_npc <- function(label, x, y, ...)
{
  ggplot2::annotation_custom(grid::textGrob(
    x = unit(x, "npc"), y = unit(y, "npc"), label = label, ...))
}

# Load melted dataframe or Rates through time
# weighted_mean_per_bioregions_all_BAMM_phylo_melted <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")
# weighted_mean_per_bioregions_median <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_median.rds")
# weighted_mean_per_bioregions_quantiles <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/weighted_mean_per_bioregions_quantiles.rds")
weighted_mean_per_bioregions_all_BAMM_phylo_melted <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")
weighted_mean_per_bioregions_median <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_median.rds")
weighted_mean_per_bioregions_quantiles <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/weighted_mean_per_bioregions_quantiles.rds")

# Load STRAPP tests through evolutionary time df for OW vs NW
# STRAPP_tests_OW_vs_NW_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df.rds")
STRAPP_tests_OW_vs_NW_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df.rds")

## Set OW vs. NW color scheme
colors_for_OW_NW <- c("mediumpurple2", "peachpuff2")
OW_NW_factors <- c("OW", "NW")
OW_NW_names <- c("Old World", "New World")
names(colors_for_OW_NW) <- OW_NW_names

# Filter df to extract bioregion rates (and not aggregated rates)
weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  filter(bioregion %in% OW_NW_factors)
weighted_mean_per_bioregions_median_OW_vs_NW <- weighted_mean_per_bioregions_median %>% 
  filter(bioregion %in% OW_NW_factors)
weighted_mean_per_bioregions_quantiles_OW_vs_NW <- weighted_mean_per_bioregions_quantiles %>% 
  filter(bioregion %in% OW_NW_factors)

# Filter time to avoid failed polygons when plotting
weighted_mean_per_bioregions_quantiles_OW_vs_NW <- weighted_mean_per_bioregions_quantiles_OW_vs_NW %>% 
  filter(time <= 85)

## Interpolate p-values to thinner time scale
x_increment <- 0.1
STRAPP_tests_OW_vs_NW_df_thin <- data.frame()
for (i in 1:(nrow(STRAPP_tests_OW_vs_NW_df)-1))
{
  # i <- 90
  
  # Extract coarse time boundaries
  start_time_i <- STRAPP_tests_OW_vs_NW_df$time[i]
  end_time_i <- STRAPP_tests_OW_vs_NW_df$time[i+1]
  
  # Extract coarse time p-values
  start_pval_i <- STRAPP_tests_OW_vs_NW_df$p_value[i]
  end_pval_i <- STRAPP_tests_OW_vs_NW_df$p_value[i+1]
  
  # Generate thinner time steps
  time_scale_thinner_i <- seq(from = start_time_i,
                              to = end_time_i,
                              by = -x_increment)
  # If any NA, fill with NA
  if (any(is.na(c(start_pval_i, end_pval_i))))
  {
    pval_thinner_i <- rep(NA, length(time_scale_thinner_i))
  } else {
    # Interpolate p-values
    pval_increment_i <- (end_pval_i - start_pval_i) / (length(time_scale_thinner_i) - 1)
    pval_thinner_i <- seq(from = start_pval_i,
                          to = end_pval_i,
                          by = pval_increment_i)
  }

  # Store results
  STRAPP_tests_OW_vs_NW_df_thin_i <- data.frame(time = time_scale_thinner_i,
                                                p_value = pval_thinner_i)
  
  STRAPP_tests_OW_vs_NW_df_thin <- rbind(STRAPP_tests_OW_vs_NW_df_thin, STRAPP_tests_OW_vs_NW_df_thin_i)
}
# Remove duplicates
STRAPP_tests_OW_vs_NW_df_thin <- STRAPP_tests_OW_vs_NW_df_thin %>% 
  dplyr::distinct(time, .keep_all = TRUE)

# Save interpolated time/p-value for STRAPP tests: OW vs NW
# saveRDS(STRAPP_tests_OW_vs_NW_df_thin, "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df_thin.rds")
saveRDS(STRAPP_tests_OW_vs_NW_df_thin, "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_df_thin.rds")


## Prepare data for significance gradient
# One polygon per x-time
nb_poly <- 0
time_x_data <- c()
max_rates <- max(weighted_mean_per_bioregions_quantiles_OW_vs_NW$quant_rates, na.rm = T)
min_rates <- min(weighted_mean_per_bioregions_quantiles_OW_vs_NW$quant_rates, na.rm = T)

# Loop per time
for (i in 1:(nrow(STRAPP_tests_OW_vs_NW_df_thin-1)))
{
  time_x <- STRAPP_tests_OW_vs_NW_df_thin$time[i]
  time_x2 <- STRAPP_tests_OW_vs_NW_df_thin$time[i+1]
  
  time_x_data_i <- c(time_x, time_x2, time_x2, time_x)
  time_x_data <- c(time_x_data, time_x_data_i)
  
  nb_poly <- nb_poly + 1
}

significance_area_gradient_df <- data.frame(time = time_x_data,
                                            rates = rep(x = c(max_rates, max_rates, min_rates, min_rates), times = nb_poly), 
                                            poly_ID = rep(1:nb_poly, each = 4))
# Join p-values
significance_area_gradient_df <- left_join(x = significance_area_gradient_df,
                                           y = STRAPP_tests_OW_vs_NW_df_thin)

# Remove polygons with p-value > 0.1
# Fix values of polygons with p-value < 0.05
min_gradient <- 0.05
max_gradient <- 0.10
significance_area_gradient_df <- significance_area_gradient_df %>% 
  filter(!(p_value > max_gradient))
significance_area_gradient_df$p_value[significance_area_gradient_df$p_value < min_gradient] <- min_gradient
  
## Select times of interest to highlight BAMM test results
STRAPP_tests_OW_vs_NW_df$time[which.min(STRAPP_tests_OW_vs_NW_df$p_value)] # Time with most significant difference
STRAPP_tests_OW_vs_NW_df$time[STRAPP_tests_OW_vs_NW_df$time < 23][which.max(STRAPP_tests_OW_vs_NW_df$p_value[STRAPP_tests_OW_vs_NW_df$time < 23])] # Time with least differences during Miocene

# T = 0 My (Present)
# T = 10 My (Tests available for Trop vs. Temp)
# T = 15 My (Mid-Miocene)
# T = 23 My (Oligocene-Miocene)
focal_times <- c(0, 10, 23)
focal_times_colors <- c("brown", "red", "sienna2")
focal_times_colors <- c("grey20", "grey50", "grey80")
  
# Get closest test time and ID to the selected focal time
focal_times_test_OW_vs_NW_df <- expand.grid(focal_times, time_scale)
names(focal_times_test_OW_vs_NW_df) <- c("focal_time", "test_time")
focal_times_test_OW_vs_NW_df$abs_diff <- abs(focal_times_test_OW_vs_NW_df$focal_time - focal_times_test_OW_vs_NW_df$test_time)
focal_times_test_OW_vs_NW_df <- focal_times_test_OW_vs_NW_df %>% 
  arrange(focal_time) %>% 
  group_by(focal_time) %>% 
  mutate(diff_ranks = rank(abs_diff)) %>%
  arrange(focal_time, diff_ranks) %>%
  filter(diff_ranks == 1) %>%
  select(-diff_ranks) %>%
  mutate(test_time_ID = which(test_time == time_scale)) %>%
  ungroup()

# Extract test results for those focal times
# Test is comparing Mann-Whitney's U test stat between observed and randomize data for 1000 of BAMM posteriors
# Thus, the null distribution is made of the 1000 delta_U stats
# Test is significant if 95% of posteriors shows higher correlation than at random
# Thus, Q5% of the null distribution needs to be positive
focal_times_test_results <- lapply(X = STRAPP_tests_OW_vs_NW_all_time[focal_times_test_OW_vs_NW_df$test_time_ID],
                                   FUN = function (x) { test_time <- x$focal_time ; p_value <- x$p.value ; Q5 <- unname(x$stat_Q5) ; return(c(test_time = test_time, p_value = p_value, Q5 = Q5)) })
focal_times_test_results <- data.frame(do.call(rbind, focal_times_test_results))
focal_times_test_OW_vs_NW_df <- left_join(focal_times_test_OW_vs_NW_df, focal_times_test_results)
focal_times_test_OW_vs_NW_df$color <- focal_times_colors
focal_times_test_OW_vs_NW_df$y_coords <- c(0.14, 0.13, 0.12)
  

## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines_with_signif.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines_with_signif.pdf", height = 8, width = 12)

ratesT_OW_vs_NW_ggplot_lines_with_signif <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW) +
  
  # Plot significance gradient
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = rates, x = time,
                             group = poly_ID, alpha = p_value),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = T) +
  
  # Plot 1000 line replicates
  geom_line(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW,
            mapping = aes(y = weighted_mean, x = time, group = interaction(bioregion, sample), col = bioregion),
            alpha = 0.01,
            linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_OW_vs_NW,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Prevent rate Y-scale to expand
  scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.155)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust alpha scale
  scale_alpha_continuous("P-values", range = c(0.8,0), breaks = c(0.05, 0.075, 0.09), labels = c("< 0.05", "< 0.075", "< 0.1")) +
  
  # Remove fill legend
  guides(fill = "none",
         color = guide_legend(order = 1,
                              override.aes = list(fill = NA))) +
  
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Add significance tests for focal times
for (i in 1:nrow(focal_times_test_OW_vs_NW_df))
{
  # i <- 1
  
  ratesT_OW_vs_NW_ggplot_lines_with_signif <- ratesT_OW_vs_NW_ggplot_lines_with_signif +
    
    # Add vertical line
    geom_vline(xintercept = focal_times_test_OW_vs_NW_df$focal_time[i],
               linewidth = 1.0, linetype = "dotted",
               # color = focal_times_test_OW_vs_NW_df$color[i],
               color = "black") +
    
    # Add test results in upper right corner
    annotate(geom = "text", x = 84, y = focal_times_test_OW_vs_NW_df$y_coords[1] + 0.01,
             hjust = 0, fontface = "bold", size = 5,
             label = "STRAPP tests") +
    annotate(geom = "text", x = 84,
             y = focal_times_test_OW_vs_NW_df$y_coords[i],
             hjust = 0, fontface = "plain", size = 5,
             # color = focal_times_test_OW_vs_NW_df$color[i],
             color = "black",
             label = paste0("T = ",focal_times_test_OW_vs_NW_df$focal_time[i]," My: ",
                            "Q5 = ", round(focal_times_test_OW_vs_NW_df$Q5[i]), ", ",
                            "p = ",focal_times_test_OW_vs_NW_df$p_value[i]))
}

# Plot
print(ratesT_OW_vs_NW_ggplot_lines_with_signif)

dev.off()

## Add symbols to the vertical line and STRAPP test legend


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_fuzzy_quantiles_with_signif.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_OW_vs_NW_fuzzy_quantiles_with_signif.pdf", height = 8, width = 12)

ratesT_OW_vs_NW_ggplot_lines_with_signif <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_OW_vs_NW) +
  
  # Plot significance gradient
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = rates, x = time,
                             group = poly_ID, alpha = p_value),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = T) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = weighted_mean_per_bioregions_quantiles_OW_vs_NW,
               mapping = aes(y = quant_rates, x = time, group = interaction(bioregion, quantile_ID), fill = bioregion),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_OW_vs_NW,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Prevent rate Y-scale to expand
  # scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.155)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust alpha scale
  scale_alpha_continuous("P-values", range = c(0.8,0), breaks = c(0.05, 0.075, 0.09), labels = c("< 0.05", "< 0.075", "< 0.1")) +
  
  # Remove fill legend
  guides(fill = "none",
         color = guide_legend(order = 1,
                              override.aes = list(fill = NA))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Add significance tests for focal times
for (i in 1:nrow(focal_times_test_OW_vs_NW_df))
{
  # i <- 1
  
  ratesT_OW_vs_NW_ggplot_lines_with_signif <- ratesT_OW_vs_NW_ggplot_lines_with_signif +
    
    # Add vertical line
    geom_vline(xintercept = focal_times_test_OW_vs_NW_df$focal_time[i],
               linewidth = 1.0, linetype = "dotted",
               # color = focal_times_test_OW_vs_NW_df$color[i]
               color = "black"
               ) +
    
    # Add test results in upper right corner
    annotate(geom = "text", x = 84, y = focal_times_test_OW_vs_NW_df$y_coords[1] + 0.01,
             hjust = 0, fontface = "bold", size = 5,
             label = "STRAPP tests") +
    annotate(geom = "text", x = 84,
             y = focal_times_test_OW_vs_NW_df$y_coords[i],
             hjust = 0, fontface = "plain", size = 5,
             # color = focal_times_test_OW_vs_NW_df$color[i],
             color = "black",
             label = paste0("T = ",focal_times_test_OW_vs_NW_df$focal_time[i]," My: ",
                            "Q5 = ", round(focal_times_test_OW_vs_NW_df$Q5[i]), ", ",
                            "p = ",focal_times_test_OW_vs_NW_df$p_value[i]))
}


# Plot
print(ratesT_OW_vs_NW_ggplot_lines_with_signif)

dev.off()

## Add symbols to the vertical line and STRAPP test legend


### 12.4/ STRAPP test along evolutionary time for Trop_vs_Temp #### 

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

## Load the binary density map for Tropics vs Temperate
# DEC_J_density_map_Trop_vs_Temp <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_rough_phylogeny_1534t/DEC_J_density_map_Trop_vs_Temp.rds"))
DEC_J_density_map_Trop_vs_Temp <- readRDS(file = paste0("./outputs/Density_maps/Ponerinae_MCC_phylogeny_1534t/DEC_J_density_map_Trop_vs_Temp.rds"))

# Load the custom STRAPP test function
source("./functions/run_STRAPP_test.R")

# Extract time scale from LTT data
# LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/LTT_per_clades_melted_df.rds")
LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/LTT_per_clades_melted_df.rds")
time_scale <- unique(LTT_per_clades_melted_df$time)
time_scale <- time_scale[order(time_scale, decreasing = T)]

## May need to unload 'raster' package to avoid conflict
detach("package:raster", unload = TRUE)

## Set seed for reproducibility of permutation tests
set.seed(seed = 12345)

## Loop per time scale
STRAPP_data_Trop_vs_Temp_all_time <- list()
STRAPP_tests_Trop_vs_Temp_all_time <- list()
for (i in seq_along(time_scale))
# for (i in 6:length(time_scale))
{
  # i <- 1
  
  # Extract focal time
  time_i <- time_scale[i]
  
  ## 12.4.1/ Get STRAPP data ####
  
  # Get tip regimes and tip rates for focal time i
  BAMM_data_updated_i <- update_tipStates_and_tipRates_for_focal_time(BAMM_object = BAMM_posterior_samples_data, time = time_i, update_rates = T, verbose = T)
  
  # Extract only $tipStates, $tipLambda and $tipMu in the loop per time to save place
  BAMM_data_updated_i <- BAMM_data_updated_i[c("tipStates", "tipLambda", "tipMu", "type", "tip.label")]
  
  # Get tip states for focal time i
  bioregion_data_Trop_vs_Temp_i <- get_most_likely_binary_states_for_focal_time(density_map = DEC_J_density_map_Trop_vs_Temp, time = time_i)
  
  ## Store input data for STRAPP test
  STRAPP_data_Trop_vs_Temp_all_time[[i]] <- list(BAMM_data = BAMM_data_updated_i,
                                             Bioregion_data_Trop_vs_Temp = bioregion_data_Trop_vs_Temp_i)
  
  # Save input data for STRAPP test
  # saveRDS(STRAPP_data_Trop_vs_Temp_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_data_Trop_vs_Temp_all_time.rds")
  saveRDS(STRAPP_data_Trop_vs_Temp_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_data_Trop_vs_Temp_all_time.rds")
  
  ## 12.4.2/ Run STRAPP test ####
  
  if ((is.null(bioregion_data_Trop_vs_Temp_i) | (any(!c("Temperate", "Tropics") %in% bioregion_data_Trop_vs_Temp_i))))
  { # If not edge data, or only one level, provide NA
    
    STRAPP_test_Trop_vs_Temp_i <- list()
    STRAPP_test_Trop_vs_Temp_i$estimate <- NA
    STRAPP_test_Trop_vs_Temp_i$p.value <- NA
    STRAPP_test_Trop_vs_Temp_i$obs.corr <- NA
    STRAPP_test_Trop_vs_Temp_i$gen <- NA
    STRAPP_test_Trop_vs_Temp_i$null <- NA
    STRAPP_test_Trop_vs_Temp_i$test_stats <- NA
    STRAPP_test_Trop_vs_Temp_i$stat_median <- NA
    STRAPP_test_Trop_vs_Temp_i$stat_Q5 <- NA
    
  } else {
    
    # Run STRAPP test if edge data is present
    STRAPP_test_Trop_vs_Temp_i <- run_STRAPP_test(BAMM_data = BAMM_data_updated_i,
                                              trait_data = bioregion_data_Trop_vs_Temp_i,
                                              reps = 1000, rate = "net diversification",
                                              return.full = T,
                                              method = "mann-whitney", # For categorical binomial data (G = 2)
                                              logrates = F, # Do not use log for Mann-Whitney as it is a rank test, so there is no need, and it prevents removing the negative rates
                                              two.tailed = T,
                                              replace = FALSE,
                                              nthreads = 1)
    
    # # Explore output
    # STRAPP_test_Trop_vs_Temp_i$estimate # Mean tip rates per categories
    # STRAPP_test_Trop_vs_Temp_i$p.value
    # STRAPP_test_Trop_vs_Temp_i$obs.corr # Observed statistic for each posterior sample
    # STRAPP_test_Trop_vs_Temp_i$gen # Generation ID of the selected posterior sample
    # STRAPP_test_Trop_vs_Temp_i$null # Null statistic for each posterior sample
    
    # Compute info for histogram 
    STRAPP_test_Trop_vs_Temp_i$test_stats <- STRAPP_test_Trop_vs_Temp_i$null - STRAPP_test_Trop_vs_Temp_i$obs.corr
    # summary(STRAPP_test_Trop_vs_Temp_i$test_stats)
    # table(STRAPP_test_Trop_vs_Temp_i$test_stats > 0)
    STRAPP_test_Trop_vs_Temp_i$stat_median <- median(STRAPP_test_Trop_vs_Temp_i$test_stats)
    STRAPP_test_Trop_vs_Temp_i$stat_Q5 <- quantile(STRAPP_test_Trop_vs_Temp_i$test_stats, p = 0.05)
  }
  
  # Store focal time information
  STRAPP_test_Trop_vs_Temp_i$focal_time <- time_i
  
  ## Store test output
  STRAPP_tests_Trop_vs_Temp_all_time[[i]] <- STRAPP_test_Trop_vs_Temp_i
  
  # Save STRAPP test outputs
  # saveRDS(STRAPP_tests_Trop_vs_Temp_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_all_time.rds")
  saveRDS(STRAPP_tests_Trop_vs_Temp_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_all_time.rds")
  
  ## Print progress
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - STRAPP test ran for time = ",time_i," My - n°", i, "/", length(time_scale),"\n"))
  }
}


## 12.4.3/ Plot evolution of p-value in time ####

# # Load input data for STRAPP test
# STRAPP_data_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_data_Trop_vs_Temp_all_time.rds")
# STRAPP_data_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_data_Trop_vs_Temp_all_time.rds")

# Load STRAPP test outputs
# STRAPP_tests_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_all_time.rds")
STRAPP_tests_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_all_time.rds")

# Extract p-value and time
STRAPP_tests_Trop_vs_Temp_df <- data.frame(time = unlist(lapply(STRAPP_tests_Trop_vs_Temp_all_time, FUN = function (x) { x$focal_time }) ),
                                       p_value = unlist(lapply(STRAPP_tests_Trop_vs_Temp_all_time, FUN = function (x) { x$p.value } )))
STRAPP_tests_Trop_vs_Temp_df

# Save STRAPP tests through evolutionary time df for OW vs NW
# saveRDS(STRAPP_tests_Trop_vs_Temp_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_df.rds")
saveRDS(STRAPP_tests_Trop_vs_Temp_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_df.rds")

# Extract root age
root_age <- time_scale[1]
# Set plot limits to Miocene
max_age <- 23

# Prepare data for significance area
significance_area_poly_df <- data.frame(p_value = c(0.00, 0.00, 0.05, 0.05),
                                        # time = c(root_age, 0, 0, root_age), 
                                        time = c(max_age, 0, 0, max_age), 
                                        poly_ID = rep(1, 4))

# Prepare data for significance gradient
# One polygon per Y/p-value
p_value_y <- 0.05
p_value_y_data <- c()
nb_poly <- 0
max_gradient <- 0.10
y_increment <- 0.001
while (p_value_y < max_gradient)
{
  p_value_y_data_i <- c(p_value_y, p_value_y, p_value_y + y_increment, p_value_y + y_increment)
  p_value_y_data <- c(p_value_y_data, p_value_y_data_i)
  
  p_value_y <- p_value_y + y_increment
  nb_poly <- nb_poly + 1
}

significance_area_gradient_df <- data.frame(p_value = p_value_y_data,
                                            # time = rep(x = c(root_age, 0, 0, root_age), times = nb_poly),
                                            time = rep(x = c(max_age, 0, 0, max_age), times = nb_poly),
                                            poly_ID = rep(1:nb_poly, each = 4))

# No period of significance! Min p-value = 0.221 at 2.5 My
STRAPP_tests_Trop_vs_Temp_df[which.min(STRAPP_tests_Trop_vs_Temp_df$p_value), ]


## GGplot
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_pvalTT.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_pvalTT.pdf", height = 6, width = 8)

STRAPP_tests_Trop_vs_Temp_plot <- ggplot(# data = STRAPP_tests_Trop_vs_Temp_df,
                                         data = STRAPP_tests_Trop_vs_Temp_df[!is.na(STRAPP_tests_Trop_vs_Temp_df$p_value), ],
                                         mapping = aes(y = p_value, x = time)) +
  
  # Plot significance area above p-value < 0.05
  geom_polygon(data = significance_area_poly_df,
               mapping = aes(y = p_value, x = time, group = poly_ID),
               fill = "limegreen", col = NA,
               alpha = 1.00,
               linewidth = 1.0) +
  
  # Plot significance gradient (0.05 < p-value < 0.10)
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = p_value, x = time,
                             group = poly_ID, alpha = poly_ID),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = F) +
  
  # Plot mean lines
  geom_line(col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add special times to investigate differences
  geom_vline(xintercept = 10, linewidth = 1.0, linetype = "dashed") + # T = 10 My (Mid Miocene) => Biggest difference
  geom_vline(xintercept = 0, linewidth = 1.0, linetype = "dashed") + # T = 0 My (Present) => Classic BAMM test
  
  # Set plot title +
  ggtitle(label = paste0("STRAPP test\nDifference in net div. rates through time\nTropics vs. Temperate")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("P-value") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     expand = c(0,0),
                     limits = c(23, 0) # Set limits to Miocene
  ) + 
  
  # Reverse p-value scale
  scale_y_continuous(transform = "reverse",
                     limits = c(1, 0) # Set limits
  ) +

  # Allow drawing outside plot inner margins
  coord_cartesian(clip = "off") +
  
  # Adjust alpha scale
  scale_alpha_continuous(range = c(1,0)) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Plot
print(STRAPP_tests_Trop_vs_Temp_plot)

dev.off()

## 12.4.4/ Add significance period to the RatesTT plot ####

# Not any significant period!

# Only add tests for focal times

# Load STRAPP tests through evolutionary time df for Trop vs Temp
# STRAPP_tests_Trop_vs_Temp_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_df.rds")
STRAPP_tests_Trop_vs_Temp_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_Trop_vs_Temp_df.rds")

# T = 0 My (Present)
# T = 10 My (Tests available for Trop vs. Temp)
focal_times <- c(0, 10)
focal_times_colors <- c("grey20", "grey50")

# Get closest test time and ID to the selected focal time
focal_times_test_Trop_vs_Temp_df <- expand.grid(focal_times, time_scale)
names(focal_times_test_Trop_vs_Temp_df) <- c("focal_time", "test_time")
focal_times_test_Trop_vs_Temp_df$abs_diff <- abs(focal_times_test_Trop_vs_Temp_df$focal_time - focal_times_test_Trop_vs_Temp_df$test_time)
focal_times_test_Trop_vs_Temp_df <- focal_times_test_Trop_vs_Temp_df %>% 
  arrange(focal_time) %>% 
  group_by(focal_time) %>% 
  mutate(diff_ranks = rank(abs_diff)) %>%
  arrange(focal_time, diff_ranks) %>%
  filter(diff_ranks == 1) %>%
  select(-diff_ranks) %>%
  mutate(test_time_ID = which(test_time == time_scale)) %>%
  ungroup()

# Extract test results for those focal times
# Test is comparing Mann-Whitney's U test stat between observed and randomize data for 1000 of BAMM posteriors
# Thus, the null distribution is made of the 1000 delta_U stats
# Test is significant if 95% of posteriors shows higher correlation than at random
# Thus, Q5% of the null distribution needs to be positive
focal_times_test_results <- lapply(X = STRAPP_tests_Trop_vs_Temp_all_time[focal_times_test_Trop_vs_Temp_df$test_time_ID],
                                   FUN = function (x) { test_time <- x$focal_time ; p_value <- x$p.value ; Q5 <- unname(x$stat_Q5) ; return(c(test_time = test_time, p_value = p_value, Q5 = Q5)) })
focal_times_test_results <- data.frame(do.call(rbind, focal_times_test_results))
focal_times_test_Trop_vs_Temp_df <- left_join(focal_times_test_Trop_vs_Temp_df, focal_times_test_results)
focal_times_test_Trop_vs_Temp_df$color <- focal_times_colors
focal_times_test_Trop_vs_Temp_df$y_coords <- c(0.14, 0.13)


## Generate plot using lines
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_ggplot_fuzzy_lines_with_signif.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_ggplot_fuzzy_lines_with_signif.pdf", height = 8, width = 12)

ratesT_Trop_vs_Temp_ggplot_lines_with_signif <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp) +
  
  # # Plot significance gradient
  # geom_polygon(data = significance_area_gradient_df,
  #              mapping = aes(y = rates, x = time,
  #                            group = poly_ID, alpha = p_value),
  #              fill = "limegreen", col = NA,
  #              linewidth = 1.0, show.legend = T) +
  
  # Plot 1000 line replicates
  geom_line(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp,
            mapping = aes(y = weighted_mean, x = time, group = interaction(bioregion, sample), col = bioregion),
            alpha = 0.01,
            linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_Trop_vs_Temp,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Prevent rate Y-scale to expand
  scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.155)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = Trop_Temp_names, labels = Trop_Temp_names, values = colors_for_Trop_Temp) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = Trop_Temp_names, labels = Trop_Temp_names, values = colors_for_Trop_Temp) +
  
  # # Adjust alpha scale
  # scale_alpha_continuous("P-values", range = c(0.8,0), breaks = c(0.05, 0.075, 0.09), labels = c("< 0.05", "< 0.075", "< 0.1")) +
  
  # # Remove fill legend
  # guides(fill = "none",
  #        color = guide_legend(order = 1,
  #                             override.aes = list(fill = NA))) +
  
  # Remove fill legend
  guides(fill = "none") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Add significance tests for focal times
for (i in 1:nrow(focal_times_test_Trop_vs_Temp_df))
{
  # i <- 1
  
  ratesT_Trop_vs_Temp_ggplot_lines_with_signif <- ratesT_Trop_vs_Temp_ggplot_lines_with_signif +
    
    # Add vertical line
    geom_vline(xintercept = focal_times_test_Trop_vs_Temp_df$focal_time[i],
               linewidth = 1.0, linetype = "dotted",
               # color = focal_times_test_Trop_vs_Temp_df$color[i],
               color = "black") +
    
    # Add test results in upper right corner
    annotate(geom = "text", x = 84, y = focal_times_test_Trop_vs_Temp_df$y_coords[1] + 0.01,
             hjust = 0, fontface = "bold", size = 5,
             label = "STRAPP tests") +
    annotate(geom = "text", x = 84,
             y = focal_times_test_Trop_vs_Temp_df$y_coords[i],
             hjust = 0, fontface = "plain", size = 5,
             # color = focal_times_test_Trop_vs_Temp_df$color[i],
             color = "black",
             label = paste0("T = ",focal_times_test_Trop_vs_Temp_df$focal_time[i]," My: ",
                            "Q5 = ", round(focal_times_test_Trop_vs_Temp_df$Q5[i]), ", ",
                            "p = ",focal_times_test_Trop_vs_Temp_df$p_value[i]))
}

# Plot
print(ratesT_Trop_vs_Temp_ggplot_lines_with_signif)

dev.off()

## Add symbols to the vertical line and STRAPP test legend


## Generate plot using quantiles
# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_fuzzy_quantiles_with_signif.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_ratesTT_Trop_vs_Temp_fuzzy_quantiles_with_signif.pdf", height = 8, width = 12)

ratesT_Trop_vs_Temp_ggplot_lines_with_signif <- ggplot(data = weighted_mean_per_bioregions_all_BAMM_phylo_melted_Trop_vs_Temp) +
  
  # Plot significance gradient
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = rates, x = time,
                             group = poly_ID, alpha = p_value),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = T) +
  
  # Plot 50 quantile polygons
  geom_polygon(data = weighted_mean_per_bioregions_quantiles_Trop_vs_Temp,
               mapping = aes(y = quant_rates, x = time, group = interaction(bioregion, quantile_ID), fill = bioregion),
               alpha = 0.02,
               linewidth = 5.0) +
  
  # Plot mean lines
  geom_line(data = weighted_mean_per_bioregions_median_Trop_vs_Temp,
            mapping = aes(y = median_rates, x = time, group = bioregion, col = bioregion),
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot div = 0 line
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +
  
  # Set plot title +
  ggtitle(label = paste0("Net diversification rates per Bioregions\nacross 1000 BAMM posterior samples")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("Net diversification rates\n[Events / lineage / My]") +
  
  # Set y-axis limits
  # ylim(c(0, 0.15)) +
  # ylim(c(0, y_max)) +
  
  # Prevent rate Y-scale to expand
  # scale_y_continuous(expand = c(0,0), limits = c(-0.02, 0.155)) +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     # limits = c(root_age, 0) # Set limits
                     limits = c(85, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust color scheme and legend
  scale_color_manual("Bioregions", breaks = OW_NW_factors, labels = OW_NW_names, values = unname(colors_for_OW_NW)) +
  
  # Adjust alpha scale
  scale_alpha_continuous("P-values", range = c(0.8,0), breaks = c(0.05, 0.075, 0.09), labels = c("< 0.05", "< 0.075", "< 0.1")) +
  
  # Remove fill legend
  guides(fill = "none",
         color = guide_legend(order = 1,
                              override.aes = list(fill = NA))) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
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
        axis.ticks.length = unit(8, "pt"),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

# Add significance tests for focal times
for (i in 1:nrow(focal_times_test_Trop_vs_Temp_df))
{
  # i <- 1
  
  ratesT_Trop_vs_Temp_ggplot_lines_with_signif <- ratesT_Trop_vs_Temp_ggplot_lines_with_signif +
    
    # Add vertical line
    geom_vline(xintercept = focal_times_test_Trop_vs_Temp_df$focal_time[i],
               linewidth = 1.0, linetype = "dotted",
               # color = focal_times_test_Trop_vs_Temp_df$color[i]
               color = "black"
               ) +
    
    # Add test results in upper right corner
    annotate(geom = "text", x = 84, y = focal_times_test_Trop_vs_Temp_df$y_coords[1] + 0.01,
             hjust = 0, fontface = "bold", size = 5,
             label = "STRAPP tests") +
    annotate(geom = "text", x = 84,
             y = focal_times_test_Trop_vs_Temp_df$y_coords[i],
             hjust = 0, fontface = "plain", size = 5,
             # color = focal_times_test_Trop_vs_Temp_df$color[i],
             color = "black",
             label = paste0("T = ",focal_times_test_Trop_vs_Temp_df$focal_time[i]," My: ",
                            "Q5 = ", round(focal_times_test_Trop_vs_Temp_df$Q5[i]), ", ",
                            "p = ",focal_times_test_Trop_vs_Temp_df$p_value[i]))
}


# Plot
print(ratesT_Trop_vs_Temp_ggplot_lines_with_signif)

dev.off()

## Add symbols to the vertical line and STRAPP test legend




##### 13/ Map of current net dispersal rates based on tip rates ####

# May need to unload the tidyverse to avoid conflicts...
library(raster)

## Load binary alpha-hull range raster stack
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

## Load bioregion sf
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
Bioregions_sf_Bioregions_level_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregions_sf_Bioregions_level_Mollweide.rds")

## Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

### 13.1/ Compute current median net diversification rates ####

# Speciation rates
lambda_tip_rates <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipLambda)
lambda_tip_rates_median <- apply(X = lambda_tip_rates, MARGIN = 2, FUN = median, na.rm = T)
names(lambda_tip_rates_median) <- BAMM_posterior_samples_data$tip.label

# Extinction rates
mu_tip_rates <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipMu)
mu_tip_rates_median <- apply(X = mu_tip_rates, MARGIN = 2, FUN = median, na.rm = T)
names(mu_tip_rates_median) <- BAMM_posterior_samples_data$tip.label

# Store rates
current_median_tip_rates_df <- data.frame(taxa = BAMM_posterior_samples_data$tip.label,
                                          lambda = lambda_tip_rates_median,
                                          mu = mu_tip_rates_median)

# Compute net diversification rates
current_median_tip_rates_df$net_div <- current_median_tip_rates_df$lambda - current_median_tip_rates_df$mu

# Save current median tip rates
#saveRDS(current_median_tip_rates_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/current_median_tip_rates_df.rds")
saveRDS(current_median_tip_rates_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/current_median_tip_rates_df.rds")


### 13.2/ Compute mean current rates raster ####

# Load terrestrial background
terrestrial_bg_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")

# Load current median tip rates
# current_median_tip_rates_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/current_median_tip_rates_df.rds")
current_median_tip_rates_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/current_median_tip_rates_df.rds")

# Reorder as in raster stack
current_median_tip_rates_df_reordered <- current_median_tip_rates_df[match(names(Ponerinae_alpha_hull_stack_WGS84), table = current_median_tip_rates_df$taxa), ]
cbind(current_median_tip_rates_df_reordered$taxa, names(Ponerinae_alpha_hull_stack_WGS84))

length(current_median_tip_rates_df$taxa)
length(current_median_tip_rates_df_reordered$taxa)
length(names(Ponerinae_alpha_hull_stack_WGS84))

# Assign rates to binary ranges 
Ponerinae_net_div_rates_stack_WGS84 <- Ponerinae_alpha_hull_stack_WGS84 * current_median_tip_rates_df_reordered$net_div
plot(Ponerinae_net_div_rates_stack_WGS84[[1:9]])
rm(Ponerinae_alpha_hull_stack_WGS84) # Remove alpha hull ranges to save RAM space
Ponerinae_net_div_rates_stack_WGS84 <- readAll(Ponerinae_net_div_rates_stack_WGS84)

## Save stack of current net div rates
# saveRDS(Ponerinae_net_div_rates_stack_WGS84, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_net_div_rates_stack_WGS84.rds")
saveRDS(Ponerinae_net_div_rates_stack_WGS84, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_net_div_rates_stack_WGS84.rds")
  

# Compute mean current rates raster (do not account for zeros as they are absences)
Ponerinae_mean_net_div_rates_WGS84 <- raster::calc(x = Ponerinae_net_div_rates_stack_WGS84,
                                                   fun = function (x) { mean(x[x > 0], na.rm = T) })
Ponerinae_mean_net_div_rates_WGS84 <- readAll(Ponerinae_mean_net_div_rates_WGS84)

# Add terrestrial boundaries
temp <- terrestrial_bg_WGS84
temp[!is.na(Ponerinae_mean_net_div_rates_WGS84@data@values)] <- Ponerinae_mean_net_div_rates_WGS84@data@values[!is.na(Ponerinae_mean_net_div_rates_WGS84@data@values)]
Ponerinae_mean_net_div_rates_WGS84 <- temp

plot(Ponerinae_mean_net_div_rates_WGS84)
rm(Ponerinae_net_div_rates_stack_WGS84) # Remove stack of current net div rates to save RAM space

## Save raster of mean current net div rates
# saveRDS(Ponerinae_mean_net_div_rates_WGS84, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_net_div_rates_WGS84.rds")
saveRDS(Ponerinae_mean_net_div_rates_WGS84, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_net_div_rates_WGS84.rds")

### 13.3/ Filter out pixels with low richness ####

## Load raster of mean current net div rates
# Ponerinae_mean_net_div_rates_WGS84 <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_net_div_rates_WGS84.rds")
Ponerinae_mean_net_div_rates_WGS84 <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_net_div_rates_WGS84.rds")

## Load raw richness raster
Ponerinae_species_richness_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/Ponerinae_species_richness_WGS84.rds")

plot(Ponerinae_mean_net_div_rates_WGS84, col = pal_bl_red_brewer)
plot(Ponerinae_species_richness_WGS84, col = pal_bl_red_brewer)

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
Ponerinae_mean_net_div_rates_WGS84_thresholded <- Ponerinae_mean_net_div_rates_WGS84
Ponerinae_mean_net_div_rates_WGS84_thresholded@data@values[Ponerinae_species_richness_WGS84@data@values < min_richness_threshold] <- 0

plot(Ponerinae_mean_net_div_rates_WGS84, col = pal_bl_red_brewer)
plot(Ponerinae_mean_net_div_rates_WGS84_thresholded, col = pal_bl_red_brewer)

pdf(file = paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_net_div_rates_maps_multiple_thresholds.pdf"),
    width = 20, height = 10)

par(mfrow = c(2,2))
plot(Ponerinae_mean_net_div_rates_WGS84, col = pal_bl_red_brewer)
# Repeat with different thresholds
for (i in 1:3)
{
  min_richness_threshold <- i 
  Ponerinae_mean_net_div_rates_WGS84_thresholded <- Ponerinae_mean_net_div_rates_WGS84
  Ponerinae_mean_net_div_rates_WGS84_thresholded@data@values[Ponerinae_species_richness_WGS84@data@values < min_richness_threshold] <- 0
  plot(Ponerinae_mean_net_div_rates_WGS84_thresholded, col = pal_bl_red_brewer, main = paste0("Min theshold = ", min_richness_threshold))
}
par(mfrow = c(1,1))

dev.off()


### 13.4/ Plot mean current rate geographic map ####

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")
pal_bl_red_brewer <- rev(RColorBrewer::brewer.pal(n = 11, "RdYlBu"))
pal_bl_red_brewer_fn <- colorRampPalette(colors = pal_bl_red_brewer)
pal_bl_red_brewer <- pal_bl_red_brewer_fn(n = 400)
# blues <- round(seq(from = 1, to = 260, length.out = 99), 0)
# reds <- round(seq(from = 261, to = 400, length.out = 100), 0)
blues <- round(seq(from = 1, to = 210, length.out = 99), 0)
reds <- round(seq(from = 211, to = 400, length.out = 100), 0) # Strong reds
reds <- round(seq(from = 211, to = 350, length.out = 100), 0) # Lighter reds
pal_bl_red_brewer <- pal_bl_red_brewer[c(blues, reds)]
pal_bl_red_brewer <- c("grey90", pal_bl_red_brewer)

# Contrast raster
source("./functions/contrasting_raster.R")
hist(Ponerinae_mean_net_div_rates_WGS84[])
table(Ponerinae_mean_net_div_rates_WGS84[])
hist(Ponerinae_mean_net_div_rates_WGS84_thresholded[])
table(Ponerinae_mean_net_div_rates_WGS84_thresholded[])
Ponerinae_mean_net_div_rates_WGS84_contrasted <- contrasting_raster(# x = Ponerinae_mean_net_div_rates_WGS84, # Use the non-thresholded full version
                                                                    x = Ponerinae_mean_net_div_rates_WGS84_thresholded, # Use the thresholded version
                                                                    # zmin = 0.06, zmax = 0.14)
                                                                    zmin = 0.04, zmax = 0.16)

# Convert raster of mean current net div rates to Mollweide
Ponerinae_mean_net_div_rates_Mollweide_contrasted <- raster::projectRaster(from = Ponerinae_mean_net_div_rates_WGS84_contrasted,
                                                                           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
                                                                           method = "ngb")
plot(Ponerinae_mean_net_div_rates_Mollweide_contrasted, col = pal_bl_red_brewer)

# Convert raster to SpatialPixelsDataFrame
Ponerinae_mean_net_div_rates_Mollweide_spdf <- as(Ponerinae_mean_net_div_rates_Mollweide_contrasted, "SpatialPixelsDataFrame")
Ponerinae_mean_net_div_rates_Mollweide_spdf <- as.data.frame(Ponerinae_mean_net_div_rates_Mollweide_spdf)
colnames(Ponerinae_mean_net_div_rates_Mollweide_spdf) <- c("value", "x", "y")

# GGplot
Ponerinae_mean_net_div_rates_ggplot <- ggplot() +
  
  # Plot species richness as raster background
  geom_tile(data = Ponerinae_mean_net_div_rates_Mollweide_spdf,
            aes(x = x, y = y, fill = value), alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Net div. rates",
                       # colors = pal_bl_red_Mannion,
                       colors = pal_bl_red_brewer) +
  
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
  ggtitle(label =  paste0("Mean current net diversification rates")) +
  
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
# # pdf(file = paste0("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_mean_net_div_rates_map.pdf"),
# pdf(file = paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_net_div_rates_map_no_threshold.pdf"),
#     width = 10, height = 5)
# 
# print(Ponerinae_mean_net_div_rates_ggplot)
# 
# dev.off()

## With threshold
pdf(file = paste0("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_mean_net_div_rates_map_threshold_",min_richness_threshold,".pdf"),
    width = 10, height = 5)

print(Ponerinae_mean_net_div_rates_ggplot)

dev.off()


##### 14/ Plot and test for a longitudinal gradient ####

# Test for longitudinal gradient vs. null model of random rates 

longitude_scale <- seq(from = -180, to = 180, by = 1)

### 14.1/ Compute the mean current net. div rates along longitudinal bands ####

## Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

## Load longitudinal bands df
# longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

## 14.1.1/ Extract Net diversification rates across BAMM posterior samples

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

# Speciation rates
lambda_tip_rates_across_BAMM_posterior_samples_df <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipLambda)
names(lambda_tip_rates_across_BAMM_posterior_samples_df) <- BAMM_posterior_samples_data$tip.label

# Extinction rates
mu_tip_rates_across_BAMM_posterior_samples_df <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipMu)
names(mu_tip_rates_across_BAMM_posterior_samples_df) <- BAMM_posterior_samples_data$tip.label

# Net diversification rates = Lambda - Mu
net_div_tip_rates_across_BAMM_posterior_samples_df <- lambda_tip_rates_across_BAMM_posterior_samples_df - mu_tip_rates_across_BAMM_posterior_samples_df

# Save df of Net diversification rates across BAMM posterior samples
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")

## 14.1.2/ Compute mean current net. div rates from taxa present in each longitudinal band

## Loop per longitudinal bands
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df <- data.frame()
for (i in 1:nrow(longitudinal_bands_df))
{
  # i <- 1
  
  # Extract list of taxa present in band
  longitude_binary_presence_i <- unlist(longitude_binary_presence_df[i,, drop = T])
  taxa_in_band_i <- names(longitude_binary_presence_df)[longitude_binary_presence_i]
  taxa_in_band_i_ID <- match(x = taxa_in_band_i, table = names(net_div_tip_rates_across_BAMM_posterior_samples_df))
  taxa_in_band_i_ID <- taxa_in_band_i_ID[!is.na(taxa_in_band_i_ID)]
  
  # Extract current net div. rates for focal taxa along BAMM posterior samples
  net_div_tip_rates_in_band_i <- net_div_tip_rates_across_BAMM_posterior_samples_df[, taxa_in_band_i_ID, drop = F]
  
  # Compute mean across taxa
  net_div_tip_rates_across_BAMM_posterior_samples_i <- apply(X = net_div_tip_rates_in_band_i, MARGIN = 1, FUN = mean, na.rm = T)
  net_div_tip_rates_across_BAMM_posterior_samples_i[is.nan(net_div_tip_rates_across_BAMM_posterior_samples_i)] <- NA
  
  # Store result
  net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df <- rbind(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df, net_div_tip_rates_across_BAMM_posterior_samples_i)
}
names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df) <- paste0("BAMM_post_", 1:ncol(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df))
row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df) <- longitudinal_bands_df$longitude_dec

View(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df)

# Save df of net div. rates along BAMM posterior samples for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df.rds")

# Melt output to use in ggplot
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df$longitude <- row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df)
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df <- 
  reshape2::melt(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df) 
names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df) <- c("longitude_dec", "BAMM_post_ID", "net_div_rate")
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df %>% 
  mutate(longitude_dec = as.numeric(as.character(longitude_dec)))

# Save melted df of net div. rates along BAMM posterior samples for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df.rds")

## Compute median net div. rates across BAMM posterior samples
net_div_tip_rates_all_longitudes_median_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df %>% 
  group_by(longitude_dec) %>% 
  summarise(median_net_div_rate = median(net_div_rate)) %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  arrange(longitude_dec)

# Inform longitudinal bands df
longitudinal_bands_df <- left_join(longitudinal_bands_df, net_div_tip_rates_all_longitudes_median_df)

# Save longitudinal bands df
# saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
saveRDS(object = longitudinal_bands_df, "./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# ### 14.2/ Compute null distribution of median net div rates with BAMM-blocked strategy ####
# 
# ## ABANDONNED has not a proper null model for our question
# # Here, the goal is to test if the rates values observed in space are more extreme than expected if species distribution was random
# # Thus, must be based on permutation of species identity => all rates should be randomly permuted
# # Blocked-permutation only useful to compare tip values with tip rates, but not tip rates in space
# 
# # Make a null model with species identity randomize (i.e., randomize net div. rates across tips)
# # Use the block-strategy of STRAPP tests to permute div. rates according to regimes
# # Compute null distribution of mean net div. rates along longitudinal bands
# 
# # Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")
# 
# # Load df for taxa presence along longitudinal bands 
# longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")
# 
# # Load taxa Net diversification rates across BAMM posterior samples
# net_div_tip_rates_across_BAMM_posterior_samples_df <- readRDS(file = "./outputs/BAMM/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")
# 
# # Load longitudinal bands df
# longitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/longitudinal_bands_df.rds")
# 
# 
# length(BAMM_posterior_samples_data$tipStates[[1]])
# dim(longitude_binary_presence_df)
# 
# ## Loop per BAMM posterior samples
# set.seed(seed = 1234)
# net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df <- data.frame()
# for (i in seq_along(BAMM_posterior_samples_data$tipStates))
# {
#   # i <- 1
#   
#   # Extract tip rates
#   tip_net_div_rates_i <- BAMM_posterior_samples_data$tipLambda[[i]] - BAMM_posterior_samples_data$tipMu[[i]]
#   names(tip_net_div_rates_i) <- BAMM_posterior_samples_data$tip.label
#   
#   # Extract tip regime rates
#   tip_regimes_rates_i <- unique(tip_net_div_rates_i)
#   # Extract tip regime ID
#   tip_regimes_i <- BAMM_posterior_samples_data$tipStates[[i]]
#   unique_tip_regimes_i <- unique(tip_regimes_i)
#   names(tip_regimes_rates_i) <- unique_tip_regimes_i
#   
#   # Randomize tip regime ID
#   new_tip_regimes_i <- sample(x = unique_tip_regimes_i, size = length(unique_tip_regimes_i), replace = F)
#   new_tip_regimes_rates_i <- tip_regimes_rates_i
#   names(new_tip_regimes_rates_i) <- new_tip_regimes_i
#   
#   # Reattribute block-randomized tip rates
#   new_tip_net_div_rates_i <- new_tip_regimes_rates_i[match(x = tip_regimes_i, table = names(new_tip_regimes_rates_i))]
#   names(new_tip_net_div_rates_i) <- BAMM_posterior_samples_data$tip.label
#   
#   ## Loop per longitudinal bands
#   net_div_tip_rates_BAMM_posterior_sample_i <- c()
#   for (j in 1:nrow(longitudinal_bands_df))
#   {
#     # j <- 1
#     
#     # Extract list of taxa present in band
#     longitude_binary_presence_j <- unlist(longitude_binary_presence_df[j,, drop = T])
#     taxa_in_band_j <- names(longitude_binary_presence_df)[longitude_binary_presence_j]
#     taxa_in_band_j_ID <- match(x = taxa_in_band_j, table = names(taxa_immigration_ages_per_BSM_simulations_BAMM_perm_i))
#     taxa_in_band_j_ID <- taxa_in_band_j_ID[!is.na(taxa_in_band_j_ID)]
#     
#     # Extract null net div. rates for focal taxa
#     taxa_net_div_rates_in_band_j <- new_tip_net_div_rates_i[taxa_in_band_j_ID]
#     
#     # Compute mean across taxa
#     net_div_rates_j <- mean(taxa_net_div_rates_in_band_j)
#     net_div_rates_j[is.nan(net_div_rates_j)] <- NA
#     
#     # Store result
#     net_div_tip_rates_BAMM_posterior_sample_i <- c(net_div_tip_rates_BAMM_posterior_sample_i, net_div_rates_j)
#   }
#   names(net_div_tip_rates_BAMM_posterior_sample_i) <- longitudinal_bands_df$longitude_dec
#   
#   # Store results from BAMM posterior sample i
#   net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df <- rbind(net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df, net_div_tip_rates_BAMM_posterior_sample_i)
#   names(net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df) <- longitudinal_bands_df$longitude_dec
#   
#   ## Print progress
#   if (i %% 10 == 0)
#   {
#     cat(paste0(Sys.time(), " - Mean net div. rates along Longitudinal bands computed for BAMM posterior sample n°", i, "/", length(BAMM_posterior_samples_data$tipStates),"\n"))
#     
#     # Save null data for mean net div. rates along Longitudinal bands
#     saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df, file = "./outputs/BAMM/net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df.rds")
#   }
# }
# 
# # Save null data for mean net div. rates along Longitudinal bands
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df, file = "./outputs/BAMM/net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df.rds")
# 
# 
# # Melt output to use in ggplot
# net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df$BAMM_post_ID <- row.names(net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df)
# net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df <- 
#   reshape2::melt(net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_df)
# names(net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df) <- c("BAMM_post_ID", "longitude_dec", "net_div_rate")
# net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df <- net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df %>% 
#   mutate(longitude_dec = as.numeric(as.character(longitude_dec)))
# 
# # Save melted df of null data for net div. rates along BAMM posterior samples for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df, file = "./outputs/BAMM/net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df.rds")
# 
# 
# ## Compute median net div. rates across BAMM posterior samples
# net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_median_df <- net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_melted_df %>% 
#   group_by(longitude_dec) %>% 
#   summarise(median_net_div_rate = median(net_div_rate)) %>% 
#   mutate(longitude_dec = as.numeric(longitude_dec)) %>%
#   arrange(longitude_dec)
# 
# # Save median values of null data for net div. rates for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_median_df, file = "./outputs/BAMM/net_div_tip_rates_across_BAMM_posterior_samples_BAMM_perm_median_df.rds")


### 14.3/ Compute null distribution of median net div rates with full permutation strategy ####

# Make a null model with species identity randomize (i.e., randomize net div. rates across tips)
# Apply each permutation to all BAMM posterior samples
# Compute null distribution of mean net div. rates along longitudinal bands

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

# Load df for taxa presence along longitudinal bands 
longitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/longitude_binary_presence_df.rds")

# Load taxa Net diversification rates across BAMM posterior samples
# net_div_tip_rates_across_BAMM_posterior_samples_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")
net_div_tip_rates_across_BAMM_posterior_samples_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")

# Load longitudinal bands df
# longitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
longitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Set nb of permutations
nb_perm <- 1000
nb_taxa <- dim(net_div_tip_rates_across_BAMM_posterior_samples_df)[2]

## Loop per permutations
set.seed(seed = 1234)
net_div_tip_rates_median_all_longitudes_perm_df <- data.frame()
for (i in 1:nb_perm)
{
  # i <- 1
  
  # Randomize taxa identity
  taxa_perm_ID <- sample(x = 1:nb_taxa, size = nb_taxa, replace = F)
  
  # Randomize immigration ages in BAMM posterior samples data
  net_div_tip_rates_per_BAMM_posterior_samples_perm_i <- net_div_tip_rates_across_BAMM_posterior_samples_df[, taxa_perm_ID]
  names(net_div_tip_rates_per_BAMM_posterior_samples_perm_i) <- names(net_div_tip_rates_across_BAMM_posterior_samples_df)
  
  ## Loop per longitudinal bands
  net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i <- data.frame()
  for (j in 1:nrow(longitudinal_bands_df))
  {
    # j <- 1
    
    # Extract list of taxa present in band
    longitude_binary_presence_j <- unlist(longitude_binary_presence_df[j,, drop = T])
    taxa_in_band_j <- names(longitude_binary_presence_df)[longitude_binary_presence_j]
    taxa_in_band_j_ID <- match(x = taxa_in_band_j, table = names(net_div_tip_rates_per_BAMM_posterior_samples_perm_i))
    taxa_in_band_j_ID <- taxa_in_band_j_ID[!is.na(taxa_in_band_j_ID)]
    
    # Extract net div rates for focal taxa along BSM simulations
    taxa_net_div_rates_df_in_band_j <- net_div_tip_rates_per_BAMM_posterior_samples_perm_i[, taxa_in_band_j_ID, drop = F]
    
    # Compute mean across taxa
    net_div_rates_across_BAMM_posterior_samples_j <- apply(X = taxa_net_div_rates_df_in_band_j, MARGIN = 1, FUN = mean, na.rm = T)
    net_div_rates_across_BAMM_posterior_samples_j[is.nan(net_div_rates_across_BAMM_posterior_samples_j)] <- NA
    
    # Store result
    net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i <- rbind(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i, net_div_rates_across_BAMM_posterior_samples_j)
  }
  names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i) <- paste0("BAMM_post_", 1:ncol(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i))
  row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i) <- longitudinal_bands_df$longitude_dec
  
  # Melt output to be able to compute median ages
  net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i$longitude <- row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i)
  net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_i <- reshape2::melt(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_df_i)
  names(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_i) <- c("longitude_dec", "BAMM_post_ID", "net_div_rate")
  
  # Compute median across all BAMM posterior samples for perm i
  
  net_div_rates_median_all_longitudes_df_i <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_i %>% 
    group_by(longitude_dec) %>% 
    summarise(median_net_div_rate = median(net_div_rate)) %>% 
    mutate(longitude_dec = as.numeric(longitude_dec)) %>%
    arrange(longitude_dec)
  
  # Store results from perm i
  net_div_tip_rates_median_all_longitudes_perm_df <- rbind(net_div_tip_rates_median_all_longitudes_perm_df, net_div_rates_median_all_longitudes_df_i$median_net_div_rate)
  names(net_div_tip_rates_median_all_longitudes_perm_df) <- net_div_rates_median_all_longitudes_df_i$longitude_dec
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Median net div. rates along Longitudinal bands computed for permutation n°", i, "/", nb_perm,"\n"))
    
    # Save null data for median net div. rates along Longitudinal bands
    # saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_df.rds")
    saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_df.rds")
  }
}

# Save null data for median net div. rates along Longitudinal bands
# saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_df.rds")
saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_df.rds")

# Melt output to use in ggplot
net_div_tip_rates_median_all_longitudes_perm_df$perm_ID <- row.names(net_div_tip_rates_median_all_longitudes_perm_df)
net_div_tip_rates_median_all_longitudes_perm_melted_df <- reshape2::melt(net_div_tip_rates_median_all_longitudes_perm_df)
names(net_div_tip_rates_median_all_longitudes_perm_melted_df) <- c("perm_ID", "longitude_dec", "net_div_rate")
net_div_tip_rates_median_all_longitudes_perm_melted_df <- net_div_tip_rates_median_all_longitudes_perm_melted_df %>% 
  mutate(longitude_dec = as.numeric(as.character(longitude_dec)))

# Save melted df of null data for net div. rates along BAMM posterior samples for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_melted_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_melted_df.rds")
saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_melted_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_melted_df.rds")

## Compute median net div. rates across Permutations
net_div_tip_rates_median_all_longitudes_perm_median_df <- net_div_tip_rates_median_all_longitudes_perm_melted_df %>% 
  group_by(longitude_dec) %>% 
  summarise(median_net_div_rate = median(net_div_rate)) %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  arrange(longitude_dec)

# Save median values of null data for net div. rates for all longitudinal bands 
# saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_median_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_median_df.rds")
saveRDS(object = net_div_tip_rates_median_all_longitudes_perm_median_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_median_df.rds")


### 14.4/ Plot net div rates along longitudinal bands ####

# Load longitudinal bands df with median net div rates
# longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_rough_phylogeny_1534t.rds")
longitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/longitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Load melted df of net div. rates along BAMM posterior samples for all longitudinal bands 
# net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df.rds")
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df.rds")

# Load melted df of null data for net div. rates for all longitudinal bands 
# net_div_tip_rates_median_all_longitudes_perm_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_melted_df.rds")
net_div_tip_rates_median_all_longitudes_perm_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_melted_df.rds")

# Load median values of null data for net div. rates for all longitudinal bands 
# net_div_tip_rates_median_all_longitudes_perm_median_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_median_df.rds")
net_div_tip_rates_median_all_longitudes_perm_median_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_longitudes_perm_median_df.rds")

## 14.4.1/ Detect NA bands to cut down polygons in sections

NA_bands <- net_div_tip_rates_median_all_longitudes_perm_melted_df %>% 
  group_by(longitude_dec) %>%
  summarize(NA_band = any(is.na(net_div_rate)))

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
net_div_tip_rates_median_all_longitudes_perm_melted_df <- net_div_tip_rates_median_all_longitudes_perm_melted_df %>% 
  left_join(y = NA_bands)
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df %>% 
  left_join(y = NA_bands)

## 14.4.2/ Compute quantiles of net div rates used to draw polygons

net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  group_by(longitude_dec, poly_ID) %>% 
  # Compute quantiles
  reframe(quant_rates = quantile(net_div_rate, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
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
max_rates <- 0.15
min_rates <- 0.05
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles$quant_rates[net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles$quant_rates > max_rates] <- max_rates
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles$quant_rates[net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles$quant_rates < min_rates] <- min_rates

# Remove most extreme quantiles to keep 95% CI (strictly a 94% CI)
net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles <- net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles %>% 
  filter(!(quantile %in% c(0, 0.01, 0.02, 0.98, 0.99, 1)))

View(net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles)

## 14.4.3/ Compute CI95% quantiles of null values for net div rates used to draw polygon
net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles <- net_div_tip_rates_median_all_longitudes_perm_melted_df %>% 
  mutate(longitude_dec = as.numeric(longitude_dec)) %>%
  group_by(longitude_dec, poly_ID) %>% 
  # Compute CI95% quantiles
  reframe(CI95_rates = quantile(net_div_rate, probs = c(0.025, 0.975), na.rm = T)) %>%
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

# Adjust values of polygons to fit y-axis limits
max_rates <- 0.15
min_rates <- 0.05
net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_rates[net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$quantile == 0.975 & net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_rates > max_rates] <- max_rates
net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_rates[net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$quantile == 0.025 & net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles$CI95_rates < min_rates] <- min_rates

## 14.4.4/ Create df for legend

net_div_rates_per_longitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                  y = c(0.05, 0.05),
                                                  data_type = c("Observed", "Null"))

## 14.4.5/ Generate plot using lines

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_rates_per_longitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_rates_per_longitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

net_div_rates_per_longitudinal_bands_ggplot_fuzzy_lines <- ggplot(data = net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_rates, x = longitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = net_div_tip_rates_median_all_longitudes_perm_median_df,
            mapping = aes(y = median_net_div_rate, x = longitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 1000 line replicates of observed data
  geom_line(data = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df,
            mapping = aes(y = net_div_rate, x = longitude_dec, group = BAMM_post_ID),
            col = "red",
            alpha = 0.01,
            linewidth = 3.0) +
  
  # Plot median lines of observed data
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = median_net_div_rate, x = longitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = net_div_rates_per_longitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            linewidth = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Set Y-axis limits
  # scale_y_continuous(limits = c(0, 0.2)) +
  scale_y_continuous(limits = c(min_rates, max_rates)) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Longitude", breaks = c(-120, -60, 0, 60, 120), labels = c("120°W", "60°W", "0°", "60°E", "120°E")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Net div. rates across longitudes")) +
  
  # Set axes labels
  xlab("Longitude") +
  ylab("Net div. rates  [Events / Lineage / My]") +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        # legend.position.inside = c(0.8, 0.80),
        legend.position.inside = c(0.8, 0.15),
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
print(net_div_rates_per_longitudinal_bands_ggplot_fuzzy_lines)

dev.off()

## 14.4.6/ Generate plot using quantiles

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_rates_per_longitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_rates_per_longitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

net_div_rates_per_longitudinal_bands_ggplot_fuzzy_quantiles <- ggplot(data = net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = net_div_tip_rates_median_all_longitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_rates, x = longitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = net_div_tip_rates_median_all_longitudes_perm_median_df,
            mapping = aes(y = median_net_div_rate, x = longitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 50 quantile polygons of observed data
  geom_polygon(data = net_div_tip_rates_across_BAMM_posterior_samples_all_longitudes_melted_df_quantiles,
               mapping = aes(y = quant_rates, x = longitude_dec, group = interaction(quantile_ID, poly_ID)),
               fill = "red",
               alpha = 0.03,
               # alpha = 0.02,
               linewidth = 2.0) +
  
  # Plot median lines of observed data
  geom_line(data = longitudinal_bands_df,
            mapping = aes(y = median_net_div_rate, x = longitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Add fake data for legend
  geom_line(data = net_div_rates_per_longitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            size = 0) +
  
  # Set Y-axis limits
  # scale_y_continuous(limits = c(0, 0.2)) +
  scale_y_continuous(limits = c(min_rates, max_rates)) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Longitude", breaks = c(-120, -60, 0, 60, 120), labels = c("120°W", "60°W", "0°", "60°E", "120°E")) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Net div. rates across longitudes")) +
  
  # Set axes labels
  xlab("Longitude") +
  ylab("Net div. rates  [Events / Lineage / My]") +

  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
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
print(net_div_rates_per_longitudinal_bands_ggplot_fuzzy_quantiles)

dev.off()


##### 15/ Plot and test for a latitudinal gradient ####

# Test for latitudinal gradient vs. null model of random rates

latitude_scale <- seq(from = -60, to = 60, by = 1)

### 15.1/ Compute the mean current net. div rates along latitudinal bands ####

## Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

## Load latitudinal bands df
# latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

## 15.1.1/ Extract Net diversification rates across BAMM posterior samples

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

# Speciation rates
lambda_tip_rates_across_BAMM_posterior_samples_df <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipLambda)
names(lambda_tip_rates_across_BAMM_posterior_samples_df) <- BAMM_posterior_samples_data$tip.label

# Extinction rates
mu_tip_rates_across_BAMM_posterior_samples_df <- do.call(rbind.data.frame, BAMM_posterior_samples_data$tipMu)
names(mu_tip_rates_across_BAMM_posterior_samples_df) <- BAMM_posterior_samples_data$tip.label

# Net diversification rates = Lambda - Mu
net_div_tip_rates_across_BAMM_posterior_samples_df <- lambda_tip_rates_across_BAMM_posterior_samples_df - mu_tip_rates_across_BAMM_posterior_samples_df

# Save df of Net diversification rates across BAMM posterior samples
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")

## 15.1.2/ Compute mean current net. div rates from taxa present in each longitudinal band

## Loop per latitudinal bands
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df <- data.frame()
for (i in 1:nrow(latitudinal_bands_df))
{
  # i <- 1
  
  # Extract list of taxa present in band
  latitude_binary_presence_i <- unlist(latitude_binary_presence_df[i,, drop = T])
  taxa_in_band_i <- names(latitude_binary_presence_df)[latitude_binary_presence_i]
  taxa_in_band_i_ID <- match(x = taxa_in_band_i, table = names(net_div_tip_rates_across_BAMM_posterior_samples_df))
  taxa_in_band_i_ID <- taxa_in_band_i_ID[!is.na(taxa_in_band_i_ID)]
  
  # Extract current net div. rates for focal taxa along BAMM posterior samples
  net_div_tip_rates_in_band_i <- net_div_tip_rates_across_BAMM_posterior_samples_df[, taxa_in_band_i_ID, drop = F]
  
  # Compute mean across taxa
  net_div_tip_rates_across_BAMM_posterior_samples_i <- apply(X = net_div_tip_rates_in_band_i, MARGIN = 1, FUN = mean, na.rm = T)
  net_div_tip_rates_across_BAMM_posterior_samples_i[is.nan(net_div_tip_rates_across_BAMM_posterior_samples_i)] <- NA
  
  # Store result
  net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df <- rbind(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df, net_div_tip_rates_across_BAMM_posterior_samples_i)
}
names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df) <- paste0("BAMM_post_", 1:ncol(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df))
row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df) <- latitudinal_bands_df$latitude_dec

View(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df)

# Save df of net div. rates along BAMM posterior samples for all latitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df.rds")

# Melt output to use in ggplot
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df$latitude <- row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df)
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df <- 
  reshape2::melt(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df) 
names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df) <- c("latitude_dec", "BAMM_post_ID", "net_div_rate")
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df %>% 
  mutate(latitude_dec = as.numeric(as.character(latitude_dec)))

# Save melted df of net div. rates along BAMM posterior samples for all latitudinal bands 
# saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df.rds")
saveRDS(object = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df.rds")

## Compute median net div. rates across BAMM posterior samples
net_div_tip_rates_all_latitudes_median_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df %>% 
  group_by(latitude_dec) %>% 
  summarise(median_net_div_rate = median(net_div_rate)) %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  arrange(latitude_dec)

# Inform latitudinal bands df
latitudinal_bands_df <- left_join(latitudinal_bands_df, net_div_tip_rates_all_latitudes_median_df)

# Save latitudinal bands df
# saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
saveRDS(object = latitudinal_bands_df, "./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

### 15.2/ Compute null distribution of median net div rates with full permutation strategy ####

# Make a null model with species identity randomize (i.e., randomize net div. rates across tips)
# Apply each permutation to all BAMM posterior samples
# Compute null distribution of mean net div. rates along latitudinal bands

# Load the BAMM posterior samples object
# BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/BAMM_posterior_samples_data.rds")
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

# Load df for taxa presence along latitudinal bands 
latitude_binary_presence_df <- readRDS(file = "./outputs/Species_richness_maps/latitude_binary_presence_df.rds")

# Load taxa Net diversification rates across BAMM posterior samples
# net_div_tip_rates_across_BAMM_posterior_samples_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")
net_div_tip_rates_across_BAMM_posterior_samples_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_df.rds")

# Load latitudinal bands df
# latitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
latitudinal_bands_df <- readRDS("./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Set nb of permutations
nb_perm <- 1000
nb_taxa <- dim(net_div_tip_rates_across_BAMM_posterior_samples_df)[2]

## Loop per permutations
set.seed(seed = 1234)
net_div_tip_rates_median_all_latitudes_perm_df <- data.frame()
for (i in 1:nb_perm)
{
  # i <- 1
  
  # Randomize taxa identity
  taxa_perm_ID <- sample(x = 1:nb_taxa, size = nb_taxa, replace = F)
  
  # Randomize immigration ages in BAMM posterior samples data
  net_div_tip_rates_per_BAMM_posterior_samples_perm_i <- net_div_tip_rates_across_BAMM_posterior_samples_df[, taxa_perm_ID]
  names(net_div_tip_rates_per_BAMM_posterior_samples_perm_i) <- names(net_div_tip_rates_across_BAMM_posterior_samples_df)
  
  ## Loop per latitudinal bands
  net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i <- data.frame()
  for (j in 1:nrow(latitudinal_bands_df))
  {
    # j <- 1
    
    # Extract list of taxa present in band
    latitude_binary_presence_j <- unlist(latitude_binary_presence_df[j,, drop = T])
    taxa_in_band_j <- names(latitude_binary_presence_df)[latitude_binary_presence_j]
    taxa_in_band_j_ID <- match(x = taxa_in_band_j, table = names(net_div_tip_rates_per_BAMM_posterior_samples_perm_i))
    taxa_in_band_j_ID <- taxa_in_band_j_ID[!is.na(taxa_in_band_j_ID)]
    
    # Extract net div rates for focal taxa along BSM simulations
    taxa_net_div_rates_df_in_band_j <- net_div_tip_rates_per_BAMM_posterior_samples_perm_i[, taxa_in_band_j_ID, drop = F]
    
    # Compute mean across taxa
    net_div_rates_across_BAMM_posterior_samples_j <- apply(X = taxa_net_div_rates_df_in_band_j, MARGIN = 1, FUN = mean, na.rm = T)
    net_div_rates_across_BAMM_posterior_samples_j[is.nan(net_div_rates_across_BAMM_posterior_samples_j)] <- NA
    
    # Store result
    net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i <- rbind(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i, net_div_rates_across_BAMM_posterior_samples_j)
  }
  names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i) <- paste0("BAMM_post_", 1:ncol(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i))
  row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i) <- latitudinal_bands_df$latitude_dec
  
  # Melt output to be able to compute median ages
  net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i$latitude <- row.names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i)
  net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_i <- reshape2::melt(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_df_i)
  names(net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_i) <- c("latitude_dec", "BAMM_post_ID", "net_div_rate")
  
  # Compute median across all BAMM posterior samples for perm i
  
  net_div_rates_median_all_latitudes_df_i <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_i %>% 
    group_by(latitude_dec) %>% 
    summarise(median_net_div_rate = median(net_div_rate)) %>% 
    mutate(latitude_dec = as.numeric(latitude_dec)) %>%
    arrange(latitude_dec)
  
  # Store results from perm i
  net_div_tip_rates_median_all_latitudes_perm_df <- rbind(net_div_tip_rates_median_all_latitudes_perm_df, net_div_rates_median_all_latitudes_df_i$median_net_div_rate)
  names(net_div_tip_rates_median_all_latitudes_perm_df) <- net_div_rates_median_all_latitudes_df_i$latitude_dec
  
  ## Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Median net div. rates along latitudinal bands computed for permutation n°", i, "/", nb_perm,"\n"))
    
    # Save null data for median net div. rates along latitudinal bands
    # saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_df.rds")
    saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_df.rds")
  }
}

# Save null data for median net div. rates along latitudinal bands
# saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_df.rds")
saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_df.rds")

# Melt output to use in ggplot
net_div_tip_rates_median_all_latitudes_perm_df$perm_ID <- row.names(net_div_tip_rates_median_all_latitudes_perm_df)
net_div_tip_rates_median_all_latitudes_perm_melted_df <- reshape2::melt(net_div_tip_rates_median_all_latitudes_perm_df)
names(net_div_tip_rates_median_all_latitudes_perm_melted_df) <- c("perm_ID", "latitude_dec", "net_div_rate")
net_div_tip_rates_median_all_latitudes_perm_melted_df <- net_div_tip_rates_median_all_latitudes_perm_melted_df %>% 
  mutate(latitude_dec = as.numeric(as.character(latitude_dec)))

# Save melted df of null data for net div. rates along BAMM posterior samples for all latitudinal bands 
# saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_melted_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_melted_df.rds")
saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_melted_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_melted_df.rds")

## Compute median net div. rates across Permutations
net_div_tip_rates_median_all_latitudes_perm_median_df <- net_div_tip_rates_median_all_latitudes_perm_melted_df %>% 
  group_by(latitude_dec) %>% 
  summarise(median_net_div_rate = median(net_div_rate)) %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  arrange(latitude_dec)

# Save median values of null data for net div. rates for all latitudinal bands 
# saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_median_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_median_df.rds")
saveRDS(object = net_div_tip_rates_median_all_latitudes_perm_median_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_median_df.rds")


### 15.3/ Plot net div rates along latitudinal bands ####

# Load latitudinal bands df with median net div rates
# latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_rough_phylogeny_1534t.rds")
latitudinal_bands_df <- readRDS(file = "./outputs/Species_richness_maps/latitudinal_bands_df_for_MCC_phylogeny_1534t.rds")

# Load melted df of net div. rates along BAMM posterior samples for all latitudinal bands 
# net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df.rds")
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df.rds")

# Load melted df of null data for net div. rates for all latitudinal bands 
# net_div_tip_rates_median_all_latitudes_perm_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_melted_df.rds")
net_div_tip_rates_median_all_latitudes_perm_melted_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_melted_df.rds")

# Load median values of null data for net div. rates for all latitudinal bands 
# net_div_tip_rates_median_all_latitudes_perm_median_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_median_df.rds")
net_div_tip_rates_median_all_latitudes_perm_median_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_tip_rates_median_all_latitudes_perm_median_df.rds")

## 15.3.1/ Detect NA bands to cut down polygons in sections

NA_bands <- net_div_tip_rates_median_all_latitudes_perm_melted_df %>% 
  group_by(latitude_dec) %>%
  summarize(NA_band = any(is.na(net_div_rate)))

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
net_div_tip_rates_median_all_latitudes_perm_melted_df <- net_div_tip_rates_median_all_latitudes_perm_melted_df %>% 
  left_join(y = NA_bands)
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df %>% 
  left_join(y = NA_bands)

## 15.3.2/ Compute quantiles of net div rates used to draw polygons

net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  group_by(latitude_dec, poly_ID) %>% 
  # Compute quantiles
  reframe(quant_rates = quantile(net_div_rate, probs = seq(from = 0, to = 1, by = 0.01), na.rm = T)) %>%
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

# Adjust min/max values to the plot limits to avoid artefct in polygon drawings
min_rates <- 0.05
max_rates <- 0.15
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles$quant_rates[net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles$quant_rates > max_rates] <- max_rates
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles$quant_rates[net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles$quant_rates < min_rates] <- min_rates

# Remove most extreme quantiles to keep 95% CI (strictly a 94% CI)
net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles <- net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles %>% 
  filter(!(quantile %in% c(0, 0.01, 0.02, 0.98, 0.99, 1)))

## 15.3.3/ Compute CI95% quantiles of null values for net div rates used to draw polygon
net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles <- net_div_tip_rates_median_all_latitudes_perm_melted_df %>% 
  mutate(latitude_dec = as.numeric(latitude_dec)) %>%
  group_by(latitude_dec, poly_ID) %>% 
  # Compute CI95% quantiles
  reframe(CI95_rates = quantile(net_div_rate, probs = c(0.025, 0.975), na.rm = T)) %>%
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

# Adjust min/max values to the plot limits to avoid artefct in polygon drawings
net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_rates[net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_rates > max_rates] <- max_rates
net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_rates[net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles$CI95_rates < min_rates] <- min_rates


## 15.3.4/ Create df for legend

net_div_rates_per_latitudinal_bands_legend_df <- data.frame(x = c(0, 0),
                                                             y = c(0.05, 0.05),
                                                             data_type = c("Observed", "Null"))

## 15.3.5/ Generate plot using lines

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_rates_per_latitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_rates_per_latitudinal_bands_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

net_div_rates_per_latitudinal_bands_ggplot_fuzzy_lines <- ggplot(data = net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_rates, x = latitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = net_div_tip_rates_median_all_latitudes_perm_median_df,
            mapping = aes(y = median_net_div_rate, x = latitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 1000 line replicates of observed data
  geom_line(data = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df,
            mapping = aes(y = net_div_rate, x = latitude_dec, group = BAMM_post_ID),
            col = "red",
            alpha = 0.01,
            linewidth = 3.0) +
  
  # Plot median lines of observed data
  geom_line(data = latitudinal_bands_df,
            mapping = aes(y = median_net_div_rate, x = latitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Flip coordinates to have vertical latitude
  coord_flip() +
  
  # Add fake data for legend
  geom_line(data = net_div_rates_per_latitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            linewidth = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Net div. rates\nacross latitudes")) +
  
  # Set axes labels
  xlab("Latitude") +
  ylab("Net div. rates  [Events / Lineage / My]") +
  
  # Set y-axis limits
  ylim(c(min_rates, max_rates)) +
  # ylim(c(0, y_max)) +
  
  # Adjust label on Latitude axis
  scale_x_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.50),
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
print(net_div_rates_per_latitudinal_bands_ggplot_fuzzy_lines)

dev.off()

## 15.3.6/ Generate plot using quantiles

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/net_div_rates_per_latitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/net_div_rates_per_latitudinal_bands_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

net_div_rates_per_latitudinal_bands_ggplot_fuzzy_quantiles <- ggplot(data = net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles) +
  
  # Plot CI95% polygon for null data
  geom_polygon(data = net_div_tip_rates_median_all_latitudes_perm_melted_df_CI95_quantiles,
               mapping = aes(y = CI95_rates, x = latitude_dec, group = poly_ID),
               fill = "grey80",
               alpha = 0.5,
               linewidth = 1.0) +
  
  # Plot median line of null data
  geom_line(data = net_div_tip_rates_median_all_latitudes_perm_median_df,
            mapping = aes(y = median_net_div_rate, x = latitude_dec),
            col = "black",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Plot 50 quantile polygons of observed data
  geom_polygon(data = net_div_tip_rates_across_BAMM_posterior_samples_all_latitudes_melted_df_quantiles,
               mapping = aes(y = quant_rates, x = latitude_dec, group = interaction(quantile_ID, poly_ID)),
               fill = "red",
               alpha = 0.03,
               linewidth = 2.0) +
  
  # Plot median lines of observed data
  geom_line(data = latitudinal_bands_df,
            mapping = aes(y = median_net_div_rate, x = latitude_dec),
            col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Flip coordinates to have vertical latitude
  coord_flip() +
  
  # Add fake data for legend
  geom_line(data = net_div_rates_per_latitudinal_bands_legend_df,
            mapping = aes(x = x, y = y, col = data_type),
            size = 0) +
  
  # Adjust color scheme and legend
  scale_color_manual("Data", breaks = c("Observed", "Null"), labels = c("Observed", "Null"),
                     values = c("red", "grey90")) +
  
  # Adjust legend aesthetics
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  
  # Set plot title +
  ggtitle(label = paste0("Net div. rates\nacross latitudes")) +
  
  # Set axes labels
  xlab("Latitude") +
  ylab("Net div. rates  [Events / Lineage / My]") +
  
  # Set y-axis limits
  ylim(c(min_rates, max_rates)) +
  # ylim(c(0, y_max)) +
  
  # Adjust label on Latitude axis
  # scale_x_continuous("Latitude", breaks = c(-60, -30, 0, 30, 60), labels = c("60°S", "30°S", "0°", "30°N", "60°N")) +
  scale_x_continuous("Latitude", limits = c(-60, 80), breaks = c(-40, 0, 40, 80), labels = c("40°S", "0°", "40°N", "80°N")) +
  
  # Adjust aesthetics
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        panel.grid.major = element_line(color = NA, linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(size  = 20, margin = margin(b = 8)),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.50),
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
print(net_div_rates_per_latitudinal_bands_ggplot_fuzzy_quantiles)

dev.off()

### 15.4/ Check effect of species range on current diversification rates ####

# Load taxa current median tip rates
# current_median_tip_rates_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/current_median_tip_rates_df.rds")
current_median_tip_rates_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/current_median_tip_rates_df.rds")

# Load alpha-hull ranges stack
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

## Loop per layers to convert raster stack into sf
Ponerinae_alpha_hull_WGS84_sf <- data.frame()
for (i in 1:nlayers(Ponerinae_alpha_hull_stack_WGS84))
{
  # i <- 1
  
  # Check number of presence pixels
  pixel_area <- sum(Ponerinae_alpha_hull_stack_WGS84[[i]]@data@values == 1, na.rm = T)
  
  if (pixel_area > 1)
  {
    # Convert to sp, then sf
    Ponerinae_alpha_hull_WGS84_sp_i <- raster::rasterToPolygons(x = Ponerinae_alpha_hull_stack_WGS84[[i]], fun = function (x) { x == 1}, dissolve = TRUE)
    if (!is.null(Ponerinae_alpha_hull_WGS84_sp_i))
    {
      Ponerinae_alpha_hull_WGS84_sf_i <- st_as_sf(Ponerinae_alpha_hull_WGS84_sp_i)
      Ponerinae_alpha_hull_WGS84_sf_i$taxa <-  names(Ponerinae_alpha_hull_stack_WGS84[[i]])
      
      # Store geometry in final sf object
      Ponerinae_alpha_hull_WGS84_sf <- rbind(Ponerinae_alpha_hull_WGS84_sf, Ponerinae_alpha_hull_WGS84_sf_i[, "taxa"])
    }
  }

  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Range raster converted to sf for taxa n°", i, "/", nlayers(Ponerinae_alpha_hull_stack_WGS84),"\n"))
  }
}

plot(Ponerinae_alpha_hull_WGS84_sf)

## Compute taxa range area
Ponerinae_alpha_hull_WGS84_sf$area <- units::set_units(st_area(Ponerinae_alpha_hull_WGS84_sf), km^2)

## Save alpha-hull ranges as sf
saveRDS(object = Ponerinae_alpha_hull_WGS84_sf, file = "./outputs/Species_richness_maps/Ponerinae_alpha_hull_WGS84_sf.rds")

## Inform range area data
min(Ponerinae_alpha_hull_WGS84_sf$area)
Taxa_ranges_rates_df <- left_join(current_median_tip_rates_df, st_drop_geometry(Ponerinae_alpha_hull_WGS84_sf))

# Set missing taxa to 1000km²
Taxa_ranges_rates_df$area[is.na(Taxa_rates_ranges_df$area)] <- 1000

## Save Taxa ranges/rates df
# saveRDS(object = Taxa_ranges_rates_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Taxa_ranges_rates_df.rds")
saveRDS(object = Taxa_ranges_rates_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Taxa_ranges_rates_df.rds")


## Run LM to get predict and test

# Check if log-transformation and sqrt transformation are needed
hist(Taxa_ranges_rates_df$net_div)
hist(log1p(Taxa_ranges_rates_df$net_div))
hist(Taxa_ranges_rates_df$area)
hist(sqrt(Taxa_ranges_rates_df$area))

plot(Taxa_ranges_rates_df$net_div ~ Taxa_ranges_rates_df$area)
plot(log(Taxa_ranges_rates_df$net_div) ~ Taxa_ranges_rates_df$area)
plot(Taxa_ranges_rates_df$net_div ~ sqrt(Taxa_ranges_rates_df$area))
plot(log(Taxa_ranges_rates_df$net_div) ~ sqrt(Taxa_ranges_rates_df$area))

# Run LM
Taxa_ranges_rates_lm <- lm(data = Taxa_ranges_rates_df, formula = net_div ~ sqrt(area))

# Extract evaluation
summary(Taxa_ranges_rates_lm)
Taxa_ranges_rates_lm_anova <- anova(Taxa_ranges_rates_lm)
R_sq <- Taxa_ranges_rates_lm_anova$`Sum Sq`[1] / sum(Taxa_ranges_rates_lm_anova$`Sum Sq`)

## True test should be a PGLS to correct for phylogenetic non-independence.
# But here we are not looking for an evolutionary relationship, just for an effect on range size on rates randomization

# Extract coefficients
Taxa_ranges_rates_lm_coef <- coef(Taxa_ranges_rates_lm)
# Extract standardized coef
Taxa_ranges_rates_lm_std <- lm(data = Taxa_ranges_rates_df, formula = scale(net_div) ~ scale(sqrt(area)))
Taxa_ranges_rates_lm_coef_std <- coef(Taxa_ranges_rates_lm_std)

## Plot current net div. rates ~ range area 

# pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Taxa_rates_ranges_ggplot.pdf", height = 6, width = 8)
pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Taxa_rates_ranges_ggplot.pdf", height = 6, width = 8)

Taxa_rates_ranges_ggplot <- ggplot(data = Taxa_ranges_rates_df) +
  
  geom_point(mapping = aes(y = net_div, x = as.numeric(area))) +
              
  geom_abline(intercept = Taxa_ranges_rates_lm_coef[1],
              slope = Taxa_ranges_rates_lm_coef[2],
              col = "red",
              linewidth = 1.2) +
  
  scale_x_continuous(transform = "sqrt") +
  
  annotate(geom = "text",
           x = 21000000, y = 0.225,
           hjust = 0,
           fontface = "bold", size = 5,
           label = paste0("R² = ", round(R_sq, 2), "\n",
                          # "Beta = ", round(Taxa_ranges_rates_lm_coef_std[2], 2), "\n",
                          "\u03B2 = ", round(Taxa_ranges_rates_lm_coef_std[2], 2), "\n",
                          # bquote(beta~round(Taxa_ranges_rates_lm_coef_std[2])), "\n",
                          "p < 0.001")) +
  
  # Add title
  ggtitle(label = "Net div. rates per taxa range") +
  
  # Add axis labels
  xlab("Area  [Km²]") +
  ylab("Net diversification rates\n[Events / Lineage / My]") +
  
  theme(panel.grid.major = element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 20, hjust = 0.5, color = "black", margin = margin(b = 15, t = 5)),
        axis.title = element_text(size = 20, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)))

print(Taxa_rates_ranges_ggplot)

dev.off()



##### 16/ Investigate relationship between PP of rate shift vs. PP of dispersal event ####

### 16.1/ Compute PP of dispersal events on branch ####

# Include both anagenetic events (d: range extension) and cladogenetic events (j: jump dispersal) where the parental area differ from the descendant area
# Get counts on each BSM
# Get frequencies across BSM to obtain Posterior probabilities

# Load phylogeny
# Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_treedata.rds")
phy <- Ponerinae_phylogeny_1534t_treedata@phylo

# Load BSM outputs
# DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_rough_phylogeny_1534t/DEC_J_BSM_output.rds")
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/Ponerinae_MCC_phylogeny_1534t/DEC_J_BSM_output.rds")

# Extract records of events
DEC_J_clado_events_tables <- DEC_J_BSM_output$RES_clado_events_tables
DEC_J_ana_events_tables <- DEC_J_BSM_output$RES_ana_events_tables

View(DEC_J_ana_events_tables[[1]])

DEC_J_clado_events_tables[[1]]$node # Node ID as in ape
DEC_J_clado_events_tables[[1]]$parent_br # ID of the edge leading to the node = gives you the row of tree$edge and tree$edge.length)

# States description. States are sampled using the MLE scaled marginal likelihoods
DEC_J_clado_events_tables[[1]]$sampled_states_AT_nodes # Node state in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$sampled_states_AT_brbots # State at the start of the branch leading to the node in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$samp_LEFT_dcorner # State at the start of the left descendant edge of the node in this particular history / stochastic map
DEC_J_clado_events_tables[[1]]$samp_RIGHT_dcorner # State at the start of the right descendant edge of the node in this particular history / stochastic map

## 16.1.1/ Counts jump-dispersal events per branches ####

## Loop per BSM
jumping_edges_ID_list <- list() # Initiate list to record ID of branches with jump dispersal (j) events
for (i in seq_along(DEC_J_clado_events_tables))
{
  # i <- 1
  
  # Extract anagenetic event table
  clado_events_table_i <- DEC_J_clado_events_tables[[i]]
  
  # Extract the founder events ($clado_event_type == "founder (j)")
  clado_events_table_i <- clado_events_table_i %>% 
    filter(clado_event_type == "founder (j)")
  
  # Identify which descendant (LEFT or RIGHT) dispersed.
  clado_events_table_i$LEFT_j <- clado_events_table_i$sampled_states_AT_nodes != clado_events_table_i$samp_LEFT_dcorner
  clado_events_table_i$RIGHT_j <- clado_events_table_i$sampled_states_AT_nodes != clado_events_table_i$samp_RIGHT_dcorner
  # table(clado_events_table_i$LEFT_j, clado_events_table_i$RIGHT_j)
  
  # Get jumping nodes ID from daughter_nds
  table(clado_events_table_i$LEFT_j)
  left_jumping_nodes_ID_i <- unlist(lapply(X = clado_events_table_i$daughter_nds[clado_events_table_i$LEFT_j], FUN = function (x) { x[1] } ))
  right_jumping_nodes_ID_i <- unlist(lapply(X = clado_events_table_i$daughter_nds[clado_events_table_i$RIGHT_j], FUN = function (x) { x[2] } ))
  jumping_nodes_ID_i <- unique(c(left_jumping_nodes_ID_i, right_jumping_nodes_ID_i))
  jumping_nodes_ID_i <- jumping_nodes_ID_i[order(jumping_nodes_ID_i)]
  
  # Get edge ID from the tipward_node_ID
  jumping_edges_ID_i <- which(phy$edge[, 2] %in% jumping_nodes_ID_i)
  
  # Store output
  jumping_edges_ID_list[[i]] <- jumping_edges_ID_i
  
  # Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Branches with jump-dispersal detected for BSM n°", i, "/", length(DEC_J_clado_events_tables),"\n"))
  }
}

## 16.1.2/ Counts range extension events per branches ####

## Loop per BSM
range_extension_edges_ID_list <- list() # Initiate list to record ID of branches with jump dispersal (j) events
for (i in seq_along(DEC_J_ana_events_tables))
{
  # i <- 1
  
  # Extract anagenetic event table
  ana_events_table_i <- DEC_J_ana_events_tables[[i]]
  
  # Extract the range extension events ($event_type == "d")
  ana_events_table_i <- ana_events_table_i %>% 
    filter(event_type == "d")
  
  # Get range extension nodes ID
  range_extension_nodes_ID_i <- ana_events_table_i$node
  range_extension_nodes_ID_i <- range_extension_nodes_ID_i[order(range_extension_nodes_ID_i)]
  
  # Get edge ID from the tipward_node_ID
  range_extension_edges_ID_i <- which(phy$edge[, 2] %in% range_extension_nodes_ID_i)
  
  # Store output
  range_extension_edges_ID_list[[i]] <- range_extension_edges_ID_i
  
  # Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Branches with range extension detected for BSM n°", i, "/", length(DEC_J_ana_events_tables),"\n"))
  }
}

## 16.1.3/ Aggregate counts of jump dispersal and range extension = dispersal events

# Initiate template vector to record binary dispersal events per edges
all_edges_ID <- as.character(1:nrow(phy$edge))
dispersal_edges_binary_template <- rep(NA, length(all_edges_ID))
names(dispersal_edges_binary_template) <- all_edges_ID

# Aggregate dispersal events per BSM in binary vectors
dispersal_edges_binary_list <- list()
for (i in seq_along(jumping_edges_ID_list))
{
  # i <- 1
  
  # Extract edge per events
  jumping_edges_ID_i <- jumping_edges_ID_list[[i]]
  range_extension_edges_ID_i <- range_extension_edges_ID_list[[i]]
  
  # Aggregate dispersal events
  dispersal_edges_ID_i <- unique(c(jumping_edges_ID_i, range_extension_edges_ID_i))
  
  # Fill binary vector
  dispersal_edges_binary_i <- dispersal_edges_binary_template
  dispersal_edges_binary_i <- names(dispersal_edges_binary_i) %in% dispersal_edges_ID_i
  # table(dispersal_edges_binary_i)
  
  # Store binary vector
  dispersal_edges_binary_list[[i]] <- dispersal_edges_binary_i
}

# Convert to df
dispersal_edges_binary_df <- do.call(rbind.data.frame, dispersal_edges_binary_list)
names(dispersal_edges_binary_df) <- all_edges_ID 
  # rows = BSM
  # cols = edges

# Save df of binary dispersal events of branches per BSM
# saveRDS(object = dispersal_edges_binary_df, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/dispersal_edges_binary_df.rds")
saveRDS(object = dispersal_edges_binary_df, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/dispersal_edges_binary_df.rds")

## 16.1.4/ Compute posterior probability of dispersal as relative frequencies ####

dispersal_PP_per_edges <- apply(X = dispersal_edges_binary_df, MARGIN = 2, FUN = sum)
dispersal_PP_per_edges <- dispersal_PP_per_edges / nrow(dispersal_edges_binary_df)

hist(dispersal_PP_per_edges)
# Bimodal distribution !

# Save posterior probability of dispersal per branches
# saveRDS(object = dispersal_PP_per_edges, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/dispersal_PP_per_edges.rds")
saveRDS(object = dispersal_PP_per_edges, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/dispersal_PP_per_edges.rds")


### 16.2/ Test for relationship between PP of diversification rate shift and PP of dispersal ####

# No valid test to account for the non-independence of edge data
# Degrees of freedom should be adjusted in the test, but not sure how...
# Only compute a stat and display a plot. Do not trust p-values 

## 16.2.1/ Convert PP to binary variables ####

## Load posterior probability of dispersal per branches
# dispersal_PP_per_edges <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/dispersal_PP_per_edges.rds")
dispersal_PP_per_edges <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/dispersal_PP_per_edges.rds")

# Distribution of dispersal PP is bimodal (or at least strongly right-skewed) so convert into a binary variable for the test
hist(dispersal_PP_per_edges)

# Apply threshold
dispersal_threshold <- 0.5
# dispersal_threshold <- 0.95
binary_dispersal_PP_per_edges <- (dispersal_PP_per_edges >= dispersal_threshold)
table(binary_dispersal_PP_per_edges)

## Load df of edge prior probs, posterior probs, and odd-ratio of rate shifts
# edge_rate_shift_probs_df <- readRDS(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/edge_rate_shift_probs_df.rds")
edge_rate_shift_probs_df <- readRDS(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/edge_rate_shift_probs_df.rds")

# Distribution of posterior probabilities and odd-ratios are also highly right-skewed
# Use a threshold to identify credible shift locations
# Does not mean a shift is likely to have occurred on all these edges as odd-ratio are independent, but in reality rate shift events alongside successive branches are not
hist(edge_rate_shift_probs_df$marg_posterior_probs)
hist(edge_rate_shift_probs_df$marginal_odd_ratios)

# Apply odd-ratio > 5 as a threshold to detect credible locations for rate shifts

binary_rate_shifts_PP_per_edges <- (edge_rate_shift_probs_df$marginal_odd_ratios > 5)
table(binary_rate_shifts_PP_per_edges)

## 16.2.2/ Apply Fisher's exact test ####

binary_rate_shifts_vs_dispersal_fisher_test <- fisher.test(x = binary_rate_shifts_PP_per_edges,
                                                           y = binary_dispersal_PP_per_edges,
                                                           alternative = "two.sided")

binary_rate_shifts_vs_dispersal_fisher_test
binary_rate_shifts_vs_dispersal_fisher_test$p.value  # p = 0.65
binary_rate_shifts_vs_dispersal_fisher_test$estimate # Odds-ratio = 0.689
binary_rate_shifts_vs_dispersal_fisher_test$conf.int # CI 95% = 0.180 - 1.891 (encompass 1 so not significant)

## 16.2.3/ Apply Independence Khi² test ####

binary_rate_shifts_vs_dispersal_chisq_test <- chisq.test(x = binary_rate_shifts_PP_per_edges, y = binary_dispersal_PP_per_edges)

binary_rate_shifts_vs_dispersal_chisq_test
binary_rate_shifts_vs_dispersal_chisq_test$statistic  # Khi² = 0.24
binary_rate_shifts_vs_dispersal_chisq_test$p.value    # p = 0.62
khi_Q95 <- qchisq(p = 0.95, df = binary_rate_shifts_vs_dispersal_chisq_test$parameter) ; khi_Q95  # Q95% = 3.841


### 16.3/ Association plot for Independence of rate shifts vs. dispersal events ####

# Compute contingency table
contingency_table <- table(binary_rate_shifts_PP_per_edges, binary_dispersal_PP_per_edges)
dimnames(contingency_table) <- list(`Rate shift events` = c("No shift", "Rate shift"),
                                    `Dispersal events` = c("No dispersal", "Dispersal"))

# Plot association plot
# pdf("./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Association_plot_rate_shifts_vs_dispersal.pdf", height = 8, width = 10)
pdf("./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Association_plot_rate_shifts_vs_dispersal.pdf", height = 8, width = 10)

par(mar = c(5.1, 4.1, 7.1, 2.1))
graphics::mosaicplot(x = contingency_table, shade = TRUE, las = 0, cex.axis = 1.5,
                     xlab = "", ylab = "",
                     main = "")
# Add title
mtext(text = "Dispersal events", side = 2, cex = 2, line = 1)
# Add axis titles
mtext(text = "Rate shifts", side = 1, cex = 2, line = 2)
mtext(text = "Independence test\nRate shifts vs. dispersal events", side = 3, cex = 2.5, line = 1.5)
# Add counts
text(x = 0.42, y = 0.54, labels = contingency_table[1,1], adj = 0.5, cex = 5)
text(x = 0.42, y = 0.03, labels = contingency_table[1,2], adj = 0.5, cex = 2)
text(x = 0.804, y = 0.54, labels = contingency_table[2,1], adj = 0.5, cex = 1)
text(x = 0.804, y = 0.01, labels = contingency_table[2,2], adj = 0.5, cex = 1)

# Add Fisher's exact test results
text(x = 0.05, y = 0.90, labels = substitute(paste(bold("Fisher's exact test"))), adj = 0, cex = 1.5)
text(x = 0.05, y = 0.83, labels = paste0("OR = ", round(binary_rate_shifts_vs_dispersal_fisher_test$estimate, 3)),
     adj = 0, cex = 1.5)
text(x = 0.05, y = 0.77, labels = paste0("CI 95% = [", round(binary_rate_shifts_vs_dispersal_fisher_test$conf.int[1], 3),
                                         " - ", round(binary_rate_shifts_vs_dispersal_fisher_test$conf.int[2], 3), "]"),
     adj = 0, cex = 1.5)
text(x = 0.05, y = 0.71, labels = paste0("p = ", round(binary_rate_shifts_vs_dispersal_fisher_test$p.value, 3)),
     adj = 0, cex = 1.5)

# Add Khi² test results
text(x = 0.05, y = 0.36, labels = substitute(paste(bold("Khi² test"))), adj = 0, cex = 1.5)
text(x = 0.05, y = 0.29, labels = paste0("Khi² = ", round(binary_rate_shifts_vs_dispersal_chisq_test$statistic, 3)),
     adj = 0, cex = 1.5)
text(x = 0.05, y = 0.23, labels = paste0("Q95% = ", round(khi_Q95, 3)),
     adj = 0, cex = 1.5)
text(x = 0.05, y = 0.17, labels = paste0("p = ", round(binary_rate_shifts_vs_dispersal_chisq_test$p.value, 3)),
     adj = 0, cex = 1.5)


dev.off()


## Test correlation between PP of rate shift vs. PP of dispersal event?
  # Which test stat to use to account for the covariation in the data due to the phylogenetic structure?
  # Use distance matrices with MRM ? rate_dist ~ disp_dist + phylo_dist 
  # Does not work for PP, but maybe for trait data ? 




