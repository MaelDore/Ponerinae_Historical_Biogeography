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
library(RColorBrewer)
library(ggtree)
library(parallel)
library(raster)
library(vcd)  # For association plots


##### 1/ Load input files ####

### 1.1/ Inform path to BAMM folder ####

BAMM_path <- "./software/bamm-2.5.0/"

### 1.2/ Load phylogeny and check if it is valid

Ponerinae_phylogeny_1534t <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_short_names.rds")

# Check if ultrametic
is.ultrametric(Ponerinae_phylogeny_1534t)
# Check if fully resolved
is.binary(Ponerinae_phylogeny_1534t)
# Check that all branch have positive length
min(Ponerinae_phylogeny_1534t$edge.length) > 0 ; min(Ponerinae_phylogeny_1534t$edge.length)

### 1.3/ Export phylogeny as .tree file for the analyses
write.tree(phy = Ponerinae_phylogeny_1534t, file = paste0(BAMM_path, "Ponerinae_phylogeny.tree"))
phy_path <- paste0(BAMM_path, "Ponerinae_phylogeny.tree")

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
numberOfGenerations = format(10000000, scientific = F)
numberOfGenerations = format(100000, scientific = F) # For the test run
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
mcmcWriteFreq <- as.numeric(numberOfGenerations) / 2000 
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
printFreq <- 1000 # Print status every 10^3 generations for short test runs
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
run_prefix_for_output_files <- "BAMM_Ponerinae_test_run"
run_prefix_for_output_files <- "BAMM_Ponerinae"
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
BD_fit <- phytools::fit.bd(tree = Ponerinae_phylogeny_1534t)

# Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
  # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
lambdaInit0 <- 0.032
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
writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt"))
writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt"))


##### 3/ Run BAMM with calls to command lines from within r using system() #####

### 3.1/ Run BAMM ####

?system

# Version
system(paste0(BAMM_path,"bamm --version"))

# Test run
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt"))

# Full run
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt"))

### Outputs
# The run_info.txt file, containing a summary of your parameters/settings
# An mcmc_log.txt containing raw MCMC information useful in diagnosing convergence
# An event_data.txt file containing all of evolutionary rate parameters and their topological mappings
# A chain_swap.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).

### 3.2/ Clean outputs (move them to a dedicated folder) ####

BAMM_output_folder_path <- "./outputs/BAMM/BAMM_Ponerinae_test_run/"
BAMM_output_folder_path <- "./outputs/BAMM/BAMM_Ponerinae/"

# Detect output files
output_files_path <- list.files(path = "./", pattern = "BAMM_")
# Move output files to dedicated folder
file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "my_phy.tree")) # Phylo file
file.rename(from = paste0(BAMM_path, "my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "my_div_control_file.txt")) # Control file


### 3.3/ Run multiple runs with different prior for the LAMBDA controlling the expected number of shifts ####

expectedNumberOfShifts_range <- c(1, 5, 10, 20)

## Load control files
# my_div_control_file <- readLines(con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_my_div_control_file.txt")) # Test run
my_div_control_file <- readLines(con = paste0(BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt")) # Full run

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
  run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i) # Full run
  outName_line <- which(str_detect(string = my_div_control_file_i, pattern = "outName = "))[1]
  my_div_control_file_i[outName_line] <- paste0("outName = ", run_prefix_for_output_files)
  
  # Write control file
  # writeLines(text = my_div_control_file_i, con = paste0(BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Test run
  writeLines(text = my_div_control_file_i, con = paste0(BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full run
  
  # Lauch BAMM run
  # system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Test run
  system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Full run
  
  # Set output folder
  # BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_test_runs/BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i, "/") # Test run
  BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i, "/") # Full run
  dir.create(path = BAMM_output_folder_path, recursive = T, showWarnings = F)
  
  # Detect output files
  output_files_path <- list.files(path = "./", pattern = run_prefix_for_output_files)
  # Move output files to dedicated folder
  file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
  file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "Ponerinae_phylogeny.tree")) # Phylo file
  # file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Control file for Test run
  file.rename(from = paste0(BAMM_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i,"_my_div_control_file.txt")) # Control file for Full run
  
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
  BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i, "/") # Full run
  
  # Extract prefix of output files
  # run_prefix_for_output_files <- paste0("BAMM_Ponerinae_test_run_nbshifts",expectedNumberOfShifts_i) # Test run
  run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",expectedNumberOfShifts_i) # Full run
  
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
pdf(file = "./outputs/BAMM/MCMC_traces_per_expected_nb_of_shifts.pdf", width = 10, height = 6)

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
BAMM_output_folder_path <- paste0("./outputs/BAMM/BAMM_Ponerinae_full_runs/BAMM_Ponerinae_nbshifts",selected_expectedNumberofShifts, "/") # Full run

# Extract prefix of output files
# run_prefix_for_output_files <- paste0("BAMM_Ponerinae_test_run_nbshifts",selected_expectedNumberofShifts) # Test run
run_prefix_for_output_files <- paste0("BAMM_Ponerinae_nbshifts",selected_expectedNumberofShifts) # Full run

# Load the phylogeny
phy <- read.tree(paste0(BAMM_output_folder_path, "Ponerinae_phylogeny.tree"))

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
saveRDS(BAMM_posterior_samples_data, file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")

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
BAMM_posterior_samples_data$eventBranchSegs[[1]] # Same but with matrix including tipward node ID and begin/end ages of the branches
BAMM_posterior_samples_data$tipStates[[1]] # Integer vectors of regime membership per tips
BAMM_posterior_samples_data$tipLambda[[1]] # Integer vectors of final speciation rates at tips = current speciation rates
BAMM_posterior_samples_data$meanTipLambda # Mean current tip speciation rates across all posterior configurations. Better to use the median and 95% HPD
BAMM_posterior_samples_data$tipMu[[1]] # Integer vectors of final extinction rates at tips = current extinction rates
BAMM_posterior_samples_data$tipLambda[[1]] - BAMM_posterior_samples_data$tipMu[[1]] # Integer vectors of final diversification  rates at tips = current diversification rates
BAMM_posterior_samples_data$meanTipMu # Mean current tip extinction rates across all posterior configurations. Better to use the median and 95% HPD


##### 6/ Plot regime shifts #####

# Marginal probabilities. If high for a sequence of branches, good practice could be to look for in the posterior samples for the joint probability of a shift happening on 0/1/2/… of these branches in a single sample.
# Goal: avoid misinterpreting those independent marginal probabilities has evidence for a sequence of shifts as what it actually supports is the likely presence of one shift among those sequence of branches, whit a bit of uncertainty on the exact location among those successive branches

# Caution: showing all marginal PP of shift location may be biases since shift location may not be independent
# The presence of one shift may be associated with the absence of another. Be careful with the interpretation of marginal PP.
# Try to investigate joint probabilities of some events.

# Load the BAMM posterior samples object
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")

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
# 930 configurations needed to reach 95% in cumulative probability!
# MAP with only PP of only 0.3% 

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
table(branch_marginal_odd_ratios$edge.length > 5)
table(branch_marginal_odd_ratios$edge.length > 100)

plot(branch_marg_posterior_probs) # Posterior probabilities
plot(branch_prior_probs) # Prior probabilities
plot(branch_marginal_odd_ratios) # Odd-ratios = Posterior probabilities / Prior probabilities

# Save prior, posterior, and odd-ratio in a single df
edge_rate_shift_probs_df <- data.frame(edge_ID = 1:nrow(branch_marg_posterior_probs$edge),
                                       edge_length = phy$edge.length,
                                       marg_posterior_probs = branch_marg_posterior_probs$edge.length,
                                       prior_probs = branch_prior_probs$edge.length,
                                       marginal_odd_ratios = branch_marginal_odd_ratios$edge.length)
saveRDS(edge_rate_shift_probs_df, file = "./outputs/BAMM/edge_rate_shift_probs_df.rds")

### 6.4/ Map bioregion membership on phylogeny scaled by shift probability ####

# Load treedata with bioregion membership information
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")

# # Assign a unique bioregion per edge
# bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
# Bioregion_binary_df <- Ponerinae_phylogeny_1534t_treedata@extraInfo[, c("mean_state_Afrotropics", "mean_state_Australasia", "mean_state_Eastern Palearctic", "mean_state_Indomalaya", "mean_state_Nearctic", "mean_state_Neotropics", "mean_state_Western Palearctic")]
# Bioregion_binary_df <- bioregion_names[unlist(apply(X = Bioregion_binary_df, MARGIN = 1, FUN = which.max))]
# Bioregion_binary_df <- as.data.frame(cbind(Ponerinae_phylogeny_1534t_treedata@extraInfo$node, Bioregion_binary_df))
# names(Bioregion_binary_df) <- c("node", "Bioregion")
# Bioregion_binary_df$node <- as.numeric(Bioregion_binary_df$node)
# Bioregion_binary_df$Bioregion <- factor(x = Bioregion_binary_df$Bioregion, levels = bioregion_names, labels = bioregion_names)
# 
# Ponerinae_phylogeny_1534t_treedata@extraInfo <- left_join(Ponerinae_phylogeny_1534t_treedata@extraInfo, y = Bioregion_binary_df)
# 
# # Save updated treedata with unique Bioregion membership (the most likely one)
# saveRDS(object = Ponerinae_phylogeny_1534t_treedata, file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")

# Load Genus-groups metadata
all_genus_groups_metadata_df <- readRDS(file = "./outputs/Grafting_missing_taxa/all_genus_groups_metadata_df.rds")

# Set angles and offset of genus-groups labels
genus_groups_angles <- c(275, 300, 345, 2, 48, 12, 357)
genus_groups_offsets <- c(15.0, 15.0, 15.0, 28.0, 15.0, 15.0, 15.0)
names(genus_groups_angles) <- all_genus_groups_metadata_df$group_name
names(genus_groups_offsets) <- all_genus_groups_metadata_df$group_name

# Set color scheme for Genus groups
genus_groups_colors <- tmaptools::get_brewer_pal("Spectral", n = nrow(all_genus_groups_metadata_df) + 2, plot = F)[-c(4:5)]

# Set color scheme for areas/bioregions (Use the BSM color scheme)
# colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
# areas_list <- c("A", "U", "E", "I", "R", "N", "W")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
bioregion_names <- c("Afrotropics", "Australasia", "Eastern Palearctic", "Indomalaya", "Nearctic", "Neotropics",  "Western Palearctic")
names(colors_list_for_areas) <- bioregion_names

# Set color scheme for tips
nb_tips <- length(Ponerinae_phylogeny_scaled_PP_shifts@phylo$tip.label)
tips_color <- c("black", "red")[Ponerinae_phylogeny_scaled_PP_shifts@extraInfo$missing_node + 1]
tips_color <- tips_color[1:nb_tips]

# Initiate new treedata
Ponerinae_phylogeny_scaled_PP_shifts <- Ponerinae_phylogeny_1534t_treedata 

# Update branch length for marginal posterior probability of rate shift
Ponerinae_phylogeny_scaled_PP_shifts@phylo$edge.length <- branch_marg_posterior_probs$edge.length

# Plot PDF to aggregate all bioregions
pdf(file = paste0("./outputs/BAMM/Ponerinae_phylogeny_scaled_PP_shifts.pdf"), height = 16, width = 17)


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
                    offset = c(-2.6, -2.6, -1.6, -1.1, -0.9, -1.1, 0.0)[i],
                    offset.text = c(0.10, 0.03, 0.05, 0.03, 0.05, 0.10, 0.05)[i],
                    barsize = 4)    # Width of the bar
}
  
# Add transparent tree with Bioregion colors (needed for the legend)
Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
  geom_tree(data = Ponerinae_phylogeny_scaled_PP_shifts,
            mapping = aes(col = Bioregion),
            alpha = 0.0, linewidth = 5.0,
            layout = "rectangular") +

  # Add Bioregion legend
  scale_color_manual("Bioregion", breaks = bioregion_names, values = colors_list_for_areas) +
  guides(color = guide_legend(override.aes = list(alpha = 1.0)))

# Adapt margins
Ponerinae_phylogeny_scaled_PP_shifts_plot <- Ponerinae_phylogeny_scaled_PP_shifts_plot +
  theme(plot.margin = unit(c(0, 0, 20, 0), "mm"), # trbl?
        legend.position = c(0.85, 0.2),
        legend.title = element_text(size = 24, face = "bold", margin = margin(b = 25)),
        legend.text = element_text(size = 18, margin = margin(l = 15)),
        legend.key.spacing.y = unit(1.5, 'line'),
        plot.title = element_text(face = "bold", colour = "black", size = 30, hjust = 0.5,
                                  angle = 0, vjust = -10),
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

# Only found in 3 posterior samples!
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
saveRDS(object = MSC_tree, file = "./outputs/BAMM/MSC_tree.rds")

### 7.3/ Extract metadata on rate shifts in the MSC ####

root_age <- max(phytools::nodeHeights(phy))

# Extract event table
MSC_shifts_df <- BAMM_posterior_samples_data$eventData[[MSC_detection$sampleindex]][, -6]
names(MSC_shifts_df) <- c("tipward_node", "start_time_since_root", "lambda0", "alpha", "mu", "regime_ID")

# Compute start age
MSC_shifts_df$start_age <- root_age - MSC_shifts_df$start_time_since_root

# Get associated branch ID
MSC_shifts_df$edge_ID <- match(x = MSC_shifts_df$tipward_node, table = phy$edge[, 2])

# Get the associated probabilities/odd-ratios
MSC_shifts_df$branch_prior_prob <- branch_prior_probs$edge.length[MSC_shifts_df$edge_ID]
MSC_shifts_df$branch_marg_posterior_prob <- branch_marg_posterior_probs$edge.length[MSC_shifts_df$edge_ID]
MSC_shifts_df$branch_odd_ratio <- branch_marginal_odd_ratios$edge.length[MSC_shifts_df$edge_ID]

# Identify dynamics
MSC_shifts_df$trend <- c("decrease", "increase")[as.numeric(MSC_shifts_df$alpha > 0)+1]
MSC_shifts_df$trend_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$alpha > 0)+1]

# Extract mean branch rates
div_rates_tree <- getMeanBranchLengthTree(BAMM_posterior_samples_data,
                                          rate = "ndr") # Net Diversification Rates

# Identify shift types
MSC_shifts_df$mean_branch_rate <- div_rates_tree$phy$edge.length[MSC_shifts_df$edge_ID]
MSC_shifts_df$rootward_node_ID <- phy$edge[MSC_shifts_df$edge_ID, 1]
MSC_shifts_df$previous_branch_ID <- match(x = MSC_shifts_df$rootward_node_ID, table = phy$edge[, 2])
MSC_shifts_df$previous_branch_mean_rate <- div_rates_tree$phy$edge.length[MSC_shifts_df$previous_branch_ID]
MSC_shifts_df$rate_shift <- MSC_shifts_df$mean_branch_rate - MSC_shifts_df$previous_branch_mean_rate
MSC_shifts_df$rate_shift_type <- c("decrease", "increase")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]
MSC_shifts_df$rate_shift_type_col <- c("dodgerblue", "brown1")[as.numeric(MSC_shifts_df$rate_shift > 0)+1]
  
# Identify clades/paraphyletic groups associated with regimes

# Plot phylogeny with rate shift location
pdf("./outputs/BAMM/MSC_tree_with_labels_rect.pdf", width = 20, height = 200)
plot.phylo(phy)
#nodelabels(node = MSC_shift_nodes, pch = 21, col = "black", bg = "red", cex = 1.5)
edgelabels(edge = MSC_shifts_df$edge_ID, pch = 21, col = "black", bg = MSC_shifts_df$rate_shift_type_col[-1], cex = MSC_shifts_df$branch_marg_posterior_prob * 3)
dev.off()

# Provide name to the regimes
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1535] <- "Root_process"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2668] <- "Ponera_group"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2615] <- "Hypoponera_fast_subgroup1"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2491] <- "Hypoponera_fast_subgroup2"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2426] <- "Hypoponera_fast_subgroup3"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 2352] <- "Plectroctena_Genus"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1537] <- "PlOHH"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1763] <- "Myopias_Leptogenys"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1767] <- "Leptogenys_subgroup"
MSC_shifts_df$regime_group_name[MSC_shifts_df$tipward_node == 1705] <- "Odontomachus_subgroup"


# Save metadata on rate shifts in the MSC tree
saveRDS(object = MSC_shifts_df, file = "./outputs/BAMM/MSC_shifts_df.rds")




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
saveRDS(div_rates_tree, file = "./outputs/BAMM/div_rates_tree.rds")


### 8.3/ Plot mean sliding rates + credible shift locations ####

## Good practice wood be to use median rather than mean
# Median can be computed from BAMM_output$eventData or using BAMMtools::dtRates to compute discretize rates

# Use the MSC configuration for plotting the shifts
# read metadata on rate shifts in the MSC tree
MSC_shifts_df <- readRDS(file = "./outputs/BAMM/MSC_shifts_df.rds")

## Plot histogram of different color palette settings for selection

pdf(file = "./outputs/BAMM/Histograms_branch_rates.pdf", height = 12, width = 6)

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

pdf(file = "./outputs/BAMM/Ponerinae_diversity_dynamics_rect.pdf", height = 10, width = 6)

par(mfrow = c(1,1), xpd = TRUE)

# tau = length of segments used to discretize the continuous rates
# tau = 0.001 => each segment is 1/1000 of the tree height (root age)

# Quantile method to discretize rates in categories of equal frequencies
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       method = "phylogram",
                                       tau = 0.001, breaksmethod = 'quantile',
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

pdf(file = "./outputs/BAMM/Ponerinae_diversity_dynamics_polar.pdf", height = 7, width = 6)

par(mfrow = c(1,1), xpd = TRUE)

# tau = length of segments used to discretize the continuous rates
# tau = 0.001 => each segment is 1/1000 of the tree height (root age)

# Quantile method to discretize rates in categories of equal frequencies
branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
                                       method = "polar",
                                       tau = 0.001, breaksmethod = 'quantile',
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


## Groups to focus on: 
 # Pachycondyla/Platythyrea: excluded: 1537, 2668
 # Ponera-group: included: 2668
 # Myopias/Leptogenys: included: 1763
 # Other taxa: included: 1537 + excluded: 1763 
 # Overall

max(div_rates_tree$phy$edge.length)

### 9.1/ Plot of rates through time with fuzzy intervals: median + all trajectories ####

pdf(file = "./outputs/BAMM/BAMM_LLT_per_clades_original_plot_fuzzy.pdf", height = 10, width = 12)

par(mfrow = c(2, 2))

root_age <- max(phytools::nodeHeights(phy))

# A/ All branches
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01),
                    intervalCol = "black", avgCol = "black",
                    start.time = root_age,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "All Ponerinae", col = "black", font = 4, cex = 2.0, pos = 4)

# B/ Only Pachycondyla/Platythyrea
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01),
                    intervalCol = "red", avgCol = "red", 
                    start.time = root_age, 
                    node = c(1537, 2668), nodetype = "exclude",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "Pachycondyla/Platythyrea", col = "red", font = 4, cex = 2.0, pos = 4)

# C/ Only Ponera
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "orange", avgCol = "orange",
                    start.time = root_age, 
                    node = 2668,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "Ponera group", col = "orange", font = 4, cex = 2.0, pos = 4)

# D/ Only Myopias/Leptogenys
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = seq(from = 0, to = 1, by = 0.01), 
                    intervalCol = "dodgerblue3", avgCol = "dodgerblue3",
                    start.time = root_age, 
                    node = 1763,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "Myopias/Leptogenys", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)

dev.off()


### 9.2/ Plot of rates through time with strict intervals: mean + 95% CI polygon ####

# Can be partitioned per clades (and per bioregions?)

pdf(file = "./outputs/BAMM/BAMM_LLT_per_clades_original_plot_CI95.pdf", height = 10, width = 12)

par(mfrow = c(2, 2))

root_age <- max(phytools::nodeHeights(phy))

# A/ All branches
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "black", avgCol = "black",
                    start.time = root_age,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "All Ponerinae", col = "black", font = 4, cex = 2.0, pos = 4)

# B/ Only Pachycondyla/Platythyrea
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2,
                    intervalCol = "red", avgCol = "red", 
                    start.time = root_age, 
                    node = c(1537, 2668), nodetype = "exclude",
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 100, y = 0.4, label = "Pachycondyla/Platythyrea", col = "red", font = 4, cex = 2.0, pos = 4)

# C/ Only Ponera
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2, 
                    intervalCol = "orange", avgCol = "orange",
                    start.time = root_age, 
                    node = 2668,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "Ponera group", col = "orange", font = 4, cex = 2.0, pos = 4)

# D/ Only Myopias/Leptogenys
plotRateThroughTime(BAMM_posterior_samples_data,
                    intervals = c(0.025, 0.975),
                    opacity = 0.2, 
                    intervalCol = "dodgerblue3", avgCol = "dodgerblue3",
                    start.time = root_age, 
                    node = 1763,
                    ylim = c(0, 0.5),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 100, y = 0.4, label = "Myopias/Leptogenys", col = "dodgerblue3", font = 4, cex = 2.0, pos = 4)

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

# A/ Extract for the whole tree
all_Ponerinae_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100)

all_Ponerinae_rates_through_time_matrix$net_div <- all_Ponerinae_rates_through_time_matrix$lambda - all_Ponerinae_rates_through_time_matrix$mu

str(all_Ponerinae_rates_through_time_matrix, max.level = 1)

dim(all_Ponerinae_rates_through_time_matrix$lambda) # Speciation rates per posterior sample (raw) x time window (columns)
dim(all_Ponerinae_rates_through_time_matrix$mu) # Extinction rates per posterior sample (raw) x time window (columns)
all_Ponerinae_rates_through_time_matrix$times # Mean time of the time windows. Names = age. Values = time since root.

all_Ponerinae_times <- round(as.numeric(names(all_Ponerinae_rates_through_time_matrix$times)), 1)

# B/ Extract for Only Pachycondyla/Platythyrea
Pachycondyla_Platythyrea_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                               node = c(1537, 2668), nodetype = "exclude")

Pachycondyla_Platythyrea_rates_through_time_matrix$net_div <- Pachycondyla_Platythyrea_rates_through_time_matrix$lambda - Pachycondyla_Platythyrea_rates_through_time_matrix$mu


# C/ Extract for Only Ponera
Ponera_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                                               node = 2668, nodetype = "include")

Ponera_rates_through_time_matrix$net_div <- Ponera_rates_through_time_matrix$lambda - Ponera_rates_through_time_matrix$mu

# Adjust to overall time scale
Ponera_times <- round(as.numeric(names(Ponera_rates_through_time_matrix$times)), 1)
Ponera_times_updated <- sapply(X = Ponera_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Ponera_rates_through_time_matrix_updated <- Ponera_rates_through_time_matrix
Ponera_update_indices <- match(all_Ponerinae_times, Ponera_times_updated)
Ponera_rates_through_time_matrix_updated$lambda <- Ponera_rates_through_time_matrix$lambda[, Ponera_update_indices]
Ponera_rates_through_time_matrix_updated$mu <- Ponera_rates_through_time_matrix$mu[, Ponera_update_indices]
Ponera_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Ponera_rates_through_time_matrix_updated$net_div <- Ponera_rates_through_time_matrix$net_div[, Ponera_update_indices]

# D/ Extract for Only Myopias/Leptogenys
Myopias_Leptogenys_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_posterior_samples_data, nslices = 100,
                                                             node = 1763, nodetype = "include")

Myopias_Leptogenys_rates_through_time_matrix$net_div <- Myopias_Leptogenys_rates_through_time_matrix$lambda - Myopias_Leptogenys_rates_through_time_matrix$mu

# Adjust to overall time scale
Myopias_Leptogenys_times <- round(as.numeric(names(Myopias_Leptogenys_rates_through_time_matrix$times)), 1)
Myopias_Leptogenys_times_updated <- sapply(X = Myopias_Leptogenys_times, FUN = find_closest_value, y = all_Ponerinae_times)

# Update matrix data
Myopias_Leptogenys_rates_through_time_matrix_updated <- Myopias_Leptogenys_rates_through_time_matrix
Myopias_Leptogenys_update_indices <- match(all_Ponerinae_times, Myopias_Leptogenys_times_updated)
Myopias_Leptogenys_rates_through_time_matrix_updated$lambda <- Myopias_Leptogenys_rates_through_time_matrix$lambda[, Myopias_Leptogenys_update_indices]
Myopias_Leptogenys_rates_through_time_matrix_updated$mu <- Myopias_Leptogenys_rates_through_time_matrix$mu[, Myopias_Leptogenys_update_indices]
Myopias_Leptogenys_rates_through_time_matrix_updated$times <- all_Ponerinae_times
Myopias_Leptogenys_rates_through_time_matrix_updated$net_div <- Myopias_Leptogenys_rates_through_time_matrix$net_div[, Myopias_Leptogenys_update_indices]


## Aggregate LTT data in an array
LTT_per_clades_array <- array(data = NA,
                              dim = c(dim(all_Ponerinae_rates_through_time_matrix$net_div), 4),
                              dimnames = list(sample = paste0("sample_",1:dim(all_Ponerinae_rates_through_time_matrix$net_div)[1]),
                                              time = names(all_Ponerinae_rates_through_time_matrix$times),
                                              group = c("All_Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys")))

LTT_per_clades_array[,,1] <- all_Ponerinae_rates_through_time_matrix$net_div                           
LTT_per_clades_array[,,2] <- Pachycondyla_Platythyrea_rates_through_time_matrix$net_div
LTT_per_clades_array[,,3] <- Ponera_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale
LTT_per_clades_array[,,4] <- Myopias_Leptogenys_rates_through_time_matrix_updated$net_div # Use the updated version with matching time scale

# Aggregate LTT data in a unique melted data.frame
LTT_per_clades_melted_df <- reshape2::melt(LTT_per_clades_array)
LTT_per_clades_melted_df$time <- round(as.numeric(LTT_per_clades_melted_df$time), 1)

# Rename rates data
LTT_per_clades_melted_df <- LTT_per_clades_melted_df %>% 
  rename(net_div_rates = value)

# Order groups
LTT_per_clades_melted_df$group <- factor(LTT_per_clades_melted_df$group, levels = c("All_Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys"), labels = c("All Ponerinae", "Pachycondyla/Platythyrea", "Ponera", "Myopias/Leptogenys"))

## Save LTT data per clades
saveRDS(LTT_per_clades_melted_df, file = "./outputs/BAMM/LTT_per_clades_melted_df.rds")


### 9.4/ Plot LTT per clades: ggplot version ####

## Compute medians
LTT_per_clades_melted_df_median <- LTT_per_clades_melted_df %>% 
  group_by(time, group) %>% 
  summarize(median_rates = median(net_div_rates)) %>%
  ungroup()

## Compute quantiles
LTT_per_clades_melted_df_quantiles <- LTT_per_clades_melted_df %>% 
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
  mutate(points_ID = c(1:100, 200:101)) %>%
  # Reorder by points ID
  arrange(group, quantile_ID, points_ID) %>%
  filter(!is.na(quant_rates)) %>%
  mutate(points_ID = row_number()) %>%
  ungroup()


## Generate plot using lines
pdf(file = "./outputs/BAMM/BAMM_LLT_per_clades_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
  ) + 
  
  # Adjust color scheme and legend
  scale_color_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) +
  
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
pdf(file = "./outputs/BAMM/BAMM_LLT_per_clades_ggplot_fuzzy_quantiles.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
  ) + 
  
  # Adjust fill scheme and legend
  scale_fill_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) +
  
  # Adjust color scheme and legend
  scale_color_manual("Groups", values = c("black", "red", "orange", "dodgerblue3")) +
  
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
saveRDS(BAMM_div_rates_trees, file = "./outputs/BAMM/BAMM_div_rates_trees.rds")

### 10.2/ Compute weighted mean rates along phylogenies ####

# Easy version = use mean rates and membership weight per branches
# Tricky version = get rates and membership per time segment as in densityMap


## Load BAMM trees scaled by div rates
BAMM_div_rates_trees <- readRDS(file = "./outputs/BAMM/BAMM_div_rates_trees.rds")

## Load treedata with bioregion membership information
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")

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
saveRDS(weighted_mean_per_bioregions_all_BAMM_phylo_array, file = "./outputs/BAMM/weighted_mean_per_bioregions_all_BAMM_phylo_array.rds")


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
saveRDS(object = weighted_mean_per_bioregions_all_BAMM_phylo_melted, file = "./outputs/BAMM/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")


### 10.4/ Compute median and quantiles rates across posterior samples ####

## Compute medians
weighted_mean_per_bioregions_median <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
  group_by(time, bioregion) %>% 
  summarize(median_rates = median(weighted_mean, na.rm = T)) %>%
  ungroup()

## Compute quantiles
weighted_mean_per_bioregions_quantiles <- weighted_mean_per_bioregions_all_BAMM_phylo_melted %>% 
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
  mutate(points_ID = c(1:100, 200:101)) %>%
  # Reorder by points ID
  arrange(bioregion, quantile_ID, points_ID) %>%
  filter(!is.na(quant_rates)) %>%
  mutate(points_ID = row_number()) %>%
  ungroup()

## Save medians and quantiles
saveRDS(object = weighted_mean_per_bioregions_median, "./outputs/BAMM/weighted_mean_per_bioregions_median.rds")
saveRDS(object = weighted_mean_per_bioregions_quantiles, "./outputs/BAMM/weighted_mean_per_bioregions_quantiles.rds")

  
### 10.5/ Plot rates through time per bioregions ####

## Set bioregion color scheme
# colors_list_for_states <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_states.rds")
# areas_list <- c("A", "U", "E", "I", "R", "N", "W")
# colors_list_for_areas <- colors_list_for_states[areas_list]
colors_list_for_areas <- readRDS(file = "./outputs/BSM/BSM_maps/colors_list_for_areas_light.rds")
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_per_bioregions_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_per_bioregions_fuzzy_quantiles.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_OW_vs_NW_fuzzy_quantiles.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_Trop_vs_Temp_ggplot_fuzzy_lines.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
pdf(file = "./outputs/BAMM/BAMM_ratesTT_Trop_vs_Temp_fuzzy_quantiles.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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
  # 4/ Compare observed statistic to its null distribution on the permutated tip values

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
  select(Taxa, Bioregion) %>%
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
  select(Taxa, Bioregion)
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

# Assign tip names to tip data
Ponerinae_bioregion_tip_data_all_bioregions <- Taxa_bioregions_tips_data$Bioregion
names(Ponerinae_bioregion_tip_data_all_bioregions) <- Taxa_bioregions_tips_data$Taxa

# Reorder as in phylogeny
Ponerinae_bioregion_tip_data_all_bioregions <- Ponerinae_bioregion_tip_data_all_bioregions[match(phy$tip.label, names(Ponerinae_bioregion_tip_data_all_bioregions))]
identical(names(Ponerinae_bioregion_tip_data_all_bioregions), phy$tip.label)

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

## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_all_bioregions$test_stats <- STRAPP_test_all_bioregions$obs.corr - STRAPP_test_all_bioregions$null
summary(STRAPP_test_all_bioregions$test_stats)
table(STRAPP_test_all_bioregions$test_stats > 0)
STRAPP_test_all_bioregions$stat_median <- median(STRAPP_test_all_bioregions$test_stats)
STRAPP_test_all_bioregions$stat_Q5 <- quantile(STRAPP_test_all_bioregions$test_stats, p = 0.05)

pdf(file = "./outputs/BAMM/STRAPP_test_all_bioregions_hist.pdf", height = 6, width = 8)

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

## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_OW_vs_NW$test_stats <- STRAPP_test_OW_vs_NW$null - STRAPP_test_OW_vs_NW$obs.corr
summary(STRAPP_test_OW_vs_NW$test_stats)
table(STRAPP_test_OW_vs_NW$test_stats > 0)
STRAPP_test_OW_vs_NW$stat_median <- median(STRAPP_test_OW_vs_NW$test_stats)
STRAPP_test_OW_vs_NW$stat_Q5 <- quantile(STRAPP_test_OW_vs_NW$test_stats, p = 0.05)

pdf(file = "./outputs/BAMM/STRAPP_test_OW_vs_NW_hist.pdf", height = 6, width = 8)

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

## Plot histogram of the test

# Compute info for histogram 
STRAPP_test_Trop_vs_Temp$test_stats <- STRAPP_test_Trop_vs_Temp$null - STRAPP_test_Trop_vs_Temp$obs.corr
summary(STRAPP_test_Trop_vs_Temp$test_stats)
table(STRAPP_test_Trop_vs_Temp$test_stats > 0)
STRAPP_test_Trop_vs_Temp$stat_median <- median(STRAPP_test_Trop_vs_Temp$test_stats)
STRAPP_test_Trop_vs_Temp$stat_Q5 <- quantile(STRAPP_test_Trop_vs_Temp$test_stats, p = 0.05)

pdf(file = "./outputs/BAMM/STRAPP_test_Trop_vs_Temp_hist.pdf", height = 6, width = 8)

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
rect(xleft = x_min_max[1], x_min_max[1] + x_range * 0.35,
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

### 12.1/ Compute density maps of aggregated bioregions

# Load simmaps with only unique areas
DEC_J_simmaps_unique_areas <- readRDS(file = "./outputs/BSM/DEC_J_simmaps_unique_areas.rds")

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
DEC_J_density_map_OW_vs_NW <- setMap(DEC_J_density_map, col_scale)
DEC_J_density_map_OW_vs_NW$states <- c("Old World", "New World")

## Plot result
plot(DEC_J_density_map_OW_vs_NW)

## Save the resulting Density map with updated color gradient
saveRDS(object = DEC_J_density_map_OW_vs_NW, file = paste0("./outputs/Density_maps/DEC_J_density_map_OW_vs_NW.rds"))


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
saveRDS(object = DEC_J_density_map_Trop_vs_Temp, file = paste0("./outputs/Density_maps/DEC_J_density_map_Trop_vs_Temp.rds"))


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
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")

## Load the binary density map for OW vs. NW
DEC_J_density_map_OW_vs_NW <- readRDS(file = paste0("./outputs/Density_maps/DEC_J_density_map_OW_vs_NW.rds"))


# Load the custom STRAPP test function
source("./functions/run_STRAPP_test.R")

# Extract time scale from LTT data
LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/LTT_per_clades_melted_df.rds")
time_scale <- unique(LTT_per_clades_melted_df$time)
time_scale <- time_scale[order(time_scale, decreasing = T)]

## Set seed for reproductibility of permutation tests
set.seed(seed = 12345)

## Loop per time scale
STRAPP_data_OW_vs_NW_all_time <- list()
STRAPP_tests_OW_vs_NW_all_time <- list()
for (i in seq_along(time_scale))
# for (i in 6:length(time_scale))
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
  saveRDS(STRAPP_data_OW_vs_NW_all_time, file = "./outputs/BAMM/STRAPP_data_OW_vs_NW_all_time.rds")

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
  saveRDS(STRAPP_tests_OW_vs_NW_all_time, file = "./outputs/BAMM/STRAPP_tests_OW_vs_NW_all_time.rds")
  
  ## Print progress
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - STRAPP test ran for time = ",time_i," My - n°", i, "/", length(time_scale),"\n"))
  }
}


## 12.3.3/ Plot evolution of p-value in time ####

# # Load input data for STRAPP test
# STRAPP_data_OW_vs_NW_all_time <- readRDS(file = "./outputs/BAMM/STRAPP_data_OW_vs_NW_all_time.rds")

# Load STRAPP test outputs
STRAPP_tests_OW_vs_NW_all_time <- readRDS(file = "./outputs/BAMM/STRAPP_tests_OW_vs_NW_all_time.rds")

# Extract p-value and time
STRAPP_tests_OW_vs_NW_df <- data.frame(time = unlist(lapply(STRAPP_tests_OW_vs_NW_all_time, FUN = function (x) { x$focal_time }) ),
                                       p_value = unlist(lapply(STRAPP_tests_OW_vs_NW_all_time, FUN = function (x) { x$p.value } )))
STRAPP_tests_OW_vs_NW_df

# Save STRAPP tests through evolutionary time df for OW vs NW
saveRDS(STRAPP_tests_OW_vs_NW_df, file = "./outputs/BAMM/STRAPP_tests_OW_vs_NW_df.rds")

# Extract root age
root_age <- time_scale[1]

# Prepare data for significance area
significance_area_poly_df <- data.frame(p_value = c(0.00, 0.00, 0.05, 0.05),
                                        time = c(root_age, 0, 0, root_age), 
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
                                            time = rep(x = c(root_age, 0, 0, root_age), times = nb_poly), 
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
pdf(file = "./outputs/BAMM/STRAPP_tests_OW_vs_NW_pvalTT.pdf", height = 6, width = 8)

STRAPP_tests_OW_vs_NW_plot <- ggplot(data = STRAPP_tests_OW_vs_NW_df,
                                     mapping = aes(y = p_value, x = time)) +
  
  # Plot significance line = 0.05
  geom_hline(yintercept = 0.05, linewidth = 1.0, linetype = "dashed") +
  
  # Plot significance area
  geom_polygon(data = significance_area_poly_df,
               mapping = aes(y = p_value, x = time, group = poly_ID),
               fill = "limegreen", col = NA,
               alpha = 1.00,
               linewidth = 1.0) +
  
  # Plot significance gradient
  geom_polygon(data = significance_area_gradient_df,
               mapping = aes(y = p_value, x = time,
                             group = poly_ID, alpha = poly_ID),
               fill = "limegreen", col = NA,
               linewidth = 1.0, show.legend = F) +
  
  # Plot vertical lines of significance
  geom_vline(xintercept = time_signif,
             col = "black", lty = 2,
             linewidth = 1.0) +
  
  # Add significance legend
  annotate(geom = "text", x = 21.5, y = 0.5,
           label = "Significant\nperiod",
           size = 5.0, hjust = 0.5,
           fontface = "bold",
           color = "limegreen") +
  
  # Plot mean lines
  geom_line(col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("STRAPP test\nDifference in net div. rates through time\nOld World vs. New World")) +

  # Set axes labels
  xlab("Time  [My]") +
  ylab("P-value") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     expand = c(0,0),
                     limits = c(root_age, 0) # Set limits
  ) + 
  
  # Reverse p-value scale
  scale_y_continuous(transform = "reverse",
                     limits = c(1, 0) # Set limits
  ) +
  
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

# Load melted dataframe or Rates through time
weighted_mean_per_bioregions_all_BAMM_phylo_melted <- readRDS(file = "./outputs/BAMM/weighted_mean_per_bioregions_all_BAMM_phylo_melted.rds")
weighted_mean_per_bioregions_median <- readRDS(file = "./outputs/BAMM/weighted_mean_per_bioregions_median.rds")
weighted_mean_per_bioregions_quantiles <- readRDS(file = "./outputs/BAMM/weighted_mean_per_bioregions_quantiles.rds")

# Load STRAPP tests through evolutionary time df for OW vs NW
STRAPP_tests_OW_vs_NW_df <- readRDS(file = "./outputs/BAMM/STRAPP_tests_OW_vs_NW_df.rds")

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
saveRDS(STRAPP_tests_OW_vs_NW_df_thin, "./outputs/BAMM/STRAPP_tests_OW_vs_NW_df_thin.rds")


## Prepare data for significance gradient
# One polygon per x-time
nb_poly <- 0
time_x_data <- c()
max_rates <- 0.155
min_rates <- -0.02

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
  

## Generate plot using lines
pdf(file = "./outputs/BAMM/BAMM_ratesTT_OW_vs_NW_ggplot_fuzzy_lines_with_signif.pdf", height = 8, width = 12)

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
                     limits = c(root_age, 0) # Set limits
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

# Plot
print(ratesT_OW_vs_NW_ggplot_lines_with_signif)

dev.off()


### 12.4/ STRAPP test along evolutionary time for Trop_vs_Temp #### 

# Load the BAMM posterior samples object
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")

## Load the binary density map for Tropics vs Temperate
DEC_J_density_map_Trop_vs_Temp <- readRDS(file = paste0("./outputs/Density_maps/DEC_J_density_map_Trop_vs_Temp.rds"))

# Load the custom STRAPP test function
source("./functions/run_STRAPP_test.R")

# Extract time scale from LTT data
LTT_per_clades_melted_df <- readRDS(file = "./outputs/BAMM/LTT_per_clades_melted_df.rds")
time_scale <- unique(LTT_per_clades_melted_df$time)
time_scale <- time_scale[order(time_scale, decreasing = T)]

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
  saveRDS(STRAPP_data_Trop_vs_Temp_all_time, file = "./outputs/BAMM/STRAPP_data_Trop_vs_Temp_all_time.rds")
  
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
  saveRDS(STRAPP_tests_Trop_vs_Temp_all_time, file = "./outputs/BAMM/STRAPP_tests_Trop_vs_Temp_all_time.rds")
  
  ## Print progress
  if (i %% 1 == 0)
  {
    cat(paste0(Sys.time(), " - STRAPP test ran for time = ",time_i," My - n°", i, "/", length(time_scale),"\n"))
  }
}


## 12.4.3/ Plot evolution of p-value in time ####

# # Load input data for STRAPP test
# STRAPP_data_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/STRAPP_data_Trop_vs_Temp_all_time.rds")

# Load STRAPP test outputs
STRAPP_tests_Trop_vs_Temp_all_time <- readRDS(file = "./outputs/BAMM/STRAPP_tests_Trop_vs_Temp_all_time.rds")

# Extract p-value and time
STRAPP_tests_Trop_vs_Temp_df <- data.frame(time = unlist(lapply(STRAPP_tests_Trop_vs_Temp_all_time, FUN = function (x) { x$focal_time }) ),
                                       p_value = unlist(lapply(STRAPP_tests_Trop_vs_Temp_all_time, FUN = function (x) { x$p.value } )))
STRAPP_tests_Trop_vs_Temp_df

# Save STRAPP tests through evolutionary time df for OW vs NW
saveRDS(STRAPP_tests_Trop_vs_Temp_df, file = "./outputs/BAMM/STRAPP_tests_Trop_vs_Temp_df.rds")

# No period of significance! Min p-value = 0.321 at 13.7 My
STRAPP_tests_Trop_vs_Temp_df[which.min(STRAPP_tests_Trop_vs_Temp_df$p_value), ]

# Extract root age
root_age <- time_scale[1]

## GGplot
pdf(file = "./outputs/BAMM/STRAPP_tests_Trop_vs_Temp_pvalTT.pdf", height = 6, width = 8)

STRAPP_tests_Trop_vs_Temp_plot <- ggplot(data = STRAPP_tests_Trop_vs_Temp_df,
                                     mapping = aes(y = p_value, x = time)) +
  
  # Plot significance line = 0.05
  geom_hline(yintercept = 0.05, linewidth = 1.0, linetype = "dashed") +
  
  # Plot mean lines
  geom_line(col = "red",
            alpha = 1.0,
            linewidth = 2.0) +
  
  # Set plot title +
  ggtitle(label = paste0("STRAPP test\nDifference in net div. rates through time\nTropics vs. Temperate")) +
  
  # Set axes labels
  xlab("Time  [My]") +
  ylab("P-value") +
  
  # Reverse time scale
  scale_x_continuous(transform = "reverse",
                     expand = c(0,0),
                     limits = c(root_age, 0) # Set limits
  ) + 
  
  # Reverse p-value scale
  scale_y_continuous(transform = "reverse",
                     limits = c(1, 0) # Set limits
  ) +

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



##### 13/ Map of current net dispersal rates based on tip rates ####

## Load binary alpha-hull range raster stack
Ponerinae_alpha_hull_stack_WGS84 <- readRDS(file = paste0("./outputs/Species_richness_maps/Ponerinae_species_alpha_hull_stack_WGS84.rds"))

## Load bioregion sf
Bioregions_sf_Bioregions_level <- readRDS(file = "./input_data/geoBoundaries/Bioregions_sf_Bioregions_level.rds")
Bioregions_sf_Bioregions_level_Mollweide <- readRDS(file = "./outputs/Species_richness_maps/Bioregions_sf_Bioregions_level_Mollweide.rds")

## Load the BAMM posterior samples object
BAMM_posterior_samples_data <- readRDS(file = "./outputs/BAMM/BAMM_posterior_samples_data.rds")

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
saveRDS(current_median_tip_rates_df, file = "./outputs/BAMM/current_median_tip_rates_df.rds")


### 13.2/ Compute mean current rates raster ####

# Load terrestrial background
terrestrial_bg_WGS84 <- readRDS(file = "./outputs/Species_richness_maps/terrestrial_bg_WGS84.rds")

# Load current median tip rates
current_median_tip_rates_df <- readRDS(file = "./outputs/BAMM/current_median_tip_rates_df.rds")

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
saveRDS(Ponerinae_net_div_rates_stack_WGS84, file = "./outputs/BAMM/Ponerinae_net_div_rates_stack_WGS84.rds")

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
saveRDS(Ponerinae_mean_net_div_rates_WGS84, file = "./outputs/BAMM/Ponerinae_mean_net_div_rates_WGS84.rds")



### 13.3/ Plot mean current rate geographic map ####

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./outputs/Species_richness_maps/pal_bl_red_Mannion.rds")
pal_bl_red_brewer <- rev(RColorBrewer::brewer.pal(n = 11, "RdYlBu"))
pal_bl_red_brewer_fn <- colorRampPalette(colors = pal_bl_red_brewer)
pal_bl_red_brewer <- pal_bl_red_brewer_fn(n = 400)
blues <- round(seq(from = 1, to = 260, length.out = 99), 0)
reds <- round(seq(from = 261, to = 400, length.out = 100), 0)
pal_bl_red_brewer <- pal_bl_red_brewer[c(blues, reds)]
pal_bl_red_brewer <- c("grey90", pal_bl_red_brewer)

# Contrast raster
source("./functions/contrasting_raster.R")
hist(Ponerinae_mean_net_div_rates_WGS84[])
table(Ponerinae_mean_net_div_rates_WGS84[])
Ponerinae_mean_net_div_rates_WGS84_contrasted <- contrasting_raster(x = Ponerinae_mean_net_div_rates_WGS84, zmin = 0.02, zmax = 0.17)

# Convert raster of mean current net div rates to Mollweide
Ponerinae_mean_net_div_rates_Mollweide_contrasted <- raster::projectRaster(from = Ponerinae_mean_net_div_rates_WGS84_contrasted,
                                                                           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
                                                                           method = "ngb")
plot(Ponerinae_mean_net_div_rates_Mollweide_contrasted)

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

pdf(file = paste0("./outputs/BAMM/Ponerinae_mean_net_div_rates_map.pdf"),
    width = 10, height = 5)

print(Ponerinae_mean_net_div_rates_ggplot)

dev.off()


##### 14/ Investigate relationship between PP of rate shift vs. PP of dispersal event ####

### 14.1/ Compute PP of dispersal events on branch ####

# Include both anagenetic events (d: range extension) and cladogenetic events (j: jump dispersal) where the parental area differ from the descendant area
# Get counts on each BSM
# Get frequencies across BSM to obtain Posterior probabilities

# Load phylogeny
Ponerinae_phylogeny_1534t_treedata <- readRDS(file = "./outputs/Grafting_missing_taxa/Ponerinae_phylogeny_1534t_treedata.rds")
phy <- Ponerinae_phylogeny_1534t_treedata@phylo

# Load BSM outputs
DEC_J_BSM_output <- readRDS(file = "./outputs/BSM/DEC_J_BSM_output.rds")

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

## 14.1.1/ Counts jump-dispersal events per branches ####

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

## 14.1.2/ Counts range extension events per branches ####

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
  range_extension_edges_ID_list[[i]] <- range_extension_ID_i
  
  # Print progress
  if (i %% 100 == 0)
  {
    cat(paste0(Sys.time(), " - Branches with range extension detected for BSM n°", i, "/", length(DEC_J_ana_events_tables),"\n"))
  }
}

## 14.1.3/ Aggregate counts of jump dispersal and range extension = dispersal events

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
saveRDS(object = dispersal_edges_binary_df, file = "./outputs/BAMM/dispersal_edges_binary_df.rds")

## 14.1.4/ Compute posterior probability of dispersal as relative frequencies ####

dispersal_PP_per_edges <- apply(X = dispersal_edges_binary_df, MARGIN = 2, FUN = sum)
dispersal_PP_per_edges <- dispersal_PP_per_edges / nrow(dispersal_edges_binary_df)

hist(dispersal_PP_per_edges)
# Bimodal distribution !

# Save posterior probability of dispersal per branches
saveRDS(object = dispersal_PP_per_edges, file = "./outputs/BAMM/dispersal_PP_per_edges.rds")


### 14.2/ Test for relationship between PP of diversification rate shift and PP of dispersal ####

# No valid test to account for the non-independence of edge data
# Degrees of freedom should be adjusted in the test, but not sure how...
# Only compute a stat and display a plot. Do not trust p-values 

## 14.2.1/ Convert PP to binary variables ####

## Load posterior probability of dispersal per branches
dispersal_PP_per_edges <- readRDS(file = "./outputs/BAMM/dispersal_PP_per_edges.rds")

# Distribution of dispersal PP is bimodal so convert into a binary variable for the test
hist(dispersal_PP_per_edges)

# Apply 0.5 threshold
binary_dispersal_PP_per_edges <- (dispersal_PP_per_edges >= 0.5)
table(binary_dispersal_PP_per_edges)

## Load df of edge prior probs, posterior probs, and odd-ratio of rate shifts
edge_rate_shift_probs_df <- readRDS(file = "./outputs/BAMM/edge_rate_shift_probs_df.rds")

# Distribution of posterior probabilities and odd-ratios are also highly right-skewed
# Use a threshold to identify credible shift locations
# Does not mean a shift is likely to have occurred on all these edges as odd-ratio are independent, but in reality rate shift events alongside successive branches are not
hist(edge_rate_shift_probs_df$marg_posterior_probs)
hist(edge_rate_shift_probs_df$marginal_odd_ratios)

# Apply odd-ratio > 5 as a threshold to detect credible locations for rate shifts

binary_rate_shifts_PP_per_edges <- (edge_rate_shift_probs_df$marginal_odd_ratios > 5)
table(binary_rate_shifts_PP_per_edges)

## 14.2.2/ Apply Fisher's exact test ####

binary_rate_shifts_vs_dispersal_fisher_test <- fisher.test(x = binary_rate_shifts_PP_per_edges,
                                                           y = binary_dispersal_PP_per_edges,
                                                           alternative = "two.sided")

binary_rate_shifts_vs_dispersal_fisher_test
binary_rate_shifts_vs_dispersal_fisher_test$p.value  # p = 0.216
binary_rate_shifts_vs_dispersal_fisher_test$estimate # Odds-ratio = 0.592
binary_rate_shifts_vs_dispersal_fisher_test$conf.int # CI 95% = 0.246 - 1.234

## 14.2.3/ Apply Independence Khi² test ####

binary_rate_shifts_vs_dispersal_chisq_test <- chisq.test(x = binary_rate_shifts_PP_per_edges, y = binary_dispersal_PP_per_edges)

binary_rate_shifts_vs_dispersal_chisq_test
binary_rate_shifts_vs_dispersal_chisq_test$statistic  # Khi² = 1.596
binary_rate_shifts_vs_dispersal_chisq_test$p.value    # p = 0.207
khi_Q95 <- qchisq(p = 0.95, df = binary_rate_shifts_vs_dispersal_chisq_test$parameter) ; khi_Q95  # Q95% = 3.841


?qchisq

binary_rate_shifts_vs_dispersal_chisq_test$parameter

### 14.3/ Association plot for Independence of rate shifts vs. dispersal events ####

# Compute contingency table
contingency_table <- table(binary_rate_shifts_PP_per_edges, binary_dispersal_PP_per_edges)
dimnames(contingency_table) <- list(`Rate shift events` = c("No shift", "Rate shift"),
                                    `Dispersal events` = c("No dispersal", "Dispersal"))

# Plot association plot
pdf("./outputs/BAMM/Association_plot_rate_shifts_vs_dispersal.pdf", height = 8, width = 10)

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



### Bonus in other scripts (?)

# For a more smooth model with many small cladogenetic changes, see CLaDS
# Compare Regional in situ radiations (RISR)
# GeoSSE, HiSSE, BiSSE for Old World vs. New world



