##### Script for BAMM diversification analyses  #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Estimate diversification rates per branches and number and location of rate shift along branches
# Plot mean rates per branches
# Plot rates through time: overall and per clades
# Compare submodel fits (Posterior probabilities of nb of shifts + Bayes Factors)
# Display set of credible rate shift configurations
# Identify location of core-shifts: odd-ratios of posterior/prior proabilities
# Identify the most likely configuration: MAP and MSC
# Identify macroevolutionary cohorts of taxa
# Test for differences in rates between trait data
  # STRAPP Test: traitDependentBAMM function

###

### Inputs

# (Set of) Time-calibrated phylogeny(ies)
# Tip trait data

###

### Outputs

# BAMM models
# Posterior Probabilities of nb of rates
# Bayes Factors between pairwise submodels with different number of rate shifts
# Plot of rates of diversification along branches
# Plot rates of diversification through time: overall and per clades
# Plot the set of credible rate shift configurations
# Plot phylogeny scaled by posterior probability of rate shift to observe relationship between regime shift and any mapped trait
# Identify the most likely configuration: MAP and MSC
# Plot the macroevolutionary cohorts matrix of taxa
# Test for differences in rates between trait data
  # STRAPP Test: traitDependentBAMM function

###

### Alternative analyses

# For a more smooth model with many small cladogenetic changes, see CLaDS
# GeoSSE, HiSSE, BiSSE, QuSSE for Trait dependant models of diversification
# Canonical Phylogenetic Ordination (CPO) for an event/node focus approach of trait/diversification relationships


# Clean environment
rm(list = ls())

library(ape)
library(BAMMtools)
library(coda)


##### 1/ Load input files ####

### 1.1/ Inform path to BAMM folder ####

BAMM_path <- "./software/bamm-2.5.0/"

### 1.2/ Load phylogeny and check if it is valid

phy <- readRDS(file = "./outputs/Grafting_missing_taxa/Platythyrea_phylogeny_44t.rds")

# Check if ultrametic
is.ultrametric(phy)
# Check if fully resolved
is.binary(phy)
# Check that all branch have positive length
min(phy$edge.length) > 0 ; min(phy$edge.length)

### 1.3/ Export phylogeny as .tree file for the analyses
write.tree(phy = phy, file = paste0(BAMM_path, "my_phy.tree"))
phy_path <- paste0(BAMM_path, "my_phy.tree")

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
lambdaInitPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaInitPrior = "))[1]
lambdaInitPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaInitPrior_default_line], pattern = "lambdaInitPrior = "))
lambdaInitPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaInitPrior = "))[1]
my_div_control_file[lambdaInitPrior_line] <- paste0("lambdaInitPrior = ", lambdaInitPrior)

# Set the standard deviation of the normal distribution prior(s) of rate variation parameters (alpha) of speciation rate regimes
# alpha in lambda(t) = lamba0 x exp(alpha*t)
# Mean of this prior(s) are fixed to zero such as a constant rate diversification process is the most probable a priori
lambdaShiftPrior <- 0.05
lambdaShiftPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaShiftPrior = "))[1]
lambdaShiftPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaShiftPrior_default_line], pattern = "lambdaShiftPrior = "))
lambdaShiftPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaShiftPrior = "))[1]
my_div_control_file[lambdaShiftPrior_line] <- paste0("lambdaShiftPrior = ", lambdaShiftPrior)

# Set the rate parameter of the exponential prior(s) of initial lambda parameters (mu0) of extinction rate regimes
# mu0 in mu(t) = mu0 x exp(alpha*t)
# As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
muInitPrior <- 1.0
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
numberOfGenerations = 1000000 # 10^7
# numberOfGenerations = 1000 # 10^3
numberOfGenerations_line <- which(str_detect(string = my_div_control_file, pattern = "numberOfGenerations = "))[1]
my_div_control_file[numberOfGenerations_line] <- paste0("numberOfGenerations = ", numberOfGenerations)

# Set the path to the MCMC output file
# Includes only summary information about MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
mcmcOutfile_path <- "mcmc_log.txt"
mcmcOutfile_path_line <- which(str_detect(string = my_div_control_file, pattern = "mcmcOutfile = "))[1]
my_div_control_file[mcmcOutfile_path_line] <- paste0("mcmcOutfile = ", mcmcOutfile_path)

# Set the frequency in which to write the MCMC output to the log file
# Aim for 5000 posterior samples ideally
mcmcWriteFreq <- 1000000 / 5000 
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
eventDataWriteFreq <- 1000 # Sample every 10^3 generations
# eventDataWriteFreq <- 100 # Sample every 100 generations
eventDataWriteFreq_line <- which(str_detect(string = my_div_control_file, pattern = "eventDataWriteFreq = "))[1]
my_div_control_file[eventDataWriteFreq_line] <- paste0("eventDataWriteFreq = ", eventDataWriteFreq)

# Set frequency in which to print MCMC status to the screen
printFreq <- 1000 # Print status every 10^3 generations
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
run_prefix_for_output_files <- "BAMM"
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

# Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
# lambda0 in lambda(t) = lamba0 x exp(alpha*t)
lambdaInit0 <- 0.032
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
minCladeSizeForShift <- 1
# minCladeSizeForShift <- 3
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
writeLines(text = my_div_control_file, con = paste0(BAMM_path, "my_div_control_file.txt"))


##### 3/ Run BAMM with calls to command lines from within r using system() #####

### 3.1/ Run BAMM ####

?system

# Version
system(paste0(BAMM_path,"bamm --version"))

# Example run
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "my_div_control_file.txt"))

### Outputs
# The run_info.txt file, containing a summary of your parameters/settings
# An mcmc_log.txt containing raw MCMC information useful in diagnosing convergence
# An event_data.txt file containing all of evolutionary rate parameters and their topological mappings
# A chain_swap.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).

### 3.2/ Clean outputs (move them to a dedicated folder) ####

BAMM_output_folder_path <- "./outputs/BAMM/BAMM_test_run/"

# Detect output files
output_files_path <- list.files(path = "./", pattern = "BAMM_")
# Move output files to dedicated folder
file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "my_phy.tree")) # Phylo file
file.rename(from = paste0(BAMM_path, "my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "my_div_control_file.txt")) # Control file


##### 4/ Check chain convergence #####

# Load the MCMC log file
MCMC_log <- read.csv(paste0(BAMM_output_folder_path, run_prefix_for_output_files, "_mcmc_log.txt"), header = T)

## Plot trace of the sampled generations (mixed of chains using swapping)
plot(MCMC_log$logLik ~ MCMC_log$generation, type = "o")

## Remove burn-in
burn_in <- 0.1
burn_in <- 0.25

burn_in_threshold <- ceiling(burn_in * MCMC_log$generation[nrow(MCMC_log)])
post_burn_MCMC <- MCMC_log[MCMC_log$generation >= burn_in_threshold, ]

## Explore effective sample sizes

coda::effectiveSize(post_burn_MCMC$N_shifts) # Nb of shifts in selected submodel
coda::effectiveSize(post_burn_MCMC$eventRate) # LAMBDA parameter for the number of shifts
coda::effectiveSize(post_burn_MCMC$logLik) # Overall likelihood


##### 5/ Load and format output files ####

### 5.1/ Load example data files ####

data(primates, events.primates)
data(whales, events.whales)

phy <- read.tree(paste0(BAMM_output_folder_path, "my_phy.tree"))

### 5.2/ Format BAMM output into R object ####

?BAMMtools::getEventData

# Create the bammdata summarizing BAMM outputs
BAMM_output_test <- getEventData(phy = phy, 
                                 eventdata = paste0(BAMM_output_folder_path, run_prefix_for_output_files, "_", eventDataOutfile_path),
                                 burnin = burn_in,
                                 type = "diversification")


# Create the bammdata summarizing BAMM outputs
BAMM_output_primates <- getEventData(phy = primates, eventdata = events.primates,
                                     burnin = 0.25,
                                     type = 'trait') # Here trait because the analysis also involved modeling of quantative trait evolution

# Create the bammdata summarizing BAMM outputs
BAMM_output_whales <- getEventData(phy = whales, eventdata = events.whales,
                                   burnin = 0.1,
                                   type = "diversification")

# Basically an extended ape phylogeny object storing BAMM output data
str(BAMM_output_primates, max.level = 1)
BAMM_output_primates$begin # Absolute time since root of edge/branch start
BAMM_output_primates$end # Absolute time since root of edge/branch end
BAMM_output_primates$numberEvents # Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts
BAMM_output_primates$eventData # List of dataframes recording shift events in each posterior configuration
BAMM_output_primates$eventData[[1]] # Dataframe recording shift events and macroevolutionary regimes in the focal posterior configuration. 1st line = Background root regime
BAMM_output_primates$eventData[[1]]$node # Tipward node ID of the branch where the shift occured
BAMM_output_primates$eventData[[1]]$time # Absolute time since root when the shift occurred on the branch
BAMM_output_primates$eventData[[1]]$lam1 # Initial rate of speciation (lambda0) of the regime
BAMM_output_primates$eventData[[1]]$lam2 # Speciation rate change parameter (alpha/z) of the regime. If = 0, rate is constant. If < 0 = decay. If > 0 = growth
BAMM_output_primates$eventData[[1]]$mu1 # Initial rate of extinction (mu0) of the regime. Should be the constant extinction rate of the regime.
BAMM_output_primates$eventData[[1]]$mu2 # Extinction rate change parameter (alpha/z) of the regime. Should be fixed to 0 as extinction rates are constant in BAMM.
BAMM_output_primates$eventVectors # List of integer vectors of regime membership per branches in each posterior configuration
BAMM_output_primates$eventVectors[[1]] # Integer vectors of regime membership per branches
BAMM_output_primates$eventBranchSegs[[1]] # Same but with matrix including tipward node ID and begin/end ages of the branches
BAMM_output_primates$tipStates[[1]] # Integer vectors of regime membership per tips
BAMM_output_primates$tipLambda[[1]] # Integer vectors of final speciation rates at tips = current speciation rates
BAMM_output_primates$meanTipLambda # Mean current tip speciation rates across all posterior configurations. Better to use the median and 95% HPD
BAMM_output_primates$tipMu[[1]] # Integer vectors of final extinction rates at tips = current extinction rates
BAMM_output_primates$tipLambda[[1]] - BAMM_output_primates$tipMu[[1]] # Integer vectors of final diversification  rates at tips = current diversification rates
BAMM_output_primates$meanTipMu # Mean current tip extinction rates across all posterior configurations. Better to use the median and 95% HPD


##### 6/ Explore rate heterogeneity = number of shifts #####

### 6.1/ Using Posterior Probabilities ####

## Quick, but sensitive to priors

## From the MCMC log

# Get posterior samples of number of shifts in sampled submodels
PP_N_shifts <- table(post_burn_MCMC$N_shifts) / nrow(post_burn_MCMC)
PP_N_shifts

# Summary stats
summary(post_burn_MCMC$N_shifts)

# Compute posterior odds ratio between different submodels
# M0 vs. submodels with rate shifts

# If no submodel with no-shift was sampled. Artificially add a half-weight as the conservative hypothesis
if (is.na(PP_N_shifts["0"])) 
{
  PP_N_shifts <- c(0.5, PP_N_shifts)
  names(PP_N_shifts)[1] <- "0"
}

PP_N_shifts["0"] / sum(PP_N_shifts[names(PP_N_shifts) != "0"])
# If odds-ratio > 5 => substantial evidence / "significance"

## From the Event data log

# Get posterior samples of number of shifts in sampled submodels
PP_N_shifts <- table(BAMM_output_test$numberEvents) / length(BAMM_output_test$numberEvents)
PP_N_shifts

# Summary stats
summary(BAMM_output_test$numberEvents)

# Compute posterior odds ratio between different submodels
# M0 vs. submodels with rate shifts

PP_N_shifts["0"] / PP_N_shifts[names(PP_N_shifts) != "0"]
# If odds-ratio > 5 => substantial evidence / "significance"

### 6.2/ Using Bayes Factors to compare models  ####

?computeBayesFactors

Bayes_factors_pairwise <- computeBayesFactors(postdata = post_burn_MCMC,
                                              expectedNumberOfShifts = expectedNumberOfShifts,
                                              burnin = burn_in)
# Matrix of Bayes factors for pairwise comparison of submodels
Bayes_factors_pairwise

# Bayes factors comparing all submodels vs. null model of no-rate shifts
Bayes_factors_pairwise[1, ]

## Bayes Factors interpretation scale
# BF > 12 = weak evidence (p < 0.05)
# BF > 20 = strong evidence (p < 0.01)
# BF > 50 = Very strong evidence  (p < 0.001)


### 6.3/ Compare prior and posterior distribution of LAMBDA (parameter controlling the nb of rates) ####

# Prior distribution is a Poisson, so we use a bar graph for a discrete distribution

data(mcmc.whales)

plotPrior(mcmc = mcmc.whales, expectedNumberOfShifts=1)
plotPrior(mcmc = post_burn_MCMC,
          expectedNumberOfShifts = expectedNumberOfShifts,
          burnin = burn_in)

##### 7/ Plot branch diversification rates on the phylogeny ####

?plot.bammdata
?BAMMtools::dtRates # To compute discretized rates on small segments along branches

### 7.1/ Extract mean branch rates ####

# Get mean PP branch speciation/extinction rates
all_branch_div_rates <- getCladeRates(BAMM_output_whales)

# Compute overall rate as the mean of mean branch rates
mean(all_branch_div_rates$lambda - all_branch_div_rates$mu)
median(all_branch_div_rates$lambda - all_branch_div_rates$mu)
quantile(all_branch_div_rates$lambda - all_branch_div_rates$mu, c(0.05, 0.95)) # 95% CI interval of the mean branch rate across the whole phylogeny

### 7.2/ Scale branch length to diversification rates ####

div_rates_tree <- getMeanBranchLengthTree(BAMM_output_whales)
plot(div_rates_tree$phy)

# Check if it is using diversification or speciation rates...

### 7.3/ Plot mean sliding rates + credible shift locations ####

## Good practice wood be to use median rather than mean
# Median can be computed from BAMM_output_primates$eventData or using BAMMtools::dtRates to compute discretize rates

# With different color palette settings

par(mfrow=c(1,3), mar=c(1, 0.5, 0.5, 0.5), xpd=TRUE)

# tau = length of segments used to discretize the continuous rates
# tau = 0.001 => each segment is 1/1000 of the tree height (root age)

# Fully linear
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'linear',
                                       lwd = 2)
addBAMMshifts(BAMM_output_whales, par.reset = FALSE, cex = 2)
title(sub = 'Linear', cex.sub = 2, line = -1)
addBAMMlegend(div_rates_plot_whales, location=c(0, 1, 140, 220))

# Linear with a max bound
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'linear', color.interval = c(NA,0.12),
                                       lwd = 2)
addBAMMshifts(BAMM_output_whales, par.reset = FALSE, cex = 2)
title(sub = 'Linear with max interval', cex.sub = 2, line = -1)
addBAMMlegend(div_rates_plot_whales, location=c(0, 1, 140, 220))

# Jenks method to discretize rates in categories min/max variance within/between groups
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'jenks',
                                       lwd = 2)
addBAMMshifts(BAMM_output_whales, par.reset = FALSE, cex = 2)

title(sub = 'Jenks method', cex.sub = 2, line = -1)
addBAMMlegend(div_rates_plot_whales, location = c(0, 1, 140, 220))


## Get the histogram used to discretize color bins

par(mfrow=c(3,1), mar=c(1, 0.5, 0.5, 0.5), xpd=TRUE)

# Fully linear
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'linear',
                                       show = F)
ratesHistogram(div_rates_plot_whales, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Linear', cex.main=1)

# Linear with a max bound
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'linear',
                                       color.interval = c(NA, 0.12),
                                       show = F)
ratesHistogram(div_rates_plot_whales, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Linear with max interval', cex.main=1)

# Jenks method to discretize rates in categories min/max variance within/between groups
div_rates_plot_whales <- plot.bammdata(BAMM_output_whales,
                                       tau = 0.001, breaksmethod = 'jenks',
                                       show = F)
ratesHistogram(div_rates_plot_whales, plotBrks = TRUE, xlab = 'Diversification rates')
title(main = 'Jenks method', cex.main=1)



### 7.4/ Plot uncertainty in diversification rates ####

# Need to define the length of segments used to discretize time
# Compute continuous rates on branches from the regime parameters
# Extract summary stats on rates per branches for each time slice (mean, median, CI, HPD)
# Map those summary stats on phylo

?BAMMtools::dtRates # To compute discretized rates on small segments along branches

# Median can be computed from BAMM_output_primates$eventData or using BAMMtools::dtRates to compute discretize rates

# To get a subset of posterior samples to fasten computation
BAMM_output_whales_subset <- subsetEventData(BAMM_output_whales, index = index)


##### 8/ Get clade-specific diversification rates #####

?getCladeRates

# Get mean PP speciation/extinction rates for branch within a clade
dolphins_branch_div_rates <- getCladeRates(BAMM_output_whales, node = 141)

mean(dolphins_branch_div_rates$lambda)
median(dolphins_branch_div_rates$lambda)
quantile(dolphins_branch_div_rates$lambda, c(0.05, 0.95))

# Get background PP speciation/extinction rates for all branches, except the focal clade
non_dolphins_branch_div_rates <- getCladeRates(BAMM_output_whales, node = 141, nodetype = "exclude")
mean(non_dolphins_branch_div_rates$lambda)
quantile(non_dolphins_branch_div_rates$lambda, c(0.05, 0.95))


##### 9/ Plot diversification rates through time #####

?branching.times # To extract divergence/branching times
?plotRateThroughTime

### 9.1/ Plot of rates through time with fuzzy intervals: mean + all trajectories ####

# Can be partitioned per clades (and per bioregions?)

plot.new()

par(mfrow = c(1, 3))

root_age <- max(branching.times(whales))

# A/ All branches
plotRateThroughTime(BAMM_output_whales,
                    intervalCol = "red", avgCol = "red",
                    start.time = root_age,
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 2)
text(x = 30, y = 0.8, label = "All whales", font = 4, cex = 2.0, pos = 4)

# B/ Only dolphins
plotRateThroughTime(BAMM_output_whales,
                    intervalCol = "blue", avgCol = "blue", 
                    start.time = root_age, 
                    node = 140,
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 30, y = 0.8, label = "Dolphins only", font = 4, cex = 2.0, pos = 4)

# C/ Only non-dolphins (background regime)
plotRateThroughTime(BAMM_output_whales, 
                    intervalCol = "darkgreen", avgCol = "darkgreen",
                    start.time = root_age, 
                    node = 140, nodetype = "exclude",
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5)
text(x = 30, y = 0.8, label = "Non-dolphins", font = 4, cex = 2.0, pos = 4)


### 9.2/ Plot of rates through time with strict intervals: mean + 95% CI polygon ####

# Can be partitioned per clades (and per bioregions?)

plot.new()

par(mfrow = c(1, 3))

root_age <- max(branching.times(whales))

# A/ All branches
plotRateThroughTime(BAMM_output_whales,
                    avgCol = "red",
                    start.time = root_age, 
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 2,
                    intervalCol = "red",
                    intervals = c(0.05, 0.95),
                    opacity = 0.2)
text(x = 30, y = 0.8, label = "All whales", font = 4, cex = 2.0, pos = 4)

# B/ Only dolphins
plotRateThroughTime(BAMM_output_whales,
                    avgCol = "blue",
                    start.time = root_age,
                    node = 140,
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 2,
                    intervalCol = "blue",
                    intervals = c(0.05, 0.95),
                    opacity = 0.2)
text(x = 30, y = 0.8, label = "Dolphins only", font = 4, cex = 2.0, pos = 4)

# C/ Only non-dolphins (background regime)
plotRateThroughTime(BAMM_output_whales, 
                    avgCol = "darkgreen",
                    start.time = root_age, 
                    node = 140, nodetype = "exclude",
                    ylim = c(0, 1),
                    xlim = c(root_age, 0),
                    cex.axis = 1.5,
                    intervalCol = "darkgreen", 
                    intervals = c(0.05, 0.95),
                    opacity = 0.2)
text(x = 30, y = 0.8, label = "Non-dolphins", font = 4, cex = 2.0, pos = 4)


### 9.3/ Extract rates through time data ####

## Can set a time window and the size of the sliding window
# Use start.time, end.time, and nslices

# Extract for the whole tree
rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_output_whales)

str(rates_through_time_matrix, max.level = 1)

dim(rates_through_time_matrix$lambda) # Speciation rates per posterior sample (raw) x time window (columns)
dim(rates_through_time_matrix$mu) # Extinction rates per posterior sample (raw) x time window (columns)
rates_through_time_matrix$times # Mean time of the time windows. Names = age. Values = time since root.

# Extract for a subclade
dolphins_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_output_whales, node = 140)

# Extract for background regime (excluding a subclade)
non_dolphins_rates_through_time_matrix <- getRateThroughTimeMatrix(BAMM_output_whales, node = 140, nodetype = "exclude")


##### 10/ Plot regime shifts #####

# Marginal probabilities. If high for a sequence of branches, good practice could be to look for in the posterior samples for the joint probability of a shift happening on 0/1/2/… of these branches in a single sample.
# Goal: avoid misinterpreting those independent marginal probabilities has evidence for a sequence of shifts as what it actually supports is the likely presence of one shift among those sequence of branches, whit a bit of uncertainty on the exact location among those successive branches

# Caution: showing all marginal PP of shift location may be biases since shift location may not be independent
# The presence of one shift may be associated with the absence of another. Be careful with the interpretation of marginal PP.
# Try to investigate joint probabilities of some events.

### 10.1/ Plot the set of credible shift configurations ####

# To extract the set of credible shift configuration found in the posterior samples
# Extract configurations until cum prob reach 95%
# Probability of shift as proportion in the posterior samples
# Use a threshold to select only shifts with a minimal marginal odd ratio that depict substantially higher PP than prior probability = core-shifts.
# If not ignored, many configurations will contain non-core shifts that are not more expected than according to the prior, thus meaningless to model our data

?credibleShiftSet
?plot.credibleshiftset

cset <- credibleShiftSet(BAMM_output_whales,
                         expectedNumberOfShifts = 1, # Initial parameter used in the config file to set the exponential prior of the lambda parameter for number of shifts
                         threshold = 3) # Significance threshold of odd-ratios between prior and posterior probabilities of rate shift on branches
# Display summary with PP of each configuration and associated number of core shifts (with odd-ratio > threshold)
summary(cset)

# Plot the set of credible shift configurations
plot.credibleshiftset(cset, lwd = 2.5)

# All PP of the set of configurations sum > 95%


### 10.2/ Plot random samples of the best configurations ####

# Provide rank (as from the PP) of the selected configuration
# Differences are in rate parameter, but shift locations are the same for a given rank

?distinctShiftConfigurations

dsc <- distinctShiftConfigurations(BAMM_output_whales,
                                   expectedNumberOfShifts = 1,
                                   threshold = 5)

# Here is one random sample with the BEST shift configuration (rank = 1)
plot.bammshifts(x = dsc, ephy = BAMM_output_whales,
                rank = 1, legend = F)
# Here is another (read the label text):
plot.bammshifts(x = dsc, ephy = BAMM_output_whales,
                rank = 1, legend = F)

# Here is one random sample with the second best shift configuration (rank = 2)
plot.bammshifts(x = dsc, ephy = BAMM_output_whales,
                rank = 2, legend = F)
# Here is another
plot.bammshifts(x = dsc, ephy = BAMM_output_whales,
                rank = 2, legend = F)


### 10.3/ Extract node ID of shift locations for a given sample ####

selected_sample <- sample(x = dsc$samplesets[[1]], size = 1) 

shift_nodes <- getShiftNodesFromIndex(BAMM_output_whales, index = selected_sample)

plot.phylo(whales)
nodelabels(node = shift_nodes, pch = 21, col = "black", bg = "red", cex = 1.5)


### 10.4/ Scale phylogeny according to the Marginal Probability of shift ####

?marginalShiftProbsTree

marg_probs <- marginalShiftProbsTree(BAMM_output_whales)
marg_probs$edge.length
plot.phylo(marg_probs)

# Highlight credible location of shifts
# Not independent! Successive high values do not mean that shifts are likely on all branches
# Rather they indicate uncertainty in the location of one single shift


### 10.5/ Compute marginal posterior odd-ratios ####

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
branch_priors <- getBranchShiftPriors(BAMM_output_whales, expectedNumberOfShifts = 1)
branch_priors$edge.length
# Should be directly proportional to branch length

# Branch length are scaled to prior probabilities of shift per branches
sum(BAMM_output_whales$edge.length)
sum(branch_priors$edge.length)

# Compute marginal odds ratio to account for prior probabilities
# MO =  posterior proba of rate shift / prior proba of rate shift. 
# The highest, the more substantial is the support for a shift occurring on that branch.
mo <- marginalOddsRatioBranches(ephy = BAMM_output_whales, expectedNumberOfShifts = 1)

# Branch length are scaled to odd-ratios of rate shift per branch
mo$edge.length
sum(mo$edge.length)

plot(marg_probs) # Posterior probabilities
plot(branch_priors) # Prior probabilities
plot(mo) # Odd-ratios = Posterior probabilities / Prior probabilities


##### 11/ Identify the most likely configuration ####

### 11.1/ Maximum A posteriori Probability (MAP) configuration ###

## Get the most frequent configuration = maximum a posteriori probability (MAP) shift configuration = the single configuration of shift location showing up the most in the posterior sample

# Idea: find the posterior configuration that is the most frequent in the posterior samples

MAP <- getBestShiftConfiguration(BAMM_output_whales,
                                 expectedNumberOfShifts = 1,
                                 threshold = 5) # Odd-ratio threshold used to select credible shifts
plot.bammdata(MAP, lwd = 1.25)
addBAMMshifts(MAP, cex=2)

### 11.2/ Maximum Shift Credibility configuration (MSC) ####

## Get the most likely configuration = maximum shift credibility configuration (MSC)

# Idea: find the posterior configuration that has the highest probability according to the marginal probability of the branches and number of shifts

# Useful for big phylogenies when there are no clear MAP because of the many possible configurations

# Identify the MSC configs
msc_detection <- maximumShiftCredibility(BAMM_output_whales)
# Extract the MSC config
msc_tree <- subsetEventData(BAMM_output_whales, index = msc_detection$sampleindex)
plot.bammdata(msc_tree, lwd = 1.25)
addBAMMshifts(msc_tree, par.reset = FALSE, cex = 2)


##### 12/ Identify clades with different regimes from the background #####

## Get the cumulative shift probability of a branch 
# Quantifies the proportion of posterior configurations where the focal branch belong to a different regime than the background regime at the root
# = Quantifies the proportion of posterior configurations where the branch i witness at least one shift transition on the path from the root to the branch i

# Get a tree with cumulative shift probability labeled on each branch 
cum_shift_prob <- cumulativeShiftProbsTree(ephy = BAMM_output_whales)
cum_shift_prob$edge.length

plot(cum_shift_prob)

# Highlight edge with cumulative shift probability > 0.95
is_cumprobshift_signif <- cum_shift_prob$edge.length >= 0.95
edgecols <- rep('black', length(cum_shift_prob$edge.length))
edgecols[is_cumprobshift_signif] <- "red"
plot.phylo(cum_shift_prob, edge.color = edgecols)

##### 13/ Identify pairs of tips sharing their macroevolutionary regime = Macroevolutionary cohort matrix ####

# Probability that two tips/taxa share the same macroevolutoinary regime

cmat <- getCohortMatrix(BAMM_output_whales)
cohorts(cmat, BAMM_output_whales, lwd = 3,
        pal = "temperature", use.plot.bammdata = TRUE)


##### 14/ Test for differences in rates associated with trait values: STRAPP test #####

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

# Load biogeographic data for Platythyrea
Platythyrea_bioregions_binary_table <- readRDS(file = "./input_data/Biogeographic_data/Platythyrea_bioregions_binary_table.rds")
# Transform arbitrarily in categorical monomorphic data
Platythyrea_bioregions_tip_data <- names(Platythyrea_bioregions_binary_table[, -1])[apply(X = Platythyrea_bioregions_binary_table[, -1], MARGIN = 1, FUN = which.max)]
Platythyrea_bioregions_tip_data <- as.factor(Platythyrea_bioregions_tip_data)

# Assign tip names to tip data
names(Platythyrea_bioregions_tip_data) <- phy$tip.label

STRAPP_test_output <- traitDependentBAMM(ephy = BAMM_output_test, traits = Platythyrea_bioregions_tip_data,
                                         reps = 1000, rate = "net diversification",
                                         return.full = T,
                                         method = "kruskal", # For categorical multinomial data (G > 2)
                                         logrates = T, 
                                         two.tailed = T)

STRAPP_test_output$estimate # Mean tip rates per categories
STRAPP_test_output$p.value
STRAPP_test_output$rate # Type of rates
STRAPP_test_output$obs.corr # Observed statistic for each posterior sample
STRAPP_test_output$gen # Generation ID of the selected posterior sample
STRAPP_test_output$null # Null statistic for each posterior sample



