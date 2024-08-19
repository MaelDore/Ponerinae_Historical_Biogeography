##### Script: Design root age prior #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Empirically design the 95% HPD interval from the root prior based on previous time estimates #
# Based on Table 1 from Borowiec et al., 2021

### Source
 # Borowiec, M. L., Moreau, C. S., & Rabeling, C. (2021). Ants: phylogeny and classification. Encyclopedia of social insects, 52-69.


# Clean environment
rm(list = ls())

##### 1/ Load packages ####

library(BayesTwin)


##### 2/ Build df of Poneroids age estimates from previous study ####

# Get min/max range from Table 1 in Borowiec et al., 2021

Poneroids_age_estimates_df <- data.frame(source = c("Brady et al., 2006", "Moreau et al., 2006", "Schmidt et al., 2013", "Moreau & Bell, 2013", "Blanchard & Moreau, 2017", "Economo et al., 2018", "Borowiec et al., 2019"),
                                         range_min = c(100, 122, 93, 69, 105, 127, 81),
                                         range_max = c(115, 134, 128, 107, 147, 149, 102))

# Compute mean
Poneroids_age_estimates_df$mean <- (Poneroids_age_estimates_df$range_min + Poneroids_age_estimates_df$range_max) / 2
Poneroids_age_estimates_df$sd <- (Poneroids_age_estimates_df$range_max - Poneroids_age_estimates_df$mean) / 1.96

##### 3/ Obtain and aggregate empirical samples from each study #####

nb_samples_per_study <- 1000000

set.seed(seed = 1234)

age_samples <- c()
for (i in 1:nrow(Poneroids_age_estimates_df))
{
  samples_study_i <- rnorm(n = nb_samples_per_study, mean = Poneroids_age_estimates_df$mean[i], sd = Poneroids_age_estimates_df$sd[i])
  age_samples <- c(age_samples, samples_study_i)
}

##### 4/ Compute 95% HPD from empirical data #####

length(age_samples)
mean(Poneroids_age_estimates_df$mean) ; mean(age_samples)
aggregated_interval <- BayesTwin::HPD(sample = age_samples, cred_int = 0.95)
aggregated_interval

# New interval: Min = 78.6 My ; Max = 144.3 My

min(Poneroids_age_estimates_df$range_min)
max(Poneroids_age_estimates_df$range_max)

# Previous interval: Min = 69 My ; Max = 149 My

