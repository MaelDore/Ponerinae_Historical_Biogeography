
##### Functions to compute residence times per window time #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compute residence time of each state/area across time-windows

###

### Inputs

# Simmaps of state evolution from Biogeographic Stochastic Maps

###

### Outputs

# Array of residence times x bioregions x time-windows x BSM maps 

###


##### 1/ Function to compute residence times per window time from a unique Simmap ####

compute_residence_times_per_sliding_windows_on_Simmap <- function (simmap,
                                                                   time_boundaries, # Data.frame of start/end time of sliding windows
                                                                   remove_empty_states = F, # Should states with no residence times be removed?
                                                                   melted = F, # Should the final array be melted into a df?
                                                                   verbose = T)  
{
  
  # Extract root age
  root_age <- max(phytools::nodeHeights(simmap))
  
  # Extract node ages per egdes
  nodes_age_df <- as.data.frame(phytools::nodeHeights(simmap))
  nodes_age_df <- round(root_age - nodes_age_df, 5)
  names(nodes_age_df) <- c("rootward_age", "tipward_age")
  nodes_age_df$branch_length <- nodes_age_df$rootward_age - nodes_age_df$tipward_age
  nodes_age_df$node_ID <- 1:nrow(nodes_age_df)
  
  # Extract ranges/states
  ranges_list <- colnames(simmap$mapped.edge)
  
  # Paste time_boundaries
  time_boundaries_labels <- paste0(time_boundaries$rootward_times, "_", time_boundaries$tipward_times)
  
  # Initiate final array
  residence_times_per_sliding_windows_array <- array(data = NA,
                                                     dim = c(length(ranges_list), length(time_boundaries_labels)),
                                                     dimnames = list(ranges_list, time_boundaries_labels))
  
  
  ### Loop per sliding window
  for (i in 1:nrow(time_boundaries))
  {
    # i <- 91
    
    # Extract window time boundaries
    rootward_time_i <- time_boundaries$rootward_times[i]
    tipward_time_i <- time_boundaries$tipward_times[i]
    
    # cat(paste0(Sys.time(), " - Extract tree section\n"))
    
    ## Old solution with master.table from BioGeoBEARS
    # cat(paste0(Sys.time(), " - Start sectionning the tree\n"))
    # tree_section_outputs <- section_the_tree_for_time_window(tree = simmap, tipward_time = 80, rootward_time = 90)
    # 
    # tree_section_outputs$master_table
    # tree_section_outputs$tree_sections_list
    
    ### Detect edge status
    
    nodes_age_df_i <- nodes_age_df
    
    # Detect edges within time window
    nodes_age_df_i$edge_within <- (nodes_age_df_i$rootward_age <= rootward_time_i) & (nodes_age_df_i$tipward_age >= tipward_time_i)
    # Detect edges outside time window
    nodes_age_df_i$edge_outside <- (nodes_age_df_i$tipward_age >= rootward_time_i) | (nodes_age_df_i$rootward_age <= tipward_time_i)
    # Detect edges with time to remove rootward
    nodes_age_df_i$edge_on_rootward_time <- (nodes_age_df_i$rootward_age > rootward_time_i) & (nodes_age_df_i$tipward_age < rootward_time_i)
    nodes_age_df_i$rootward_diff <- sapply(X = nodes_age_df_i$rootward_age, FUN = function (x) { max(0, x - rootward_time_i) } )
    # Detect edges with time to remove tipward
    nodes_age_df_i$edge_on_tipward_time <- (nodes_age_df_i$rootward_age > tipward_time_i) & (nodes_age_df_i$tipward_age < tipward_time_i)
    nodes_age_df_i$tipward_diff <- sapply(X = nodes_age_df_i$tipward_age, FUN = function (x) { max(0, tipward_time_i - x) } )
    # Detect branches to keep
    nodes_age_df_i$to_keep <- nodes_age_df_i$edge_within | nodes_age_df_i$edge_on_rootward_time | nodes_age_df_i$edge_on_tipward_time
    
    # table(nodes_age_df_i$to_keep)
    
    ### Adjust mapping to remove states outside of time window
    
    state_maps_i <- simmap$maps
    
    # Keep only branch from the focal time window in the state maps
    state_maps_i <- state_maps_i[nodes_age_df_i$to_keep]
    names(state_maps_i) <- which(nodes_age_df_i$to_keep)
    state_maps_i
    
    ## Remove time in edges overlapping rootward time boundary
    
    if (sum(nodes_age_df_i$edge_on_rootward_time) > 0)
    {
      # Compute cumulative time from root
      states_maps_cumtime_from_root_i <- lapply(X = state_maps_i, FUN = cumsum)
      
      # Remove time before 
      for (j in which(nodes_age_df_i$edge_on_rootward_time))
      {
        # j <- 1
        
        branch_ID <- as.character(j)
        
        # How much time to remove
        rootward_diff_j <- nodes_age_df_i$rootward_diff[j]
        
        # Extract nb of states
        nb_states <- length(state_maps_i[[branch_ID]])
        
        # Find where to cut
        time_cut_location <- which.max(states_maps_cumtime_from_root_i[[branch_ID]] - rootward_diff_j > 0)
        # Remove residence time before cut (set to zero)
        state_maps_i[[branch_ID]][1:(time_cut_location-1)] <- 0
        # Use cum time since root for the cut
        state_maps_i[[branch_ID]][time_cut_location] <- states_maps_cumtime_from_root_i[[branch_ID]][time_cut_location] - rootward_diff_j
        # Keep all after the cut
        state_maps_i[[branch_ID]]
      }
      state_maps_i
    }
    
    # Reverse state order (from tip to root)
    state_maps_i <- lapply(X = state_maps_i, FUN = rev) 
    
    ## Remove time in next/younger strata
    
    if(sum(nodes_age_df_i$edge_on_tipward_time) > 0)
    {
      # Compute cumulative time from tips
      states_maps_cumtime_from_tips_i <- lapply(X = state_maps_i, FUN = cumsum)
      
      # Remove time after 
      for (j in which(nodes_age_df_i$edge_on_tipward_time))
      {
        # j <- 5
        
        branch_ID <- as.character(j)
        
        # How much time to remove
        tipward_diff_j <- nodes_age_df_i$tipward_diff[j]
        
        # Extract nb of states
        nb_states <- length(state_maps_i[[branch_ID]])
        
        # Find where to cut
        time_cut_location <- which.max(states_maps_cumtime_from_tips_i[[branch_ID]] - tipward_diff_j > 0)
        # Remove residence time after cut (set to zero)
        state_maps_i[[branch_ID]][1:(time_cut_location-1)] <- 0
        # Use cum time since root for the cut
        state_maps_i[[branch_ID]][time_cut_location] <- states_maps_cumtime_from_tips_i[[branch_ID]][time_cut_location] - tipward_diff_j
        # Keep all before the cut
        state_maps_i[[branch_ID]]
      }
    }
    
    ## Reverse again state order (from root to tip)
    state_maps_i <- lapply(X = state_maps_i, FUN = rev) 
    
    # print(state_maps_i)
    
    ##  From maps to residence times per states
    
    # cat(paste0(Sys.time(), " - Compute residence times within the time window\n"))
    
    # Compute residence times per states
    edge_states_matrix <- plyr::ldply(.data = state_maps_i, .fun = base::rbind)
    edge_states_matrix <- edge_states_matrix[,-1, drop = F] # Remove node ID
    state_times <- apply(X = edge_states_matrix, MARGIN = 2, FUN = sum, na.rm = T)
    state_times_reordered <- setNames(object = rep(0, length(ranges_list)), nm = ranges_list)
    state_times_reordered[names(state_times)] <- state_times
    
    ## Store results in final array
    residence_times_per_sliding_windows_array[,i] <- state_times_reordered
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Residence times computed for ",rootward_time_i,"-",tipward_time_i," My - Window n°", i, "/", nrow(time_boundaries),"\n"))
    }
  }
  
  ### Remove states with no residence times if needed
  if (remove_empty_states)
  {
    non_empty_states <- rowSums(residence_times_per_sliding_windows_array) != 0
    residence_times_per_sliding_windows_array <- residence_times_per_sliding_windows_array[non_empty_states, ]
  }
  
  output <- residence_times_per_sliding_windows_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(residence_times_per_sliding_windows_array)
    names(output) <- c("state", "window", "residence_time")
  }
  
  # Export result
  return(output)
  
  # cat(paste0(Sys.time(), " - Done\n"))
}


##### 2/ Function to compute residence times per window time from a list of Simmaps ####

compute_residence_times_per_sliding_windows_on_MultiSimmap <- function (multiSimmap,
                                                                        time_boundaries,  # Data.frame of start/end time of sliding windows)
                                                                        remove_empty_states = T, # Should states with no residence times be removed?
                                                                        melted = F, # Should the final array be melted into a df?
                                                                        verbose = T)  
{
  # Initiate final array
  
  # Extract ranges/states
  ranges_list <- unique(unlist(lapply(X = multiSimmap, FUN = function (x) { colnames(x$mapped.edge) })))
  
  # Paste time_boundaries
  time_boundaries_labels <- paste0(time_boundaries$rootward_times, "_", time_boundaries$tipward_times)
  
  # Paste maps
  nb_maps <- length(multiSimmap)
  map_labels <- paste0("Map_", 1:nb_maps)
  
  # Initiate final array
  residence_times_per_sliding_windows_all_maps_array <- array(data = NA,
                                                              dim = c(length(ranges_list), length(time_boundaries_labels), nb_maps),
                                                              dimnames = list(ranges_list, time_boundaries_labels, map_labels))
  
  # Loop per maps
  for (i in seq_along(multiSimmap))
  {
    # i <- 1
    
    # Compute residence times for window i
    array_i <- compute_residence_times_per_sliding_windows_on_Simmap(simmap = multiSimmap[[i]],
                                                                     time_boundaries = time_boundaries, # Data.frame of start/end time of sliding windows
                                                                     remove_empty_states = F, # Should states with no residence times be removed?
                                                                     melted = F, # Should the final array be melted into a df?
                                                                     verbose = F)  
    # Store residence times
    ranges_i <- dimnames(array_i)[[1]]
    residence_times_per_sliding_windows_all_maps_array[ranges_i, ,i] <- array_i
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Residence times computed for Map n°", i, "/", length(multiSimmap),"\n"))
    }
  }
  
  ### Remove states with no residence times if needed
  if (remove_empty_states)
  {
    non_empty_states <- rowSums(residence_times_per_sliding_windows_all_maps_array) != 0
    residence_times_per_sliding_windows_all_maps_array <- residence_times_per_sliding_windows_all_maps_array[non_empty_states, ,]
  }
  
  output <- residence_times_per_sliding_windows_all_maps_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(residence_times_per_sliding_windows_all_maps_array)
    names(output) <- c("state", "window", "map", "residence_time")
  }
  
  # Export result
  return(output)
}


##### 3/ Function to aggregated residence times per unique areas ####

# Aggregate residence times per unique areas instead of ranges
# Time for ranges encompassing multiple areas is equally split across those areas

aggregate_residence_times_per_sliding_windows_in_unique_areas <- function (residence_times_per_sliding_windows_in_ranges_array,
                                                                           remove_empty_areas = T, # Should areas with no residence times be removed?
                                                                           melted = F) # Should the final array be melted into a df?
{
  # Extract ranges
  all_ranges <- dimnames(residence_times_per_sliding_windows_in_ranges_array)[[1]]
  
  # Extract unique areas
  all_areas <- all_ranges[nchar(all_ranges) == 1]
  all_areas <- all_areas[!(all_areas == "_")] # Remove empty area
  
  # Split ranges in unique areas
  ranges_split_list <- str_split(string = all_ranges, pattern = "")
  
  ## Divide residence times of ranges equally across areas
  
  residence_times_per_sliding_windows_in_ranges_split_array <- residence_times_per_sliding_windows_in_ranges_array
  for (i in seq_along(all_ranges))
  {
    # i <- 10
    
    # Extract focal range
    range_i <- all_ranges[i]
    
    # Extract number of areas composing the range
    nb_areas_i <- nchar(range_i)
    
    # Divide residence times accordingly
    residence_times_per_sliding_windows_in_ranges_split_array[range_i,,] <- residence_times_per_sliding_windows_in_ranges_array[range_i,,] / nb_areas_i
    
  }
  
  # residence_times_per_sliding_windows_in_ranges_array[,1,"Map_1"]
  # residence_times_per_sliding_windows_in_ranges_split_array[,1,"Map_1"]
  
  ## Loop per areas to aggregate data
  residence_times_per_sliding_windows_in_areas_array <- residence_times_per_sliding_windows_in_ranges_array[all_areas,,]
  residence_times_per_sliding_windows_in_areas_array[,,] <- NA
  for (i in seq_along(all_areas))
  {
    # i <- 1
    
    # Extract focal area
    areas_i <- all_areas[i]
    
    # Find matching ranges
    matching_ranges_i <- suppressWarnings(all_ranges[str_detect(string = ranges_split_list, pattern = areas_i)])
    
    # Sum counts across all matching ranges
    residence_times_per_sliding_windows_in_areas_array[areas_i,,] <- apply(X = residence_times_per_sliding_windows_in_ranges_split_array[matching_ranges_i,,], MARGIN = c(2,3), FUN = sum)
  }
  
  # residence_times_per_sliding_windows_in_ranges_array[,1,"Map_1"]
  # residence_times_per_sliding_windows_in_ranges_split_array[,1,"Map_1"]
  # residence_times_per_sliding_windows_in_areas_array[,1,"Map_1"]
  
  ### Remove areas with no residence times if needed
  if (remove_empty_areas)
  {
    non_empty_areas <- rowSums(residence_times_per_sliding_windows_in_areas_array) != 0
    residence_times_per_sliding_windows_in_areas_array <- residence_times_per_sliding_windows_in_areas_array[non_empty_areas, ,]
  }
  
  output <- residence_times_per_sliding_windows_in_areas_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(residence_times_per_sliding_windows_in_areas_array)
    names(output) <- c("area", "window", "map", "residence_time")
  }
  
  # Export result
  return(output)
  
}
