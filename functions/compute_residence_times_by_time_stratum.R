# Compute residence times per time stratum

compute_residence_times_by_time_stratum <-  function (simmap, model_fit, stratum_ID)
{
  # Use the master table to filter out times from the state_maps that do not belong to a time stratum
  
  state_maps <- simmap$maps
  master_table <- model_fit$inputs$master_table
  
  # Get ID of branches found in Stratum i
  branches_ID <- master_table$parent_br[master_table$stratum == stratum_ID]
  branches_ID <- branches_ID[!is.na(branches_ID)] # Remove root as NA  
  
  # Remove stratum 0 from the master table (origin tips)
  master_table_filtered <- master_table[master_table$stratum != 0, ]
  
  # Filter master table to get time to remove before focal stratum per branch
  time_to_remove_before_df <- master_table_filtered[(master_table_filtered$parent_br %in% branches_ID & master_table_filtered$stratum == (stratum_ID+1)), ]
  # Filter master table to get time to remove after focal stratum per branch
  time_to_remove_after_df <- master_table_filtered[(master_table_filtered$parent_br %in% branches_ID & master_table_filtered$stratum == (stratum_ID-1)), ]
  
  # View(time_to_remove_before_df)
  # View(time_to_remove_after_df)
  
  # Keep only branch from the focal stratum in the state maps
  state_maps_i <- state_maps[branches_ID]
  names(state_maps_i) <- branches_ID
  state_maps_i
  
  ## Remove time in previous/older strata
  
  if(nrow(time_to_remove_before_df) > 0)
  {
    # Compute cumulative time from root
    states_maps_cumtime_from_root_i <- lapply(X = state_maps_i, FUN = cumsum)
    
    # Remove time before 
    for (i in 1:nrow(time_to_remove_before_df))
    {
      # i <- 7
      branch_i <- time_to_remove_before_df[i,]
      branch_original_ID <- branch_i$parent_br # Which branch
      
      # How much time to remove
      branch_time_before <- branch_i$SUBedge.length # Age of branch if tree stops after this strata
      
      branch_ID <- which(names(states_maps_cumtime_from_root_i) == as.character(branch_original_ID))
      
      # Extract nb of states
      nb_states <- length(state_maps_i[[branch_ID]])
      
      # Find where to cut
      time_cut_location <- which.max(states_maps_cumtime_from_root_i[[branch_ID]] - branch_time_before > 0)
      # Remove residence time before cut (set to zero)
      state_maps_i[[branch_ID]][1:(time_cut_location-1)] <- 0
      # Use cum time since root for the cut
      state_maps_i[[branch_ID]][time_cut_location] <- states_maps_cumtime_from_root_i[[branch_ID]][time_cut_location] - branch_time_before
      # Keep all after the cut
      state_maps_i[[branch_ID]]
    }
    state_maps_i
  }
  
  # Reverse state order (from tip to root)
  state_maps_i <- lapply(X = state_maps_i, FUN = rev) 
  
  ## Remove time in next/younger strata
  
  if(nrow(time_to_remove_after_df) > 0)
  {
    # Compute cumulative time from tips
    states_maps_cumtime_from_tips_i <- lapply(X = state_maps_i, FUN = cumsum)
    
    # Remove time after 
    for (i in 1:nrow(time_to_remove_after_df))
    {
      # i <- 1
      branch_i <- time_to_remove_after_df[i,]
      branch_original_ID <- branch_i$parent_br # Which branch
      
      # How much time to remove
      full_branch_time <- branch_i$edge.length # Full branch length
      focal_stratum_branch_time <- master_table_filtered$SUBedge.length[master_table_filtered$parent_br == branch_original_ID & master_table_filtered$stratum == stratum_ID] # Age of branch if tree stops at the focal strata
      focal_stratum_branch_time <- focal_stratum_branch_time[!is.na(focal_stratum_branch_time)] # Remove root as NA
      branch_time_after <- full_branch_time - focal_stratum_branch_time # Branch length after the focal stratum
      
      branch_ID <- which(names(states_maps_cumtime_from_tips_i) == as.character(branch_original_ID))
      
      # Extract nb of states
      nb_states <- length(state_maps_i[[branch_ID]])
      
      # Find where to cut
      time_cut_location <- which.max(states_maps_cumtime_from_tips_i[[branch_ID]] - branch_time_after > 0)
      # Remove residence time after cut (set to zero)
      state_maps_i[[branch_ID]][1:(time_cut_location-1)] <- 0
      # Use cum time since root for the cut
      state_maps_i[[branch_ID]][time_cut_location] <- states_maps_cumtime_from_tips_i[[branch_ID]][time_cut_location] - branch_time_after
      # Keep all before the cut
      state_maps_i[[branch_ID]]
    }
  }
  
  # Reverse again state order (from root to tip)
  state_maps_i <- lapply(X = state_maps_i, FUN = rev) 
  
  # print(state_maps_i)
  
  # Extract full range list
  returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = model_fit$inputs, include_null_range = model_fit$inputs$include_null_range)
  ranges_list <- returned_mats$ranges_list
  
  # Compute residence times per states
  edge_states_matrix <- plyr::ldply(.data = state_maps_i, rbind)
  row.names(edge_states_matrix) <- edge_states_matrix$.id
  edge_states_matrix <- edge_states_matrix[, -1, drop = FALSE]
  state_times <- apply(X = edge_states_matrix, MARGIN = 2, FUN = sum, na.rm = T)
  state_times_reordered <- setNames(object = rep(0, times = length(ranges_list)), nm = ranges_list)
  state_times_reordered[ranges_list] <- state_times[ranges_list]
  names(state_times_reordered) <- ranges_list
  state_times_reordered[is.na(state_times_reordered)] <- 0
  
  # Add total
  state_times_df <- as.data.frame(t(state_times_reordered))
  state_times_df$total <- sum(state_times_reordered)
  
  # Add percentage of times
  state_times_perc <- state_times_df / state_times_df$total * 100
  state_times_df <- rbind(state_times_df, state_times_perc)
  row.names(state_times_df) <- c("raw", "perc")
  
  # Return residence times x states for the focal time stratum
  return(state_times_df)
}

compute_residence_times_for_all_time_strata <-  function (simmap, model_fit)
{
  # Extract number of time-strata
  nb_strata <- length(model_fit$inputs$timeperiods)
  
  # Extract full range list
  returned_mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = model_fit$inputs, include_null_range = model_fit$inputs$include_null_range)
  ranges_list <- returned_mats$ranges_list
  
  # Initiate final array
  state_times_per_strata <- array(data = NA,
                                  dim = c(2, length(ranges_list) + 1, nb_strata),
                                  dimnames = list(c("raw", "perc"), c(ranges_list, "total"), paste0("Stratum_", 1:nb_strata)))
  dim(state_times_per_strata)
  
  # Loop per strata
  for (i in 1:nb_strata)
  {
    state_times_df_i <- compute_residence_times_by_time_stratum(simmap = simmap, model_fit = model_fit, stratum_ID = i)
    
    state_times_per_strata[,,i] <- as.matrix(state_times_df_i)
  }
  
  state_times_per_strata
  
  # Return final array
  return(state_times_per_strata)
}


