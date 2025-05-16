##### Functions to compute origins of lineages along time steps #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Compute origins of lineages along time steps: dispersal vs. in situ speciation

###

### Inputs

# Simmaps of state evolution from Biogeographic Stochastic Maps

###

### Outputs

# Array of endemicity scores x bioregions x time steps x BSM maps 

###


##### 1/ Function to extract current states at a given focal time from a unique simmap  ####

find_edge_states_at_focal_time_on_Simmap <- function (simmap, focal_time)
{
  # Check validity of focal time
  root_age <- max(phytools::nodeHeights(tree = simmap))
  if (focal_time > root_age)
  {
    stop("Focal time should be equal or lower than root age. Root age = ", root_age,". Focal time = ", focal_time,".")
  }
  if (focal_time < 0)
  {
    stop("Focal time should be equal or higher than 0. Focal time = ", focal_time,".")
  }
  
  # Extract state maps
  initial_state_maps <- simmap$maps
  state_maps <- initial_state_maps
  
  # Extract edge ages
  edges_data_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = simmap) + root_age, 5))
  names(edges_data_df) <- c("rootward_age", "tipward_age")
  edges_data_df$edge_ID <- 1:nrow(edges_data_df)
  edges_data_df <- edges_data_df[, c("edge_ID", "rootward_age", "tipward_age")]
  
  # Detect edges crossing the focal_time
  edges_data_df$focal_edge <- (edges_data_df$rootward_age > focal_time) & (edges_data_df$tipward_age <= focal_time)
  # table(edges_data_df$focal_edge)
  
  # Record time to remove tipward to reach the focal time
  edges_data_df$tipward_diff <- sapply(X = edges_data_df$tipward_age, FUN = function (x) { max(0, focal_time - x) } )
  
  # Reverse state order (from tip to root)
  state_maps <- lapply(X = state_maps, FUN = rev) 
  
  ## Remove tip ward time beyond focal time from state maps
  edges_data_df$current_state <- NA
  
  # Convert state maps as compute cumulative time from tips
  state_maps_cumtime_from_tips <- lapply(X = state_maps, FUN = cumsum)
  
  # Remove time after focal time
  for (i in which(edges_data_df$focal_edge))
  {
    # i <- 128
    
    edge_ID <- as.character(i)
    
    # How much time to remove
    tipward_diff_i <- edges_data_df$tipward_diff[i]
    
    # Extract nb of states
    nb_states <- length(state_maps[[i]])
    
    # Find current state at focal time
    time_cut_location <- which.max(state_maps_cumtime_from_tips[[i]] - tipward_diff_i > 0)
    current_state <- names(time_cut_location)
    
    # Store current state
    edges_data_df$current_state[i] <- current_state
  }
  
  ## Filter to keep only current egdes
  current_edges_data_df <- edges_data_df %>% 
    filter(!is.na(current_state)) %>% 
    dplyr::select(edge_ID, current_state)
  
  # Add time
  current_edges_data_df$time <- focal_time
  
  # Output = df with edge ID and current states
  return(current_edges_data_df)
}


##### 2/ Function to extract state origin per areas at a given focal time from a unique simmap ####

find_egde_origins_at_focal_time_on_Simmap <- function (simmap, focal_time, areas_list = NULL)
{
  # Check validity of focal time
  root_age <- max(phytools::nodeHeights(tree = simmap))
  if (focal_time > root_age)
  {
    stop("Focal time should be equal or lower than root age. Root age = ", root_age,". Focal time = ", focal_time,".")
  }
  if (focal_time < 0)
  {
    stop("Focal time should be equal or higher than 0. Focal time = ", focal_time,".")
  }
  
  # Extract state maps
  initial_state_maps <- simmap$maps
  state_maps <- initial_state_maps
  
  # Extract list of areas recorded in state maps
  ranges_list <- unique(names(unlist(state_maps)))
  ranges_list <- ranges_list[order(ranges_list)]
  areas_recorded_list <- unique(unlist(str_split(string = ranges_list, pattern = "")))
  areas_recorded_list <- areas_list[order(areas_recorded_list)]
  
  # Replace areas_list if not provided
  if (is.null(areas_list))
  {
    areas_list <- areas_recorded_list
  } else {
    # If provided, check adequacy with recorded areas
    missing_areas <- areas_recorded_list[!(areas_recorded_list %in% areas_list)]
    empty_areas <- areas_list[!(areas_list %in% areas_recorded_list)]
    if (length(missing_areas) > 0)
    {
      warning(paste0(length(missing_areas), " area(s) recorded in the state maps are missing from areas_list: ", paste(missing_areas, collapse = ", "),"."))
    }
    if (length(empty_areas) > 0)
    {
      warning(paste0(length(empty_areas), " area(s) provided in areas_list are not recorded in the state maps: ", paste(empty_areas, collapse = ", "),"."))
    }
  }
  
  # Extract root state and ID
  root_node_ID <- phangorn::getRoot(simmap)
  first_edges_ID <- which(simmap$edge[,1] == root_node_ID)
  first_edges_states <- names(unlist(state_maps[first_edges_ID]))
  root_areas <- unique(unlist(str_split(string = first_edges_states, pattern = ""))) # Unique areas
  root_areas <- root_areas[rank(match(x = root_areas, table = areas_list))] # Order them as in areas_list
  root_state <- paste(root_areas, collapse = "") # Collapse to unique range
  
  # Extract edge ages
  edges_data_df <- as.data.frame(round(-1 * phytools::nodeHeights(tree = simmap) + root_age, 5))
  names(edges_data_df) <- c("rootward_age", "tipward_age")
  edges_data_df$edge_ID <- 1:nrow(edges_data_df)
  edges_data_df <- edges_data_df[, c("edge_ID", "rootward_age", "tipward_age")]
  
  # Detect edges crossing the focal_time
  edges_data_df$focal_edge <- (edges_data_df$rootward_age > focal_time) & (edges_data_df$tipward_age <= focal_time)
  # table(edges_data_df$focal_edge)
  
  # Regular case, with edges crossing the focal_time
  if (any(edges_data_df$focal_edge))
  {
    # Record time to remove tipward to reach the focal time
    edges_data_df$tipward_diff <- sapply(X = edges_data_df$tipward_age, FUN = function (x) { max(0, focal_time - x) } )
    
    # Reverse state order (from tip to root)
    state_maps <- lapply(X = state_maps, FUN = rev) 
    
    ## Find current state
    edges_data_df$current_state <- NA
    
    # Convert state maps as compute cumulative time from tips
    state_maps_cumtime_from_tips <- lapply(X = state_maps, FUN = cumsum)
    
    # Remove tipward time beyond focal time from state maps
    for (i in which(edges_data_df$focal_edge))
    {
      # i <- 128
      
      edge_ID <- as.character(i)
      
      # How much time to remove
      tipward_diff_i <- edges_data_df$tipward_diff[i]
      
      # Extract nb of states
      nb_states <- length(state_maps[[i]])
      
      # Find current state at focal time
      time_cut_location <- which.max(state_maps_cumtime_from_tips[[i]] - tipward_diff_i > 0)
      current_state <- names(time_cut_location)
      
      # Store current state
      edges_data_df$current_state[i] <- current_state
      
      # Remove residence time after cut
      state_maps[[i]] <- state_maps[[i]][time_cut_location:length(state_maps[[i]])]
      # Use cum time since root to update residence time of current state
      state_maps[[i]][1] <- state_maps_cumtime_from_tips[[i]][time_cut_location] - tipward_diff_i
    }
    
    ## Find parental state 
    edges_data_df$parental_state <- NA
    # edges_data_df$event_type <- NA
    for (i in which(edges_data_df$focal_edge))
    {
      # i <- 13
      # i <- 128
      # i <- 2980
      
      # # Case of lineage with anagenetic events => retrieve previous state from the state maps of egde i
      # if (length(state_maps[[i]]) > 1)
      # {
      #   previous_state_i <- names(state_maps[[i]][2])
      #   edges_data_df$event_type[i] <- "ana"
      #   
      # } else { # Case of lineage with no anagenetic events => retrieve previous state from the state maps of the parental egde
      #   
      #   # Identify parental edge
      #   rootward_node_ID_i <- simmap$edge[i,1]
      #   
      #   # Special case of edge branching from the root
      #   if (rootward_node_ID_i == root_node_ID)
      #   {
      #     previous_state_i <- root_state
      #   } else {
      #     # Retrieve parental node ID
      #     parental_egde_ID_i <- which(simmap$edge[,2] == rootward_node_ID_i)
      #     
      #     # Retrieve last state of parental edge
      #     previous_state_i <- names(state_maps[[parental_egde_ID_i]][1])
      #     
      #   }
      #   
      #   # Store type of event
      #   edges_data_df$event_type[i] <- "clado"
      # }
      
      # Identify parental edge
      rootward_node_ID_i <- simmap$edge[i,1]
      
      # Special case of edge branching from the root
      if (rootward_node_ID_i == root_node_ID)
      {
        parental_state_i <- root_state
      } else {
        # Retrieve parental node ID
        parental_egde_ID_i <- which(simmap$edge[,2] == rootward_node_ID_i)
        
        # Retrieve last state of parental edge
        parental_state_i <- names(state_maps[[parental_egde_ID_i]][1])
        
      }
      
      # Store parental state
      edges_data_df$parental_state[i] <- parental_state_i
    }
    
    ## Filter to keep only current egdes
    current_edges_data_df <- edges_data_df %>% 
      filter(!is.na(current_state)) %>% 
      dplyr::select(edge_ID, rootward_age, tipward_age, current_state, parental_state)
    
    ## Compute endemicity_score per lineage
    
    current_edges_data_df$endemicity_score <- NA
    # Proportion of current unique areas found among the previous unique areas
    for (i in 1:nrow(current_edges_data_df))
    {
      # i <- 43
      
      # Extract current and previous state
      current_state_i <- current_edges_data_df$current_state[i]
      parental_state_i <- current_edges_data_df$parental_state[i]
      
      # Split in unique areas
      current_unique_areas_i <- unlist(str_split(string = current_state_i, pattern = ""))
      previous_unique_areas_i <- unlist(str_split(string = parental_state_i, pattern = ""))
      
      # Compute endemicity score
      current_match_i <- current_unique_areas_i %in% previous_unique_areas_i
      endemicity_score_i <- sum(current_match_i) / length(current_match_i)
      
      # Store endemicity score
      current_edges_data_df$endemicity_score[i] <- endemicity_score_i
    }
    
    ## Compute endemicity scores per current states
    endemicity_score_df <- current_edges_data_df %>% 
      group_by(current_state) %>% 
      summarize(endemicity_raw = mean(endemicity_score),
                nb_lineages = dplyr::n())
    
    ## Aggregate per unique bioregions
    endemicity_score_df$nb_areas <- nchar(endemicity_score_df$current_state)
    endemicity_score_df$weights <- endemicity_score_df$nb_lineages / endemicity_score_df$nb_areas
    
    # Split ranges in unique areas
    recorded_areas_endemicity_score_df <- data.frame(area = areas_recorded_list)
    recorded_areas_endemicity_score_df$endemicity <- NA
    recorded_areas_endemicity_score_df$area_richness <- 0
    for (i in seq_along(areas_recorded_list))
    {
      # i <- 1
      
      # Extract focal area
      focal_area_i <- areas_recorded_list[i]
      
      # Find matching ranges/states
      focal_states_i <- endemicity_score_df$current_state[str_detect(string = endemicity_score_df$current_state, pattern = focal_area_i)]
      focal_ID_i <- which(endemicity_score_df$current_state %in% focal_states_i)
      
      # Only if at least 1 lineage is present
      if (length(focal_states_i) > 0)
      {
        # Compute weighted mean using endemicity and the weights
        weighted_endemicity_i <- weighted.mean(x = endemicity_score_df$endemicity_raw[focal_ID_i], w = endemicity_score_df$weights[focal_ID_i])
        # Compute current richness
        area_richness_i <- sum(endemicity_score_df$weights[focal_ID_i])
        
        # Store area endemicity and richness
        recorded_areas_endemicity_score_df$endemicity[i] <- weighted_endemicity_i
        recorded_areas_endemicity_score_df$area_richness[i] <- area_richness_i
      }
    }
    
    # Reorder as in areas_list
    areas_endemicity_score_df <- data.frame(area = areas_list)
    areas_endemicity_score_df <- left_join(x = areas_endemicity_score_df, y = recorded_areas_endemicity_score_df, join_by(area == area))
    
    # Add total endemicity
    total_endemicity_score <- mean(current_edges_data_df$endemicity_score)
    total_richness <- nrow(current_edges_data_df)
    areas_endemicity_score_df <- rbind(areas_endemicity_score_df, data.frame(area = "total", endemicity = total_endemicity_score, area_richness = total_richness))
    
    # Add time
    areas_endemicity_score_df$time <- focal_time
    
    # Output = Endemicity score + edge df with edge ID, ages, current state, previous state, endemicity
    output <- list(areas_endemicity_score_df = areas_endemicity_score_df, edges_data_df = current_edges_data_df)
    
  } else {  # Special case when focal_time is the root_time
    
    # Extract root state/areas
    root_areas <- unlist(str_split(string = root_state, pattern = ""))
    nb_root_areas <- length(root_areas)
    root_areas_ID <- which(areas_list %in% root_areas)
    
    # Create areas_endemicity_score_df 
    areas_endemicity_score_df <- data.frame(area = areas_list)
    # Assume endemicity of root state
    areas_endemicity_score_df$endemicity <- NA
    areas_endemicity_score_df$endemicity[root_areas_ID] <- 1
    # Share weights among unique areas
    areas_endemicity_score_df$area_richness <- 0
    areas_endemicity_score_df$area_richness[root_areas_ID] <- 1/nb_root_areas
    
    # Add total endemicity
    areas_endemicity_score_df <- rbind(areas_endemicity_score_df, data.frame(area = "total", endemicity = 1, area_richness = 1))
    
    # Add time
    areas_endemicity_score_df$time <- focal_time
    
    # Do not provide an edge_df
    output <- list(areas_endemicity_score_df = areas_endemicity_score_df, edges_data_df = NULL)
    
  }
  
  # Output = Endemicity score + edges_data_df
  return(output)
}


##### 3/ Function to extract state origin per areas for multiple time steps from a unique simmap ####

find_egde_origins_on_time_steps_on_Simmap <- function (simmap,
                                                       time_steps,
                                                       areas_list = NULL,
                                                       melted = F,  # Should the final array be melted into a df?
                                                       verbose = F)  
{
  
  ## Check validity of time steps
  root_age <- max(phytools::nodeHeights(tree = simmap))
  oldest_time_step <- max(time_steps)
  lowest_focal_time <- min(time_steps)
  if (oldest_time_step > root_age)
  {
    stop("Focal time should be equal or lower than root age. Root age = ", root_age,". Oldest focal time = ", oldest_time_step,".")
  }
  if (lowest_focal_time < 0)
  {
    stop("Focal time should be equal or higher than 0. Lowest time = ", lowest_focal_time,".")
  }
  
  ## Initiate final array
  
  # Extract state maps
  initial_state_maps <- simmap$maps
  
  # Extract list of areas recorded in state maps
  ranges_list <- unique(names(unlist(initial_state_maps)))
  ranges_list <- ranges_list[order(ranges_list)]
  areas_recorded_list <- unique(unlist(str_split(string = ranges_list, pattern = "")))
  areas_recorded_list <- areas_list[order(areas_recorded_list)]
  
  # Replace areas_list if not provided
  if (is.null(areas_list))
  {
    areas_list <- areas_recorded_list
  } else {
    # If provided, check adequacy with recorded areas
    missing_areas <- areas_recorded_list[!(areas_recorded_list %in% areas_list)]
    empty_areas <- areas_list[!(areas_list %in% areas_recorded_list)]
    if (length(missing_areas) > 0)
    {
      warning(paste0(length(missing_areas), " area(s) recorded in the state maps are missing from areas_list: ", paste(missing_areas, collapse = ", "),"."))
    }
    if (length(empty_areas) > 0)
    {
      warning(paste0(length(empty_areas), " area(s) provided in areas_list are not recorded in the state maps: ", paste(empty_areas, collapse = ", "),"."))
    }
  }
  
  # Initiate final array
  areas_endemicity_all_time_steps_array <- array(data = NA,
                                                 dim = c(length(areas_list) + 1, length(time_steps)),
                                                 dimnames = list(c(areas_list, "total"), time_steps))
  
  # Loop per time steps
  for (i in seq_along(time_steps))
  {
    # i <- 114
    
    # Extract focal time
    focal_time_i <- time_steps[i]
    
    # Compute areas endemicity for focal time i
    areas_endemicity_time_i <- find_egde_origins_at_focal_time_on_Simmap(simmap = simmap,
                                                                         focal_time = focal_time_i,
                                                                         areas_list = areas_list)
    areas_endemicity_time_i <- areas_endemicity_time_i$areas_endemicity_score_df
    
    # Store endemicity scores
    areas_i <- areas_endemicity_time_i$area
    areas_endemicity_all_time_steps_array[areas_i, as.character(focal_time_i)] <- areas_endemicity_time_i$endemicity
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Endemicity per bioregions computed for Time = ",focal_time_i," - n°", i, "/", length(time_steps),"\n"))
    }
  }
  
  output <- areas_endemicity_all_time_steps_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(areas_endemicity_all_time_steps_array)
    names(output) <- c("area", "time", "endemicity")
  }
  
  # Export result
  return(output)
}


##### 4/ Function to extract state origin per areas for a given focal time from multiple Simmaps ####

find_egde_origins_at_focal_time_on_MultiSimmap <- function (multiSimmap,
                                                            focal_time,
                                                            areas_list = NULL,
                                                            melted = F,  # Should the final array be melted into a df?
                                                            verbose = T)  
{
  
  ## Check validity of focal time
  root_age <- max(phytools::nodeHeights(tree = multiSimmap[[1]]))
  if (focal_time > root_age)
  {
    stop("Focal time should be equal or lower than root age. Root age = ", root_age,". Focal time = ", focal_time,".")
  }
  if (focal_time < 0)
  {
    stop("Focal time should be equal or higher than 0. Focal time = ", focal_time,".")
  }
  
  ## Initiate final array
  
  # Extract ranges/states recorded in state maps
  ranges_list <- unique(unlist(lapply(X = multiSimmap, FUN = function (x) { colnames(x$mapped.edge) })))
  ranges_list <- ranges_list[order(ranges_list)]
  ranges_list <- ranges_list[ranges_list != "_"]
  areas_recorded_list <- unique(unlist(str_split(string = ranges_list, pattern = "")))
  areas_recorded_list <- areas_list[order(areas_recorded_list)]
  
  # Replace areas_list if not provided
  if (is.null(areas_list))
  {
    areas_list <- areas_recorded_list
  } else {
    # If provided, check adequacy with recorded areas
    missing_areas <- areas_recorded_list[!(areas_recorded_list %in% areas_list)]
    empty_areas <- areas_list[!(areas_list %in% areas_recorded_list)]
    if (length(missing_areas) > 0)
    {
      warning(paste0(length(missing_areas), " area(s) recorded in the state maps are missing from areas_list: ", paste(missing_areas, collapse = ", "),"."))
    }
    if (length(empty_areas) > 0)
    {
      warning(paste0(length(empty_areas), " area(s) provided in areas_list are not recorded in the state maps: ", paste(empty_areas, collapse = ", "),"."))
    }
  }
  
  # Paste maps
  nb_maps <- length(multiSimmap)
  map_labels <- paste0("Map_", 1:nb_maps)
  
  # Initiate final array
  areas_endemicity_all_maps_array <- array(data = NA,
                                           dim = c(length(areas_list) + 1, nb_maps),
                                           dimnames = list(c(areas_list, "total"), map_labels))
  
  # Loop per maps
  for (i in seq_along(multiSimmap))
  {
    # i <- 1
    
    # Compute areas endemicity for map i
    areas_endemicity_map_i <- find_egde_origins_at_focal_time_on_Simmap(simmap = multiSimmap[[i]],
                                                                        focal_time = focal_time,
                                                                        areas_list = areas_list) 
    
    areas_endemicity_map_i <- areas_endemicity_map_i$areas_endemicity_score_df
    
    # Store endemicity scores
    areas_i <- areas_endemicity_map_i$area
    areas_endemicity_all_maps_array[areas_i, i] <- areas_endemicity_map_i$endemicity
    
    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Endemicity per bioregions computed for Map n°", i, "/", length(multiSimmap),"\n"))
    }
  }
  
  output <- areas_endemicity_all_maps_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(areas_endemicity_all_maps_array)
    names(output) <- c("area", "map", "endemicity")
  }
  
  # Export result
  return(output)
}


##### 5/ Function to extract state origin per areas for multiple time steps from multiple Simmaps ####

find_egde_origins_on_time_steps_on_MultiSimmap <- function (multiSimmap,
                                                            time_steps,
                                                            areas_list = NULL,
                                                            melted = F,  # Should the final array be melted into a df?
                                                            verbose = T)  
{
  
  ## Check validity of focal time
  root_age <- max(phytools::nodeHeights(tree = multiSimmap[[1]]))
  oldest_time_step <- max(time_steps)
  lowest_focal_time <- min(time_steps)
  if (oldest_time_step > root_age)
  {
    stop("Focal time should be equal or lower than root age. Root age = ", root_age,". Oldest focal time = ", oldest_time_step,".")
  }
  if (lowest_focal_time < 0)
  {
    stop("Focal time should be equal or higher than 0. Lowest time = ", lowest_focal_time,".")
  }
  
  ## Initiate final array
  
  # Extract ranges/states recorded in state maps
  ranges_list <- unique(unlist(lapply(X = multiSimmap, FUN = function (x) { colnames(x$mapped.edge) })))
  ranges_list <- ranges_list[order(ranges_list)]
  ranges_list <- ranges_list[ranges_list != "_"]
  areas_recorded_list <- unique(unlist(str_split(string = ranges_list, pattern = "")))
  areas_recorded_list <- areas_list[order(areas_recorded_list)]
  
  # Replace areas_list if not provided
  if (is.null(areas_list))
  {
    areas_list <- areas_recorded_list
  } else {
    # If provided, check adequacy with recorded areas
    missing_areas <- areas_recorded_list[!(areas_recorded_list %in% areas_list)]
    empty_areas <- areas_list[!(areas_list %in% areas_recorded_list)]
    if (length(missing_areas) > 0)
    {
      warning(paste0(length(missing_areas), " area(s) recorded in the state maps are missing from areas_list: ", paste(missing_areas, collapse = ", "),"."))
    }
    if (length(empty_areas) > 0)
    {
      warning(paste0(length(empty_areas), " area(s) provided in areas_list are not recorded in the state maps: ", paste(empty_areas, collapse = ", "),"."))
    }
  }
  
  # Paste maps
  nb_maps <- length(multiSimmap)
  map_labels <- paste0("Map_", 1:nb_maps)
  
  # Initiate final array
  areas_endemicity_all_time_all_maps_array <- array(data = NA,
                                                    dim = c(length(areas_list) + 1, length(time_steps), nb_maps),
                                                    dimnames = list(c(areas_list, "total"), time_steps, map_labels))
  
  # Loop per maps
  for (i in seq_along(multiSimmap))
  {
    # i <- 1
    
    # Compute areas endemicity for map i, for all time steps
    areas_endemicity_map_array_i <- find_egde_origins_on_time_steps_on_Simmap(simmap = multiSimmap[[i]],
                                                                              time_steps = time_steps,
                                                                              areas_list = areas_list,
                                                                              melted = F,  # Should the final array be melted into a df?
                                                                              verbose = F) 
    
    # areas_endemicity_map_array_i
    
    # Store endemicity scores
    areas_endemicity_all_time_all_maps_array[, , i] <- areas_endemicity_map_array_i
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Endemicity per bioregions computed for Map n°", i, "/", length(multiSimmap),"\n"))
    }
  }
  
  output <- areas_endemicity_all_time_all_maps_array
  
  ### Melt to df if needed
  if (melted)
  {
    output <- reshape2::melt(areas_endemicity_all_time_all_maps_array)
    names(output) <- c("area", "time", "map", "endemicity")
  }
  
  # Export result
  return(output)
}


