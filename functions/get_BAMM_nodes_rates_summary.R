
### Modified version of the getMarginalBranchRateMatrix from BAMMtools package used to retrieve mean rates at nodes rather than along branches

# Original function from BAMMtools package. Reference: 10.1111/2041-210X.12199

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Original function to compute rates per edges from a BAMM output ####

# Input = BAMM output object

# Output = list of matrices
  # Matrix of speciation rates as edges (rows) x BAMM posteriors (columns)
  # Matrix of extinction rates as edges (rows) x BAMM posteriors (columns) 

getMarginalBranchRateMatrix <- function (ephy, verbose = FALSE) 
{
  if (!inherits(ephy, "bammdata")) {
    stop("Object must be of class bammdata\n")
  }
  
  # Build matrix of speciation rates
  lammat <- matrix(0, ncol = length(ephy$eventBranchSegs), 
                   nrow = nrow(ephy$edge))
  
  # Build matrix of extinction rates
  mumat <- matrix(0, ncol = length(ephy$eventBranchSegs), nrow = nrow(ephy$edge))
  
  # Loop per posterior configurations
  for (i in 1:length(ephy$eventBranchSegs))
  {
    if (verbose)
    {
      cat("Processing sample ", i, "\n")
    }
    
    # Extract edge metadata with vectors of regime membership per branches, tipward node ID and begin/end ages of the branches
      # Strictly, not just edges, but segments if an edge is split among multiple regimes by a shift
      # First column is NOT edge ID, but the ID of the tipward node!
    esegs <- ephy$eventBranchSegs[[i]]
    # Extract regimes with parameters
    events <- ephy$eventData[[i]]
    events <- events[order(events$index), ]
    
    # Extract time since start of the segment regime and the beginning of the segment = relative start time = use to compute start rate
    relsegmentstart <- esegs[, 2] - ephy$eventData[[i]]$time[esegs[, 4]]
    # Extract time since start of the segment regime and the end of the segment = relative end time = use to compute final rate
    relsegmentend <- esegs[, 3] - ephy$eventData[[i]]$time[esegs[, 4]]
    
    # Extract initial rate of speciation (lambda0) of the associated regime
    lam1 <- ephy$eventData[[i]]$lam1[esegs[, 4]]
    # Extract speciation rate change parameter (alpha/z) of the associated regime. If = 0, rate is constant. If < 0 = decay. If > 0 = growth
    lam2 <- ephy$eventData[[i]]$lam2[esegs[, 4]]
    # Extract initial rate of extinction (mu0) of the regime. Should be the constant extinction rate of the regime.
    mu1 <- ephy$eventData[[i]]$mu1[esegs[, 4]]
    # Extract extinction rate change parameter (alpha/z) of the regime. Should be fixed to 0 as extinction rates are constant in BAMM.
    mu2 <- ephy$eventData[[i]]$mu2[esegs[, 4]]
    
    # Compute the integral of speciation rates over the segment time
    lamint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, 
                                       lam1, lam2)
    # Compute the integral of extinction rates over the segment time
    muint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, 
                                      mu1, mu2)
    # Extract segment length
    seglengths <- esegs[, 3] - esegs[, 2]
    
    # Compute mean rates as the integral divided by the edge length
    # Loop per edge to sum integrals and lengths of segments belonging to the same edge
    for (k in 1:nrow(ephy$edge)) 
    {
      # Find segments matching the tipward node ID
      isRightBranch <- esegs[, 1] == ephy$edge[k, 2]
      # Compute edge mean speciation rate as the sum of segment integral / sum of segment lengths
      lammat[k, i] <- sum(lamint[isRightBranch])/sum(seglengths[isRightBranch])
      # Compute edge mean extinction rate as the sum of segment integral / sum of segment lengths
      mumat[k, i] <- sum(muint[isRightBranch])/sum(seglengths[isRightBranch])
    }
  }
  
  # Export matrix of speciation and extinction rates
  if (ephy$type == "diversification") {
    return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat))
  }
  # Export matrix of phenotypic rate of trait evolution
  if (ephy$type == "trait") {
    return(list(beta_branch_matrix = lammat))
  }
}


### 2/ Internal function to compute the integral of rates on a branch segment with respect to time ####

# This is not the average. It is the mean rates based on the integral as Intergal = Area, Time = X-axis, Rates = Y-axis.
# Thus mean rates = mean-Y = Area / X-axis = Integral / Time range

timeIntegratedBranchRate <- function(t1, t2, # Starting and ending time (boudaries for the integral)
                                     p1, # Initial rate parameter in the exponential function
                                     p2) # Rate variation parameter in the exponential function 
{
  tol <- 0.00001;
  res <- vector(mode = 'numeric', length = length(t1));
  
  # Constant rate
  zero <- which(abs(p2) < tol); # Variation rate is close to zero
  p1s <- p1[zero];
  t1s <- t1[zero];
  t2s <- t2[zero];
  res[zero] <- p1s * (t2s - t1s); # Area/Integral is Fixed rate x Time range
  
  # Declining rate
  nonzero <- which(p2 < -tol); # Variation rate is negative
  p1s <- p1[nonzero];
  p2s <- p2[nonzero];
  t1s <- t1[nonzero];
  t2s <- t2[nonzero];
  res[nonzero] <- (p1s/p2s)*(exp(p2s*t2s) - exp(p2s*t1s)); # Area/Integral is based on the exponential function
  
  # Increasing rate
  nonzero <- which(p2 > tol); # Variation rate is positive
  p1s <- p1[nonzero];
  p2s <- p2[nonzero];
  t1s <- t1[nonzero];
  t2s <- t2[nonzero];
  res[nonzero] <- (p1s/p2s)*(2*p2s*(t2s-t1s) + exp(-p2s*t2s) - exp(-p2s*t1s)); # Area/Integral is based on the exponential function
  
  # Return mean rates
  return(res);
}

### 3/ New function to compute rates per nodes from a BAMM output ####

# Input = BAMM output object

# Output = list of matrices
  # Matrix of speciation rates as Nodes (rows) x BAMM posteriors (columns)
  # Matrix of extinction rates as Nodes (rows) x BAMM posteriors (columns) 

getMarginalNodeRateMatrix <- function (BAMM_output, verbose = FALSE) 
{
  if (!inherits(BAMM_output, "bammdata")) {
    stop("Object must be of class bammdata\n")
  }
  
  # Build matrix of speciation rates
  lambda_mat <- matrix(data = NA, ncol = length(BAMM_output$eventBranchSegs), 
                   nrow = nrow(BAMM_output$edge) + 1) # Number of nodes (internal + tips, including root = nb of edges + 1)
  
  # Build matrix of extinction rates
  mu_mat <- matrix(data = NA, ncol = length(BAMM_output$eventBranchSegs),
                  nrow = nrow(BAMM_output$edge) + 1) # Number of nodes (internal + tips, including root = nb of edges + 1)
  
  # Loop per posterior configurations
  for (i in 1:length(BAMM_output$eventBranchSegs))
  {
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Processing BAMM posterior sample n°", i, "/", length(BAMM_output$eventBranchSegs),"\n"))
    }
    
    # Extract regimes with parameters
    regimes_metadata_df <- BAMM_output$eventData[[i]]
    regimes_metadata_df <- regimes_metadata_df[order(regimes_metadata_df$index), ]
    
    # Extract edge metadata with vectors of regime membership per branches, tipward node ID and begin/end ages of the branches
      # Strictly, not just edges, but segments if an edge is split among multiple regimes by a shift
    segs_metadata_df <- as.data.frame(BAMM_output$eventBranchSegs[[i]])
    names(segs_metadata_df) <- c("tipward_node_ID", "edge_rootward_time", "edge_tipward_time", "regime_ID")
    # Assign edge ID
    segs_metadata_df$edge_ID <- match(x = segs_metadata_df$tipward_node_ID, table = BAMM_output$edge[, 2])
    # Filter only the terminal segment of each edge
    edges_metadata_df <- segs_metadata_df %>% 
      group_by(edge_ID) %>%
      arrange(edge_ID, desc(edge_tipward_time)) %>%
      mutate(segment_rank = row_number()) %>%
      filter(segment_rank == 1) %>%
      dplyr::select(edge_ID, tipward_node_ID, edge_tipward_time, regime_ID) %>%
      ungroup()
    
    # Convert to node metadata
    nodes_metadata_df <- edges_metadata_df %>% 
      rename(node_ID = tipward_node_ID,
             parental_edge_ID = edge_ID,
             time_since_root = edge_tipward_time) %>%
      dplyr::select(node_ID, time_since_root, regime_ID, parental_edge_ID) %>%
      arrange(node_ID)
    
    # Add the root node
    root_node_ID <- regimes_metadata_df$node[1]
    root_metadata_df <- data.frame(node_ID = root_node_ID, time_since_root = 0, regime_ID = 1, parental_edge_ID = NA)
    nodes_metadata_df <- nodes_metadata_df %>% 
      rbind(., root_metadata_df) %>% 
      arrange(node_ID)
    
    # Extract time since node and the beginning of the segment = relative time = use to compute node rate
    nodes_metadata_df$time_since_regime_start <- nodes_metadata_df$time_since_root - BAMM_output$eventData[[i]]$time[nodes_metadata_df$regime_ID]

    # Extract initial rate of speciation (lambda0) of the associated regime
    nodes_metadata_df$lambda0 <- BAMM_output$eventData[[i]]$lam1[nodes_metadata_df$regime_ID]
    # Extract speciation rate change parameter (alpha/z) of the associated regime. If = 0, rate is constant. If < 0 = decay. If > 0 = growth
    nodes_metadata_df$lambda_var <- BAMM_output$eventData[[i]]$lam2[nodes_metadata_df$regime_ID]
    # Extract initial rate of extinction (mu0) of the regime. Should be the constant extinction rate of the regime.
    nodes_metadata_df$mu0 <- BAMM_output$eventData[[i]]$mu1[nodes_metadata_df$regime_ID]
    # Extract extinction rate change parameter (alpha/z) of the regime. Should be fixed to 0 as extinction rates are constant in BAMM.
    nodes_metadata_df$mu_var <- BAMM_output$eventData[[i]]$mu2[nodes_metadata_df$regime_ID]
    
    # Compute the node speciation rate
    nodes_metadata_df$lambda <- nodes_metadata_df$lambda0 * exp(nodes_metadata_df$lambda_var * nodes_metadata_df$time_since_regime_start)
    # Compute the node extinction rate
    nodes_metadata_df$mu <- nodes_metadata_df$mu0 * exp(nodes_metadata_df$mu_var * nodes_metadata_df$time_since_regime_start)
    
    # Inform matrices
    lambda_mat[, i] <- nodes_metadata_df$lambda
    mu_mat[, i] <- nodes_metadata_df$mu
    
  }
  
  # Add col.names to output matrices
  colnames(lambda_mat) <- paste0("BAMM_posterior_", 1:length(BAMM_output$eventBranchSegs))
  colnames(mu_mat) <- paste0("BAMM_posterior_", 1:length(BAMM_output$eventBranchSegs))
  
  # Export matrix of speciation and extinction rates
  if (BAMM_output$type == "diversification")
  {
    return(list(lambda_nodes_matrix = lambda_mat, mu_nodes_matrix = mu_mat))
  }
  # Export matrix of phenotypic rate of trait evolution
  if (BAMM_output$type == "trait")
  {
    return(list(beta_nodes_matrix = lambda_mat))
  }
}


### 4/ New function to compute mean/median rates per nodes from a BAMM output ####

# Input = BAMM output object

# Output = Array of speciation, extinction, and net diversification rates per Nodes (rows) x types of rates (columns) x type of summary stats

get_BAMM_nodes_rates_summary <- function (x = BAMM_output, verbose = FALSE)
{
  # Get speciation and extinction rates per nodes across all BAMM posteriors
  nodes_matrix_list <- getMarginalNodeRateMatrix(BAMM_output = BAMM_output, verbose = verbose) 
  nb_nodes <- dim(nodes_matrix_list$lambda_nodes_matrix)[1]
  nb_BAMM_posteriors <- dim(nodes_matrix_list$lambda_nodes_matrix)[2]
  
  # Compute net diversification rates
  nodes_matrix_list[[3]] <- nodes_matrix_list[[1]] - nodes_matrix_list[[2]]
  names(nodes_matrix_list)[3] <- "net_div_nodes_matrix"
  
  dim(nodes_matrix_list[[1]])
  
  # Convert list of marginal node rates to array
  nodes_rates_array <- array(data = NA,
                             dim = c(nb_nodes, nb_BAMM_posteriors, 3),
                             dimnames = list(node_ID = 1:nb_nodes, paste0("BAMM_posterior_", BAMM_posterior_ID = 1:nb_BAMM_posteriors), stats = c("speciation", "extinction", "net diversification")))
  nodes_rates_array[ , ,1] <- nodes_matrix_list$lambda_nodes_matrix
  nodes_rates_array[ , ,2] <- nodes_matrix_list$mu_nodes_matrix
  nodes_rates_array[ , ,3] <- nodes_matrix_list$net_div_nodes_matrix
  
  # Initiate final summary array
  nodes_rates_summary_array <- array(data = NA,
                                     dim = c(nb_nodes, 3, 6),
                                     dimnames = list(node_ID = 1:nb_nodes, rate_types = c("speciation", "extinction", "net diversification"),  stats = c("mean", "median", "Q_2.5", "Q_97.5", "HPD_2.5", "HPD_97.5")))
  
  # Compute summary stats across BAMM posteriors
  nodes_rates_summary_array[ , , "mean"] <- apply(X = nodes_rates_array, MARGIN = c(1,3), FUN = mean)
  nodes_rates_summary_array[ , , "median"] <- apply(X = nodes_rates_array, MARGIN = c(1,3), FUN = median)
  Q95 <- apply(X = nodes_rates_array, MARGIN = c(1,3), FUN = quantile, probs = c(0.025, 0.975))
  nodes_rates_summary_array[ , , c("Q_2.5", "Q_97.5")] <- aperm(Q95, c(2,3,1))
  HPD95 <- apply(X = nodes_rates_array, MARGIN = c(1,3), FUN = BayesTwin::HPD, cred_int = 0.95)
  nodes_rates_summary_array[ , , c("HPD_2.5", "HPD_97.5")] <- aperm(HPD95, c(2,3,1))
  
  # Export output
  return(nodes_rates_summary_array)
}


