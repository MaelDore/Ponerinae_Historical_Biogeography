
built_posterior_phylogenies_from_MCMC_outputs <- function (MCMC_node_ages, phy, time_conversion = 1, keep_metadata = T, verbose = T)
{
  # Get the internal node labels
  internal_nodes_ID <- names(MCMC_node_ages)[str_detect(string = names(MCMC_node_ages), pattern = "t_n")]
  internal_nodes_ID <- str_remove(string = internal_nodes_ID, pattern = "t_n")
  
  # Extract node_ages
  row.names(MCMC_node_ages) <- MCMC_node_ages$Gen
  node_ages_df <- MCMC_node_ages[, str_detect(string = names(MCMC_node_ages), pattern = "t_n")]
  names(node_ages_df) <- internal_nodes_ID
  
  # Convert node ages
  node_ages_df <- node_ages_df * time_conversion
  
  # Extract metadata
  metadata_df <- MCMC_node_ages[, !str_detect(string = names(MCMC_node_ages), pattern = "t_n")]
  
  # Initiate final multiPhylo
  PP_phylo <- list()
  
  ## Loop per posterior samples
  for (i in 1:nrow(node_ages_df))
  {
    # i <- 1
    
    # Initiate new phylogeny
    PP_phy_i <- phy
    
    # Extract node ages for this posterior sample
    node_ages_df_i <- data.frame(node_ID = as.numeric(names(node_ages_df)), node_age = unlist(node_ages_df[i, , drop = T]))
    
    # Compute new branch length per edge
    branch_ages_df <- as.data.frame(phy$edge)
    names(branch_ages_df) <- c("root_node_ID", "tip_node_ID")
    # Extract root age
    root_age <- max(node_ages_df_i$node_age)
    # Extract rootward node ages
    branch_ages_df <- branch_ages_df %>% 
      left_join(y = node_ages_df_i,
                by = join_by("root_node_ID" == "node_ID")) %>%
      rename(root_node_age = node_age)
    # Extract tipward node ages
    branch_ages_df <- branch_ages_df %>% 
      left_join(y = node_ages_df_i,
                by = join_by("tip_node_ID" == "node_ID")) %>%
      rename(tip_node_age = node_age) %>%
      mutate(tip_node_age = tidyr::replace_na(tip_node_age, replace = 0))
    # Compute new branch length per edge
    branch_ages_df$branch_length <- branch_ages_df$root_node_age - branch_ages_df$tip_node_age
    branch_ages_df$branch_length[branch_ages_df$branch_length  < 0] <- 0 # Fix negative branches
    
    # Update branch length
    PP_phy_i$edge.length <- branch_ages_df$branch_length
    
    # Plot initial and new phylogeny
    # plot(phy)
    # plot(PP_phy_i)
    
    # Inform metadata
    if (keep_metadata)
    {
      PP_phy_i$Generation <- metadata_df$Gen[i]
      PP_phy_i$mu <- metadata_df$mu[i]
      PP_phy_i$sigma2 <- metadata_df$sigma2[i]
      PP_phy_i$lnLk <- metadata_df$lnL[i]
    }

    # Store result
    PP_phylo[[i]] <- PP_phy_i
    
    ## Print progress
    if ((i %% 1000 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Phylogeny created from posterior sample n°", i, "/", nrow(node_ages_df),"\n"))
    }
  }
  
  # Assign ape multiPhylo class
  class(PP_phylo) <- "multiPhylo"
  
  # Export output
  return(PP_phylo)
}

