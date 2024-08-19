

# Original function written by Liam Revell, 2012 in the phytools package

# Additions 
  # Can modify manually the tolerance to handle issue with mismatch between node ages and time steps used
  # Can print progress across egdes
  # Can provide a manual color scale (N = 1001)


densityMap_custom <- function (trees, res = 100, fsize = NULL, ftype = NULL, lwd = 3, 
                               tol = 1e-5, verbose = T, col_scale = NULL,
                               check = FALSE, legend = NULL, outline = FALSE, type = "phylogram", 
                               direction = "rightwards", plot = TRUE, ...) 
{
  require(phytools)
  
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(0.3, 4)
  if (hasArg(offset)) 
    offset <- list(...)$offset
  else offset <- NULL
  if (hasArg(states)) 
    states <- list(...)$states
  else states <- NULL
  if (hasArg(hold)) 
    hold <- list(...)$hold
  else hold <- TRUE
  if (length(lwd) == 1) 
    lwd <- rep(lwd, 2)
  else if (length(lwd) > 2) 
    lwd <- lwd[1:2]
  
  # Adjust tolerance
  tol <- tol
  
  if (!inherits(trees, "multiPhylo") && inherits(trees, "phylo")) 
    stop("trees not \"multiPhylo\" object; just use plotSimmap.")
  if (!inherits(trees, "multiPhylo")) 
    stop("trees should be an object of class \"multiPhylo\".")
  
  # Extract root age
  h <- sapply(unclass(trees), function(x) max(phytools::nodeHeights(x)))
  
  # Define time steps
  steps <- 0:res/res * max(h)
  
  # Rescale trees to ensure they all have the same root age
  trees <- rescaleSimmap(trees, totalDepth = max(h))
  
  # Check that phylogeny topology and branch length are equal
  if (check) {
    X <- matrix(FALSE, length(trees), length(trees))
    for (i in 1:length(trees)) X[i, ] <- sapply(trees, all.equal.phylo, 
                                                current = trees[[i]])
    if (!all(X)) 
      stop("some of the trees don't match in topology or relative branch lengths")
  }
  
  # Extract first tree as reference
  tree <- trees[[1]]
  
  # Remove class
  trees <- unclass(trees)
  
  # Extract all states from the first tree (dangerous if some states are not present in this tree but in other!)
  if (is.null(states)) 
    ss <- sort(unique(c(getStates(tree, "nodes"), getStates(tree, 
                                                            "tips"))))
  else ss <- states
  
  # If states are not binary, rename the first two states as "0" and "1" (dangerous as if there are more states, will lead to errors)
  if (!all(ss == c("0", "1"))) 
  {
    c1 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    c2 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    trees <- lapply(trees, mergeMappedStates, ss[1], c1)
    trees <- lapply(trees, mergeMappedStates, ss[2], c2)
    trees <- lapply(trees, mergeMappedStates, c1, "0")
    trees <- lapply(trees, mergeMappedStates, c2, "1")
  }
  
  # Extract all node ages per edge
  H <- phytools::nodeHeights(tree)
  message("sorry - this might take a while; please be patient")
  
  # Reinitiate the map of the reference tree with NULL data
  tree$maps <- vector(mode = "list", length = nrow(tree$edge))
  
  # Loop per edge/item in tree$maps
  for (i in 1:nrow(tree$edge))
  {
    # i <- 5
    # i <- 983
    
    # YY = Matrix of ages to use as time step along edge i
    # Include start (0) and end (edge length)
    # One raw = one interval
    # Columns = Start and End relative age
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]),
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) - H[i, 1]
    
    # Initiate vector of final time step values for edge i
    ZZ <- rep(0, nrow(YY))
    
    # Loop per trees/simmaps
    for (j in 1:length(trees))
    {
      # j <- 1
      
      # XX = a matrix of the states detected for edge i on simmap j
      # One row = one state
      # Colmuns = start and end relative age
      XX <- matrix(data = 0, 
                   nrow = length(trees[[j]]$maps[[i]]), 
                   ncol = 2, 
                   dimnames = list(names(trees[[j]]$maps[[i]]), 
                                   c("start", "end")))
      # Fill the first raw with information on start and end relative age of the first state
      XX[1, 2] <- trees[[j]]$maps[[i]][1]
      
      # Case with multiple states: fill information for other states
      if (length(trees[[j]]$maps[[i]]) > 1)
      {
        for (k in 2:length(trees[[j]]$maps[[i]]))
        {
          XX[k, 1] <- XX[k - 1, 2]
          XX[k, 2] <- XX[k, 1] + trees[[j]]$maps[[i]][k]
        }
      }
      
      # Loop per time interval wanted for the density mapping
      for (k in 1:nrow(YY))
      {
        # k <- 1
        
        # Detect which state start before the k time step
        lower <- which(XX[, 1] <= YY[k, 1])
        lower <- lower[length(lower)] # Take the last one as the last state recorded before the beginning of the time step
        
        # Detect which state end after the k time step
        upper <- which(XX[, 2] >= (YY[k, 2] - tol))[1]
        
        
        AA <- 0
        names(lower) <- names(upper) <- NULL
        if (!all(XX == 0)) 
        {
          # Case for internal edge (end time > 0)
          
          # Loop per states on the time interval
          for (l in lower:upper)
          {
            # Compute weighted mean of the time interval
            AA <- AA + (min(XX[l, 2], YY[k, 2]) - max(XX[l, 1], YY[k, 1]))/(YY[k, 2] - YY[k, 1]) * as.numeric(rownames(XX)[l])
          }
        } else {
          # Case for tips (or null branches) (start and end time = 0)
          AA <- as.numeric(rownames(XX)[1]) # Use the tip state
        }
        # Increment the final value of the edge i
        ZZ[k] <- ZZ[k] + AA/length(trees)
      }
    }
    # Record time steps length
    tree$maps[[i]] <- YY[, 2] - YY[, 1]
    
    # Record time steps continuous value as names of the maps of edge i
    names(tree$maps[[i]]) <- round(ZZ * 1000) # Convert proportion in a scale from 0 to 1000
    
    # Print progress every 100 edges
    if ((i %% 100 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probability computed for edge n°", i, "/", nrow(tree$edge),"\n"))
    }
  }
  
  # Create color scale
  if (is.null(col_scale))
  {
    cols <- rainbow(1001, start = 0.7, end = 0)
    names(cols) <- 0:1000
  } else {
    cols <- col_scale
  }

  # Recreate map using continuous values
  tree$mapped.edge <- makeMappedEdge(tree$edge, tree$maps)
  tree$mapped.edge <- tree$mapped.edge[, order(as.numeric(colnames(tree$mapped.edge)))]
  
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  attr(tree, "map.order") <- "right-to-left"
  x <- list(tree = tree, cols = cols, states = ss)
  class(x) <- "densityMap"
  
  # Plot
  if (plot) 
    plot.densityMap(x, fsize = fsize, ftype = ftype, lwd = lwd, 
                    legend = legend, outline = outline, type = type, 
                    mar = mar, direction = direction, offset = offset, 
                    hold = hold)
  
  # Return final simmap
  invisible(x)
}


makeMappedEdge <- function (edge, maps) 
{
  st <- sort(unique(unlist(sapply(maps, function(x) names(x)))))
  mapped.edge <- matrix(0, nrow(edge), length(st))
  rownames(mapped.edge) <- apply(edge, 1, function(x) paste(x, 
                                                            collapse = ","))
  colnames(mapped.edge) <- st
  for (i in 1:length(maps)) for (j in 1:length(maps[[i]])) mapped.edge[i, 
                                                                       names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + 
    maps[[i]][j]
  return(mapped.edge)
}
