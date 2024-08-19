
# Modified from phytools::ltt (by Liam Revell)

ltt_simmap_custom <- function (tree, plot = TRUE, log.lineages = FALSE, gamma = TRUE, ...) 
{
  if (hasArg(tol)) 
    tol <- list(...)$tol
  else tol <- .Machine$double.eps^0.5
  is_ultrametric <- is.ultrametric(tree)
  if (!inherits(tree, "simmap")) {
    stop("tree must be an object of class \"simmap\".")
  }
  else {
    
    # Extract states
    levs <- sort(unique(c(getStates(tree, "tips"), getStates(tree, "nodes"))))
    
    # Tree with unique branch status
    tt <- map.to.singleton(tree)
    
    # Extract node height
    H <- nodeHeights(tt)
    h <- c(0, max(H) - BRANCHING(tt, is_ultrametric), TIPHEIGHTS(tt, is_ultrametric))
    
    # Extract branch states as factor
    ss <- setNames(as.factor(names(tt$edge.length)), tt$edge[, 2])
    
    # Bonary matrix of states per branch
    lineages <- matrix(0, length(h), length(levs), dimnames = list(names(h), levs))
    # Fill root state
    lineages[1, getStates(tree, "nodes")[1]] <- 1
    
    # Fill binary states of internal branches
    ROOT <- Ntip(tt) + 1
    for (i in 2:length(h))
    {
      if (h[i] < max(h))
      {
        ii <- intersect(which(h[i] >= H[, 1]), which(h[i] < H[, 2]))
      } else {
        ii <- intersect(which(h[i] >= H[, 1]), which(h[i] <= H[, 2]))
      } 
      lineages[i, ] <- summary(ss[ii])[names(summary(ss[ii])) %in% levs]
    }
    ii <- order(h)
    times <- h[ii]
    lineages <- lineages[ii, , drop = FALSE]
    lineages <- cbind(lineages, total = rowSums(lineages))
    obj <- list(times = times, ltt = lineages)
    if (gamma == FALSE)
    {
      obj <- list(ltt = lineages, times = times, tree = tree)
      class(obj) <- "ltt.simmap"
    } else {
      gam <- gammatest(ltt(as.phylo(tree), plot = FALSE))
      obj <- list(ltt = lineages, times = times, gamma = gam$gamma, 
                  p = gam$p, tree = tree)
      class(obj) <- "ltt.simmap"
    }
  }
  if (plot) 
    plot(obj, log.lineages = log.lineages, ...)
  obj
}


## For multiSimmap

ltt_multiSimmap_custom <-function (tree, gamma = TRUE, ...) 
{
  if (!inherits(tree, "multiSimmap")) {
    stop("tree must be object of class \"multiSimmap\".")
  }
  else {
    obj <- lapply(tree, ltt_simmap_custom, plot = FALSE, log.lineages = FALSE, 
                  gamma = gamma)
    class(obj) <- "ltt.multiSimmap"
    return(obj)
  }
}

## Internally used functions

BRANCHING <- function(phy,is_ultrametric)
{
  x <- if(is_ultrametric) branching.times(phy)
  else { 
    sort(setNames(max(nodeHeights(phy))-sapply(1:phy$Nnode+Ntip(phy),
                                               nodeheight,tree=phy),1:phy$Nnode+Ntip(phy)))
  }
  x
}

TIPHEIGHTS <- function(phy,is_ultrametric)
{
  x <- if(is_ultrametric) {
    min(setNames(sapply(1:Ntip(phy),nodeheight,tree=phy),1:Ntip(phy)))
  } else setNames(sapply(1:Ntip(phy),nodeheight,tree=phy),1:Ntip(phy))
  x
}
