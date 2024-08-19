# Function to apply the branch lengths of a reference tree to a second tree, including mappings
# written by Liam J. Revell 2012

applyBranchLengths <- function(tree, edge.length)
{
  class_tree <- class(tree)
  if("multiPhylo" %in% class_tree)
  {
    trees <- lapply(X = tree, FUN = applyBranchLengths, edge.length = edge.length)
    class(trees) <- class_tree
    return(trees)
    
  } else {
    
    tree$edge.length <- edge.length
    if(!is.null(tree$maps))
    {
      for(i in 1:nrow(tree$edge))
      {
        temp <- tree$maps[[i]] / sum(tree$maps[[i]])
        tree$maps[[i]] <- temp*tree$edge.length[i]
      }
    }
    
    if(!is.null(tree$mapped.edge)){
      a <- vector()
      
      for (i in 1:nrow(tree$edge))
      {
        a <- c(a, names(tree$maps[[i]]))
      } 
      a <- unique(a)
      tree$mapped.edge <- matrix(data = 0, length(tree$edge.length), length(a), dimnames = list(apply(X = tree$edge, MARGIN = 1, FUN = function(x) {paste(x,collapse=",") } ), state = a))
      
      for(i in 1:length(tree$maps))
      {
        for(j in 1:length(tree$maps[[i]]))
        {
          tree$mapped.edge[i, names(tree$maps[[i]])[j]] <- tree$mapped.edge[i, names(tree$maps[[i]])[j]] + tree$maps[[i]][j]
        }
      }
    }
    
    class(tree) <- class_tree
    return(tree)
  }
}
