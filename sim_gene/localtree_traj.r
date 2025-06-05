tree_trajectory2 <- function(child, tree) {
  #' input: child node, tree with tree$edge not changing order.
  if (!inherits(tree, "localtree")) {
    stop("Object must be of class 'localtree'")
  }

  trajectory <- c(child)
  edge_index <- c()
  current_node <- child
  # complexity: O(n)
  for (i in 1:nrow(tree$edge)) {
    if (tree$edge[i, 2] == current_node) {
      current_node <- tree$edge[i, 1]
      trajectory <- c(trajectory, current_node)
      edge_index <- c(edge_index, i)
    }
  }
  
  traj <- list(trajectory=trajectory, edge_index=edge_index)
  return(traj) # return a list including trajectory and indices for tree$edge
}
