tree_trajectory <- function(child, tree) {
  trajectory <- c(child)
  edge_index <- c()
  current_node <- child
  # complexity: O(n^2) in worst cases
  while (current_node != (tree$n + 1)) {
    # Find the parent of the given node
    parent_to_child <- which(tree$edge[, 2] == current_node)
    parent <- tree$edge[parent_to_child, 1]
    trajectory <- c(trajectory, parent)
    edge_index <- c(edge_index, parent_to_child)
    current_node <- parent
  }
  
  traj <- list(trajectory=trajectory, edge_index=edge_index)
  return(traj) # return a list including trajectory and indices for tree$edge
}


tree_trajectory2 <- function(child, tree) {
  #' input: child node, tree with tree$edge not changing order.
  if (!inherits(tree, "simtree")) {
    stop("Object must be of class 'simtree'")
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
