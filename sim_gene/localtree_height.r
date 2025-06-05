tree_height <- function(tree) {
  if (!inherits(tree, "localtree")) {
    stop("Object must be of class 'localtree'")
  }

  # provide the height (time to MRCA) of the tree
  traj <- tree_trajectory2(1, tree)
  return(sum(tree$edge[traj$edge_index, 3]))
}
