source("ape_like_tree/tree_trajectory.r")


tree_length <- function(tree) {
  if (!inherits(tree, "simtree")) {
    stop("Object must be of class 'simtree'")
  }

  # provide the total length of the tree
  return(sum(tree$edge.length))
}


tree_height <- function(tree) {
  if (!inherits(tree, "simtree")) {
    stop("Object must be of class 'simtree'")
  }

  # provide the height (time to MRCA) of the tree
  traj <- tree_trajectory(1, tree)
  return(sum(tree$edge.length[traj$edge_index]))
}