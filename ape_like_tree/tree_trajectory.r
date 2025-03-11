library(ape)


# Generate a random coalescent tree with 10 tips
tree <- rcoal(10)

# Function to find the trajectory from a tip to the root
tree_trajectory <- function(tip, edge) {
  trajectory <- c(tip)
  current_node <- tip
  
  # Traverse the tree upward until the root is reached
  while (current_node != (Ntip(tree)+1)) {
    # Find the parent of the current node
    parent <- edge[edge[,2] == current_node, 1]
    trajectory <- c(trajectory, parent)
    current_node <- parent
  }
  
  return(trajectory)
}

tree_trajectory(1, edge=tree$edge)
tree$tip.label
plot(tree)
nodelabels()
