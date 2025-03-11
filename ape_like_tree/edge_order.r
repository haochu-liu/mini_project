library(ape)


# Generate a random coalescent tree with 10 tips
tree <- rcoal(10)

# Function to perform postorder traversal and reorder edge matrix
reorder_cladewise <- function(tree) {
  # Initialize a stack for postorder traversal
  postorder <- c()
  
  # Recursive function for postorder traversal
  postorder_traversal <- function(node) {
    children <- tree$edge[tree$edge[,1] == node, 2]
    for (child in children) {
      # Recursively visit children
      postorder_traversal(child)
    }
    # Append node after visiting all children
    postorder <<- c(postorder, node)
  }
  
  # Start traversal from the root node
  root <- tree$n + 1
  postorder_traversal(root)
  
  # Create a new edge matrix in postorder
  new_edge <- matrix(NA, nrow = nrow(tree$edge), ncol = 2)
  edge_index <- 1
  for (node in postorder) {
    children <- tree$edge[tree$edge[,1] == node, 2]
    for (child in children) {
      new_edge[edge_index, ] <- c(node, child)
      edge_index <- edge_index + 1
    }
  }
  
  # Update the tree with the new edge matrix
  tree$edge <- new_edge
  return(tree)
}

# Apply the manual reordering function
tree_cladewise_manual <- reorder_cladewise(tree)

# Compare the result with the built-in reorder()
tree_cladewise_builtin <- reorder(tree, order = "cladewise")

# Check if both edge matrices are identical
identical(tree_cladewise_manual$edge, tree_cladewise_builtin$edge)