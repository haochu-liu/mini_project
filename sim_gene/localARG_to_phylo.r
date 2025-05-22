localARG_to_phylo <- function(tree, label=FALSE) {
  if (!inherits(tree, "localARG")) {
    stop("Object must be of class 'localARG'")
  }

  # find the node1 that only occurs once
  edge_df <- tree$edge[, 1:3]
  counts_node1 <- table(edge_df$node1)
  single_node1 <- as.numeric(names(counts_node1[counts_node1 == 1]))
  single_node_edge <- which(edge_df$node1 %in% single_node1)

  if (length(single_node_edge)) {
    # for loop to modify the edges
    for (i in single_node_edge) {

    }
  }

}


vec <- local_tree_matrix[, 1]
counts <- table(vec)
single_occurrence_numbers <- as.numeric(names(counts[counts == 1]))
indices_of_single_occurrence <- which(vec %in% single_occurrence_numbers)
