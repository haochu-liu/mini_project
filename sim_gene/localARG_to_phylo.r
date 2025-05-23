localARG_to_phylo <- function(tree, label=FALSE) {
  if (!inherits(tree, "localARG")) {
    stop("Object must be of class 'localARG'")
  }

  # find the node1 that only occurs once
  edge_df <- tree$edge[, 1:3]
  counts_node1 <- table(edge_df$node1)
  single_node1 <- as.numeric(names(counts_node1[counts_node1 == 1]))

  if (length(single_node1)) {
    # for loop to modify the edges
    for (i in single_node1) {
      delete_edge <- which(edge_df$node1 %in% i)
      target_edge <- which(edge_df$node2 %in% i)
      edge_df$node2[target_edge] <- edge_df$node2[delete_edge]
      edge_df$length[target_edge] <- edge_df$length[target_edge] +
                                     edge_df$length[delete_edge]
      edge_df[delete_edge, ] <- NA
    }
    edge_df <- na.omit(edge_df)
  }

  # update the node index for phylo object
  unique_node1 <- sort(unique(edge_df$node1))
  for (i in 1:length(unique_node1)) {
    edge_df[, 1:2] <- replace(edge_df[, 1:2],
                              edge_df[, 1:2] == unique_node1[i],
                              2*tree$n - i)
  }

  if (label & !is.null(tree$node$gene_str)) {
    leaf_labels <- c()
    for (i in 1:tree$n) {
      node_str <- paste(i, tree$node$gene_str[i], sep=":")
      leaf_labels <- c(leaf_labels, node_str)
    }
  } else {
    leaf_labels <- as.character(1:tree$n)
  }

  # convert the localARG object to phylo object for ape::plot.phylo
  tree_phylo <- list(edge=as.matrix(edge_df[, 1:2]),
                     edge.length=edge_df$length,
                     tip.label=leaf_labels,
                     Nnode=as.integer(tree$n - 1))
  class(tree_phylo) <- "phylo"

  return(tree_phylo)
}
