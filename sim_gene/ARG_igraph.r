ARG_igraph <- function(ARG) {
  if (!inherits(ARG, "simARG")) {
    stop("Object must be of class 'simARG'")
  }

  n <- nrow(ARG$node)
  layout_coord <- matrix(data=NA, nrow=n, ncol=2)
  leaf_nodes <- unique(ARG$edge$node2[ARG$edge$node2 %in% 1:ARG$n])
  layout_coord[leaf_nodes, 1] <- 1:ARG$n * 3
  layout_coord[leaf_nodes, 2] <- 0
  current_level <- 1

  recomb_nodes <- ARG$edge$node2[duplicated(ARG$edge$node2) | duplicated(ARG$edge$node2,
    fromLast=TRUE)]
  recomb_nodes <- unique(recomb_nodes)
  recomb_parents <- ARG$edge$node1[ARG$edge$node2 %in% recomb_nodes]

  for (i in (ARG$n+1):n) {
    target_node <- ARG$node$index[i]
    if (target_node %in% recomb_parents) {
      if (is.na(layout_coord[i, 1])) {
        base_node <- ARG$edge$node2[which(target_node == ARG$edge$node1)]
        base_node_index <- which(base_node == ARG$node$index)
        target_node <- ARG$edge$node1[ARG$edge$node2 %in% base_node]
        target_node_index <- which(ARG$node$index %in% target_node)
        layout_coord[target_node_index, 1] <- c(-1, 1) + layout_coord[base_node_index, 1]
        layout_coord[target_node_index, 2] <- current_level
        current_level <- current_level + 1
      }
    } else {
      base_node <- ARG$edge$node2[which(target_node == ARG$edge$node1)]
      base_node_index <- which(ARG$node$index %in% base_node)
      target_node_index <- which(ARG$node$index %in% target_node)
      layout_coord[target_node_index, 1] <- mean(layout_coord[base_node_index, 1])
      layout_coord[target_node_index, 2] <- current_level
      current_level <- current_level + 1
    }
  }
  return(layout_coord)
}
