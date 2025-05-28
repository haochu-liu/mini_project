#' Input: ARG (simARG)
#' Construct a layout matrix for tree-like ARG
#' Output: layout for igraph::plot.igraph
ARG_igraph <- function(ARG) {
  if (!inherits(ARG, "sim_ISM_ARG") & !inherits(ARG, "sim_FSM_ARG")) {
    stop("Object must be of class 'sim_ISM_ARG' or 'sim_FSM_ARG'")
  }

  # create layout matrix and leaf nodes
  n <- nrow(ARG$node)
  layout_coord <- matrix(data=NA, nrow=n, ncol=2)
  leaf_nodes <- unique(ARG$edge$node2[ARG$edge$node2 %in% 1:ARG$n])
  layout_coord[leaf_nodes, 1] <- 1:ARG$n * 3
  layout_coord[leaf_nodes, 2] <- 0
  current_level <- 1

  # store the recombination nodes
  recomb_nodes <- ARG$edge$node2[duplicated(ARG$edge$node2) | duplicated(ARG$edge$node2,
    fromLast=TRUE)]
  recomb_nodes <- unique(recomb_nodes)
  recomb_parents <- ARG$edge$node1[ARG$edge$node2 %in% recomb_nodes]

  # add node coordinate to layout matrix
  for (i in (ARG$n+1):n) {
    target_node <- i
    if (target_node %in% recomb_parents) {
      # if the target node is from recombination
      if (is.na(layout_coord[i, 1])) {
        base_node <- ARG$edge$node2[which(target_node == ARG$edge$node1)]
        target_node <- ARG$edge$node1[ARG$edge$node2 %in% base_node] # get all two nodes
        layout_coord[target_node, 1] <- c(-1, 1) + layout_coord[base_node, 1]
        layout_coord[target_node, 2] <- current_level
        current_level <- current_level + 1
      }
    } else {
      # if the target node is from coalescent
      base_node <- ARG$edge$node2[which(target_node == ARG$edge$node1)]
      layout_coord[target_node, 1] <- mean(layout_coord[base_node, 1])
      layout_coord[target_node, 2] <- current_level
      current_level <- current_level + 1
    }
  }
  return(layout_coord)
}
