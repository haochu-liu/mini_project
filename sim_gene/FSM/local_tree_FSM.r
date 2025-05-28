#' Input: ARG (sim_FSM_ARG), site location for local tree
#' Select edges for local tree graph
#' Output: local tree
local_tree_FSM <- function(ARG, location) {
  keep_edge <- c()
  for (i in 1:nrow(ARG$edge_mat)) {
    if (ARG$edge_mat[i, location]) {
      keep_edge <- c(keep_edge, i)
    }
  }
  ARG$edge_matrix <- ARG$edge_matrix[keep_edge, ]
  ARG$edge_length <- ARG$edge_length[keep_edge]
  ARG$edge_mat <- ARG$edge_mat[keep_edge, ]

  class(ARG) <- "localtree"
  return(ARG)
}
