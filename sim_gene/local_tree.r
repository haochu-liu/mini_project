source("sim_gene/localtree_tools/localtree_height.r")
source("sim_gene/localtree_tools/localtree_traj.r")


#' Input: ARG (sim_ISM_ARG or sim_FSM_ARG), site location for local tree
#' Call function local_tree_ISM or local_tree_FSM
#' Output: local tree
local_tree <- function(ARG, location) {
  if (!inherits(ARG, "sim_ISM_ARG") & !inherits(ARG, "sim_FSM_ARG")) {
    stop("Object must be of class 'sim_ISM_ARG' or 'sim_FSM_ARG'")
  }

  if (inherits(ARG, "sim_ISM_ARG")) {
    source("sim_gene/ISM/local_tree_ISM.r")
    local_tree <- local_tree_ISM(ARG, location)
  } else {
    source("sim_gene/FSM/local_tree_FSM.r")
    local_tree <- local_tree_FSM(ARG, location)
  }

  local_tree$waiting_time <- NULL
  local_tree$k <- NULL
  local_tree$sum_time <- tree_height(local_tree)
  local_tree$mutation[local_tree$mutation[, 1] %in% local_tree$edge_index, ]

  return(local_tree)
}
