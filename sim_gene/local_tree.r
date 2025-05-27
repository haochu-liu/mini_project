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
  }
  return(local_tree)
}
