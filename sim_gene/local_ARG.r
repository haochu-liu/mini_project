#' Input: tree (simARG), site location for local ARG
#' Select edges for local tree graph
#' Output: local ARG
local_ARG <- function(tree, location) {
  if (!inherits(tree, "simARG")) {
    stop("Object must be of class 'simARG'")
  }

  interval <- iv(location, location+.Machine$double.eps)
  keep_edge <- c()
  for (i in 1:nrow(tree$edge)) {
    if (iv_count_overlaps(interval, tree$edge$material[[i]])) {
      keep_edge <- c(keep_edge, i)
    }
  }
  local_tree <- tree
  local_tree$edge <- tree$edge[keep_edge, ]

  class(local_tree) <- "localARG"
  return(local_tree)
}
