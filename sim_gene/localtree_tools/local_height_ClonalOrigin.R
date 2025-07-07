local_height_ClonalOrigin <- function(ARG, location) {
  if (!inherits(ARG, "ClonalOrigin")) {
    cli::cli_abort("Object must be of class 'ClonalOrigin'!")
  }
  
  if (is.null(ARG$recomb_edge)) {
    return(ARG$sum_time)
  }
  edge_index <- which(ARG$recomb_edge[, 3] == -1 &
                      ARG$recomb_edge[, 5] <= location &
                      ARG$recomb_edge[, 6] >= location)
  if (!length(edge_index)) {
    return(ARG$sum_time)
  } else {
    return(max(ARG$recomb_edge[edge_index, 4]))
  }
}
