coale_model <- function(sample_size) {
  coale_list <- list(n=sample_size, coale_event=list())
  class(coale_list) <- "human_genealogy"
  return(coale_list)
}

add_recombination <- function(obj, recombination) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  } else {
    coale_list <- append(append(obj, list(recomb=list(recomb_rate=recombination,
                                                      recomb_event=list()))))
    return(obj)
  }
}
