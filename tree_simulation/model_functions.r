coale_model <- function(sample_size) {
  if (sample_size!=as.integer(sample_size)) {
    stop("Sample size must be an integer")
  }
  coale_list <- list(n=sample_size, items=list(coale_event=list()))
  class(coale_list) <- "human_genealogy"
  return(coale_list)
}

add_recombination <- function(obj, recombination) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  }
  obj$items <- append(obj$items, list(recomb_rate=recombination,
                                      recomb_event=list()))
  return(obj)
}

add_mutation <- function(obj, mutation) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  }
  obj$items <- append(obj$items, list(mutation_rate=mutation,
                                      mutation_event=list()))
  return(obj)
}
