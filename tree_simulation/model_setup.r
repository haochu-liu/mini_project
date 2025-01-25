coale_model <- function(sample_size, pop_size) {
  if (sample_size!=as.integer(sample_size) | pop_size!=as.integer(pop_size)) {
    stop("Sample size and population size must be integers")
  }
  coale_list <- list(n=sample_size, N=pop_size,
                     items=list(coale_event=list(),
                     tree=as.list(1:sample_size)))
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
