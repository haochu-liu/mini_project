add_bacteria_recomb <- function(obj, recombination) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  } else if (!is.null(obj$items$recomb_rate)) {
    stop("Human recombination exists")
  }
  obj$items <- append(obj$items, list(recomb_rate=recombination,
                                      recomb_event=list()))
  class(obj) <- "bacteria_genealogy"
  return(obj)
}


