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

bacteria_recomb_simulator <- function(obj) {
  if (!inherits(obj, "bacteria_genealogy")) {
    stop("Object must be of class 'bacteria_genealogy'")
  }
  T <- c()
  for (i in 1:length(obj$items$coale_event)) {
    T <- c(T, obj$items$coale_event[[i]][[1]])
  }
  T <- cumsum(T)
  rho <- 2 * obj$N * obj$items$mutation_rate
  R <- rpois(1, rho * tail(T, 1) / 2) # R ~ Poi(rho * T / 2)
  b <- runif(R, min=0, max=tail(T, 1)) # b_i|R ~ U(0, T)
  t <- rexp(R, rate=1) # t_i|b_i ~ Exp(1)
  a <- b + t # a_i = b_i + t_i
}
