add_bacteria_recomb <- function(obj, recombination, delta) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  } else if (!is.null(obj$items$recomb_rate)) {
    stop("Human recombination exists")
  }
  obj$items <- append(obj$items, list(recomb_rate=recombination,
                                      recomb_length=delta,
                                      recomb_event=list()))
  class(obj) <- "bacteria_genealogy"
  return(obj)
}

bacteria_recomb_simulator <- function(obj) {
  if (!inherits(obj, "bacteria_genealogy")) {
    stop("Object must be of class 'bacteria_genealogy'")
  }
  T <- c()
  L <- c()
  for (i in 1:length(obj$items$coale_event)) {
    L <- c(L, obj$items$coale_event[[i]][[1]])
    T <- c(T, l * length(obj$items$tree[[i]]))
  }
  T_sum <- sum(T)
  L_sum <- sum(L)
  rho <- 2 * obj$N * obj$items$recomb_rate

  R <- rpois(1, rho * T_sum / 2) # R ~ Poi(rho * T / 2)
  if (R == 0) {return(obj)}
  for (i in 1:R) {
    # b <- runif(1, min=0, max=tail(T, 1)) # b_i|R ~ U(0, T)
    t <- rexp(1, rate=1) # t_i|b_i ~ Exp(1)
    # a <- b + t # a_i = b_i + t_i
    while (TRUE) {
      l <- rexp(1, rate=1/obj$items$recomb_length) # l ~ Exp(1/delta)
      if (l <= 1) {break}
    }
    x <- runif(1, min=0, max=1-l) # x ~ U(0, 1-l)
    y <- x + l
  }

  
  
}
