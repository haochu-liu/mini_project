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
  T <- c() # T[i] = sum of branch lengths when n+1-i lineages
  L <- c() # L[i] = time / length when n+1-i lineages
  for (i in 1:length(obj$items$coale_event)) {
    l <- obj$items$coale_event[[i]][[1]]
    L <- c(L, l)
    T <- c(T, l * length(obj$items$tree[[i]]))
  }
  T_sum <- sum(T)
  L_sum <- sum(L)
  rho <- 2 * obj$N * obj$items$recomb_rate

  R <- rpois(1, rho * T_sum / 2) # R ~ Poi(rho * T / 2)
  if (R == 0) {return(obj)}
  for (i in 1:R) {
    # b_i|R ~ U(0, T)
    b_i_t <- sample(1:length(T), 1, prob = T)
    b_i_lineage <- sample(1:length(obj$items$tree[[b_i_t]]), 1)
    b_i_position <- runif(1, min=0, max=L[b_i_t])
    t_i <- rexp(1, rate=1) # t_i|b_i ~ Exp(1)
    # a_i = b_i + t_i
    target_length <- b_i_position + t_i
    accumulated_length <- 0
    for (j in b_i_t:length(T)) {
      accumulated_length <- accumulated_length + L[j]
      if (target_length <= accumulated_length) {
        a_i_t <- j
        a_i_position <- target_length - accumulated_length + L[j]
        a_i_lineage <- sample(1:length(obj$items$tree[[a_i_t]]), 1)
        break
      }
    }
    if (target_length > accumulated_length & j == length(T)) {
      a_i_t <- length(T) + 1
      a_i_position <- target_length - accumulated_length
      a_i_lineage <- 1
    }
    
    while (TRUE) {
      l <- rexp(1, rate=1/obj$items$recomb_length) # l ~ Exp(1/delta)
      if (l <= 1) {break}
    }
    x <- runif(1, min=0, max=1-l) # x ~ U(0, 1-l)
    y <- x + l

    recomb_list <- list(b=list(n=b_i_t,
                               lineage=b_i_lineage,
                               position=b_i_position),
                        a=list(n=a_i_t,
                               lineage=a_i_lineage,
                               position=a_i_position),
                        x=x, y=y)
    obj$items$recomb_event <- append(obj$items$recomb_event, list(recomb_list))
  }
}
