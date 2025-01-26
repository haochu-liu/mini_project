model_simulator <- function(obj) {
  if (!inherits(obj, "human_genealogy")) {
    stop("Object must be of class 'human_genealogy'")
  }
  recomb_rate <- ifelse(is.null(obj$items$recomb_rate), 0, obj$items$recomb_rate)
  mutation_rate <- ifelse(is.null(obj$items$mutation_rate), 0, obj$items$mutation_rate)
  rho <- 2*obj$N*recomb_rate
  theta <- 2*obj$N*mutation_rate
  j <- obj$n
  tree_i <- 1
  while (j > 1) {
    exp_rate <- j * ((j-1) + rho + theta) / 2
    t <- rexp(1, rate = exp_rate)
    event <- sample(1:3, 1, prob = c((j-1), rho, theta))
    if (event==1) {
      pick_two <- sample(obj$items$tree[[tree_i]], 2)
      coale_lineage <- as.list(c(as.numeric(pick_two[[1]]),
                                 as.numeric(pick_two[[2]])))
      tree_c <- setdiff(obj$items$tree[[tree_i]], list(pick_two[[1]]))
      tree_c <- setdiff(tree_c, list(pick_two[[2]]))
      obj$items$tree[[tree_i+1]] <- append(tree_c, list(coale_lineage))
      tree_i <- tree_i + 1
      j <- j - 1
      obj$items$coale_event <- append(obj$items$coale_event,
                                      list(list(t, coale_lineage)))
    } else if (event==2) {
      recomb_position <- runif(1)
      pick_one <- sample(obj$items$tree[[tree_i]], 1)
      recomb_lineage <- as.list(c(obj$n + 1 + length(obj$items$recomb_event)))
      obj$items$tree[[tree_i+1]] <- append(obj$items$tree[[tree_i]],
                                           list(recomb_lineage))
      tree_i <- tree_i + 1
      j <- j + 1
      obj$items$recomb_event <- append(obj$items$recomb_event,
                              list(list(t, pick_one, recomb_position)))
    } else {
      mutation_position <- runif(1)
      pick_one <- sample(obj$items$tree[[tree_i]], 1)
      obj$items$mutation_event <- append(obj$items$mutation_event,
                                list(list(t, pick_one, mutation_position)))
    }
  }
  return(obj)
}
