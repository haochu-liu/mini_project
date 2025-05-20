sim_coale_recomb <- function(n, rho) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  }

  k = n
  t <- vector("numeric", length = 0) # vector of event times
  edge <- tibble(
    node1 = integer(), # root node
    node2 = integer(), # leaf node
    length = numeric(), # edge length
    material = list() # edge material interval
  )
  node <- tibble(
    index = as.integer(1:n), # node number
    height = rep(0, n), # node height to recent time
    material = list(iv(0, 1))
  )
  pool <- as.integer(1:n)
  next_node <- as.integer(2*n)

  while (k > 1) {
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t <- c(t, event_time)
    p_coale <- rbinom(n=1, size=1, prob=(k-1)/(k-1+rho))
    if (p_coale == 1) {
      # coalescent event
      leaf_node <- sample(pool, size=2, replace=FALSE)
      append_edge <- matrix(NA, nrow = 2, ncol = 2)
      append_edge[, 2] <- leaf_node
      append_edge[, 1] <- next_node
      edge <- rbind(edge, append_edge)
      edge.material

      edge.length[index] <- t_sum[i] - node.length[offspring_node]
      node.length[next_node] <- t_sum[i]
      pool <- c(setdiff(pool, leaf_node), next_node)
      next_node <- next_node + 1L
    }
  }
}