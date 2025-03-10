simcoale <- function(n) {
  n_edge <- 2 * (n-1) # number of edges in a coale tree
  edge <- matrix(NA, nrow=n_edge, ncol=2)
  edge.length <- rep(NA, n_edge)
  t <- rexp(n-1, rate=1) / (n:2 * (n-1):1 / 2) # vector of coalescent times
  t_sum <- cumsum(t) # cumulative sum of coalescent times
  node.length <- rep(0, 2*n-1)
  pool <- as.integer(1:n)
  next_node <- as.integer(2*n-1)
  # assign node and edges to coalescent times
  for (i in 1:(n-1)) {
    offspring_node <- sample(pool, size=2, replace=FALSE)
    index <- i*2 - 1:0
    edge[index, 2] <- offspring_node
    edge[index, 1] <- next_node
  }
}