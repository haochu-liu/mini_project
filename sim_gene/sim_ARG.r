sim_ARG <- function(n, rho) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  }

  k = n
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0
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
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t <- c(t, event_time)
    t_sum <- t_sum + event_time
    # sample whether the event is a coalescent
    p_coale <- rbinom(n=1, size=1, prob=(k-1)/(k-1+rho))
    if (p_coale == 1) {
      # coalescent event
      leaf_node <- sample(pool, size=2, replace=FALSE)
      leaf1_index <- which(node$index == leaf_node[1])
      leaf2_index <- which(node$index == leaf_node[2])
      # append edges
      append_edge <- tibble(
        node1 = rep(next_node, 2),
        node2 = leaf_node,
        length = c(t_sum-node$height[leaf1_index],
                   t_sum-node$height[leaf2_index]),
        material = list(node$material[[leaf1_index]],
                        node$material[[leaf2_index]])
      )
      edge <- bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble(
        index = next_node,
        height = t_sum,
        material = iv_set_intersect(append_edge$material[[1]],
                                    append_edge$material[[2]])
      )
      node <- bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node)
      next_node <- next_node + 1L
      k <- k - 1
    } else {
      # recombination event
      u <- runif(1, min=0, max=1)
      leaf_node <- sample(pool, size=1, replace=FALSE)
      leaf_index <- which(node$index == leaf_node)
      # append edges
      append_edge <- tibble(
        node1 = c(next_node, next_node + 1L),
        node2 = rep(leaf_node, 2),
        length = c(t_sum-node$height[leaf_index],
                   t_sum-node$height[leaf_index]),
        material = node$material[[leaf_index]]
      )
      edge <- bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble(
        index = c(next_node, next_node + 1L),
        height = t_sum,
        material = list(iv_set_intersect(iv(0, u), node$material[[leaf_index]]),
                        iv_set_intersect(iv(u, 1), node$material[[leaf_index]]))
      )
      node <- bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node, next_node+1L)
      next_node <- next_node + 2L
      k <- k + 1
    }
  }
  ARG = list(edge=edge, node=node, waiting_time=t, sum_time=t_sum, n=n, rho=rho)
  class(ARG) <- "simARG"
  return(ARG)
}
