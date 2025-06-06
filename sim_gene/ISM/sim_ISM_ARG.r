#' Input: number of leaf alleles, recombination parameter,
#' bacteria recombination or not,
#' if yes, delta is mean of the length of recombinant segment
#' Create a full ARG with coalescent and recombination
#' Edge dataframe: root node, leaf node, edge length, edge material interval
#' Node dataframe: node index, node height, node material interval
#' Output: edge dataframe, node dataframe, waiting time for each event,
#' total time, number of lineages at each event time, number of leaf alleles,
#' recombination parameter, bacteria recombination or not, and parameter delta
sim_ISM_ARG <- function(n, rho, bacteria=FALSE, delta=NULL) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (bacteria & is.null(delta)) {
    stop("Must provide parameter for the length of recombinant segment")
  }

  k = n
  k_vector <- c(k)
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0
  edge <- tibble(
    node1 = integer(), # root node
    node2 = integer(), # leaf node
    length = numeric(), # edge length
    material = list() # edge material interval
  )
  node <- tibble(
    height = rep(0, n), # node height to recent time
    material = list(iv(0, 1)) # node material interval
  )
  pool <- as.integer(1:n)
  next_node <- as.integer(n+1)

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
      # append edges
      append_edge <- tibble(
        node1 = rep(next_node, 2),
        node2 = leaf_node,
        length = c(t_sum-node$height[leaf_node[1]],
                   t_sum-node$height[leaf_node[2]]),
        material = list(node$material[[leaf_node[1]]],
                        node$material[[leaf_node[2]]])
      )
      edge <- bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble(
        height = t_sum,
        material = list(iv_set_union(append_edge$material[[1]],
                                     append_edge$material[[2]]))
      )
      node <- bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node)
      next_node <- next_node + 1L
      k <- k - 1
    } else {
      # recombination event
      leaf_node <- sample(pool, size=1, replace=FALSE)
      if (bacteria) {
        x <- runif(1, min=0, max=1)
        y <- min(1, x + rexp(1, rate=1/delta))
        iv1 <- iv(x, y)
        iv2 <- iv_set_complement(iv1, lower = 0, upper = 1)
        df_mat <- list(iv_set_intersect(iv1, node$material[[leaf_node]]),
                       iv_set_intersect(iv2, node$material[[leaf_node]]))
      } else {
        u <- runif(1, min=0, max=1)
        df_mat <- list(iv_set_intersect(iv(0, u), node$material[[leaf_node]]),
                       iv_set_intersect(iv(u, 1), node$material[[leaf_node]]))
      }
      # append edges
      append_edge <- tibble(
        node1 = c(next_node, next_node + 1L),
        node2 = rep(leaf_node, 2),
        length = c(t_sum-node$height[leaf_node],
                   t_sum-node$height[leaf_node]),
        material = df_mat
      )
      edge <- bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble(
        height = t_sum,
        material = df_mat
      )
      node <- bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node, next_node+1L)
      next_node <- next_node + 2L
      k <- k + 1
    }
    k_vector <- c(k_vector, k)
  }
  ARG = list(edge=edge, node=node, waiting_time=t, sum_time=t_sum, k=k_vector,
             n=n, rho=rho, bacteria=bacteria, delta=delta)
  class(ARG) <- "sim_ISM_ARG"
  return(ARG)
}
