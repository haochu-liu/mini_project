#' Input: number of leaf alleles, recombination parameter, number of sites
#' bacteria recombination or not,
#' if yes, delta is mean of the length of recombinant segment,
#' initial maximal node size (default = 1000)
#' Create a full ARG with coalescent and recombination
#' Edge dataframe: root node, leaf node, edge length, edge material interval
#' Node dataframe: node index, node height, node material interval
#' Output: edge dataframe, node dataframe, waiting time for each event,
#' total time, number of lineages at each event time, number of leaf alleles,
#' recombination parameter, bacteria recombination or not, and parameter delta
sim_FSM_ARG <- function(n, rho, L, bacteria=FALSE, delta=NULL, node_max=1000) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (L!=as.integer(L)) {
    stop("Number of sites must be an integer")
  } else if (bacteria & is.null(delta)) {
    stop("Must provide parameter for the length of recombinant segment")
  }

  k = n
  k_vector <- c(k)
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0

  edge_matrix <- matrix(NA, nrow=node_max, ncol=3) # root and leaf nodes, length
  colnames(edge_matrix) <- c("node1", "node2", "length")
  edge_mat <- matrix(NA, nrow=node_max, ncol=L)    # edge material
  node_height <- rep(NA, node_max)                 # node height to recent time
  node_mat <- matrix(NA, nrow=node_max, ncol=L)    # node material
  node_height[1:n] <- 0                            # initialize first n nodes
  node_mat[1:n, ] <- 1

  # Probability of starting recombination at each site
  probstart <- rep(1, L)
  if (bacteria) {probstart[1] <- delta}
  probstart <- probstart / sum(probstart)
  probstartcum <- cumsum(probstart)

  edge_index <- 1
  node_index <- n + 1
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
      edge_matrix[c(edge_index, edge_index+1), 1] <- next_node
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
      edge_mat[c(edge_index, edge_index+1), ] <- node_mat[leaf_node, ]
      # append root node
      node_height[node_index] <- t_sum
      node_mat[node_index, ] <- as.integer(node_mat[leaf_node[1], ] |
                                           node_mat[leaf_node[2], ])
      # updates for iteration
      edge_index <- edge_index + 2
      node_index <- node_index + 1
      pool <- c(setdiff(pool, leaf_node), next_node)
      next_node <- next_node + 1L
      k <- k - 1
    } else {
      # recombination event
      leaf_node <- sample(pool, size=1, replace=FALSE)
      if (bacteria) {
        x <- which(runif(1) < probstartcum)[1]
        y <- min(x + rgeom(1, 1/delta), L)

        if (sum(node_mat[leaf_node, x:y])==0 |
            sum(node_mat[leaf_node, -(x:y)])==0) {
          next
        }

        edge_mat[c(edge_index, edge_index+1), ] <- 0
        edge_mat[edge_index, x:y] <- node_mat[leaf_node, x:y]
        edge_mat[edge_index+1, -(x:y)] <- node_mat[leaf_node, -(x:y)]

        node_mat[c(node_index, node_index+1), ] <- 0
        node_mat[node_index, x:y] <- node_mat[leaf_node, x:y]
        node_mat[node_index+1, -(x:y)] <- node_mat[leaf_node, -(x:y)]
      } else {
        x <- which(runif(1) < probstartcum)[1]

        if (sum(node_mat[leaf_node, 1:(x-1)])==0 |
            sum(node_mat[leaf_node, x:L])==0) {
          next
        }

        edge_mat[c(edge_index, edge_index+1), ] <- 0
        edge_mat[edge_index, 1:(x-1)] <- node_mat[leaf_node, 1:(x-1)]
        edge_mat[edge_index+1, x:L] <- node_mat[leaf_node, x:L]

        node_mat[c(node_index, node_index+1), ] <- 0
        node_mat[node_index, 1:(x-1)] <- node_mat[leaf_node, 1:(x-1)]
        node_mat[node_index+1, x:L] <- node_mat[leaf_node, x:L]
      }
      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- c(next_node, next_node + 1L)
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
      # append root node
      node_height[c(node_index, node_index+1)] <- t_sum
      # updates for iteration
      edge_index <- edge_index + 2
      node_index <- node_index + 2
      pool <- c(setdiff(pool, leaf_node), next_node, next_node+1L)
      next_node <- next_node + 2L
      k <- k + 1
    }
    k_vector <- c(k_vector, k)
    if (edge_index >= node_max) {
      # add empty rows or elements if more edges than expected
      edge_matrix <- rbind(edge_matrix, matrix(NA, nrow=node_max, ncol=3))
      edge_mat <- rbind(edge_mat, matrix(NA, nrow=node_max, ncol=L))
      node_height <- c(node_height, rep(NA, node_max))
      node_mat <- rbind(node_mat, matrix(NA, nrow=node_max, ncol=L))
    }
  }

  ARG = list(edge=edge_matrix[complete.cases(edge_matrix), ],
             edge_mat=edge_mat[complete.cases(edge_mat), ],
             node_height=node_height[!is.na(node_height)],
             node_mat=node_mat[complete.cases(node_mat), ],
             waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho,
             bacteria=bacteria, delta=delta)
  class(ARG) <- "sim_FSM_ARG"
  return(ARG)
}
