#' Input: number of leaf lineages, recombination parameter, number of sites,
#' delta is mean of the length of recombinant segment,
#' initial maximal node size (default = 1000),
#' optimise recombination edges or not
#' Create a full ARG using ClonalOrigin
#' Edge dataframe: root node, leaf node, edge length, edge material interval
#' Node dataframe: node index, node height, node material interval
#' Output: edge dataframe, node dataframe, waiting time for each event,
#' total time, number of lineages at each event time, number of leaf alleles,
#' recombination parameter, and parameter delta
sim_ClonalOrigin_treetoARG <- function(n, rho, L, delta, node_max=1000, edgemat=TRUE) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (L!=as.integer(L)) {
    stop("Number of sites must be an integer")
  } else if (n >= node_max) {
    stop("Maximal node size must greater than the number of leaf lineages")
  }
  
  k = n
  t_sum <- 0
  
  edge_matrix <- matrix(NA, nrow=node_max, ncol=3) # root and leaf nodes, length
  colnames(edge_matrix) <- c("node1", "node2", "length")
  node_height <- rep(NA, node_max)                 # node height to recent time
  node_height[1:n] <- 0                            # initialize first n nodes
  
  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)
  
  # Clonal Tree by coalescent only
  while (k > 1) {
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t_sum <- t_sum + event_time
    # coalescent event
    leaf_node <- sample(pool, size=2, replace=FALSE)
    
    # append edges
    edge_matrix[c(edge_index, edge_index+1), 1] <- node_index
    edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
    edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
    
    # append root node
    node_height[node_index] <- t_sum
    
    # updates for iteration
    pool <- c(setdiff(pool, leaf_node), node_index)
    edge_index <- edge_index + 2L
    node_index <- node_index + 1L
    k <- k - 1

    if (max(edge_index, node_index) >= node_max - 1) {
      # add empty rows or elements if more edges than expected
      edge_matrix <- rbind(edge_matrix, matrix(NA, nrow=node_max, ncol=3))
      node_height <- c(node_height, rep(NA, node_max))
      node_max <- 2 * node_max
    }
  }
  
  if (edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  } else {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  }
  class(ARG) <- "sim_FSM_ARG"
  return(ARG)
}
