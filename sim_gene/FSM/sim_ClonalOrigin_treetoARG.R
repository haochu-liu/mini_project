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
sim_ClonalOrigin_treetoARG <- function(n, rho, L, delta, edgemat=TRUE) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (L!=as.integer(L)) {
    stop("Number of sites must be an integer")
  }
  
  k = n
  t_sum <- 0
  
  clonal_edge <- matrix(NA, nrow=2*(n-1), ncol=3) # root and leaf nodes, length
  colnames(clonal_edge) <- c("node1", "node2", "length")
  clonal_node_height <- rep(NA, 2*n-1)            # node height to recent time
  clonal_node_height[1:n] <- 0                    # initialize first n nodes
  
  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)
  
  # clonal tree by coalescent only
  while (k > 1) {
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t_sum <- t_sum + event_time
    # coalescent event
    leaf_node <- sample(pool, size=2, replace=FALSE)
    
    # append edges
    clonal_edge[c(edge_index, edge_index+1), 1] <- node_index
    clonal_edge[c(edge_index, edge_index+1), 2] <- leaf_node
    clonal_edge[c(edge_index, edge_index+1), 3] <- t_sum-clonal_node_height[leaf_node]
    
    # append root node
    clonal_node_height[node_index] <- t_sum
    
    # updates for iteration
    pool <- c(setdiff(pool, leaf_node), node_index)
    edge_index <- edge_index + 2L
    node_index <- node_index + 1L
    k <- k - 1
  }
  
  # number of recombination edges
  l <- sum(clonal_edge[, 3])
  n_recomb <- rpois(1, rho*l/2) # num of recombs | l ~ Poisson(rho*l/2)
  
  
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
