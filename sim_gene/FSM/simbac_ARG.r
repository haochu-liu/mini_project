source("sim_gene/effective_R.r")


#' Input: number of leaf lineages, recombination parameter, number of sites,
#' delta is mean of the length of recombinant segment,
#' initial maximal node size (default = 1000),
#' optimise recombination edges or not
#' Create a full ARG with coalescent and recombination
#' Edge dataframe: root node, leaf node, edge length, edge material interval
#' Node dataframe: node index, node height, node material interval
#' Output: edge dataframe, node dataframe, waiting time for each event,
#' total time, number of lineages at each event time, number of leaf alleles,
#' recombination parameter, bacteria recombination or not, and parameter delta
simbac_ARG <- function(n, rho, L, delta, node_max=1000, output_eff_R=FALSE,
                       optimise_site=FALSE, edgemat=TRUE) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (L!=as.integer(L)) {
    stop("Number of sites must be an integer")
  } else if (n >= node_max) {
    stop("Maximal node size must greater than the number of leaf lineages")
  }

  k = n
  k_vector <- c(k)
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0

  edge_matrix <- matrix(NA, nrow=node_max, ncol=3)    # root and leaf nodes, length
  colnames(edge_matrix) <- c("node1", "node2", "length")
  edge_mat_index <- rep(NA, node_max)                 # edge material index
  node_height <- rep(NA, node_max)                    # node height to recent time
  node_mat <- matrix(NA, nrow=node_max, ncol=L)       # node material
  node_eff_R <- matrix(NA, nrow=node_max, ncol=2)     # node effective R
  colnames(node_eff_R) <- c("effective_R", "clonal")
  node_probstart <- matrix(NA, nrow=node_max, ncol=L) # node recomb starting prob
  node_height[1:n] <- 0                               # initialize nodes height
  node_mat[1:n, ] <- 1                                # initialize nodes material
  node_eff_R[1:n, 1] <- rho * (1-(1-1/delta)^(L-1))   # initialize effective R
  node_eff_R[1:n, 2] <- TRUE
  node_probstart[1:n, ] <- matrix(cumsum(rep(1/L, L)), nrow=n, ncol=L, byrow=TRUE)

  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)
  if (optimise_site) {
    include_site <- rep(TRUE, L)
  } else {
    include_site <- NULL
  }

  while (k > 1) {
    # sample a new event time
    rho_eff <- sum(node_eff_R[pool, 1])
    event_time <- rexp(1, rate=(k*(k-1)+rho_eff)/2)
    t <- c(t, event_time)
    t_sum <- t_sum + event_time
    # sample whether the event is a coalescent
    p_coale <- rbinom(n=1, size=1, prob=k*(k-1)/(k*(k-1)+rho_eff))
    if (p_coale == 1) {
      # coalescent event
      leaf_node <- sample(pool, size=2, replace=FALSE)

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- node_index
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
      edge_mat_index[c(edge_index, edge_index+1)] <- leaf_node

      # append root node
      node_height[node_index] <- t_sum
      node_mat[node_index, ] <- as.integer(node_mat[leaf_node[1], ] |
                                           node_mat[leaf_node[2], ])
      
      # update effective R
      node_eff_R[node_index, 2] <- any(as.logical(node_eff_R[leaf_node, 2]))
      if (optimise_site) {
        list_eff_R <- effective_R(node_mat[node_index, ], delta, rho,
                                  node_eff_R[node_index, 2], include_site=include_site)
      } else {
        list_eff_R <- effective_R(node_mat[node_index, ], delta, rho,
                                  node_eff_R[node_index, 2])
      }
      node_eff_R[node_index, 1] <- list_eff_R$R_eff
      node_probstart[node_index, ] <- list_eff_R$probstartcum
      
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), node_index)
      edge_index <- edge_index + 2L
      node_index <- node_index + 1L
      k <- k - 1

      if (optimise_site & any(include_site)) {
        # check if contain exclude sites
        include_site_index <- which(include_site)
        for (i in 1:sum(include_site)) {
          if (sum(node_mat[pool, include_site_index[i]]) == 1) {
            include_site[include_site_index[i]] <- FALSE
            for (j in pool) {
              list_eff_R <- effective_R(node_mat[j, ], delta, rho,
                                        node_eff_R[j, 2],
                                        include_site=include_site)
              node_eff_R[j, 1] <- list_eff_R$R_eff
              node_probstart[j, ] <- list_eff_R$probstartcum
            }
          }
        }
      }

    } else {
      # recombination event
      leaf_node <- sample(pool, size=1, replace=FALSE, prob=node_eff_R[pool, 1])

      # sample recombination segment
      x <- which(runif(1) < node_probstart[leaf_node, ])[1]
      # pmf for conditional geom distribution
      v_s <- which(node_mat[leaf_node, ] &
        (node_mat[leaf_node, ] != c(0, node_mat[leaf_node, ][1:(L-1)])))
      site_index <- which(v_s %in% x)
      if (length(site_index)) {
        v_e <- which(node_mat[leaf_node, ] &
          (node_mat[leaf_node, ] != c(node_mat[leaf_node, ][2:L], 0)))
        v_e <- c(v_e[length(v_e)] - L, v_e)
        s <- min(L - (v_s[site_index] - v_e[site_index]), L - x + 1)
      } else {s <- L - x + 1}
      r_pmf <- (1 - 1/delta)^(1:s - 1) / (delta * (1 - (1 - 1/delta)^s))
      y <- which(runif(1) < cumsum(r_pmf))[1] + x - 1

      # repeat {
      #   x <- which(runif(1) < node_probstart[leaf_node, ])[1]
      #   y <- min(x + rgeom(1, 1/delta), L)
      #   if ((!(sum(node_mat[leaf_node, x:y])==0 |
      #          sum(node_mat[leaf_node, -(x:y)])==0)) &
      #       (!node_eff_R[leaf_node, 2])) {
      #     break
      #   } else if ((!(sum(node_mat[leaf_node, x:y])==0)) &
      #              node_eff_R[leaf_node, 2]) {
      #     break
      #   }
      # }
      
      edge_mat_index[c(edge_index, edge_index+1)] <- c(node_index, node_index+1)

      node_mat[c(node_index, node_index+1), ] <- 0
      node_mat[node_index, x:y] <- node_mat[leaf_node, x:y]
      node_mat[node_index+1, -(x:y)] <- node_mat[leaf_node, -(x:y)]
      
      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- c(node_index, node_index+1L)
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]

      # append root node
      node_height[c(node_index, node_index+1)] <- t_sum

      # update effective R
      node_eff_R[node_index, 2] <- FALSE
      if (optimise_site) {
        list_eff_R <- effective_R(node_mat[node_index, ], delta, rho,
                                  node_eff_R[node_index, 2], include_site=include_site)
      } else {
        list_eff_R <- effective_R(node_mat[node_index, ], delta, rho,
                                  node_eff_R[node_index, 2])
      }
      node_eff_R[node_index, 1] <- list_eff_R$R_eff
      node_probstart[node_index, ] <- list_eff_R$probstartcum

      node_eff_R[node_index+1, 2] <- node_eff_R[leaf_node, 2]
      if (optimise_site) {
        list_eff_R <- effective_R(node_mat[node_index+1, ], delta, rho,
                                  node_eff_R[node_index+1, 2], include_site=include_site)
      } else {
        list_eff_R <- effective_R(node_mat[node_index+1, ], delta, rho,
                                  node_eff_R[node_index+1, 2])
      }
      node_eff_R[node_index+1, 1] <- list_eff_R$R_eff
      node_probstart[node_index+1, ] <- list_eff_R$probstartcum

      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), node_index, node_index+1L)
      edge_index <- edge_index + 2L
      node_index <- node_index + 2L
      k <- k + 1
    }
    k_vector <- c(k_vector, k)
    if (max(edge_index, node_index) >= node_max - 1) {
      # add empty rows or elements if more edges than expected
      edge_matrix <- rbind(edge_matrix, matrix(NA, nrow=node_max, ncol=3))
      edge_mat_index <- c(edge_mat_index, rep(NA, node_max))
      node_height <- c(node_height, rep(NA, node_max))
      node_mat <- rbind(node_mat, matrix(NA, nrow=node_max, ncol=L))
      node_eff_R <- rbind(node_eff_R, matrix(NA, nrow=node_max, ncol=2))
      node_probstart <- rbind(node_probstart, matrix(NA, nrow=node_max, ncol=L))
      node_max <- 2 * node_max
    }
  }

  if (output_eff_R & edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta, node_eff_R=node_eff_R[1:(node_index-1), ],
               node_probstart=node_probstart[1:(node_index-1), ])
  } else if (edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta)
  } else if (output_eff_R) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta, node_eff_R=node_eff_R[1:(node_index-1), ],
               node_probstart=node_probstart[1:(node_index-1), ])
  } else {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta)
  }
  class(ARG) <- "sim_FSM_ARG"
  return(ARG)
}
