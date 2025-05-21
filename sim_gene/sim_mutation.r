#' Input: tree (simARG), mutation rate
#' Add mutation uniformly to edges
#' Store in a dataframe: edge index, position(time), mutation
#' Output: add mutation dataframe and a new column in node dataframe of genotype 
sim_mutation <- function(tree, theta) {
  if (!inherits(tree, "simARG")) {
    stop("Object must be of class 'simARG'")
  }

  l <- sum(tree$edge$length)
  n <- rpois(1, theta*l/2) # num of mutations | l ~ Poisson(theta*l/2)

  # if there is no mutation
  if (n == 0) {
    mutations <- data.frame(edge_index = NA,
        pos = NA,
        site = NA)
    node_seq <- data.frame(node = NA)
    return(list(mutations=mutations, node_seq=node_seq))
  }

  # if there are mutations
  mutate_edges <- sample(1:nrow(tree$edge), n,
                         replace=TRUE, prob=tree$edge.length)
  mutate_site <- runif(n)
  # dataframe mutations to store information of every mutation
  mutations <- data.frame(edge_index = mutate_edges,
                          pos = rep(NA, n),
                          site = mutate_site)
  for (i in 1:n) {
      mutations[i, 2] <- runif(1, max=tree$edge.length[mutate_edges[i]])
  }
  
  # simulate the mutations before every node
  node_seq <- data.frame(node=rep(NA, 2*tree$n-1))
  node_seq$node <- as.list(node_seq$node)
      node_seq$node[[tree$n+1]] <- numeric(0)
      for (i in nrow(tree$edge):1) {
          edge_mutation <- mutations$site[mutations$edge_index==i]
          parent_seq <- node_seq$node[[tree$edge[i, 1]]]
          node_seq$node[[tree$edge[i, 2]]] <- sort(c(parent_seq, edge_mutation))
      }
  return(list(mutations=mutations, node_seq=node_seq))
}
