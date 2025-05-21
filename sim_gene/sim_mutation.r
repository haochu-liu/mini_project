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
  new_tree <- tree
  new_tree$mutation <- tibble(
    edge_index = NA,
    pos = NA,
    site = NA
  )
  new_tree$node$gene <- numeric()

  # if there is no mutation
  if (n == 0) {
    new_tree$node$gene_str <- "[]"
    return(new_tree)
  }

  # if there are mutations
  mutate_edge <- sample(1:nrow(tree$edge), n,
                        replace=TRUE, prob=tree$edge$length)
  mutate_site <- runif(n)
  drop_mutation <- c()
  for (i in 1:n) {
    edge_iv <- tree$edge$material[[mutate_edge[i]]]
    site_iv <- iv(mutate_site[i], mutate_site[i]+.Machine$double.eps)
    if (!iv_count_overlaps(site_iv, edge_iv)) {
      drop_mutation <- c(drop_mutation, i)
    }
  }
  mutate_edge <- mutate_edge[-drop_mutation]
  mutate_site <- mutate_site[-drop_mutation]
  # dataframe mutations to store information of mutations
  new_tree$mutation <- tibble(edge_index = mutate_edge,
                              pos = NA,
                              site = mutate_site)
  for (i in 1:length(mutate_edge)) {
    new_tree$mutation$pos <- runif(1, max=tree$edge$length[mutate_edge[i]])
  }

  # simulate the mutations at every node
  node_seq$node <- as.list(node_seq$node)
      node_seq$node[[tree$n+1]] <- numeric(0)
      for (i in nrow(tree$edge):1) {
          edge_mutation <- mutations$site[mutations$edge_index==i]
          parent_seq <- node_seq$node[[tree$edge[i, 1]]]
          node_seq$node[[tree$edge[i, 2]]] <- sort(c(parent_seq, edge_mutation))
      }
  return(list(mutations=mutations, node_seq=node_seq))
}
