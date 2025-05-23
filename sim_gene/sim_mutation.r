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
  new_tree$node$gene <- list(numeric())

  # if there is no mutation
  if (n == 0) {
    new_tree$node$gene_str <- "[]"
    return(new_tree)
  }

  # if there are mutations
  mutate_edge <- sample(1:nrow(tree$edge), n,
                        replace=TRUE, prob=tree$edge$length)
  mutate_site <- runif(n)
  keep_mutation <- c()
  for (i in 1:n) {
    edge_iv <- tree$edge$material[[mutate_edge[i]]]
    site_iv <- iv(mutate_site[i], mutate_site[i]+.Machine$double.eps)
    if (iv_count_overlaps(site_iv, edge_iv)) {
      keep_mutation <- c(keep_mutation, i)
    }
  }
  mutate_edge <- mutate_edge[keep_mutation]
  mutate_site <- mutate_site[keep_mutation]
  if (!length(keep_mutation)) {
    new_tree$node$gene_str <- "[]"
    return(new_tree)
  }
  # dataframe mutations to store information of mutations
  new_tree$mutation <- tibble(
    edge_index = mutate_edge,
    pos = rep(NA, length(mutate_edge)),
    site = mutate_site
  )
  for (i in 1:length(mutate_edge)) {
    new_tree$mutation$pos[i] <- runif(1, max=tree$edge$length[mutate_edge[i]])
  }

  # simulate the mutations at every node
  for (i in nrow(tree$edge):1) {
    edge_mutation <- new_tree$mutation$site[new_tree$mutation$edge_index==i]
    parent_seq <- new_tree$node$gene[[which(tree$node$index==tree$edge$node1[i])]]
    new_tree$node$gene[[which(tree$node$index==tree$edge$node2[i])]] <- unique(sort(c(parent_seq,
      edge_mutation, new_tree$node$gene[[which(tree$node$index==tree$edge$node2[i])]])))
  }
  # convert to string
  new_tree$node$gene_str <- NA
  for (i in 1:nrow(new_tree$node)) {
    new_tree$node$gene_str[i] <- paste0("[", paste(round(new_tree$node$gene[[i]], 3),
                                        collapse = ", "), "]")
  }
  return(new_tree)
}
