#' Input: ARG (sim_ISM_ARG), mutation rate
#' Add mutation uniformly to edges
#' Store in a dataframe: edge index, position(time), mutation
#' Output: add mutation dataframe and a new column in node dataframe of genotype
sim_mutation <- function(ARG, theta) {
  if (!inherits(ARG, "sim_ISM_ARG")) {
    stop("Object must be of class 'sim_ISM_ARG'")
  }

  l <- sum(ARG$edge$length)
  n <- rpois(1, theta*l/2) # num of mutations | l ~ Poisson(theta*l/2)
  new_ARG <- ARG
  new_ARG$mutation <- tibble(
    edge_index = NA,
    pos = NA,
    site = NA
  )
  new_ARG$node$gene <- list(numeric())

  # if there is no mutation
  if (n == 0) {
    new_ARG$node$gene_str <- "[]"
    return(new_ARG)
  }

  # if there are mutations
  mutate_edge <- sample(1:nrow(ARG$edge), n,
                        replace=TRUE, prob=ARG$edge$length)
  mutate_site <- runif(n)
  keep_mutation <- c()
  for (i in 1:n) {
    edge_iv <- ARG$edge$material[[mutate_edge[i]]]
    site_iv <- iv(mutate_site[i], mutate_site[i]+.Machine$double.eps)
    if (iv_count_overlaps(site_iv, edge_iv)) {
      keep_mutation <- c(keep_mutation, i)
    }
  }
  mutate_edge <- mutate_edge[keep_mutation]
  mutate_site <- mutate_site[keep_mutation]
  if (!length(keep_mutation)) {
    new_ARG$node$gene_str <- "[]"
    return(new_ARG)
  }
  # dataframe mutations to store information of mutations
  new_ARG$mutation <- tibble(
    edge_index = mutate_edge,
    pos = rep(NA, length(mutate_edge)),
    site = mutate_site
  )
  for (i in 1:length(mutate_edge)) {
    new_ARG$mutation$pos[i] <- runif(1, max=ARG$edge$length[mutate_edge[i]])
  }

  # simulate the mutations at every node
  for (i in nrow(ARG$edge):1) {
    edge_mutation <- new_ARG$mutation$site[new_ARG$mutation$edge_index==i]
    parent_seq <- new_ARG$node$gene[[ARG$edge$node1[i]]]
    new_ARG$node$gene[[ARG$edge$node2[i]]] <- unique(sort(c(parent_seq,
      edge_mutation, new_ARG$node$gene[[ARG$edge$node2[i]]])))
  }
  # convert to string
  new_ARG$node$gene_str <- NA
  for (i in 1:nrow(new_ARG$node)) {
    new_ARG$node$gene_str[i] <- paste0("[", paste(round(new_ARG$node$gene[[i]], 3),
                                        collapse = ", "), "]")
  }
  return(new_ARG)
}
