#' Input: ARG (sim_FSM_ARG), mutation rate
#' Add mutation uniformly to edges
#' Store in a dataframe: edge index, position(time), mutation
#' Output: add mutation dataframe and a new column in node dataframe of genotype
sim_FSM_mutation <- function(ARG, theta) {
  if (!inherits(ARG, "sim_FSM_ARG")) {
    stop("Object must be of class 'sim_FSM_ARG'")
  }

  l <- sum(ARG$edge[, 3])
  n <- rpois(1, theta*l/2) # num of mutations | l ~ Poisson(theta*l/2)
  ARG$mutation <- matrix(NA, nrow=n, ncol=3)
  colnames(ARG$mutation) <- c("edge_index", "pos", "site")
  ARG$node <- data.frame(gene_str=rep("[]", length(ARG$node_height)))
  ARG$node$gene <- list(numeric())

  # if there is no mutation
  if (n == 0) {return(ARG)}

  # if there are mutations
  mutate_edge <- sample(1:nrow(ARG$edge), n,
                        replace=TRUE, prob=ARG$edge[, 3])
  mutate_site <- sample(1:ARG$L, n, replace=TRUE)
  keep_mutation <- c()

  # ignore mutations not in the edge material
  for (i in 1:n) {
    if (ARG$edge_mat[mutate_edge[i], mutate_site[i]]) {
      keep_mutation <- c(keep_mutation, i)
    }
  }
  mutate_edge <- mutate_edge[keep_mutation]
  mutate_site <- mutate_site[keep_mutation]
  ARG$mutation <- ARG$mutation[keep_mutation, ]
  if (!length(keep_mutation)) {return(ARG)}

  # dataframe mutations to store information of mutations
  ARG$mutation[, 1] <- mutate_edge
  ARG$mutation[, 3] <- mutate_site
  for (i in 1:length(mutate_edge)) {
    ARG$mutation[i, 2] <- runif(1, max=ARG$edge[mutate_edge[i], 3])
  }

  # simulate the mutations at every node
  for (i in nrow(ARG$edge):1) {
    edge_mutation <- ARG$mutation[ARG$mutation[, 1]==i, 3]
    parent_seq <- ARG$node$gene[[ARG$edge[i, 1]]]
    ARG$node$gene[[ARG$edge[i, 2]]] <- unique(sort(c(parent_seq,
      edge_mutation, ARG$node$gene[[ARG$edge[i, 2]]])))
  }
  # convert to string
  for (i in 1:nrow(ARG$node)) {
    ARG$node$gene_str[i] <- paste0("[", paste(round(ARG$node$gene[[i]], 3),
                                   collapse = ", "), "]")
  }
  return(ARG)
}
