# library(ape)


simtree_to_phylo <- function(tree, label=NULL) {
  if (!inherits(tree, "simtree")) {
    stop("Object must be of class 'simtree'")
  }

  if (is.null(label)) {label <- paste0("t", 1:tree$n)}
  # convert the simtree object to phylo object for ape::plot.phylo
  tree_phylo <- list(edge=tree$edge, edge.length=tree$edge.length,
                     tip.label=label,
                     Nnode=as.integer(tree$n - 1))
  class(tree_phylo) <- "phylo"
  # tree_phylo <- ape::reorder.phylo(tree_phylo, order="cladewise")

  return(tree_phylo)
}
