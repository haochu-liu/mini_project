library(ivs)
library(tibble)
library(dplyr)
library(igraph)
source("sim_gene/sim_ARG.r")
source("sim_gene/sim_mutation.r")


set.seed(100)
tree <- sim_ARG(5, 1)
tree_matrix <- as.matrix(tree$edge[, c(1, 2)])
g <- graph_from_edgelist(tree_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
plot(g)

tree_mutation <- sim_mutation(tree, 1)
