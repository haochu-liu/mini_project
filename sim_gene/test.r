library(ivs)
library(tibble)
library(dplyr)
library(igraph)
source("sim_gene/sim_ARG.r")


set.seed(100)
tree1 <- sim_ARG(5, 1)
tree1_matrix <- as.matrix(tree1$edge[, c(1, 2)])
g <- graph_from_edgelist(tree1_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
plot(g)
