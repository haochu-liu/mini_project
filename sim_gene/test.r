library(ivs)
library(tibble)
library(dplyr)
library(igraph)
source("sim_gene/sim_ARG.r")
source("sim_gene/sim_mutation.r")
source("sim_gene/local_ARG.r")


set.seed(10)
tree <- sim_ARG(5, 1)
tree_matrix <- as.matrix(tree$edge[, c(1, 2)])
g <- graph_from_edgelist(tree_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
plot(g)

tree_mutation <- sim_mutation(tree, 2)

local_tree <- local_ARG(tree, 0.5)
local_tree_matrix <- as.matrix(local_tree$edge[, c(1, 2)])
local_g <- graph_from_edgelist(local_tree_matrix, directed = FALSE)
local_g <- delete_vertices(local_g, V(local_g)[degree(local_g) == 0])
plot(local_g)
