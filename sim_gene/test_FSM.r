library(igraph)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/sim_FSM_mutation.r")
source("sim_gene/local_tree.r")
source("sim_gene/localtree_to_phylo.r")
source("sim_gene/ARG_igraph.r")


set.seed(11)
ARG <- sim_FSM_ARG(5, 1, 10)
ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
g <- graph_from_edgelist(ARG_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(ARG)
plot.igraph(g)
plot.igraph(g, layout=layout_coord)
ARG_mutation <- sim_FSM_mutation(ARG, 2)

local_tree1 <- local_tree(ARG_mutation, 2)
local_tree1_matrix <- as.matrix(local_tree1$edge[, c(1, 2)])
local_g1 <- graph_from_edgelist(local_tree1_matrix, directed = FALSE)
local_g1 <- delete_vertices(local_g1, V(local_g1)[degree(local_g1) == 0])
plot(local_g1)

local_tree2 <- local_tree(ARG_mutation, 9)
local_tree2_matrix <- as.matrix(local_tree2$edge[, c(1, 2)])
local_g2 <- graph_from_edgelist(local_tree2_matrix, directed = FALSE)
local_g2 <- delete_vertices(local_g2, V(local_g2)[degree(local_g2) == 0])
plot(local_g2)

library(ape)

phylo_tree1 <- localtree_to_phylo(local_tree1, label=TRUE)
plot(phylo_tree1)

phylo_tree2 <- localtree_to_phylo(local_tree2, label=TRUE)
plot(phylo_tree2)
