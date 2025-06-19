library(igraph)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/sim_FSM_mutation.r")
source("sim_gene/local_tree.r")
source("sim_gene/localtree_tools/localtree_to_phylo.r")
source("sim_gene/ARG_igraph.r")


set.seed(11)
ARG <- sim_FSM_ARG(5, 1, 10, bacteria = TRUE, delta = 5, optimise_recomb = TRUE, clonal = TRUE)
ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
g <- graph_from_edgelist(ARG_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(ARG)
plot.igraph(g)
plot.igraph(g, layout=layout_coord)
ARG_mutation <- sim_FSM_mutation(ARG, 2)

local_tree1 <- local_tree(ARG_mutation, 2)
local_tree2 <- local_tree(ARG_mutation, 5)
local_tree3 <- local_tree(ARG_mutation, 8)

library(ape)

phylo_tree1 <- localtree_to_phylo(local_tree1, label=TRUE)
plot(phylo_tree1)

phylo_tree2 <- localtree_to_phylo(local_tree2, label=TRUE)
plot(phylo_tree2)

phylo_tree3 <- localtree_to_phylo(local_tree3, label=TRUE)
plot(phylo_tree3)
