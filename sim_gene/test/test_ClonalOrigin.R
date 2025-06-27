library(igraph)
source("sim_gene/FSM/sim_ClonalOrigin_ARG.r")
source("sim_gene/FSM/sim_FSM_mutation.r")
source("sim_gene/local_tree.r")
source("sim_gene/ARG_igraph.r")


set.seed(11)
ARG <- sim_ClonalOrigin_ARG(5, 1, 10, 5, optimise_recomb = FALSE)
ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
g <- graph_from_edgelist(ARG_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(ARG)
plot.igraph(g)
plot.igraph(g, layout=layout_coord)
ARG_mutation <- sim_FSM_mutation(ARG, 2)

local_tree1 <- local_tree(ARG_mutation, 2)




