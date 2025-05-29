library(igraph)
source("sim_gene/FSM/sim_FSM_ARG.r")


set.seed(11)
ARG <- sim_ISM_ARG(5, 1, 10)
ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
g <- graph_from_edgelist(ARG_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(ARG)
plot.igraph(g)
plot.igraph(g, layout=layout_coord)
ARG_mutation <- sim_ISM_mutation(ARG, 2)