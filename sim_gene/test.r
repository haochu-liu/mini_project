library(ivs)
library(tibble)
library(dplyr)
library(igraph)
library(ggplot2)
library(purrr)
library(deSolve)
source("sim_gene/sim_ARG.r")
source("sim_gene/sim_mutation.r")
source("sim_gene/local_ARG.r")
source("sim_gene/localARG_to_phylo.r")


set.seed(11)
tree <- sim_ARG(5, 1)
tree_matrix <- as.matrix(tree$edge[, c(1, 2)])
g <- graph_from_edgelist(tree_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(tree)
plot.igraph(g, layout=layout_coord)
plot.igraph(g)
tree_mutation <- sim_mutation(tree, 2)

local_tree1 <- local_ARG(tree_mutation, 0.1)
local_tree1_matrix <- as.matrix(local_tree1$edge[, c(1, 2)])
local_g1 <- graph_from_edgelist(local_tree1_matrix, directed = FALSE)
local_g1 <- delete_vertices(local_g1, V(local_g1)[degree(local_g1) == 0])
plot(local_g1)

local_tree2 <- local_ARG(tree_mutation, 0.2)
local_tree2_matrix <- as.matrix(local_tree2$edge[, c(1, 2)])
local_g2 <- graph_from_edgelist(local_tree2_matrix, directed = FALSE)
local_g2 <- delete_vertices(local_g2, V(local_g2)[degree(local_g2) == 0])
plot(local_g2)

library(ape)

phylo_tree1 <- localARG_to_phylo(local_tree1, label=TRUE)
plot(phylo_tree1)

phylo_tree2 <- localARG_to_phylo(local_tree2, label=TRUE)
plot(phylo_tree2)

n_trees <- 10
n <- 20
rho <- 1
tree_list <- map(1:n_trees, function(i) {
  tree <- sim_ARG(n, rho)
  x_values <- cumsum(c(0, tree$waiting_time))
  y_values <- tree$k
  tibble(
    x = x_values,
    y = y_values,
    id = paste0("ARG_", i) # Identifier for each step function
  )
})
step_df <- bind_rows(tree_list)

ODE_process <- function(x, state, parameters) {
  with(as.list(c(state, parameters)), {
    # ODE: dy/dx = (rho * y - y^2 + y) / 2
    dy_dx <- (rho * y - y^2 + y) / 2
    return(list(c(dy_dx)))
  })
}

initial_y <- c(y = n)
x_values <- seq(0, 2*max(step_df$x), by = 0.1)
ode_parameters <- c(rho = rho)

ode_solution <- ode(
  y = initial_y,
  times = x_values,
  func = ODE_process,
  parms = ode_parameters
)

ode_df <- as_tibble(ode_solution)
colnames(ode_df)[1] <- "x"

ggplot() +
  # 10 step functions 
  geom_step(data = step_df, aes(x = x, y = y, group = id, color = id),
            linewidth = 0.7, alpha = 0.7) +
  # one continuous function
  geom_line(data = ode_df, aes(x = x, y = y), color = "black",
            linetype = "solid", linewidth = 1.2) +
  
  labs(
    title = "Birth and Death Process",
    x = "Time",
    y = "Number of alleles",
    color = "ARG ID"
  ) +

  theme_minimal() +
  guides(color = guide_legend(ncol = 2))
