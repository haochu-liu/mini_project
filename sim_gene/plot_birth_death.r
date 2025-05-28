library(ivs)
library(tibble)
library(dplyr)
library(ggplot2)
library(purrr)
library(deSolve)
source("sim_gene/ISM/sim_ISM_ARG.r")
source("sim_gene/ISM/sim_ISM_mutation.r")
source("sim_gene/local_tree.r")
source("sim_gene/localtree_to_phylo.r")
source("sim_gene/ARG_igraph.r")

n_ARG <- 10
n <- 20
rho <- 1
ARG_list <- map(1:n_ARG, function(i) {
  ARG <- sim_ISM_ARG(n, rho)
  x_values <- cumsum(c(0, ARG$waiting_time))
  y_values <- ARG$k
  tibble(
    x = x_values,
    y = y_values,
    id = paste0("ARG_", i) # Identifier for each step function
  )
})
step_df <- bind_rows(ARG_list)

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
