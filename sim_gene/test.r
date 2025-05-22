library(ivs)
library(tibble)
library(dplyr)
library(igraph)
library(ape)
library(ggplot2)
library(purrr)
library(deSolve)
source("sim_gene/sim_ARG.r")
source("sim_gene/sim_mutation.r")
source("sim_gene/local_ARG.r")
source("sim_gene/localARG_to_phylo.r")


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

phylo_tree <- localARG_to_phylo(local_tree)
plot(phylo_tree)



n_functions <- 10
n_points_per_function <- 20

step_data_list <- map(1:n_functions, function(i) {
  x_values <- sort(runif(n_points_per_function, 0, 10)) # Random x values
  y_values <- cumsum(rep(1/n_points_per_function, n_points_per_function)) # CDF-like y values
  tibble(
    x = x_values,
    y = y_values,
    id = paste0("Step_", i) # Identifier for each step function
  )
})

# Combine all step function data into a single tibble
step_df <- bind_rows(step_data_list)

# --- 2. Prepare data for the continuous function ---
# Let's use a simple sine wave or a polynomial
continuous_data <- tibble(
  x = seq(0, 10, length.out = 100), # Densely sampled x values
  y = 0.5 * sin(x / 2) + 0.5 # A sine wave for demonstration
)

# --- 3. Plotting with ggplot2 ---
ggplot() +
  # Add the 10 step functions
  geom_step(data = step_df, aes(x = x, y = y, group = id, color = id), size = 0.7, alpha = 0.7) +
  # Add the one continuous function
  geom_line(data = continuous_data, aes(x = x, y = y), color = "black", linetype = "solid", size = 1.2) +
  
  # Customize labels and title
  labs(
    title = "Ten Step Functions and One Continuous Function",
    x = "X-axis",
    y = "Y-axis",
    color = "Step Function ID" # Legend title for step functions
  ) +
  # Set limits if desired
  # xlim(0, 10) +
  # ylim(0, 1) +
  
  # Choose a theme
  theme_minimal() +
  # Customize legend (optional, but good for many lines)
  guides(color = guide_legend(ncol = 2)) # Adjust legend columns if many entries

my_ode_function <- function(x, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Extract y from state and rho from parameters
    # The ODE is dy/dx = (rho * y - y^2 + y) / 2
    dy_dx <- (rho * y - y^2 + y) / 2
    return(list(c(dy_dx))) # Must return a list
  })
}

# Initial condition: y at x = 0
initial_y <- c(y = 5) # Let's start y at 0.5 when x is 0

# Range for the independent variable x
x_values <- seq(0, 10, by = 0.1) # Solve from x=0 to x=10, with steps of 0.1

# Parameter value
# The behavior of the ODE will strongly depend on the value of rho.
# Let's try a few different values of rho to see different behaviors.
# For example: rho = 0, rho = 0.5, rho = 1.5, rho = -0.5
rho_value <- 1

# Combine parameters into a named vector
ode_parameters <- c(rho = rho_value)

# Solve the ODE
ode_solution <- ode(
  y = initial_y,        # Initial state
  times = x_values,     # Values of the independent variable (x) to solve for
  func = my_ode_function, # Your ODE function
  parms = ode_parameters  # ODE parameters
)

# Convert the solution matrix to a data frame for easier plotting
ode_df <- as_tibble(ode_solution) # Using tibble for consistency with ggplot2
colnames(ode_df)[1] <- "x" # Rename the first column from 'time' to 'x'

print(head(ode_df))

# Plotting with ggplot2
ggplot(ode_df, aes(x = x, y = y)) +
  geom_line(color = "darkblue", linewidth = 1) +
  labs(
    title = paste0("Solution of dy/dx = (rho * y + y^2 - y) / 2 (rho = ", rho_value, ")"),
    x = "x",
    y = "y"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title
