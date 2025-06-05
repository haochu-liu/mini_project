library(ggplot2)
library(patchwork)
library(sdprisk)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/local_tree.r")
source("sim_gene/FSM/sim_FSM_mutation.r")
source("sim_gene/localtree_height.r")
source("sim_gene/localtree_traj.r")


set.seed(11)
ARG <- sim_FSM_ARG(5, 1, 10, bacteria = TRUE, delta = 3)
ARG_mutation <- sim_FSM_mutation(ARG, 2)
local_tree1 <- local_tree(ARG_mutation, 2)
local_tree2 <- local_tree(ARG_mutation, 5)
local_tree3 <- local_tree(ARG_mutation, 8)
tree_height(local_tree1)
tree_height(local_tree2)
tree_height(local_tree3)


height_t_df <- data.frame(s1=rep(NA, 1000),
                          s50=rep(NA, 1000),
                          s80=rep(NA, 1000))
set.seed(11)
for (i in 1:1000) {
  ARG <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                     node_max = 100000, optimise_recomb = TRUE)
  tree1 <- local_tree(ARG, 1)
  tree50 <- local_tree(ARG, 50)
  tree80 <- local_tree(ARG, 80)
  height_t_df$s1[i] <- tree_height(tree1)
  height_t_df$s50[i] <- tree_height(tree50)
  height_t_df$s80[i] <- tree_height(tree80)
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

n <- 100
height_rate <- n:2 * (n-1):1 / 2
# x <- seq(0, 10, length.out = 500)
# height_density <- dhypoexp(x, rate=height_rate)
height_density <- function(x) {
  dhypoexp(x, rate=height_rate)
}

hist1 <- ggplot(height_t_df, aes(x = s1)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 1",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist2 <- ggplot(height_t_df, aes(x = s50)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 50",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist3 <- ggplot(height_t_df, aes(x = s80)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 80",
       x = "Height",
       y = "Density") +
  theme_minimal()
combined_hist <- hist1 + hist2 + hist3
combined_hist <- combined_hist +
  plot_annotation(
    title = "sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1, node_max = 100000, optimise_recomb = TRUE)")
print(combined_hist)


height_f_df <- data.frame(s1=rep(NA, 1000),
                          s50=rep(NA, 1000),
                          s80=rep(NA, 1000))
set.seed(11)
for (i in 1:1000) {
  ARG <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                     node_max = 100000, optimise_recomb = FALSE)
  tree1 <- local_tree(ARG, 1)
  tree50 <- local_tree(ARG, 50)
  tree80 <- local_tree(ARG, 80)
  height_f_df$s1[i] <- tree_height(tree1)
  height_f_df$s50[i] <- tree_height(tree50)
  height_f_df$s80[i] <- tree_height(tree80)
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

hist1 <- ggplot(height_f_df, aes(x = s1)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 1",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist2 <- ggplot(height_f_df, aes(x = s50)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 50",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist3 <- ggplot(height_f_df, aes(x = s80)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 80",
       x = "Height",
       y = "Density") +
  theme_minimal()
combined_hist <- hist1 + hist2 + hist3
combined_hist <- combined_hist +
  plot_annotation(
    title = "sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1, node_max = 100000, optimise_recomb = FALSE)")
print(combined_hist)
