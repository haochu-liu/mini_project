library(igraph)
library(ggplot2)
library(microbenchmark)
library(patchwork)
library(sdprisk)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/simbac_ARG.r")
source("sim_gene/FSM/sim_FSM_mutation.r")
source("sim_gene/ref/simu.R")
source("sim_gene/sim_birth_death.r")
source("sim_gene/local_tree.r")
source("sim_gene/ARG_igraph.r")


set.seed(100)
ARG <- simbac_ARG(5, 1, 10, 5, output_eff_R = TRUE)
ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
g <- graph_from_edgelist(ARG_matrix, directed = FALSE)
g <- delete_vertices(g, V(g)[degree(g) == 0])
layout_coord <- ARG_igraph(ARG)
plot.igraph(g)
plot.igraph(g, layout=layout_coord)
ARG_mutation <- sim_FSM_mutation(ARG, 2)


time_df <- data.frame(time=rep(NA, 4000),
                      func=c(rep("Birth-death Process", 1000),
                             rep("sim_FSM_ARG with optimisation", 1000),
                             rep("sim_FSM_ARG without optimisation", 1000),
                             rep("simbac", 1000)))

set.seed(10)
for (i in 1:1000) {
  t <- sim_birth_death(100, 5)
  time_df$time[i] <- t
  
  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 5,
                   node_max = 100000, optimise_recomb = TRUE)
  time_df$time[1000+i] <- r$sum_time
  
  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 5,
                   node_max = 100000, optimise_recomb = FALSE)
  time_df$time[2000+i] <- r$sum_time
  
  r <- simbac_ARG(100, 5, 100, delta = 5, node_max = 100000, output_eff_R = TRUE)
  time_df$time[3000+i] <- r$sum_time
  
  if (i%%10 == 0) {print(paste("Complete", i, "iterations"))}
}

ggplot(time_df, aes(x=func, y=time)) +
  geom_violin() +
  labs(title = "# of leaf lineages = 100, rho = 5",
       x = "Functions",
       y = "ARG height") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )


height_t_df <- data.frame(s1=rep(NA, 1000),
                          s50=rep(NA, 1000),
                          s80=rep(NA, 1000))
set.seed(11)
for (i in 1:1000) {
  ARG <- simbac_ARG(100, 5, 100, delta = 5, node_max = 100000, output_eff_R = TRUE)
  tree1 <- local_tree(ARG, 1)
  tree50 <- local_tree(ARG, 50)
  tree80 <- local_tree(ARG, 80)
  height_t_df$s1[i] <- tree1$sum_time
  height_t_df$s50[i] <- tree50$sum_time
  height_t_df$s80[i] <- tree80$sum_time
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
    title = "simbac_ARG(100, 5, 100, delta = 5, node_max = 100000, output_eff_R = TRUE)")
print(combined_hist)


benchmark_with_opt <- microbenchmark(
  sim_FSM_ARG(100, 10, 100, bacteria = TRUE, delta = 10,
              node_max = 100000, optimise_recomb = TRUE),
  times = 10,
  setup=set.seed(10)
)

benchmark_without_opt <- microbenchmark(
  sim_FSM_ARG(100, 10, 100, bacteria = TRUE, delta = 10,
              node_max = 100000, optimise_recomb = FALSE),
  times = 10,
  setup=set.seed(10)
)

benchmark_simbac <- microbenchmark(
  simbac_ARG(100, 10, 100, delta = 10, node_max = 100000),
  times = 10,
  setup=set.seed(10)
)

benchmark_with_opt
benchmark_without_opt
benchmark_simbac
