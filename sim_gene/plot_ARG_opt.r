library(ggplot2)
library(microbenchmark)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/ref/simu.R")
source("sim_gene/sim_birth_death.r")


time_df <- data.frame(time=rep(NA, 5000),
                      optimise=c(rep("Simulated Hitting Time", 1000),
                                 c(rep("TRUE", 1000), rep("FALSE", 1000)),
                                 c(rep("TRUE", 1000), rep("FALSE", 1000))),
                      func=c(rep("Birth-death Process", 1000),
                             rep("sim_FSM_ARG", 2000),
                             rep("simu", 2000)))

for (i in 1:1000) {
  t <- sim_birth_death(100, 5)
  time_df$time[i] <- t

  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                   node_max = 100000, optimise_recomb = TRUE)
  time_df$time[1000+i] <- r$sum_time

  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                   node_max = 100000, optimise_recomb = FALSE)
  time_df$time[2000+i] <- r$sum_time

  r = simu(n=100, rho = 5, delta = 1, blocks = c(100), optimise = T)
  time_df$time[3000+i] <- tail(r$ages, n=1)

  r = simu(n=100, rho = 5, delta = 1, blocks = c(100), optimise = F)
  time_df$time[4000+i] <- tail(r$ages, n=1)

  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

ggplot(time_df, aes(x=func, y=time, fill=optimise)) +
  geom_boxplot() +
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



node_df <- data.frame(n_node=rep(0, 200),
                      n=c(10*c(1:100), 10*c(1:100)),
                      optimise=c(c(rep("TRUE", 100), rep("FALSE", 100))))

for (i in 1:100) {
  for (j in 1:10) {
    r <- sim_FSM_ARG(node_df$n[i], 5, 100, bacteria = TRUE, delta = 10,
                     node_max = 100000, optimise_recomb = TRUE)
    node_df$n_node[i] <- node_df$n_node[i] + length(r$node_height)

    r <- sim_FSM_ARG(node_df$n[i], 5, 100, bacteria = TRUE, delta = 10,
                     node_max = 100000, optimise_recomb = FALSE)
    node_df$n_node[100+i] <- node_df$n_node[100+i] + length(r$node_height)
  }
  if (i%%10 == 0) {print(paste("Complete", i, "iterations"))}
}
node_df$n_node <- node_df$n_node/10

ggplot(node_df, aes(x=n, y=n_node, color=optimise)) +
  geom_line(size = 1.2) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  scale_color_manual(values=c("TRUE"="darkblue", "FALSE"="darkred")) +
  labs(
    title = "Number of Nodes from sim_FSM_ARG()",
    x = "# of Leaf Lineages",
    y = "# of Nodes",
    color = "Optimise"
  ) +
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



time_df <- data.frame(mean_time=rep(0, 20),
                      n=c(100*c(1:10), 100*c(1:10)),
                      optimise=c(c(rep("TRUE", 10), rep("FALSE", 10))))

for (i in 1:10) {
  benchmark_with_opt <- microbenchmark(
    sim_FSM_ARG(time_df$n[i], 5, 100, bacteria = TRUE, delta = 10,
                node_max = 100000, optimise_recomb = TRUE),
    times = 100,
    setup=set.seed(10)
  )
  
  benchmark_without_opt <- microbenchmark(
    sim_FSM_ARG(time_df$n[i], 5, 100, bacteria = TRUE, delta = 10,
                node_max = 100000, optimise_recomb = FALSE),
    times = 100,
    setup=set.seed(10)
  )

  time_df$mean_time[i] <- summary(benchmark_with_opt)$mean
  time_df$mean_time[10+i] <- summary(benchmark_without_opt)$mean

  print(paste("Complete", i, "iteration"))
}

ggplot(time_df, aes(x=n, y=mean_time, color=optimise)) +
  geom_line(size = 1.2) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  scale_color_manual(values=c("TRUE"="darkblue", "FALSE"="darkred")) +
  labs(
    title = "Average time of running sim_FSM_ARG()",
    x = "# of Leaf Lineages",
    y = "Time",
    color = "Optimise"
  ) +
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
