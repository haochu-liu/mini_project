library(ggplot2)
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
