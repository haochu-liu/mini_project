library(ggplot2)
library(patchwork)
library(sdprisk)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/simbac_ARG.r")


n_values <- c(20, 100, 200)
rho_values <- seq(2, 10, by = 2)
delta_values <- c(1, 10, 20)
func_vec <- c("sim_FSM_ARG(optimise=T)", "sim_FSM_ARG(optimise=T, clonal=T)",
              "sim_FSM_ARG(optimise=F)", "simbac_ARG(optimise_site=F)",
              "simbac_ARG(optimise_site=T)")
time_df <- expand.grid(
  t = NA,
  n = n_values,
  rho = rho_values,
  delta = delta_values,
  func = func_vec
)

set.seed(100)
for (i in 1:nrow(time_df)) {
  n <- time_df$n[i]
  rho <- time_df$rho[i]
  delta <- time_df$delta[i]
  func_str <- time_df$func[i]
  time_vec <- rep(NA, 10)
  if (func_str == "sim_FSM_ARG(optimise=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        sim_FSM_ARG(n, rho, 100, bacteria = TRUE, delta = delta,
                    node_max = 100000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "sim_FSM_ARG(optimise=T, clonal=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        sim_FSM_ARG(n, rho, 100, bacteria = TRUE, delta = delta,
                    node_max = 100000, optimise_recomb = TRUE, clonal = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "sim_FSM_ARG(optimise=F)") {
    for (j in 1:10) {
      time_result <- system.time(
        sim_FSM_ARG(n, rho, 100, bacteria = TRUE, delta = delta,
                    node_max = 100000, optimise_recomb = FALSE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "simbac_ARG(optimise_site=F)") {
    for (j in 1:10) {
      time_result <- system.time(
        simbac_ARG(n, rho, 100, delta, node_max = 100000)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "simbac_ARG(optimise_site=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        simbac_ARG(n, rho, 100, delta, node_max = 100000, optimise_site = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  }
  time_df$t[i] <- median(time_vec)
  print(paste("Complete", i, "iterations"))
}


time_plot <- time_df
time_plot$n <- paste0("n = ", time_plot$n)
time_plot$delta <- paste0("delta = ", time_plot$delta)

ggplot(time_plot, aes(x=rho, y=t, color=func)) +
  geom_line(aes(group = interaction(n, delta, func)), linewidth = 1.2, alpha=0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1.2, alpha=0.7) +
  scale_color_manual(values=c("sim_FSM_ARG(optimise=T)"="darkblue",
                              "sim_FSM_ARG(optimise=T, clonal=T)"="darkgreen",
                              "sim_FSM_ARG(optimise=F)"="darkred",
                              "simbac_ARG(optimise_site=F)"="orange",
                              "simbac_ARG(optimise_site=T)"="yellow")) +
  facet_grid(n ~ delta) +
  labs(
    title = "Time by Recombination Rates",
    x = "rho",
    y = "time",
    color = "function"
  ) +
  scale_y_continuous(trans='log10') +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
