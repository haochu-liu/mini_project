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
    title = "Running time for different ARG simulation methods",
    x = "rho",
    y = "log time",
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



library(devtools)
install("C:/Users/u2008181/simARG")


rho_values <-  seq(5, 30, by = 5)
func_vec <- c("FSM_ARG(optimise=T)", "FSM_ARG(optimise=T, clonal=T)",
              "simbac_ARG(optimise_site=F)", "simbac_ARG(optimise_site=T)",
              "ClonalOrigin_ARG(optimise=T)")
time_df <- expand.grid(
  t = NA,
  rho = rho_values,
  func = func_vec
)

set.seed(100)
for (i in 1:nrow(time_df)) {
  n <- 20L
  rho <- time_df$rho[i]
  delta <- 10
  func_str <- time_df$func[i]
  time_vec <- rep(NA, 10)
  if (func_str == "FSM_ARG(optimise=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        simARG::FSM_ARG(n, rho, 100L, bacteria = TRUE, delta = delta,
                        node_max = 100000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "FSM_ARG(optimise=T, clonal=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        simARG::FSM_ARG(n, rho, 100L, bacteria = TRUE, delta = delta,
                        node_max = 100000, optimise_recomb = TRUE, clonal = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "simbac_ARG(optimise_site=F)") {
    for (j in 1:10) {
      time_result <- system.time(
        simARG::simbac_ARG(n, rho, 100L, delta, node_max = 100000)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "simbac_ARG(optimise_site=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        simARG::simbac_ARG(n, rho, 100L, delta, node_max = 100000, optimise_site = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (func_str == "ClonalOrigin_ARG(optimise=T)") {
    for (j in 1:10) {
      time_result <- system.time(
        simARG::ClonalOrigin_ARG(n, rho, 100L, delta, node_max = 100000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  }
  time_df$t[i] <- median(time_vec)
  print(paste("Complete", i, "iterations"))
}

time_plot <- time_df

ggplot(time_plot, aes(x=rho, y=t, color=func)) +
  geom_line(, linewidth = 1.2, alpha=0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1.2, alpha=0.7) +
  scale_color_manual(values=c("FSM_ARG(optimise=T)"="darkblue",
                              "FSM_ARG(optimise=T, clonal=T)"="darkgreen",
                              "ClonalOrigin_ARG(optimise=F)"="darkred",
                              "simbac_ARG(optimise_site=F)"="orange",
                              "simbac_ARG(optimise_site=T)"="yellow")) +
  labs(
    title = "Running time for different ARG simulation methods",
    x = "rho",
    y = "log time",
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
