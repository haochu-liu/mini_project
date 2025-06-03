source("sim_gene/FSM/sim_FSM_ARG.r")

fsm_arg_df <- data.frame(ages=rep(NA, 2000),
                         optimise=c(rep("TRUE", 1000), rep("FALSE", 1000)))
for (i in 1:1000) {
  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                   node_max = 100000, optimise_recomb = TRUE)
  fsm_arg_df$ages[i] <- r$sum_time
  r <- sim_FSM_ARG(100, 5, 100, bacteria = TRUE, delta = 1,
                   node_max = 100000, optimise_recomb = FALSE)
  fsm_arg_df$ages[i+1000] <- r$sum_time
}
boxplot(ages ~ optimise, data = fsm_arg_df,
        main = "sim_FSM_ARG",
        xlab = "Optimise",
        ylab = "ARG height")
mean(fsm_arg_df$ages[1:1000])
mean(fsm_arg_df$ages[1001:2000])
