source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/simbac_ARG.r")


arg1 <- sim_FSM_ARG(100, 10, 1e6, bacteria = TRUE, delta = 10,
                    optimise_recomb = TRUE)
dim(arg1$node_mat)
View(arg1$node_mat)
