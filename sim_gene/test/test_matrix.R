library(Matrix)
source("sim_gene/FSM/sim_FSM_ARG.r")
source("sim_gene/FSM/simbac_ARG.r")


arg1 <- sim_FSM_ARG(100, 10, 1e6, bacteria = TRUE, delta = 10,
                    optimise_recomb = TRUE)
dim(arg1$node_mat)
View(arg1$node_mat)
object.size(arg1$node_mat)
sparseMat <- as(arg1$node_mat, "sparseMatrix")
object.size(sparseMat)

rows <- c(1, 3, 2, 5)
cols <- c(2, 4, 1, 3)
values <- c(10, 20, 30, 40)
sparse_mat_1 <- sparseMatrix(i = rows, j = cols, x = values, dims = c(5, 4))
print(sparse_mat_1)
full_mat_1 <- as.matrix(sparse_mat_1)
print(full_mat_1)

object.size(sparse_mat_1)
object.size(full_mat_1)



