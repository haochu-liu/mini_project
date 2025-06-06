library(coala)
library(abc)
source("abc_functions/abc_kernel.r")
source("abc_functions/abc_basic.r")
source("abc_functions/abc_IEnKI.r")


sfs <- c(112, 57, 24, 34, 16, 29, 8, 10, 15)

model <- coal_model(10, 50) +
  feat_mutation(par_prior("theta", runif(1, 1, 5))) +
  sumstat_sfs()
sim_data <- simulate(model, nsim = 2000, seed = 17)

# Getting the parameters
sim_param <- create_abc_param(sim_data, model)
# Getting the summary statistics
sim_sumstat <- create_abc_sumstat(sim_data, model)

posterior <- abc(sfs, sim_param, sim_sumstat, 0.05, method = "rejection")
print(posterior)
hist(posterior, breaks = 20)

abc_basic_list <- abc_basic(sfs, matrix(sim_param$theta), sim_sumstat, 50,
                            kernel='uniform')
plot(abc_basic_list$param[, 1], abc_basic_list$weights)
abc_basic_list <- abc_basic(sfs, matrix(sim_param$theta), sim_sumstat, 50,
                            kernel='Gaussian', sigma=diag(rep(1, 9)))
plot(abc_basic_list$param[, 1], exp(abc_basic_list$log_weights))
plot(abc_basic_list$sumstat[, 1], exp(abc_basic_list$log_weights))

epsilon <- c(50, 20, 10)
Sigma <- diag(rep(1, length(sfs)))
abc_IEnKI_list <- abc_IEnKI(sfs, matrix(sim_param$theta),
                            as.matrix(sim_sumstat), epsilon, Sigma)
pdf_values <- dnorm(abc_IEnKI_list$sumstat[, 1],
                    mean=abc_IEnKI_list$mu_T[1], sd=abc_IEnKI_list$C_T[1, 1])
plot(abc_IEnKI_list$sumstat[, 1], pdf_values)
