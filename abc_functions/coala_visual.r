library(coala)
library(ggplot2)
source("abc_functions/abc_kernel.r")
source("abc_functions/abc_basic.r")
source("abc_functions/abc_IEnKI.r")


# Set the observation parameters
sfs <- c(112, 57, 24, 34, 16, 29, 8, 10, 15)

# Run model 2000 times for simulation
model <- coal_model(10, 50) +
  feat_mutation(par_prior("theta", runif(1, 1, 5))) +
  sumstat_sfs()
sim_data <- simulate(model, nsim = 2000, seed = 17)

# Getting the parameters
sim_param <- create_abc_param(sim_data, model)
# Getting the summary statistics
sim_sumstat <- create_abc_sumstat(sim_data, model)

# Basic ABC method with uniform or Gaussian kernel
basic_uniform <- abc_basic(sfs, matrix(sim_param$theta), sim_sumstat, 50,
                           kernel='uniform')
basic_gaussian <- abc_basic(sfs, matrix(sim_param$theta), sim_sumstat, 50,
kernel='Gaussian', sigma=rep(1, 9))

# ABC-IEnKI method
epsilon <- 10^(seq(from=10, to=-2, by=-1))
Sigma <- diag(rep(1, length(sfs)))
abc_IEnKI_list <- abc_IEnKI(sfs, matrix(sim_param$theta),
                            as.matrix(sim_sumstat), epsilon, Sigma)
pdf_values <- dnorm(abc_IEnKI_list$sumstat[, 1],
                    mean=abc_IEnKI_list$mu_T[1], sd=abc_IEnKI_list$C_T[1, 1])

# Plot and compare ABC methods on distributions for summary statistics
plot(abc_basic_list$sumstat[, 1], abc_basic_list$weights)
plot(abc_IEnKI_list$sumstat[, 1], pdf_values)
