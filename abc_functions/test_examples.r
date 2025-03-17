library(ggplot2)
source("abc_functions/abc_rejection.r")
source("abc_functions/abc_mcmc.r")
source("abc_functions/abc_pmc.r")


# abc_rejection and abc_knn
mu <- 0
# observation
s_obs <- mean(rnorm(10, mean=mu, sd=1))
# simulations
mu_vec <- runif(100, min=-1, max=2)
s_vec <- c()
for (i in 1:length(mu_vec)) {
  s_i <- mean(rnorm(10, mean=mu_vec[i]), sd=1)
  s_vec <- c(s_vec, s_i)
}
weights1 <- abc_rejection(c(s_obs), matrix(mu_vec), matrix(s_vec),
                          tol=0.1, kernel = "uniform")
weights2 <- abc_rejection(c(s_obs), matrix(mu_vec), matrix(s_vec),
                          tol=0.1, kernel = "epanechnikov")
weights3 <- abc_rejection(c(s_obs), matrix(mu_vec), matrix(s_vec),
                          tol=0.1, kernel = "gaussian")
abc_df <- data.frame(
    x = rep(mu_vec, 3),
    y = c(weights1, weights2, weights3),
    type = rep(c("uniform kernel", "epanechnikov kernel", "Gaussian kernel"),
                 each = length(mu_vec))
  )

ggplot(abc_df, aes(x = x, y = y, color = type)) +
  geom_point(data=abc_df, size=3, alpha=0.5) +
  scale_color_manual(values=c("red", "blue", "green")) +
  labs(x="mu", y="pi(mu|s_obs)") +
  theme_minimal()

w_k <- abc_knn(c(s_obs), matrix(mu_vec), matrix(s_vec),
               k=10, kernel = "uniform")

# abc-mcmc
mu <- 0
# observation
s_obs <- mean(rnorm(10, mean=mu, sd=1))
# s_obs <- 0
# proposal function
p_theta <- function(theta_0) {
  runif(1, min=-2, max=2)
}
d_theta <- function(theta_1, theta_0) {
  log(1/4)
}
# model sample
p_s <- function(theta) {
  mean(rnorm(10, mean=theta), sd=1)
}
# initial theta
mu_0 <- runif(1, min=-2, max=2)
# M-H
matrix_list <- abc_mcmc(c(s_obs), 0.01, "gaussian", p_theta, d_theta, p_s, d_theta, mu_0, 10000)
# hist of posterior
hist(matrix_list$theta_matrix, probability = TRUE, main = "pi(mu|s_obs)",
     breaks = 20, col = "gray", border = "black")
abline(v = mu, col = "red", lwd = 2, lty = 2)

# abc-pmc
mu <- 0
# observation
s_obs <- mean(rnorm(10, mean=mu, sd=1))
# s_obs <- 0
# proposal function
p_theta <- function() {
  runif(1, min=-2, max=2)
}
d_theta <- function(theta) {
  log(1/4)
}
# model sample
p_s <- function(theta) {
  mean(rnorm(10, mean=theta), sd=1)
}
# list of tolerance
tol <- c(10, 1, 0.1, 0.01)
matrix_list <- abc_pmc(c(s_obs), tol, "gaussian", p_theta, d_theta, p_s, 1000)
# hist of posterior
hist(matrix_list$theta_matrix, probability = TRUE, main = "pi(mu|s_obs)",
     breaks = 20, col = "gray", border = "black")
abline(v = mu, col = "red", lwd = 2, lty = 2)
