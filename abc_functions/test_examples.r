library(ggplot2)
source("abc_functions/abc_rejection.r")


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
