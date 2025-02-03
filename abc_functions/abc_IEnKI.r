library(matlib)


abc_IEnKI <- function(obs, param, sumstat, epsilon, sigma) {
  # param should have sample points in row-wise, and provide a matrix
  M <- nrow(param)
  N <- ncol(param)
  d <- length(obs)
  T <- length(epsilon)
  e <- tail(epsilon, n=1)
  S_obs <- matrix(obs, length(obs), 1)
  S <- t(param)
  mu_t <- rowMeans(S)
  C_t <- 0
  for (i in 1:M) {
    C_t <- C_t + (S[, i]-mu_t) %*% t(S[, i]-mu_t)
  }
  C_t <- C_t / (M-1)
  Z_t <- 1
  for (t in 1:T) {
    gamma_t <- e^(-2) / (epsilon[t+1]^(-2) - epsilon[t]^(-2))
    K_t <- C_t %*% inv(C_t + gamma_t*sigma)
    for (j in 1:M) {
      s_t <- rmvnorm(1, mean=S[, j], sigma=gamma_t*sigma)
      S[, j] <- S[, j] + K_t %*% (S_obs - s_t)
    }

    c_t <- gamma_t^(d/2) * (2*pi)^((1-gamma_t)*d/2) * det(sigma)^((1-gamma_t)/2)
    Z_t <- Z_t * c_t * dmvnorm(obs, mean=mu_t, sigma=(C_t + gamma_t*sigma))

    mu_t <- rowMeans(S)
    C_t <- 0
    for (i in 1:M) {
      C_t <- C_t + (S[, i]-mu_t) %*% t(S[, i]-mu_t)
    }
    C_t <- C_t / (M-1)
  }
  
  abc_list <- list(param=param, sumstat=sumstat, obs=obs, epsilon=epsilon,
                   sigma=sigma, Z_T=Z_t, mu_T=mu_t, C_T=C_t)
  class(abc_list) <- "abc_IEnKI_list"
  return(abc_list)
}