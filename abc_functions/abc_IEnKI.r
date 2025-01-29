library(matlib)


abc_IEnKI <- function(obs, param, sumstat, epsilon, sigma) {
  # param should have sample points in row-wise, and provide a matrix
  M <- nrow(param)
  N <- ncol(param)
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
  for (t in 1:T) {
    gamma_t <- e^(-2) / (epsilon[t+1]^(-2) - epsilon[t]^(-2))
    K_t <- C_t %*% inv(C_t + gamma_t*sigma)
    for (j in 1:M) {
      s_t <- rnorm(1, mean=S[, j], sd=sqrt(gamma_t*sigma))
      S[, j] <- S[, j] + K_t %*% (S_obs - s_t)
    }
    mu_t <- rowMeans(S)
    C_t <- 0
    for (i in 1:M) {
      C_t <- C_t + (S[, i]-mu_t) %*% t(S[, i]-mu_t)
    }
    C_t <- C_t / (M-1)
  }
  return()
}