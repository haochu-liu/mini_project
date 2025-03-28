abc_IEnKI <- function(obs, param, sumstat, epsilon, sigma) {
  # sumstat should have sample points in row-wise, and provide a matrix
  M <- nrow(sumstat)
  N <- ncol(sumstat)
  T <- length(epsilon)
  e <- tail(epsilon, n=1)
  S_obs <- matrix(obs, length(obs), 1)
  S <- t(sumstat)
  mu_t <- rowMeans(S)
  C_t <- 0
  for (i in 1:M) {
    C_t <- C_t + (S[, i]-mu_t) %*% t(S[, i]-mu_t)
  }
  C_t <- C_t / (M-1)
  log_Z_t <- 0
  for (t in 2:T) {
    gamma_t <- e^(-2) / (epsilon[t]^(-2) - epsilon[t-1]^(-2))
    K_t <- C_t %*% solve(C_t + gamma_t*sigma)
    for (j in 1:M) {
      s_t <- rmvnorm(1, mean=S[, j], sigma=gamma_t*sigma)
      s_t <- t(s_t)
      S[, j] <- S[, j] + K_t %*% (S_obs - s_t)
    }

    log_c_t <- log(gamma_t) * N / 2 + log(2*pi) * ((1-gamma_t)*N/2) +
      log(det(sigma)) * ((1-gamma_t)/2)
    # c_t <- gamma_t^(N/2)*(2*pi)^((1-gamma_t)*N/2)*det(sigma)^((1-gamma_t)/2)
    log_Z_t <- log_Z_t + log_c_t + dmvnorm(obs,
                                            mean=mu_t,
                                            sigma=(C_t + gamma_t*sigma),
                                            log=TRUE)
    print(c(t, gamma_t, log_c_t, exp(log_c_t), log_Z_t))
    mu_t <- rowMeans(S)
    C_t <- 0
    for (i in 1:M) {
      C_t <- C_t + (S[, i]-mu_t) %*% t(S[, i]-mu_t)
    }
    C_t <- C_t / (M-1)
  }

  abc_list <- list(param=param, sumstat=sumstat, obs=obs, epsilon=epsilon,
                   sigma=sigma, log_Z_T=log_Z_t, Z_T=exp(log_Z_t),
                   mu_T=mu_t, C_T=C_t)
  class(abc_list) <- "abc_IEnKI_list"
  return(abc_list)
}
