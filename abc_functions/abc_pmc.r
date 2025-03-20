library(mvtnorm)
library(matrixStats)
source("abc_functions/kernel_functions.r")


abc_pmc <- function(obs, tol, kernel_func, p_theta, d_theta, p_s, n_iter, sigma=NULL) {
#' input: obs, sigma, tol(vec in descendant order), kernel,
#' prior p_theta(), log-density d_theta(theta), model p_s(theta), n_iter
#' output: matrix (n_iter x length(theta)), matrix (n_iter x length(s))

  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  theta_matrix <- matrix(NA, nrow=n_iter, ncol=length(p_theta()))
  s_matrix <- matrix(NA, nrow=n_iter, ncol=length(obs))
  T_tol <- length(tol)
  # t = 1
  for (i in 1:n_iter) {
    while (TRUE) {
      theta_0 <- p_theta()
      s_0 <- p_s(theta_0)
      k_0 <- kernel_func(obs, s_0, 1, sigma)
      if (k_0<log(tol[1])) {break}
    }
    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0
  }
  w <- rep(log(1/n_iter), n_iter)
  tau <- 2*var(theta_matrix)

  for (t in 2:T_tol) {
    theta_1_matrix <- matrix(NA, nrow=n_iter, ncol=length(p_theta()))
    w_1 <- rep(NA, n_iter)
    for (i in 1:n_iter) {
      theta_index <- sample(1:n_iter, 1, prob=exp(w))
      while (TRUE) {
        theta_0 <- rmvnorm(n=1, mean=theta_matrix[theta_index, ], sigma=tau)
        s_0 <- p_s(theta_0)
        k_0 <- kernel_func(obs, s_0, 1, sigma)
        if (k_0<log(tol[t])) {break}
      }
      w_1[i] <- d_theta(theta_0) -
        logSumExp(w+dmvnorm((as.numeric(theta_0)-theta_matrix), sigma=tau,
                                 log=TRUE)+c(log(sqrt(tau))))
      theta_1_matrix[i, ] <- theta_0
      s_matrix[i, ] <- s_0
    }
    w <- w_1 - logSumExp(w_1)
    theta_matrix <- theta_1_matrix
    tau <- 2*var(theta_matrix)
  }

  return(list(theta_matrix=theta_matrix, s_matrix=s_matrix))
}
