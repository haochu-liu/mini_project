library(mvtnorm)
library(matrixStats)
source("abc_functions/kernel_functions.r")


abc_pmc <- function(obs, tol, kernel, p_theta, d_theta, p_s, n_iter, sigma=NULL) {
#' input: obs, sigma, tol(vec in descendant order), kernel,
#' prior p_theta(theta), model p_s(theta), n_iter
#' output: matrix (n_iter x length(theta)), matrix (n_iter x length(s))

  if (kernel!="uniform" & kernel!="triangular" & kernel!="epanechnikov" &
    kernel!="biweight" & kernel!="gaussian") {stop("No valid kernel.")}
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  theta_matrix <- matrix(NA, nrow=n_iter, ncol=length(p_theta()))
  s_matrix <- matrix(NA, nrow=n_iter, ncol=length(obs))
  T_tol <- length(tol)
  # t = 1
  for (i in 1:n_iter) {
    while (TRUE) {
      theta_0 <- p_theta()
      s_0 <- p_s(theta_0)
      if (kernel=="uniform") {
        k_0 <- uniform_kernel(obs, s_0, 1, sigma)
        if (k_0<tol[1]) {break}
      } else if (kernel=="triangular") {
        k_0 <- triangular_kernel(obs, s_0, 1, sigma)
        if (k_0<tol[1]) {break}
      } else if (kernel=="epanechnikov") {
        k_0 <- epanechnikov_kernel(obs, s_0, 1, sigma)
        if (k_0<tol[1]) {break}
      } else if (kernel=="biweight") {
        k_0 <- biweight_kernel(obs, s_0, 1, sigma)
        if (k_0<tol[1]) {break}
      } else if (kernel=="gaussian") {
        k_0 <- gaussian_kernel(obs, s_0, 1, sigma)
        if (k_0<log(tol[1])) {break}
      }
    }
    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0
  }
  w <- rep(1/n_iter, n_iter)
  tau <- 2*var(theta_matrix)

  for (t in 2:T_tol) {
    theta_1_matrix <- matrix(NA, nrow=n_iter, ncol=length(p_theta()))
    for (i in 1:n_iter) {
      theta_index <- sample(1:n_iter, 1, prob=w)
      w_1 <- rep(NA, n_iter)
      while (TRUE) {
        theta_0 <- rmvnorm(n=1, mean=theta_matrix[theta_index, ], sigma=tau)
        s_0 <- p_s(theta_0)
        if (kernel=="uniform") {
          k_0 <- uniform_kernel(obs, s_0, 1, sigma)
          if (k_0<tol[1]) {break}
        } else if (kernel=="triangular") {
          k_0 <- triangular_kernel(obs, s_0, 1, sigma)
          if (k_0<tol[1]) {break}
        } else if (kernel=="epanechnikov") {
          k_0 <- epanechnikov_kernel(obs, s_0, 1, sigma)
          if (k_0<tol[1]) {break}
        } else if (kernel=="biweight") {
          k_0 <- biweight_kernel(obs, s_0, 1, sigma)
          if (k_0<tol[1]) {break}
        } else if (kernel=="gaussian") {
          k_0 <- gaussian_kernel(obs, s_0, 1, sigma)
          if (k_0<log(tol[1])) {break}
        }
      }
      w_1[i] <- d_theta(theta_0) - logSumExp(log(w)+dmvnorm((theta_0-theta_matrix), sigma=tau))
      theta_1_matrix[i, ] <- theta_0
      s_matrix[i, ] <- s_0
    }
    w <- w_1
    theta_matrix <- theta_1_matrix
    tau <- 2*var(theta_matrix)
  }

  return(list(theta_matrix=theta_matrix, s_matrix=s_matrix))
}
