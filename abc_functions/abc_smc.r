library(matrixStats)
source("abc_functions/kernel_functions.r")


ess <- function(w) {
#' input: normalized weights(vector)
#' output: effective sample size

  return(1/sum(w^2))
}


abc_smc <- function(obs, tol, kernel_func, p_theta, d_theta, p_s, n_iter, sigma=NULL) {
#' input: obs, sigma, tol(vec in descendant order), kernel,
#' prior p_theta(), log-density d_theta(), model p_s(theta), n_iter
#' output: matrix (n_iter x length(theta)), matrix (n_iter x length(s))

  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  theta_matrix <- matrix(NA, nrow=n_iter, ncol=length(p_theta()))
  s_matrix <- matrix(NA, nrow=n_iter, ncol=length(obs))
  T_tol <- length(tol)
  t <- 1
  
}
