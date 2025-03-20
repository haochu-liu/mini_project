source("abc_functions/kernel_functions.r")


abc_mcmc <- function(obs, tol, kernel_func, p_theta, d_theta, p_s, prior,
                     theta_0, n_iter, sigma=NULL) {
#' input: kernel, tol, sigma, obs,
#' proposal p_theta(theta'), log-density d_theta(theta, theta') (q(theta|theta')),
#' model p_s(theta), prior log-density prior(theta), initial theta, n_iter
#' output: matrix (n_iter x length(theta)), matrix (n_iter x length(s))

  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  theta_matrix <- matrix(NA, nrow=(n_iter+1), ncol=length(theta_0))
  s_matrix <- matrix(NA, nrow=(n_iter+1), ncol=length(obs))
  while (TRUE) {
    s_0 <- p_s(theta_0)
    k_0 <- kernel_func(obs, s_0, tol, sigma)
    if (exp(k_0)) {break}
  }
  theta_matrix[1, ] <- theta_0
  s_matrix[1, ] <- s_0

  for (i in 2:(n_iter+1)) {
    theta_1 <- p_theta(theta_0)
    s_1 <- p_s(theta_1)
    k_1 <- kernel_func(obs, s_1, tol, sigma)
    if (k_1) {
      log_alpha <- k_1+prior(theta_1)+d_theta(theta_0, theta_1)-
        k_0-prior(theta_0)-d_theta(theta_1, theta_0)
      if (log(runif(1)) < log_alpha) {
        theta_0 <- theta_1
        s_0 <- s_1
      }
    }
    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0
  }

  return(list(theta_matrix=theta_matrix, s_matrix=s_matrix))
}
