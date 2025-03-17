source("abc_functions/kernel_functions.r")


abc_mcmc <- function(obs, tol, kernel, p_theta, d_theta, p_s, prior,
                     theta_0, n_iter, sigma=NULL) {
#' input: kernel, tol, sigma, obs,
#' proposal p_theta(theta'), log-density d_theta(theta'), model p_s(theta),
#' prior log-density prior(theta), initial theta, n_iter
#' output: matrix (n_iter x length(theta)), matrix (n_iter x length(s))

  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  theta_matrix <- matrix(NA, nrow=n_iter, ncol=length(theta_0))
  s_matrix <- matrix(NA, nrow=n_iter, ncol=length(obs))
  while (TRUE) {
    s_0 <- p_s(theta_0)
    if (kernel=="uniform") {
        k <- uniform_kernel(obs, s_0, tol, sigma)
    }
  }


}