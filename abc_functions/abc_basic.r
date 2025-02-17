library(matrixStats)


abc_basic <- function(obs, param, sumstat, tol, kernel='uniform', sigma=NULL) {
  weights <- c()
  n <- dim(param)[1]
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  if (kernel=='uniform') {
    for (i in 1:n) {
      k <- uniform_kernel(obs, as.numeric(sumstat[i, ]), tol, sigma)
      weights <- c(weights, k)
    }
    weights <- weights / sum(weights)
    log_weights <- NULL
  } else if (kernel=='Gaussian') {
    for (i in 1:n) {
      k <- Gaussian_kernel(obs, as.numeric(sumstat[i, ]), tol, sigma)
      weights <- c(weights, k)
    }
    log_weights <- weights - logSumExp(weights)
    weights <- exp(log_weights)
  }

  abc_list <- list(param=param, sumstat=sumstat, obs=obs,
                   weights=weights, log_weights=log_weights,
                   tol=tol, kernel=kernel, sigma=sigma)
  class(abc_list) <- "abc_list"
  return(abc_list)
}
