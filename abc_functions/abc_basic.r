abc_basic <- function(obs, param, sumstat, tol, kernel='uniform', sigma=NULL) {
  weights <- c()
  n <- dim(param)[1]
  if (kernel=='uniform') {
    for (i in 1:n) {
      k <- uniform_kernel(obs, as.numeric(sumstat[i, ]), tol)
      weights <- c(weights, k)
    }
  } else if (kernel=='Gaussian') {
    for (i in 1:n) {
      k <- Gaussian_kernel(obs, as.numeric(sumstat[i, ]), tol, sigma)
      weights <- c(weights, k)
    }
  }

  if (kernel=='Gaussian') {
    weights <- weights / sum(weights)
  }
  abc_list <- list(param=param, sumstat=sumstat, weights=weights, obs=obs,
                   tol=tol, kernel=kernel, sigma=sigma)
  class(abc_list) <- "abc_list"
  return(abc_list)
}
