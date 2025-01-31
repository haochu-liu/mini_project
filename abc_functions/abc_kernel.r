library(mvtnorm)


uniform_kernel <- function(obs, sumstat, tol) {
  d <- norm((obs - sumstat), type="2")
  if (d < tol) {
    return(1)
  } else{
    return(0)
  }
}

Gaussian_kernel <- function(obs, sumstat, tol, sigma) {
  sigma <- diag(sigma)
  p <- dmvnorm(sumstat, mean=obs, sigma=tol*sigma)
  return(p)
}
