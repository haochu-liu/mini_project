library(mvtnorm)


uniform_kernel <- function(obs, sumstat, tol, sigma) {
  distance_m <- as.matrix(obs - sumstat)
  d <- t(distance_m) %*% solve(sigma) %*% distance_m
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
