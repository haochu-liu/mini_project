distance <- function(y, z, tol=1, sigma) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: distance between y and z with cov matrix sigma
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  return(u)
}


uniform_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: log(K_tol(y, z))
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(log(as.numeric(u <= 1) / 2))
  } else {return(as.numeric(u <= 1) / 2)}
}


triangular_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: log(K_tol(y, z))
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(log((1-u) * as.numeric(u <= 1)))
  } else {return((1-u) * as.numeric(u <= 1))}
}


epanechnikov_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: log(K_tol(y, z))
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(log((1-u^2) * as.numeric(u <= 1) * 3 / 4))
  } else {return((1-u^2) * as.numeric(u <= 1) * 3 / 4)}
}


biweight_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: log(K_tol(y, z))
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(log((1-u^2)^3 * as.numeric(u <= 1) * 15 / 16))
  } else {return((1-u^2)^3 * as.numeric(u <= 1) * 15 / 16)}
}


gaussian_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
#' input: y(vector), z(vector), tol(numeric), sigma(matrix)
#' output: log(K_tol(y, z))
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(dnorm(u, log=TRUE))
  } else {dnorm(u, log=FALSE)}
}
