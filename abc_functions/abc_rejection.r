library(matrixStats)
source("abc_functions/kernel_functions.r")


abc_rejection <- function(obs, param, sumstat, tol, kernel,
                          sigma=NULL, log.weight=FALSE) {
#' input: obs(vector),
#' param(matrix: n x length(theta)), sumstat(matrix: n x length(s)),
#' tol(numeric), kernel(string), sigma=NULL(matrix), log.weight=FALSE(bool)
#' output: weight or log weight given by the abc rejection sampling

  weights <- c()
  n <- dim(param)[1]
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  if (kernel=="uniform") {
    weights <- vapply(1:n,
                      function(i) uniform_kernel(obs, sumstat[i, ], tol, sigma),
                      numeric(1))
    weights <- weights / sum(weights)
    if (log.weight) {weights <- log(weights)}
  } else if (kernel=="triangular") {
    weights <- vapply(1:n,
                      function(i) triangular_kernel(obs, sumstat[i, ], tol, sigma),
                      numeric(1))
    weights <- weights / sum(weights)
    if (log.weight) {weights <- log(weights)}
  } else if (kernel=="epanechnikov") {
    weights <- vapply(1:n,
                      function(i) epanechnikov_kernel(obs, sumstat[i, ], tol, sigma),
                      numeric(1))
    weights <- weights / sum(weights)
    if (log.weight) {weights <- log(weights)}
  } else if (kernel=="biweight") {
    weights <- vapply(1:n,
                      function(i) biweight_kernel(obs, sumstat[i, ], tol, sigma),
                      numeric(1))
    weights <- weights / sum(weights)
    if (log.weight) {weights <- log(weights)}
  } else if (kernel=="gaussian") {
    weights <- vapply(1:n,
                      function(i) gaussian_kernel(obs, sumstat[i, ], tol, sigma),
                      numeric(1))
    weights <- weights - logSumExp(weights)
    if (!log.weight) {weights <- exp(weights)}
  } else {stop("No valid kernel")}

  return(weights)
}
