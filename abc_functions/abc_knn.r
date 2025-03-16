source("abc_functions/abc_rejection.r")
source("abc_functions/kernel_functions.r")


abc_knn <- function(obs, param, sumstat, k, kernel,
  sigma=NULL, log.weight=FALSE) {
#' input: obs(vector),
#' param(matrix: n x length(theta)), sumstat(matrix: n x length(s)),
#' k(numeric), kernel(string), sigma=NULL(matrix), log.weight=FALSE(bool)
#' output: weight or log weight given by the abc rejection sampling

  if (k!=as.integer(k)) {
    stop("k must be an integer")
  }
  
  n <- dim(param)[1]
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  d <- vapply(1:n, function(i) distance(obs, sumstat[i, ], 1, sigma), numeric(1))
  tol <- sort(d)[k] # set the k-th smallest distance as tol
  weights <- abc_rejection(obs, param, sumstat, tol, kernel, sigma, log.weight)

  return(weights)
}
