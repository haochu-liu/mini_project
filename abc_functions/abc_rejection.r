library(matrixStats)
source("abc_functions/kernel_functions.r")


abc_rejection <- function(obs, param, sumstat, tol, kernel_func,
                          sigma=NULL, log.weight=FALSE) {
#' input: obs(vector),
#' param(matrix: n x length(theta)), sumstat(matrix: n x length(s)),
#' tol(numeric), kernel(string), sigma=NULL(matrix), log.weight=FALSE(bool)
#' output: weight or log weight given by the abc rejection sampling

  weights <- c()
  n <- dim(param)[1]
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
  # get log weighs
  weights <- vapply(1:n,
                    function(i) kernel_func(obs, sumstat[i, ], tol, sigma),
                    numeric(1))
  weights <- weights - logSumExp(weights)
  if (!log.weight) {weights <- exp(weights)}
  
  return(weights)
}
