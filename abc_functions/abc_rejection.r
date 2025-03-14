source("abc_functions/kernel_functions.r")


abc_rejection <- function(obs, param, sumstat, tol, kernel,
                          sigma=NULL, log.weight=FALSE) {
#' input: obs(vector),
#' param(matrix: n x length(theta)), sumstat(matrix: n x length(s)),
#' tol(numeric), kernel(string), sigma=NULL(matrix), log=FALSE(bool)
#' output: weight or log weight given by the abc rejection sampling

  weights <- c()
  n <- dim(param)[1]
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}
}