abc_llr <- function(obs, param, sumstat, weights) {
#' input: obs(vector), param(matrix: n x length(theta)),
#' sumstat(matrix: n x length(s)), weights(vector)
#' output: corrected param(matrix: n x length(theta))

  param_hat <- param
  sumstat_d <- sumstat - obs
  for (i in 1:ncol(param)) {
    loc_model <- lm(param[, i] ~ sumstat_d, weights=weights)
    param_hat[, i] <- predict(loc_model, newdata=as.data.frame(sumstat_d))
  }

  return(param_hat)
}
