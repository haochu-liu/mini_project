abc_llr <- function(obs, param, sumstat, weights) {
#' input: obs(vector), param(matrix: n x length(theta)),
#' sumstat(matrix: n x length(s)), weights(vector)
#' output: adjusted param(matrix: n x length(theta))

  param_hat <- param
  sumstat_d <- sumstat - obs
  for (i in 1:ncol(param)) {
    loc_model <- lm(param[, i] ~ sumstat_d, weights=weights)
    beta_hat <- matrix(coef(loc_model)[2:length(coef(loc_model))], ncol=1)
    param_hat[, i] <- param[, i] - sumstat_d %*% beta_hat
  }

  return(param_hat)
}
