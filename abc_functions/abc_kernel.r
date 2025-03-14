# library(mvtnorm)


# uniform_kernel <- function(obs, sumstat, tol, sigma) {
#   distance_m <- as.matrix(obs - sumstat)
#   d <- sqrt(t(distance_m) %*% solve(sigma) %*% distance_m)
#   if (d < tol) {
#     return(1)
#   } else{
#     return(0)
#   }
# }

# Gaussian_kernel <- function(obs, sumstat, tol, sigma) {
#   p <- dmvnorm(sumstat, mean=obs, sigma=tol*sigma, log=TRUE)
#   return(p) # return log density
# }
