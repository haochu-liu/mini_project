#' Input: n and rho
#' Simulate a birth death process until it hits the boundary
#' Output: the time when the process first hit k=1
sim_birth_death <- function(n, rho) {
  t <- 0
  k <- n
  while (k > 1) {
    t <- t + rexp(1, rate=k*(k-1+rho)/2)
    if (runif(1) < (k - 1) / (k - 1 + rho)) {
      k <- k - 1
    } else {
      k <- k + 1
    }
  }
  return(t)
}
