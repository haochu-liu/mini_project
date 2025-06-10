#' Input: node ancestral material, mean of the length of recombinant segment,
#' number of sites, recombination rate per genome
#' Compute the effective recombination rate,
#' and provide the probability of recombination initiating sites
#' Output: Effective recombination rate, cum sum of initiating sites probability
effective_R <- function(mat, delta, L, rho) {
  if (L < 2) {
    stop("Number of sites must be larger than one")
  }

  if (sum(mat) == L) {
    # if contain full material
    R_eff <- rho * (1-(1-1/delta)^(L-1)) / 2
    return(list(R_eff = R_eff,
                probstartcum=cumsum(rep(1/L, L))))
  } else {
    R <- rho / L
    v_s <- c(1:L)[mat & (mat != c(0, mat[1:(L-1)]))] # compute s1, ..., sb
    v_e <- c(1:L)[mat & (mat != c(mat[2:L], 0))]     # compute e1, ..., eb
    v_e <- c(0, v_e)                                 # add e0
  }
}