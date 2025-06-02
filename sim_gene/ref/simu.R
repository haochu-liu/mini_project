simu <- function(n = 100, rho = 1, delta = 100, blocks = rep(500, 7), optimise = T) {
  # This function simulates the coalescent with recombination
  # It is optimized using the concept of ancestral material
  
  L <- sum(blocks)                        # Total length of the data
  delta <- 1 / delta                      # Inverse of average tract length
  blockscum <- cumsum(blocks)            # Cumulative block lengths
  
  # Probability of starting recombination at each site
  probstart <- rep(1, L)
  loc <- 1
  for (i in seq_along(blocks)) {
    probstart[loc] <- 1 / delta
    loc <- loc + blocks[i]
  }
  probstart <- probstart / sum(probstart)
  probstartcum <- cumsum(probstart)
  
  # Initialize ARG data structures
  max_nodes <- 1000
  s <- matrix(0L, nrow = max_nodes, ncol = 6)      # s: nodes with [child1, child2, clonal parent, imported parent, start, end]
  ages <- rep(0, max_nodes)                        # Ages of nodes
  toCoal <- rep(0L, 10 * n)                        # Remaining lineages
  toCoalAncMat <- matrix(0L, nrow = L, ncol = 10 * n)  # Which sites are ancestral material
  toCoal[1:n] <- 1:n
  toCoalAncMat[, 1:n] <- 1                         # Initially, everything is ancestral material
  ancmat <- matrix(0L, nrow = L + 1, ncol = max_nodes)
  ancmat[, 1:n] <- 1                               # Full ancestry matrix
  nodes <- n                                       # Initial number of nodes
  k <- n                                           # Number of lineages
  time <- 0                                        # Time tracker
  eq <- 1                                          # Flag for recombination across segments
  
  while (k > 1) {
    # Time until next event
    time <- time - log(runif(1)) / (k * (k - 1) / 2 + k * rho / 2)
    
    if (runif(1) < (k - 1) / (k - 1 + rho)) {
      # Coalescent event
      i <- sample(1:k, 1)
      j <- sample(setdiff(1:k, i), 1)
      
      nodes <- nodes + 1
      s[nodes, 1:2] <- c(toCoal[i], toCoal[j])
      s[toCoal[i], 3] <- nodes
      s[toCoal[j], 3] <- nodes
      ages[nodes] <- time
      
      if (sum(toCoalAncMat[, i] & toCoalAncMat[, j]) > 0 &&
          sum(toCoalAncMat[, i] & toCoalAncMat[, j]) < L) eq <- 0
      
      # Merge ancestral material
      toCoalAncMat[, i] <- toCoalAncMat[, i] | toCoalAncMat[, j]
      ancmat[1:L, nodes] <- toCoalAncMat[, i]
      ancmat[L + 1, nodes] <- ancmat[L + 1, toCoal[i]] | ancmat[L + 1, toCoal[j]]
      
      # Update lineage list
      toCoal[i] <- nodes
      toCoal[j] <- toCoal[k]
      toCoal[k] <- 0
      toCoalAncMat[, j] <- toCoalAncMat[, k]
      k <- k - 1
    } else {
      # Recombination event
      i <- sample(1:k, 1)
      beg <- which(runif(1) < probstartcum)[1]
      b <- which(beg <= blockscum)[1]
      nd <- min(beg + rgeom(1, delta), blockscum[b])
      
      # Skip if recombination doesn't affect ancestral material
      if (optimise &&
          all(toCoalAncMat[beg:nd, i] == 0) ||
          all(toCoalAncMat[-(beg:nd), i] == 0)) next
      
      nodes <- nodes + 2
      s[nodes - 1, 1] <- toCoal[i]
      s[nodes, 1] <- toCoal[i]
      s[toCoal[i], 3:4] <- c(nodes - 1, nodes)
      s[nodes, 5:6] <- c(beg, nd)
      
      # Split ancestry
      toCoalAncMatSet <- rep(0L, L)
      toCoalAncMatSet[beg:nd] <- toCoalAncMat[beg:nd, i]
      toCoalAncMat[beg:nd, i] <- 0
      toCoalAncMat[, k + 1] <- toCoalAncMatSet
      
      ancmat[1:L, nodes - 1] <- toCoalAncMat[, i]
      ancmat[1:L, nodes] <- toCoalAncMatSet
      ancmat[L + 1, nodes - 1] <- ancmat[L + 1, toCoal[i]]
      
      toCoal[i] <- nodes - 1
      toCoal[k + 1] <- nodes
      k <- k + 1
    }
  }
  
  # Finalize
  s <- s[1:nodes, , drop = FALSE]
  ages <- ages[1:nodes]
  ancmat <- ancmat[, 1:nodes, drop = FALSE]
  
  return(list(s = s, n=n, ages = ages, ancmat=ancmat, blocks=blocks, eq = eq))
}
