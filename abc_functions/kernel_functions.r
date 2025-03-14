# input: y(vector), z(vector), tol(numeric)
# output: K_tol(y, z)

uniform_kernel <- function(y, z, tol=1) {
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- norm(y-z, type="2")
  return(as.numeric(d/tol < 1))
}


