library(ape)


tree1 <- rcoal(5)
str(tree1)
class(tree1)
tree1$edge
tree1$tip.label
plot(tree1)

tree2 <- simcoale(3)
str(tree2)
class(tree2)
tree2$edge


rcoal <- function(n, tip.label = NULL, br = "coalescent", ...)
{
    n <- as.integer(n)
    nbr <- 2*n - 2
    edge <- matrix(NA, nbr, 2)
    ## coalescence times by default:
    x <- if (is.character(br)) 2*rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
    else if (is.numeric(br)) rep(br, length.out = n - 1) else br(n - 1, ...)
    if (n == 2) {
        edge[] <- c(3L, 3L, 1:2)
        edge.length <- rep(x, 2)
    } else if (n == 3) {
        edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
        edge.length <- c(x[c(2, 1, 1)], sum(x))
    } else {
        edge.length <- numeric(nbr)
        h <- numeric(2*n - 1)
        node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 2L*n - 1L
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            ind <- (i - 1)*2 + 1:2
            edge[ind, 2] <- y
            edge[ind, 1] <- nextnode
            edge.length[ind] <- node.height[i] - h[y]
            h[nextnode] <- node.height[i]
            pool <- c(pool[! pool %in% y], nextnode)
            nextnode <- nextnode - 1L
        }
    }
    phy <- list(edge = edge, edge.length = edge.length)
    phy$tip.label <- sample(tip.label)
    phy$Nnode <- n - 1L
    phy <- reorder(phy)
    ## to avoid crossings when converting with as.hclust:
    phy$edge[phy$edge[, 2] <= n, 2] <- 1:n
    phy
}