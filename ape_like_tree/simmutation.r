source("ape_like_tree/length_and_height.r")


# input: tree, mutation rate, finite or infinite alleles
# add mutation uniformly to edges
# store in a dataframe: edge index, position(time), mutation
# infinite alleles and finite alleles by tracing lineages from root to leaf
# output: genotype for each leaf allele
simmutation <- function(tree, rate, model="infinite sites") {
    if (!inherits(tree, "simtree")) {
        stop("Object must be of class 'simtree'")
    }

    l <- tree_length(tree=tree)
    n <- rpois(1, rate*l/2) # num of mutations | l ~ Poisson(rate*l/2)
    mutate_edges <- sample(1:nrow(tree$edge), n,
                           replace=TRUE, prob=tree$edge.length)
    # dataframe mutations to store information of every mutation
    mutations <- data.frame(edge_index = mutate_edges,
                            mutate_pos = rep(NA, n),
                            mutate_site = rep(NA, n))
    for (i in 1:n) {
        mutations[i, 2] <- runif(1, max=tree$edge.length[mutate_edges[i]])
    }
    
}