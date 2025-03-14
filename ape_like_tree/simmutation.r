source("ape_like_tree/length_and_height.r")
source("ape_like_tree/tree_trajectory.r")


#' input: tree, mutation rate, finite or infinite alleles
#' add mutation uniformly to edges
#' store in a dataframe: edge index, position(time), mutation
#' infinite alleles and finite alleles by tracing lineages from root to leaf
#' output: genotype for each leaf allele
simmutation <- function(tree, rate, model="infinite sites", l_seq=NULL) {
    if (!inherits(tree, "simtree")) {
        stop("Object must be of class 'simtree'")
    }

    l <- tree_length(tree=tree)
    n <- rpois(1, rate*l/2) # num of mutations | l ~ Poisson(rate*l/2)
    if (n == 0) {return(list(mutations=NA, node_seq=NA))}

    mutate_edges <- sample(1:nrow(tree$edge), n,
                           replace=TRUE, prob=tree$edge.length)
    if (model=="infinite sites") {
        mutate_site <- runif(n)
    } else if (model=="finite sites" & (!is.null(l_seq))) {
        mutate_site <- sample(1:l_seq, size=n, replace=TRUE)
    }
    # dataframe mutations to store information of every mutation
    mutations <- data.frame(edge_index = mutate_edges,
                            pos = rep(NA, n),
                            site = mutate_site)
    for (i in 1:n) {
        mutations[i, 2] <- runif(1, max=tree$edge.length[mutate_edges[i]])
    }
    
    # simulate the mutations before every node
    node_seq <- data.frame(node=rep(NA, 2*tree$n-1))
    node_seq$node <- as.list(node_seq$node)
    if (model=="infinite sites") {
        node_seq$node[[tree$n+1]] <- numeric(0)
        for (i in c((tree$n+2):(2*tree$n-1), 1:tree$n)) {
            traj <- tree_trajectory(i, tree=tree)
            parent_node <- traj$trajectory[2]
            i_edge <- traj$edge_index[1]
            edge_mutation <- mutations$site[mutations$edge_index==i_edge]
            node_seq$node[[i]] <- sort(c(node_seq$node[[parent_node]],
                                         edge_mutation))
        }
    } else if (model=="finite sites" & (!is.null(l_seq))) {
        node_seq$node[[tree$n+1]] <- rep(0, l_seq)
        for (i in c((tree$n+2):(2*tree$n-1), 1:tree$n)) {
            traj <- tree_trajectory(i, tree=tree)
            parent_node <- traj$trajectory[2]
            i_edge <- traj$edge_index[1]
            edge_mutation <- mutations$site[mutations$edge_index==i_edge]
            parent_seq <- node_seq$node[[parent_node]]
            parent_seq[edge_mutation] <- 1 - parent_seq[edge_mutation]
            node_seq$node[[i]] <- parent_seq
        }
    }
    return(list(mutations=mutations, node_seq=node_seq))
}
