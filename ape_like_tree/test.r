library(ape)


tree1 <- rcoal(5)
str(tree1)
class(tree1)
tree1$edge
tree1$tip.label
plot(tree1)

tree2 <- simcoale(5)
str(tree2)
class(tree2)
tree2$edge

tree2_phylo <- simtree_to_phylo(tree2)
str(tree2_phylo)
class(tree2_phylo)
tree2_phylo$edge
plot(tree2_phylo)
