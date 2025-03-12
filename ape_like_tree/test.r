library(ape)
library(phangorn)
source("ape_like_tree/simcoale.r")
source("ape_like_tree/simtree_to_phylo.r")
source("ape_like_tree/tree_trajectory.r")


tree1 <- rcoal(5)
str(tree1)
class(tree1)
tree1$edge
plot(tree1)
nodelabels()

tree2 <- simcoale(10)
str(tree2)
class(tree2)
tree2$edge

tree2_phylo <- simtree_to_phylo(tree2)
str(tree2_phylo)
class(tree2_phylo)
tree2_phylo$edge
plot(tree2_phylo)
nodelabels()

tree_trajectory(1, tree=tree2)

df <- simmutation(tree=tree2, rate=1, model="infinite sites")
df <- simmutation(tree=tree2, rate=1, model="finite sites", l_seq = 10)
