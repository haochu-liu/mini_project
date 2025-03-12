library(ape)
library(ggplot2)
source("ape_like_tree/simcoale.r")
source("ape_like_tree/simtree_to_phylo.r")
source("ape_like_tree/simmutation.r")


tree2 <- simcoale(5)
str(tree2)
class(tree2)
tree2$edge

df1 <- simmutation(tree=tree2, rate=1, model="infinite sites")
df2 <- simmutation(tree=tree2, rate=1, model="finite sites", l_seq = 4)

df1_node <- df1$node_seq$node[1:5]
df1_label <- c()
for (i in 1:5) {
  string_node <- paste0("[", paste(round(df1_node[[i]], 3), collapse = ", "), "]")
  df1_label <- c(df1_label, paste(i, string_node, sep=":"))
}
tree2_phylo <- simtree_to_phylo(tree2, label=df1_label)
str(tree2_phylo)
class(tree2_phylo)
tree2_phylo$edge
plot(tree2_phylo)
nodelabels()
# edgelabels()

df2_node <- df2$node_seq$node[1:5]
df2_label <- c()
for (i in 1:5) {
  string_node <- paste0("[", paste(df2_node[[i]], collapse = ", "), "]")
  df2_label <- c(df2_label, paste(i, string_node, sep=":"))
}
tree2_phylo <- simtree_to_phylo(tree2, label=df2_label)
str(tree2_phylo)
class(tree2_phylo)
tree2_phylo$edge
plot(tree2_phylo)
nodelabels()
# edgelabels()

