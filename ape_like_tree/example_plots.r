library(ape)
library(sdprisk)
source("ape_like_tree/simcoale.r")
source("ape_like_tree/simtree_to_phylo.r")
source("ape_like_tree/simmutation.r")
source("ape_like_tree/length_and_height.r")


n <- 10
tree2 <- simcoale(n)
str(tree2)
class(tree2)
tree2$edge

df1 <- simmutation(tree=tree2, rate=1, model="infinite sites")
df2 <- simmutation(tree=tree2, rate=1, model="finite sites", l_seq = 5)

df1_node <- df1$node_seq$node[1:n]
df1_label <- c()
for (i in 1:n) {
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

df2_node <- df2$node_seq$node[1:n]
df2_label <- c()
for (i in 1:n) {
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

height_vec <- c()
length_vec <- c()
n <- 10
for (i in 1:1000) {
  tree <- simcoale(n)
  height_vec <- c(height_vec, tree_height(tree))
  length_vec <- c(length_vec, tree_length(tree))
}
height_rate <- n:2 * (n-1):1 / 2
x <- seq(0, 10, length.out = 500)
height_density <- dhypoexp(x, rate=height_rate)
hist(height_vec, probability = TRUE, col="gray", breaks=20,
     main = "Plot of TMRCA and hypoexponential distribution",
     xlab = "TMRCA", ylab = "Density")
lines(x, height_density)

length_rate <- (n-1):1 / 2
x <- seq(0, 20, length.out = 500)
length_density <- dhypoexp(x, rate=length_rate)
hist(length_vec, probability = TRUE, col="gray", breaks=20,
     main = "Plot of total tree length and hypoexponential distribution",
     xlab = "total tree length", ylab = "Density")
lines(x, length_density)
