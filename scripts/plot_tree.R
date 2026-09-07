# plot_tree.R - Phylogenomic tree (in English) with ggtree
library(ggtree)
library(ggplot2)
library(ape)

# 1) READ the tree (Newick format)
tree <- read.tree("/home/oussama/stage/results/arbre.tree")

# 2) ROOT the tree on the outgroup (most divergent genome)
outgroup <- "Acetilactobacillus-jinshanensis-HSLZ-75"
if (outgroup %in% tree$tip.label) {
  tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
}

# 3) DRAW the tree
p <- ggtree(tree) +
  geom_tiplab(size = 3) +
  geom_treescale() +
  ggtitle("Phylogenomic tree (120 core genes, bac120)") +
  xlim(0, 1.2)

# 4) SAVE as PNG
ggsave("/home/oussama/stage/results/tree.png", p, width = 11, height = 7, dpi = 200)
cat("Tree saved in results/tree.png\n")
