# plot_all.R - Toutes les figures, version moderne (viridis)
suppressMessages({
  library(ggplot2); library(reshape2); library(ggtree)
  library(ape); library(viridisLite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_all.R <results_dir>")
}
RES <- args[1]

# thème moderne commun a toutes les heatmaps
theme_modern <- theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey30", size = 10),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y   = element_text(size = 7),
        panel.grid    = element_blank(),
        legend.position = "right")

clean <- function(x) gsub(".fas$", "", gsub(".*/", "", x))

# fonction generique : dessine une heatmap moderne
# direction = -1 inverse l'echelle de couleur (utile pour une DISTANCE,
# ou une valeur BASSE = genomes PROCHES, contrairement a ANI/AAI/POCP
# ou une valeur HAUTE = genomes proches)
heatmap_fig <- function(data, valcol, titre, sous, legende, fichier, decimales = 1, direction = 1) {
  p <- ggplot(data, aes(GenomeA, GenomeB, fill = .data[[valcol]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = round(.data[[valcol]], decimales)), size = 2, color = "white") +
    scale_fill_viridis_c(option = "D", direction = direction, name = legende) +
    coord_fixed() +
    labs(title = titre, subtitle = sous, x = NULL, y = NULL) +
    theme_modern
  ggsave(fichier, p, width = 9.5, height = 8.5, dpi = 300)
  cat("  ->", fichier, "\n")
}

# ---------- 1) ANI ----------
d <- read.table(file.path(RES, "ani_resultats.txt"), sep = "\t", header = FALSE)[, 1:3]
colnames(d) <- c("GenomeA", "GenomeB", "ANI")
d$GenomeA <- clean(d$GenomeA); d$GenomeB <- clean(d$GenomeB)
heatmap_fig(d, "ANI", "Average Nucleotide Identity (ANI)",
            "Species threshold = 95%", "ANI (%)", file.path(RES, "heatmap_ani.png"))

# ---------- 2) AAI ----------
raw <- read.table(file.path(RES, "aai_resultats.tsv"), sep = "\t", header = TRUE, check.names = FALSE)
d <- raw[, c(3, 4, 5)]; colnames(d) <- c("GenomeA", "GenomeB", "AAI")
heatmap_fig(d, "AAI", "Average Amino acid Identity (AAI)",
            "Approximate genus boundary around 65%", "AAI (%)", file.path(RES, "heatmap_aai.png"))

# ---------- 3) dDDH (distance GBDP, formule 2) ----------
# On affiche la DISTANCE (0 = identique), pas le %similarite estime :
# ce dernier n'est fiable que pres du seuil espece et devient absurde
# (parfois negatif) pour des genomes deja tres divergents -- voir dddh.py.
mat <- read.table(file.path(RES, "dddh_out/dddh_distance.tsv"), sep = "\t",
                  header = TRUE, row.names = 1, check.names = FALSE)
d <- melt(as.matrix(mat)); colnames(d) <- c("GenomeA", "GenomeB", "dDDH")
heatmap_fig(d, "dDDH", "digital DNA-DNA Hybridization (dDDH, formule GBDP 2)",
            "Seuil espece : distance > 0.0412 (~<70% DDH) = especes differentes",
            "distance dDDH", file.path(RES, "heatmap_dddh.png"), decimales = 3, direction = -1)

# ---------- 4) POCP ----------
mat <- read.table(file.path(RES, "pocp_out/pocp_matrice.tsv"), sep = "\t",
                  header = TRUE, row.names = 1, check.names = FALSE)
d <- melt(as.matrix(mat)); colnames(d) <- c("GenomeA", "GenomeB", "POCP")
heatmap_fig(d, "POCP", "Percentage of Conserved Proteins (POCP)",
            "Genus threshold = 50%", "POCP (%)", file.path(RES, "heatmap_pocp.png"))

# ---------- 5) TREE (colored by genus) ----------
tree <- read.tree(file.path(RES, "arbre.tree"))
og <- "Acetilactobacillus-jinshanensis-HSLZ-75"
if (og %in% tree$tip.label) tree <- root(tree, outgroup = og, resolve.root = TRUE)
genus <- sub("-.*", "", tree$tip.label)   # le genre = 1er mot du nom
dd <- data.frame(label = tree$tip.label, Genus = genus)
pt <- ggtree(tree) %<+% dd +
  geom_tippoint(aes(color = Genus), size = 3) +
  geom_tiplab(aes(color = Genus), size = 3, hjust = -0.05) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  ggtitle("Phylogenomic tree (120 core genes, bac120)") +
  xlim(0, 1.3) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))
ggsave(file.path(RES, "tree.png"), pt, width = 12, height = 7.5, dpi = 300)
cat("  -> tree.png\n")

cat("Toutes les figures modernes ont ete generees.\n")
