# plot_all_web.R - Parameterized version of plot_all.R for the web interface.
# Same figures, same styling; only the results directory comes from argv[1]
# instead of being hardcoded, so each web job renders into its own folder.
suppressMessages({
  library(ggplot2); library(reshape2); library(ggtree)
  library(ape); library(viridisLite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_all_web.R <results_dir>")
}
RES <- args[1]

# theme "planche scientifique d'epoque" -- inspire des gravures naturalistes
# du XIXe siecle (type Haeckel, "Pedigree of Man") : papier vieilli, encre
# sepia, typographie serif. Coherent avec le fait qu'on trace aussi un arbre
# genomique, comme les vieilles planches d'arbres evolutifs.
BG_PAPER   <- "#f1e7d0"   # papier vieilli
INK_DARK   <- "#2b1d0e"   # encre sepia fonce
INK_MED    <- "#6b4a2a"   # encre sepia moyen
GRID_LINE  <- "#b9a06f"   # liseret papier/encre clair
ACCENT_RED <- "#7a2419"   # rouge bordeaux d'epoque (sous-titres/seuils)

theme_vintage <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(plot.background    = element_rect(fill = BG_PAPER, color = INK_DARK, linewidth = 1.4),
        panel.background   = element_rect(fill = BG_PAPER, color = NA),
        panel.border       = element_rect(color = INK_DARK, fill = NA, linewidth = 0.9),
        plot.title         = element_text(face = "bold", size = 16, color = INK_DARK,
                                           family = "serif", hjust = 0.5),
        plot.subtitle      = element_text(size = 10, color = ACCENT_RED, family = "serif",
                                           hjust = 0.5, face = "italic"),
        plot.caption       = element_text(size = 7, color = INK_MED, family = "serif", hjust = 0.5),
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 7, color = INK_DARK,
                                           family = "serif", face = "italic"),
        axis.text.y        = element_text(size = 7, color = INK_DARK, family = "serif", face = "italic"),
        panel.grid         = element_blank(),
        legend.background  = element_rect(fill = BG_PAPER, color = NA),
        legend.key         = element_rect(fill = BG_PAPER, color = NA),
        legend.text        = element_text(color = INK_MED, family = "serif", size = 7),
        legend.title       = element_text(color = INK_DARK, family = "serif", size = 9, face = "bold"),
        legend.position    = "right",
        plot.margin        = margin(16, 20, 12, 12))

clean <- function(x) gsub("\\.(fas|fasta|fna)$", "", gsub(".*/", "", x))

# fonction generique : dessine une heatmap "planche scientifique d'epoque"
heatmap_fig <- function(data, valcol, titre, sous, legende, fichier) {
  rng <- range(data[[valcol]], na.rm = TRUE)
  mid <- mean(rng)
  data$label_col <- ifelse(data[[valcol]] > mid, BG_PAPER, INK_DARK)

  p <- ggplot(data, aes(GenomeA, GenomeB, fill = .data[[valcol]])) +
    geom_tile(color = GRID_LINE, linewidth = 0.6) +
    geom_text(aes(label = round(.data[[valcol]], 1), color = label_col), size = 2.1,
              family = "serif", fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradientn(colors = c(BG_PAPER, "#c9a876", "#8b5a2b", INK_DARK), name = legende) +
    coord_fixed() +
    labs(title = titre, subtitle = sous, x = NULL, y = NULL,
         caption = "Genome Taxonomy Explorer  —  Local Analysis") +
    theme_vintage
  ggsave(fichier, p, width = 9.5, height = 8.5, dpi = 300, bg = BG_PAPER)
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

# ---------- 3) POCP ----------
mat <- read.table(file.path(RES, "pocp_out/pocp_matrice.tsv"), sep = "\t",
                  header = TRUE, row.names = 1, check.names = FALSE)
d <- melt(as.matrix(mat)); colnames(d) <- c("GenomeA", "GenomeB", "POCP")
heatmap_fig(d, "POCP", "Percentage of Conserved Proteins (POCP)",
            "Genus threshold = 50%", "POCP (%)", file.path(RES, "heatmap_pocp.png"))

# ---------- 4) TREE (colored by genus) ----------
tree <- read.tree(file.path(RES, "arbre.tree"))
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

cat("All figures generated.\n")
