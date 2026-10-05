# plot_all_web.R - Figures of the pipeline (heatmaps + tree), used both by the
# web interface and the command line (pipeline.py runs pipeline_web.py).
# The results directory comes from argv[1], so each web job renders into its
# own folder.
suppressMessages({
  library(ggplot2); library(reshape2); library(ggtree)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_all_web.R <results_dir>")
}
RES <- args[1]

# theme "planche d'atlas" -- les memes encres que l'interface web (style
# Haeckel, Kunstformen der Natur) : papier de planche, encre brune, bleu de
# Prusse du cote "meme taxon", terre d'ombre de l'autre, carmin pour le seuil.
# Les couleurs doivent rester identiques a webapp/frontend/js/plates.js.
BG_PAPER   <- "#F6EEE3"   # papier de planche
INK_DARK   <- "#2A211B"   # encre
INK_MED    <- "#5C4A3D"   # encre secondaire
GRID_LINE  <- "#A8957D"   # filet fin
ACCENT_RED <- "#9B2F3A"   # carmin (seuils)
NA_FILL    <- "#E6DCCB"   # case sans valeur (NA), hors de l'echelle
SAME_RAMP  <- c("#DCE2E4", "#93A6B6", "#4F6A82", "#2F4A63")   # du seuil vers "meme taxon"
DIFF_RAMP  <- c("#EFE1CA", "#D2B48D", "#A47E57", "#7A5636")   # du seuil vers "autre taxon"

theme_plate <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(plot.background    = element_rect(fill = BG_PAPER, color = INK_MED, linewidth = 0.8),
        panel.background   = element_rect(fill = BG_PAPER, color = NA),
        panel.border       = element_blank(),
        plot.title         = element_text(size = 16, color = INK_DARK, family = "serif", hjust = 0),
        plot.subtitle      = element_text(size = 9.5, color = INK_MED, family = "serif",
                                           hjust = 0, face = "italic", margin = margin(b = 10)),
        plot.caption       = element_text(size = 7, color = INK_MED, family = "serif", hjust = 0),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 7, color = INK_DARK,
                                           family = "serif", face = "italic"),
        axis.text.y        = element_text(size = 7, color = INK_DARK, family = "serif", face = "italic"),
        panel.grid         = element_blank(),
        legend.background  = element_rect(fill = BG_PAPER, color = NA),
        legend.key         = element_rect(fill = BG_PAPER, color = NA),
        legend.text        = element_text(color = INK_MED, family = "serif", size = 7),
        legend.title       = element_text(color = INK_DARK, family = "serif", size = 9),
        legend.position    = "right",
        plot.margin        = margin(18, 20, 14, 16))

clean <- function(x) gsub("\\.(fas|fasta|fna)$", "", gsub(".*/", "", x))
# "Bacillus-subtilis-C3-41" -> "Bacillus subtilis C3-41" (le tiret de la souche reste)
affiche_nom <- function(x) sub("^([^-_ ]+)[-_ ]([^-_ ]+)[-_ ]?", "\\1 \\2 ", x)

# Ordre des genomes : celui de l'arbre s'il existe (les clades se lisent en
# blocs, comme dans l'interface), sinon un regroupement UPGMA sur la mesure.
tree_file <- file.path(RES, "arbre.tree")
tree_order <- if (file.exists(tree_file)) read.tree(tree_file)$tip.label else NULL

genome_order <- function(data, valcol, distance) {
  noms <- sort(unique(c(as.character(data$GenomeA), as.character(data$GenomeB))))
  if (!is.null(tree_order) && setequal(noms, tree_order)) return(tree_order)
  m <- matrix(NA_real_, length(noms), length(noms), dimnames = list(noms, noms))
  m[cbind(as.character(data$GenomeA), as.character(data$GenomeB))] <- data[[valcol]]
  d <- if (distance) m else 100 - m
  d <- (d + t(d)) / 2                      # ANI est asymetrique : moyenne des deux sens
  d[is.na(d)] <- max(d, na.rm = TRUE)      # paire sans valeur = la plus eloignee
  diag(d) <- 0
  noms[hclust(as.dist(d), method = "average")$order]
}

# Couleur d'une valeur : deux encres qui divergent exactement au seuil, plus
# claires pres du seuil, plus sombres loin de lui (meme regle que plates.js).
couleur <- function(v, seuil, lo, hi, distance) {
  rampe <- function(cols, t) rgb(colorRamp(cols)(pmin(1, pmax(0, t))^0.8), maxColorValue = 255)
  out <- rep(NA_character_, length(v))
  ok <- !is.na(v)
  meme <- ok & (if (distance) v <= seuil else v >= seuil)
  autre <- ok & !meme
  if (distance) {
    out[meme]  <- rampe(SAME_RAMP, (seuil - v[meme]) / (seuil - lo))
    out[autre] <- rampe(DIFF_RAMP, (v[autre] - seuil) / (hi - seuil))
  } else {
    out[meme]  <- rampe(SAME_RAMP, (v[meme] - seuil) / (hi - seuil))
    out[autre] <- rampe(DIFF_RAMP, (seuil - v[autre]) / (seuil - lo))
  }
  out
}
sombre <- function(cols) {
  rgbv <- col2rgb(cols) / 255
  lin <- ifelse(rgbv <= 0.03928, rgbv / 12.92, ((rgbv + 0.055) / 1.055)^2.4)
  (0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]) < 0.3
}

# fonction generique : dessine une heatmap "planche d'atlas"
# distance = TRUE pour une DISTANCE (valeur BASSE = genomes PROCHES, donc du
# cote "meme taxon" sous le seuil) au lieu d'un pourcentage de similarite.
heatmap_fig <- function(data, valcol, titre, sous, legende, fichier, seuil,
                        decimales = 1, distance = FALSE) {
  vals <- data[[valcol]][as.character(data$GenomeA) != as.character(data$GenomeB)]
  if (distance) { lo <- 0; hi <- max(c(vals, seuil * 1.5), na.rm = TRUE) }
  else { lo <- min(c(vals, seuil - 1), na.rm = TRUE); hi <- 100 }

  ordre <- genome_order(data, valcol, distance)
  data$GenomeA <- factor(data$GenomeA, levels = rev(ordre))   # premier genome en haut
  data$GenomeB <- factor(data$GenomeB, levels = ordre)

  # Une valeur manquante (ex. paire dDDH sans alignement significatif) est
  # affichee "NA" sur un fond neutre, jamais confondue avec une vraie valeur.
  manquant <- is.na(data[[valcol]])
  fond <- couleur(data[[valcol]], seuil, lo, hi, distance)
  meme <- !manquant & (if (distance) data[[valcol]] <= seuil else data[[valcol]] >= seuil)
  data$label_txt <- ifelse(manquant, "NA", formatC(data[[valcol]], format = "f", digits = decimales))
  data$label_col <- ifelse(manquant, INK_MED, ifelse(!manquant & sombre(ifelse(manquant, NA_FILL, fond)),
                                                     BG_PAPER, INK_DARK))
  data$label_face <- ifelse(meme, "bold", "plain")             # chiffres gras = meme taxon

  # la legende reprend les deux rampes, jointes au seuil
  # Echelle divergente : le seuil au milieu de la legende, chaque cote sur sa
  # propre plage (sinon le cote "meme espece" du dDDH n'est qu'un trait).
  cols <- if (distance) c(rev(SAME_RAMP), DIFF_RAMP) else c(rev(DIFF_RAMP), SAME_RAMP)
  data$t <- ifelse(data[[valcol]] <= seuil, 0.5 * (data[[valcol]] - lo) / (seuil - lo),
                   0.5 + 0.5 * (data[[valcol]] - seuil) / (hi - seuil))
  pos <- c(seq(0, 0.5, length.out = 4), seq(0.5001, 1, length.out = 4))

  p <- ggplot(data, aes(GenomeB, GenomeA, fill = t)) +
    geom_tile(color = BG_PAPER, linewidth = 0.8) +
    geom_text(aes(label = label_txt, color = label_col, fontface = label_face), size = 2.1,
              family = "serif") +
    scale_color_identity() +
    scale_fill_gradientn(colors = cols, values = pos, limits = c(0, 1), name = legende, na.value = NA_FILL,
                         breaks = c(0, 0.5, 1), labels = formatC(c(lo, seuil, hi), format = "f", digits = decimales),
                         guide = guide_colourbar(barheight = grid::unit(5, "cm"), ticks.colour = ACCENT_RED)) +
    scale_x_discrete(labels = affiche_nom) +
    scale_y_discrete(labels = affiche_nom) +
    coord_fixed() +
    labs(title = titre, subtitle = sous, x = NULL, y = NULL,
         caption = "Taxo Genomics Analyses  ·  bold figures: same taxon  ·  NA: no value") +
    theme_plate
  ggsave(fichier, p, width = 9.5, height = 8.5, dpi = 300, bg = BG_PAPER)
  cat("  ->", fichier, "\n")
}

# Chaque etape du pipeline est optionnelle : on ne dessine une figure que si
# l'etape qui produit son fichier d'entree a tourne, sinon on la saute
# (au lieu de faire echouer toute l'etape Figures sur un fichier absent).
n_figures <- 0
present <- function(chemin, figure) {
  if (file.exists(chemin)) return(TRUE)
  cat("  -- skipped", figure, "(", basename(chemin), "not found: step not run )\n")
  FALSE
}

# ---------- 1) ANI ----------
f <- file.path(RES, "ani_resultats.txt")
if (present(f, "ANI heatmap")) {
  d <- read.table(f, sep = "\t", header = FALSE)[, 1:3]
  colnames(d) <- c("GenomeA", "GenomeB", "ANI")
  d$GenomeA <- clean(d$GenomeA); d$GenomeB <- clean(d$GenomeB)
  # fastANI n'ecrit pas les paires trop eloignees : on les ajoute en NA
  noms <- unique(c(d$GenomeA, d$GenomeB))
  d <- merge(expand.grid(GenomeA = noms, GenomeB = noms, stringsAsFactors = FALSE), d, all.x = TRUE)
  heatmap_fig(d, "ANI", "Average Nucleotide Identity (ANI)",
              "Species threshold 95 %. Rows: query genome. NA: too distant for FastANI.",
              "ANI (%)", file.path(RES, "heatmap_ani.png"), seuil = 95)
  n_figures <- n_figures + 1
}

# ---------- 2) AAI ----------
f <- file.path(RES, "aai_resultats.tsv")
if (present(f, "AAI heatmap")) {
  raw <- read.table(f, sep = "\t", header = TRUE, check.names = FALSE)
  d <- raw[, c(3, 4, 5)]; colnames(d) <- c("GenomeA", "GenomeB", "AAI")
  heatmap_fig(d, "AAI", "Average Amino acid Identity (AAI)",
              "Approximate genus boundary around 65 %", "AAI (%)", file.path(RES, "heatmap_aai.png"),
              seuil = 65)
  n_figures <- n_figures + 1
}

# ---------- 3) dDDH (distance GBDP, formule 2) ----------
# On affiche la DISTANCE (0 = identique), pas le %similarite estime :
# ce dernier n'est fiable que pres du seuil espece et devient absurde
# (parfois negatif) pour des genomes deja tres divergents -- voir dddh.py.
f <- file.path(RES, "dddh_out/dddh_distance.tsv")
if (present(f, "dDDH heatmap")) {
  mat <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  d <- melt(as.matrix(mat)); colnames(d) <- c("GenomeA", "GenomeB", "dDDH")
  d$GenomeA <- as.character(d$GenomeA); d$GenomeB <- as.character(d$GenomeB)
  heatmap_fig(d, "dDDH", "digital DNA-DNA Hybridization (dDDH, GBDP formula 2)",
              "Auch et al. 2010 method · species threshold: distance > 0.0412 (~<70% DDH) = distinct species",
              "dDDH distance", file.path(RES, "heatmap_dddh.png"), seuil = 0.0412,
              decimales = 3, distance = TRUE)
  n_figures <- n_figures + 1
}

# ---------- 4) POCP ----------
f <- file.path(RES, "pocp_out/pocp_matrice.tsv")
if (present(f, "POCP heatmap")) {
  mat <- read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  d <- melt(as.matrix(mat)); colnames(d) <- c("GenomeA", "GenomeB", "POCP")
  d$GenomeA <- as.character(d$GenomeA); d$GenomeB <- as.character(d$GenomeB)
  heatmap_fig(d, "POCP", "Percentage of Conserved Proteins (POCP)",
              "Genus threshold 50 %", "POCP (%)", file.path(RES, "heatmap_pocp.png"), seuil = 50)
  n_figures <- n_figures + 1
}

# ---------- 5) TREE (engraved: ink branches, genus brackets) ----------
f <- tree_file
if (present(f, "tree figure")) {
  tree <- read.tree(f)
  # ggtree dessine la premiere feuille en bas : on fait pivoter l'arbre pour
  # qu'elle soit en haut, dans le meme ordre que l'interface et les heatmaps
  tree <- rotateConstr(tree, rev(tree$tip.label))
  genus <- sub("[-_].*", "", tree$tip.label)   # le genre = 1er mot du nom
  dd <- data.frame(label = tree$tip.label, Genus = genus,
                   affiche = affiche_nom(tree$tip.label))
  # ladderize = FALSE : meme ordre des feuilles que l'interface et les heatmaps
  pt <- ggtree(tree, ladderize = FALSE, color = INK_DARK, linewidth = 0.45) %<+% dd +
    geom_tiplab(aes(label = affiche), size = 3.1, family = "serif", fontface = "italic",
                color = INK_DARK, align = TRUE, linetype = "dotted", linesize = 0.25,
                offset = 0.005)
  # une accolade par genre, sur la marge droite (ordre des feuilles de l'arbre)
  pos <- pt$data[pt$data$isTip, c("label", "y")]
  pos <- pos[order(-pos$y), ]
  pos$Genus <- sub("[-_].*", "", pos$label)
  blocs <- rle(pos$Genus)
  fin <- cumsum(blocs$lengths); debut <- fin - blocs$lengths + 1
  for (k in seq_along(blocs$values)) {
    pt <- pt + geom_strip(pos$label[debut[k]], pos$label[fin[k]], label = blocs$values[k],
                          offset = max(pt$data$x) * 0.42, offset.text = max(pt$data$x) * 0.02,
                          barsize = 0.4, color = INK_MED, fontsize = 3, family = "serif",
                          extend = 0.25)
  }
  pt <- pt + geom_treescale(fontsize = 2.6, linesize = 0.4, color = INK_MED, family = "serif") +
    xlim(0, max(pt$data$x) * 1.75) +
    labs(title = "Phylogenomic tree (120 core genes, bac120)",
         caption = "Branch lengths in substitutions per site  ·  Taxo Genomics Analyses") +
    theme_tree() +
    theme(plot.background = element_rect(fill = BG_PAPER, color = INK_MED, linewidth = 0.8),
          panel.background = element_rect(fill = BG_PAPER, color = NA),
          plot.title = element_text(size = 15, color = INK_DARK, family = "serif"),
          plot.caption = element_text(size = 7, color = INK_MED, family = "serif", hjust = 0),
          plot.margin = margin(18, 20, 14, 16))
  ggsave(file.path(RES, "tree.png"), pt, width = 12, height = 7.5, dpi = 300, bg = BG_PAPER)
  cat("  -> tree.png\n")
  n_figures <- n_figures + 1
}

if (n_figures == 0) {
  cat("No figure generated: none of ANI, dDDH, AAI, POCP or the tree was run.\n")
} else {
  cat(n_figures, "figure(s) generated.\n")
}
