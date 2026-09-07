# plot_ani.R - Premiere heatmap en R : la matrice ANI
library(ggplot2)

# 1) LIRE le fichier ANI produit par FastANI (5 colonnes, separees par tabulation)
data <- read.table("/home/oussama/stage/results/ani_resultats.txt",
                    sep = "\t", header = FALSE)

# 2) PREPARER : garder les 3 premieres colonnes et les nommer
data <- data[, 1:3]
colnames(data) <- c("GenomeA", "GenomeB", "ANI")

# nettoyer les noms (enlever le chemin et l'extension .fas)
data$GenomeA <- gsub(".*/", "", data$GenomeA)   # enleve le chemin
data$GenomeA <- gsub(".fas", "", data$GenomeA)  # enleve .fas
data$GenomeB <- gsub(".*/", "", data$GenomeB)
data$GenomeB <- gsub(".fas", "", data$GenomeB)

# 3) DESSINER la heatmap (on empile les couches avec +)
p <- ggplot(data, aes(x = GenomeA, y = GenomeB, fill = ANI)) +
  geom_tile(color = "white") +                                  # les carreaux
  geom_text(aes(label = round(ANI, 1)), size = 2) +             # la valeur dans chaque carreau
  scale_fill_gradient(low = "blue", high = "red") +             # bleu = faible, rouge = eleve
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Average Nucleotide Identity (ANI) heatmap",
       x = "", y = "", fill = "ANI (%)")

# 4) SAUVER en image PNG
ggsave("/home/oussama/stage/results/heatmap_ani.png", p, width = 9, height = 8, dpi = 200)
cat("ANI heatmap saved in results/heatmap_ani.png\n")
