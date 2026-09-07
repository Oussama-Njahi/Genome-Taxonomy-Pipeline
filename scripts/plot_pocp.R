# plot_pocp.R - POCP heatmap (in English)
library(ggplot2)
library(reshape2)

# 1) READ the POCP matrix (square matrix, genomes as rows and columns)
mat <- read.table("/home/oussama/stage/results/pocp_out/pocp_matrice.tsv",
                   sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# 2) MELT: turn the square matrix into a long table (one row per pair)
data <- melt(as.matrix(mat))
colnames(data) <- c("GenomeA", "GenomeB", "POCP")

# 3) DRAW the heatmap
p <- ggplot(data, aes(x = GenomeA, y = GenomeB, fill = POCP)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(POCP, 1)), size = 2) +
  scale_fill_gradient(low = "blue", high = "red", name = "POCP (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(face = "bold")) +
  labs(title = "Percentage of Conserved Proteins (POCP) heatmap",
       subtitle = "Genus threshold = 50%",
       x = "Genome", y = "Genome")

# 4) SAVE as PNG
ggsave("/home/oussama/stage/results/heatmap_pocp.png", p, width = 9, height = 8, dpi = 200)
cat("POCP heatmap saved in results/heatmap_pocp.png\n")
