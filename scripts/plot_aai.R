# plot_aai.R - AAI heatmap (in English)
library(ggplot2)
library(reshape2)

# 1) READ the AAI matrix produced by EzAAI
#    EzAAI output has several columns; we read it and build a matrix.
raw <- read.table("/home/oussama/stage/results/aai_resultats.tsv",
                   sep = "\t", header = TRUE, check.names = FALSE)

# EzAAI 'calculate' output columns: ID 1, ID 2, Label 1, Label 2, AAI, ...
# We keep the two labels and the AAI value (columns 3, 4, 5).
data <- raw[, c(3, 4, 5)]
colnames(data) <- c("GenomeA", "GenomeB", "AAI")

# 2) DRAW the heatmap
p <- ggplot(data, aes(x = GenomeA, y = GenomeB, fill = AAI)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(AAI, 1)), size = 2) +
  scale_fill_gradient(low = "blue", high = "red", name = "AAI (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(face = "bold")) +
  labs(title = "Average Amino acid Identity (AAI) heatmap",
       subtitle = "Approximate genus boundary around 65%",
       x = "Genome", y = "Genome")

# 3) SAVE as PNG
ggsave("/home/oussama/stage/results/heatmap_aai.png", p, width = 9, height = 8, dpi = 200)
cat("AAI heatmap saved in results/heatmap_aai.png\n")
