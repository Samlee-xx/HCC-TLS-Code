library(tidyverse)
library(pheatmap)
library(ggpubr)


# Plot the cholangiocyte-like cell and neutrophil spatial co-occurrence results

input_dir <- "F:/TLS Material/06-Code/result/Figure2_Chol_Neu_spatial_cooccurrence"
figure_dir <- "F:/TLS Material/06-Code/result/Figure2_Chol_Neu_spatial_cooccurrence/figure"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# 1. Radial-profile heatmap
score_z <- read.csv(
  file.path(input_dir, "chol_neu_score_row_z.csv.gz"),
  row.names = 1,
  check.names = FALSE
)

cell_info <- read.csv(
  file.path(input_dir, "chol_neu_cell_assignment.csv.gz"),
  check.names = FALSE
)

table(cell_info$cluster)
table(cell_info$chol_group)
dim(score_z)

cell_info$cluster <- factor(
  cell_info$cluster,
  levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
)
cell_info <- cell_info[order(cell_info$cluster), ]
score_z <- score_z[cell_info$cellid, ]

ann_row <- data.frame(cluster = cell_info$cluster)
rownames(ann_row) <- cell_info$cellid

cluster_color <- c(
  Cluster1 = "#B2182B",
  Cluster2 = "#EF8A62",
  Cluster3 = "#67A9CF",
  Cluster4 = "#2166AC"
)

gap_row <- cumsum(table(cell_info$cluster))
gap_row <- gap_row[-length(gap_row)]

pdf(file.path(figure_dir, "chol_neu_radial_heatmap.pdf"), width = 6, height = 8)
pheatmap(
  as.matrix(score_z),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = ann_row,
  annotation_colors = list(cluster = cluster_color),
  gaps_row = gap_row,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA
)
dev.off()


# 2. Mean raw co-occurrence score of each cluster
cluster_mean <- read.csv(
  file.path(input_dir, "chol_neu_cluster_mean_raw.csv"),
  check.names = FALSE
)

cluster_mean_long <- cluster_mean %>%
  pivot_longer(
    cols = starts_with("r"),
    names_to = "radius",
    values_to = "score"
  ) %>%
  mutate(
    radius = as.numeric(sub("r", "", radius)),
    cluster = factor(cluster, levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4"))
  )

p_curve <- ggplot(
  cluster_mean_long,
  aes(x = radius, y = score, color = cluster)
) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey60") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cluster_color) +
  scale_x_continuous(breaks = seq(50, 300, 50)) +
  labs(x = "Radius (pixels)", y = "Co-occurrence score", color = NULL) +
  theme_classic()

ggsave(
  file.path(figure_dir, "chol_neu_cluster_curve.pdf"),
  p_curve,
  width = 5,
  height = 4
)


# 3. TLS-level NeuProx-Chol percentage
tls_result <- read.csv(
  file.path(input_dir, "chol_neu_TLS_percentage.csv"),
  check.names = FALSE
)

table(tls_result$tls_subtype)
summary(tls_result$n_clusterable_chol)

tls_result$tls_subtype <- factor(
  tls_result$tls_subtype,
  levels = c("TNC-TLS", "BF-TLS")
)

p_box <- ggplot(
  tls_result,
  aes(x = tls_subtype, y = neu_prox_percentage, fill = tls_subtype)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.65) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("TNC-TLS" = "#74ADD1", "BF-TLS" = "#FDAE61")) +
  labs(x = NULL, y = "NeuProx-Chols (%)") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "NeuProx_Chol_percentage.pdf"),
  p_box,
  width = 3.6,
  height = 4
)

