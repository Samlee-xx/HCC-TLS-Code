###### CODEX and Visium HD integration plots ######

library(pheatmap)


###### files ######

input_dir <- paste0(
  "F:/TLS Material/06-Code/result/",
  "Figure4_CODEX_HD_integration"
)

matching_file <- file.path(
  input_dir,
  "CODEX_HD_candidate_TLS_pairs.csv"
)
correlation_file <- file.path(
  input_dir,
  "TNC_TLS_CODEX_HD_Pearson_r_matrix.csv"
)
fdr_file <- file.path(
  input_dir,
  "TNC_TLS_CODEX_HD_BH_FDR_matrix.csv"
)


###### TLS-pair composition concordance ######

matching <- read.csv(matching_file, check.names = FALSE)
matching <- matching[order(matching$pearson_r), ]

bar_color <- ifelse(matching$matched, "#C84E57", "#BDBDBD")

pdf(
  file.path(input_dir, "CODEX_HD_TLS_pair_Pearson_matching.pdf"),
  width = 8,
  height = 4
)

barplot(
  matching$pearson_r,
  col = bar_color,
  border = NA,
  space = 0,
  xlab = "Candidate TLS pairs",
  ylab = "Pearson r"
)
abline(h = 0, col = "black", lwd = 0.8)
box(bty = "l")

legend(
  "topleft",
  legend = c("Matched", "Not matched"),
  fill = c("#C84E57", "#BDBDBD"),
  border = NA,
  bty = "n"
)

dev.off()


###### correlations in the 30 matched TNC-TLSs ######

pearson_r <- as.matrix(
  read.csv(correlation_file, row.names = 1, check.names = FALSE)
)
bh_fdr <- as.matrix(
  read.csv(fdr_file, row.names = 1, check.names = FALSE)
)

stars <- matrix(
  "",
  nrow = nrow(pearson_r),
  ncol = ncol(pearson_r),
  dimnames = dimnames(pearson_r)
)
stars[bh_fdr < 0.05] <- "*"

heatmap_color <- colorRampPalette(
  c("#739559", "white", "#F50538")
)(100)

pdf(
  file.path(input_dir, "TNC_TLS_CODEX_HD_function_Pearson_heatmap.pdf"),
  width = 9,
  height = 6.5
)

pheatmap(
  pearson_r,
  color = heatmap_color,
  breaks = seq(-0.5, 0.5, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  display_numbers = stars,
  number_color = "black",
  fontsize_number = 14,
  angle_col = 45,
  border_color = "grey80",
  main = "Matched TNC-TLSs (n = 30); * BH-FDR < 0.05"
)

dev.off()


