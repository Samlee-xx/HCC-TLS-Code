###### CosMx CCL5 and CXCL1/6/8 expression across cell types ######

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})


input_rds <- "J:/Y25-0292_napari/Y25_0292_20260812.RDS"
annotation_file <- paste0(
  "F:/TLS Material/06-Code/result/Figure5_CosMx/",
  "01_celltype_definition/CosMx_celltype_annotation.csv.gz"
)
result_dir <- paste0(
  "F:/TLS Material/06-Code/result/",
  "Figure5_CosMx/03_chemokine_expression"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

celltype_order <- c(
  "T cell", "B cell", "GC-B", "Plasma cell", "Macrophage",
  "Neutrophil", "Fibroblast", "Endothelial", "HSC", "Mural cells",
  "Cholangiocyte", "Hepatocyte", "Tumor", "Unassigned"
)
genes <- c("CCL5", "CXCL1", "CXCL6", "CXCL8")



obj <- readRDS(input_rds)
meta <- obj[[]]
meta$cell_barcode <- rownames(meta)
counts <- LayerData(obj[["RNA"]], layer = "counts")

annotation <- read.csv(gzfile(annotation_file), check.names = FALSE)
index <- match(colnames(counts), annotation$cell_barcode)
stopifnot(!anyNA(index))
celltype <- annotation$broad_celltype[index]

stopifnot(all(genes %in% rownames(counts)))
counts <- counts[genes, , drop = FALSE]

summaries <- list()

for (group in celltype_order) {
  cell_index <- which(celltype == group)
  if (length(cell_index) == 0L) next

  one <- counts[, cell_index, drop = FALSE]
  total_rna <- as.numeric(meta$nCount_RNA[cell_index])
  total_rna[!is.finite(total_rna) | total_rna <= 0] <- 1

  percent_expressing <- 100 * Matrix::rowMeans(one > 0)
  one@x <- log1p(one@x * rep(10000 / total_rna, diff(one@p)))

  summaries[[group]] <- data.frame(
    celltype = group,
    gene = genes,
    n_cells = length(cell_index),
    percent_expressing = as.numeric(percent_expressing),
    mean_log1p_10k = as.numeric(Matrix::rowMeans(one))
  )
}

expression <- do.call(rbind, summaries)
expression$expression_z <- ave(
  expression$mean_log1p_10k,
  expression$gene,
  FUN = function(x) {
    if (!is.finite(sd(x)) || sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x))
  }
)


ccl5 <- expression[
  expression$gene == "CCL5",
  c("celltype", "expression_z")
]
names(ccl5)[2] <- "score"
ccl5$signature <- "CCL5"

cxcl <- aggregate(
  expression_z ~ celltype,
  data = expression[expression$gene %in% c("CXCL1", "CXCL6", "CXCL8"), ],
  FUN = mean
)
names(cxcl)[2] <- "score"
cxcl$signature <- "CXCL1/6/8"

scores <- rbind(ccl5, cxcl)


expression$celltype <- factor(expression$celltype, levels = rev(celltype_order))
expression$gene <- factor(expression$gene, levels = genes)

p1 <- ggplot(
  expression,
  aes(x = gene, y = celltype, size = percent_expressing, color = expression_z)
) +
  geom_point() +
  scale_color_gradient2(
    low = "#2C6AA6", mid = "white", high = "#BF3F3A", midpoint = 0
  ) +
  labs(
    x = NULL, y = NULL,
    size = "Cells expressing (%)",
    color = "Mean expression\nz score"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(result_dir, "CosMx_CCL5_CXCL168_expression_dotplot.pdf"),
  p1, width = 7, height = 7
)

scores$celltype <- factor(scores$celltype, levels = rev(celltype_order))

p2 <- ggplot(scores, aes(x = score, y = celltype)) +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_col(fill = "#9EA7AD", width = 0.72) +
  facet_grid(. ~ signature, scales = "free_x") +
  labs(x = "Mean gene-wise expression z score", y = NULL) +
  theme_classic(base_size = 12)

ggsave(
  file.path(result_dir, "CosMx_CCL5_CXCL168_expression_ranking.pdf"),
  p2, width = 8, height = 6
)
