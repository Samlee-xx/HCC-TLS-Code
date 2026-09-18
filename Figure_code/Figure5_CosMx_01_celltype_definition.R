###### CosMx cell-type definition ######

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})


###### files ######

input_rds <- "J:/Y25-0292_napari/Y25_0292_20260812.RDS"
neutrophil_file <- paste0(
  "F:/TLS DATA/07-Cosmx/02-analysis-2/04-updated_celltype_dotplot/",
  "data/01_neutrophil_candidate_cells.csv.gz"
)

result_dir <- paste0(
  "F:/TLS Material/06-Code/result/",
  "Figure5_CosMx/01_celltype_definition"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)


###### canonical markers ######

marker_sets <- list(
  `T cell` = c("CD3E", "CD3D", "CCL5", "NKG7", "IL7R"),
  `B cell` = c("MS4A1", "CD79A", "CD19", "CD74"),
  `GC-B` = c("STMN1", "MKI67", "PCLAF"),
  `Plasma cell` = c("JCHAIN", "IGHG1", "MZB1", "XBP1", "IGHA1"),
  Macrophage = c("CD163", "CD68", "C1QA", "C1QB", "C1QC"),
  Neutrophil = c("MPO", "S100A8", "S100A9", "CSF3R"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Endothelial = c("ACKR1", "VWF", "PLVAP", "PECAM1", "SELE", "SELP"),
  HSC = c("RELN", "GDF2", "RBP1", "LRAT", "COLEC10", "ANGPTL6"),
  `Mural cells` = c("MYH11", "ACTA2", "CNN1", "TAGLN", "RGS5", "PDGFRB"),
  Cholangiocyte = c("KRT7", "MUC6", "MMP7", "SPP1", "CFTR"),
  Hepatocyte = c("TTR", "CYP3A4", "ALB", "CYP2C8", "HAMP"),
  Tumor = c("EPCAM", "KRT8", "KRT18", "KRT19")
)

broad_celltype_order <- c(
  "T cell", "B cell", "GC-B", "Plasma cell", "Macrophage",
  "Neutrophil", "Fibroblast", "Endothelial", "HSC", "Mural cells",
  "Cholangiocyte", "Hepatocyte", "Tumor", "Unassigned"
)


###### apply the final cell-type definitions ######

obj <- readRDS(input_rds)
meta <- obj[[]]
meta$cell_barcode <- rownames(meta)
counts <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(identical(meta$cell_barcode, colnames(counts)))

original_celltype <- as.character(meta$celltype_20260812)
updated_celltype <- original_celltype
updated_celltype[
  updated_celltype %in% c("unassgin", "unassign", "Unassigned")
] <- "Unassigned"

# Leiden cluster 8 is not used as the neutrophil definition.
updated_celltype[as.character(meta$leiden) == "8"] <- "Unassigned"

# Neutrophils are defined by the final marker-supported cell-barcode table.
neutrophil_cells <- read.csv(gzfile(neutrophil_file), check.names = FALSE)
neutrophil_index <- match(neutrophil_cells$cell_barcode, meta$cell_barcode)
stopifnot(!anyNA(neutrophil_index), !anyDuplicated(neutrophil_index))
updated_celltype[neutrophil_index] <- "Neutrophil"


###### broad labels used in marker and expression summaries ######

broad_celltype <- updated_celltype
broad_celltype[
  broad_celltype %in% c("Fibroblast", "FRC-like CAF", "FDC-like", "FDC-like cell")
] <- "Fibroblast"
broad_celltype[
  broad_celltype %in% c(
    "Endothelial", "Endothelial cell", "LSEC",
    "lymphatic endothelial", "lymphatic endothelial cell"
  )
] <- "Endothelial"
broad_celltype[
  broad_celltype %in% c("mural cells", "Mural cell", "Mural cells")
] <- "Mural cells"
broad_celltype[
  tolower(broad_celltype) %in% c("unassigned", "unassign", "unassgin")
] <- "Unassigned"


###### section-specific neutrophil diagnostic QC ######

neutrophil_qc <- rep(NA, nrow(meta))

for (section in unique(meta$sample)) {
  index <- which(meta$sample == section & updated_celltype == "Neutrophil")
  if (length(index) == 0L) next

  count_cutoff <- quantile(meta$nCount_RNA[index], 0.10, na.rm = TRUE)
  area_cutoff <- quantile(meta$Area.um2[index], c(0.01, 0.99), na.rm = TRUE)

  neutrophil_qc[index] <- (
    meta$nCount_RNA[index] >= count_cutoff &
    meta$Area.um2[index] >= area_cutoff[1] &
    meta$Area.um2[index] <= area_cutoff[2]
  )
}


###### save the cell-level annotation ######

annotation <- data.frame(
  cell_barcode = meta$cell_barcode,
  sample = meta$sample,
  fov = meta$fov,
  leiden = as.character(meta$leiden),
  original_celltype = original_celltype,
  updated_celltype = updated_celltype,
  broad_celltype = broad_celltype,
  neutrophil_qc = neutrophil_qc
)

annotation_connection <- gzfile(
  file.path(result_dir, "CosMx_celltype_annotation.csv.gz"),
  open = "wt"
)
write.csv(annotation, annotation_connection, row.names = FALSE)
close(annotation_connection)

celltype_counts <- as.data.frame(table(updated_celltype), stringsAsFactors = FALSE)
names(celltype_counts) <- c("celltype", "n_cells")
celltype_counts <- celltype_counts[order(celltype_counts$n_cells, decreasing = TRUE), ]
write.csv(
  celltype_counts,
  file.path(result_dir, "CosMx_celltype_counts.csv"),
  row.names = FALSE
)


###### verify cell identities with canonical markers ######

marker_manifest <- data.frame(
  gene = unlist(marker_sets, use.names = FALSE),
  marker_for = rep(names(marker_sets), lengths(marker_sets))
)
marker_manifest <- marker_manifest[!duplicated(marker_manifest$gene), ]
marker_manifest$present <- marker_manifest$gene %in% rownames(counts)
write.csv(
  marker_manifest,
  file.path(result_dir, "CosMx_canonical_marker_manifest.csv"),
  row.names = FALSE
)

genes <- marker_manifest$gene[marker_manifest$present]
marker_counts <- counts[genes, , drop = FALSE]
marker_summary <- list()

for (group in broad_celltype_order) {
  index <- which(broad_celltype == group)
  if (length(index) == 0L) next

  one <- marker_counts[, index, drop = FALSE]
  total_rna <- as.numeric(meta$nCount_RNA[index])
  total_rna[!is.finite(total_rna) | total_rna <= 0] <- 1

  percent_expressing <- 100 * Matrix::rowMeans(one > 0)
  one@x <- log1p(one@x * rep(10000 / total_rna, diff(one@p)))

  marker_summary[[group]] <- data.frame(
    celltype = group,
    gene = rownames(one),
    marker_for = marker_manifest$marker_for[match(rownames(one), marker_manifest$gene)],
    n_cells = length(index),
    percent_expressing = as.numeric(percent_expressing),
    mean_log1p_10k = as.numeric(Matrix::rowMeans(one))
  )
}

marker_summary <- do.call(rbind, marker_summary)
marker_summary$expression_z <- ave(
  marker_summary$mean_log1p_10k,
  marker_summary$gene,
  FUN = function(x) {
    if (!is.finite(sd(x)) || sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x))
  }
)

write.csv(
  marker_summary,
  file.path(result_dir, "CosMx_canonical_marker_expression.csv"),
  row.names = FALSE
)

marker_summary$celltype <- factor(
  marker_summary$celltype,
  levels = rev(broad_celltype_order)
)
marker_summary$gene <- factor(marker_summary$gene, levels = genes)

p <- ggplot(
  marker_summary,
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
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  file.path(result_dir, "CosMx_celltype_canonical_marker_dotplot.pdf"),
  p, width = 14, height = 7
)
