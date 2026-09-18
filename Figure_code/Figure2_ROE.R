###### Ro/e analysis of cell types across TLS spatial layers ######

library(pheatmap)

# Delaunay layers were assigned before this step.
# This script compares cell-type distributions across Edge, Transition and Core.

input.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/03-TLS Delaunay/03-result/01-tables/',
  'tls_delaunay_layer_cells_final.csv'
)

result.path <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/03-TLS Delaunay/03-result/01-tables'
)

figure.path <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/03-TLS Delaunay/03-result/02-figures/',
  'tls_layer_celltype_ROE'
)

dir.create(figure.path, recursive = TRUE, showWarnings = FALSE)

tls.layer <- read.csv(input.file, check.names = FALSE)
tls.layer <- subset(tls.layer, !is.na(delaunay_layer))

tls.layer$layer <- tls.layer$delaunay_layer
tls.layer$layer[tls.layer$layer == 'Periphery'] <- 'Edge'
tls.layer$layer <- factor(
  tls.layer$layer,
  levels = c('Edge', 'Transition', 'Core')
)

tls.layer$celltype <- tls.layer$celltype_l2

celltype.order <- c(
  'B', 'T-B', 'CD4T', 'CD8T', 'DC', 'Neutrophil', 'Macrophage',
  'Endothelial', 'HEV', 'Lymphatic', 'Fibroblast',
  'Cholangiocyte', 'Unassigned'
)

celltype.label <- c(
  'B', 'T-B', 'CD4+T', 'CD8+T', 'DC', 'Neutrophil', 'Macrophage',
  'Endothelium', 'HEV', 'LEC', 'Fibroblast',
  'Cholangiocyte', 'Unassigned'
)

tls.subtype.order <- c('TNC-TLS', 'BF-TLS')
layer.order <- c('Edge', 'Transition', 'Core')

roe.source <- data.frame()
roe.matrix <- NULL

###### calculate observed counts, expected counts and Ro/e ######

for (tls.subtype in tls.subtype.order) {

  one.subtype <- subset(tls.layer, tls_subtype == tls.subtype)

  observed.matrix <- as.matrix(
    table(one.subtype$celltype, one.subtype$layer)
  )

  observed.matrix <- observed.matrix[celltype.order, layer.order]

  roe <- observed.matrix

  for (i in 1:nrow(observed.matrix)) {
    for (j in 1:ncol(observed.matrix)) {

      observed <- observed.matrix[i, j]

      target.celltype.other.layers <-
        sum(observed.matrix[i, ]) - observed

      other.celltypes.target.layer <-
        sum(observed.matrix[, j]) - observed

      other.celltypes.other.layers <-
        sum(observed.matrix) - observed -
        target.celltype.other.layers -
        other.celltypes.target.layer

      fisher.matrix <- matrix(
        c(
          observed,
          target.celltype.other.layers,
          other.celltypes.target.layer,
          other.celltypes.other.layers
        ),
        nrow = 2,
        byrow = TRUE
      )

      expected <-
        sum(observed.matrix[i, ]) *
        sum(observed.matrix[, j]) /
        sum(observed.matrix)

      roe[i, j] <- observed / expected
      p.value <- fisher.test(fisher.matrix)$p.value

      roe.source <- rbind(
        roe.source,
        data.frame(
          tls_subtype = tls.subtype,
          celltype = celltype.order[i],
          celltype_label = celltype.label[i],
          layer = layer.order[j],
          observed = observed,
          expected = expected,
          roe = observed / expected,
          p_value = p.value
        )
      )
    }
  }

  colnames(roe) <- paste(tls.subtype, layer.order, sep = '_')
  rownames(roe) <- celltype.label
  roe.matrix <- cbind(roe.matrix, roe)
}

roe.matrix <- as.matrix(roe.matrix)
storage.mode(roe.matrix) <- 'numeric'

write.csv(
  roe.source,
  file.path(result.path, 'tls_layer_celltype_ROE_source_data_final.csv'),
  row.names = FALSE
)

write.csv(
  roe.matrix,
  file.path(result.path, 'tls_layer_celltype_ROE_matrix_final.csv')
)

###### draw Ro/e heatmap ######

# The complete Ro/e values are retained in the CSV file.

roe.plot <- roe.matrix
breaks <- c(
  seq(0.5, 1, by = 0.01),
  seq(1.001, 1.6, by = 0.01)
)

heatmap.color <- c(
  colorRampPalette(c('white', '#F7E6D7'))(50),
  colorRampPalette(c('#F7E6D7', '#FA151D'))(60)
)

pdf(
  file.path(figure.path, 'tls_layer_celltype_ROE_heatmap.pdf'),
  height = 5.8,
  width = 5.4
)

pheatmap(
  roe.plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = 3,
  gaps_row = 7,
  cellwidth = 25,
  cellheight = 20,
  color = heatmap.color,
  breaks = breaks,
  labels_col = rep(layer.order, 2),
  legend_breaks = c(0.8, 1.2, 1.6),
  main = 'TNC-TLS            BF-TLS',
  angle_col = 90
)

dev.off()

head(roe.source)
roe.matrix
