###### CODEX TLS subtype classification and tree plot ######

library(TreeAndLeaf)
library(igraph)

set.seed(777)

input.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/02-TLS clustering/01-input/',
  'tls_cell_composition_for_clustering_final.csv'
)

result.path <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/02-TLS clustering/03-result'
)

figure.path <- file.path(result.path, '02-figures/TLS_treeAndLeaf_BF_TNC')
dir.create(figure.path, recursive = TRUE, showWarnings = FALSE)

tls <- read.csv(input.file, check.names = FALSE)

# These ten cell-type percentages were used for TLS clustering.
features <- c(
  'pct_B',
  'pct_T_combined',
  'pct_Cholangiocyte',
  'pct_DC',
  'pct_Endothelial',
  'pct_Fibroblast',
  'pct_HEV',
  'pct_Lymphatic',
  'pct_Macrophage',
  'pct_Neutrophil'
)

celltype.matrix <- as.matrix(tls[, features])
rownames(celltype.matrix) <- tls$tls_id

###### standardize cell-type percentages ######

celltype.zscore <- scale(celltype.matrix)
celltype.zscore[is.na(celltype.zscore)] <- 0

###### cosine distance and complete-linkage clustering ######

row.norm <- sqrt(rowSums(celltype.zscore ^ 2))
celltype.unit <- celltype.zscore

celltype.unit[row.norm > 0, ] <-
  celltype.zscore[row.norm > 0, ] / row.norm[row.norm > 0]

celltype.unit[row.norm == 0, ] <- 0

cosine.distance <- 1 - tcrossprod(celltype.unit)
diag(cosine.distance) <- 0
cosine.distance[cosine.distance < 0] <- 0

tls.tree <- hclust(
  as.dist(cosine.distance),
  method = 'complete'
)

tls$cluster <- cutree(tls.tree, k = 2)

###### name the two TLS subtypes ######

mean.composition <- aggregate(
  tls[, features],
  by = list(cluster = tls$cluster),
  FUN = mean
)

bf.cluster <- mean.composition$cluster[
  which.max(mean.composition$pct_B)
]

tls$tls_subtype <- ifelse(
  tls$cluster == bf.cluster,
  'BF-TLS',
  'TNC-TLS'
)

table(tls$tls_subtype)

subtype.mean <- aggregate(
  tls[, features],
  by = list(tls_subtype = tls$tls_subtype),
  FUN = mean
)

subtype.mean

# BF-TLS has higher mean B-cell and fibroblast percentages.
# TNC-TLS has higher mean T-cell, neutrophil and cholangiocyte percentages.

###### save TLS subtype assignments ######

assignment <- tls[, c(
  'tls_id', 'slide', 'fov', 'tls_number',
  'n_cells', 'fov_region', 'cluster', 'tls_subtype'
)]

write.csv(
  assignment,
  file.path(result.path, 'tls_subtype_assignment_final.csv'),
  row.names = FALSE
)

write.csv(
  tls,
  file.path(result.path, 'tls_subtype_celltype_composition_final.csv'),
  row.names = FALSE
)

###### Figure J: TLS hierarchical tree ######

tls.graph <- treeAndLeaf(tls.tree)
leaf.node <- V(tls.graph)$isLeaf

leaf.subtype <- tls$tls_subtype[
  match(V(tls.graph)$name, tls$tls_id)
]

subtype.color <- c(
  'BF-TLS' = '#C095E4',
  'TNC-TLS' = '#72B043'
)

node.color <- rep('black', vcount(tls.graph))
node.color[leaf.node] <- subtype.color[leaf.subtype[leaf.node]]

# The only internal node with degree 2 is the root of the hclust tree.
root.node <- which(!leaf.node & degree(tls.graph) == 2)[1]
node.color[root.node] <- '#E31A1C'

set.seed(777)
tree.layout <- layout_with_fr(
  tls.graph,
  niter = 3500,
  grid = 'nogrid'
)

pdf(
  file.path(figure.path, 'TLS_treeAndLeaf_BF_TNC.pdf'),
  width = 8,
  height = 7
)

plot(
  tls.graph,
  layout = tree.layout,
  vertex.label = NA,
  vertex.size = ifelse(leaf.node, 2.2, 0.6),
  vertex.color = node.color,
  vertex.frame.color = NA,
  edge.color = 'grey60',
  edge.width = 0.45,
  margin = 0.02
)

dev.off()

