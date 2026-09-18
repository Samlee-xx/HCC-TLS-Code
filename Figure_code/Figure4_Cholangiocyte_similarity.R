###### Cholangiocyte-like transcriptional similarity along the tumor axis ######

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)


###### input and output files ######

base.path <- 'F:/TLS DATA/02-HD'

reference.file <- paste0(
  base.path,
  '/03-analysis/03-BileDuct similarity/previous data and code/',
  'Merge Normal HCC bile duct.rds'
)

# This table contains the cholangiocyte annotation, sample identity and
# array_row/array_col coordinates for each 8-um Visium HD bin.
bin.file <- paste0(
  base.path,
  '/03-analysis/03-BileDuct similarity/01-new-sample-axis-qc/',
  'all_new_sample_bins_with_celltype_for_axis_qc.csv'
)

result.path <- paste0(
  base.path,
  '/03-analysis/03-BileDuct similarity/',
  '02-new-sample-cholangiocyte-similarity/submission_results'
)

dir.create(result.path, recursive = TRUE, showWarnings = FALSE)


###### read the single-cell cholangiocyte reference ######

reference <- readRDS(reference.file)
reference.expression <- reference@assays$RNA@data

hcc.cells <- reference$group == 'hcc'
normal.cells <- reference$group == 'health'

table(reference$group)


###### read the Visium HD cholangiocyte-like bins ######

bin.table <- read.csv(bin.file, check.names = FALSE)
bin.table <- subset(bin.table, cell_type == 'Cholangiocyte')

table(bin.table$sample)


###### cosine similarity ######

cosine.similarity <- function(expression, reference.profile) {

  numerator <- as.numeric(
    Matrix::crossprod(expression, reference.profile)
  )
  bin.norm <- sqrt(Matrix::colSums(expression ^ 2))
  reference.norm <- sqrt(sum(reference.profile ^ 2))

  similarity <- numerator / (bin.norm * reference.norm)
  similarity[!is.finite(similarity)] <- NA
  similarity
}


###### calculate similarity for the four Visium HD slides ######

all.similarity <- data.frame()

for (slide.id in c('LS-1', 'LS-2', 'LS-3', 'LS-4')) {

  message('Processing ', slide.id)

  matrix.file <- file.path(
    base.path,
    '01-New data',
    slide.id,
    'outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5'
  )

  hd.expression <- Read10X_h5(matrix.file)
  if (is.list(hd.expression)) {
    hd.expression <- hd.expression[[1]]
  }

  one.slide <- subset(bin.table, ls_id == slide.id)

  matched.column <- match(one.slide$Barcode, colnames(hd.expression))
  one.slide <- one.slide[!is.na(matched.column), ]
  matched.column <- matched.column[!is.na(matched.column)]

  hd.cholangiocyte <- hd.expression[, matched.column, drop = FALSE]
  colnames(hd.cholangiocyte) <- one.slide$Barcode

  # Project the HD bins and scRNA-seq reference into the same gene space.
  common.genes <- intersect(
    rownames(reference.expression),
    rownames(hd.cholangiocyte)
  )

  hd.cholangiocyte <- hd.cholangiocyte[
    common.genes,
    ,
    drop = FALSE
  ]

  hcc.average <- Matrix::rowMeans(
    reference.expression[common.genes, hcc.cells, drop = FALSE]
  )
  normal.average <- Matrix::rowMeans(
    reference.expression[common.genes, normal.cells, drop = FALSE]
  )

  similarity.hcc <- cosine.similarity(
    hd.cholangiocyte,
    hcc.average
  )
  similarity.normal <- cosine.similarity(
    hd.cholangiocyte,
    normal.average
  )

  one.result <- data.frame(
    barcode = one.slide$Barcode,
    sample = one.slide$sample,
    slide = one.slide$slide,
    ls_id = one.slide$ls_id,
    array_row = one.slide$array_row,
    array_col = one.slide$array_col,
    similarity_hcc = similarity.hcc,
    similarity_normal = similarity.normal
  )

  all.similarity <- rbind(all.similarity, one.result)
}


###### standardize similarity within each specimen ######

all.similarity <- all.similarity %>%
  group_by(sample) %>%
  mutate(
    similarity_hcc_z = as.numeric(scale(similarity_hcc)),
    similarity_normal_z = as.numeric(scale(similarity_normal)),
    delta_similarity = similarity_hcc_z - similarity_normal_z
  ) %>%
  ungroup() %>%
  filter(
    is.finite(delta_similarity),
    is.finite(array_row),
    is.finite(array_col)
  )

write.csv(
  all.similarity,
  file.path(result.path, 'cholangiocyte_bin_cosine_similarity.csv'),
  row.names = FALSE
)


###### select the spatial axis oriented toward the tumor ######

# The appropriate local axis was chosen after inspecting the spatial
# distributions of tumor and cholangiocyte-like bins in each specimen.
select.axis <- function(data, sample.id) {

  data <- subset(data, sample == sample.id)

  if (sample.id == 'S1-1') {
    data <- subset(data, array_col >= 700 & array_col <= 1000)
    data$axis_value <- data$array_row
    data$axis_name <- 'y = array_row, with x 700-1000'
  }

  if (sample.id == 'S1-4') {
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col'
  }

  if (sample.id == 'S2-1') {
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col'
  }

  if (sample.id == 'S2-3') {
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col'
  }

  if (sample.id == 'S2-4') {
    data <- subset(data, array_row >= 1000 & array_row <= 1400)
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col, with y 1000-1400'
  }

  if (sample.id == 'S3-1') {
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col'
  }

  if (sample.id == 'S3-2') {
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col'
  }

  if (sample.id == 'S3-3') {
    data <- subset(data, array_col >= 1000 & array_col <= 1400)
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col, with x 1000-1400'
  }

  if (sample.id == 'S3-4') {
    data <- subset(data, array_row >= 900 & array_row <= 1400)
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col, with y 900-1400'
  }

  if (sample.id == 'S4-1') {
    data <- subset(data, array_col >= 900 & array_col <= 1200)
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col, with x 900-1200'
  }

  if (sample.id == 'S4-2') {
    x1 <- 700
    y1 <- 800
    x2 <- 100
    y2 <- 200

    dx <- x2 - x1
    dy <- y2 - y1
    line.length <- sqrt(dx ^ 2 + dy ^ 2)

    data$axis_value <- (
      (data$array_col - x1) * dx +
        (data$array_row - y1) * dy
    ) / line.length

    data <- subset(
      data,
      axis_value >= 0 & axis_value <= line.length
    )
    data$axis_name <- 'diagonal: (700,800) to (100,200)'
  }

  if (sample.id == 'S4-3') {
    data <- subset(data, array_col >= 1000 & array_col <= 1250)
    data$axis_value <- data$array_col
    data$axis_name <- 'x = array_col, with x 1000-1250'
  }

  if (sample.id == 'S4-4') {
    x1 <- 100
    y1 <- 800
    x2 <- 700
    y2 <- 1400

    dx <- x2 - x1
    dy <- y2 - y1
    line.length <- sqrt(dx ^ 2 + dy ^ 2)

    data$axis_value <- (
      (data$array_col - x1) * dx +
        (data$array_row - y1) * dy
    ) / line.length

    data <- subset(
      data,
      axis_value >= 0 & axis_value <= line.length
    )
    data$axis_name <- 'diagonal: (100,800) to (700,1400)'
  }

  data
}


###### prepare the selected-axis data ######

# S1-2 and S1-3 were not used for the final local-axis curves due to no tumor boundary.
# S2-2 contained only 42 cholangiocyte-like bins and was excluded.
selected.samples <- c(
  'S1-1', 'S1-4',
  'S2-1', 'S2-3', 'S2-4',
  'S3-1', 'S3-2', 'S3-3', 'S3-4',
  'S4-1', 'S4-2', 'S4-3', 'S4-4'
)

selected.data <- bind_rows(
  lapply(
    selected.samples,
    function(x) select.axis(all.similarity, x)
  )
)

selected.data <- subset(
  selected.data,
  is.finite(axis_value) & is.finite(delta_similarity)
)

write.csv(
  selected.data,
  file.path(result.path, 'selected_axis_cholangiocyte_similarity.csv'),
  row.names = FALSE
)


###### correlation between spatial position and delta similarity ######

correlation.result <- selected.data %>%
  group_by(sample, axis_name) %>%
  summarise(
    cholangiocyte_bins = n(),
    axis_min = min(axis_value),
    axis_max = max(axis_value),
    spearman_rho = cor(
      axis_value,
      delta_similarity,
      method = 'spearman'
    ),
    spearman_p = cor.test(
      axis_value,
      delta_similarity,
      method = 'spearman',
      exact = FALSE
    )$p.value,
    .groups = 'drop'
  )

write.csv(
  correlation.result,
  file.path(result.path, 'selected_axis_similarity_correlation.csv'),
  row.names = FALSE
)


plot.list <- list()
for (sample.id in selected.samples) {

  one.sample <- subset(selected.data, sample == sample.id)

  p <- ggplot(
    one.sample,
    aes(x = axis_value, y = delta_similarity)
  ) +
    geom_smooth(
      method = 'loess',
      se = TRUE,
      color = 'darkred',
      fill = 'grey75',
      linewidth = 1.1
    ) +
    theme_classic() +
    labs(
      title = sample.id,
      x = unique(one.sample$axis_name),
      y = expression(Delta~'cosine similarity score')
    )

  plot.list[[sample.id]] <- p

  ggsave(
    file.path(
      result.path,
      paste0(sample.id, '_cholangiocyte_similarity.pdf')
    ),
    p,
    width = 3.5,
    height = 3
  )
}

overview <- wrap_plots(plot.list, ncol = 4)

ggsave(
  file.path(result.path, 'all_cholangiocyte_similarity_curves.pdf'),
  overview,
  width = 14,
  height = 10
)

print(correlation.result)
