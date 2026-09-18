###### Copy-number inference using CopyKAT ######

library(data.table)
library(Matrix)
library(Seurat)
library(copykat)
library(ggplot2)

set.seed(20260708)


###### paths ######

base.path <- 'F:/TLS DATA/02-HD'
annotation.path <- file.path(
  base.path,
  '03-analysis/02-Slide and cluster barcode'
)
region.path <- file.path(
  base.path,
  '03-analysis/04-Cholangiocyte copykat/Tumor region barcodes'
)
result.path <- file.path(
  base.path,
  '03-analysis/04-Cholangiocyte copykat/submission_results'
)

dir.create(result.path, recursive = TRUE, showWarnings = FALSE)


###### analysis settings ######

minimum.detected.genes <- 200
minimum.normal.references <- 20

# These are the final CopyKAT settings used for the sparse Visium HD data.
copykat.low.dr <- 0
copykat.up.dr <- 0.01
copykat.ks.cut <- 0.01

# Visium HD contains many bins. The same fixed limits are applied to every
# specimen before CopyKAT to keep the calculation manageable.
maximum.bile.duct.bins <- 800
maximum.tumor.bins <- 700
maximum.each.normal.type <- 500
maximum.other.bins <- 200

normal.celltypes <- c(
  'Heptocyte',
  'Hepatocyte',
  'Endothelial',
  'Fibroblast'
)


###### sample information ######

old.samples <- data.frame(
  sample_id = c('238966', '239117', '176125', '135121'),
  data_folder = c(
    '02-Old data/H1-WCPQQ8F-A1-238966-2',
    '02-Old data/H1-WCPQQ8F-D1-239117-2',
    '02-Old data/H1-98MVJHD-D1-176125-2-jiace',
    '02-Old data/H1-JCNDPTV-A1-135121-3'
  ),
  cluster_file = c(
    '238966 cluster barcode.csv',
    '239117 cluster barcode.csv',
    '176125 Graph.based.csv',
    '135121 Graph.based.csv'
  ),
  sample_file = NA,
  slide = NA,
  region_file = c(
    '238966.csv',
    '239117.csv',
    '176125 ta-tumor region.csv',
    '135121 Bile duct region.csv'
  ),
  region_rule = c('tumor_text', 'tumor_text', 'non_empty', 'non_empty')
)

new.samples <- data.frame(
  sample_id = c(
    paste0('S1-', 1:4),
    'S2-1', 'S2-3', 'S2-4',
    paste0('S3-', 1:4),
    paste0('S4-', 1:4)
  ),
  data_folder = c(
    rep('01-New data/LS-1', 4),
    rep('01-New data/LS-2', 3),
    rep('01-New data/LS-3', 4),
    rep('01-New data/LS-4', 4)
  ),
  cluster_file = c(
    rep('slide1 cluster barcode.csv', 4),
    rep('sldie2 cluster barcode.csv', 3),
    rep('slide3 cluster barcode.csv', 4),
    rep('slide4 cluster barcode.csv', 4)
  ),
  sample_file = c(
    rep('slide1 sample barcode.csv', 4),
    rep('slide2 sample barcode.csv', 3),
    rep('slide3 sample barcode.csv', 4),
    rep('slide4 sample barcode.csv', 4)
  ),
  slide = c(
    rep('S1', 4),
    rep('S2', 3),
    rep('S3', 4),
    rep('S4', 4)
  ),
  region_file = c(
    rep('S1.csv', 4),
    rep('S2.csv', 3),
    rep('S3.csv', 4),
    rep('S4.csv', 4)
  ),
  region_rule = 'tumor_text'
)

sample.information <- rbind(old.samples, new.samples)


###### graph-based cluster annotation for the four new slides ######

cluster.annotation <- list(
  S1 = c(
    '1' = 'Immune-Stroma', '2' = 'Hepatocyte',
    '3' = 'Hepatocyte', '4' = 'Immune-Stroma',
    '5' = 'Tumor', '6' = 'Tumor', '7' = 'Bile Duct',
    '8' = 'Tumor', '9' = 'Endothelial',
    '10' = 'Immune-Stroma', '11' = 'unassign',
    '12' = 'unassign'
  ),
  S2 = c(
    '1' = 'Hepatocyte', '2' = 'Immune-Stroma',
    '3' = 'Immune-Stroma', '4' = 'Tumor', '5' = 'Tumor',
    '6' = 'Immune-Stroma', '7' = 'Tumor', '8' = 'Tumor',
    '9' = 'Tumor', '10' = 'Endothelial', '11' = 'Bile Duct'
  ),
  S3 = c(
    '1' = 'Hepatocyte', '2' = 'Fibroblast',
    '3' = 'Immune-Stroma', '4' = 'Immune-Stroma',
    '5' = 'Tumor', '6' = 'Tumor', '7' = 'Fibroblast',
    '8' = 'Tumor', '9' = 'Tumor', '10' = 'Fibroblast',
    '11' = 'Tumor', '12' = 'Bile Duct',
    '13' = 'Endothelial', '14' = 'Immune-Stroma',
    '15' = 'unassign'
  ),
  S4 = c(
    '1' = 'Hepatocyte', '2' = 'Immune-Stroma', '3' = 'Tumor',
    '4' = 'Immune-Stroma', '5' = 'Tumor', '6' = 'Tumor',
    '7' = 'Immune-Stroma', '8' = 'Tumor',
    '9' = 'Immune-Stroma', '10' = 'Bile Duct',
    '11' = 'unassign', '12' = 'Immune-Stroma',
    '13' = 'unassign', '14' = 'Bile Duct'
  )
)


###### choose a fixed number of bins from each annotation ######

select.bins <- function(metadata) {

  bile.duct <- metadata$barcode[metadata$celltype == 'Bile Duct']
  if (length(bile.duct) > maximum.bile.duct.bins) {
    bile.duct <- sample(bile.duct, maximum.bile.duct.bins)
  }

  tumor <- metadata$barcode[metadata$celltype == 'Tumor']
  if (length(tumor) > maximum.tumor.bins) {
    tumor <- sample(tumor, maximum.tumor.bins)
  }

  normal <- c()
  for (one.celltype in normal.celltypes) {
    one.normal <- metadata$barcode[metadata$celltype == one.celltype]
    if (length(one.normal) > maximum.each.normal.type) {
      one.normal <- sample(one.normal, maximum.each.normal.type)
    }
    normal <- c(normal, one.normal)
  }

  other <- metadata$barcode[
    !(metadata$celltype %in% c('Bile Duct', 'Tumor', normal.celltypes))
  ]
  if (length(other) > maximum.other.bins) {
    other <- sample(other, maximum.other.bins)
  }

  unique(c(bile.duct, tumor, normal, other))
}


###### run CopyKAT separately for every specimen ######

all.proportions <- data.frame()

for (i in 1:nrow(sample.information)) {

  one.sample <- sample.information[i, ]
  sample.id <- one.sample$sample_id
  message('CopyKAT: ', sample.id)

  one.result.path <- file.path(result.path, sample.id)
  dir.create(one.result.path, recursive = TRUE, showWarnings = FALSE)

  matrix.file <- file.path(
    base.path,
    one.sample$data_folder,
    'outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5'
  )

  counts <- Read10X_h5(matrix.file)
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  counts <- as(counts, 'dgCMatrix')
  rownames(counts) <- make.unique(rownames(counts))

  cluster <- fread(file.path(annotation.path, one.sample$cluster_file))
  setnames(cluster, names(cluster)[1:2], c('barcode', 'celltype'))

  if (is.na(one.sample$sample_file)) {

    metadata <- unique(cluster[, .(barcode, celltype)])

  } else {

    sample.table <- fread(
      file.path(annotation.path, one.sample$sample_file)
    )
    setnames(sample.table, names(sample.table)[1:2], c('barcode', 'sample_id'))
    sample.table[sample_id == 'S-1', sample_id := 'S1-1']

    cluster[, cluster_number := gsub('[^0-9]', '', celltype)]
    cluster[, celltype := cluster.annotation[[one.sample$slide]][cluster_number]]

    metadata <- merge(
      sample.table[, .(barcode, sample_id)],
      cluster[, .(barcode, celltype)],
      by = 'barcode'
    )
    metadata <- metadata[sample_id == sample.id]
  }

  common.barcodes <- intersect(colnames(counts), metadata$barcode)
  counts <- counts[, common.barcodes, drop = FALSE]
  metadata <- metadata[match(common.barcodes, barcode)]

  metadata[, detected_genes := as.numeric(
    Matrix::colSums(counts > 0)[barcode]
  )]

  # Remove bins with fewer than 200 detected genes.
  metadata <- metadata[detected_genes >= minimum.detected.genes]

  selected.barcodes <- select.bins(metadata)
  counts.copykat <- counts[, selected.barcodes, drop = FALSE]
  metadata.copykat <- metadata[match(selected.barcodes, barcode)]

  # Use annotated normal bins only when at least 20 are available.
  normal.references <- metadata.copykat$barcode[
    metadata.copykat$celltype %in% normal.celltypes
  ]
  if (length(normal.references) < minimum.normal.references) {
    normal.references <- ''
  }

  write.csv(
    metadata.copykat,
    file.path(one.result.path, paste0(sample.id, '_CopyKAT_input.csv')),
    row.names = FALSE
  )

  old.working.directory <- getwd()
  setwd(one.result.path)

  copykat.result <- copykat(
    rawmat = as.matrix(counts.copykat),
    id.type = 'S',
    cell.line = 'no',
    ngene.chr = 1,
    LOW.DR = copykat.low.dr,
    UP.DR = copykat.up.dr,
    win.size = 25,
    norm.cell.names = normal.references,
    KS.cut = copykat.ks.cut,
    sam.name = gsub('-', '_', sample.id),
    distance = 'euclidean',
    n.cores = 1,
    genome = 'hg20',
    plot.genes = FALSE
  )

  setwd(old.working.directory)

  saveRDS(
    copykat.result,
    file.path(one.result.path, paste0(sample.id, '_CopyKAT_result.rds'))
  )

  prediction <- as.data.frame(copykat.result$prediction)
  colnames(prediction)[1] <- 'barcode'
  prediction <- merge(
    metadata.copykat,
    prediction,
    by = 'barcode',
    all.x = TRUE
  )

  write.csv(
    prediction,
    file.path(one.result.path, paste0(sample.id, '_CopyKAT_prediction.csv')),
    row.names = FALSE
  )


  ###### classify the three spatial groups ######

  region <- fread(file.path(region.path, one.sample$region_file))
  setnames(region, names(region)[1:2], c('barcode', 'region_label'))
  region[, barcode := trimws(barcode)]

  if (one.sample$region_rule == 'non_empty') {
    tumor.region.barcodes <- region[
      !is.na(region_label) & trimws(region_label) != '',
      barcode
    ]
  } else {
    tumor.region.barcodes <- region[
      grepl('tumor', region_label, ignore.case = TRUE),
      barcode
    ]
  }

  plot.data <- as.data.table(prediction)
  plot.data <- plot.data[
    celltype %in% c('Bile Duct', 'Tumor') &
      copykat.pred %in% c('diploid', 'aneuploid')
  ]

  plot.data[, group := ifelse(
    celltype == 'Tumor',
    'Tumor',
    ifelse(
      barcode %in% tumor.region.barcodes,
      'Tumor-proximal cholangiocyte-like',
      'Non-tumor cholangiocyte-like'
    )
  )]

  group.order <- c(
    'Tumor-proximal cholangiocyte-like',
    'Non-tumor cholangiocyte-like',
    'Tumor'
  )

  proportion <- plot.data[, .(n = .N), by = .(group, copykat.pred)]

  complete.table <- CJ(
    group = group.order,
    copykat.pred = c('diploid', 'aneuploid'),
    unique = TRUE
  )
  proportion <- merge(
    complete.table,
    proportion,
    by = c('group', 'copykat.pred'),
    all.x = TRUE
  )
  proportion[is.na(n), n := 0]
  proportion[, total := sum(n), by = group]
  proportion[, proportion := ifelse(total == 0, NA, n / total)]
  proportion[, sample_id := sample.id]

  proportion$group <- factor(proportion$group, levels = group.order)
  proportion$copykat.pred <- factor(
    proportion$copykat.pred,
    levels = c('diploid', 'aneuploid')
  )

  write.csv(
    proportion,
    file.path(one.result.path, paste0(sample.id, '_CopyKAT_proportion.csv')),
    row.names = FALSE
  )

  p <- ggplot(
    proportion,
    aes(x = group, y = proportion, fill = copykat.pred)
  ) +
    geom_col(width = 0.7, color = 'grey30', linewidth = 0.2) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_fill_manual(
      values = c('diploid' = '#4E79A7', 'aneuploid' = '#E15759')
    ) +
    labs(
      title = sample.id,
      x = NULL,
      y = 'Proportion',
      fill = 'CopyKAT call'
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(
    file.path(one.result.path, paste0(sample.id, '_CopyKAT_proportion.pdf')),
    p,
    width = 5.2,
    height = 4.2
  )

  all.proportions <- rbind(all.proportions, proportion)
}


###### combined proportion table ######

write.csv(
  all.proportions,
  file.path(result.path, 'all_samples_CopyKAT_proportions.csv'),
  row.names = FALSE
)
