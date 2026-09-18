.libPaths(c("D:/software/Rlibs/R_copy", .libPaths()))

library(data.table)
library(Matrix)
library(Seurat)


# Cholangiocyte-like genes changing along the tumor-proximity axis

base_dir <- "F:/TLS DATA/02-HD"
out_dir <- "F:/TLS Material/06-Code/result/Figure4_HD_along_axis_gene"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "per_sample"), showWarnings = FALSE)

metadata_file <- file.path(
  base_dir,
  "03-analysis/03-BileDuct similarity/01-new-sample-axis-qc/all_new_sample_bins_with_celltype_for_axis_qc.csv"
)

sample_id <- c(
  "S1-1", "S1-4",
  "S2-1", "S2-3", "S2-4",
  "S3-1", "S3-2", "S3-4",
  "S4-1", "S4-2", "S4-3"
)


# 1. Read the bin annotation and retain cholangiocyte-like bins
metadata <- fread(metadata_file)
metadata <- metadata[
  sample %in% sample_id &
    cell_type %in% c("Bile Duct", "Cholangiocyte")
]

dim(metadata)
table(metadata$sample)
table(metadata$cell_type)


# 2. Define the tumor-facing axis separately for each sample

axis_list <- list()
axis_summary <- list()

for (sid in sample_id) {
  temp <- copy(metadata[sample == sid])
  n_before <- nrow(temp)

  if (sid == "S1-1") {
    temp <- temp[array_col >= 700 & array_col <= 1000]
    temp[, axis_value := array_row]
    temp[, axis_name := "y = array_row, x 700-1000"]
  }

  if (sid %in% c("S1-4", "S2-1", "S2-3", "S3-1", "S3-2")) {
    temp[, axis_value := array_col]
    temp[, axis_name := "x = array_col"]
  }

  if (sid == "S2-4") {
    temp <- temp[array_row >= 1000 & array_row <= 1400]
    temp[, axis_value := array_col]
    temp[, axis_name := "x = array_col, y 1000-1400"]
  }

  if (sid == "S3-4") {
    temp <- temp[array_row >= 900 & array_row <= 1400]
    temp[, axis_value := array_col]
    temp[, axis_name := "x = array_col, y 900-1400"]
  }

  if (sid == "S4-1") {
    temp <- temp[array_col >= 900 & array_col <= 1200]
    temp[, axis_value := array_col]
    temp[, axis_name := "x = array_col, x 900-1200"]
  }

  if (sid == "S4-2") {
    x1 <- 700
    y1 <- 800
    x2 <- 100
    y2 <- 200
    dx <- x2 - x1
    dy <- y2 - y1
    axis_length <- sqrt(dx^2 + dy^2)

    temp[, axis_value :=
      ((array_col - x1) * dx + (array_row - y1) * dy) / axis_length]
    temp <- temp[axis_value >= 0 & axis_value <= axis_length]
    temp[, axis_name := "diagonal: (700,800) to (100,200)"]
  }

  if (sid == "S4-3") {
    temp <- temp[array_col >= 1000 & array_col <= 1250]
    temp[, axis_value := array_col]
    temp[, axis_name := "x = array_col, x 1000-1250"]
  }

  temp[, tumor_proximity := max(axis_value) - axis_value]

  axis_list[[sid]] <- temp
  axis_summary[[sid]] <- data.table(
    sample = sid,
    n_chol_before_axis_selection = n_before,
    n_chol_used = nrow(temp),
    axis_name = unique(temp$axis_name),
    axis_min = min(temp$axis_value),
    axis_max = max(temp$axis_value)
  )
}

axis_data <- rbindlist(axis_list)
axis_summary <- rbindlist(axis_summary)

axis_summary
table(axis_data$sample)
dim(axis_data)

fwrite(
  axis_data,
  file.path(out_dir, "cholangiocyte_bins_with_tumor_proximity.csv")
)
fwrite(
  axis_summary,
  file.path(out_dir, "sample_axis_summary.csv")
)


# 3. Read the raw 8-um count matrices
library_table <- data.table(
  ls_id = c("LS-1", "LS-2", "LS-3", "LS-4"),
  h5_file = file.path(
    base_dir,
    c(
      "01-New data/LS-1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
      "01-New data/LS-2/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
      "01-New data/LS-3/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
      "01-New data/LS-4/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"
    )
  )
)


# 4. For each gene, calculate its Spearman correlation with tumor proximity
all_result <- list()

for (one_ls in library_table$ls_id) {
  message("Reading ", one_ls)

  count <- Read10X_h5(library_table[ls_id == one_ls, h5_file])
  if (is.list(count)) count <- count[[1]]
  count <- as(count, "dgCMatrix")
  rownames(count) <- make.unique(rownames(count))

  samples_in_library <- unique(axis_data[ls_id == one_ls, sample])

  for (sid in samples_in_library) {
    temp_axis <- axis_data[ls_id == one_ls & sample == sid]
    common_barcode <- intersect(colnames(count), temp_axis$Barcode)
    temp_axis <- temp_axis[match(common_barcode, Barcode)]
    temp_count <- count[, common_barcode, drop = FALSE]

    print(sid)
    print(dim(temp_count))
    print(summary(temp_axis$tumor_proximity))

    n_bin <- ncol(temp_count)
    positive_percentage <- Matrix::rowSums(temp_count > 1) / n_bin * 100
    mean_expression <- Matrix::rowMeans(temp_count)

    rho <- rep(NA_real_, nrow(temp_count))

    # The genes are handled in small blocks to avoid converting the full
    # sparse expression matrix into one very large dense matrix.
    for (start in seq(1, nrow(temp_count), by = 500)) {
      end <- min(start + 499, nrow(temp_count))
      temp_matrix <- as.matrix(temp_count[start:end, , drop = FALSE])

      rho[start:end] <- apply(temp_matrix, 1, function(gene_expression) {
        if (sd(gene_expression) == 0) return(NA_real_)
        suppressWarnings(cor(
          gene_expression,
          temp_axis$tumor_proximity,
          method = "spearman"
        ))
      })
    }

    rho_for_test <- pmin(pmax(rho, -0.999999), 0.999999)
    test_statistic <- rho_for_test * sqrt((n_bin - 2) / (1 - rho_for_test^2))
    p_value <- 2 * pt(-abs(test_statistic), df = n_bin - 2)
    p_value[is.na(rho)] <- NA_real_

    result <- data.table(
      sample = sid,
      Gene = rownames(temp_count),
      Correlation = rho,
      PValue = p_value,
      AdjustedPValue = p.adjust(p_value, method = "BH"),
      PositivePercentage = as.numeric(positive_percentage),
      MeanExpression = as.numeric(mean_expression),
      cholangiocyte_bins = n_bin,
      axis_name = unique(temp_axis$axis_name)
    )

    all_result[[sid]] <- result

    fwrite(
      result,
      file.path(out_dir, "per_sample", paste0(sid, "_all_gene_correlation.csv"))
    )
    fwrite(
      result[Correlation > 0 & PValue < 0.05 & PositivePercentage > 0],
      file.path(out_dir, "per_sample", paste0(sid, "_positive_genes_P005.csv"))
    )
  }
}

new_result <- rbindlist(all_result, fill = TRUE)

dim(new_result)
table(new_result$sample)

fwrite(
  new_result,
  file.path(out_dir, "all_new_sample_gene_correlation.csv")
)


# 5. Add the two previous samples

old_238966 <- fread(file.path(base_dir, "03-analysis/238966 significant gene.csv"))
old_238966 <- old_238966[, .(
  sample = "238966",
  Gene,
  Correlation,
  PValue,
  AdjustedPValue = NA_real_,
  PositivePercentage,
  MeanExpression = NA_real_,
  cholangiocyte_bins = NA_integer_,
  axis_name = NA_character_
)]

old_239117 <- fread(file.path(base_dir, "03-analysis/239117 significant gene.csv"))
old_239117 <- old_239117[, .(
  sample = "239117",
  Gene,
  Correlation,
  PValue,
  AdjustedPValue = NA_real_,
  PositivePercentage,
  MeanExpression = NA_real_,
  cholangiocyte_bins = NA_integer_,
  axis_name = NA_character_
)]

new_positive <- new_result[
  Correlation > 0 & PValue < 0.05 & PositivePercentage > 0
]

all_positive <- rbindlist(
  list(new_positive, old_238966, old_239117),
  fill = TRUE
)

dim(all_positive)
table(all_positive$sample)

fwrite(
  all_positive,
  file.path(out_dir, "all_old_new_positive_genes_P005.csv")
)


# 6. Count in how many samples each positive gene is observed
gene_recurrence <- all_positive[, .(
  n_samples = uniqueN(sample),
  samples = paste(sort(unique(sample)), collapse = ";"),
  mean_correlation = mean(Correlation, na.rm = TRUE),
  minimum_p_value = min(PValue, na.rm = TRUE),
  mean_positive_percentage = mean(PositivePercentage, na.rm = TRUE)
), by = Gene]

gene_recurrence <- gene_recurrence[order(-n_samples, -mean_correlation)]

table(gene_recurrence$n_samples)
head(gene_recurrence)

fwrite(
  gene_recurrence,
  file.path(out_dir, "tumor_proximity_gene_recurrence.csv")
)
fwrite(
  gene_recurrence[n_samples >= 2],
  file.path(out_dir, "TA_Chol_genes_P005_at_least_2_samples.csv")
)
