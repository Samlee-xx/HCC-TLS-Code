.libPaths(c("D:/software/Rlibs/R_copy", .libPaths()))

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


input_dir <- "F:/TLS Material/06-Code/result/Figure4_HD_along_axis_gene"
out_dir <- file.path(input_dir, "GO_BP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_result <- fread(
  file.path(input_dir, "all_old_new_positive_genes_P005.csv")
)

dim(gene_result)
table(gene_result$sample)

go_gene_p_cutoff <- 0.01

go_summary <- list()

for (sid in unique(gene_result$sample)) {
  significant_gene <- gene_result[
    sample == sid &
      Correlation > 0 &
      PValue <= go_gene_p_cutoff &
      PositivePercentage > 0,
    unique(Gene)
  ]

  print(sid)
  print(length(significant_gene))

  gene_id <- mapIds(
    org.Hs.eg.db,
    keys = significant_gene,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  gene_id <- unique(na.omit(gene_id))

  go_result <- enrichGO(
    gene = gene_id,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  go_table <- as.data.table(go_result@result)

  if (nrow(go_table) > 0) {
    ratio_part <- strsplit(go_table$GeneRatio, "/")
    go_table[, GeneRatioNumeric := vapply(
      ratio_part,
      function(x) as.numeric(x[1]) / as.numeric(x[2]),
      numeric(1)
    )]
  }

  fwrite(
    go_table,
    file.path(out_dir, paste0(sid, "_GO_BP_all.csv"))
  )

  go_summary[[sid]] <- data.table(
    sample = sid,
    input_gene_number = length(significant_gene),
    mapped_gene_number = length(gene_id),
    significant_GO_number = sum(go_table$p.adjust <= 0.05, na.rm = TRUE)
  )
}

go_summary <- rbindlist(go_summary)

go_summary
fwrite(go_summary, file.path(out_dir, "GO_BP_summary.csv"))


