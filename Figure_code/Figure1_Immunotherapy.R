###### TLS gene signatures in immunotherapy bulk RNA-seq cohorts ######

library(IOBR)
library(dplyr)
library(ggpubr)

###### TLS subtype gene signatures ######

marker <- readRDS(
  'E:/01-TLS/01-Omics/01-ST/03-analysis final/02-DEGs/DEGs.rds'
)

marker <- marker[marker$p_val_adj <= 0.05, ]

tls.signature <- list(
  Type1_B_TLS = marker$gene[marker$cluster == 'TLS.2'],
  Type2_T_TLS = marker$gene[marker$cluster == 'TLS.0']
)

length(tls.signature$Type1_B_TLS)
length(tls.signature$Type2_T_TLS)

###### EGAS00001005503: IMbrave150 ######

imbrave.info <- read.delim(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/03-Imbrave150/IMbrave150.cli.txt',
  check.names = FALSE
)

imbrave.tpm <- read.delim(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/03-Imbrave150/IMbrave150.exp.txt',
  check.names = FALSE
)

imbrave.tpm <- imbrave.tpm[!duplicated(imbrave.tpm$Symbol), ]
rownames(imbrave.tpm) <- imbrave.tpm$Symbol
imbrave.tpm <- imbrave.tpm[, -c(1, 2)]

dim(imbrave.tpm)
table(imbrave.info$Visit)
table(imbrave.info$Treatment)
sum(!imbrave.info$anon_sampleId %in% colnames(imbrave.tpm))

imbrave.tpm <- imbrave.tpm[, imbrave.info$anon_sampleId]

imbrave.score <- calculate_sig_score(
  pdata = NULL,
  eset = imbrave.tpm,
  signature = tls.signature,
  method = 'zscore',
  mini_gene_count = 0,
  adjust_eset = TRUE
)

imbrave.score <- as.data.frame(imbrave.score)
imbrave.score$Index <- NULL

imbrave <- left_join(
  imbrave.score,
  imbrave.info,
  by = c('ID' = 'anon_sampleId')
)

# Baseline samples from the atezolizumab plus bevacizumab arm.
imbrave <- imbrave[
    imbrave$Treatment == 'Atezolizumab+Bevacizumab' &
    imbrave$Confirmed.Response_IRF %in% c('CR', 'PR', 'SD', 'PD'),
]

imbrave$response <- factor(
  imbrave$Confirmed.Response_IRF,
  levels = c('PD', 'SD', 'PR', 'CR')
)

dim(imbrave)
table(imbrave$response)

imbrave.comparisons <- list(
  c('PD', 'SD'),
  c('PD', 'PR'),
  c('PD', 'CR'),
  c('SD', 'PR'),
  c('SD', 'CR'),
  c('PR', 'CR')
)

imbrave.b.test <- compare_means(
  Type1_B_TLS ~ response,
  data = imbrave,
  method = 'wilcox.test',
  comparisons = imbrave.comparisons,
  p.adjust.method = 'holm'
)

imbrave.t.test <- compare_means(
  Type2_T_TLS ~ response,
  data = imbrave,
  method = 'wilcox.test',
  comparisons = imbrave.comparisons,
  p.adjust.method = 'holm'
)

imbrave.b.test
imbrave.t.test

ggboxplot(
  imbrave,
  x = 'response',
  y = 'Type1_B_TLS',
  fill = 'response',
  add = 'jitter',
  outlier.shape = NA
) + theme_classic()

ggboxplot(
  imbrave,
  x = 'response',
  y = 'Type2_T_TLS',
  fill = 'response',
  add = 'jitter',
  outlier.shape = NA
) + theme_classic()

###### GSE235863 ######

gse.tpm <- read.delim(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/04-ICI bulk/GSE235863_bulk_rna_seq_tpm.txt',
  row.names = 1,
  check.names = FALSE
)

gse.tpm <- gse.tpm[!duplicated(gse.tpm$geneName), ]
rownames(gse.tpm) <- gse.tpm$geneName
gse.tpm <- gse.tpm[, -1]

dim(gse.tpm)
colnames(gse.tpm)

gse.score <- calculate_sig_score(
  pdata = NULL,
  eset = gse.tpm,
  signature = tls.signature,
  method = 'zscore',
  mini_gene_count = 0,
  adjust_eset = TRUE
)

gse.score <- as.data.frame(gse.score)
gse.score$Index <- NULL

gse.response <- c(
  P1 = 'NR',
  P3 = 'PR',
  P5 = 'PR',
  P16 = 'NR',
  P18 = 'CR',
  P26 = 'NR',
  P27 = 'PR',
  P31 = 'PR',
  P32 = 'PR',
  P34 = 'NR',
  P35 = 'CR',
  P36 = 'PR',
  P38 = 'CR',
  P40 = 'CR',
  P41 = 'CR'
)

gse.score$response <- gse.response[gse.score$ID]
gse.score$response <- factor(gse.score$response, levels = c('NR', 'PR', 'CR'))

dim(gse.score)
table(gse.score$response, useNA = 'ifany')

gse.comparisons <- list(
  c('NR', 'PR'),
  c('NR', 'CR'),
  c('PR', 'CR')
)

gse.b.test <- compare_means(
  Type1_B_TLS ~ response,
  data = gse.score,
  method = 'wilcox.test',
  comparisons = gse.comparisons,
  p.adjust.method = 'holm'
)

gse.t.test <- compare_means(
  Type2_T_TLS ~ response,
  data = gse.score,
  method = 'wilcox.test',
  comparisons = gse.comparisons,
  p.adjust.method = 'holm'
)

gse.b.test
gse.t.test

ggboxplot(
  gse.score,
  x = 'response',
  y = 'Type1_B_TLS',
  color = 'response',
  add = 'jitter',
  outlier.shape = NA
) + theme_classic()

ggboxplot(
  gse.score,
  x = 'response',
  y = 'Type2_T_TLS',
  color = 'response',
  add = 'jitter',
  outlier.shape = NA
) + theme_classic()
