###### TLS deconvolution by CIBERSORTx ######
# ?????????????????????????????????,???????????????????????????
count.file <- paste0(
  'E:/01-TLS/01-Omics/01-ST/04-analysis final/',
  '08-Alternative_deconvolution/inputs/tls_counts_461_gene_symbol.tsv'
)

reference.file <- paste0(
  'E:/01-TLS/01-Omics/01-ST/04-analysis final/',
  '08-Alternative_deconvolution/inputs/',
  'cibersortx_reference_sc_11types_2200cells.txt'
)

mixture.file <- paste0(
  'E:/01-TLS/01-Omics/01-ST/04-analysis final/',
  '08-Alternative_deconvolution/inputs/',
  'cibersortx_mixture_461TLS_CPM.txt'
)

result.file <- paste0(
  'E:/01-TLS/01-Omics/01-ST/04-analysis final/',
  '08-Alternative_deconvolution/CIBERSORTx_web_upload/',
  '03_CIBERSORTx_Results.csv'
)

###### prepare the TLS mixture ######

tls.count <- read.delim(count.file,
                        row.names = 1,
                        check.names = FALSE)

dim(tls.count)

tls.cpm <- sweep(tls.count, 2, colSums(tls.count), '/') * 1000000

dim(tls.cpm)
range(colSums(tls.cpm))

# The prepared reference contains 11 cell types and 200 cells per cell type.
reference.header <- strsplit(readLines(reference.file, n = 1), '\t')[[1]]
table(reference.header[-1])

# ???????????????????????????
# Files used on the CIBERSORTx website:
# reference.file: single-cell reference matrix
# mixture.file: TLS CPM mixture matrix
#
# Website settings:
# Custom signature matrix; scRNA-seq; minimum expression = 0.5
# Relative fractions; S-mode batch correction
# Quantile normalization disabled; permutations = 100

###### read the CIBERSORTx result ######

cibersortx <- read.csv(result.file,
                       check.names = FALSE)

colnames(cibersortx)[1] <- 'TLS_ID'

cell.type <- c(
  'B Cell',
  'Plasma Cell',
  'CD4 T Cell',
  'CD8 T Cell',
  'Dendritic Cell',
  'Macrophage',
  'Neutrophil',
  'Fibroblast',
  'Endothelial',
  'Epithelial',
  'Bile_Duct'
)

cibersortx.fraction <- as.matrix(cibersortx[, cell.type])
cibersortx.fraction <- cibersortx.fraction / rowSums(cibersortx.fraction)
rownames(cibersortx.fraction) <- cibersortx$TLS_ID

dim(cibersortx.fraction)
range(rowSums(cibersortx.fraction))
head(cibersortx.fraction)

summary(cibersortx$Correlation)
summary(cibersortx$RMSE)
