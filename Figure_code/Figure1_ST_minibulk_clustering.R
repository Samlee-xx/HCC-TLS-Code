###### ST TLS analysis from mini-bulk counts ######

library(sva)
library(IOBR)
library(Seurat)

set.seed(777) # I like the 7 number

# TLS regions were identified before this step using H&E, TLS scores and SPATA2.
# This script starts from the raw mini-bulk count matrix of the TLS regions.

load('E:/01-TLS/01-Omics/01-ST/02-New data/02-analysis/02-TLS Subtype/Version1.rds')

tls.info <- read.csv(
  'E:/01-TLS/01-Omics/01-ST/04-analysis final/08-Alternative_deconvolution/inputs/tls_metadata_461.csv',
  check.names = FALSE
)

tls.info <- tls.info[, c('TLS_ID', 'sample', 'location')]

names(batch.level) <- colnames(tls.total)
tls.count <- as.matrix(tls.total[, tls.info$TLS_ID])
batch <- batch.level[tls.info$TLS_ID]

dim(tls.count)
length(unique(tls.info$sample))
length(unique(batch))
table(tls.info$sample)
table(batch)

tls.info$TLS_location <- ifelse(tls.info$location == 'N',
                                'N-TLS', 'TA-TLS')
table(tls.info$TLS_location)

###### remove batch effect ######

adjusted <- ComBat_seq(tls.count, batch = batch)
dim(adjusted)

###### convert counts to TPM ######

# This gene-length table was prepared from GENCODE v41.
gene.length <- read.csv(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/01-TLS subtypes/gene_length_20220807.csv'
)

gene.length <- gene.length[!duplicated(gene.length$gene_name), ]
rownames(gene.length) <- gene.length$gene_name

common.gene <- intersect(rownames(adjusted), rownames(gene.length))
adjusted <- adjusted[common.gene, ]
gene.length <- gene.length[common.gene, ]

TPM <- count2tpm(adjusted,
                 effLength = gene.length,
                 id = 'gene_name',
                 length = 'gene_length',
                 gene_symbol = 'gene_name')

dim(TPM)

###### PCA and clustering ######

rownames(tls.info) <- tls.info$TLS_ID
tls.info <- tls.info[colnames(TPM), ]

tls <- CreateSeuratObject(counts = TPM,
                          meta.data = tls.info,
                          project = 'TLS',
                          min.cells = 0,
                          min.features = 0)

tls@assays[['RNA']]@scale.data <- as.matrix(TPM)

tls <- RunPCA(tls, features = rownames(tls), npcs = 50)
ElbowPlot(tls)

tls <- FindNeighbors(tls, dims = 1:2, k.param = 25)
tls <- FindClusters(tls, resolution = 0.3, random.seed = 777)

table(tls$seurat_clusters)
table(tls$TLS_location, tls$seurat_clusters)

tls <- RunUMAP(tls, dims = 1:2, seed.use = 777)
DimPlot(tls, reduction = 'umap', label = TRUE)
