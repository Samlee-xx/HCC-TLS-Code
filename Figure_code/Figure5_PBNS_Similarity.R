
library(dplyr)
library(Seurat)
neu.hcc.single <- readRDS('01 Neu seurat.rds')
marker.neu <- FindAllMarkers(neu.hcc.single,only.pos = TRUE)
counts_matrix <- GetAssayData(neu.hcc.single, slot = "counts")
metadata <- neu.hcc.single@meta.data
counts_df <- t(counts_matrix) %>% as.data.frame(counts_matrix) 
sum(!colnames(counts_matrix) == rownames(metadata))
rownames(counts_df)
counts_df$celltype <- metadata$celltype

bulk_counts_df <- counts_df %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
bulk_counts_df = bulk_counts_df %>% as.data.frame()
rownames(bulk_counts_df) = bulk_counts_df$celltype
bulk_counts_df <- bulk_counts_df[,-1] %>% t() %>% as.data.frame()
bulk_counts_df <- bulk_counts_df[rowMeans(bulk_counts_df)>1,]
single.bulk <- bulk_counts_df

bulk.data <- read.csv('gene_sample_count_with_symbol.csv',
                     sep = ',',row.names = 1)
bulk.data$Symbol <- ifelse(bulk.data$Symbol == '---',NA,bulk.data$Symbol)
bulk.data <- bulk.data[complete.cases(bulk.data$Symbol),]
gene.common <- intersect(rownames(bulk.data),rownames(bulk_counts_df))
# Remove duplicate rows based on specific columns 
bulk.data <- bulk.data %>% distinct(Symbol, .keep_all = TRUE) 
rownames(bulk.data) <- bulk.data$Symbol
colnames(bulk.data)
bulk.data <- bulk.data[,-c(1,2)]
bulk.data <- bulk.data[rowMeans(bulk.data)>1,]
colnames(bulk.data)

gene.single <- VariableFeatures(neu.hcc.single)
common.gene <- intersect(gene.single,rownames(bulk.data)) # 1483

bulk.data <- bulk.data[common.gene,]
single.bulk <- single.bulk[common.gene,]

merge.data <- cbind(bulk.data,single.bulk)
merge.data <- na.omit(merge.data)
pheatmap(cor(merge.data))

library_sizes <- colSums(merge.data)
scaling_factors <- library_sizes / mean(library_sizes)
normalized_counts <- sweep(merge.data, 2, scaling_factors, FUN = "/")
normalized_counts <- normalized_counts * 1e6
pheatmap(cor(normalized_counts))

cor.data <- cor(normalized_counts)
cor.data.sub <- cor.data[c(16:21),
                        c(1:15)] %>% t() %>%as.data.frame()
rownames(cor.data.sub)
cor.data.sub$group <- c('IFNG','IFNG','IFNG','TAN','TAN','TAN','LPS','LPS','LPS','TAN','TAN','TAN',
                       'Control','Control','Control')

cor.data.sub <- cor.data.sub %>% group_by(group) %>% summarise_each(funs = median) %>%as.data.frame()
rownames(cor.data.sub) <- cor.data.sub$group
cor.data.sub <- cor.data.sub[,-1]
cor.data.sub <- scale(cor.data.sub)

colnames(cor.data.sub)
data <- data.frame(
  ContionType = rownames(cor.data.sub),
  ScaledCor = cor.data.sub[,6]
)

data$ContionType <- factor(data$ContionType,levels = rownames(data))
a = ggplot(data, aes(x = reorder(ContionType, ScaledCor), y = ScaledCor)) +
  geom_segment(aes(x = ContionType, xend = ContionType, y = 0, yend = ScaledCor), color = "grey",size =1) + 
  geom_point(aes(color = ContionType), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  labs(x = NULL, y = "Neu_C3_ISG15") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    "Control" = "#E41A1C", "IFNG" = "#377EB8", "LPS" = "#4DAF4A", "TAN" = "#984EA3"
  ))
a

