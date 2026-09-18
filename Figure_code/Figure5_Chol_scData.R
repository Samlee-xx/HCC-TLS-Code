
###### Bile Duct single cell ######

normal.data = readRDS('...\\Normal Liver Buli Duct.rds')
hcc.data = readRDS('...\\Bilt Duct from HCC sample.rds')
DimPlot(hcc.data)
FeaturePlot(hcc.data, features = c('MUC6','EPCAM','KAT7'))+DimPlot(hcc.data,label = TRUE)
hcc.data = subset(hcc.data,seurat_clusters %in% c(0,7))

normal.data$group = 'Normal'
hcc.data$group = 'HCC'

normal.data <- NormalizeData(normal.data)
hcc.data <- NormalizeData(hcc.data)
anchors <- FindIntegrationAnchors(object.list = list(normal.data, hcc.data), dims = 1:30)
expr <- IntegrateData(anchorset = anchors)
DefaultAssay(expr) <- "integrated" 
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(expr, verbose = FALSE)
PCS <- 30
k.param <- 20
expr <- RunPCA(expr, features = VariableFeatures(expr), verbose = FALSE)
expr <- FindNeighbors(expr, dims = 1:PCS, verbose = FALSE, k.param = k.param)
expr <- RunUMAP(expr, dims = 1:PCS, verbose = FALSE) 
expr <- FindClusters(expr, resolution = 0.1, verbose = FALSE)
DimPlot(expr,cols = c("#DA8DC0", '#AFD766', "#92A0CD",'#EA966B'))

saveRDS(expr,'E:\\01-TLS\\04-SingleCell\\06-Nature Subtypes\\03-Analysis\\Bile Duct.rds')
expr = readRDS('E:\\01-TLS\\04-SingleCell\\06-Nature Subtypes\\03-Analysis\\Bile Duct.rds')
a = FeaturePlot(expr,features = c('MUC1','EPCAM','KRT19','SOX9','CFTR'),min.cutoff = 0)
pdf(file ='E:\\01-TLS\\04-SingleCell\\06-Nature Subtypes\\03-Analysis\\Bile Duct marker.pdf',height = 8,width = 7.5)
print(a)
dev.off()

a = DimPlot(expr,cols = c("#DA8DC0", '#AFD766', "#92A0CD",'#EA966B'))
pdf(file ='E:\\01-TLS\\04-SingleCell\\06-Nature Subtypes\\03-Analysis\\UMAP Bile Duct.pdf',height = 4,width = 5)
print(a)
dev.off()

DimPlot(expr, group.by = 'group')

Idents(expr) = expr$seurat_clusters
VlnPlot(expr,features = c('MUC6','KRT7','EPCAM','CFTR'))
FeaturePlot(expr,features = c('MUC6','KRT7','EPCAM','CFTR'),min.cutoff = 2)+DimPlot(expr, group.by = 'group')
FeaturePlot(expr,features = c('CXCL6','CXCL8','CXCL5','CXCL1'))+DimPlot(expr, group.by = 'group')
marker = FindAllMarkers(expr,only.pos = TRUE)
write.csv(marker,'...\\BileDuct Marker.csv')

library(ggplot2)
library(patchwork)

p <- VlnPlot(expr, 
             features = c('CXCL6', 'CXCL8', 'CXCL1', 'CCL2'), 
             ncol = 2, 
             cols = c("#DA8DC0", '#AFD766', "#92A0CD", '#EA966B'),
             pt.size = 0, combine = FALSE)
p <- lapply(p, function(g) g + geom_jitter(size = 0.5, width = 0.2, alpha = 0.3))

a = wrap_plots(p)
pdf(file ='...\\CXCL combined.pdf',height = 4,width =6)
print(a)
dev.off()

FeaturePlot(expr,features = c('CXCL8'),max.cutoff = 3,min.cutoff = 0)
a = FeaturePlot(expr, features = c('CXCL8'), max.cutoff = 3, min.cutoff = 0, cols = c("white", "#9362ff"))

library(org.Hs.eg.db)
library(clusterProfiler)

marker$enterz <- mapIds(org.Hs.eg.db, keys = marker$gene, keytype = "SYMBOL", column="ENTREZID")
enrich <- compareCluster(enterz ~cluster, data=marker,
                         fun="enrichGO",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP" ,
                         pAdjustMethod = "BH",
                         readabl = TRUE)
dotplot(enrich,showCategory=10)
result=enrich@compareClusterResult

marker = subset(marker,cluster==0)
GO <- enrichGO(
  gene = marker$enterz,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  ont = "BP",
  readable = TRUE
)
result = GO@result 
Go.sub <- GO@result %>%
  filter(p.adjust <= 0.05) %>%
  filter(Description %in% c(
    'epithelial cell proliferation',
    'negative regulation of immune system process',
    'neutrophil chemotaxis',
    'neutrophil migration',
    'neutrophil activation'
  )) %>%
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

Go.sub$Description = factor(Go.sub$Description,levels = c(
  'epithelial cell proliferation',
  'negative regulation of immune system process',
  'neutrophil chemotaxis',
  'neutrophil migration',
  'neutrophil activation'
))

a = ggplot(data = Go.sub, aes(x = GeneRatio, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_gradientn(colors = c("#d7e3fc", "#c1d3fe", "#93b5ff", "#7a71f8"))+
  #scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "GeneRatio", y = "", title = "") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5),
    axis.text.y = element_blank()
  ) +
  geom_text(aes(x = 0.000, label = Description), hjust = 0) 

result = GO@result 
Go.sub <- GO@result %>%
  filter(p.adjust <= 0.05)

