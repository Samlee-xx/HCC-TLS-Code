
library(Seurat)
setwd('...\\14 HCC Neutrophil')
neu <- readRDS('...\\nature_neu_addname_after_specific_patient.rds')
marker <- FindAllMarkers(neu,only.pos = TRUE)
DimPlot(neu,group.by = 'seurat_clusters',label = TRUE)
neu <- subset(neu,seurat_clusters %in% c(0,1,2,3,4,5,6,8,9)) # remove unknown clusters 
neu <- FindVariableFeatures(neu)
neu <- ScaleData(neu,features = VariableFeatures(neu))
neu <- RunPCA(object = neu, npcs = 30, verbose = TRUE)
neu <- RunUMAP(object = neu, dims = 1:30, verbose = TRUE)
neu <- FindNeighbors(object = neu, dims = 1:30, verbose = TRUE)
neu <- FindClusters(object = neu, verbose = TRUE, resolution = 0.2)

DimPlot(neu,label = TRUE)
FeaturePlot(neu,features = c('CCR1'))
marker = FindAllMarkers(neu,only.pos = TRUE)

new.cluster.ids <- c('Neu_S0_S100A8',
                     'Neu_S1_APOC3',
                     'Neu_S2_MIF',
                     'Neu_S3_ISG15',
                     'Neu_S4_S100A12',
                     'Neu_S5_CD74')

table(new.cluster.ids)
names(new.cluster.ids) <- levels(neu)
neu<- RenameIdents(neu, new.cluster.ids)
neu$celltype <- neu@active.ident
DimPlot(neu,group.by = 'celltype')
saveRDS(neu,'01 Neu seurat.rds')






