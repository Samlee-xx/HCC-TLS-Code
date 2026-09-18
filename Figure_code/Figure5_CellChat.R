
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

combined <- readRDS('01 Combined Immunecell.rds')
data.input <- combined[["RNA"]]
labels <- combined$celltype
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
table(meta$labels)
cellchat <- createCellChat(object = data.input@data, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human
db_int <- CellChatDB[["interaction"]]

new_pairs <- data.frame(
  interaction_name = "LGALS3_LAG3",
  ligand = "LGALS3",
  receptor = "LAG3",
  agonist = NA,
  antagonist = NA,
  cofactor = NA,
  evidence = "literature",
  annotation = "Cell-Cell Contact",   
  pathway_name = "LGALS3",          
  stringsAsFactors = FALSE
)

miss_cols <- setdiff(colnames(db_int), colnames(new_pairs))
for (m in miss_cols) new_pairs[[m]] <- NA
new_pairs <- new_pairs[, colnames(db_int)]

CellChatDB2 <- CellChatDB
CellChatDB2[["interaction"]] <- rbind(db_int, new_pairs)

CellChatDB.use <- subsetDB(CellChatDB2,
                           search = c("Secreted Signaling", "Cell-Cell Contact"),
                           key    = "annotation")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE, trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
cellchat@LR$LRsig$interaction_name_2
ix <- which(is.na(cellchat@LR$LRsig$interaction_name_2) & 
              cellchat@LR$LRsig$interaction_name == "LGALS3_LAG3")
cellchat@LR$LRsig$interaction_name_2[ix] <- "LGALS3 - LAG3"

a = netVisual_bubble(cellchat, sources.use = c(4,5,6,7,8,9), targets.use = c(1,2,3), remove.isolate = FALSE)+
  scale_color_gradient(low = "#D9D9D9", high = "#C25E72")

a = netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(4,5,6,7,8,9), remove.isolate = FALSE)+
  scale_color_gradient(low = "#D9D9D9", high = "#C25E72")


saveRDS(cellchat,'Cellchat Analysis 0822.rds')
cellchat <- readRDS('Cellchat Analysis 0822.rds')

a <- netVisual_aggregate(
  cellchat,
  sources.use = 6,
  targets.use = c(1, 2, 3),
  signaling = cellchat@netP$pathways
)


