
rm(list = ls())
set.seed(234)
options(stringsAsFactors = F)

setwd("")
cts <- read.csv("cts.csv", row.names = 1, header = T) %>% as.matrix()
cts <- cts[,c(13:15,1:12)]
saveRDS(cts, "cts.rds")

### IFNG vs quanpei 
cts <- readRDS("cts.rds")
cts <- cts[rowMeans(cts)>1,]
cts <- cts[,c(1:3,4:6)]

# create metadata
condition <- factor(c("Neu_IFNG","Neu_IFNG","Neu_IFNG",
                      "Neu_quanpei","Neu_quanpei","Neu_quanpei"))
colData <- data.frame(row.names=colnames(cts), condition)

# create dds
dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds1 <- DESeq(dds)

resultsNames(dds1)
res <- results(dds1)
summary(res)

# extract degs matrix
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up <- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]
res1_up$ENSEMBL <- rownames(res1_up)
res1_up <- res1_up[,c(7,1:6)]

s2e <- bitr(rownames(res1_up),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)
res1_up <- left_join(res1_up, s2e, by = "ENSEMBL")
res1_up <- res1_up[!duplicated(res1_up[[1]]), ]
write.csv(res1_up, "Neu_IFNG_upRegulated_Genes.csv")

res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]
res1_down$ENSEMBL <- rownames(res1_down)
res1_down <- res1_down[,c(7,1:6)]

s2e <- bitr(rownames(res1_down),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)

res1_down <- left_join(res1_down, s2e, by = "ENSEMBL")
res1_down <- res1_down[!duplicated(res1_down[[1]]), ]

write.csv(res1_down, "Neu_IFNG_downRegulated_Genes.csv")

### LPS vs quanpei 
cts <- readRDS("cts.rds")
cts <- cts[rowMeans(cts)>1,]
cts <- cts[,c(1:3,10:12)]

condition <- factor(c("Neu_LPS","Neu_LPS","Neu_LPS",
                      "Neu_quanpei","Neu_quanpei","Neu_quanpei"))
colData <- data.frame(row.names=colnames(cts), condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds1 <- DESeq(dds)

resultsNames(dds1)
res <- results(dds1)
summary(res)

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up <- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]
res1_up$ENSEMBL <- rownames(res1_up)
res1_up <- res1_up[,c(7,1:6)]

s2e <- bitr(rownames(res1_up),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)

res1_up <- left_join(res1_up, s2e, by = "ENSEMBL")
res1_up <- res1_up[!duplicated(res1_up[[1]]), ]
write.csv(res1_up, "Neu_LPS_upRegulated_Genes.csv")

res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]
res1_down$ENSEMBL <- rownames(res1_down)
res1_down <- res1_down[,c(7,1:6)]

s2e <- bitr(rownames(res1_down),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)
res1_down <- left_join(res1_down, s2e, by = "ENSEMBL")
res1_down <- res1_down[!duplicated(res1_down[[1]]), ]

write.csv(res1_down, "Neu_LPS_downRegulated_Genes.csv")


#### Merge LM3 PLC as TAN-Neu

cts <- readRDS("cts.rds")
cts <- cts[rowMeans(cts)>1,]
colnames(cts)
cts <- cts[,c(1:3,7,8,9,13,14,15)]

condition <- factor(c("Neu_control","Neu_control","Neu_control",
                      "Neu_TAN","Neu_TAN","Neu_TAN",'Neu_TAN','Neu_TAN','Neu_TAN'))
colData <- data.frame(row.names=colnames(cts), condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds1 <- DESeq(dds)

resultsNames(dds1)
res <- results(dds1)
summary(res)

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1_up <- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]
res1_up$ENSEMBL <- rownames(res1_up)
res1_up <- res1_up[,c(7,1:6)]

s2e <- bitr(rownames(res1_up),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)

res1_up <- left_join(res1_up, s2e, by = "ENSEMBL")
res1_up <- res1_up[!duplicated(res1_up[[1]]), ]

write.csv(res1_up, "Neu_TAN_upRegulated_Genes.csv")

res1_down <- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]
res1_down$ENSEMBL <- rownames(res1_down)
res1_down <- res1_down[,c(7,1:6)]

s2e <- bitr(rownames(res1_down),
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)

res1_down <- left_join(res1_down, s2e, by = "ENSEMBL")
res1_down <- res1_down[!duplicated(res1_down[[1]]), ]

write.csv(res1_down, "Neu_TAN_downRegulated_Genes.csv")

