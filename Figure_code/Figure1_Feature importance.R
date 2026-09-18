###### Cross-platform TLS subtype feature importance ######

library(dplyr)
library(ggplot2)
library(patchwork)
library(randomForest)

set.seed(777)


st.file <- paste0(
  'E:/01-TLS/01-Omics/01-ST/03-analysis final/',
  '05-cell2location/04-tables/',
  'cell2location_percent_wide_all_models.csv'
)

codex.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/02-TLS clustering/03-result/',
  'tls_subtype_celltype_composition_final.csv'
)

result.path <- paste0(
  'E:/01-TLS/01-Omics/01-ST/03-analysis final/',
  '06-feature importance/03-tables'
)

figure.path <- paste0(
  'E:/01-TLS/01-Omics/01-ST/03-analysis final/',
  '06-feature importance/01-main_figures'
)

dir.create(result.path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure.path, recursive = TRUE, showWarnings = FALSE)

###### prepare 10x Visium CytAssist data ######

st <- read.csv(st.file, check.names = FALSE)

st <- subset(
  st,
  model == 'Version2' &
    posterior == 'q05' &
    TLS_subtype %in% c('B-TLS', 'T-TLS')
)

st.data <- data.frame(
  group = factor(
    ifelse(st$TLS_subtype == 'B-TLS', 'BF-TLS', 'TNC-TLS'),
    levels = c('TNC-TLS', 'BF-TLS')
  ),
  `B-lineage` = st$`B Cell` + st$`Plasma Cell`,
  `T-lineage` = st$`CD4 T Cell` + st$`CD8 T Cell`,
  Neutrophil = st$Neutrophil,
  Fibroblast = st$Fibroblast,
  Macrophage = st$Macrophage,
  DC = st$`Dendritic Cell`,
  Endothelium = st$Endothelial,
  check.names = FALSE
)

table(st.data$group)

###### prepare CODEX data ######

codex <- read.csv(codex.file, check.names = FALSE)

codex.data <- data.frame(
  group = factor(
    codex$tls_subtype,
    levels = c('TNC-TLS', 'BF-TLS')
  ),
  `B-lineage` = codex$pct_B,
  `T-lineage` = codex$pct_CD4T + codex$pct_CD8T,
  Neutrophil = codex$pct_Neutrophil,
  DC = codex$pct_DC,
  `T-B` = codex$`pct_T-B`,
  Fibroblast = codex$pct_Fibroblast,
  Macrophage = codex$pct_Macrophage,
  Cholangiocyte = codex$pct_Cholangiocyte,
  HEV = codex$pct_HEV,
  Endothelium = codex$pct_Endothelial,
  LEC = codex$pct_Lymphatic,
  check.names = FALSE
)

table(codex.data$group)

###### random forest: 10x Visium CytAssist ######

st.features <- setdiff(colnames(st.data), 'group')

set.seed(777)
st.model <- randomForest(
  x = st.data[, st.features],
  y = st.data$group,
  ntree = 500,
  mtry = 2,
  nodesize = 20,
  maxnodes = 5,
  sampsize = round(nrow(st.data) * 0.4),
  replace = TRUE,
  importance = TRUE
)

st.importance <- importance(st.model)
st.result <- data.frame(
  Feature = rownames(st.importance),
  MeanDecreaseGini = st.importance[, 'MeanDecreaseGini'],
  check.names = FALSE
)

st.mean.TNC <- colMeans(
  st.data[st.data$group == 'TNC-TLS', st.features],
  na.rm = TRUE
)

st.mean.BF <- colMeans(
  st.data[st.data$group == 'BF-TLS', st.features],
  na.rm = TRUE
)

st.result$mean_TNC_TLS <- st.mean.TNC[st.result$Feature]
st.result$mean_BF_TLS <- st.mean.BF[st.result$Feature]
st.result$log2FC_BF_over_TNC <- log2(
  (st.result$mean_BF_TLS + 0.001) /
    (st.result$mean_TNC_TLS + 0.001)
)

st.result <- st.result[order(st.result$MeanDecreaseGini, decreasing = TRUE), ]
st.result$Gini_rank <- 1:nrow(st.result)
st.result$platform <- '10x Visium CytAssist'

###### random forest: CODEX ######

codex.features <- setdiff(colnames(codex.data), 'group')

set.seed(777)
codex.model <- randomForest(
  x = codex.data[, codex.features],
  y = codex.data$group,
  ntree = 500,
  mtry = 2,
  nodesize = 20,
  maxnodes = 5,
  sampsize = round(nrow(codex.data) * 0.4),
  replace = TRUE,
  importance = TRUE
)

codex.importance <- importance(codex.model)
codex.result <- data.frame(
  Feature = rownames(codex.importance),
  MeanDecreaseGini = codex.importance[, 'MeanDecreaseGini'],
  check.names = FALSE
)

codex.mean.TNC <- colMeans(
  codex.data[codex.data$group == 'TNC-TLS', codex.features],
  na.rm = TRUE
)

codex.mean.BF <- colMeans(
  codex.data[codex.data$group == 'BF-TLS', codex.features],
  na.rm = TRUE
)

codex.result$mean_TNC_TLS <- codex.mean.TNC[codex.result$Feature]
codex.result$mean_BF_TLS <- codex.mean.BF[codex.result$Feature]
codex.result$log2FC_BF_over_TNC <- log2(
  (codex.result$mean_BF_TLS + 0.001) /
    (codex.result$mean_TNC_TLS + 0.001)
)

codex.result <- codex.result[
  order(codex.result$MeanDecreaseGini, decreasing = TRUE),
]
codex.result$Gini_rank <- 1:nrow(codex.result)
codex.result$platform <- 'CODEX'

###### save feature importance and log2FC ######

feature.importance <- rbind(st.result, codex.result)

write.csv(
  feature.importance,
  file.path(result.path, 'feature_importance_primary_common_RF.csv'),
  row.names = FALSE
)

st.result
codex.result

###### Figure K ######

st.top3 <- st.result$Feature[1:3]
codex.top3 <- codex.result$Feature[1:3]

st.top3
codex.top3

st.result$Feature_label <- ifelse(
  st.result$Feature %in% st.top3,
  paste0(st.result$Feature, ' *'),
  st.result$Feature
)

codex.result$Feature_label <- ifelse(
  codex.result$Feature %in% codex.top3,
  paste0(codex.result$Feature, ' *'),
  codex.result$Feature
)

st.result$Feature_label <- factor(
  st.result$Feature_label,
  levels = rev(st.result$Feature_label)
)

codex.result$Feature_label <- factor(
  codex.result$Feature_label,
  levels = rev(codex.result$Feature_label)
)

st.color.limit <- max(abs(st.result$log2FC_BF_over_TNC))
codex.color.limit <- max(abs(codex.result$log2FC_BF_over_TNC))

plot.st <- ggplot(
  st.result,
  aes(x = MeanDecreaseGini, y = Feature_label, fill = log2FC_BF_over_TNC)
) +
  geom_col(color = 'grey30', width = 0.78) +
  geom_text(
    aes(label = sprintf('%.2f', MeanDecreaseGini)),
    hjust = -0.12,
    size = 3.8
  ) +
  scale_fill_gradient2(
    low = '#43A9E6',
    mid = 'white',
    high = '#E52B32',
    midpoint = 0,
    limits = c(-st.color.limit, st.color.limit),
    name = expression(log[2] * FC~'(BF-TLS / TNC-TLS)')
  ) +
  coord_cartesian(
    xlim = c(0, max(st.result$MeanDecreaseGini) * 1.3),
    clip = 'off'
  ) +
  labs(
    title = 'Gini (10x Visium CytAssist)',
    subtitle = '* Within-platform top 3',
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'bottom'
  )

plot.codex <- ggplot(
  codex.result,
  aes(x = MeanDecreaseGini, y = Feature_label, fill = log2FC_BF_over_TNC)
) +
  geom_col(color = 'grey30', width = 0.78) +
  geom_text(
    aes(label = sprintf('%.2f', MeanDecreaseGini)),
    hjust = -0.12,
    size = 3.8
  ) +
  scale_fill_gradient2(
    low = '#43A9E6',
    mid = 'white',
    high = '#E52B32',
    midpoint = 0,
    limits = c(-codex.color.limit, codex.color.limit),
    name = expression(log[2] * FC~'(BF-TLS / TNC-TLS)')
  ) +
  coord_cartesian(
    xlim = c(0, max(codex.result$MeanDecreaseGini) * 1.3),
    clip = 'off'
  ) +
  labs(
    title = 'Gini (CODEX)',
    subtitle = '* Within-platform top 3',
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'bottom'
  )

plot.combined <- plot.st + plot.codex

ggsave(
  file.path(figure.path, 'Feature_importance_10x_and_CODEX_common_RF.pdf'),
  plot.combined,
  width = 10.2,
  height = 6.1
)

ggsave(
  'F:/TLS Material/06-Code/Feature_importance_10x_and_CODEX_common_RF.png',
  plot.combined,
  width = 10.2,
  height = 6.1,
  dpi = 400
)
