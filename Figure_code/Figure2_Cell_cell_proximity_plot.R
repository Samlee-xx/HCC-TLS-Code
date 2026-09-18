###### cell-cell proximity plots ######

.libPaths('D:/software/Rlibs/R_copy')

library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)


###### read results ######

pair = read.csv(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/03-BileDuct analysis/03-result/01-tables/cell_cell_proximity_by_TLS_simple.csv',
  check.names = FALSE
)

chol.neu = read.csv(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/03-BileDuct analysis/03-result/01-tables/chol_neu_proximity_by_TLS_simple.csv',
  check.names = FALSE
)

dim(pair)
dim(chol.neu)
table(pair$tls_subtype)
table(chol.neu$tls_subtype)


###### panel E: cholangiocyte-like cell proximity network ######

celltype.order = c(
  'Endothelial', 'Fibroblast', 'HEV', 'Lymphatic', 'Macrophage',
  'Neutrophil', 'T-B', 'B', 'CD8T', 'CD4T', 'DC'
)

celltype.label = c(
  'Endothelial' = 'Endo',
  'Fibroblast' = 'Fibro',
  'HEV' = 'HEV',
  'Lymphatic' = 'LEC',
  'Macrophage' = 'Macro',
  'Neutrophil' = 'Neu',
  'T-B' = 'T-B',
  'B' = 'B',
  'CD8T' = 'CD8+T',
  'CD4T' = 'CD4+T',
  'DC' = 'DC'
)

network = pair %>%
  filter(neighbor == 'Cholangiocyte') %>%
  filter(center %in% celltype.order) %>%
  group_by(center) %>%
  summarise(
    n_TLS = sum(is.finite(z_score)),
    z_score = mean(z_score, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    center = factor(center, levels = celltype.order),
    label = celltype.label[as.character(center)]
  ) %>%
  arrange(center)

network

angles = seq(-pi / 2, pi / 2, length.out = nrow(network))

nodes = data.frame(
  name = c('Cholangiocyte', as.character(network$center)),
  label = c('Chol-like', network$label),
  x = c(0, cos(angles)),
  y = c(0, sin(angles))
)

edges = data.frame(
  from = 'Cholangiocyte',
  to = as.character(network$center),
  z_score = network$z_score
) %>%
  left_join(nodes[, c('name', 'x', 'y')], by = c('from' = 'name')) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(nodes[, c('name', 'x', 'y')], by = c('to' = 'name')) %>%
  rename(x_to = x, y_to = y)

node.color = nodes %>%
  left_join(
    network %>%
      transmute(name = as.character(center), z_score),
    by = 'name'
  ) %>%
  mutate(z_score = ifelse(is.na(z_score), 0, z_score))

p.network = ggplot() +
  geom_curve(
    data = edges,
    aes(
      x = x_from, y = y_from,
      xend = x_to, yend = y_to,
      colour = z_score
    ),
    curvature = 0.3,
    linewidth = 1.15,
    lineend = 'round'
  ) +
  geom_point(
    data = subset(node.color, name != 'Cholangiocyte'),
    aes(x = x, y = y, fill = z_score),
    shape = 21,
    size = 5.4,
    color = 'black',
    stroke = 0.45
  ) +
  geom_point(
    data = subset(node.color, name == 'Cholangiocyte'),
    aes(x = x, y = y),
    shape = 21,
    size = 4.2,
    fill = 'white',
    color = 'black'
  ) +
  geom_text(
    data = nodes,
    aes(
      x = x,
      y = y,
      label = label,
      hjust = ifelse(name == 'Cholangiocyte', 0.5, -0.25)
    ),
    vjust = ifelse(nodes$name == 'Cholangiocyte', 1.8, 0.5),
    size = 4
  ) +
  scale_colour_gradient2(
    low = '#43B0F1',
    mid = 'white',
    high = '#CA0020',
    midpoint = 0,
    limits = c(-0.1, 0.1),
    oob = squish,
    name = 'Z score'
  ) +
  scale_fill_gradient2(
    low = '#43B0F1',
    mid = 'white',
    high = '#CA0020',
    midpoint = 0,
    limits = c(-0.1, 0.1),
    oob = squish,
    guide = 'none'
  ) +
  coord_equal(xlim = c(-0.35, 1.35), ylim = c(-1.15, 1.15), clip = 'off') +
  labs(title = 'Cell-cell proximity (TA-TLSs)') +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
    legend.position = 'bottom',
    plot.margin = margin(10, 45, 10, 10)
  )

ggsave(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/03-BileDuct analysis/03-result/02-figures/04-cell-pairing/cell_cell_proximity_all_TA_TLS_simple.pdf',
  p.network,
  width = 4.5,
  height = 4.3
)


###### panel F: neutrophils close to cholangiocyte-like cells ######

chol.neu$tls_subtype = factor(
  chol.neu$tls_subtype,
  levels = c('TNC-TLS', 'BF-TLS')
)

group.color = c('TNC-TLS' = '#0073C2', 'BF-TLS' = '#EFC000')

p1 = ggboxplot(
  chol.neu,
  x = 'tls_subtype',
  y = 'percentage_of_all_cells',
  color = 'tls_subtype',
  palette = group.color,
  add = 'jitter',
  add.params = list(size = 0.9, alpha = 0.55),
  outlier.shape = NA
) +
  stat_compare_means(
    method = 'wilcox.test',
    method.args = list(exact = FALSE),
    label = 'p.format'
  ) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(
    title = 'Normalized to\nall cells',
    x = NULL,
    y = 'Percentage (%)',
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'none')

p2 = ggboxplot(
  chol.neu,
  x = 'tls_subtype',
  y = 'percentage_of_all_neutrophils',
  color = 'tls_subtype',
  palette = group.color,
  add = 'jitter',
  add.params = list(size = 0.9, alpha = 0.55),
  outlier.shape = NA
) +
  stat_compare_means(
    method = 'wilcox.test',
    method.args = list(exact = FALSE),
    label = 'p.format'
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = 'Normalized to\nall neutrophils',
    x = NULL,
    y = 'Percentage (%)',
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = 'none')

p.chol.neu = ggarrange(p1, p2, ncol = 2, widths = c(1, 1))

ggsave(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/03-BileDuct analysis/03-result/02-figures/04-cell-pairing/chol_neu_proximity_BF_TNC_TLS_simple.pdf',
  p.chol.neu,
  width = 5.6,
  height = 3.8
)
