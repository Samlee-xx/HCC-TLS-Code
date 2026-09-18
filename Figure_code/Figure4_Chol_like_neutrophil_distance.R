###### Neutrophil-to-cholangiocyte-like nearest-neighbor distance ###

library(data.table)
library(ggplot2)


##### paths ####

data.path <- paste0(
  'F:/TLS DATA/03-CODEX-ref_to_HD/03-analysis-final/',
  '06-Revision/01-close far tumor/01-data'
)

result.path <- paste0(
  'F:/TLS Material/06-Code/result/',
  'Figure4_Chol_like_neutrophil_distance'
)

dir.create(result.path, recursive = TRUE, showWarnings = FALSE)


#### fluorescence thresholds ###

# Read the PanCK and MPO thresholds used for each CODEX chip.
threshold <- fread(file.path(data.path, 'threshodd.txt'))
setnames(threshold, c('chip', 'marker', 'threshold'))

threshold$chip <- tolower(threshold$chip)
threshold$marker <- toupper(threshold$marker)
threshold <- threshold[chip %in% paste0('chip', 1:4)]

panck.threshold <- setNames(
  threshold[marker == 'PANCK']$threshold,
  threshold[marker == 'PANCK']$chip
)

mpo.threshold <- setNames(
  threshold[marker == 'MPO']$threshold,
  threshold[marker == 'MPO']$chip
)


###### read CODEX cell data ######

chip.files <- paste0('Chip', 1:4, '.txt')
cell.list <- vector('list', length(chip.files))

for (i in seq_along(chip.files)) {

  input.file <- file.path(data.path, chip.files[i])
  chip.name <- tolower(tools::file_path_sans_ext(chip.files[i]))

  column.names <- names(fread(input.file, nrows = 0))
  x.column <- grep('^Centroid X', column.names, value = TRUE)[1]
  y.column <- grep('^Centroid Y', column.names, value = TRUE)[1]

  required.columns <- c(
    'Parent',
    x.column,
    y.column,
    'Cell: Pan-Cytokeratin: Mean',
    'Cell: MPO: Mean'
  )

  one.chip <- fread(
    input.file,
    select = required.columns,
    na.strings = c('', 'NA', 'NaN')
  )

  setnames(
    one.chip,
    required.columns,
    c('parent', 'x', 'y', 'panck', 'mpo')
  )

  one.chip <- one.chip[
    !is.na(parent) & parent != '' &
      is.finite(x) & is.finite(y) &
      is.finite(panck) & is.finite(mpo)
  ]

  one.chip[, chip := chip.name]
  one.chip[, sample := sub('-(Close|Far)-.*$', '', parent, ignore.case = TRUE)]
  one.chip[, condition := fifelse(
    grepl('-Close-', parent, ignore.case = TRUE),
    'Close',
    fifelse(grepl('-Far-', parent, ignore.case = TRUE), 'Far', NA_character_)
  )]

  if (any(is.na(one.chip$condition))) {
    stop('Close/Far information could not be parsed from Parent in ', chip.files[i])
  }

  one.chip[, panck.positive := panck >= panck.threshold[chip.name]]
  one.chip[, mpo.positive := mpo >= mpo.threshold[chip.name]]
  one.chip[, chol.like := panck.positive & !mpo.positive]
  one.chip[, neutrophil := mpo.positive & !panck.positive]

  cell.list[[i]] <- one.chip
}

cell.data <- rbindlist(cell.list, use.names = TRUE)


###### nearest-neighbor distance within each FOV ######

nearest.distance <- function(neutrophil.xy, chol.like.xy) {

  vapply(
    seq_len(nrow(neutrophil.xy)),
    function(i) {
      min(sqrt(
        (neutrophil.xy[i, 1] - chol.like.xy[, 1])^2 +
          (neutrophil.xy[i, 2] - chol.like.xy[, 2])^2
      ))
    },
    numeric(1)
  )
}

fov.list <- split(
  cell.data,
  paste(cell.data$chip, cell.data$parent, sep = '__')
)

fov.result <- rbindlist(lapply(fov.list, function(one.fov) {

  chol.like.xy <- as.matrix(one.fov[chol.like == TRUE, .(x, y)])
  neutrophil.xy <- as.matrix(one.fov[neutrophil == TRUE, .(x, y)])

  distance.um <- nearest.distance(neutrophil.xy, chol.like.xy)

  data.table(
    chip = one.fov$chip[1],
    sample = one.fov$sample[1],
    condition = one.fov$condition[1],
    fov = one.fov$parent[1],
    n_chol_like = nrow(chol.like.xy),
    n_neutrophil = nrow(neutrophil.xy),
    median_nearest_distance_um = median(distance.um)
  )
}))


###### average FOV values within each biological sample ######

sample.result <- fov.result[
  , .(
    n_fovs = .N,
    nearest_neighbor_distance_um = mean(median_nearest_distance_um)
  ),
  by = .(chip, sample, condition)
]

paired.result <- dcast(
  sample.result,
  chip + sample ~ condition,
  value.var = 'nearest_neighbor_distance_um'
)

paired.result <- paired.result[complete.cases(Far, Close)]


###### paired statistical test ######

wilcox.result <- wilcox.test(
  paired.result$Close,
  paired.result$Far,
  paired = TRUE,
  alternative = 'two.sided',
  exact = TRUE,
  correct = FALSE
)

statistics <- data.frame(
  test = 'Two-sided paired Wilcoxon signed-rank test',
  n_pairs = nrow(paired.result),
  distal_mean_um = mean(paired.result$Far),
  distal_median_um = median(paired.result$Far),
  proximal_mean_um = mean(paired.result$Close),
  proximal_median_um = median(paired.result$Close),
  mean_difference_proximal_minus_distal_um = mean(
    paired.result$Close - paired.result$Far
  ),
  p_value = wilcox.result$p.value
)


###### paired plot ######

plot.data <- melt(
  paired.result,
  id.vars = c('chip', 'sample'),
  measure.vars = c('Far', 'Close'),
  variable.name = 'condition',
  value.name = 'distance_um'
)

plot.data[, pair_id := paste(chip, sample, sep = '_')]
plot.data[, group := factor(
  condition,
  levels = c('Far', 'Close'),
  labels = c('Distal', 'Proximal')
)]

p.label <- if (wilcox.result$p.value < 0.0001) {
  'P < 0.0001'
} else {
  paste0('P = ', formatC(wilcox.result$p.value, format = 'f', digits = 4))
}

y.maximum <- max(plot.data$distance_um) * 1.30

p <- ggplot(
  plot.data,
  aes(x = group, y = distance_um, group = pair_id)
) +
  geom_line(color = 'grey65', linewidth = 0.7) +
  geom_point(aes(color = group), size = 3) +
  stat_summary(
    aes(group = group),
    fun = median,
    geom = 'crossbar',
    width = 0.28,
    linewidth = 0.8,
    color = 'black'
  ) +
  annotate(
    'text',
    x = 1.02,
    y = y.maximum * 0.94,
    label = paste('Paired Wilcoxon', p.label, sep = '\n'),
    hjust = 0,
    vjust = 1,
    size = 5.2
  ) +
  scale_color_manual(values = c('Distal' = '#4C78A8', 'Proximal' = '#E45756')) +
  scale_y_continuous(
    limits = c(0, y.maximum),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = 'Neutrophil close to Chol-like',
    x = NULL,
    y = expression(paste('Nearest-neighbor distance (', mu, 'm)'))
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.y = element_text(size = 17),
    axis.text = element_text(color = 'black'),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

ggsave(
  file.path(result.path, 'Neutrophil_Chol_like_nearest_distance.pdf'),
  p,
  width = 4.2,
  height = 5.0
)


