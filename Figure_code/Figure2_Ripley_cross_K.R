###### Polygon-boundary Ripley's cross-K analysis ######

library(data.table)
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)

# This analysis asks whether neutrophils are spatially enriched around
# cholangiocyte-like cells within each TA-TLS.


###### input and output files ######

cell.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '03-BileDuct analysis/01-input/',
  'dbscan_tls_cells_final.csv.gz'
)

subtype.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '03-BileDuct analysis/01-input/',
  'tls_subtype_assignment_final.csv'
)

# The TLS polygons were generated during DBSCAN TLS identification.
# They are alpha-shape polygons with alpha radius = 120 pixels.
polygon.file <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '01-TLS subtype/01-DBSCAN/03-result/01-tables/',
  'tls_polygon_final.csv'
)

result.path <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '03-BileDuct analysis/03-result/01-tables/',
  '05-ripley k analysis'
)

figure.path <- paste0(
  'F:/TLS DATA/01-CODEX/03-Analysis-final/',
  '03-BileDuct analysis/03-result/02-figures/',
  '05-ripley k analysis'
)

dir.create(result.path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure.path, recursive = TRUE, showWarnings = FALSE)

target.radius <- 100
minimum.cells <- 3


###### read data ######

subtype <- read.csv(subtype.file, check.names = FALSE)
subtype <- subset(
  subtype,
  tls_subtype %in% c('BF-TLS', 'TNC-TLS'),
  select = c(tls_id, slide, fov, tls_subtype)
)
subtype <- subtype[!duplicated(subtype$tls_id), ]

polygon <- read.csv(polygon.file, check.names = FALSE)
polygon <- subset(polygon, tls_id %in% subtype$tls_id)

cells <- fread(
  cell.file,
  select = c('tls_id', 'centroid.0', 'centroid.1', 'celltype_l2')
)
cells <- subset(
  cells,
  tls_id %in% subtype$tls_id &
    celltype_l2 %in% c('Cholangiocyte', 'Neutrophil') &
    is.finite(centroid.0) &
    is.finite(centroid.1)
)

# Confirm that the upstream TLS boundaries match the manuscript method.
stopifnot(nrow(polygon) == nrow(subtype))
stopifnot(!anyDuplicated(polygon$tls_id))
stopifnot(all(polygon$boundary_method == 'alpha_shape'))
stopifnot(all(polygon$alpha_radius == 120))


###### convert one WKT polygon to spatstat coordinates ######

read.polygon <- function(wkt) {

  text <- sub('POLYGON ((', '', wkt, fixed = TRUE)
  text <- substr(text, 1, nchar(text) - 2)
  pieces <- strsplit(text, ',', fixed = TRUE)[[1]]

  xy <- do.call(
    rbind,
    lapply(pieces, function(one.piece) {
      as.numeric(strsplit(trimws(one.piece), '[[:space:]]+')[[1]][1:2])
    })
  )

  if (!all(xy[1, ] == xy[nrow(xy), ])) {
    xy <- rbind(xy, xy[1, ])
  }

  # spatstat requires the outer polygon boundary in anticlockwise order.
  signed.area <- sum(
    xy[-nrow(xy), 1] * xy[-1, 2] -
      xy[-1, 1] * xy[-nrow(xy), 2]
  )

  if (signed.area < 0) {
    xy <- xy[nrow(xy):1, ]
  }

  list(x = xy[, 1], y = xy[, 2])
}


###### calculate cross-K for each TA-TLS ######

result <- data.frame()

for (tls.id in subtype$tls_id) {

  one.subtype <- subtype[subtype$tls_id == tls.id, ]
  one.polygon <- polygon[polygon$tls_id == tls.id, ]
  one.cells <- cells[cells$tls_id == tls.id, ]

  n.cholangiocyte <- sum(one.cells$celltype_l2 == 'Cholangiocyte')
  n.neutrophil <- sum(one.cells$celltype_l2 == 'Neutrophil')

  if (nrow(one.cells) == 0) {
    result <- rbind(
      result,
      data.frame(
        tls_id = tls.id,
        slide = one.subtype$slide,
        fov = one.subtype$fov,
        tls_subtype = one.subtype$tls_subtype,
        n_cholangiocyte = 0,
        n_neutrophil = 0,
        radius_pixel = NA,
        k_theoretical = NA,
        k_observed = NA,
        normalized_difference = NA,
        status = 'no_target_cells'
      )
    )
    next
  }

  polygon.coordinates <- read.polygon(one.polygon$polygon_wkt)
  tls.window <- owin(poly = polygon.coordinates)

  inside <- inside.owin(
    one.cells$centroid.0,
    one.cells$centroid.1,
    tls.window
  )
  one.cells <- one.cells[inside, ]

  n.cholangiocyte <- sum(one.cells$celltype_l2 == 'Cholangiocyte')
  n.neutrophil <- sum(one.cells$celltype_l2 == 'Neutrophil')

  # Exclude TLSs containing fewer than three cells of either type.
  if (n.cholangiocyte < minimum.cells || n.neutrophil < minimum.cells) {
    result <- rbind(
      result,
      data.frame(
        tls_id = tls.id,
        slide = one.subtype$slide,
        fov = one.subtype$fov,
        tls_subtype = one.subtype$tls_subtype,
        n_cholangiocyte = n.cholangiocyte,
        n_neutrophil = n.neutrophil,
        radius_pixel = NA,
        k_theoretical = NA,
        k_observed = NA,
        normalized_difference = NA,
        status = 'insufficient_cells'
      )
    )
    next
  }

  cell.mark <- factor(
    one.cells$celltype_l2,
    levels = c('Cholangiocyte', 'Neutrophil')
  )

  point.pattern <- ppp(
    x = one.cells$centroid.0,
    y = one.cells$centroid.1,
    window = tls.window,
    marks = cell.mark
  )

  cross.k <- Kcross(
    point.pattern,
    i = 'Cholangiocyte',
    j = 'Neutrophil',
    correction = 'isotropic',
    rmax = target.radius
  )

  radius.index <- which.min(abs(cross.k$r - target.radius))
  k.observed <- cross.k$iso[radius.index]
  k.theoretical <- cross.k$theo[radius.index]
  normalized.difference <-
    (k.observed - k.theoretical) / k.theoretical

  if (!is.finite(normalized.difference)) {
    status <- 'non_finite_k'
  } else {
    status <- 'ok'
  }

  result <- rbind(
    result,
    data.frame(
      tls_id = tls.id,
      slide = one.subtype$slide,
      fov = one.subtype$fov,
      tls_subtype = one.subtype$tls_subtype,
      n_cholangiocyte = n.cholangiocyte,
      n_neutrophil = n.neutrophil,
      radius_pixel = cross.k$r[radius.index],
      k_theoretical = k.theoretical,
      k_observed = k.observed,
      normalized_difference = normalized.difference,
      status = status
    )
  )
}

write.csv(
  result,
  file.path(
    result.path,
    'ripley_cross_k_chol_neu_polygon_all_TA_TLS_100px.csv'
  ),
  row.names = FALSE
)


###### TLS-level statistical test ######

valid.result <- subset(
  result,
  status == 'ok' & is.finite(normalized_difference)
)

# H0: the TLS-level normalized deviations are centred at zero.
wilcoxon.result <- wilcox.test(
  valid.result$normalized_difference,
  mu = 0,
  alternative = 'two.sided',
  exact = FALSE
)

summary.result <- data.frame(
  n_tls_input = nrow(subtype),
  n_tls_valid = nrow(valid.result),
  n_tls_excluded = nrow(subtype) - nrow(valid.result),
  n_bf_tls_valid = sum(valid.result$tls_subtype == 'BF-TLS'),
  n_tnc_tls_valid = sum(valid.result$tls_subtype == 'TNC-TLS'),
  target_radius_pixel = target.radius,
  minimum_cells_per_type = minimum.cells,
  median_normalized_difference = median(valid.result$normalized_difference),
  mean_normalized_difference = mean(valid.result$normalized_difference),
  q25_normalized_difference = quantile(
    valid.result$normalized_difference,
    0.25
  ),
  q75_normalized_difference = quantile(
    valid.result$normalized_difference,
    0.75
  ),
  n_positive = sum(valid.result$normalized_difference > 0),
  percentage_positive = mean(valid.result$normalized_difference > 0) * 100,
  wilcoxon_statistic = unname(wilcoxon.result$statistic),
  wilcoxon_p_value = wilcoxon.result$p.value
)

write.csv(
  summary.result,
  file.path(
    result.path,
    'ripley_cross_k_chol_neu_polygon_all_TA_TLS_100px_summary.csv'
  ),
  row.names = FALSE
)

status.result <- as.data.frame(table(result$status))
colnames(status.result) <- c('status', 'n_tls')

write.csv(
  status.result,
  file.path(
    result.path,
    'ripley_cross_k_chol_neu_polygon_all_TA_TLS_100px_status.csv'
  ),
  row.names = FALSE
)


###### plot the normalized difference for each valid TLS ######

plot.data <- valid.result[order(valid.result$normalized_difference), ]
plot.data$tls_id <- factor(plot.data$tls_id, levels = plot.data$tls_id)
plot.data$direction <- ifelse(
  plot.data$normalized_difference > 0,
  'Enrichment',
  'Avoidance'
)

p <- ggplot(
  plot.data,
  aes(x = tls_id, y = normalized_difference, color = direction)
) +
  geom_segment(
    aes(xend = tls_id, y = 0, yend = normalized_difference),
    color = 'grey80',
    linewidth = 0.3
  ) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.4) +
  scale_color_manual(
    values = c('Enrichment' = '#adc178', 'Avoidance' = '#ffa0c5')
  ) +
  labs(
    x = 'TA-TLS',
    y = expression((K[observed] - K[theoretical]) / K[theoretical]),
    color = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  file.path(
    figure.path,
    'ripley_cross_k_chol_neu_polygon_all_TA_TLS_100px.pdf'
  ),
  p,
  width = 6,
  height = 3
)

print(summary.result)
print(wilcoxon.result)
