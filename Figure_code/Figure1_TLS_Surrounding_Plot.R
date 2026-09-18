###### Figure 1I: cell composition and enrichment of TA-TLSs ######

library(dplyr)
library(ggplot2)
library(ggpubr)


###### paths and data ######

data_dir <- "F:/TLS Material/06-Code/Figure1I_source_data"

percentage <- read.csv(
  file.path(data_dir, "Figure1I_percentage_by_FOV.csv"),
  check.names = FALSE
)

percentage_stat <- read.csv(
  file.path(data_dir, "Figure1I_percentage_statistics.csv"),
  check.names = FALSE
)

enrichment <- read.csv(
  file.path(data_dir, "Figure1I_enrichment_by_FOV.csv"),
  check.names = FALSE
)

enrichment_stat <- read.csv(
  file.path(data_dir, "Figure1I_enrichment_statistics.csv"),
  check.names = FALSE
)


###### cell types displayed in Figure 1I ######

selected_celltypes <- c(
  "B", "Neutrophil", "Cholangiocyte", "Fibroblast"
)

display_names <- c(
  B = "B",
  Neutrophil = "Neutrophil",
  Cholangiocyte = "Cholangiocyte-like",
  Fibroblast = "Fibroblast"
)

group_colors <- c(
  TLS = "#46A9E0",
  `TLS-surrounding` = "#D94B55"
)

format_p <- function(x) {
  if (x < 0.001) {
    formatC(x, format = "e", digits = 2)
  } else {
    formatC(x, format = "f", digits = 4)
  }
}


###### paired cell-type percentages ######

percentage_plots <- list()

for (one_celltype in selected_celltypes) {
  one_data <- percentage %>%
    filter(celltype == one_celltype) %>%
    mutate(region = factor(region, levels = c("TLS", "TLS-surrounding")))

  one_stat <- percentage_stat %>%
    filter(celltype == one_celltype)

  stat_label <- paste0(
    "P = ", format_p(one_stat$p_value),
    "\nFDR = ", format_p(one_stat$p_adj_BH),
    "\nr_rb = ", format(round(one_stat$r_rb, 2), trim = TRUE)
  )

  one_plot <- ggplot(
    one_data,
    aes(x = region, y = percentage, fill = region)
  ) +
    geom_line(
      aes(group = fov_key),
      color = "grey75",
      linewidth = 0.25,
      alpha = 0.15
    ) +
    geom_boxplot(
      width = 0.55,
      color = "black",
      linewidth = 0.55,
      outlier.size = 0.7
    ) +
    annotate(
      "text",
      x = 1,
      y = Inf,
      label = stat_label,
      hjust = 0,
      vjust = 1.05,
      size = 3.5
    ) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.23))) +
    labs(
      title = display_names[[one_celltype]],
      x = NULL,
      y = "Percentage (%)",
      fill = "Groups"
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "top"
    )

  percentage_plots[[one_celltype]] <- one_plot
}


###### enrichment scores ######

enrichment_plots <- list()

for (one_celltype in selected_celltypes) {
  one_data <- enrichment %>%
    filter(celltype == one_celltype)

  one_stat <- enrichment_stat %>%
    filter(celltype == one_celltype)

  stat_label <- paste0(
    "P = ", format_p(one_stat$p_value_against_1),
    "\nFDR = ", format_p(one_stat$p_adj_BH),
    "\nMedian = ", formatC(
      one_stat$median_enrichment,
      format = "f",
      digits = 2
    )
  )

  one_plot <- ggplot(one_data, aes(x = "", y = enrichment_score)) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_violin(
      fill = "#B8AED1",
      color = "black",
      width = 0.9,
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.18,
      fill = "grey45",
      color = "black",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      size = 0.55,
      alpha = 0.25,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0.58,
      y = Inf,
      label = stat_label,
      hjust = 0,
      vjust = 1.05,
      size = 3.5
    ) +
    coord_cartesian(ylim = c(0, 4), clip = "off") +
    labs(
      title = display_names[[one_celltype]],
      x = NULL,
      y = "Enrichment score"
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  enrichment_plots[[one_celltype]] <- one_plot
}



