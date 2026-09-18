###### Figure 2K: predicted TLS subtype across serial sections ######

library(ggplot2)


###### paths and data ######

data_dir <- "F:/TLS Material/06-Code/Serial_TLS_source_data"

data <- read.csv(
  file.path(data_dir, "serial_TLS_predictions_TLS.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


###### TLS names and section order ######

tls_names <- paste0("TLS", 1:7)

for (one_tls in tls_names) {
  one_data <- data[data$TLS == one_tls, ]
  if (
    nrow(one_data) != 17 ||
    !all(sort(one_data$section_index) == 1:17)
  ) {
    stop(paste(one_tls, "does not contain all 17 serial sections."))
  }
}

data$TLS <- factor(data$TLS, levels = rev(tls_names))
data$section <- factor(
  data$section_number,
  levels = seq(1, 49, by = 3)
)


###### probability heatmap ######

figure_2K <- ggplot(
  data,
  aes(x = section, y = TLS, fill = prob_TNC_logistic)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#F2F2F2",
    high = "#D7191C",
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("BF-TLS\n0", "0.5", "TNC-TLS\n1"),
    name = "TNC-TLS probability"
  ) +
  labs(
    title = "Stability analysis of TLS subtypes",
    x = "Section number",
    y = NULL
  ) +
  coord_fixed(ratio = 0.75) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.key.width = grid::unit(5, "cm")
  )

ggsave(
  file.path(data_dir, "Figure2K_serial_TLS_probability.pdf"),
  figure_2K,
  width = 10,
  height = 4.5
)

ggsave(
  file.path(data_dir, "Figure2K_serial_TLS_probability.png"),
  figure_2K,
  width = 10,
  height = 4.5,
  dpi = 600
)
