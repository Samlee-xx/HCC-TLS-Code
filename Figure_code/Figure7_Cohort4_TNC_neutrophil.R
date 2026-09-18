###### Neutrophil percentage in TNC-TLS: responder versus non-responder ######

library(dplyr)
library(ggplot2)


###### files ######

tls.file <- paste0(
  'F:/TLS Material/06-Code/result/Figure6_Cohort4/',
  '01-TLS_clustering/Cohort4_TLS_subtype_assignment.csv'
)

response.file <- paste0(
  'F:/TLS DATA/05-ICI-CODEX/03-analysis/04-ICI response/',
  'ICI response and TLS subtype.csv'
)

output.dir <- paste0(
  'F:/TLS Material/06-Code/result/Figure6_Cohort4/',
  '03-TNC_TLS_neutrophil'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### read data ######

tls <- read.csv(tls.file, check.names = FALSE)
response <- read.csv(response.file, check.names = FALSE)

dim(tls)
table(tls$TLS_subtype)
table(response$Response)


tnc <- tls[tls$TLS_subtype == 'TNC-TLS', ]

response.map <- response[, c(1, 3, 4)]
colnames(response.map) <- c('FOV', 'patient_id', 'Response')
response.map$patient_id <- as.character(response.map$patient_id)
response.map$fov_key <- paste0('TLS-', response.map$FOV)

tnc$fov_key <- tnc$fov
not.exact <- !(tnc$fov_key %in% response.map$fov_key)
tnc$fov_key[not.exact] <- sub('-[^-]+$', '', tnc$fov_key[not.exact])

tnc <- left_join(
  tnc,
  response.map[, c('fov_key', 'patient_id', 'Response')],
  by = 'fov_key'
)

tnc$response_group <- NA
tnc$response_group[tnc$Response == 'PR'] <- 'R'
tnc$response_group[tnc$Response %in% c('SD', 'PD')] <- 'NR'
tnc$response_group <- factor(tnc$response_group, levels = c('R', 'NR'))

dim(tnc)
table(tnc$response_group, useNA = 'ifany')
table(tnc$Response, useNA = 'ifany')

write.csv(
  tnc,
  file.path(output.dir, 'Cohort4_TNC_TLS_neutrophil_TLS_level.csv'),
  row.names = FALSE
)


###### TLS-level comparison, kept as a sensitivity analysis ######

tls.t.test <- t.test(
  dpct_Neutrophil ~ response_group,
  data = tnc
)

tls.wilcox <- wilcox.test(
  dpct_Neutrophil ~ response_group,
  data = tnc,
  exact = FALSE
)

tls.t.test
tls.wilcox


patient <- tnc %>%
  group_by(patient_id, response_group) %>%
  summarise(
    n_TNC_TLS = n(),
    neutrophil_percent = mean(dpct_Neutrophil),
    .groups = 'drop'
  )

dim(patient)
table(patient$response_group)
summary(patient$neutrophil_percent)

patient.t.test <- t.test(
  neutrophil_percent ~ response_group,
  data = patient
)

patient.wilcox <- wilcox.test(
  neutrophil_percent ~ response_group,
  data = patient,
  exact = FALSE
)

patient.t.test
patient.wilcox

write.csv(
  patient,
  file.path(output.dir, 'Cohort4_TNC_TLS_neutrophil_patient_level.csv'),
  row.names = FALSE
)

test.result <- data.frame(
  level = c('TLS', 'TLS', 'Patient', 'Patient'),
  test = c('Welch t test', 'Wilcoxon test', 'Welch t test', 'Wilcoxon test'),
  p_value = c(
    tls.t.test$p.value,
    tls.wilcox$p.value,
    patient.t.test$p.value,
    patient.wilcox$p.value
  )
)

write.csv(
  test.result,
  file.path(output.dir, 'Cohort4_TNC_TLS_neutrophil_group_tests.csv'),
  row.names = FALSE
)




plot.patient <- ggplot(
  patient,
  aes(x = response_group, y = neutrophil_percent, fill = response_group)
) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, colour = 'grey40') +
  scale_fill_manual(values = c('R' = '#B78AD6', 'NR' = '#6BAA45')) +
  labs(
    x = NULL,
    y = 'Neutrophil percentage in TNC-TLS',
    title = paste0('Patient-level P = ', signif(patient.wilcox$p.value, 3))
  ) +
  theme_classic() +
  theme(legend.position = 'none')

ggsave(
  file.path(output.dir, 'Cohort4_TNC_TLS_neutrophil_patient_level_sensitivity.pdf'),
  plot.patient,
  width = 3,
  height = 3.5
)


plot.tls <- ggplot(
  tnc,
  aes(x = response_group, y = dpct_Neutrophil, fill = response_group)
) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, colour = 'grey40') +
  scale_fill_manual(values = c('R' = '#B78AD6', 'NR' = '#6BAA45')) +
  labs(
    x = NULL,
    y = 'Neutrophil percentage in TNC-TLS',
    title = paste0('TLS-level P = ', signif(tls.t.test$p.value, 3))
  ) +
  theme_classic() +
  theme(legend.position = 'none')

ggsave(
  file.path(output.dir, 'Cohort4_TNC_TLS_neutrophil_by_response.pdf'),
  plot.tls,
  width = 3,
  height = 3.5
)
