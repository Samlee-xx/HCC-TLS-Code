###### Cohort 4 univariable and multivariable logistic regression ######

library(dplyr)
library(ggplot2)
library(broom)


###### files ######

clinical.file <- paste0(
  'F:/TLS DATA/05-ICI-CODEX/03-analysis/05-multi-factor analysis/',
  '02_NR_R_multivariable_analysis/08_TLS_burden_total_immune_model/',
  '01_tables/NR_R_TLS_burden_total_immune_model_input_ascii.csv'
)

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
  '02-logistic_regression'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### read the 103-patient cohort ######

dat <- read.csv(clinical.file, check.names = FALSE)
tls <- read.csv(tls.file, check.names = FALSE)
response <- read.csv(response.file, check.names = FALSE)

dim(dat)
table(dat$Response, useNA = 'ifany')
table(dat$tls_subtype, useNA = 'ifany')


response.map <- response[, c(1, 3)]
colnames(response.map) <- c('FOV', 'patient_id')
response.map$patient_id <- as.character(response.map$patient_id)
response.map$fov_key <- paste0('TLS-', response.map$FOV)
  
tls$fov_key <- tls$fov # check and confirm again

tls <- left_join(
  tls,
  response.map[, c('fov_key', 'patient_id')],
  by = 'fov_key'
)

table(is.na(tls$patient_id))
table(tls$TLS_subtype)

patient.tls <- tls %>%
  group_by(patient_id) %>%
  summarise(
    tls_count = n(),
    bf_tls_count = sum(TLS_subtype == 'BF-TLS'),
    tnc_tls_count = sum(TLS_subtype == 'TNC-TLS'),
    .groups = 'drop'
  )

dat$patient_id <- as.character(dat$patient_id)
dat <- dat[, !(colnames(dat) %in% c(
  'tls_subtype', 'tls_count', 'bf_tls_count', 'tnc_tls_count'
))]
dat <- left_join(dat, patient.tls, by = 'patient_id')

dat$tls_count[is.na(dat$tls_count)] <- 0
dat$bf_tls_count[is.na(dat$bf_tls_count)] <- 0
dat$tnc_tls_count[is.na(dat$tnc_tls_count)] <- 0

dat$tls_subtype <- 'Non-TLS'
dat$tls_subtype[dat$tnc_tls_count > 0] <- 'TNC-TLS'
dat$tls_subtype[dat$bf_tls_count > 0] <- 'BF-TLS'

table(dat$tls_subtype, dat$Response)
table(dat$tls_count)


###### variables ######

dat$responder <- as.numeric(dat$responder)
dat$tls_subtype <- factor(
  dat$tls_subtype,
  levels = c('Non-TLS', 'TNC-TLS', 'BF-TLS')
)
dat$sex2 <- factor(dat$sex2, levels = c('Male', 'Female'))
dat$treatment2 <- factor(
  dat$treatment2,
  levels = c('Single ICI', 'Combination')
)
dat$bclc_c <- factor(dat$bclc_c, levels = c('A/B', 'C'))
dat$vascular2 <- factor(dat$vascular2, levels = c('No', 'Yes'))
dat$multifocal_primary <- factor(
  dat$multifocal_primary,
  levels = c('Single/non-multiple', 'Multiple')
)

dat$additional_tls_count <- pmax(dat$tls_count - 1, 0)

write.csv(
  dat,
  file.path(output.dir, 'Cohort4_103_patient_analysis_table.csv'),
  row.names = FALSE
)


###### univariable logistic regression ######

fit.tls <- glm(responder ~ tls_subtype, data = dat, family = binomial())
fit.burden <- glm(responder ~ tls_count, data = dat, family = binomial())
fit.age <- glm(responder ~ age_10, data = dat, family = binomial())
fit.sex <- glm(responder ~ sex2, data = dat, family = binomial())
fit.treatment <- glm(responder ~ treatment2, data = dat, family = binomial())
fit.bclc <- glm(responder ~ bclc_c, data = dat, family = binomial())
fit.vascular <- glm(responder ~ vascular2, data = dat, family = binomial())
fit.diameter <- glm(
  responder ~ tumor_diameter_10mm,
  data = dat,
  family = binomial()
)
fit.multifocal <- glm(
  responder ~ multifocal_primary,
  data = dat,
  family = binomial()
)
fit.immune <- glm(
  responder ~ total_immune_10,
  data = dat,
  family = binomial()
)

summary(fit.tls)
summary(fit.burden)
summary(fit.age)
summary(fit.sex)
summary(fit.treatment)
summary(fit.bclc)
summary(fit.vascular)
summary(fit.diameter)
summary(fit.multifocal)
summary(fit.immune)

uni.tls <- tidy(fit.tls)
uni.tls$model <- 'TLS subtype'
uni.tls$n_model <- nobs(fit.tls)

uni.burden <- tidy(fit.burden)
uni.burden$model <- 'TLS burden'
uni.burden$n_model <- nobs(fit.burden)

uni.age <- tidy(fit.age)
uni.age$model <- 'Age'
uni.age$n_model <- nobs(fit.age)

uni.sex <- tidy(fit.sex)
uni.sex$model <- 'Sex'
uni.sex$n_model <- nobs(fit.sex)

uni.treatment <- tidy(fit.treatment)
uni.treatment$model <- 'Treatment'
uni.treatment$n_model <- nobs(fit.treatment)

uni.bclc <- tidy(fit.bclc)
uni.bclc$model <- 'BCLC stage'
uni.bclc$n_model <- nobs(fit.bclc)

uni.vascular <- tidy(fit.vascular)
uni.vascular$model <- 'Vascular invasion'
uni.vascular$n_model <- nobs(fit.vascular)

uni.diameter <- tidy(fit.diameter)
uni.diameter$model <- 'Tumor diameter'
uni.diameter$n_model <- nobs(fit.diameter)

uni.multifocal <- tidy(fit.multifocal)
uni.multifocal$model <- 'Multifocal tumor'
uni.multifocal$n_model <- nobs(fit.multifocal)

uni.immune <- tidy(fit.immune)
uni.immune$model <- 'Total immune infiltration'
uni.immune$n_model <- nobs(fit.immune)

univariable <- bind_rows(
  uni.tls,
  uni.burden,
  uni.age,
  uni.sex,
  uni.treatment,
  uni.bclc,
  uni.vascular,
  uni.diameter,
  uni.multifocal,
  uni.immune
)

univariable <- univariable %>%
  filter(term != '(Intercept)') %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  )

write.csv(
  univariable,
  file.path(output.dir, 'Cohort4_univariable_logistic_regression.csv'),
  row.names = FALSE
)


###### multivariable logistic regression ######

model.variables <- c(
  'responder', 'tls_subtype', 'additional_tls_count', 'age_10',
  'sex2', 'treatment2', 'bclc_c', 'vascular2',
  'tumor_diameter_10mm', 'multifocal_primary', 'total_immune_10'
)

dat.complete <- dat[complete.cases(dat[, model.variables]), ]

dim(dat.complete)
table(dat.complete$responder)
table(dat.complete$tls_subtype)

fit.multi <- glm(
  responder ~ tls_subtype + additional_tls_count + age_10 +
    sex2 + treatment2 + bclc_c + vascular2 +
    tumor_diameter_10mm + multifocal_primary + total_immune_10,
  data = dat.complete,
  family = binomial()
)

summary(fit.multi)

multivariable <- tidy(fit.multi) %>%
  filter(term != '(Intercept)') %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    n_model = nobs(fit.multi)
  )

write.csv(
  multivariable,
  file.path(output.dir, 'Cohort4_multivariable_logistic_regression.csv'),
  row.names = FALSE
)


###### forest plot ######

term.label <- c(
  'tls_subtypeBF-TLS' = 'BF-TLS vs Non-TLS',
  'tls_subtypeTNC-TLS' = 'TNC-TLS vs Non-TLS',
  'additional_tls_count' = 'TLS burden, per additional TLS',
  'age_10' = 'Age, per 10 years',
  'sex2Female' = 'Female vs Male',
  'treatment2Combination' = 'Combination vs Single ICI',
  'bclc_cC' = 'BCLC C vs A/B',
  'vascular2Yes' = 'Vascular invasion: Yes vs No',
  'tumor_diameter_10mm' = 'Tumor diameter, per 10 mm',
  'multifocal_primaryMultiple' = 'Multifocal tumor: Yes vs No',
  'total_immune_10' = 'Total immune infiltration, per 10%'
)

plot.order <- c(
  'BF-TLS vs Non-TLS',
  'TNC-TLS vs Non-TLS',
  'TLS burden, per additional TLS',
  'Age, per 10 years',
  'Female vs Male',
  'Combination vs Single ICI',
  'BCLC C vs A/B',
  'Vascular invasion: Yes vs No',
  'Tumor diameter, per 10 mm',
  'Multifocal tumor: Yes vs No',
  'Total immune infiltration, per 10%'
)

multivariable$label <- term.label[multivariable$term]
multivariable$label <- factor(
  multivariable$label,
  levels = rev(plot.order)
)

forest <- ggplot(
  multivariable,
  aes(x = OR, y = label)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = 'grey60') +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = 'y',
    width = 0.15,
    colour = 'grey40'
  ) +
  geom_point(size = 2, colour = '#377EB8') +
  scale_x_log10() +
  labs(
    x = 'Covariate-adjusted odds ratio for response',
    y = NULL
  ) +
  theme_classic()

ggsave(
  file.path(output.dir, 'Cohort4_multivariable_forest.pdf'),
  forest,
  width = 7.5,
  height = 4.5
)
