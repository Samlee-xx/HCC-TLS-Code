###### Establishment of the final TLSProfiler model ######

library(pROC)


###### files ######

input.dir <- paste0(
  'F:/TLS Material/05-table_methods/05-Figure7/',
  '06-reproducibility_bundle_20260913/01-runnable_workflow/',
  'results/feature_selection_slide_loso'
)

cohort2.file <- file.path(
  input.dir,
  'Cohort1_training_final_with_slide.csv'
)

cohort3.file <- file.path(
  input.dir,
  'Cohort2_validation_final_with_slide.csv'
)

output.dir <- paste0(
  'F:/TLS Material/06-Code/result/Figure7_TLSProfiler/',
  '02-final_model'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### read data ######

cohort2 <- read.csv(cohort2.file, check.names = FALSE)
cohort3 <- read.csv(cohort3.file, check.names = FALSE)

required.columns <- c(
  'tls_id', 'slide', 'B', 'tcell', 'Neutrophil', 'group'
)

stopifnot(all(complete.cases(cohort2[, required.columns])))
stopifnot(all(complete.cases(cohort3[, required.columns])))
stopifnot(all(cohort2$group %in% c('BF-TLS', 'TNC-TLS')))
stopifnot(all(cohort3$group %in% c('BF-TLS', 'TNC-TLS')))

cohort2$y <- ifelse(cohort2$group == 'TNC-TLS', 1, 0)
cohort3$y <- ifelse(cohort3$group == 'TNC-TLS', 1, 0)

dim(cohort2)
dim(cohort3)
table(cohort2$group)
table(cohort3$group)
length(unique(cohort2$slide))
length(unique(cohort3$slide))
head(cohort2[, c('tls_id', 'slide', 'B', 'tcell', 'Neutrophil', 'group')])


###### log1p transformation of the selected percentages ######

# B, tcell and Neutrophil are percentages on the 0-100 scale.
cohort2$log_B <- log1p(cohort2$B)
cohort2$log_T <- log1p(cohort2$tcell)
cohort2$log_Neutrophil <- log1p(cohort2$Neutrophil)

cohort3$log_B <- log1p(cohort3$B)
cohort3$log_T <- log1p(cohort3$tcell)
cohort3$log_Neutrophil <- log1p(cohort3$Neutrophil)


###### Cohort 2 slide-level leave-one-out prediction ######

model.formula <- y ~ log_B + log_T + log_Neutrophil
cohort2$probability_TNC <- NA_real_

for (held.slide in unique(cohort2$slide)) {

  test.index <- cohort2$slide == held.slide

  fit.one.fold <- glm(
    model.formula,
    data = cohort2[!test.index, ],
    family = binomial()
  )

  cohort2$probability_TNC[test.index] <- predict(
    fit.one.fold,
    newdata = cohort2[test.index, ],
    type = 'response'
  )
}

summary(cohort2$probability_TNC)


###### train the locked model with all Cohort 2 TLSs ######

tlsprofiler <- glm(
  model.formula,
  data = cohort2,
  family = binomial()
)

summary(tlsprofiler)
coef(tlsprofiler)


###### determine the cutoff from Cohort 2 LOSO predictions ######

# The cutoff is the observed probability giving the largest Youden index.
candidate.cutoff <- sort(unique(cohort2$probability_TNC))
youden.index <- rep(NA_real_, length(candidate.cutoff))

for (i in seq_along(candidate.cutoff)) {

  predicted <- ifelse(
    cohort2$probability_TNC >= candidate.cutoff[i],
    1,
    0
  )

  sensitivity <- mean(predicted[cohort2$y == 1] == 1)
  specificity <- mean(predicted[cohort2$y == 0] == 0)

  youden.index[i] <- sensitivity + specificity - 1
}

decision.cutoff <- candidate.cutoff[which.max(youden.index)]
decision.cutoff

cohort2$predicted_TLS_subtype <- ifelse(
  cohort2$probability_TNC >= decision.cutoff,
  'TNC-TLS',
  'BF-TLS'
)

table(cohort2$group, cohort2$predicted_TLS_subtype)


###### apply the locked model to Cohort 3 without refitting ######

cohort3$probability_TNC <- predict(
  tlsprofiler,
  newdata = cohort3,
  type = 'response'
)

cohort3$predicted_TLS_subtype <- ifelse(
  cohort3$probability_TNC >= decision.cutoff,
  'TNC-TLS',
  'BF-TLS'
)

table(cohort3$group, cohort3$predicted_TLS_subtype)


###### ROC analysis ######

roc.cohort2 <- roc(
  cohort2$y,
  cohort2$probability_TNC,
  levels = c(0, 1),
  direction = '<',
  quiet = TRUE
)

roc.cohort3 <- roc(
  cohort3$y,
  cohort3$probability_TNC,
  levels = c(0, 1),
  direction = '<',
  quiet = TRUE
)

auc.cohort2 <- as.numeric(auc(roc.cohort2))
auc.cohort3 <- as.numeric(auc(roc.cohort3))

auc.cohort2
auc.cohort3


###### slide-cluster bootstrap confidence intervals ######

bootstrap.number <- 5000

bootstrap.auc.cohort2 <- rep(NA_real_, bootstrap.number)
slides.cohort2 <- unique(cohort2$slide)

for (i in 1:bootstrap.number) {
  sampled.slides <- sample(
    slides.cohort2,
    length(slides.cohort2),
    replace = TRUE
  )

  sampled.index <- unlist(
    lapply(sampled.slides, function(x) which(cohort2$slide == x))
  )

  sampled.y <- cohort2$y[sampled.index]
  sampled.probability <- cohort2$probability_TNC[sampled.index]

  if (length(unique(sampled.y)) == 2) {
    bootstrap.auc.cohort2[i] <- as.numeric(
      auc(
        sampled.y,
        sampled.probability,
        levels = c(0, 1),
        direction = '<',
        quiet = TRUE
      )
    )
  }
}

bootstrap.auc.cohort3 <- rep(NA_real_, bootstrap.number)
slides.cohort3 <- unique(cohort3$slide)

for (i in 1:bootstrap.number) {
  sampled.slides <- sample(
    slides.cohort3,
    length(slides.cohort3),
    replace = TRUE
  )

  sampled.index <- unlist(
    lapply(sampled.slides, function(x) which(cohort3$slide == x))
  )

  sampled.y <- cohort3$y[sampled.index]
  sampled.probability <- cohort3$probability_TNC[sampled.index]

  if (length(unique(sampled.y)) == 2) {
    bootstrap.auc.cohort3[i] <- as.numeric(
      auc(
        sampled.y,
        sampled.probability,
        levels = c(0, 1),
        direction = '<',
        quiet = TRUE
      )
    )
  }
}

ci.cohort2 <- quantile(
  bootstrap.auc.cohort2,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci.cohort3 <- quantile(
  bootstrap.auc.cohort3,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci.cohort2
ci.cohort3


###### save the final model and source data ######

model.object <- list(
  model = tlsprofiler,
  input_features = c('B_pct', 'T_pct', 'Neutrophil_pct'),
  input_scale = 'percentage, 0-100',
  transformation = 'log1p',
  positive_class = 'TNC-TLS',
  probability_cutoff = decision.cutoff,
  clinical_minimum_cells = 300
)

saveRDS(
  model.object,
  file.path(output.dir, 'TLSProfiler_log1p_BTN_model.rds')
)

coefficient.table <- data.frame(
  term = names(coef(tlsprofiler)),
  coefficient = as.numeric(coef(tlsprofiler))
)

roc.summary <- data.frame(
  cohort = c('Cohort 2', 'Cohort 3'),
  evaluation = c('slide-LOSO', 'cross-cohort'),
  n_TLS = c(nrow(cohort2), nrow(cohort3)),
  n_slides = c(length(unique(cohort2$slide)), length(unique(cohort3$slide))),
  AUC = c(auc.cohort2, auc.cohort3),
  CI_low = c(ci.cohort2[1], ci.cohort3[1]),
  CI_high = c(ci.cohort2[2], ci.cohort3[2])
)

coefficient.table
roc.summary

write.csv(
  coefficient.table,
  file.path(output.dir, 'TLSProfiler_coefficients.csv'),
  row.names = FALSE
)

write.csv(
  data.frame(probability_cutoff = decision.cutoff),
  file.path(output.dir, 'TLSProfiler_cutoff.csv'),
  row.names = FALSE
)

write.csv(
  roc.summary,
  file.path(output.dir, 'TLSProfiler_ROC_summary.csv'),
  row.names = FALSE
)

write.csv(
  cohort2[, c(
    'tls_id', 'slide', 'group', 'B', 'tcell', 'Neutrophil',
    'probability_TNC', 'predicted_TLS_subtype'
  )],
  file.path(output.dir, 'Cohort2_slide_LOSO_predictions.csv'),
  row.names = FALSE
)

write.csv(
  cohort3[, c(
    'tls_id', 'slide', 'group', 'B', 'tcell', 'Neutrophil',
    'probability_TNC', 'predicted_TLS_subtype'
  )],
  file.path(output.dir, 'Cohort3_predictions.csv'),
  row.names = FALSE
)


###### ROC plot ######

pdf(
  file.path(output.dir, 'Figure7K_TLSProfiler_ROC.pdf'),
  width = 4.5,
  height = 4.5
)

plot(
  roc.cohort3,
  col = '#2C7FB8',
  lwd = 2,
  legacy.axes = TRUE,
  main = 'TLSProfiler cross-cohort evaluation'
)
abline(a = 0, b = 1, col = 'grey70')
legend(
  'bottomright',
  legend = paste0(
    'Cohort 3 AUC = ',
    sprintf('%.3f', auc.cohort3),
    ' (', sprintf('%.3f', ci.cohort3[1]),
    '-', sprintf('%.3f', ci.cohort3[2]), ')'
  ),
  col = '#2C7FB8',
  lwd = 2,
  bty = 'n'
)

dev.off()


