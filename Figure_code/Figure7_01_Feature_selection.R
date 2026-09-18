###### Selection of a small and stable TLS feature set ######

library(pROC)
library(ggplot2)


###### files ######

input.dir <- paste0(
  'F:/TLS Material/05-table_methods/05-Figure7/feature_selection_slide_loso'
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
  '01-feature_selection'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### read the two CODEX cohorts ######

cohort2 <- read.csv(cohort2.file, check.names = FALSE)
cohort3 <- read.csv(cohort3.file, check.names = FALSE)

dim(cohort2)
dim(cohort3)
head(cohort2)
table(cohort2$group)
table(cohort3$group)

length(unique(cohort2$slide))
length(unique(cohort3$slide))


###### seven candidate cell-type features ######

feature.names <- c(
  'B', 'Bile_Duct', 'Endothelial', 'Macrophage',
  'Neutrophil', 'Fibroblast', 'tcell'
)

use.columns <- c('tls_id', 'slide', feature.names, 'group')

stopifnot(all(complete.cases(cohort2[, use.columns])))
stopifnot(all(complete.cases(cohort3[, use.columns])))
stopifnot(all(cohort2$group %in% c('BF-TLS', 'TNC-TLS')))
stopifnot(all(cohort3$group %in% c('BF-TLS', 'TNC-TLS')))

cohort2 <- cohort2[, use.columns]
cohort3 <- cohort3[, use.columns]

cohort2$y <- ifelse(cohort2$group == 'TNC-TLS', 1, 0)
cohort3$y <- ifelse(cohort3$group == 'TNC-TLS', 1, 0)

dim(cohort2)
dim(cohort3)
table(cohort2$group)
table(cohort3$group)


###### make all 127 non-empty feature combinations ######

feature.combinations <- list()

for (n.feature in seq_along(feature.names)) {
  feature.combinations <- c(
    feature.combinations,
    combn(feature.names, n.feature, simplify = FALSE)
  )
}

length(feature.combinations)


###### logistic regression for every combination ######

result <- data.frame()

for (i in seq_along(feature.combinations)) {

  one.combination <- feature.combinations[[i]]
  model.formula <- reformulate(one.combination, response = 'y')

  cohort2.probability <- rep(NA_real_, nrow(cohort2))

  for (held.slide in unique(cohort2$slide)) {
    test.index <- cohort2$slide == held.slide

    fit <- glm(
      model.formula,
      data = cohort2[!test.index, ],
      family = binomial()
    )

    cohort2.probability[test.index] <- predict(
      fit,
      newdata = cohort2[test.index, ],
      type = 'response'
    )
  }

  auc.cohort2 <- as.numeric(
    auc(cohort2$y, cohort2.probability, levels = c(0, 1), direction = '<')
  )

  fit.all <- glm(
    model.formula,
    data = cohort2,
    family = binomial()
  )

  cohort3.probability <- predict(
    fit.all,
    newdata = cohort3,
    type = 'response'
  )

  auc.cohort3 <- as.numeric(
    auc(cohort3$y, cohort3.probability, levels = c(0, 1), direction = '<')
  )

  result <- rbind(
    result,
    data.frame(
      feature_combination = paste(one.combination, collapse = '+'),
      number_of_features = length(one.combination),
      Cohort2_slide_LOSO_AUC = auc.cohort2,
      Cohort3_AUC = auc.cohort3,
      robust_score = min(auc.cohort2, auc.cohort3),
      mean_AUC = mean(c(auc.cohort2, auc.cohort3)),
      AUC_gap = abs(auc.cohort2 - auc.cohort3)
    )
  )
}

dim(result)
head(result)
summary(result$robust_score)


###### choose the smallest model with robust score >= 0.85 ######

performance.floor <- 0.85

eligible <- result[result$robust_score >= performance.floor, ]
eligible <- eligible[
  order(
    eligible$number_of_features,
    -eligible$robust_score,
    -eligible$mean_AUC,
    eligible$AUC_gap
  ),
]

selected <- eligible[1, ]
selected

result$selected <- result$feature_combination == selected$feature_combination

write.csv(
  result,
  file.path(output.dir, 'all_127_feature_combinations.csv'),
  row.names = FALSE
)

write.csv(
  selected,
  file.path(output.dir, 'selected_feature_combination.csv'),
  row.names = FALSE
)


###### best robust score at each model size ######

best.by.number <- data.frame()

for (n.feature in 1:length(feature.names)) {
  one.size <- result[result$number_of_features == n.feature, ]
  one.size <- one.size[order(-one.size$robust_score, -one.size$mean_AUC), ]
  best.by.number <- rbind(best.by.number, one.size[1, ])
}

best.by.number

write.csv(
  best.by.number,
  file.path(output.dir, 'best_combination_by_number_of_features.csv'),
  row.names = FALSE
)


### plOTS

set.seed(20260917)

p1 <- ggplot(
  result,
  aes(x = number_of_features, y = robust_score,
      colour = factor(number_of_features))
) +
  geom_hline(yintercept = performance.floor, linetype = 2, colour = 'grey50') +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.8) +
  geom_point(
    data = result[result$selected, ],
    colour = '#D62728', size = 3
  ) +
  scale_x_continuous(breaks = 1:7) +
  labs(x = 'Number of features', y = 'Robust score') +
  theme_classic() +
  theme(legend.position = 'none')

p2 <- ggplot(
  best.by.number,
  aes(x = number_of_features, y = robust_score)
) +
  geom_hline(yintercept = performance.floor, linetype = 2, colour = 'grey50') +
  geom_line(colour = '#2C7FB8') +
  geom_point(colour = '#2C7FB8', size = 2) +
  geom_point(
    data = best.by.number[best.by.number$selected, ],
    colour = '#D62728', size = 3
  ) +
  scale_x_continuous(breaks = 1:7) +
  labs(x = 'Number of features', y = 'Robust score') +
  theme_classic()

ggsave(
  file.path(output.dir, 'Figure7H_feature_combinations.pdf'),
  p1, width = 5.2, height = 3.8
)

ggsave(
  file.path(output.dir, 'Figure7I_best_model_by_feature_number.pdf'),
  p2, width = 4.5, height = 3.8
)


