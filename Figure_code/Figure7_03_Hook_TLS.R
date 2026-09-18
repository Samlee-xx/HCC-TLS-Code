###### HookNet-TLS analysis in Cohort 7 ######

library(ggplot2)


###### files ######

hook.file <- paste0(
  'F:/TLS DATA/06-ICI-TLSProfiler/05-HookTLS/',
  'batch_output/beijing/cohort7_hooknet_response/',
  'final_supplementary_table/',
  'Supplementary_Table_HookNet_TLS_Cohort7_final.csv'
)

output.dir <- paste0(
  'F:/TLS Material/06-Code/result/Figure7_TLSProfiler/',
  '03-HookNet_TLS'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### read the finalized 120-patient table ######

hook <- read.csv(hook.file, check.names = FALSE)

dim(hook)
head(hook)
length(unique(hook$`Patient ID`))
table(hook$`Response group`, useNA = 'ifany')
table(hook$`HookNet-TLS category`, useNA = 'ifany')

# This table is already the finalized Cohort 7 table.
# No patient is removed in this script.


###### patient-level HookNet-TLS grouping ######

hook$HookNet_group <- hook$`HookNet-TLS category`
hook$HookNet_group[hook$HookNet_group == 'GC-negative TLS'] <- 'GC-negative'
hook$HookNet_group[hook$HookNet_group == 'GC-positive TLS'] <- 'GC-positive'

hook$HookNet_group <- factor(
  hook$HookNet_group,
  levels = c('non-TLS', 'GC-negative', 'GC-positive')
)

hook$`Response group` <- factor(
  hook$`Response group`,
  levels = c('R', 'NR')
)

count.table <- table(hook$HookNet_group, hook$`Response group`)
count.table

count.result <- as.data.frame.matrix(count.table)
count.result$HookNet_group <- rownames(count.result)
rownames(count.result) <- NULL
count.result <- count.result[, c('HookNet_group', 'R', 'NR')]

write.csv(
  count.result,
  file.path(output.dir, 'HookNet_TLS_response_counts.csv'),
  row.names = FALSE
)


###### GC-negative versus non-TLS ######

test.table <- matrix(
  c(
    count.table['GC-negative', 'R'],
    count.table['non-TLS', 'R'],
    count.table['GC-negative', 'NR'],
    count.table['non-TLS', 'NR']
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Response = c('R', 'NR'),
    HookNet_group = c('GC-negative', 'non-TLS')
  )
)
test.table

fisher.result <- fisher.test(test.table, alternative = 'two.sided')
fisher.result

test.result <- data.frame(
  comparison = 'GC-negative versus non-TLS',
  non_TLS_R = count.table['non-TLS', 'R'],
  non_TLS_NR = count.table['non-TLS', 'NR'],
  GC_negative_R = count.table['GC-negative', 'R'],
  GC_negative_NR = count.table['GC-negative', 'NR'],
  odds_ratio_for_R_GC_negative_vs_non_TLS = unname(fisher.result$estimate),
  fisher_P_two_sided = fisher.result$p.value
)

test.result

write.csv(
  test.result,
  file.path(output.dir, 'HookNet_TLS_Fisher_test.csv'),
  row.names = FALSE
)


