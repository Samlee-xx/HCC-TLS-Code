###### Prognostic analysis of TLS gene signatures in TCGA-LIHC ######

library(IOBR)
library(dplyr)
library(survival)
library(survminer)

###### TLS subtype gene signatures ######

marker <- readRDS(
  'E:/01-TLS/01-Omics/01-ST/03-analysis final/02-DEGs/DEGs.rds'
)

marker <- marker[marker$p_val_adj <= 0.05, ]

tls.signature <- list(
  Type1_B_TLS = marker$gene[marker$cluster == 'TLS.2'],
  Type2_T_TLS = marker$gene[marker$cluster == 'TLS.0']
)

length(tls.signature$Type1_B_TLS)
length(tls.signature$Type2_T_TLS)

###### calculate TLS scores ######

tcga.tpm <- readRDS(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/02-TCGA/total_tpm_408.rds'
)

rownames(tcga.tpm) <- tcga.tpm$gene_name
tcga.tpm <- tcga.tpm[, -1]

dim(tcga.tpm)

tcga.score <- calculate_sig_score(
  pdata = NULL,
  eset = tcga.tpm,
  signature = tls.signature,
  method = 'ssgsea',
  mini_gene_count = 0,
  adjust_eset = TRUE
)

tcga.score <- as.data.frame(tcga.score)
tcga.score$Index <- NULL
colnames(tcga.score)[colnames(tcga.score) == 'ID'] <- 'Sample.ID'

dim(tcga.score)
head(tcga.score)

###### clinical information ######

sample.sheet <- read.delim(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/02-TCGA/gdc_sample_sheet.2022-12-03.tsv',
  quote = ''
)

clinical <- read.delim(
  'E:/01-TLS/01-Omics/01-ST/01-Previous version/02-TCGA/clinical.tsv',
  quote = ''
)

sample.sheet <- sample.sheet[, c(
  'Sample.ID', 'Sample.Type', 'Case.ID'
)]

clinical <- clinical[, c(
  'case_submitter_id',
  'days_to_death',
  'days_to_last_follow_up',
  'vital_status',
  'ajcc_pathologic_stage'
)]

# The downloaded clinical table contains two identical rows for each patient.
clinical <- unique(clinical)

colnames(clinical)[1] <- 'Case.ID'

dim(clinical)

tcga <- left_join(tcga.score, sample.sheet, by = 'Sample.ID')
tcga <- left_join(tcga, clinical, by = 'Case.ID')

dim(tcga)
table(tcga$Sample.Type, useNA = 'ifany')
table(tcga$vital_status, useNA = 'ifany')

tcga$days_to_death[tcga$days_to_death == "'--"] <- NA
tcga$days_to_last_follow_up[tcga$days_to_last_follow_up == "'--"] <- NA

tcga$days_to_death <- as.numeric(tcga$days_to_death)
tcga$days_to_last_follow_up <- as.numeric(tcga$days_to_last_follow_up)

tcga$event <- NA
tcga$event[tcga$vital_status == 'Alive'] <- 0
tcga$event[tcga$vital_status == 'Dead'] <- 1

tcga$OS_time <- tcga$days_to_last_follow_up
dead <- which(tcga$event == 1)
tcga$OS_time[dead] <- tcga$days_to_death[dead]

tcga <- tcga[
  !is.na(tcga$Sample.Type) &
    tcga$Sample.Type == 'Primary Tumor' &
    !is.na(tcga$OS_time) &
    !is.na(tcga$event) &
    !is.na(tcga$Type1_B_TLS) &
    !is.na(tcga$Type2_T_TLS),
]

dim(tcga)
table(tcga$event)
table(tcga$ajcc_pathologic_stage, useNA = 'ifany')

###### Type 1 B-TLS ######

cut.b <- surv_cutpoint(
  tcga,
  time = 'OS_time',
  event = 'event',
  variables = 'Type1_B_TLS',
  minprop = 0.45
)

cut.b

b.value <- cut.b$cutpoint$cutpoint[1]
tcga$B_TLS_group <- ifelse(tcga$Type1_B_TLS <= b.value, 'Low', 'High')
tcga$B_TLS_group <- factor(tcga$B_TLS_group, levels = c('Low', 'High'))

table(tcga$B_TLS_group)

fit.b <- survfit(Surv(OS_time, event) ~ B_TLS_group, data = tcga)
cox.b <- coxph(Surv(OS_time, event) ~ B_TLS_group, data = tcga)

summary(cox.b)

plot.b <- ggsurvplot(
  fit.b,
  data = tcga,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c('#1F77B4', '#D62728'),
  ggtheme = theme_classic(),
  title = 'Type 1 B-TLS'
)

print(plot.b)

###### Type 2 T-TLS ######

cut.t <- surv_cutpoint(
  tcga,
  time = 'OS_time',
  event = 'event',
  variables = 'Type2_T_TLS',
  minprop = 0.45
)

cut.t

t.value <- cut.t$cutpoint$cutpoint[1]
tcga$T_TLS_group <- ifelse(tcga$Type2_T_TLS <= t.value, 'Low', 'High')
tcga$T_TLS_group <- factor(tcga$T_TLS_group, levels = c('Low', 'High'))

table(tcga$T_TLS_group)

fit.t <- survfit(Surv(OS_time, event) ~ T_TLS_group, data = tcga)
cox.t <- coxph(Surv(OS_time, event) ~ T_TLS_group, data = tcga)

summary(cox.t)

plot.t <- ggsurvplot(
  fit.t,
  data = tcga,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c('#1F77B4', '#D62728'),
  ggtheme = theme_classic(),
  title = 'Type 2 T-TLS'
)

print(plot.t)
