
###### Identify real tls ######
dbscan.out <- read.csv('...\\TLS FOV DBSCAN.csv')
dbscan.out$tls.label <- paste0(dbscan.out$fov,'_',dbscan.out$cluster_label)
unique(dbscan.out$tls.label)

tls.label = c('ICI_02_ROI3_full_cluster1',
             'ICI_04_ROI1_bottom_cluster1',
             'ICI_04_ROI1_top_cluster1',
             'ICI_04_ROI1_top_cluster2',
             'ICI_05_ROI1_full_cluster1',
             'ICI_11_ROI1_bottom_cluster1',
             'ICI_11_ROI1_bottom_cluster2',
             'ICI_13_ROI1_full_cluster7',
             'ICI_13_ROI1_full_cluster9',
             'ICI_13_ROI2_full_cluster2',
             'ICI_13_ROI2_full_cluster8',
             'ICI_17_ROI2_full_cluster2',
             'ICI_19_ROI1_full_cluster1',
             'ICI_19_ROI2_full_cluster1',
             'ICI_19_ROI2_full_cluster2',
             'ICI_19_ROI2_full_cluster3',
             'ICI_21_ROI2_full_cluster1',
             'ICI_21_ROI2_full_cluster4',
             'ICI_23_ROI1_left_cluster5',
             'ICI_23_ROI1_right_cluster1', 
             'ICI_23_ROI1_right_cluster2',
             'ICI_23_ROI1_right_cluster3',
             'ICI_25_ROI2_full_cluster2',
             'ICI_25_ROI2_full_cluster6',
             'ICI_26_ROI1_left_cluster1',
             'ICI_27_ROI1_full_cluster1',
             'ICI_31_ROI1_bottom_cluster3',
             'ICI_31_ROI1_bottom_cluster4',
             'ICI_31_ROI1_top_cluster3',
             'ICI_37_ROI2_full_cluster1',
             'ICI_39_ROI1_full_cluster2',
             'ICI_39_ROI1_full_cluster3',
             'ICI_39_ROI1_full_cluster4',
             'ICI_39_ROI2_full_cluster4',
             'ICI_40_ROI3_full_cluster3',
             'ICI_41_ROI1_full_cluster2')

tls.meta <- subset(dbscan.out,tls.label %in% tls.label)

library(dplyr)

celltype_count <- tls.meta %>%
  group_by(tls.label, celltype) %>%
  summarise(count = n(), .groups = 'drop')

celltype_count <- celltype_count %>%
  group_by(tls.label) %>%
  mutate(percent = count / sum(count)*100)

wide_percent <- dcast(celltype_count, tls.label ~ celltype, value.var = "percent", fill = 0)
wide_percent <- wide_percent %>% as.data.frame() 
rownames(wide_percent) <- wide_percent$tls.label
wide_percent <- wide_percent[,-1]
colnames(wide_percent)
write.csv(wide_percent,'TLS Celltype add group.csv')


###### Identify T/B-TLS ######
wide_percent <- read.csv('TLS Celltype add group.csv')
colnames(wide_percent)
rownames(wide_percent) = wide_percent$X
wide_percent <- wide_percent[,-c(1,13)]
colnames(wide_percent)
wide_percent$tcell <- wide_percent$CD4T+wide_percent$CD8T
wide_percent <- wide_percent[,-c(3,4)]

a <- pheatmap(as.matrix(wide_percent), silent = TRUE)
library(pheatmap)
row_hc <- a$tree_row
plot(row_hc, main = "Row Clustering Dendrogram from pheatmap")
groups <- cutree(row_hc, k = 2)
wide_percent$group <- groups[rownames(wide_percent)]
table(wide_percent$group)
wide_percent$group <- ifelse(wide_percent$group %in% c(1),'B-TLS','T-TLS')
table(wide_percent$group)

colnames(wide_percent)
wide_percent$tb <- log((wide_percent$B+1)/(wide_percent$tcell+1))

a = ggplot(wide_percent,aes(x=group,y=tb,fill=group))+geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(width =0.3,size=1,color = "grey50", alpha = 0.6)+ 
  theme(panel.background = element_blank(),axis.line = element_line())+ 
  stat_compare_means(label="p.format")+    
  stat_boxplot(geom = "errorbar",width=0.5)+   
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#c095e4','#72b043'))
colnames(wide_percent)

red_line <- median(wide_percent$tb)
a = ggplot(wide_percent, aes_string(x = 'group', y = 'tb')) +
  geom_boxplot(outlier.shape = 21, width = 0.5, aes(fill = ifelse(..middle.. > red_line, "#C25E72", "#D9D9D9"))) +
  theme_minimal() +
  geom_jitter(width = 0.2, alpha = 0.5, size =1, color = "#919191") +
  theme(
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = '', y = 'Percentage') +
  stat_compare_means(aes(label = "p.format"), method = "wilcox.test", label = "p.format") +
  scale_fill_identity()

###### Prognosis ######

library(survminer)
library(survival)
patient.prognosis <- read.csv('Patient Clinical Info3.csv')
colnames(patient.prognosis)
fit <- survfit(Surv(OS, Censored) ~ Group, data = patient.prognosis)
a = ggsurvplot(fit,
               data = patient.prognosis,
               pval = TRUE,
               risk.table = TRUE,
               palette = "jco",
               title = "Prognosis of Immunotherapy",
               xlab = "Days",
               ylab = "Survival Probability")

pairwise_result <- pairwise_survdiff(
  Surv(OS, Censored) ~Group,
  data = patient.prognosis,
  p.adjust.method = "BH" 
)








