library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(UCell)
library(RColorBrewer)
library(pals)
library(ComplexHeatmap)
library(stringr)

set.seed(12345)

cluster_ord <- c('HSPC', 'LMPP', 'CLP', 'Pre-pro-B', 'Pro-B', 'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B', 'T','DC-Progenitor','Mono', 'GMP','NK', 'cDC', 'pDC', 'MEP')

clustcolors <- c(
  'HSPC' = '#b33e52',
  'LMPP' = '#1b7f7f',
  'CLP' = '#dcabc5',
  'Pre-pro-B' = '#db79b1',
  'Pro-B' = '#9471b4',
  'Pre-B' = '#9a45af',
  'Immature-B' = '#6b3e97',
  'Mature-B' = '#c9b2d5',
  'Plasma-B' = '#9590ff',
  'T' = '#7c223e',
  'DC-Progenitor' = '#a6cee3',
  'Mono' = '#2eb0ce',
  'GMP' = '#1f78b3',
  'NK' = '#f7c68b',
  'cDC' = '#cc8d42',
  'pDC' = '#c66b09',
  'MEP' = '#9dbfab'
)

#Read the updated file# 
ChenData_Healthy <- readRDS('human_data/Data/chen/seurat_FiveHealthyD_updated.rds')

ChenData_Healthy$Ctype <- factor(ChenData_Healthy$Ctype, levels = cluster_ord)
Idents(ChenData_Healthy) <- 'Ctype'

ChenData <- readRDS('human_data/Data/chen/seurat_scRNA_patient_clean.rds')
ChenData <- UpdateSeuratObject(ChenData)

table(ChenData@meta.data$Ctype)
table(ChenData@meta.data$sample)
table(ChenData@meta.data$Ctype_Final)

ChenData_nonblast <- subset(ChenData, Ctype %in% c('Blasts', 'Progenitors'), invert = T)
ChenData <- subset(ChenData, Ctype %in% c('Blasts', 'Progenitors'))

LeukBStem <- list(LeukBStem = c('ARID2','BAHCC1','BANK1','CAMK2D','CDK17','CSNK1G3','DIAPH3',
                                'DLG2','ELF1','FOXP1','GEM','HMCN1','HMGA2','HOXA10','HOXA5',
                                'HOXA7','HOXA9','HSPA4','IGF1','IL9R','KANSL1L','MACROD2',
                                'MAGI1','MAN1A1','MCTP1','MECOM','MEIS1','MTMR3','NKX2-3','NRIP1',
                                'PCDH7','PDE3B','PHTF2','PLCB4','PPP1R9B','PRKG1','RANBP9','RNF220',
                                'RSBN1L','SLIT2','SOX4','THSD4','TOX','UTRN',
                                'PPIA-','COX4I1-','NDUFB10-','RTRAF-','MACROH2A1-','MYC-','FAU-','ATP6V0C-'))

ChenData <- AddModuleScore_UCell(ChenData, features = LeukBStem, name = NULL)
ChenData_nonblast <- AddModuleScore_UCell(ChenData_nonblast, features = LeukBStem, name = NULL)

table(ChenData_nonblast$Ctype)

subs <- subset(ChenData, PresumedFusion %in% c('MLL-AF4', 'MLL-AF9', 'MLL-ENL'))

tb <- table(subs$sample, subs$PresumedFusion)
fusions <- colnames(tb)[apply(tb, 1, function(x) which(x > 0))]
names(fusions) <- rownames(tb)

subs$sample <- factor(subs$sample, levels = c(names(which(fusions == 'MLL-AF4')), names(which(fusions == 'MLL-AF9')), names(which(fusions == 'MLL-ENL'))))

VlnPlot(subs, 'LeukBStem', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/chen_all3groups_sign.pdf')

VlnPlot(ChenData_nonblast, 'LeukBStem', group.by = 'Ctype_Final') + NoLegend()
ggsave('figs/march2025/chen_nonblast_sign.pdf')

VlnPlot(subs, 'HMGA2', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/chen_all3groups_HMGA2.pdf')

VlnPlot(subs, 'CD24', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/chen_all3groups_CD24.pdf')

rm(ChenData_nonblast)
rm(subs)

ChenData <- ScaleData(ChenData, features = rownames(ChenData))

table(ChenData$PresumedFusion)

VlnPlot(ChenData, 'HMGA2', group.by = 'projCtype') + geom_abline(intercept = 0.3, slope = 0) + NoLegend()

ChenData$CD24_status <- 'other'
ChenData$CD24_status[WhichCells(ChenData, expression = CD19 == 0 & CD24 > 0.5)] <- 'CD24pos_CD19lo'
ChenData$CD24_status[WhichCells(ChenData, expression = CD19 > 0.5 & CD24 > 0.5)] <- 'CD24pos_CD19pos'

ChenData$HMGA2_status <- 'HMGA2neg'
ChenData$HMGA2_status[WhichCells(ChenData, expression = HMGA2 > 0.3 )] <- 'HMGA2pos'

table(ChenData$CD24_status, ChenData$HMGA2_status)

ChenData$group <- paste0(ChenData$CD24_status, '_', ChenData$HMGA2_status)

MA4 <- subset(ChenData, PresumedFusion == 'MLL-AF4')
MA9 <- subset(ChenData, PresumedFusion == 'MLL-AF9')
MENL <- subset(ChenData, PresumedFusion == 'MLL-ENL')

rm(ChenData)

MA4$sign_pos <- ifelse(MA4$LeukBStem > 0.03, 'pos', 'neg')
sign <- table(MA4$sample, MA4$sign_pos)
write.csv(sign, 'figs/march2025/chen_MA4_sign_pos.csv')

MA9$sign_pos <- ifelse(MA9$LeukBStem > 0.03, 'pos', 'neg')
sign <- table(MA9$sample, MA9$sign_pos)
write.csv(sign, 'figs/march2025/chen_MA9_sign_pos.csv')

MENL$sign_pos <- ifelse(MENL$LeukBStem > 0.03, 'pos', 'neg')
sign <- table(MENL$sample, MENL$sign_pos)
write.csv(sign, 'figs/march2025/chen_MENL_sign_pos.csv')


MA4$projCtype <- factor(str_remove(MA4$projCtype, pattern = '-like'), levels = cluster_ord)
MA9$projCtype <- factor(str_remove(MA9$projCtype, pattern = '-like'), levels = cluster_ord)
MENL$projCtype <- factor(str_remove(MENL$projCtype, pattern = '-like'), levels = cluster_ord)

Idents(MA4) <- 'projCtype'
Idents(MA9) <- 'projCtype'
Idents(MENL) <- 'projCtype'

((VlnPlot(MA4, features = 'CD24', group.by = 'projCtype', cols = clustcolors) + NoLegend() + ggtitle('CD24 - MLL-AF4') + geom_abline(intercept = 0.5, slope = 0)) /
    (VlnPlot(MA9, features = 'CD24', group.by = 'projCtype', cols = clustcolors) + NoLegend()+ ggtitle('CD24 - MLL-AF9') + geom_abline(intercept = 0.5, slope = 0)) /
    (VlnPlot(MENL, features = 'CD24', group.by = 'projCtype', cols = clustcolors) + NoLegend()+ ggtitle('CD24 - MLL-ENL') + geom_abline(intercept = 0.5, slope = 0))) |
  ((VlnPlot(MA4, features = 'CD19', group.by = 'projCtype', cols = clustcolors) + NoLegend() + ggtitle('CD19 - MLL-AF4')) /
     (VlnPlot(MA9, features = 'CD19', group.by = 'projCtype', cols = clustcolors) + NoLegend()+ ggtitle('CD19 - MLL-AF9')) /
     (VlnPlot(MENL, features = 'CD19', group.by = 'projCtype', cols = clustcolors) + NoLegend()+ ggtitle('CD19 - MLL-ENL')))

ymax <- max(c(MA4$LeukBStem, MA9$LeukBStem, MENL$LeukBStem))

VlnPlot(MA4, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/chen_MA4_sign_vln.pdf')

VlnPlot(MA9, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/chen_MA9_sign_vln.pdf')

VlnPlot(MENL, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/chen_MENL_sign_vln.pdf')

MA4_sub <- subset(MA4, LeukBStem > 0.03)
MA4_sub$projCtype <- factor(MA4_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA4_sub$projCtype])

predictions <- table(MA4_sub$sample, MA4_sub$projCtype)
predictions <- predictions/rowSums(predictions)  # divide by number of cells for each patient
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA4_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA4_predictions_bar_sign_pos.pdf', width = 5.5, height = 4.5)


MA4_sub <- subset(MA4, CD24_status == 'CD24pos_CD19lo')
MA4_sub$projCtype <- factor(MA4_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA4_sub$projCtype])

predictions <- table(MA4_sub$sample, MA4_sub$projCtype)
predictions <- predictions/rowSums(predictions)
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA4_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA4_predictions_bar_CD24pos_CD19lo.pdf', width = 5.5, height = 4.5)


predictions <- table(MA4$sample, MA4$projCtype)
predictions <- predictions/rowSums(predictions)  
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA4$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA4_predictions_bar_allcells.pdf', width = 5.5, height = 4.5)


MA9_sub <- subset(MA9, LeukBStem > 0.03)
MA9_sub$projCtype <- factor(MA9_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA9_sub$projCtype])

predictions <- table(MA9_sub$sample, MA9_sub$projCtype)
predictions <- predictions/rowSums(predictions)  
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA9_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA9_predictions_bar_sign_pos.pdf', width = 5.5, height = 4.5)

MA9_sub <- subset(MA9, CD24_status == 'CD24pos_CD19lo')
MA9_sub$projCtype <- factor(MA9_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA9_sub$projCtype])

predictions <- table(MA9_sub$sample, MA9_sub$projCtype)
predictions <- predictions/rowSums(predictions)
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA9_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA9_predictions_bar_CD24pos_CD19lo.pdf', width = 5.5, height = 4.5)


predictions <- table(MA9$sample, MA9$projCtype)
predictions <- predictions/rowSums(predictions)  
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MA9$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MA9_predictions_bar_allcells.pdf', width = 5.5, height = 4.5)


MENL_sub <- subset(MENL, LeukBStem > 0.03)
MENL_sub$projCtype <- factor(MENL_sub$projCtype, levels = cluster_ord[cluster_ord %in% MENL_sub$projCtype])

predictions <- table(MENL_sub$sample, MENL_sub$projCtype)
predictions <- predictions/rowSums(predictions)  
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MENL_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MENL_predictions_bar_sign_pos.pdf', width = 5.5, height = 4.5)

MENL_sub <- subset(MENL, CD24_status == 'CD24pos_CD19lo')
MENL_sub$projCtype <- factor(MENL_sub$projCtype, levels = cluster_ord[cluster_ord %in% MENL_sub$projCtype])

predictions <- table(MENL_sub$sample, MENL_sub$projCtype)
predictions <- predictions/rowSums(predictions)
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MENL_sub$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MENL_predictions_bar_CD24pos_CD19lo.pdf', width = 5.5, height = 4.5)

predictions <- table(MENL$sample, MENL$projCtype)
predictions <- predictions/rowSums(predictions)  
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% MENL$projCtype]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/chen_MENL_predictions_bar_allcells.pdf', width = 5.5, height = 4.5)

hmga2 <- table(MA4$sample, MA4$group)
write.csv(hmga2, 'figs/march2025/chen_MA4_gene_groups.csv')

hmga2 <- table(MA9$sample, MA9$group)
write.csv(hmga2, 'figs/march2025/chen_MA9_gene_groups.csv')

hmga2 <- table(MENL$sample, MENL$group)
write.csv(hmga2, 'figs/march2025/chen_MENL_gene_groups.csv')


# heatmaps
goi <- c('CD24', 'HMGA2', 'MEIS1', 'PAX5', 'CD34', 'MME', 'CD19')

pal2 <- brewer.pal(n = length(unique(MA4$sample)), name = 'Set1')
names(pal2) <- unique(MA4$sample)

MA4_sub <- subset(MA4, CD24_status == 'CD24pos_CD19lo')
MA4_sub$projCtype <- factor(MA4_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA4_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MA4_sub$projCtype, 'patient' = MA4_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MA4_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MA4_CD24pos_CD19neg_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()

MA4_sub <- subset(MA4, CD24_status == 'CD24pos_CD19pos')
MA4_sub$projCtype <- factor(MA4_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA4_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MA4_sub$projCtype, 'patient' = MA4_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MA4_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MA4_CD24pos_CD19pos_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


pal2 <- brewer.pal(n = length(unique(MA9$sample)), name = 'Set1')
names(pal2) <- unique(MA9$sample)

MA9_sub <- subset(MA9, CD24_status == 'CD24pos_CD19lo')
MA9_sub$projCtype <- factor(MA9_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA9_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MA9_sub$projCtype, 'patient' = MA9_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)


mtx <- GetAssayData(MA9_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MA9_CD24pos_CD19neg_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()

MA9_sub <- subset(MA9, CD24_status == 'CD24pos_CD19pos')
MA9_sub$projCtype <- factor(MA9_sub$projCtype, levels = cluster_ord[cluster_ord %in% MA9_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MA9_sub$projCtype, 'patient' = MA9_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MA9_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MA9_CD24pos_CD19pos_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


pal2 <- brewer.pal(n = length(unique(MENL$sample)), name = 'Set1')
names(pal2) <- unique(MENL$sample)

MENL_sub <- subset(MENL, CD24_status == 'CD24pos_CD19lo')
MENL_sub$projCtype <- factor(MENL_sub$projCtype, levels = cluster_ord[cluster_ord %in% MENL_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MENL_sub$projCtype, 'patient' = MENL_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MENL_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MENL_CD24pos_CD19neg_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()

MENL_sub <- subset(MENL, CD24_status == 'CD24pos_CD19pos')
MENL_sub$projCtype <- factor(MENL_sub$projCtype, levels = cluster_ord[cluster_ord %in% MENL_sub$projCtype])

annot_df <- data.frame('projection_celltype' = MENL_sub$projCtype, 'patient' = MENL_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MENL_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/chen_MENL_CD24pos_CD19pos_heatmap.pdf', width = 11, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


# healthy BM
ChenData_Healthy$CD24_status <- 'other'
ChenData_Healthy$CD24_status[WhichCells(ChenData_Healthy, expression = CD19 == 0 & CD24 > 0.5)] <- 'CD24pos_CD19lo'
ChenData_Healthy$CD24_status[WhichCells(ChenData_Healthy, expression = CD19 > 0.5 & CD24 > 0.5)] <- 'CD24pos_CD19pos'

ChenData_Healthy$HMGA2_status <- 'HMGA2neg'
ChenData_Healthy$HMGA2_status[WhichCells(ChenData_Healthy, expression = HMGA2 > 0.3 )] <- 'HMGA2pos'

table(ChenData_Healthy$CD24_status, ChenData_Healthy$HMGA2_status)
ChenData_Healthy$group <- paste0(ChenData_Healthy$CD24_status, '_', ChenData_Healthy$HMGA2_status)

hmga2 <- table(ChenData_Healthy$sample, ChenData_Healthy$group)
write.csv(hmga2, 'figs/march2025/chen_healthy_gene_groups.csv')

ChenData_Healthy$Ctype <- factor(ChenData_Healthy$Ctype, levels = rev(cluster_ord))

ChenData_Healthy <- ScaleData(ChenData_Healthy, features = rownames(ChenData_Healthy))

goi2 <- c('CD24', 'HMGA2', 'PAX5', 'CD34', 'IL7R', 'MME', 'CD19')

DotPlot(ChenData_Healthy, features = goi2, group.by = 'Ctype') &
  scale_colour_gradientn(colours = (brewer.pal(name = 'PuRd', n = 9)[2:9]))
ggsave('figs/march2025/chen_healthyBM_dotplot.pdf', width = 8, height = 4.5)


VlnPlot(ChenData_Healthy, 'HMGA2')

pal2 <- brewer.pal(n = length(unique(MA4$sample)), name = 'Set1')
names(pal2) <- unique(MA4$sample)

MA4$projCtype <- factor(MA4$projCtype, levels = cluster_ord[cluster_ord %in% MA4$projCtype])

annot_df <- data.frame('projection_celltype' = MA4$projCtype, 'patient' = MA4$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MA4, assay = 'RNA', layer = 'scale.data')[goi, ]

pdf('figs/march2025/chen_MA4_all_heatmap.pdf', width = 14, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


pal2 <- brewer.pal(n = length(unique(MA9$sample)), name = 'Set1')
names(pal2) <- unique(MA9$sample)

MA9$projCtype <- factor(MA9$projCtype, levels = cluster_ord[cluster_ord %in% MA9$projCtype])

annot_df <- data.frame('projection_celltype' = MA9$projCtype, 'patient' = MA9$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MA9, assay = 'RNA', layer = 'scale.data')[goi, ]

pdf('figs/march2025/chen_MA9_all_heatmap.pdf', width = 14, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


pal2 <- brewer.pal(n = length(unique(MENL$sample)), name = 'Set1')
names(pal2) <- unique(MENL$sample)

MENL$projCtype <- factor(MENL$projCtype, levels = cluster_ord[cluster_ord %in% MENL$projCtype])

annot_df <- data.frame('projection_celltype' = MENL$projCtype, 'patient' = MENL$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(MENL, assay = 'RNA', layer = 'scale.data')[goi, ]

pdf('figs/march2025/chen_MENL_all_heatmap.pdf', width = 14, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()
