library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(UCell)
library(RColorBrewer)
library(pals)
library(ComplexHeatmap)

set.seed(12345)

# generate a new UMAp with the model saved to be able to project samples onto the reference

#ChenData_Healthy<-readRDS(file= "seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
#ChenData_Healthy<-UpdateSeuratObject(ChenData_Healthy)
#ChenData_Healthy<-RunUMAP(ChenData_Healthy, return.model = T, dims = 1:20)
#DimPlot(ChenData_Healthy, group.by = "Ctype")
#Save the updated object #
#saveRDS(ChenData_Healthy,file= "seurat_FiveHealthyD_updated.rds")

cluster_ord <- c('HSPC', 'LMPP', 'CLP', 'Pre-pro-B', 'Pro-B', 'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B', 'T','DC-Progenitor','Mono', 'GMP','NK', 'cDC', 'pDC', 'MEP')

clustcolors <-  c(
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

xlim <- c(min(ChenData_Healthy[['umap']]@cell.embeddings[,1]), max(ChenData_Healthy[['umap']]@cell.embeddings[,1]))
ylim <- c(min(ChenData_Healthy[['umap']]@cell.embeddings[,2]), max(ChenData_Healthy[['umap']]@cell.embeddings[,2]))

DimPlot(ChenData_Healthy, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/chen_healthyBM_umap.pdf')

###PB from iALL patients 
P636N <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/636N')
P636N <- CreateSeuratObject(counts = P636N, min.cells = 3, min.features = 200)
P636N[['percent.mt']] <- PercentageFeatureSet(P636N, pattern ='^MT-')
VlnPlot(P636N, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P636N <- subset(P636N, subset = nFeature_RNA > 200 & nFeature_RNA < 2600 & percent.mt < 10)
VlnPlot(P636N, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P636N <- NormalizeData(P636N) %>% FindVariableFeatures()
P636N <- ScaleData(P636N, features = rownames(P636N))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P636N,  reference.reduction = 'pca', dims = 1:20)
P636N <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P636N,
                  refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P636N) <- factor(P636N$predicted.id, levels = cluster_ord)
DimPlot(P636N, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P636N_projection.pdf')

FeaturePlot(P636N, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P636N_genes.pdf', width = 7, height = 5)


P1977N <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/1977N')
P1977N <- CreateSeuratObject(counts = P1977N, min.cells = 3, min.features = 200)
P1977N[['percent.mt']] <- PercentageFeatureSet(P1977N, pattern ='^MT-')
VlnPlot(P1977N, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1977N <- subset(P1977N, subset = nFeature_RNA > 200 & nFeature_RNA < 5300 & percent.mt < 10)
VlnPlot(P1977N, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1977N <- NormalizeData(P1977N) %>% FindVariableFeatures()
P1977N <- ScaleData(P1977N, features = rownames(P1977N))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P1977N,  reference.reduction = 'pca', dims = 1:20)
P1977N <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P1977N,
                   refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P1977N) <- factor(P1977N$predicted.id, levels = cluster_ord)
DimPlot(P1977N, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P1977N_projection.pdf')

FeaturePlot(P1977N, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P1977N_genes.pdf', width = 7, height = 5)


P6806R <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/6086R')
P6806R <- CreateSeuratObject(counts = P6806R, min.cells = 3, min.features = 200)
P6806R[['percent.mt']] <- PercentageFeatureSet(P6806R, pattern ='^MT-')
VlnPlot(P6806R, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P6806R <- subset(P6806R, subset = nFeature_RNA > 200 & nFeature_RNA < 4800 & percent.mt < 10 & nCount_RNA < 28000)
VlnPlot(P6806R, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P6806R <- NormalizeData(P6806R) %>% FindVariableFeatures()
P6806R <- ScaleData(P6806R, features = rownames(P6806R))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P6806R,  reference.reduction = 'pca', dims = 1:20)
P6806R <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P6806R,
                   refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P6806R) <- factor(P6806R$predicted.id, levels = cluster_ord)
DimPlot(P6806R, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P6806R_projection.pdf')

FeaturePlot(P6806R, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P6806R_genes.pdf', width = 7, height = 5)


P8010R <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/8010R')
P8010R <- CreateSeuratObject(counts = P8010R, min.cells = 3, min.features = 200)
P8010R[['percent.mt']] <- PercentageFeatureSet(P8010R, pattern ='^MT-')
VlnPlot(P8010R, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P8010R <- subset(P8010R, subset = nFeature_RNA > 200 & nFeature_RNA < 4600 & percent.mt < 10)
VlnPlot(P8010R, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P8010R <- NormalizeData(P8010R) %>% FindVariableFeatures()
P8010R <- ScaleData(P8010R, features = rownames(P8010R))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P8010R,  reference.reduction = 'pca', dims = 1:20)
P8010R <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P8010R,
                   refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P8010R) <- factor(P8010R$predicted.id, levels = cluster_ord)
DimPlot(P8010R, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P8010R_projection.pdf')

FeaturePlot(P8010R, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P8010R_genes.pdf', width = 7, height = 5)


P1776 <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/LX025_1776')
P1776 <- CreateSeuratObject(counts = P1776, min.cells = 3, min.features = 200)
P1776[['percent.mt']] <- PercentageFeatureSet(P1776, pattern ='^MT-')
VlnPlot(P1776, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1776 <- subset(P1776, subset = nFeature_RNA > 200 & nFeature_RNA < 5200 & percent.mt < 10)
VlnPlot(P1776, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1776 <- NormalizeData(P1776) %>% FindVariableFeatures()
P1776 <- ScaleData(P1776, features = rownames(P1776))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P1776,  reference.reduction = 'pca', dims = 1:20)
P1776 <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P1776,
                  refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P1776) <- factor(P1776$predicted.id, levels = cluster_ord)
DimPlot(P1776, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P1776_projection.pdf')

FeaturePlot(P1776, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P1776_genes.pdf', width = 7, height = 5)


P1175 <- Read10X(data.dir = 'human_data/Data/scRNAseq_Datarequest/LX024_1175')
P1175 <- CreateSeuratObject(counts = P1175, min.cells = 3, min.features = 200)
P1175[['percent.mt']] <- PercentageFeatureSet(P1175, pattern ='^MT-')
VlnPlot(P1175, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1175 <- subset(P1175, subset = nFeature_RNA > 200 & nFeature_RNA < 5300 & percent.mt < 10 & nCount_RNA < 28000)
VlnPlot(P1175, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
P1175 <- NormalizeData(P1175) %>% FindVariableFeatures()
P1175 <- ScaleData(P1175, features = rownames(P1175))

anchors <- FindTransferAnchors(reference = ChenData_Healthy, reference.assay = 'RNA',
                               query = P1175,  reference.reduction = 'pca', dims = 1:20)
P1175 <- MapQuery(anchorset = anchors, reference = ChenData_Healthy, query = P1175,
                  refdata = ChenData_Healthy$Ctype, reference.reduction = 'pca', reduction.model = 'umap')

Idents(P1175) <- factor(P1175$predicted.id, levels = cluster_ord)
DimPlot(P1175, cols = clustcolors) + coord_fixed() + xlim(xlim) + ylim(ylim)
ggsave('figs/march2025/P1175_projection.pdf')

FeaturePlot(P1175, features = c('CD24', 'CD19', 'HMGA2', 'IL7R'), pt.size = 1.5, order = T)
ggsave('figs/march2025/P1175_genes.pdf', width = 7, height = 5)


P636N$sample <- 'P636N'
P8010R$sample <- 'P8010R'
P1977N$sample <- 'P1977N'
P6806R$sample <- 'P6806R'
P1175$sample <- 'P1175'
P1776$sample <- 'P1776'

P636N$age <- 2.8
P8010R$age <- 0
P1977N$age <- 6.5
P6806R$age <- 11.3
P1175$age <- 6.6
P1776$age <- 0.7

# preleukemic signature
LeukBStem <- list(LeukBStem = c('ARID2','BAHCC1','BANK1','CAMK2D','CDK17','CSNK1G3','DIAPH3',
                                'DLG2','ELF1','FOXP1','GEM','HMCN1','HMGA2','HOXA10','HOXA5',
                                'HOXA7','HOXA9','HSPA4','IGF1','IL9R','KANSL1L','MACROD2',
                                'MAGI1','MAN1A1','MCTP1','MECOM','MEIS1','MTMR3','NKX2-3','NRIP1',
                                'PCDH7','PDE3B','PHTF2','PLCB4','PPP1R9B','PRKG1','RANBP9','RNF220',
                                'RSBN1L','SLIT2','SOX4','THSD4','TOX','UTRN',
                                'PPIA-','COX4I1-','NDUFB10-','RTRAF-','MACROH2A1-','MYC-','FAU-','ATP6V0C-'))

ChenData_Healthy <- AddModuleScore_UCell(ChenData_Healthy, features = LeukBStem, name = NULL)

P636N <- AddModuleScore_UCell(P636N, features = LeukBStem, name = NULL)
P1977N <- AddModuleScore_UCell(P1977N, features = LeukBStem, name = NULL)
P6806R <- AddModuleScore_UCell(P6806R, features = LeukBStem, name = NULL)
P8010R <- AddModuleScore_UCell(P8010R, features = LeukBStem, name = NULL)
P1175 <- AddModuleScore_UCell(P1175, features = LeukBStem, name = NULL)
P1776 <- AddModuleScore_UCell(P1776, features = LeukBStem, name = NULL)

ymax <- 0.35 #manually set ymax

my_palette <- brewer.pal(name = 'RdPu', n = 9)[3:9]
sc <- scale_color_gradientn(colors = my_palette, limits = c(0, ymax))

FeaturePlot(ChenData_Healthy, features = 'LeukBStem', pt.size = 1.5, order = T) + ggtitle('Healthy BM') + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/chen_healthyBM_sign_umap.pdf')

FeaturePlot(P636N, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P636N_sign_umap.pdf')

FeaturePlot(P1977N, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P1977N_sign_umap.pdf')

FeaturePlot(P8010R, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P8010R_sign_umap.pdf')

FeaturePlot(P6806R, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P6806R_sign_umap.pdf')

FeaturePlot(P1175, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P1175_sign_umap.pdf')

FeaturePlot(P1776, features = 'LeukBStem', pt.size = 1.5, order = T) + coord_fixed() + xlim(xlim) + ylim(ylim) & sc
ggsave('figs/march2025/P1776_sign_umap.pdf')

VlnPlot(P636N, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P636N_sign_vln.pdf')

VlnPlot(P1977N, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P1977N_sign_vln.pdf')

VlnPlot(P8010R, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P8010R_sign_vln.pdf')

VlnPlot(P6806R, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P6806R_sign_vln.pdf')

VlnPlot(P1175, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P1175_sign_vln.pdf')

VlnPlot(P1776, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend()
ggsave('figs/march2025/P1776_sign_vln.pdf')


((VlnPlot(P636N, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0)) | 
    (VlnPlot(P8010R, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0))) /
  ((VlnPlot(P6806R, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0)) | 
     (VlnPlot(P1977N, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0)))

VlnPlot(P1175, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0) | 
  VlnPlot(P1776, features = 'LeukBStem', cols = clustcolors, y.max = ymax) + NoLegend() + geom_abline(intercept = 0.03, slope = 0) 


merged <- merge(P636N, c(P8010R, P1977N, P6806R, P1175, P1776))
rm(P636N, P8010R, P1977N, P6806R, P1175, P1776)

VlnPlot(merged, 'nFeature_RNA', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/stam_nFeature_RNA_vln.pdf')

merged[['RNA']] <- JoinLayers(merged[['RNA']])
merged <- ScaleData(merged, assay = 'RNA', features = rownames(merged))

merged$sign_pos <- ifelse(merged$LeukBStem > 0.03, 'pos', 'neg')
sign <- table(merged$sample, merged$sign_pos)
write.csv(sign, 'figs/march2025/stam_sign_pos.csv')

merged$sample <- factor(merged$sample, levels = c('P8010R', 'P1776', 'P636N', 'P1977N', 'P1175', 'P6806R'))

pal2 <- brewer.pal(n = length(unique(merged$sample)), name = 'Set1')
names(pal2) <- unique(merged$sample)

VlnPlot(merged, features = 'LeukBStem', cols = pal2, y.max = ymax, group.by = 'sample') + NoLegend()
ggsave('figs/march2025/stam_sign_vln.pdf')

VlnPlot(merged, features = 'HMGA2', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/stam_HMGA2_vln.pdf')

VlnPlot(merged, features = 'CD24', group.by = 'sample') + NoLegend()
ggsave('figs/march2025/stam_CD24_vln.pdf')

merged_sub <- subset(merged, LeukBStem > 0.03)
merged_sub$predicted.id <- factor(merged_sub$predicted.id, levels = cluster_ord[cluster_ord %in% merged_sub$predicted.id])

predictions <- table(merged_sub$sample, merged_sub$predicted.id)
predictions <- predictions/rowSums(predictions)  # divide by number of cells for each patient
predictions <- data.frame(predictions)

table(merged$sample, merged$age)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% merged_sub$predicted.id]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/stam_predictions_bar_sign_pos.pdf', width = 7, height = 4.8)


merged_sub <- subset(merged, LeukBStem == 0)
merged_sub$predicted.id <- factor(merged_sub$predicted.id, levels = cluster_ord[cluster_ord %in% merged_sub$predicted.id])

predictions <- table(merged_sub$sample, merged_sub$predicted.id)
predictions <- predictions/rowSums(predictions)  # divide by number of cells for each patient
predictions <- data.frame(predictions)

table(merged$sample, merged$age)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% merged_sub$predicted.id]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/stam_predictions_bar_sign_neg.pdf', width = 7, height = 4.8)


VlnPlot(merged, 'CD19', group.by = 'sample') | (VlnPlot(merged, 'CD24', group.by = 'sample') + geom_abline(intercept = 0.5, slope = 0))

merged$CD24_status <- 'other'
merged$CD24_status[WhichCells(merged, expression = CD19 == 0 & CD24 > 0.5)] <- 'CD24pos_CD19lo'
merged$CD24_status[WhichCells(merged, expression = CD19 > 0.5 & CD24 > 0.5)] <- 'CD24pos_CD19pos'

VlnPlot(merged, 'HMGA2', group.by = 'sample') + geom_abline(intercept = 0.3, slope = 0)

merged$HMGA2_status <- 'HMGA2neg'
merged$HMGA2_status[WhichCells(merged, expression = HMGA2 > 0.3 )] <- 'HMGA2pos'

table(merged$CD24_status, merged$HMGA2_status)

merged$group <- paste0(merged$CD24_status, '_', merged$HMGA2_status)

hmga2 <- table(merged$sample, merged$group)
write.csv(hmga2, 'figs/march2025/stam_gene_groups.csv')

merged_sub <- subset(merged, CD24_status == 'CD24pos_CD19lo')
# merged_sub$predicted.id <- factor(merged_sub$predicted.id, levels = cluster_ord[cluster_ord %in% merged_sub$predicted.id])

predictions <- table(merged_sub$sample, merged_sub$predicted.id)
predictions <- predictions/rowSums(predictions)  # divide by number of cells for each patient
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% merged_sub$predicted.id]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/stam_predictions_bar_CD24pos_CD19lo.pdf', width = 7, height = 4.8)


predictions <- table(merged$sample, merged$predicted.id)
predictions <- predictions/rowSums(predictions)  # divide by number of cells for each patient
predictions <- data.frame(predictions)

predictions$Var2 <- factor(predictions$Var2, levels = rev(cluster_ord[cluster_ord %in% merged$predicted.id]))

ggplot(predictions, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_fill_manual(values = clustcolors, name = 'Predicted cell type') +
  xlab('Predicted cell type') +
  theme_bw() +
  ylab('Frequency')
ggsave('figs/march2025/stam_predictions_bar_allcells.pdf', width = 7, height = 4.8)


# heatmaps
merged_sub <- subset(merged, CD24_status == 'CD24pos_CD19lo')

merged_sub$predicted.id <- factor(merged_sub$predicted.id, levels = cluster_ord[cluster_ord %in% merged_sub$predicted.id])

annot_df <- data.frame('projection_celltype' = merged_sub$predicted.id, 'patient' = merged_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

goi <- c('CD24', 'HMGA2', 'MEIS1', 'PAX5', 'CD34', 'MME', 'CD19')

mtx <- GetAssayData(merged_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/stam_CD24pos_CD19neg_heatmap.pdf', width = 10, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()


merged_sub <- subset(merged, CD24_status == 'CD24pos_CD19pos')
merged_sub$predicted.id <- factor(merged_sub$predicted.id, levels = cluster_ord[cluster_ord %in% merged_sub$predicted.id])

annot_df <- data.frame('projection_celltype' = merged_sub$predicted.id, 'patient' = merged_sub$sample)
column_ha <- HeatmapAnnotation(df = annot_df, 
                               col = list(
                                 'projection_celltype' = clustcolors,
                                 'patient' = pal2)
)

mtx <- GetAssayData(merged_sub, assay = 'RNA', layer = 'scale.data')[goi, ]
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)

pdf('figs/march2025/stam_CD24pos_CD19pos_heatmap.pdf', width = 10, height = 4)
Heatmap(mtx, cluster_columns = T, show_row_names = T, show_column_names = F, top_annotation = column_ha, cluster_rows = F)
dev.off()
