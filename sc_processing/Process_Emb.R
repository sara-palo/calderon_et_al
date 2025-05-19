library(Seurat)
library(ShinyCell)
library(dplyr)
library(assertthat)


cb.data <- Read10X_h5("/mnt/VR_Project/Data_For_Shamit/embryo/filtered_feature_bc_matrix_DONT_USE_ADTS_OR_HTOS.h5")

emb <- CreateSeuratObject(counts = cb.data[[1]], project = "CB", min.cells = 3, min.features = 200)



cells.needed <- read.delim("/mnt/VR_Project/Data_For_Shamit/Cell2Sample.embryo.R1_001.fastq.gz.tsv",header=T,row.names=1,sep="\t")
nms <- paste(rownames(cells.needed),"-1",sep="")

#Keep all cells
###############


#subs <- cells.needed[grep("PreproB",cells.needed$Most.likely.name),]

#subs.nm <- paste(rownames(subs),"-1",sep="")

#emb <- emb[,which(is.na(match(colnames(emb),subs.nm))==F)]


emb[["percent.mt"]] <- PercentageFeatureSet(emb, pattern = "^mt-")

plot1 <- FeatureScatter(emb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(emb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#emb <- subset(emb, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)
emb <- subset(emb, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)


emb <- emb[-sort(c(grep("^Rps",rownames(emb)),grep("^Rpl",rownames(emb)))),]




emb <- NormalizeData(emb, normalization.method = "LogNormalize", scale.factor = 10000)
emb <- FindVariableFeatures(emb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(emb)
emb <- ScaleData(emb, features = all.genes)

emb <- RunPCA(emb, features = VariableFeatures(object = emb))

emb <- FindNeighbors(emb, dims = 1:30)
emb <- FindClusters(emb, resolution = 0.5)

emb <- RunUMAP(emb, dims = 1:30)

emb[["MLN"]] <- cells.needed$Most.likely.name[match(colnames(emb),nms)]

DimPlot(emb, reduction = "umap")


FeaturePlot(emb,"Meis1")


##

scConf = createConfig(emb)
makeShinyApp(emb, scConf, gene.mapping = TRUE,shiny.title = "Charlotta")


FeaturePlot(emb,"Hmga2")+DimPlot(emb, group.by=c("MLN","seurat_cluster"))


##
emb <- SetIdent(emb,value="MLN")

emb.markers <- FindAllMarkers(emb, only.pos = TRUE)

emb.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.1) -> embryo.diffs


write.table(embryo.diffs,"Embryo_differentials_log2FC_1_padj_0.1.txt",row.names=T,col.names=NA,quote=F,sep="\t")

emb.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
DoHeatmap(emb, features = top20$gene) + NoLegend()

pdf("Embryo_Top20_markers__log2FC_1_padj_0.1.pdf",height=8,width=13)
DoHeatmap(emb, features = top20$gene) + NoLegend()
dev.off()


##################
## Diffs


emb <- SetIdent(emb,value="MLN")

PreproB_Cre.markers <- FindMarkers(emb, ident.1 = "PreproB_Cre-", ident.2 = "PreproB_Cre+")
ProB_Cre.markers <- FindMarkers(emb, ident.1 = "ProB_Cre-", ident.2 = "ProB_Cre+")
Cre.markers <- FindMarkers(emb, ident.1 = "Cre-", ident.2 = "Cre+")

write.table(PreproB_Cre.markers,"Embryo_PreproB_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")
write.table(ProB_Cre.markers,"Embryo_ProB_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")
write.table(Cre.markers,"Embryo_Cre_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")


################
meta.all <- readRDS("/home/shamit/Boiers/Test/Scripts/shinyApp/sc1meta.rds")
meta <- readRDS("/home/shamit/Boiers/Test/Scripts/shinyApp/sc1meta.rds")[,14:15]


####
emb1 <- emb

metat <- as.matrix(meta)
rownames(metat) <- colnames(emb)
colnames(metat) <- c("umap_1","umap_2")


emb1@reductions$umap@cell.embeddings <- metat
DimPlot(emb1)


emb1[[]]$seurat_clusters <- as.factor(readRDS("/home/shamit/Boiers/Test/Scripts/shinyApp/sc1meta.rds")$seurat_clusters)
DimPlot(emb1,group.by="seurat_clusters")

emb1 <- SetIdent(emb1,value="seurat_clusters")
DimPlot(emb1,group.by="MLN")
saveRDS(emb1,"emb1.rds")

cl5 <- readLines("Cluster_5_Charlotta.txt")

MLN.wc5 <- emb1$MLN
MLN.wc5[match(cl5,colnames(emb1))] <- "Clus5"

emb1$MLNc5 <- as.factor(MLN.wc5)


emb1 <- SetIdent(emb1,value="MLNc5")
cluster5.markers <- FindMarkers(emb1, ident.1 = "Clus5", ident.2 = c("PreproB_Cre-", "PreproB_Cre+"))

write.table(cluster5.markers,"cluster5_vs_AllProproB_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")

cluster5.markers.ppbwto <- FindMarkers(emb1, ident.1 = "Clus5", ident.2 = c("PreproB_Cre-"))
write.table(cluster5.markers.ppbwto,"cluster5_vs_PreproB_Cre-_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")


cluster5.markers.all <- FindAllMarkers(emb1)


 cluster5.markers.all%>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup() -> top30

write.table(top30,"cluster5_vs_allgroups.txt",row.names=T,col.names=NA,quote=F,sep="\t")


scConf = createConfig(emb1)
makeShinyApp(emb1, scConf, gene.mapping = TRUE,shiny.title = "All_data_with_cluster5",shiny.dir = "shinyApp_All_data_with_Cluster5/")

makeShinyApp(emb1, scConf, gene.mapping = TRUE,gex.assay="RNA",shiny.title = "All_data_with_cluster5",shiny.dir = "shinyApp_All_data_with_Cluster5/")

emb1j <- JoinLayers(emb1)

makeShinyApp(emb1j, scConf, gene.mapping = TRUE,gex.assay="RNA",shiny.title = "All_data_with_cluster5",shiny.dir = "shinyApp_All_data_with_Cluster5/")


pdf("Cluster5.pdf",7,7)
DimPlot(emb1,label=T,cells.highlight=cl5)
dev.off()

pdf("Cluster5_nolabs.pdf",7,7)
DimPlot(emb1,label=F,cells.highlight=cl5)
dev.off()


### Cre- Clust 2 vs Cluster 5

#emb1

clus.mln <- paste(emb1[[]]$seurat_clusters,emb1[[]]$MLN,sep="_")
clus.mln[which(emb1[[]]$MLNc5=="Clus5")] <- "Clus5"
emb1$clus.mln.c5 <- as.factor(clus.mln)

emb1 <- SetIdent(emb1,value="clus.mln.c5")
cl2.cren.vs.c5 <- FindMarkers(emb1, ident.1 = "Clus5", ident.2 = c("2_Cre-"))
write.table(cl2.cren.vs.c5,"cluster5_vs_Clus2_Cre_neg_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")



#### cluster 2 from emb1 cre- cre +
clus.mln <- paste(emb1[[]]$seurat_clusters,emb1[[]]$MLN,sep="_")
emb1$clus.mln <- as.factor(clus.mln)
emb1 <- SetIdent(emb1,value="clus.mln")
cl2.cren.vscrep <- FindMarkers(emb1, ident.1 = "2_Cre-", ident.2 = c("2_Cre+"))

write.table(cl2.cren.vscrep,"Clus2_Cre_neg_vs_Cre_pos_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")

#### cluster 4 from emb1 cre- cre +

cl4.cren.vscrep <- FindMarkers(emb1, ident.1 = "4_Cre-", ident.2 = c("4_Cre+"))
write.table(cl4.cren.vscrep,"Clus4_Cre_neg_vs_Cre_pos_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")

#### cluster 1 from emb1 cre- cre +

cl1.cren.vscrep <- FindMarkers(emb1, ident.1 = "1_Cre-", ident.2 = c("1_Cre+"))
write.table(cl1.cren.vscrep,"Clus1_Cre_neg_vs_Cre_pos_Markers.txt",row.names=T,col.names=NA,quote=F,sep="\t")


#######################

