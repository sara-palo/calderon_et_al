library(Seurat)
library(ShinyCell)
library(dplyr)

# Load from teh MTX files at GEO into cb.data
#cb.data <- Read10X_h5("/mnt/VR_Project/Data_For_Shamit/embryo/filtered_feature_bc_matrix_DONT_USE_ADTS_OR_HTOS.h5")

pbmc <- CreateSeuratObject(counts = cb.data[[1]], project = "CB", min.cells = 3, min.features = 200)

# HTOs from GEO
cells.needed <- read.delim("/mnt/VR_Project/Data_For_Shamit/Cell2Sample.embryo.R1_001.fastq.gz.tsv",header=T,row.names=1,sep="\t")

subs <- cells.needed[grep("PreproB",cells.needed$Most.likely.name),]

subs.nm <- paste(rownames(subs),"-1",sep="")

pbmc <- pbmc[,which(is.na(match(colnames(pbmc),subs.nm))==F)]


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern =  "^mt-")

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)


pbmc <- pbmc[-sort(c(grep("^Rps",rownames(pbmc)),grep("^Rpl",rownames(pbmc)))),]


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:20)

pdf("For_CB.pdf",7,7)
DimPlot(pbmc, reduction = "umap",label=T)
dev.off()

FeaturePlot(pbmc,"Hmga2")

all.markers <- FindAllMarkers(pbmc)

write.table(all.markers,"Embryo_PreproB_Cluster_All_DiffStats.txt",row.names=T,col.names=NA,quote=F,sep="\t")

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)



all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::filter(p_val_adj < 0.1) %>%
    ungroup() -> cluster.diffs

write.table(cluster.diffs,"Embryo_PreproB_Cluster_lfc_1_fdr_0.1.txt",row.names=T,col.names=NA,quote=F,sep="\t")

all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::filter(p_val_adj < 0.1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20

pdf("Embryo_PreproB_cluster_diffs_lfc_1_fdr_0.1_Heatmap.pdf",10,18)
DoHeatmap(pbmc, features = top20$gene) + NoLegend()
dev.off()


##
pbmc[["MLN"]] <- as.factor(subs$Most.likely.name[match(colnames(pbmc),subs.nm)])
DimPlot(pbmc, group.by=c("MLN"))

scConf = createConfig(pbmc)
makeShinyApp(pbmc, scConf, gene.mapping = TRUE,shiny.title = "Charlotta")


FeaturePlot(pbmc,"Hmga2")+DimPlot(pbmc, group.by=c("MLN","seurat_clusters"))

pbmc <- SetIdent(pbmc,value="MLN")

PreproB_Cre.markers <- FindMarkers(pbmc, ident.1 = "PreproB_Cre-", ident.2 = "PreproB_Cre+")

write.table(PreproB_Cre.markers,"Embryo_PreproB_Cre-_vs_Cre+_for_Preprob_Only.txt",row.names=T,col.names=NA,quote=F,sep="\t")

##Split the data

cre.neg <- pbmc[,which(pbmc[[]]$MLN=="PreproB_Cre-")]
cre.pos <- pbmc[,which(pbmc[[]]$MLN=="PreproB_Cre+")]

scConf = createConfig(cre.neg)
makeShinyApp(cre.neg, scConf, gene.mapping = TRUE,shiny.title = "Embryo_PreProB_CreNeg",shiny.dir = "shinyApp_Embryo_PreProB_CreNeg/")

scConf = createConfig(cre.neg)
makeShinyApp(cre.pos, scConf, gene.mapping = TRUE,shiny.title = "Embryo_PreProB_CrePos",shiny.dir = "shinyApp_Embryo_PreProB_CrePos/"")

