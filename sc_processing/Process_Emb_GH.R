library(Seurat)
library(ShinyCell)
library(dplyr)

# Load from teh MTX files at GEO into cb.data
#cb.data <- Read10X_h5("/mnt/VR_Project/Data_For_Shamit/embryo/filtered_feature_bc_matrix_DONT_USE_ADTS_OR_HTOS.h5")

emb <- CreateSeuratObject(counts = cb.data[[1]], project = "CB", min.cells = 3, min.features = 200)

# HTOs from GEO
cells.needed <- read.delim("/mnt/VR_Project/Data_For_Shamit/Cell2Sample.embryo.R1_001.fastq.gz.tsv",header=T,row.names=1,sep="\t")

subs <- cells.needed[grep("PreproB",cells.needed$Most.likely.name),]

subs.nm <- paste(rownames(subs),"-1",sep="")

emb <- emb[,which(is.na(match(colnames(emb),subs.nm))==F)]


emb[["percent.mt"]] <- PercentageFeatureSet(emb, pattern =  "^mt-")

plot1 <- FeatureScatter(emb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(emb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#emb <- subset(emb, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)
emb <- subset(emb, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)


emb <- emb[-sort(c(grep("^Rps",rownames(emb)),grep("^Rpl",rownames(emb)))),]


emb <- NormalizeData(emb, normalization.method = "LogNormalize", scale.factor = 10000)
emb <- FindVariableFeatures(emb, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(emb)
emb <- ScaleData(emb, features = all.genes)

emb <- RunPCA(emb, features = VariableFeatures(object = emb))

emb <- FindNeighbors(emb, dims = 1:20)
emb <- FindClusters(emb, resolution = 0.5)

emb <- RunUMAP(emb, dims = 1:20)

pdf("For_CB.pdf",7,7)
DimPlot(emb, reduction = "umap",label=T)
dev.off()

FeaturePlot(emb,"Hmga2")

all.markers <- FindAllMarkers(emb)

write.table(all.markers,"Embryo_PreproB_Cluster_All_DiffStats.txt",row.names=T,col.names=NA,quote=F,sep="\t")

emb.markers %>%
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
DoHeatmap(emb, features = top20$gene) + NoLegend()
dev.off()


##
emb[["MLN"]] <- as.factor(subs$Most.likely.name[match(colnames(emb),subs.nm)])
DimPlot(emb, group.by=c("MLN"))

scConf = createConfig(emb)
makeShinyApp(emb, scConf, gene.mapping = TRUE,shiny.title = "Charlotta")


FeaturePlot(emb,"Hmga2")+DimPlot(emb, group.by=c("MLN","seurat_clusters"))

emb <- SetIdent(emb,value="MLN")

PreproB_Cre.markers <- FindMarkers(emb, ident.1 = "PreproB_Cre-", ident.2 = "PreproB_Cre+")

write.table(PreproB_Cre.markers,"Embryo_PreproB_Cre-_vs_Cre+_for_Preprob_Only.txt",row.names=T,col.names=NA,quote=F,sep="\t")

##Split the data

cre.neg <- emb[,which(emb[[]]$MLN=="PreproB_Cre-")]
cre.pos <- emb[,which(emb[[]]$MLN=="PreproB_Cre+")]

scConf = createConfig(cre.neg)
makeShinyApp(cre.neg, scConf, gene.mapping = TRUE,shiny.title = "Embryo_PreProB_CreNeg",shiny.dir = "shinyApp_Embryo_PreProB_CreNeg/")

scConf = createConfig(cre.neg)
makeShinyApp(cre.pos, scConf, gene.mapping = TRUE,shiny.title = "Embryo_PreProB_CrePos",shiny.dir = "shinyApp_Embryo_PreProB_CrePos/"")

