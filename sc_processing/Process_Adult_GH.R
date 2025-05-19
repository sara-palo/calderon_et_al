library(Seurat)
library(ShinyCell)
library(harmony)

# Load the MTX from GEO
#ad.spring <- Read10X_h5("/mnt/VR_Project/Data_For_Shamit/adult-spring/filtered_feature_bc_matrix_DONT_USE_ADTS_OR_HTOS.h5")[[1]]
#ad.summer <- Read10X_h5("/mnt/VR_Project/Data_For_Shamit/adult-summer/filtered_feature_bc_matrix_DONT_USE_ADTS_OR_HTOS.h5")[[1]]

spring <- CreateSeuratObject(counts = ad.spring, project = "CB1", min.cells = 3, min.features = 200)
summer <- CreateSeuratObject(counts = ad.summer, project = "CB2", min.cells = 3, min.features = 200)

# HTOs from GEO
spring.meta <- read.delim("/mnt/VR_Project/Data_For_Shamit/Cell2Sample.adult-spring.R1_001.fastq.gz.tsv",header=T,row.names=1,sep="\t")
summer.meta <- read.delim("/mnt/VR_Project/Data_For_Shamit/Cell2Sample.adult-summer.R1_001.fastq.gz.tsv",header=T,row.names=1,sep="\t")

rownames(spring.meta) <- paste(rownames(spring.meta),"-1",sep="")
rownames(summer.meta) <- paste(rownames(summer.meta),"-1",sep="")

spring[["MostLikelyName"]] <- spring.meta$Most.likely.name[match(colnames(spring),rownames(spring.meta))]
summer[["MostLikelyName"]] <- summer.meta$Most.likely.name[match(colnames(summer),rownames(summer.meta))]


adlt <- merge(spring, y = summer, add.cell.ids = c("spring", "summer"), project = "adult")

adlt[["percent.mt"]] <- PercentageFeatureSet(adlt, pattern = "^mt-")

plot1 <- FeatureScatter(adlt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adlt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#adlt <- subset(adlt, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA >1000 & percent.mt > 0.5)
adlt <- subset(adlt, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA >1000& nFeature_RNA)


adlt <- adlt[-sort(c(grep("^Rps",rownames(adlt)),grep("^Rpl",rownames(adlt)))),]


adlt <- NormalizeData(adlt, normalization.method = "LogNormalize", scale.factor = 10000)
adlt <- FindVariableFeatures(adlt, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(adlt)
adlt <- ScaleData(adlt, features = all.genes)

adlt <- RunPCA(adlt, features = VariableFeatures(object = adlt))

adlt <- RunHarmony(adlt, group.by.vars = "orig.ident")

#adlt <- FindNeighbors(adlt, dims = 1:30)
adlt <- FindNeighbors(adlt, reduction = "harmony",dims = 1:30)
adlt <- FindClusters(adlt, resolution = 0.5)

#adlt <- RunUMAP(adlt, dims = 1:30)
adlt <- RunUMAP(adlt, reduction = "harmony", dims = 1:30)
#adlt[["MLN"]] <- cells.needed$Most.likely.name[match(colnames(adlt),nms)]

DimPlot(adlt, reduction = "umap",group.by="orig.ident")



FeaturePlot(adlt,"Hmga2")


##

scConf = createConfig(adlt)
makeShinyApp(adlt, scConf, gene.mapping = TRUE,shiny.title = "Adult- with Harmony integration",shiny.dir = "shinyApp_Adult_withHarmony/")


FeaturePlot(adlt,"Hmga2")+DimPlot(adlt, group.by=c("MLN","seurat_cluster"))


#######
#Diffs

adlt <- SetIdent(adlt,value="MostLikelyName")

adlt.markers <- FindAllMarkers(adlt, only.pos = TRUE)

adlt.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.1) -> adult.diffs


write.table(adult.diffs,"Adult_differentials_log2FC_1_padj_0.1.txt",row.names=T,col.names=NA,quote=F,sep="\t")

adlt.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
DoHeatmap(adlt, features = top20$gene) + NoLegend()

pdf("Adult_Top20_markers__log2FC_1_padj_0.1.pdf",height=8,width=13)
DoHeatmap(adlt, features = top20$gene) + NoLegend()
dev.off()

####################################

PreproB_Cre.markers <- FindMarkers(adlt, ident.1 = "PreproB_Cre-", ident.2 = "PreproB_Cre+")
ProB_Cre.markers <- FindMarkers(adlt, ident.1 = "ProB_Cre-", ident.2 = "ProB_Cre+")
Cre.markers <- FindMarkers(adlt, ident.1 = "Cre-", ident.2 = "Cre+")

write.table(PreproB_Cre.markers,"Adult_PreproB_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")
write.table(ProB_Cre.markers,"Adult_ProB_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")
write.table(Cre.markers,"Adult_Cre_Cre-_vs_Cre+.txt",row.names=T,col.names=NA,quote=F,sep="\t")



VlnPlot(adlt, features =rownames(PreproB_Cre.markers)[1:10])
