
####packages####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("harmony")
BiocManager::install("UCell")
install.packages("devtools")
BiocManager::install ("limma")
install.packages('Seurat') #V5
options(SeuratData.repo.use = "http://seurat.nygenome.org")

####libraries####

library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)
library(dplyr)
#for signature
library(UCell)
library(RColorBrewer)
library("viridis") 

#for reading and managing the SCPCL data files
library("AnnotationDbi")
library("org.Hs.eg.db")
library(readxl)

####Pediatric Healthy Bone Marrow-Creating a REF-Data Series GSE132509 ####
PBM1<-Read10X(data.dir="data_dir/GSE132509_RAW/PBMMC_1")
PBM1<- CreateSeuratObject(counts= PBM1, project = "PBM1", min.cells = 3, min.features = 500)
PBM1[["percent.mt"]]<-PercentageFeatureSet(PBM1, pattern ="^MT-")
PBM1<- subset(PBM1, subset = nFeature_RNA>500 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(PBM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PBM3<-Read10X(data.dir="data_dir/GSE132509_RAW/PBMMC_3")
PBM3<- CreateSeuratObject(counts= PBM3, project = "PBM3", min.cells = 3, min.features = 500)
PBM3[["percent.mt"]]<-PercentageFeatureSet(PBM3, pattern ="^MT-")
PBM3<- subset(PBM3, subset = nFeature_RNA>500 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(PBM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Merging the BM data
PBM.Healthy<- merge(x=PBM1, y= PBM3, add.cell.ids = c("PBM1", "PBM3"))
PBM.Healthy<- NormalizeData(PBM.Healthy)
PBM.Healthy<-FindVariableFeatures(PBM.Healthy)
PBM.Healthy<-ScaleData(PBM.Healthy, features = rownames(PBM.Healthy))
PBM.Healthy<- RunPCA(PBM.Healthy, features = VariableFeatures(object = PBM.Healthy))
ElbowPlot(PBM.Healthy)
PBM.Healthy <- FindNeighbors(PBM.Healthy,reduction = "pca", dims = 1:15)
PBM.Healthy <- FindClusters(PBM.Healthy,resolution = 0.9)
PBM.Healthy <- RunUMAP(PBM.Healthy, dims = 1:15,  reduction = "pca")
DimPlot(PBM.Healthy, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = T)

##Integrating healthy BM data with harmony and creating a REF, not relevant if mapping to a ref
PBM.Healthy <- RunHarmony(PBM.Healthy,"orig.ident")
PBM.Healthy[["RNA"]] <- JoinLayers(PBM.Healthy[["RNA"]])
PBM.Healthy <- FindNeighbors(PBM.Healthy,reduction = "harmony", dims = 1:15)
ElbowPlot(PBM.Healthy)
PBM.Healthy <- FindClusters(PBM.Healthy,resolution = 0.9)
PBM.Healthy <- RunUMAP(PBM.Healthy, dims = 1:15,  reduction = "harmony", return.model = TRUE)
DimPlot(PBM.Healthy, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = T)
PBM.Healthy.embed<- Embeddings(PBM.Healthy,"harmony")
SaveSeuratRds(PBM.Healthy,file= "data_dir/PediatricBM.Healthy")
PBM.Healthy<-readRDS(file= "data_dir/Pediatric_ZM_BM.Healthy")
#For cell count per cluster
table(PBM.Healthy@meta.data$seurat_clusters)

###MAKING the REF UMAP similar to the UMAP in the original article
PBM.Healthy <- FindNeighbors(PBM.Healthy,reduction = "harmony", dims = 1:12)
ElbowPlot(PBM.Healthy)
PBM.Healthy <- FindClusters(PBM.Healthy,resolution = 0.3)
PBM.Healthy <- RunUMAP(PBM.Healthy, dims = 1:12, reduction = "harmony", return.model = TRUE, seed.use=123456789)
DimPlot(PBM.Healthy, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = T)

#PBM.Healthy.embed<- Embeddings(PBM.Healthy,"harmony")
Refmap<- DimPlot(PBM.Healthy, reduction = "umap", group.by = "seurat_clusters",
                 pt.size = 1, label = T)  + NoLegend() + ggtitle("Healthy pediatric BM-Reference annotations")

#Gene expression
FeaturePlot(PBM.Healthy, features = c("HMGA2", "CD34", "CD24", "CD19"),pt.size = 1,order = T)&
  scale_colour_gradientn(colours = (brewer.pal(name="PuBu",n=9)[2:9]))


####Adding 4 patient samples from PMID: 34304246####
#PB from P636N
P636N<-Read10X(data.dir="data_dir/scRNAseq_Datarequest/636N")
P636N<- CreateSeuratObject(counts= P636N, project = "P636N", min.cells = 3, min.features = 200)
P636N[["percent.mt"]]<-PercentageFeatureSet(P636N, pattern ="^MT-")
P636N<- subset(P636N, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(P636N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P636N <- NormalizeData(P636N)
P636N <- FindVariableFeatures(P636N)
P636N <- ScaleData(P636N, features = rownames(P636N))
#Anchors and mapping to REF
Anchors <- FindTransferAnchors(reference = PBM.Healthy, reference.assay = "RNA",
                               query = P636N,  reference.reduction = "pca", dims = 1:10)
predictions <- TransferData(anchorset = Anchors, refdata = PBM.Healthy$seurat_clusters, dims = 1:10)
P636N <- AddMetaData(P636N, metadata = predictions)
P636N <- MapQuery(anchorset = Anchors, reference = PBM.Healthy, query = P636N,
                  refdata = list(celltype = "seurat_clusters"), reference.reduction = "harmony", reduction.model = "umap")
#Count per cluster
table(P636N@meta.data$predicted.id)
ALL.P636N<- DimPlot(P636N, reduction = "ref.umap", group.by = "predicted.id",
                    pt.size = 1)  + NoLegend()+ ggtitle("P636N transferred labels")
#Gene expression
FeaturePlot(P636N, features = c("HMGA2","MECOM", "CD24", "CD19"),pt.size = 1,order = T)&
  scale_colour_gradientn(colours = (brewer.pal(name="PuBu",n=9)[2:9]))

#PB from P1977N
P1977N<-Read10X(data.dir="data_dir/scRNAseq_Datarequest/1977N")
P1977N<- CreateSeuratObject(counts= P1977N, project = "P1977N", min.cells = 3, min.features = 200)
P1977N[["percent.mt"]]<-PercentageFeatureSet(P1977N, pattern ="^MT-")
P1977N<- subset(P1977N, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(P1977N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P1977N <- NormalizeData(P1977N)
P1977N <- FindVariableFeatures(P1977N)
P1977N <- ScaleData(P1977N, features = rownames(P1977N))
#Anchors and mapping to REF
Anchors <- FindTransferAnchors(reference = PBM.Healthy, reference.assay = "RNA",
                               query = P1977N,  reference.reduction = "pca", dims = 1:10)
predictions <- TransferData(anchorset = Anchors, refdata = PBM.Healthy$seurat_clusters, dims = 1:10)
P1977N <- AddMetaData(P1977N, metadata = predictions)
P1977N <- MapQuery(anchorset = Anchors, reference = PBM.Healthy, query = P1977N,
                   refdata = list(celltype = "seurat_clusters"), reference.reduction = "harmony", reduction.model = "umap")
#Count per cluster
table(P1977N@meta.data$predicted.id)
ALL.P1977N<-DimPlot(P1977N, reduction = "ref.umap", group.by = "predicted.id",
                    pt.size = 1)  + NoLegend()+ ggtitle("P1977N transferred labels")
#Gene expression
FeaturePlot(P1977N, features = c("HMGA2","MECOM", "CD24", "CD19"),pt.size = 1,order = T)&
  scale_colour_gradientn(colours = (brewer.pal(name="PuBu",n=9)[2:9]))

#PB from P6086R
P6086R<-Read10X(data.dir="data_dir/scRNAseq_Datarequest/6086R")
P6086R<- CreateSeuratObject(counts= P6086R, project = "P6086R", min.cells = 3, min.features = 200)
P6086R[["percent.mt"]]<-PercentageFeatureSet(P6086R, pattern ="^MT-")
P6086R<- subset(P6086R, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(P6086R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P6086R <- NormalizeData(P6086R)
P6086R <- FindVariableFeatures(P6086R)
P6086R <- ScaleData(P6086R, features = rownames(P6086R))
#Anchors and mapping to REF
Anchors <- FindTransferAnchors(reference = PBM.Healthy, reference.assay = "RNA",
                               query = P6086R,  reference.reduction = "pca", dims = 1:10)
predictions <- TransferData(anchorset = Anchors, refdata = PBM.Healthy$seurat_clusters, dims = 1:10)
P6086R <- AddMetaData(P6086R, metadata = predictions)
P6086R <- MapQuery(anchorset = Anchors, reference = PBM.Healthy, query = P6086R,
                   refdata = list(celltype = "seurat_clusters"), reference.reduction = "harmony", reduction.model = "umap")
#Count per cluster
table(P6086R@meta.data$predicted.id)
ALL.P6086R<- DimPlot(P6086R, reduction = "ref.umap", group.by = "predicted.id",
                     pt.size = 1)  + NoLegend()+ ggtitle("P6086R transferred labels")
#Gene expression
FeaturePlot(P6086R, features = c("HMGA2","MECOM", "CD24", "CD19"),pt.size = 1,order = T)&
  scale_colour_gradientn(colours = (brewer.pal(name="PuBu",n=9)[2:9]))

#PB from P8010R
P8010R<-Read10X(data.dir="data_dir/scRNAseq_Datarequest/8010R")
P8010R<- CreateSeuratObject(counts= P8010R, project = "P8010R", min.cells = 3, min.features = 200)
P8010R[["percent.mt"]]<-PercentageFeatureSet(P8010R, pattern ="^MT-")
P8010R<- subset(P8010R, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(P8010R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P8010R <- NormalizeData(P8010R)
P8010R <- FindVariableFeatures(P8010R)
P8010R <- ScaleData(P8010R, features = rownames(P8010R))
#Anchors and mapping to REF
Anchors <- FindTransferAnchors(reference = PBM.Healthy, reference.assay = "RNA",
                               query = P8010R,  reference.reduction = "pca", dims = 1:10)
predictions <- TransferData(anchorset = Anchors, refdata = PBM.Healthy$seurat_clusters, dims = 1:10)
P8010R <- AddMetaData(P8010R, metadata = predictions)
P8010R <- MapQuery(anchorset = Anchors, reference = PBM.Healthy, query = P8010R,
                   refdata = list(celltype = "seurat_clusters"), reference.reduction = "harmony", reduction.model = "umap")
#Count per cluster
table(P8010R@meta.data$predicted.id)
ALL.P8010R<-DimPlot(P8010R, reduction = "ref.umap", group.by = "predicted.id",
                    pt.size = 1)  + NoLegend()+ ggtitle("P8010R transferred labels")
#Gene expression
FeaturePlot(P8010R, features = c("HMGA2","MECOM", "CD24", "CD19"),pt.size = 1,order = T)&
  scale_colour_gradientn(colours = (brewer.pal(name="PuBu",n=9)[2:9]))

####For gene Signature via UCell####
#48 up, 8 down gene signature for leukemic pre-proB: Charlotta Böiers
##"KMT2A" was removed to avoid gene-expression bias
LeukBStem<- list (LeukBStem=c("ARID2","BAHCC1","BANK1","CAMK2D","CDK17","CSNK1G3","DIAPH3",
                              "DLG2","ELF1","FOXP1","GEM","HMCN1","HMGA2","HOXA10","HOXA5",
                              "HOXA7","HOXA9","HSPA4","IGF1","IL9R","KANSL1L","MACROD2",
                              "MAGI1","MAN1A1","MCTP1","MECOM","MEIS1","MTMR3","NKX2-3","NRIP1",
                              "PCDH7","PDE3B","PHTF2","PLCB4","PPP1R9B","PRKG1","RANBP9","RNF220",
                              "RSBN1L","SLIT2","SOX4","THSD4","TOX","UTRN",
                              "PPIA-","COX4I1-","NDUFB10-","RTRAF-","MACROH2A1-","MYC-","FAU-","ATP6V0C-"))

my_palette <- brewer.pal(name="RdPu",n=9)[3:9]
sc<-scale_color_gradientn(limits = c(0,0.35),colors=my_palette )

#Different color palette for the gene expression

my_palette2 <- brewer.pal(name="BuPu",n=9)[2:9]
sc2<-scale_color_gradientn(colors= my_palette2)

#signature in HEALTHY           
PBM.Healthy<-AddModuleScore_UCell(PBM.Healthy,features = LeukBStem,name = NULL)
FeaturePlot(PBM.Healthy,features = names(LeukBStem),pt.size = 1.5,order = T)+ ggtitle("Pediatric healthy BM")&
  sc

#HMGA2 expression check for signature pathway
FeaturePlot(PBM.Healthy, features = c("HMGA2", "LIN28B"),pt.size = 2,order = T)&
  sc2

#Signature in KMT2A::AFF1 Patient Samples

#P636N
P636N<-AddModuleScore_UCell(P636N,features = LeukBStem,name = NULL)
FeaturePlot(P636N, features = names(LeukBStem),pt.size = 1.5,order = T)+ ggtitle("P636N BM")&
  sc

#HMGA2 expression check for signature pathway
FeaturePlot(P636N, features = c("HMGA2"),pt.size = 2,order = T)&
  sc2

#P1977N
P1977N<-AddModuleScore_UCell(P1977N,features = LeukBStem,name = NULL)
FeaturePlot(P1977N, features = names(LeukBStem),pt.size = 1.5,order = T)+ ggtitle("P1977N BM")&
  sc

#HMGA2 expression check for signature pathway
FeaturePlot(P1977N, features = c("HMGA2"),pt.size = 2,order = T)&
  sc2

#P6086R
P6086R<-AddModuleScore_UCell(P6086R,features = LeukBStem,name = NULL)
FeaturePlot(P6086R, features = names(LeukBStem),pt.size = 1.5,order = T)+ ggtitle("P6086R BM")&
  sc


#HMGA2 expression check for signature pathway
FeaturePlot(P6086R, features = c("HMGA2"),pt.size = 2,order = T)&
  sc2

#P8010R
P8010R<-AddModuleScore_UCell(P8010R,features = LeukBStem,name = NULL)
FeaturePlot(P8010R, features = names(LeukBStem),pt.size = 1.5,order = T)+ ggtitle("P8010R BM")&
  sc


#HMGA2 expression check for signature pathway
FeaturePlot(P8010R, features = c("HMGA2"),pt.size = 2,order = T)&
  sc2
