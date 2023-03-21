## Project: To explore the EFF induced cell state of HADF cells at single cell resolution
```
author: "Rajneesh Srivastava"
#dates: 
	#original: 04/20/2022
	#last updated: 03/01/2023
```
#### Load libraries
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(future)
library(sctransform)
plan("multisession", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
```
#### UPLOAD DATA
```
setwd("path/EFF_project/")
EFF.dir = "path/EFF_project/dataset/"
EFF.list=list("CHDAF", "EFFCSD5", "EFFTET123D5")

for (file in EFF.list){
               EFF_data <- Read10X_h5(filename =    
                                    paste0(EFF.dir, file, "/",
                                           "filtered_feature_bc_matrix.h5"))
               
               EFF_obj <- CreateSeuratObject(counts = 
                                    EFF_data,
                                    min.cells=3,
                                    min.features = 300,
                                    project = file)
               assign(file, EFF_obj)
                            }
```
#### QUALITY FILTERING AND VISUALIZATION
```
sample.list=list(CHDAF, EFFCSD5, EFFTET123D5)

for (i in 1:length(sample.list)) {
    sample.list[[i]][["percent.mt"]] <-
            PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
    sample.list[[i]] <- subset(sample.list[[i]], 
                            subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 5000 & percent.mt < 15 & 
                            nCount_RNA < 25000 & nCount_RNA > 2000)     
#SCT transformation
     sample.list[[i]] <- SCTransform(sample.list[[i]],
                                    vars.to.regress = "percent.mt", 
                                    verbose = FALSE)
                                 }
#rm(CHDAF,EFFCSD5, EFFTET123D5, EFF.dir,EFF_data,EFF_obj,EFF.list,file)
```
##### Generate violin plots
```
VD.list=list()

for (i in 1:length(sample.list)) {
    VD.list[[i]] <- VlnPlot(sample.list[[i]], 
                            features = c("nFeature_SCT",       
                                          "nCount_SCT",
                                          "percent.mt"), 
                                          ncol = 3)
                                }

VD.list                        
```
#### DATA INTEGRATION
```
for (i in 1:length(sample.list)) {
     sample.list[[i]] <- SCTransform(sample.list[[i]],
                                    vars.to.regress = "percent.mt", 
                                    verbose = FALSE)
                                 }
sample.features <- SelectIntegrationFeatures(object.list = sample.list,    
                   nfeatures = 3000)

sample.list <- lapply(X = sample.list, FUN = function(x) {
                   x <- RunPCA(x, features = sample.features)
                                                         } )

sample.list <- PrepSCTIntegration(object.list = sample.list, 
                   anchor.features = sample.features,
                   verbose = FALSE)
                   
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, 
                   normalization.method = "SCT", 
                   anchor.features = sample.features, 
                   reference = c(1, 2), 
                   reduction = "rpca", 
                   verbose = TRUE)

sample.integrated <- IntegrateData(anchorset = sample.anchors, 
                   normalization.method = "SCT",
                   verbose = TRUE)
```
#### CLUSTERING ANALYSIS
```
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 
#run tSNE/UMAP
sample.integrated <- FindNeighbors(sample.integrated) # dim=1:10
sample.integrated = RunUMAP(sample.integrated, dims = 1:30)
sample.integrated = RunTSNE(sample.integrated, dims = 1:30)

#saveRDS (sample.integrated, file = "sample.integrated.rds")
```
#### FIND CLUSTERS WITH DEFINED RESOLUTION
```
sample15 <- FindClusters(sample.integrated, resolution = 0.15)
#saveRDS (sample15, file = "sample15.rds")
#sample15=readRDS("sample15.rds")
```

#### DOWN SAMPLING
```
#sample15=readRDS("path/sample15.rds")
#Idents(sample15)="seurat_clusters"

CHDAF=subset(sample15,subset=orig.ident=="CHDAF")
nonCHDAF=subset(sample15,subset=orig.ident!="CHDAF")

sd_CHDAF = subset(CHDAF, cells = sample(Cells(CHDAF), 2000))
sample15 <- merge(nonCHDAF, y = sd_CHDAF, merge.data = TRUE, add.cell.ids = c("trns", "sd_ctrl"),merge.dr=c("pca","tsne","umap"))
```
##### tweak-in IDs
```
#meta=read.table("metadata.txt",sep="\t", header=T)
head(meta)
      orig.id              group
1       CHDAF            control
2     EFFCSD5       EFF+scramble
3 EFFTET123D5 EFF+TET1/2/3 siRNA

GSM=sample15@meta.data
GSM$group=1
#head(GSM)

for(i in 1:nrow(GSM)){
  for (j in 1:nrow(meta)){
   if (GSM[i,1] == meta[j,1])
        { GSM[i,9] = meta[j,2] }
                 }
             }

sample15@meta.data=GSM
#saveRDS (sample15, file = "sample15.rds")
```
##### Cell distribution across groups, and expression profile of Fb markers
```
DimPlot(sample15, split.by="group") #Figure 2A
VlnPlot(sample15,features=c("THY1","S100A4","VIM","COL1A1","COL1A2", "FBLN1"),flip=T,stack=T) #Figure 2B
```
#### DE analysis for all clusters (related to Figure 2C)
```
DefaultAssay(sample15)="SCT"
x=c("0","1","2","3")

for (i in (x)){
		sigma <- FindMarkers(sample15, 
									 ident.1 = i, 
									 ident.2 = dput(paste0(x[!(x %in% i)])),   #i.e. ident.2=i not in x
									 logfc.threshold = 0.10,
									 min.pct=0.05,recorrect_umi=FALSE)
		write.csv(sigma,paste0(i,"-vs-others",".csv")) 
				  }
```
##### Dot plot (Figure S1A)
```
DefaultAssay(sample25)="SCT"
TMPA_1.0=c("ANTXR1", "ANTXR2", "B2M", "CAV1", "CD320", "CD46", "CD47", "COLEC12", "ENG", "F3", "HLA−B", "HLA−E", "HMMR", "IFNAR1", "IFNGR1", "IGF1R", "IGF2R", "IL13RA1", "IL13RA2", "IL1R1", "IL6ST", "ITGA2", "ITGA4", "ITGA5", "ITGAV", "ITGB1", "KTN1", "LGALS3BP", "LRP1", "LRP10", "NRP1", "OSMR", "PGRMC1", "PLA2R1", "PLAUR", "PLXNA1", "PLXNA3", "SCARB2", "SIGMAR1", "THBD", "TNFRSF12A", "TNFRSF1A", "TSPO", "UTRN")
DotPlot(sample15, features = unique(c(TMPA_1.0)))+ 
     theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
     theme(axis.text.y=element_text(size=7.5, face="italic"))
```
#### CELL CHAT ANALYSIS
```
#load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(future)
library(sctransform)
plan("multisession", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
```
##### INPUT Seurat object
```
DefaultAssay(sample15)="SCT"
G2=subset(sample,subset=group=="EFF+scramble")
G3=subset(sample,subset=group=="EFF+TET1/2/3 siRNA")

# Create CellChat objects
G2_cell_chat=createCellChat(object=G2,meta=G2@meta.data,group.by="seurat_clusters")
G3_cell_chat=createCellChat(object=G3,meta=G3@meta.data,group.by="seurat_clusters")

# Assign Idents
G2_cell_chat <-setIdent(G2_cell_chat, ident.use="seurat_clusters")
G3_cell_chat <-setIdent(G3_cell_chat, ident.use="seurat_clusters")

#Database usage
CellChatDB<- CellChatDB.human
#showDatabaseCategory(CellChatDB)
CellChatDB.use<- CellChatDB
G2_cell_chat@DB <-CellChatDB.use
G3_cell_chat@DB <-CellChatDB.use

#subset the dataset for enhanced performance
G2_cell_chat <-subsetData(G2_cell_chat)
G3_cell_chat <-subsetData(G3_cell_chat)
```
#### Individual Analysis
```
G2_cell_chat <- identifyOverExpressedGenes(G2_cell_chat)
G2_cell_chat <- identifyOverExpressedInteractions(G2_cell_chat)
G2_cell_chat <- projectData(G2_cell_chat, PPI.human) #(optional)
G2_cell_chat <- computeCommunProb(G2_cell_chat, raw.use = TRUE)
G2_cell_chat <- filterCommunication(G2_cell_chat,min.cells = 10)
G2_cell_chat <- computeCommunProbPathway(G2_cell_chat)
G2_cell_chat <- netAnalysis_computeCentrality(G2_cell_chat, slot.name = "netP")
G2_cell_chat <- aggregateNet(G2_cell_chat)

G3_cell_chat <- identifyOverExpressedGenes(G3_cell_chat)
G3_cell_chat <- identifyOverExpressedInteractions(G3_cell_chat)
G3_cell_chat <- projectData(G3_cell_chat, PPI.human) #(optional)
G3_cell_chat <- computeCommunProb(G3_cell_chat, raw.use = TRUE)
G3_cell_chat <- filterCommunication(G3_cell_chat,min.cells = 10)
G3_cell_chat <- computeCommunProbPathway(G3_cell_chat)
G3_cell_chat <- netAnalysis_computeCentrality(G3_cell_chat, slot.name = "netP")
G3_cell_chat <- aggregateNet(G3_cell_chat)
```
##### Merging cell chats
```
object.list=c(group2=G2_cell_chat,group3=G3_cell_chat)
merged_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
```
##### save cell_chat_objects
```
saveRDS(G2_cell_chat,"group2_cell_chat.rds")
saveRDS(G3_cell_chat,"group3_cell_chat.rds")
saveRDS(merged_cellchat, "merged_cellchat.rds")
```
##### Comparison of cell chat objects
```
#Figure S1B
g1=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2))
g2=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
g1+g2

#Figure S1C
#object.list=c(group2=G2_cell_chat,group3=G3_cell_chat)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[[i]]))
}

#Figure S1D - Information_flow
rankNet(merged_cellchat, comparison=c(1,2),stacked = T, do.stat = TRUE)

#Figure 2D DE-connectome
netVisual_diffInteraction(merged_cellchat, weight.scale = T,comparison = c(1,3))
```
##### Identify signaling groups based on functional similarity
```
#Figure 2E

library(reticulate)
reticulate::py_install(packages ='umap-learn')
UMAP<-import("umap", convert=FALSE)

#merged_cellchat=readRDS("./merged_cellchat")

cellchatf <- computeNetSimilarityPairwise(merged_cellchat, type = "functional")
cellchatf <- netEmbedding(cellchatf, type = "functional")
cellchatf <- netClustering(cellchatf, type = "functional")
netVisual_embeddingPairwise(cellchatf, type = "functional", label.size = 3.5)

#cellchats <- computeNetSimilarityPairwise(merged_cellchat, type = "structural")
#cellchats <- netEmbedding(cellchats, type = "structural")
#cellchats <- netClustering(cellchats, type = "structural")
#netVisual_embeddingPairwise(cellchats, type = "structural", label.size = 3.5)
```
#### Thank you
