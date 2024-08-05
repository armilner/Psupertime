library(Seurat)
library(dplyr)
obj1K <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSeurat.rds")
obj2K <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/2KSeurat.rds")

merge <- merge(obj1K, y = obj2K)
rm(obj1K, obj2K)
gc()
merge <- SCTransform(merge, vst.flavor = "v2", vars.to.regress = "percent.mt",  verbose = F, conserve.memory = T)
VariableFeatures(merge[["SCT"]]) <- rownames(merge[["SCT"]]@scale.data)
merge <- RunPCA(merge, verbose = F)

ElbowPlot(merge)

#Harmony integration
merge <- IntegrateLayers(
  object = merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
#Unintegrated

merge <- FindNeighbors(merge, dims = 1:15, reduction = "pca", verbose = F)
merge <- FindClusters(merge, resolution = 0.5, cluster.name = "unintegrated", verbose = F)
merge <- RunUMAP(merge, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated", verbose = F)

#Harmony

merge <- FindNeighbors(merge, dims = 1:15, reduction = "harmony", verbose = F)
merge <- FindClusters(merge, resolution = 0.5, cluster.name = "harmony_clusters", verbose = F)
merge <- RunUMAP(merge, dims = 1:19, reduction = "harmony", reduction.name = "umap.harmony", verbose = F)

saveRDS(merge, file = "../../nfs_share/Milner/Psupertime/Allcombined/Objects/allsamplesmergeprocessUnannotated.rds")

#Quickly Annotate-Will need to redo DEGs

DimPlot(merge, reduction = "umap.harmony", label = T)

DimPlot(merge, reduction = "umap.unintegrated", label = T)


#FAILURE....######################################################################################

#read in DF samples, and process

library(Seurat)

#17.5 samples
S1_1K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/176C1KDF.rds")
S2_1K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/725D1KDF.rds")
S3_1K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/477D1KDF.rds")

#Add Time Label
S1_1K_17$time <- "17.5"
S2_1K_17$time <- "17.5"
S3_1K_17$time <- "17.5"

#16.5 labels
S1_1K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/167E1KDF.rds")
S2_1K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/235I1KDF.rds")
S3_1K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/723D1KDF.rds")

#Add Time Label

S1_1K_16$time <- "16.5"
S2_1K_16$time <- "16.5"
S3_1K_16$time <- "16.5"


#15.5 Labels
S1_1K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/177C895E1KDF.rds")
S2_1K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/188G233A1KDF.rds")
S3_1K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/189F824B1KDF.rds")

#Add Time Label

S1_1K_15$time <- "15.5"
S2_1K_15$time <- "15.5"
S3_1K_15$time <- "15.5"

#2K Samples

#17.5 samples
S4_2K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/176E2KDF.rds")
S5_2K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/179A2KDF.rds")
S6_2K_17 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/17.5/837A2KDF.rds")

#Add Time Label (M
S4_2K_17$time <- "17.5"
S5_2K_17$time <- "17.5"
S6_2K_17$time <- "17.5"

#16.5 labels
S4_2K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/178E2KDF.rds")
S5_2K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/200J2KDF.rds")
S6_2K_16 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/16.5/993F2KDF.rds")

#Add Time Label

S4_2K_16$time <- "16.5"
S5_2K_16$time <- "16.5"
S6_2K_16$time <- "16.5"


#15.5 Labels
S4_2K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/188D189G2KDF.rds")
S5_2K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/189G187C2KDF.rds")
S6_2K_15 <- readRDS(file = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/DoubletFinder/15.5/232I231D2KDF.rds")

#Add Time Label

S4_2K_15$time <- "15.5"
S5_2K_15$time <- "15.5"
S6_2K_15$time <- "15.5"


merge <- merge(S1_1K_15, y = c(S1_1K_16, S1_1K_17, S2_1K_15, S2_1K_16, S2_1K_17, S3_1K_15, S3_1K_16, S3_1K_17, S4_2K_15, S4_2K_16, S4_2K_17, S5_2K_15, S5_2K_16,
                                S5_2K_17, S6_2K_15, S6_2K_16, S6_2K_17))

rm(S1_1K_15, S1_1K_16, S1_1K_17, S2_1K_15, S2_1K_16, S2_1K_17, S3_1K_15, S3_1K_16, S3_1K_17, S4_2K_15, S4_2K_16, S4_2K_17, S5_2K_15, S5_2K_16,
   S5_2K_17, S6_2K_15, S6_2K_16, S6_2K_17)

saveRDS(merge, file = "../../nfs_share/Milner/Psupertime/Allcombined/Objects/mergedsingleDFsamples.rds")

#Processing 

VariableFeatures(merge[["SCT"]]) <- rownames(merge[["SCT"]]@scale.data)

merge <- RunPCA(merge)

ElbowPlot(merge)


#Harmony integration
merge <- IntegrateLayers(
  object = merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
#Unintegrated

merge <- FindNeighbors(merge, dims = 1:15, reduction = "pca", verbose = F)
merge <- FindClusters(merge, resolution = 0.6, cluster.name = "unintegrated", verbose = F)
merge <- RunUMAP(merge, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated", verbose = F)

#Harmony

merge <- FindNeighbors(merge, dims = 1:19, reduction = "harmony", verbose = F)
merge <- FindClusters(merge, resolution = 0.5, cluster.name = "harmony_clusters", verbose = F)
merge <- RunUMAP(merge, dims = 1:19, reduction = "harmony", reduction.name = "umap.harmony", verbose = F)


saveRDS(merge, file = "../../nfs_share/Milner/Psupertime/Allcombined/Objects/allsamplesmergeprocessUnannotated.rds")

#Quickly Annotate-Will need to redo DEGs

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/UMAPharmony.png")
DimPlot(merge, reduction = "umap.harmony", label = T)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/UMAPunit.png")
DimPlot(merge, reduction = "umap.unintegrated", label = T)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/allSamplesIntegratedAcrosstime.png", width = 960)
DimPlot(merge, reduction = "umap.harmony", label = T, split.by = "time")
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/allSamplesIntegratedUnlabeledGroupSample.png", width = 960)
DimPlot(merge, reduction = "umap.harmony", label = T, split.by = "time", group.by = "orig.ident")
dev.off()


#subsetting re run clustering to look at harmony int umap

merge <- subset(merge, subset = nFeature_RNA > 1000 & nFeature_RNA < 2500 & percent.mt < 5)

#Here I reran the find neighbors, find clusters...

#Remove small unnecessary clusters (15,19,16,20,18,17)

merge <- subset(merge, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))


#FeaturePlots to label cell clusters

FeaturePlot(merge, reduction = "umap.harmony", features = "Lrp2", label = T, split.by = "time")
FeaturePlot(merge, reduction = "umap.harmony", features = "Eya1", label = T, split.by = "time")
FeaturePlot(merge, reduction = "umap.harmony", features = "Six2", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Pax2", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Ret", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Gdnf", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Wnt11", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Gata3", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Pecam1", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Kdr", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Emcn", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Aqp3", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Scnn1b", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Pdgfrb", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Foxd1", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Lgr5", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Zeb2", label = T)

FeaturePlot(merge, reduction = "umap.harmony", features = "Sfrp1", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Apoe", label = T)
FeaturePlot(merge, reduction = "umap.harmony", features = "Slc34a1", label = T)

FeaturePlot(merge, reduction = "umap.harmony", features = "Gpc6", label = T)

FeaturePlot(merge, reduction = "umap.harmony", features = "Nphs1", label = T)


#REnaming clusters

merge <- readRDS(file = "../../nfs_share/Milner/Psupertime/Allcombined/Objects/allsamplesmergeprocessUnannotated.rds")



merge <- RenameIdents(merge, '0' = "MM")
merge <- RenameIdents(merge, '1' = "UB")
merge <- RenameIdents(merge, '2' = "Pdgfrb+ Stromal")
merge <- RenameIdents(merge, '3' = "Foxd1+ Stromal")
merge <- RenameIdents(merge, '4' = "MM")
merge <- RenameIdents(merge, '5' = "Tubule Progenitor")
merge <- RenameIdents(merge, '6' = "Collecting Duct")
merge <- RenameIdents(merge, '7' = "Pdgfrb+ Stromal")
merge <- RenameIdents(merge, '8' = "Foxd1+ Stromal")
merge <- RenameIdents(merge, '9' = "Pdgfrb+ Stromal")
merge <- RenameIdents(merge, '10' = "Foxd1+ Stromal")
merge <- RenameIdents(merge, '11' = "MM")
merge <- RenameIdents(merge, '12' = "Podocytes")
merge <- RenameIdents(merge, '13' = "Transitory")
merge <- RenameIdents(merge, '14' = "Endothelium")

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/AllsampHarmonyAnnotated.png")
DimPlot(merge, reduction = "umap.harmony", label = T)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/AllsampHarmonyAnnotatedsplitbytime.png", width = 960)
DimPlot(merge, reduction = "umap.harmony", label = T, split.by = "time", repel = T)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/AllsampHarmonyAnnotatedsplitbytimegroupbtsample.png", width = 960)
DimPlot(merge, reduction = "umap.harmony", label = T, split.by = "time", group.by = "orig.ident", repel = T)
dev.off()

saveRDS(merge, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/Objects/allsamplesHarmAnnotated.rds")

merge <- readRDS(file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/Objects/allsamplesHarmAnnotated.rds")

DefaultAssay(merge) <- "RNA"

merge <- NormalizeData(merge)

#Join layers before DE analysis-prior to scaling for RAM

merge <- JoinLayers(merge)

merge <- ScaleData(merge)

saveRDS(merge, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/Objects/allsamplesHarmonyAnnotatedScaledforDEA.rds")

#####PICK UP HERE WITH SCALED FILE FOR DEA, need to add extra annotation column to perform DEA.

####################################################################################################

#Differential Expression Analysis
#Compare 1K vs 2K for each celltype/timepoint


merge$celltype <- Idents(merge)
merge$celltype.group <- paste(merge$celltype, merge$orig.ident, sep = "_")
merge$celltype.group.tp <- paste(merge$celltype.group, merge$time, sep = "_")

Idents(merge) <- "celltype.group.tp"

#Start with 15.5

#Celltype: MM, Podocytes, Tubule Progenitor, UB, Collecting Duct, Foxd1+ Stromal, Pdgfrb+ Stromal, Transitory, Endothelium

MM_15 <- FindMarkers(merge, ident.1 = "MM_1K_15.5", ident.2 = "MM_2K_15.5")

#Error on Podocytes to few cells
Podo_15 <- FindMarkers(merge, ident.1 = "Podocytes_1K_15.5", ident.2 = "Podocytes_2K_15.5")

Tubprog_15 <- FindMarkers(merge, ident.1 = "Tubule Progenitor_1K_15.5", ident.2 = "Tubule Progenitor_2K_15.5")
UB_15 <- FindMarkers(merge, ident.1 = "UB_1K_15.5", ident.2 = "UB_2K_15.5")
CD_15 <- FindMarkers(merge, ident.1 = "Collecting Duct_1K_15.5", ident.2 = "Collecting Duct_2K_15.5")
FxStromal_15 <- FindMarkers(merge, ident.1 = "Foxd1+ Stromal_1K_15.5", ident.2 = "Foxd1+ Stromal_2K_15.5")
PdStromal_15 <- FindMarkers(merge, ident.1 = "Pdgfrb+ Stromal_1K_15.5", ident.2 = "Pdgfrb+ Stromal_2K_15.5")
Trans_15 <- FindMarkers(merge, ident.1 = "Transitory_1K_15.5", ident.2 = "Transitory_2K_15.5")
Endo_15 <- FindMarkers(merge, ident.1 = "Endothelium_1K_15.5", ident.2 = "Endothelium_2K_15.5")

#save as csv file

write.csv(MM_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/MM15.csv")
write.csv(Tubprog_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Tubprog15.csv")
write.csv(UB_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/UB15.csv")
write.csv(CD_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/CD15.csv")
write.csv(FxStromal_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Fxstromal15.csv")
write.csv(PdStromal_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Pdstromal15.csv")
write.csv(Trans_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Trans15.csv")
write.csv(Endo_15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Endo15.csv")

#16.5

MM_16 <- FindMarkers(merge, ident.1 = "MM_1K_16.5", ident.2 = "MM_2K_16.5")
Podo_16 <- FindMarkers(merge, ident.1 = "Podocytes_1K_16.5", ident.2 = "Podocytes_2K_16.5")
Tubprog_16 <- FindMarkers(merge, ident.1 = "Tubule Progenitor_1K_16.5", ident.2 = "Tubule Progenitor_2K_16.5")
UB_16 <- FindMarkers(merge, ident.1 = "UB_1K_16.5", ident.2 = "UB_2K_16.5")
CD_16 <- FindMarkers(merge, ident.1 = "Collecting Duct_1K_16.5", ident.2 = "Collecting Duct_2K_16.5")
FxStromal_16 <- FindMarkers(merge, ident.1 = "Foxd1+ Stromal_1K_16.5", ident.2 = "Foxd1+ Stromal_2K_16.5")
PdStromal_16 <- FindMarkers(merge, ident.1 = "Pdgfrb+ Stromal_1K_16.5", ident.2 = "Pdgfrb+ Stromal_2K_16.5")
Trans_16 <- FindMarkers(merge, ident.1 = "Transitory_1K_16.5", ident.2 = "Transitory_2K_16.5")
Endo_16 <- FindMarkers(merge, ident.1 = "Endothelium_1K_16.5", ident.2 = "Endothelium_2K_16.5")

#Write to csv

write.csv(MM_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/MM16.csv")
write.csv(Tubprog_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Tubprog16.csv")
write.csv(UB_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/UB16.csv")
write.csv(CD_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/CD16.csv")
write.csv(FxStromal_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Fxstromal16.csv")
write.csv(PdStromal_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Pdstromal16.csv")
write.csv(Trans_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Trans16.csv")
write.csv(Endo_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Endo16.csv")
write.csv(Podo_16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Podo16.csv")


#17.5

MM_17 <- FindMarkers(merge, ident.1 = "MM_1K_17.5", ident.2 = "MM_2K_17.5")
Podo_17 <- FindMarkers(merge, ident.1 = "Podocytes_1K_17.5", ident.2 = "Podocytes_2K_17.5")
Tubprog_17 <- FindMarkers(merge, ident.1 = "Tubule Progenitor_1K_17.5", ident.2 = "Tubule Progenitor_2K_17.5")
UB_17 <- FindMarkers(merge, ident.1 = "UB_1K_17.5", ident.2 = "UB_2K_17.5")
CD_17 <- FindMarkers(merge, ident.1 = "Collecting Duct_1K_17.5", ident.2 = "Collecting Duct_2K_17.5")
FxStromal_17 <- FindMarkers(merge, ident.1 = "Foxd1+ Stromal_1K_17.5", ident.2 = "Foxd1+ Stromal_2K_17.5")
PdStromal_17 <- FindMarkers(merge, ident.1 = "Pdgfrb+ Stromal_1K_17.5", ident.2 = "Pdgfrb+ Stromal_2K_17.5")
Trans_17 <- FindMarkers(merge, ident.1 = "Transitory_1K_17.5", ident.2 = "Transitory_2K_17.5")
Endo_17 <- FindMarkers(merge, ident.1 = "Endothelium_1K_17.5", ident.2 = "Endothelium_2K_17.5")

write.csv(MM_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/MM17.csv")
write.csv(Tubprog_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Tubprog17.csv")
write.csv(UB_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/UB17.csv")
write.csv(CD_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/CD17.csv")
write.csv(FxStromal_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Fxstromal17.csv")
write.csv(PdStromal_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Pdstromal17.csv")
write.csv(Trans_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Trans17.csv")
write.csv(Endo_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Endo17.csv")
write.csv(Podo_17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Podo17.csv")

#Global Degs

merge$group.tp <- paste(merge$orig.ident, merge$time, sep = "_")

Idents(merge) <- "group.tp"

global15 <- FindMarkers(merge, ident.1 = "1K_15.5", ident.2 = "2K_15.5")
global16 <- FindMarkers(merge, ident.1 = "1K_16.5", ident.2 = "2K_16.5")
global17 <- FindMarkers(merge, ident.1 = "1K_17.5", ident.2 = "2K_17.5")

write.csv(global15, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/global15.csv")
write.csv(global16, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/global16.csv")
write.csv(global17, file = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/global17.csv")


#Making Barplot of DEGs

#DEGs/celltype adj.pval 0.05

process_DEGs <- function(file_path) {
  deg_data <- read.csv(file_path)
  significant_degs <- deg_data %>%
    filter(p_val_adj < 0.05) %>%
    mutate(Direction = ifelse(avg_log2FC > 0, "Upregulated", "Downregulated")) %>%
    count(Direction)
  return(significant_degs)
}



degs <- process_DEGs("../../nfs_share/Milner/")
degs$Cluster <- ""


#17.5 Barplot


MM_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/MM17.csv")
MM_17$Cluster <- "MM"

Podo_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Podo17.csv")
Podo_17$Cluster <- "Podocyte"

Tubprog_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Tubprog17.csv")
Tubprog_17$Cluster <- "Tubule Progenitors"

UB_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/UB17.csv")
UB_17$Cluster <- "UB"

CD_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/CD17.csv")
CD_17$Cluster <- "Collecting Duct"

FxStromal_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Fxstromal17.csv")
FxStromal_17$Cluster <- "Foxd1+ Stromal"

PdStromal_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Pdstromal17.csv")
PdStromal_17$Cluster <- "Pdgfrb+ Stromal"

Trans_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Trans17.csv")
Trans_17$Cluster <- "Transitory"

Endo_17 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/Endo17.csv")
Endo_17$Cluster <- "Endothelium"

combined_degs <- bind_rows(MM_17, Podo_17, Tubprog_17, UB_17, CD_17, FxStromal_17, PdStromal_17, Trans_17, Endo_17)

level_order <- c("MM", "Podocyte", "Tubule Progenitors", "UB", "Collecting Duct", "Foxd1+ Stromal", "Pdgfrb+ Stromal", "Transitory", "Endothelium")

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/17/image/HSRAFK175DEGSbarplot.png", res = 85)
ggplot(combined_degs, aes(x = factor(Cluster, levels = level_order), y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  coord_flip() +
  labs(title = "GD 17.5 (-S vs. -C)", x = "Cluster", y = "Number of Significant DEGs", fill = "Gene Regulation") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



#16.5 Barplot


MM_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/MM16.csv")
MM_16$Cluster <- "MM"

Podo_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Podo16.csv")
Podo_16$Cluster <- "Podocyte"

Tubprog_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Tubprog16.csv")
Tubprog_16$Cluster <- "Tubule Progenitors"

UB_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/UB16.csv")
UB_16$Cluster <- "UB"

CD_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/CD16.csv")
CD_16$Cluster <- "Collecting Duct"

FxStromal_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Fxstromal16.csv")
FxStromal_16$Cluster <- "Foxd1+ Stromal"

PdStromal_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Pdstromal16.csv")
PdStromal_16$Cluster <- "Pdgfrb+ Stromal"

Trans_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Trans16.csv")
Trans_16$Cluster <- "Transitory"

Endo_16 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/Endo16.csv")
Endo_16$Cluster <- "Endothelium"

combined_degs <- bind_rows(MM_16, Podo_16, Tubprog_16, UB_16, CD_16, FxStromal_16, PdStromal_16, Trans_16, Endo_16)

level_order <- c("MM", "Podocyte", "Tubule Progenitors", "UB", "Collecting Duct", "Foxd1+ Stromal", "Pdgfrb+ Stromal", "Transitory", "Endothelium")

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/16/image/HSRAFK165DEGSbarplot.png", res = 85)
ggplot(combined_degs, aes(x = factor(Cluster, levels = level_order), y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  coord_flip() +
  labs(title = "GD 16.5 (-S vs. -C)", x = "Cluster", y = "Number of Significant DEGs", fill = "Gene Regulation") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#15.5 Barplot


MM_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/MM15.csv")
MM_15$Cluster <- "MM"

Podo_15 <- Trans_16
Podo_15$Cluster <- "Podocyte"

Tubprog_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Tubprog15.csv")
Tubprog_15$Cluster <- "Tubule Progenitors"

UB_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/UB15.csv")
UB_15$Cluster <- "UB"

CD_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/CD15.csv")
CD_15$Cluster <- "Collecting Duct"

FxStromal_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Fxstromal15.csv")
FxStromal_15$Cluster <- "Foxd1+ Stromal"

PdStromal_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Pdstromal15.csv")
PdStromal_15$Cluster <- "Pdgfrb+ Stromal"

Trans_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Trans15.csv")
Trans_15$Cluster <- "Transitory"

Endo_15 <- process_DEGs("../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/Endo15.csv")
Endo_15$Cluster <- "Endothelium"

combined_degs <- bind_rows(MM_15, Podo_15, Tubprog_15, UB_15, CD_15, FxStromal_15, PdStromal_15, Trans_15, Endo_15)

level_order <- c("MM", "Podocyte", "Tubule Progenitors", "UB", "Collecting Duct", "Foxd1+ Stromal", "Pdgfrb+ Stromal", "Transitory", "Endothelium")

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/DEA/15/image/HSRAFK155DEGSbarplot.png", res = 85)
ggplot(combined_degs, aes(x = factor(Cluster, levels = level_order), y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  coord_flip() +
  labs(title = "GD 15.5 (-S vs. -C)", x = "Cluster", y = "Number of Significant DEGs", fill = "Gene Regulation") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

Idents(merge) <- "celltype"

FeaturePlot(merge, features = "Smad2")

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/ViolinplotMarkers.png")
VlnPlot(merge, features = c("Emcn", "Gata3", "Aqp2", "Foxd1", "Pdgfrb", "Eya1", "Nphs1", "Lrp2"), stack = T, sort = 'increasing')
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Allcombined/Subsetted QC and Clusters/Annotated/ViolinplotMarkersFlip.png", width = 720)
VlnPlot(merge, features = c("Emcn", "Gata3", "Aqp2", "Foxd1", "Pdgfrb", "Eya1", "Nphs1", "Lrp2"), stack = T, sort = 'increasing', flip = T)
dev.off()

#Heatmap Top 10 unique genes per cluster
#Do conserved markers for each cluster

allmarkers <- FindConservedMarkers(merge, ident.1 = "celltype", grouping.var = "orig.ident")


#HEATMAP NEEDS TWEAKING
#Heatmap top 10 markers

levels(obj) <- c("Podocytes", "CD-PC", "MC", "CD-ICB", "PT-S", "CT", "ATL-LOH", "DT-C", "Proliferating", "GC", "CD-ICA", "DTL-LOH", "TAL-LOH", "PT-C")
all_top_genes <- c()

#Reorder to set decreasing
csv_files <- c("../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/Podo.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/PC.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/Mes.csv","../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/ICB.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/PTS.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/CNT.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/ATL.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/DCT.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/Prol.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/GlomEndo.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/ICA.csv",  "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/DTL.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/TAL.csv", "../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Markers/PTC.csv")


for (file_path in csv_files) {
  markers <- read.csv(file_path)
  top_genes <- head(markers$X, 10)
  all_top_genes <- union(all_top_genes, top_genes)
} 

png("../../nfs_share/Milner/HSRA4 Unforced/with7and8no5/Figures/heatmapHSRA.png", width = 900, height = 1350)
DoHeatmap(obj, features = all_top_genes, label = T, angle = 90, size = 7.5) + NoLegend()
dev.off()  
# + theme(axis.text.y = element_text(size = 5) #to change gene name size


