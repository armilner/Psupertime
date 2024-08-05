library(Seurat)
library(SingleCellExperiment)
library(psupertime)

#1K first

#E15.5 Data, run 1 at a time and save after
######################################################################################################################################
oneK_15_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_1K_177C895E_LK/filtered_feature_bc_matrix/")
oneK_15_S1 <- CreateSeuratObject(counts = oneK_15_S1_data, project = "1K", min.cells = 3)

o1K_15_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_1K_188G233A_LK/filtered_feature_bc_matrix/")
o1K_15_S2<- CreateSeuratObject(counts = o1K_15_S2_data, project = "1K", min.cells = 3)

o1K_15_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_1K_189F824B_RK/filtered_feature_bc_matrix/")
o1K_15_S3 <- CreateSeuratObject(counts = o1K_15_S3_data, project = "1K", min.cells = 3)

#Add Time Label

oneK_15_S1$time <- "15.5"
o1K_15_S2$time <- "15.5"
o1K_15_S3$time <- "15.5"

#E16.5
######################################################################################################################################
o1K_16_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_1K_167E_LK/filtered_feature_bc_matrix/")
o1K_16_S1 <- CreateSeuratObject(counts = o1K_16_S1_data, project = "1K", min.cells = 3)


o1K_16_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_1K_235I_LK/filtered_feature_bc_matrix/")
o1K_16_S2 <- CreateSeuratObject(counts = o1K_16_S2_data, project = "1K", min.cells = 3)


o1K_16_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_1K_723D_RK/filtered_feature_bc_matrix/")
o1K_16_S3 <- CreateSeuratObject(counts = o1K_16_S3_data, project = "1K", min.cells = 3)

#Add Time Label

o1K_16_S1$time <- "16.5"
o1K_16_S2$time <- "16.5"
o1K_16_S3$time <- "16.5"


#E17.5
######################################################################################################################################
o1K_17_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_1K_176C_RK/filtered_feature_bc_matrix/")
o1K_17_S1 <- CreateSeuratObject(counts = o1K_17_S1_data, project = "1K", min.cells = 3)

o1K_17_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_1K_477D_LK/filtered_feature_bc_matrix/")
o1K_17_S2 <- CreateSeuratObject(counts = o1K_17_S2_data, project = "1K", min.cells = 3)

o1K_17_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_1K_725D_LK/filtered_feature_bc_matrix/")
o1K_17_S3 <- CreateSeuratObject(counts = o1K_17_S3_data, project = "1K", min.cells = 3)

#Add Time Label
o1K_17_S1$time <- "17.5"
o1K_17_S2$time <- "17.5"
o1K_17_S3$time <- "17.5"




merged <- merge(o1K_15_S2, y =  c(o1K_15_S3, oneK_15_S1, o1K_16_S1, o1K_16_S2, o1K_16_S3, o1K_17_S1, o1K_17_S2, o1K_17_S3))

rm(oneK_15_S1, oneK_15_S1_data, o1K_15_S2, o1K_15_S2_data, o1K_15_S3, o1K_15_S3_data, o1K_16_S1, o1K_16_S1_data, o1K_16_S2, o1K_16_S2_data, o1K_16_S3, o1K_16_S3_data,
   o1K_17_S1, o1K_17_S1_data, o1K_17_S2, o1K_17_S2_data, o1K_17_S3, o1K_17_S3_data)



#QC

merged <- PercentageFeatureSet(merged, pattern = "^Mt-", col.name = "percent.mt")

merged

merged <- subset(merged, subset = nFeature_RNA > 1000 & nFeature_RNA < 2500 & percent.mt < 5)


merged

merged <- JoinLayers(merged)
merged <- NormalizeData(merged)



saveRDS(merged, file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSeurat.rds")



merged.sce <- as.SingleCellExperiment(merged)
rm(merged)

saveRDS(merged.sce, file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSCE.rds")

y = merged.sce$time

#Does changing min expression alter the results

psuper_obj = psupertime(merged.sce, y, sel_genes = "hvg", scale = T, smooth = F, min_expression = 0.01)

saveRDS(psuper_obj, file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/1Kpsuperobj.rds")

#Checking out the results

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Images/trainingresults1K.png")
plot_train_results(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Images/psupertimeorder1K.png")
plot_labels_over_psupertime(psuper_obj, palette = "Set1")
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Images/plotgenes1K.png")
plot_identified_gene_coefficients(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Images/plotgenesgraph1K.png", width = "720")
plot_identified_genes_over_psupertime(psuper_obj, label_name = "Gestational Day", palette = "Set1")
dev.off()

write.csv(psuper_obj$beta_dt, file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/1KpsupertimeResults.csv")




#############################################################################################################################################
#############################################################################################################################################



#2K next

#E15.5 Data, run 1 at a time and save after
######################################################################################################################################
t2K_15_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_2K_188D189G_LK/filtered_feature_bc_matrix/")
t2K_15_S1 <- CreateSeuratObject(counts = t2K_15_S1_data, project = "2K", min.cells = 3)

t2K_15_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_2K_189G187C_RK/filtered_feature_bc_matrix/")
t2K_15_S2<- CreateSeuratObject(counts = t2K_15_S2_data, project = "2K", min.cells = 3)

t2K_15_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/15.5/P1312_E15_2K_232I231D_LK/filtered_feature_bc_matrix/")
t2K_15_S3 <- CreateSeuratObject(counts = t2K_15_S3_data, project = "2K", min.cells = 3)

#Add Time Label

t2K_15_S1$time <- "15.5"
t2K_15_S2$time <- "15.5"
t2K_15_S3$time <- "15.5"

#E16.5
######################################################################################################################################
t2K_16_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_2K_178E_RK/filtered_feature_bc_matrix/")
t2K_16_S1 <- CreateSeuratObject(counts = t2K_16_S1_data, project = "2K", min.cells = 3)


t2K_16_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_2K_200J_LK/filtered_feature_bc_matrix/")
t2K_16_S2 <- CreateSeuratObject(counts = t2K_16_S2_data, project = "2K", min.cells = 3)


t2K_16_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/16.5/P1312_E16_2K_993F_LK/filtered_feature_bc_matrix/")
t2K_16_S3 <- CreateSeuratObject(counts = t2K_16_S3_data, project = "2K", min.cells = 3)

#Add Time Label

t2K_16_S1$time <- "16.5"
t2K_16_S2$time <- "16.5"
t2K_16_S3$time <- "16.5"


#E17.5
######################################################################################################################################
t2K_17_S1_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_2K_176E_LK/filtered_feature_bc_matrix/")
t2K_17_S1 <- CreateSeuratObject(counts = t2K_17_S1_data, project = "2K", min.cells = 3)

t2K_17_S2_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_2K_179A_LK/filtered_feature_bc_matrix/")
t2K_17_S2 <- CreateSeuratObject(counts = t2K_17_S2_data, project = "2K", min.cells = 3)

t2K_17_S3_data <- Read10X(data.dir = "../../nfs_share/Milner/HSRAGD 15.5 17.5/Poster Analysis/17.5/P1312_E17_2K_837A_RK/filtered_feature_bc_matrix/")
t2K_17_S3 <- CreateSeuratObject(counts = t2K_17_S3_data, project = "2K", min.cells = 3)

#Add Time Label
t2K_17_S1$time <- "17.5"
t2K_17_S2$time <- "17.5"
t2K_17_S3$time <- "17.5"




merged <- merge(t2K_15_S2, y =  c(t2K_15_S3, t2K_15_S1, t2K_16_S1, t2K_16_S2, t2K_16_S3, t2K_17_S1, t2K_17_S2, t2K_17_S3))

rm(t2K_15_S1, t2K_15_S1_data, t2K_15_S2, t2K_15_S2_data, t2K_15_S3, t2K_15_S3_data, t2K_16_S1, t2K_16_S1_data, t2K_16_S2, t2K_16_S2_data, t2K_16_S3, t2K_16_S3_data,
   t2K_17_S1, t2K_17_S1_data, t2K_17_S2, t2K_17_S2_data, t2K_17_S3, t2K_17_S3_data)



#QC

merged <- PercentageFeatureSet(merged, pattern = "^Mt-", col.name = "percent.mt")

merged

merged <- subset(merged, subset = nFeature_RNA > 1000 & nFeature_RNA < 2500 & percent.mt < 5)


merged

merged <- JoinLayers(merged)
merged <- NormalizeData(merged)



saveRDS(merged, file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/2KSeurat.rds")



merged.sce <- as.SingleCellExperiment(merged)
rm(merged)

saveRDS(merged.sce, file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/2KSCE.rds")

y = merged.sce$time

#Does changing min expression alter the results

psuper_obj = psupertime(merged.sce, y, sel_genes = "hvg", scale = T, smooth = F, min_expression = 0.01)

saveRDS(psuper_obj, file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/2Kpsuperobj.rds")

#Checking out the results

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Images/trainingresults2K.png")
plot_train_results(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Images/psupertimeorder2K.png")
plot_labels_over_psupertime(psuper_obj, palette = "Set1")
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Images/plotgenes2K.png")
plot_identified_gene_coefficients(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Images/plotgenesgraph2K.png", width = "720")
plot_identified_genes_over_psupertime(psuper_obj, label_name = "Gestational Day", palette = "Set1")
dev.off()

write.csv(psuper_obj$beta_dt, file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/2KpsupertimeResults.csv")


#Try with denoise 1K first

merged.sce <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSCE.rds")
y = merged.sce$time

psuper_obj = psupertime(merged.sce, y, sel_genes = "hvg", scale = T, smooth = T, min_expression = 0.01)

#Still aborted rsession with 2GB unprocessed object
#Denoise prior or reduced size of object further








