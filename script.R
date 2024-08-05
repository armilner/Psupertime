library(Seurat)
library(SingleCellExperiment)
library(psupertime)
library(scran)
library(scater)
library(SAVER)
library(Rmagic)
library(reticulate)

use_python("/home/amilner2/Desktop/amilner2-Data/anaconda3/bin/python", required = TRUE)




#1K denoise SCE prior to psupertime with scran?

sce <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSCE.rds")



#Rmagic???

expr_matrix <- logcounts(sce)

# Convert the dgCMatrix to a Python scipy sparse matrix and save it as an npz file
expr_matrix_py <- r_to_py(expr_matrix)

# Save the matrix as an npz file

py_run_string("
+ import scipy.sparse
+ import numpy as np
+ 
+ # Convert the R sparse matrix to a scipy sparse matrix
+ expr_matrix = r.expr_matrix_py
+ 
+ # Save as an npz file
+ scipy.sparse.save_npz('../../nfs_share/Milner/sparse_matrix.npz', expr_matrix)
+ ")


#Loading it back in

magic_results <- py_run_string("
import numpy as np
npzfile = np.load('../../nfs_share/Milner/logcounts_matrix_magic.npz')
logcounts_matrix_magic = npzfile['arr_0']
")

logcounts_matrix_magic <- py$logcounts_matrix_magic

saveRDS(logcounts_matrix_magic, file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/magiclogcount1K.rds")

# Assign the MAGIC results back to the logcounts of your SCE object
rownames(matrix) <- rownames(sce)
colnames(matrix) <- colnames(sce)

logcounts(sce) <- matrix

# Proceed with further analysis in R

y = sce$time

#Does changing min expression alter the results

#identifying highly variable genes

sce <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/1K/Objects/1KSCE.rds")

# Identify highly variable genes using scran
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 2000)
hvg <- as.character(hvg)


psuper_obj = psupertime(sce, y, sel_genes = "list", gene_list = hvg, scale = F, smooth = F, min_expression = 0.01)

saveRDS(psuper_obj, file = "../../nfs_share/Milner/Psupertime/With denoising/1K/objects/psuperobj1Kdenoise.rds")

#Checking out the results

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/1K/images/trainingresults1K.png")
plot_train_results(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/1K/images/psupertimeorder1KNogeomvline.png")
plot_labels_over_psupertime(psuper_obj, palette = "Set1")
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/1K/images/plotgenes1K.png")
plot_identified_gene_coefficients(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/1K/images/plotgenesgraph1K.png", width = "720")
plot_identified_genes_over_psupertime(psuper_obj, label_name = "Gestational Day", palette = "Set1")
dev.off()

write.csv(psuper_obj$beta_dt, file = "../../nfs_share/Milner/Psupertime/With denoising/1K/images/1KpsupertimeResults.csv")





