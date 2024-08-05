library(SingleCellExperiment)
library(psupertime)

library(scran)
library(Rmagic)
library(reticulate)

use_python("/home/amilner2/Desktop/amilner2-Data/anaconda3/bin/python", required = TRUE)
#2K denoise SCE prior to psupertime with scran?

sce <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/2KSCE.rds")

#Rmagic???

expr_matrix <- logcounts(sce)

# Convert the dgCMatrix to a Python scipy sparse matrix and save it as an npz file
expr_matrix_py <- r_to_py(expr_matrix)

# Save the matrix as an npz file

py_run_string("
import scipy.sparse
import numpy as np

# Convert the R sparse matrix to a scipy sparse matrix
expr_matrix = r.expr_matrix_py

# Save as an npz file
scipy.sparse.save_npz('../../nfs_share/Milner/2Ksparse_matrix.npz', expr_matrix)
")

#Loading it back in

magic_results <- py_run_string("
import numpy as np
npzfile = np.load('../../nfs_share/Milner/2Klogcounts_matrix_magic.npz')
logcounts_matrix_magic = npzfile['arr_0']
")

logcounts_matrix_magic <- py$logcounts_matrix_magic

saveRDS(logcounts_matrix_magic, file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/magiclogcount2K.rds")

# Identify highly variable genes using scran
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 2000)
hvg <- as.character(hvg)

# Assign the MAGIC results back to the logcounts of your SCE object
rownames(logcounts_matrix_magic) <- rownames(sce)
colnames(logcounts_matrix_magic) <- colnames(sce)

logcounts(sce) <- logcounts_matrix_magic

# Proceed with further analysis in R

y = sce$time

#Does changing min expression alter the results

#identifying highly variable genes

sce <- readRDS(file = "../../nfs_share/Milner/Psupertime/Both consistent/2K/Objects/2KSCE.rds")



psuper_obj = psupertime(sce, y, sel_genes = "list", gene_list = hvg, scale = F, smooth = F, min_expression = 0.01)

saveRDS(psuper_obj, file = "../../nfs_share/Milner/Psupertime/With denoising/2K/objects/psuperobj2Kdenoise.rds")

#Checking out the results

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/2K/images/trainingresults2K.png")
plot_train_results(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/2K/images/psupertimeorder2KNogeomvline.png")
plot_labels_over_psupertime(psuper_obj, palette = "Set1")
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/2K/images/plotgenes2K.png")
plot_identified_gene_coefficients(psuper_obj)
dev.off()

png(filename = "../../nfs_share/Milner/Psupertime/With denoising/2K/images/plotgenesgraph2K.png", width = "720")
plot_identified_genes_over_psupertime(psuper_obj, label_name = "Gestational Day", palette = "Set1")
dev.off()

write.csv(psuper_obj$beta_dt, file = "../../nfs_share/Milner/Psupertime/With denoising/2K/images/2KpsupertimeResults.csv")


