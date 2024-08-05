library(psupertime)
library(SingleCellExperiment)

# load the data
data(acinar_hvg_sce)

# run psupertime
y           = acinar_hvg_sce$donor_age
psuper_obj  = psupertime(acinar_hvg_sce, y, sel_genes='all')


psuper_obj

g = plot_train_results(psuper_obj)

g

g2 = plot_labels_over_psupertime(psuper_obj, label_name='Donor age')
g2

g3 = plot_identified_gene_coefficients(psuper_obj)
g3


g4 = plot_identified_genes_over_psupertime(psuper_obj, label_name='Donor age')
g4

g5 = plot_predictions_against_classes(psuper_obj)
g5
