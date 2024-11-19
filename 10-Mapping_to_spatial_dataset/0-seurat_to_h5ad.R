# Using seruat_to_h5ad-env environment
library(sceasy)
library(reticulate)
seurat_object <- readRDS(".data/GSE135893_ILD_annotated_fullsize.rds")
seurat_object@images <- list()
sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                       outFile='.data/GSE135893_ILD_annotated_fullsize.h5ad')