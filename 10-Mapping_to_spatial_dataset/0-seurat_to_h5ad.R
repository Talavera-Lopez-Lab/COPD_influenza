library(sceasy)
library(reticulate)
seurat_object <- readRDS(".data/GSE135893_ILD_annotated_fullsize.rds")
sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                       outFile='GSE135893_ILD_annotated_fullsize')