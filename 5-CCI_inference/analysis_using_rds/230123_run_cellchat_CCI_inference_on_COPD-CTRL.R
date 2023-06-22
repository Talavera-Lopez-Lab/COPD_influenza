### Load required modules

library(NMF)
library(dplyr)
library(igraph)
library(Matrix)
library(ggplot2)
library(CellChat) 
library(patchwork)
library(ggalluvial)
library(reticulate)
options(stringsAsFactors = FALSE)

### Read in data

ad <- import("anndata", convert = FALSE)
pd <- import("pandas", convert = FALSE)
ad_object <- ad$read_h5ad("../data/COPD-CTRL_anotated.")

### Access expression matrix

data.input <- t(py_to_r(ad_object$X))

### Add metadata

rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))

meta.data <- py_to_r(ad_object$obs)
meta <- meta.data

### Create `cellchat` object

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "leiden_states")

### Set up ligand-receptor interaction database for `cellchat`

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

### Process expression data

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)

### Export results as dataframe

df.net <- subsetCommunication(cellchat)
head(df.net)
write.table(df.net, sep = ',', row.names = FALSE, './inferences/COPD-CTRL_cellchat_net.csv')

### Infer cell-cell communication

cellchat <- computeCommunProbPathway(cellchat)

### Calculate aggregated cell-cell communication

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,3), xpd = TRUE)
options(repr.plot.width = 40, repr.plot.height = 40)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

unique(df.net$pathway_name)

options(repr.plot.width = 10, repr.plot.height = 15)
pathways.show <- c("FN1")
#pathways.show <- c("CXCL")
netAnalysis_contribution(cellchat, signaling = pathways.show)
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

n   