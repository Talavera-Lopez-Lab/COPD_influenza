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
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)
trace(netClustering, edit=TRUE)

use_python("/Users/cartalop/mambaforge/envs/scanpy/bin", required = TRUE)

### Read in data

cellchat.ctrl <- readRDS("../../data/Healthy-CTRL_annotated.rds")
cellchat.ctrl <- updateCellChat(cellchat.ctrl)
cellchat.ctrl

cellchat.iav <- readRDS("../../data/Healthy-IAV_annotated.rds")
cellchat.iav <- updateCellChat(cellchat.iav)
cellchat.iav

### Merge objects

object.list <- list(CTRL = cellchat.ctrl, IAV = cellchat.iav)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

##Predict general principles of cell-cell communication

###Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

###Differential number of interactions or interaction strength among different cell populations

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

###Differential number of interactions or interaction strength among different cell types

group.cellType <- c(rep("Mixed_11", 3), rep("Mixed_12", 3), rep("Mixed_13", 3), rep("Mixed_16", 3), rep("Mixed_17", 3), rep("Basal_3", 3))
group.cellType <- factor(group.cellType, levels = c("Mixed_11", "Mixed_12", "Mixed_13", 'Mixed_16', 'Mixed_17', 'Basal_3'))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

###Show the differential number of interactions or interaction strength between any two cell types using circle plot. Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

###Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

###Explore pathways of major souces and targets for CCIs

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mixed_13", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mixed_12", signaling.exclude = c("MIF"))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mixed_11", signaling.exclude = "MIF")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mixed_16", signaling.exclude = c("MIF"))
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mixed_17", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4,gg5))
                                                                                                                                    
###Part II: Identify the conserved and context-specific signaling pathways
###Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

###Identify signaling groups based on structure similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

###Compute and visualize the pathway distance in the learned joint manifold

rankSimilarity(cellchat, type = "functional")
rankSimilarity(cellchat, type = "structural")

###Compare the overall information flow of each signaling pathway

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

###Compare outgoing (or incoming) signaling associated with each cell population

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 16, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 16, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

###Identify the upgulated and down-regulated signaling ligand-receptor pairs

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Healthy", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Healthy", angle.x = 45, remove.isolate = T)
gg1 + gg2

###Identify dysfunctional signaling by using differential expression analysis

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "IAV"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "IAV",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(3:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(3:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(2:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(2:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

computeEnrichmentScore(net.down, species = 'human')
computeEnrichmentScore(net.up, species = 'human')

###Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("COL9A2") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("MK") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(18, 14, 17), lab.cex = 0.5, title.name = paste0("Signaling sending from Sig - ", names(object.list)[i]))
}

###Compare all the interactions sending from fibroblast to inflamatory immune cells

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(4, 16), targets.use = c(17, 14, 18), lab.cex = 0.5,  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}

###Show all the significant signaling pathways from fibroblast to immune cells

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(17,14,16,18), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}

###Compare the signaling gene expression distribution between different datasets

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "IAV")) # set factor level
plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", colors.ggplot = T)


# Save object 

saveRDS(cellchat, file = "/Users/carlos.lopez/Downloads/cellchat_objects/Healthy-CTRL_anotated.rds")
