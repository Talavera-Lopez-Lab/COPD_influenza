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
library(wordcloud)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)
#trace(netClustering, edit = TRUE)

use_python("/Users/cartalop/mambaforge/envs/scanpy/bin", required = TRUE)

### Read in data

cellchat.H_ctrl <- readRDS("../../../data/Epithelial_Healthy-CTRL_anotated.rds")
cellchat.H_ctrl <- updateCellChat(cellchat.H_ctrl)
cellchat.H_ctrl

cellchat.H_iav <- readRDS("../../../data/Epithelial_Healthy-IAV_anotated.rds")
cellchat.H_iav <- updateCellChat(cellchat.H_iav)
cellchat.H_iav

cellchat.C_ctrl <- readRDS("../../../data/Epithelial_COPD-CTRL_anotated.rds")
cellchat.C_ctrl <- updateCellChat(cellchat.C_ctrl)
cellchat.C_ctrl

cellchat.C_iav <- readRDS("../../../data/Epithelial_COPD-IAV_anotated.rds")
cellchat.C_iav <- updateCellChat(cellchat.C_iav)
cellchat.C_iav

### Lift objects to level cell numbers 

group.new = levels(cellchat.H_iav@idents)
cellchat.H_ctrl <- liftCellChat(cellchat.H_ctrl, group.new)

group.new = levels(cellchat.C_iav@idents)
cellchat.C_ctrl <- liftCellChat(cellchat.C_ctrl, group.new)

### Merge objects

object.list <- list(H_CTRL = cellchat.H_ctrl, H_IAV = cellchat.H_iav, C_CTRL = cellchat.C_ctrl, C_IAV = cellchat.C_iav)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat

df.net <- subsetCommunication(cellchat)

unique_to_Hctrl <- setdiff(unique(df.net$H_CTRL$pathway_name), unique(df.net$H_IAV$pathway_name))
unique_to_Hctrl

unique_to_Hiav <- setdiff(unique(df.net$H_IAV$pathway_name), unique(df.net$H_CTRL$pathway_name))
unique_to_Hiav

unique_to_Cctrl <- setdiff(unique(df.net$C_CTRL$pathway_name), unique(df.net$C_IAV$pathway_name))
unique_to_Cctrl

unique_to_Ciav <- setdiff(unique(df.net$C_IAV$pathway_name), unique(df.net$C_CTRL$pathway_name))
unique_to_Ciav

all_unique <- c(unique_to_Hctrl, unique_to_Hiav, unique_to_Cctrl, unique_to_Ciav)
all_unique

### Differential number of interactions or interaction strength among different cell populations

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

### Same visualisation but in heatmap mode

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

### Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

### Visualise which pathways are active in the most significant cell states

####SERPINE+Basal
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SERPINE1+Basal", signaling.exclude = "MIF", comparison = c(1,2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SERPINE1+Basal", signaling.exclude = "MIF", comparison = c(3,4))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SERPINE2+Basal", signaling.exclude = "MIF", comparison = c(1,2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SERPINE2+Basal", signaling.exclude = "MIF", comparison = c(3,4))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "iavAPC_Epi", signaling.exclude = "MIF", comparison = c(1,2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "iavAPC_Epi", signaling.exclude = "MIF", comparison = c(3,4))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MHCII+Club", signaling.exclude = "MIF", comparison = c(1,2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MHCII+Club", signaling.exclude = "MIF", comparison = c(3,4))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "TNC+Basal", signaling.exclude = "MIF", comparison = c(1,2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "TNC+Basal", signaling.exclude = "MIF", comparison = c(3,4))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "FB-like_Basal", signaling.exclude = "MIF", comparison = c(1,2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "FB-like_Basal", signaling.exclude = "MIF", comparison = c(3,4))
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))

###Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional", comparison = c(1,2,3,4))
cellchat <- netEmbedding(cellchat, type = "functional", comparison = c(1,2,3,4))
cellchat <- netClustering(cellchat, type = "functional", comparison = c(1,2,3,4))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

###Identify signaling groups based on their structural similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural", comparison = c(1,2,3,4))
cellchat <- netEmbedding(cellchat, type = "structural", comparison = c(1,2,3,4))
cellchat <- netClustering(cellchat, type = "structural", comparison = c(1,2,3,4))
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")
rankSimilarity(cellchat, type = "structural")

### Compare the overall information flow of each signaling pathway

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2,3,4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2,3,4))
gg1 + gg2

### Compare outgoing (or incoming) signaling associated with each cell population

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'RdPu')
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'RdPu')
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

i = 3
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'RdPu')
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'RdPu')
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'Blues')
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'Blues')
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

i = 3
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'Blues')
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'Blues')
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'OrRd')
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'OrRd')
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

i = 3
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 22, height = 26,   color.heatmap = 'OrRd')
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 22, height = 26,   color.heatmap = 'OrRd')
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

### Identify dysfunctional signaling by comparing the communication probabities

netVisual_bubble(cellchat, 
                 sources.use = c('SERPINE1+Basal', 'SERPINE2+Basal', 'iavAPC_Epi', 'MHCII+Club', 'TNC+Basal', 'FB-like_Basal'), 
                 targets.use = c('SERPINE1+Basal', 'SERPINE2+Basal', 'iavAPC_Epi', 'MHCII+Club', 'TNC+Basal', 'FB-like_Basal'), 
                 comparison = c(1,2,3,4), 
                 angle.x = 45)


### Identify dysfunctional signaling by using differential expression analysis

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "healthy_ctrl"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "group", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.05, 
                                       thresh.p = 1) 

#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.1, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.05, receptor.logFC = -0.05)























