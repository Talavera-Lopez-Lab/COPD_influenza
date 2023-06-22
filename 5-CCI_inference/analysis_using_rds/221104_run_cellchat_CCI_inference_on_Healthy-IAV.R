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

### Read in data

cellchat <- readRDS("../../data/Healthy-IAV_anotated.rds")
cellchat

### Set up ligand-receptor interaction database for `cellchat`

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

### Check 

df.net <- subsetCommunication(cellchat)
unique(df.net$pathway_name)

### Visualisations

options(repr.plot.width = 10, repr.plot.height = 15)
pathways.show <- c("CXCL")

netAnalysis_contribution(cellchat, signaling = pathways.show)
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

options(repr.plot.width = 10, repr.plot.height = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("FN1", "IL1"))
gg1 + gg2

options(repr.plot.width = 5, repr.plot.height = 5)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 17, height = 19, color.heatmap = "RdPu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 17, height = 19, color.heatmap = "RdPu")
ht1 + ht2

### Identify global communication patterns

selectK(cellchat, pattern = "outgoing")

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,  width = 10, height = 10)

### Sankey plot

options(repr.plot.width = 40, repr.plot.height = 22.5)
netAnalysis_river(cellchat, pattern = "outgoing", font.size = 2.5, font.size.title = 20)

netAnalysis_dot(cellchat, pattern = "outgoing")

pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathway[1,] # show one ligand-receptor pair

### Identify global communication patterns

selectK(cellchat, pattern = "incoming")

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,  width = 10, height = 10)

### Sankey plot

options(repr.plot.width = 40, repr.plot.height = 22.5)
netAnalysis_river(cellchat, pattern = "incoming", font.size = 2.5, font.size.title = 20)

netAnalysis_dot(cellchat, pattern = "incoming")

pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathway[1,] # show one ligand-receptor pair


# Hierarchy plot

vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

plotGeneExpression(cellchat, signaling = "FN1")

# Save object 

saveRDS(cellchat, file = "/Users/carlos.lopez/Downloads/cellchat_objects/Healthy-IAV_anotated.rds")
