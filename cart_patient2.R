# find all markers distinguishing cluster CR.BA from CR.CD
DefaultAssay(cart.combined) <- "integrated"

cr.basal <- c("BA112", "BA118", "BA154", "BA165")
nr.basal <- c("BA103", "BA108", "BA167")

cr.stim <- c("CD112", "CD118", "CD154", "CD165")
nr.stim <- c("CD103", "CD108", "CD167")

basal <- c("BA112", "BA118", "BA154", "BA165","BA103", "BA108", "BA167")
stim <- c("CD112", "CD118", "CD154", "CD165", "CD103", "CD108", "CD167")

# ----- all basal & all stim
basal.stim.markers <- FindMarkers(cart.combined, ident.1 = basal, ident.2 = stim, min.pct = 0.25)
head(basal.stim.markers, n = 5)
write.table(rownames(basal.stim.markers), file = "txt_files/basal_stim_markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# ----- cr.basal & cr.stim
cr_basal.stim.markers <- FindMarkers(cart.combined, ident.1 = cr.basal, ident.2 = cr.stim, min.pct = 0.25)
head(cr_basal.stim.markers, n = 5)
write.table(rownames(cr_basal.stim.markers), file = "txt_files/cr.basal_cr.stim_markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# ----- nr.basal & nr.stim
nr_basal.stim.markers <- FindMarkers(cart.combined, ident.1 = nr.basal, ident.2 = nr.stim, min.pct = 0.25)
head(nr_basal.stim.markers, n = 5)
write.table(rownames(nr_basal.stim.markers), file = "txt_files/nr.basal_nr.stim_markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# ----- cr.stim & nr.stim
cr_stim.nr_stim.markers <- FindMarkers(cart.combined, ident.1 = cr.stim, ident.2 = nr.stim, min.pct = 0.25)
head(cr_stim.nr_stim.markers, n = 5)
write.table(rownames(cr_stim.nr_stim.markers), file = "txt_files/cr.stim_nr.stim_markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# ----- cr_basal.nr_basal
cr_basal.nr_basal.markers <- FindMarkers(cart.combined, ident.1 = cr.basal, ident.2 = nr.basal, min.pct = 0.25)
head(cr_basal.nr_basal.markers, n = 5)
write.table(rownames(cr_basal.nr_basal.markers), file = "txt_files/cr.basal_nr.basal_markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# ----- ------ VolcanoPlots
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")
BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
library(airway)
library(magrittr)

saveRDS(cart.combined, "cart.combined.rds")

# ---------- basal.stim
png("cart_images/volcanoplots/basal_stim_volcanoplot.png", type = "cairo", width = 900, height = 600)
EnhancedVolcano(basal.stim.markers , 
                rownames(basal.stim.markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
dev.off()

# ---------- cr_basal.stim
png("cart_images/volcanoplots/cr_basal_stim_volcanoplot.png", type = "cairo", width = 900, height = 600)
EnhancedVolcano(cr_basal.stim.markers , 
                rownames(cr_basal.stim.markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
dev.off()

# ---------- nr_basal.stim
png("cart_images/volcanoplots/nr_basal_stim_volcanoplot.png", type = "cairo", width = 900, height = 600)
EnhancedVolcano(nr_basal.stim.markers , 
                rownames(nr_basal.stim.markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
dev.off()

# ---------- cr_stim.nr_stim
png("cart_images/volcanoplots/cr_stim_nr_stim_volcanoplot.png", type = "cairo", width = 900, height = 600)
EnhancedVolcano(cr_stim.nr_stim.markers , 
                rownames(cr_stim.nr_stim.markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
dev.off()

# ---------- cr_basal.nr_basal
png("cart_images/volcanoplots/cr_basal_nr_basal_volcanoplot2.png", type = "cairo", width = 900, height = 600)
EnhancedVolcano(cr_basal.nr_basal.markers, 
                rownames(cr_basal.nr_basal.markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
dev.off()

# Check Cells per Cluster UMAP
n_cells <- FetchData(cart.combined, 
                     vars = c("ident", "seurat_clusters")) %>%
  dplyr::count(ident, seurat_clusters) %>%
  tidyr::spread(ident, n)

View(n_cells)

# FindAllMarkers
Idents(cart.combined) <- "seurat_clusters"

cart.markers <- FindAllMarkers(cart.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cart.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# find all markers distinguishing cluster 1 from cluster 1
cluster0.markers <- FindMarkers(cart.combined, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
head(cluster0.markers, n = 5)
write.table(rownames(cluster0.markers), file = "txt_files/cluster0_cluster1.markers.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)