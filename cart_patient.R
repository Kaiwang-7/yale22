library(dplyr)
library(Seurat)
library(patchwork)

library(ggplot2)
library(sctransform)

library(cowplot)

setwd("/gpfs/ycga/project/fan/kw824")

# ----- Add sample name to the colnames for patient 103
NR_BA103 <- RenameCells(object = NR_BA103, add.cell.id = "NR_BA103")
NR_CD103 <- RenameCells(object = NR_CD103, add.cell.id = "NR_CD103")

saveRDS(NR_BA103, "NR_BA103.rds")
saveRDS(NR_CD103, "NR_CD103.rds")

file_paths <- list.files("Patient_Data")
file_names <- gsub(file_paths, pattern = "\\.rds$", replacement = "")
cart.list <- lapply(paste0("Patient_Data/", file_paths), readRDS)
names(cart.list) <- file_names

for(patient_name in file_names) {
  DefaultAssay(cart.list[[patient_name]]) <- "RNA"
}

saveRDS(cart.list, "Patient_ALL.rds")

cart.list <- readRDS("Patient_ALL.rds")

# store mitochondrial percentage in object meta data
cart.list <- lapply(cart.list, function(X) {X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-"); X} )

#SCTransform
cart.list <- lapply(X = cart.list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)

#PCA Integration
features <- SelectIntegrationFeatures(object.list = cart.list, nfeatures = 3000)
cart.list <- lapply(X = cart.list, FUN = RunPCA, features = features)
cart.list <- PrepSCTIntegration(object.list = cart.list, anchor.features = features)


saveRDS(cart.list, "cart_integration.rds")

cart.anchors <- FindIntegrationAnchors(object.list = cart.list, normalization.method = "SCT",
                                       anchor.features = features, dims = 1:50, reduction = "rpca")

saveRDS(cart.anchors, "cart_anchors.rds")

cart.combined <- IntegrateData(anchorset = cart.anchors, normalization.method = "SCT", dims = 1:50)

saveRDS(cart.combined, "cart_combined.rds")

#QC
png("cart_images/VlnPlot.png", type = "cairo", width = 900, height = 600)
levels(cart.combined@active.ident) <- restrp("anything",44)
VlnPlot(cart.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

cart.combined <- subset(cart.combined, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)

# Visualization and Clustering Steps
cart.combined <- RunPCA(cart.combined, verbose = FALSE, approx = FALSE)
cart.combined <- FindNeighbors(cart.combined, dims = 1:30)
cart.combined <- FindClusters(cart.combined, resolution = 0.35, verbose = FALSE)
cart.combined <- RunUMAP(cart.combined, dims = 1:30)

png("cart_images/Patient_DimPlot_Rerun.png", type = "cairo", width = 900, height = 600)
DimPlot(cart.combined, label = TRUE)
dev.off()

# Normalize ADT
cart.combined <- NormalizeData(cart.combined, normalization.method = "CLR", margin = 2)

# Plot ADT Markers
DefaultAssay(cart.combined) <- "ADT"

ADT <- c("ADT-CD4", "ADT-CD45RA", "ADT-CD45RO", "ADT-CD62L", "ADT-CD95", "ADT-CD8", "ADT-CD127", "ADT-CCR7", "ADT-CD27", "ADT-HLA-DR", "ADT-CD69", "ADT-CD28", "ADT-PD-1", "ADT-Tim-3", "ADT-LAG-3", "ADT-CTLA-4", "ADT-TIGIT")

png("cart_images/Patient_FeaturePlot.png", type = "cairo", width = 900, height = 600)
FeaturePlot(cart.combined, features = c("ADT"))
dev.off()

cart.combined$samples <- NA
cart.combined$samples[grepl('BA112', colnames(cart.combined))] <- 'BA112'
cart.combined$samples[grepl('BA118', colnames(cart.combined))] <- 'BA118'
cart.combined$samples[grepl('BA154', colnames(cart.combined))] <- 'BA154'
cart.combined$samples[grepl('BA165', colnames(cart.combined))] <- 'BA165'

cart.combined$samples[grepl('CD112', colnames(cart.combined))] <- 'CD112'
cart.combined$samples[grepl('CD118', colnames(cart.combined))] <- 'CD118'
cart.combined$samples[grepl('CD154', colnames(cart.combined))] <- 'CD154'
cart.combined$samples[grepl('CD165', colnames(cart.combined))] <- 'CD165'

cart.combined$samples[grepl('BA103', colnames(cart.combined))] <- 'BA103'
cart.combined$samples[grepl('BA108', colnames(cart.combined))] <- 'BA108'
cart.combined$samples[grepl('BA167', colnames(cart.combined))] <- 'BA167'

cart.combined$samples[grepl('CD103', colnames(cart.combined))] <- 'CD103'
cart.combined$samples[grepl('CD108', colnames(cart.combined))] <- 'CD108'
cart.combined$samples[grepl('CD167', colnames(cart.combined))] <- 'CD167'

# Create subset of patients
CR.BA <-cart.combined[,grepl('BA112', colnames(cart.combined)) | grepl('BA118', colnames(cart.combined)) | grepl('BA154', colnames(cart.combined)) | grepl('BA165', colnames(cart.combined))]
CR.CD <-cart.combined[,grepl('CD112', colnames(cart.combined)) | grepl('CD118', colnames(cart.combined)) | grepl('CD154', colnames(cart.combined)) | grepl('CD165', colnames(cart.combined))]
NR.BA <-cart.combined[,grepl('BA103', colnames(cart.combined)) | grepl('BA108', colnames(cart.combined)) | grepl('BA167', colnames(cart.combined))]
NR.CD <-cart.combined[,grepl('CD103', colnames(cart.combined)) | grepl('CD108', colnames(cart.combined)) | grepl('CD167', colnames(cart.combined))]

# ----- CD and BA
CD.all <- cart.combined[,grepl('CD112', colnames(cart.combined)) | grepl('CD118', colnames(cart.combined)) | grepl('CD154', colnames(cart.combined)) | grepl('CD165', colnames(cart.combined)) | grepl('CD103', colnames(cart.combined)) | grepl('CD108', colnames(cart.combined)) | grepl('CD167', colnames(cart.combined))]
BA.all <- cart.combined[,grepl('BA112', colnames(cart.combined)) | grepl('BA118', colnames(cart.combined)) | grepl('BA154', colnames(cart.combined)) | grepl('BA165', colnames(cart.combined)) | grepl('BA103', colnames(cart.combined)) | grepl('BA108', colnames(cart.combined)) | grepl('BA167', colnames(cart.combined))]

# ClusterPlot by samples
png("cart_images/Patient_DimPlot3.png", type = "cairo", width = 900, height = 600)
DimPlot(cart.combined, label = TRUE, group.by = "samples")
dev.off()

# ClusterPlot
png("cart_images/Patient_DimPlot4.png", type = "cairo", width = 900, height = 600)
DimPlot(CR.BA, label = TRUE)
dev.off()

png("cart_images/Patient_DimPlot5.png", type = "cairo", width = 900, height = 600)
DimPlot(CR.CD, label = TRUE)
dev.off()

# VlnPlot by Group
#features <- c("GZMB", "IFNG", "CCL3", "PRF1", "TNF", "IL4", "IL10", "IL13",
#              "IL22", "CDF2", "IL2", "IL5", "IL9", "IL21", "IL8", "IL17A",
#              "IL17F", "CCL4", "CCL5", "CXCL10", "CCL20")
# ----- Not Found: TGFB1, IL7, IL15
features <- c("GZMB", "IFNG", "CCL3", "PRF1", "TNF", "IL13",
              "IL2", "IL21", "CCL4", "CCL5", "CCL20",
              "CSF2", "XCL1", "XCL2", "CXCL8")

#, group.by = "samples",
png("cart_images/VlnPlot3.png", type = "cairo", width = 900, height = 600)
VlnPlot(cart.combined, features = features, stack = TRUE, flip = TRUE, fill.by = "ident",
        idents = c("BA112", "BA118", "BA154", "BA165", "CD112", "CD118", "CD154", "CD165", "BA103", "BA108", "BA167", "CD103", "CD108", "CD167"))
dev.off()

# FindAllMarkers
CR.BA.markers <- FindAllMarkers(CR.BA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CR.BA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

CR.CD.markers <- FindAllMarkers(CR.CD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CR.CD.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

NR.BA.markers <- FindAllMarkers(NR.BA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NR.BA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

NR.CD.markers <- FindAllMarkers(NR.CD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NR.CD.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# ----- Normalize
DefaultAssay(CR.BA) <- "ADT"
DefaultAssay(CR.CD) <- "ADT"
DefaultAssay(NR.BA) <- "ADT"
DefaultAssay(NR.CD) <- "ADT"

CR.BA <- NormalizeData(CR.BA, normalization.method = "CLR", margin = 2)
CR.CD <- NormalizeData(CR.CD, normalization.method = "CLR", margin = 2)
NR.BA <- NormalizeData(NR.BA, normalization.method = "CLR", margin = 2)
NR.CD <- NormalizeData(NR.CD, normalization.method = "CLR", margin = 2)

# ----- # VlnPlot find markers of cluster CR.BA
DefaultAssay(CR.BA) <- "ADT"
png("cart_images/VlnPlotC1.png", type = "cairo", width = 900, height = 600)
VlnPlot(CR.BA, features = c(row.names(CR.BA.markers)[1], row.names(CR.BA.markers)[2]))
dev.off()

# ----- ADT RidgePlots
png("cart_images/RidgePlot.png", type = "cairo", width = 900, height = 600)
RidgePlot(CR.BA, features = "ADT-LAG-3", ncol = 2)
dev.off()

# ----- ADT FeaturePlots
png("cart_images/test.png", type = "cairo", width = 900, height = 600)
FeaturePlot(object = CR.BA, 
            features = c("ADT-LAG-3", "ADT-TIGIT", "ADT-CD45RO", "ADT-LAG-31", "ADT-CD28", "ADT-CD62L", "ADT-CD95", "ADT-CD45RA", "ADT-Tim-3"),
            cols = c("grey", "blue"))
dev.off()

# ----- StackedBarPlot to show cell proportion within clusters

ggplot(df, aes(x = x, y = y, fill = group)) + 
  geom_bar(stat = "identity") 


# ----- Cluster Barplot
library(RColorBrewer)

pt <- table(Idents(CD.all), CD.all$seurat_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)


png("cart_images/barplots/CD_all_cluster_samples1.png", type = "cairo", width = 900, height = 600)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

# ----- Add metadata stim & basal
#cart.combined$data[cart.combined@meta.data$samples %in% c("CD112", "CD118", "CD154", "CD165", "CD103", "CD108", "CD167")] <- "stim"
#cart.combined$data["BA112", "BA118", "BA154", "BA165","BA103", "BA108", "BA167"] <- basal
cart.combined$stim <- NA
cart.combined$stim[grepl('CD', colnames(cart.combined))] <- 'stim'
cart.combined$stim[grepl('BA', colnames(cart.combined))] <- 'basal'

# ------ newbarplot
Idents(cart.combined) <- "stim"

df <- table(Idents(cart.combined), cart.combined$seurat_clusters)
df <- as.data.frame(df)
df$Var1 <- as.character(df$Var1)

png("cart_images/barplots/grouped/stim_basal.png", type = "cairo", width = 900, height = 600)
ggplot(df, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Clusters") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

cart.markers <- list()
for (clusterind in 0:10) {
  cart.markers[[as.character(clusterind)]] <- FindMarkers(cart.combined, ident.1 = clusterind, assay = "integrated", only.pos = TRUE)
}

for (clusterind in 0:10) {
  write.table(rownames(cart.combined[[as.character(clusterind)]])
              quote = 