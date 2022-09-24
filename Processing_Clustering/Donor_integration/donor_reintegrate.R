library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(monocle3)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggridges)
library(cowplot)

# find trajectory of d101 only organoids
organoids = readRDS("../../v2_filtered/combined_all_organoids_filtered_v3.rds")
donor = organoids[, organoids$source == 'Donor']

# find lineage on d101
reintegration = function(organoids_dx, split.by = "sample")
{
  DefaultAssay(organoids_dx) = "RNA"
  
  organoids_dx <- SplitObject(organoids_dx, split.by = split.by)
  for (i in seq_along(organoids_dx)) {
    organoids_dx[[i]] <- NormalizeData(organoids_dx[[i]]) %>% FindVariableFeatures()
  }
  features <- SelectIntegrationFeatures(organoids_dx)
  for (i in seq_along(along.with = organoids_dx)) {
    organoids_dx[[i]] <- ScaleData(organoids_dx[[i]], features = features) %>% RunPCA(features = features)
  }
  anchors <- FindIntegrationAnchors(organoids_dx, anchor.features = features, reduction = "rpca", dims = 1:30)
  integrated <- IntegrateData(anchors, dims = 1:30)
  
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated)
  integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
  integrated <- FindNeighbors(integrated, dims = 1:30)
  # integrated <- FindClusters(integrated)
  return (integrated)
}

donor_reintegrated = reintegration(donor)
saveRDS(donor_reintegrated, 'donor_reintegrated.rds')
write.table(donor_reintegrated@reductions$UMAP@cell.embeddings, 'umap_reintegrated.txt', sep = '\t')

DimPlot(donor_reintegrated, group.by = 'cellclass_draft_v2', label = T) + theme(legend.position = 'none')
DimPlot(donor_reintegrated, group.by = 'timepoint_v2') + theme(legend.position = c(0.8,0.92))

# re-annotate original clusters
original_cluster_cols = colnames(donor_reintegrated@meta.data)[grepl('snn', colnames(donor_reintegrated@meta.data))]
for (col in original_cluster_cols)
{
  names(donor_reintegrated@meta.data)[names(donor_reintegrated@meta.data) == col] = paste0(col, '-ori')
}

# do clustering
res_vec = c(0.3, 0.5, 1, 2)
for (res in res_vec)
{
  donor_reintegrated = FindClusters(donor_reintegrated, resolution = res)
}

# transfer data 
DefaultAssay(donor_reintegrated) = 'RNA'
SaveH5Seurat(donor_reintegrated, filename = 'donor_reintegrated.h5Seurat')
Convert('donor_reintegrated.h5Seurat', dest = 'h5ad')

