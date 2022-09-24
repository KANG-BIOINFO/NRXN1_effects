library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggridges)
library(cowplot)

source('../Figures/functions/_global.R')

samples = c('Cre-d160', 'Flp-d160')

organoid.list = list()
for (i in 1:2)
{
  sample = samples[i]
  organoid.counts = Read10X_h5(paste0("/data/aronow/Kang/single_cell_projects/10X-Pak/processed data/",sample,"/filtered_feature_bc_matrix.h5"))
  organoid = CreateSeuratObject(counts = organoid.counts)
  organoid = RenameCells(organoid, new.names = paste0(sample,"-",colnames(organoid)))
  organoid$sample = sample
  organoid[["percent.mt"]] = PercentageFeatureSet(organoid, pattern = "^MT-")
  organoid.list[[i]] = organoid
}

# find lineage on d101
run_integration = function(organoids_dx)
{
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
  integrated <- FindClusters(integrated, resolution = 0.5)
  integrated <- FindClusters(integrated, resolution = 1.0)
  integrated <- FindClusters(integrated, resolution = 2.0)
  return (integrated)
}

d160 = run_integration(organoid.list)
saveRDS(d160, 'Engineered_d160_integrated.rds')

DimPlot(d160, group.by = 'integrated_snn_res.0.5', label = T, label.size = 6) + theme(legend.position = 'none')
DimPlot(d160, group.by = 'sample', cols = c('#EC706B', '#A4CBF1')) 
DefaultAssay(d160) = 'RNA'
FeaturePlot(d160, features = genes_atlas, ncol = 6)
Idents(d160) = d160$sample
VlnPlot(d160, features = c('NRXN1'))
VlnPlot(d160, features = genes_in)
VlnPlot(d160, features = genes_cn)

