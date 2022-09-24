# RPCA (fast integration using reciprocal PCA) after removal of 9 low quality clusters
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(monocle3)
library(harmony)

options(future.globals.maxSize=(90*1024^3)) # set max memory size as 60G
setwd("/data/aronow/Kang/single_cell_projects/10X-Pak/integration//")

# get high quality cells
df_cells = read.table("./ctype_integrated_C2_RPCA_v2.txt", sep = "\t", header=1, row.names = 1)
high_quality_cells = rownames(df_cells[df_cells$low_quality == "False",])
# megre 10x outputs
organoid.list = list()
samples = list.files("../processed data/")
for (i in 1:length(samples))
{
  sample = samples[i]
  # load data
  organoid.counts = Read10X_h5(paste0("../processed data/",sample,"/filtered_feature_bc_matrix.h5"))
  # create seurat objects
  organoid = CreateSeuratObject(counts = organoid.counts)
  # rename cell ids
  organoid = RenameCells(organoid, new.names = paste0(sample,"-",colnames(organoid)))
  organoid$sample = sample
  # calculate mitochondrial gene percentage
  organoid[["percent.mt"]] = PercentageFeatureSet(organoid, pattern = "^MT-")
  # filter cells
  organoid = organoid[,intersect(colnames(organoid), high_quality_cells)]
  # organoid = subset(organoid, subset = nCount_RNA > 1200 & nCount_RNA < 25000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 10)
  # normalization
  organoid = NormalizeData(organoid, verbose = F)
  # find variable genes
  organoid = FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 3000, verbose = F)
  organoid.list[[i]] = organoid
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = organoid.list)
organoid.list <- lapply(X = organoid.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

organoid.anchors <- FindIntegrationAnchors(object.list = organoid.list, anchor.features = features, reduction = "rpca")

organoid.combined <- IntegrateData(anchorset = organoid.anchors)
DefaultAssay(organoid.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
organoid.combined <- ScaleData(organoid.combined, verbose = FALSE)
# run principle component analysis
organoid.combined = RunPCA(organoid.combined, npcs = 30)
# calculate the UMAP coordinates
organoid.combined = RunUMAP(organoid.combined, reduction = "pca", dims = 1:30)
organoid.combined = RunUMAP(organoid.combined, reduction = "pca", dims = 1:30, n.components = 3L)

# find neighbors and clusters (multiple resolutions)
organoid.combined = FindNeighbors(organoid.combined, dim = 1:30)
organoid.combined = FindClusters(organoid.combined, resolution = 1)
organoid.combined = FindClusters(organoid.combined, resolution = 2)
organoid.combined = FindClusters(organoid.combined, resolution = 0.5)
organoid.combined = FindClusters(organoid.combined, resolution = 0.3)

# save results
saveRDS(organoid.combined, "organoid_integrated_highQuality.RDS")
write.table(organoid.combined@meta.data, "ctype_integrated_C2_RPCA_highQuality_v2.txt",sep="\t")
write.table(organoid.combined@reductions$umap@cell.embeddings, "umap_integrated_C2_RPCA_highQuality.txt",sep="\t")
