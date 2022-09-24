# RPCA (fast integration using reciprocal PCA)
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(SeuratWrappers)
library(monocle3)

options(future.globals.maxSize=(90*1024^3)) # set max memory size as 60G

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
  print(sample)
  print(dim(organoid))
  organoid = subset(organoid, subset = nCount_RNA > 1200 & nCount_RNA < 25000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 10)
  print(dim(organoid))
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

# find neighbors and clusters
organoid.combined = FindNeighbors(organoid.combined, dim = 1:30)
organoid.combined = FindClusters(organoid.combined, resolution = 1)
organoid.combined = FindClusters(organoid.combined, resolution = 2)
organoid.combined = FindClusters(organoid.combined, resolution = 0.5)
organoid.combined = FindClusters(organoid.combined, resolution = 0.3)

saveRDS(organoid.combined, "organoid_integrated_C2_RPCA.RDS")
write.table(organoid.combined@meta.data, "ctype_integrated_C2_RPCA.txt",sep="\t")
write.table(organoid.combined@reductions$umap@cell.embeddings, "umap_integrated_C2_RPCA.txt",sep="\t")

# regress out ncount and nfeature effects
organoid.combined = readRDS("./organoid_integrated_C2_RPCA.RDS")
organoid.regressed <- ScaleData(organoid.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA"))
organoid.regressed = RunPCA(organoid.regressed, npcs = 30)
organoid.regressed = RunUMAP(organoid.regressed, reduction = "pca", dims = 1:30)
organoid.regressed = FindNeighbors(organoid.regressed, dim = 1:30)
organoid.regressed = FindClusters(organoid.regressed, resolution = 1)
organoid.regressed = FindClusters(organoid.regressed, resolution = 2)
organoid.regressed = FindClusters(organoid.regressed, resolution = 0.5)
organoid.regressed = FindClusters(organoid.regressed, resolution = 0.3)
saveRDS(organoid.regressed, "organoid_integrated_C2_RPCA_regressed.RDS")

# read previous annotations
df_cells = read.table("ctype_integrated_v2.txt", sep = "\t", header = 1, row.names = 1)
df_cells = df_cells[colnames(organoid.combined),]
for (col in colnames(df_cells))
{
  if (! col %in% colnames(organoid.combined@meta.data))
  {
    organoid.combined@meta.data[[col]] = df_cells[[col]]
  }
}
df_cellxgene = read.table("cellxgene_0126.csv", sep = ",", header = 1, row.names = 1)
organoid.combined$Cell.class_v2 = df_cellxgene[colnames(organoid.combined),"Cell.class_v2"]

# label clusters as low quality cells
organoid.combined$low_quality = "False"
organoid.combined@meta.data[organoid.combined$integrated_snn_res.2 %in% c(1,5,8,9,10,11,12,24,25),"low_quality"] = "True"
write.table(organoid.combined@meta.data, "ctype_integrated_C2_RPCA_v2.txt",sep="\t")

# pseudotime trajectory
cds = as.cell_data_set(organoid.combined)
cds = cluster_cells(cds)
cds = learn_graph(cds)
plot_cells(cds)

# plots
DimPlot(organoid.combined, group.by = "Cell.class_v1")
DimPlot(organoid.combined, group.by = "Cell.class_v2")
DimPlot(organoid.combined, group.by = "integrated_snn_res.0.5", label = T)

VlnPlot(organoid.combined, group.by = "integrated_snn_res.2", features = c("nCount_RNA"), pt.size = 0, y.max = 18000) + theme(legend.position = "none")
VlnPlot(organoid.combined, group.by = "integrated_snn_res.2", features = c("nFeature_RNA"), pt.size = 0, y.max = 5000) + theme(legend.position = "none")
VlnPlot(organoid.combined, group.by = "integrated_snn_res.2", features = c("percent.mt"), pt.size = 0, y.max = 10) + theme(legend.position = "none")

