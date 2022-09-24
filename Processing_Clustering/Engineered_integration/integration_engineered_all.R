library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(plyr)

original_engineered_samples = c("Cre-P61-d112-Live", "Cre-P61-d50-Live", "Cre-P61-d50-Methanol",
                                "Cre-P63-d23-Live", "Flp-P61-d112-Live", "Flp-P61-d50-Live",
                                "Flp-P61-d50-Methanol", "Flp-P63-d23-Live")
engineered_ipsc_samples = c(original_engineered_samples, c('Cre-d160', 'Flp-d160'))

df_engineered_selected = read.table("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/ctype_integrated_C2_RPCA_highQuality_v5.txt",
                                    sep = "\t", header = 1, row.names = 1)
original_selected_cells = rownames(df_engineered_selected)

seurat_integration = function(samples)
{
  organoid.list = list()
  # 1. integration and preprocessing
  for (i in 1:length(samples))
  {
    sample = samples[i]
    organoid.counts = Read10X_h5(paste0("/data/aronow/Kang/single_cell_projects/10X-Pak/processed data/",sample,"/filtered_feature_bc_matrix.h5"))
    organoid = CreateSeuratObject(counts = organoid.counts)
    organoid = RenameCells(organoid, new.names = paste0(sample,"-",colnames(organoid)))
    organoid$sample = sample
    organoid[["percent.mt"]] = PercentageFeatureSet(organoid, pattern = "^MT-")

    # QC
    organoid = subset(organoid, subset = nCount_RNA > 1200 & nCount_RNA < 25000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 10)
    
    # filtering
    if (sample %in% original_engineered_samples)
    {
      cells_filtered = intersect(original_selected_cells, colnames(organoid))
      organoid = organoid[, cells_filtered]
    }
    
    organoid = NormalizeData(organoid, verbose = F)
    organoid = FindVariableFeatures(organoid, selection.method = "vst", nfeatures = 3000, verbose = F)
    organoid.list[[i]] = organoid
  }
  
  # 2. integration feature selection
  features <- SelectIntegrationFeatures(object.list = organoid.list)
  organoid.list <- lapply(X = organoid.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  # 3. find anchors
  organoid.anchors <- FindIntegrationAnchors(object.list = organoid.list, anchor.features = features, reduction = "rpca")
  
  # 4. integrate the data
  organoid.combined <- IntegrateData(anchorset = organoid.anchors)
  DefaultAssay(organoid.combined) <- "integrated"
  organoid.combined <- ScaleData(organoid.combined, verbose = FALSE)
  organoid.combined = RunPCA(organoid.combined, npcs = 30)
  organoid.combined = RunUMAP(organoid.combined, reduction = "pca", dims = 1:30)
  # organoid.combined = RunUMAP(organoid.combined, reduction = "pca", dims = 1:30, n.components = 3L)
  
  # 5. find neighbors and clusters
  organoid.combined = FindNeighbors(organoid.combined, dim = 1:30)
  organoid.combined = FindClusters(organoid.combined, resolution = 1)
  organoid.combined = FindClusters(organoid.combined, resolution = 2)
  organoid.combined = FindClusters(organoid.combined, resolution = 0.5)
  organoid.combined = FindClusters(organoid.combined, resolution = 0.3)
  organoid.combined = FindClusters(organoid.combined, resolution = 3)
  
  return(organoid.combined)
  #saveRDS(organoid.combined, paste0(name, ".rds"))
}

combined_organoid = seurat_integration(engineered_ipsc_samples)
saveRDS(combined_organoid, 'engineered_combined_4timepoints.rds')

# add metadata
df_cells = read.table('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/ctype_v6.txt',
                      sep = '\t', header = 1, row.names = 1)
overlap_cells = intersect(rownames(df_cells), colnames(combined_organoid))
df_cells = df_cells[overlap_cells, ]

for (col in colnames(df_cells)){
  if (!col %in% colnames(combined_organoid@meta.data)){
    combined_organoid@meta.data[[col]] = ''
    combined_organoid@meta.data[overlap_cells, col] = df_cells[[col]]
  }
}

DimPlot(combined_organoid[, combined_organoid$timepoint_v2 != ''], 
        group.by = 'cellclass_draft_v2', label = T) + theme(legend.position = 'none')
DimPlot(combined_organoid, group.by = 'timepoint_v2', cols = c('red','grey','grey','grey')) + theme(legend.position = 'none')
combined_organoid$sample = factor(combined_organoid$sample, levels = unique(combined_organoid$sample))
DimPlot(combined_organoid, group.by = 'sample', cols = c(rep('grey',8), 'red','blue')) + theme(legend.position = 'none')

# transfer to h5ad 
DefaultAssay(combined_organoid) = "RNA"
SaveH5Seurat(combined_organoid, filename = "combined_all_eng_organoids.h5seurat")
Convert("combined_all_eng_organoids.h5seurat", dest = "h5ad")




