library(Seurat)
library(ggplot2)

# donor
donor = readRDS('donor_reintegrated.rds')

# make 3d umap
original_umap = donor@reductions$UMAP@cell.embeddings
donor = RunUMAP(donor, reduction = 'pca', dims = 1:30, n.components = 3, reduction.name = 'UMAP_3d')

df_meta = donor@meta.data
df_umap = donor@reductions$UMAP_3d@cell.embeddings
df_meta = data.frame(cbind(df_meta, df_umap))

df_umap_2d = donor@reductions$UMAP@cell.embeddings
df_meta = data.frame(cbind(df_meta, df_umap_2d))
write.table(df_meta, 'umap_3d_metadata_donor.txt', sep = '\t')

# all
all = readRDS('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_filtered_v3.rds')
all = RunUMAP(all, reduction = 'pca', dims = 1:30, n.components = 3, reduction.name = 'UMAP_3d')

df_umap_3d = all@reductions$UMAP_3d@cell.embeddings
df_meta = data.frame(cbind(all@meta.data, df_umap_3d))
write.table(df_meta, 'umap_3d_metadata_all.txt', sep = '\t')

# engineered
engineer = readRDS('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/organoid_integrated_normalized_C2_RPCA_highQuality_v5.RDS')
df_meta = all@meta.data
df_meta = df_meta[df_meta$sample_group %in% c('Engineered-Cre', 'Engineered-Flp'),]

engineer = RunUMAP(engineer, reduction = 'pca', dims = 1:30, n.components = 3, reduction.name = 'UMAP_3d')
df_umap_3d = engineer@reductions$UMAP_3d@cell.embeddings
df_meta = data.frame(cbind(df_meta, df_umap_3d))
write.table(df_meta, 'umap_3d_metadata_engineered.txt', sep = '\t')