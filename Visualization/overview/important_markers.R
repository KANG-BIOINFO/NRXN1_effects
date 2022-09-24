library(Seurat)
library(ggplot2)

source('../functions/_global.R')
source('../functions/_visualize.R')

# load data
organoid = readRDS('../../../v2_filtered/combined_all_organoids_filtered_v3.rds')
organoid = organoid[, organoid$cellclass_draft_v2 %in% cellclass_ordered_v2]
organoid$cellclass_draft_v2 = factor(organoid$cellclass_draft_v2, levels = cellclass_ordered_v2)
df_umap = data.frame(organoid@reductions$umap@cell.embeddings)
organoid = organoid[, rownames(df_umap[df_umap$UMAP_1 < 10, ])]
saveRDS(organoid, 'combined_all_organoids_filtered_v3_filtered.rds')

markers = c('S100B', 'ALDOC', 'HES5', 'CA2', 'TNC', 'HOPX', 'MOXD1',
            'FAM107A', 'HES1', 'CRYAB', 'FBXO32', 'OLIG1', 'OLIG2', 
            'EOMES', 'TBR2', 'NEUROG1', 'PPP1R17', 'BTG2', 'SLC17A7',
            'SCL17A6', 'GAD1', 'GAD2', 'CALB1', 'SST', 'CALB2', 'VGAT')
markers_v2 = c('SLC17A6','ALDH1L1','AQP4','GFAP','SOX9','SPARCL1','WNT7B','SATB2','FEZF2','BCL11B','NEUROD2','SYT1','SYT4','DLX5','MIK67','TOP2A','BIRC5',
                'PTTG1', 'PCNA', 'HMGB2', 'CDK1', 'MIK67', 'CCNB2', 'SOX2', 'HES1', 'VIM')

DefaultAssay(organoid) = 'RNA'
# genes = rownames(organoid)
# for (marker in markers)
# {
#   if (marker %in% genes)
#   {
#     p = FeaturePlot(organoid, features = marker,raster = F, cols = c('#F3F3F3', '#990000'))
#     ggsave(paste0('markers/umap_',marker,'.png'), p, height = 10, width = 10, units = 'cm')
#   }
# }

genes = rownames(organoid)
for (marker in markers_v2)
{
  if (marker %in% genes)
  {
    p = FeaturePlot(organoid, features = marker,raster = F, cols = c('#F3F3F3', '#990000'))
    ggsave(paste0('markers_v2/umap_',marker,'.png'), p, height = 10, width = 10, units = 'cm')
  }
}
