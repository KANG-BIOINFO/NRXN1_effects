library(Seurat)
library(ggplot2)

source('../functions/_global.R')
source('../functions/_visualize.R')

# umap for combined organoids can be found in /data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/visualize/visualize.R
# umap for donor
donor = readRDS('../../donor_reintegrate/donor_reintegrated.rds')
donor = donor[, donor$cellclass_draft_v2 %in% cellclass_ordered_v2]
donor$cellclass_draft_v2 = factor(donor$cellclass_draft_v2, levels = cellclass_ordered_v2)
donor$timepoint_v2 = factor(donor$timepoint_v2, levels = c('d022','d050','d101'))
donor$genotype = factor(donor$genotype, levels = c('Control', 'NRXN_del'))

# dimplot
visualize_dimplot(seu_obj = donor, key = 'cellclass_draft_v2', colors = colour_cellclass_v2,
                  output_name = 'donor_umap_all_cellclass')

# dimplot (timepoint)
visualize_dimplot(seu_obj = donor, key = 'timepoint_v2', colors = colour_timepoint_v2,
                  alpha = 0.5, output_name = 'donor_umap_all_timepoint')

# dimplot (genotype)
visualize_dimplot(seu_obj = donor, key = 'genotype', colors = colour_genotype,
                  alpha = 0.5, output_name = 'donor_umap_genotype')

# umap for combined organoids
organoid = readRDS('../../../v2_filtered/combined_all_organoids_filtered_v3.rds')
organoid = organoid[, organoid$cellclass_draft_v2 %in% cellclass_ordered_v2]
organoid$cellclass_draft_v2 = factor(organoid$cellclass_draft_v2, levels = cellclass_ordered_v2)
df_umap = data.frame(organoid@reductions$umap@cell.embeddings)
organoid = organoid[, rownames(df_umap[df_umap$UMAP_1 < 10, ])]

p = DimPlot(organoid, group.by = 'cellclass_draft_v2', cols = colour_cellclass_v2) + 
  theme(legend.position = 'none') + ggtitle('')
ggsave('donor_engineered_umap_all.png', p, width = 20, height = 20, units = 'cm')

organoid$timepoint_v2 = factor(organoid$timepoint_v2, levels = c('d022','d023','d050','d101','d112'))
p2 = DimPlot(organoid, group.by = 'timepoint_v2', 
             cols = c('#006600', '#006600', '#3399FF', '#CC0000', '#CC0000')) + 
  theme(legend.position = 'none') + ggtitle('')
ggsave('donor_engineered_umap_all_timepoints.png', p2, width = 20, height = 20, units = 'cm')

p3 = DimPlot(organoid, group.by = 'sample_group') + 
  theme(legend.position = 'none') + ggtitle('')
ggsave('donor_engineered_umap_all_sample_groups.png', p3, width = 20, height = 20, units = 'cm')
