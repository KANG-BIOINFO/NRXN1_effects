library(Seurat)
library(ggplot2)

source('../functions/_global.R')
source('../functions/_visualize.R')

# load data for integrated organoid samples
donor = readRDS('../../donor_reintegrate/donor_reintegrated.rds')
donor = donor[, donor$cellclass_draft_v2 %in% cellclass_ordered_v2]
donor$cellclass_draft_v2 = factor(donor$cellclass_draft_v2, levels = cellclass_ordered_v2)
donor$timepoint_v2 = factor(donor$timepoint_v2, levels = c('d022','d050','d101'))
donor$genotype = factor(donor$genotype, levels = c('Control', 'NRXN_del'))

# load data for each timepoint
path = '/data//aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/trajectory/'
time_points = c('d22', 'd50', 'd101')
donors = list()
for (i in 1:length(time_points)){
  time_point = time_points[i]
  donors[[i]] = readRDS(paste0(path, 'Donor_', time_point, '_reintegrated.rds'))
}

for (i in 1:length(donors))
{
  time_point = time_points[[i]]
  donor_i = donors[[i]]
  donor_i = donor_i[, donor_i$cellclass_draft_v2 %in% cellclass_ordered_v2]
  donor_i$cellclass_draft_v2 = factor(donor_i$cellclass_draft_v2, levels = cellclass_ordered_v2)
  donor_i$genotype = factor(donor_i$genotype, levels = c('Control', 'NRXN_del'))
  
  # dimplot (cell class)
  visualize_dimplot(seu_obj = donor_i, key = 'cellclass_draft_v2', colors = colour_cellclass_v2,
                    output_name = paste0('donor_umap_all_cellclass_', time_point))
  # dimplot (genotype)
  visualize_dimplot(seu_obj = donor_i, key = 'genotype', colors = colour_genotype,
                    alpha = 0.5, output_name = paste0('donor_umap_genotype', time_point))
}


# visualize marker genes 
df_umap = data.frame(donor@reductions$UMAP@cell.embeddings)
remained_cells = rownames(df_umap[df_umap$UMAP_1 < 7.5, ])
donor = donor[, remained_cells]
DefaultAssay(donor) = 'RNA'
FeaturePlot(donor, features = genes_npc, cols = c('#F7F7F7','#872729'))
FeaturePlot(donor, features = genes_glia, cols = c('#F7F7F7','#872729'))
FeaturePlot(donor, features = c('STMN2','GAP43','DCX'), cols = c('#F7F7F7','#872729'), ncol = 3)
FeaturePlot(donor, features = c('SLC32A1','GAD1','GAD2'), cols = c('#F7F7F7','#872729'), ncol = 3)
FeaturePlot(donor, features = c('NEUROD2','NEUROD6', 'FOXP2'), cols = c('#F7F7F7','#872729'), ncol = 3)
