# load environment
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(plyr)
library(readxl)

source('../../functions/_global.R')
source('../../functions/_visualize.R')

# load data
fetal = readRDS('single_cell_projects/10X-Pak/reference/fetalBrain_Bhaduri_Nature/public_h5ad/Bhaduri_raw.rds')
fetal$development_stage_v2 = gsub(' post-fertilization human stage', '', fetal$development_stage)
df_cells = fetal@meta.data

# create bar charts
draw_barchart_genePercentage = function(organoid, gene = 'NRXN1', max_color = '#6600CC'){

  df_meta = organoid@meta.data[, c('cell_type', 'development_stage')]
  df_meta$expr = organoid@assays$RNA@data[gene, ]
  df_meta$nonzero_expr = df_meta$expr > 0
  df_meta$count = 1
  labs = aggregate(count~cell_type, df_meta, sum)
  labs_v2 = aggregate(nonzero_expr~cell_type, df_meta, sum)
  labs$nonzero_expr = labs_v2$nonzero_expr
  labs$nonzero_prop = round(labs$nonzero_expr / labs$count * 100, 1)
  labs_ordered = labs[order(labs$nonzero_prop),]
  labs_ordered[['cell_type']] = factor(labs_ordered[['cell_type']],
                                                levels = labs_ordered[['cell_type']])
  
  p = ggplot(data = labs_ordered, aes(x = cell_type, y = nonzero_prop, fill = nonzero_prop)) +
      geom_bar(stat = 'identity') + coord_flip() + theme_minimal() +
    scale_fill_gradient2(low = '#FFFFFF', high = max_color) +
    ylab(paste0('percentage of ', gene, ' expressing cells')) +
    xlab('cell class')
  
  return(p)
}

# run it
gene = 'NRXN1'

for (stage in unique(df_cells$development_stage_v2))
{
  organoid = fetal[, fetal$development_stage_v2 == stage]
  # p = draw_barchart_genePercentage(organoid, gene = gene)
  # ggsave(paste0('Barplot_percentile_', stage, '_', gene, '.pdf'), p,
  #        width = 18, height = 10, units = 'cm')
  p = draw_barchart_genePercentage(organoid, gene = gene, max_color = '#006400')
  p = p + ylim(0,100)
  ggsave(paste0('v2_same_scale/Barplot_percentile_', stage, '_', gene, '.pdf'), p,
         width = 18, height = 10, units = 'cm')
}

#####################################
#####################################
#####################################
# run it for another annotation group
df_origin = read_excel('./41586_2021_3910_MOESM3_ESM.xlsx', sheet = 1)
overlap_cells = intersect(rownames(df_cells), df_origin$cell.name)

# transfer annotations
df_cells_v2 = df_cells[overlap_cells, ]

df_origin_v2 = as.data.frame(df_origin[df_origin$cell.name %in% overlap_cells, ])
rownames(df_origin_v2) = df_origin_v2$cell.name
df_origin_v2 = df_origin_v2[overlap_cells,]

df_cells_v2$cell_type_v2 = df_origin_v2$`cell type`
df_match = as.data.frame(table(df_cells_v2$cell_type, df_cells_v2$cell_type_v2))
colnames(df_match) = c('Annotation_1', 'Annotation_2', 'Number of cells')
write.table(df_match, 'Number_cells_in_two_annotations.txt', sep = '\t')

# create bar charts
draw_barchart_genePercentage = function(organoid, gene = 'NRXN1', max_color = '#6600CC'){
  
  df_meta = organoid@meta.data[, c('cell_type_v2', 'development_stage')]
  df_meta$expr = organoid@assays$RNA@data[gene, ]
  df_meta$nonzero_expr = df_meta$expr > 0
  df_meta$count = 1
  labs = aggregate(count~cell_type_v2, df_meta, sum)
  labs_v2 = aggregate(nonzero_expr~cell_type_v2, df_meta, sum)
  labs$nonzero_expr = labs_v2$nonzero_expr
  labs$nonzero_prop = round(labs$nonzero_expr / labs$count * 100, 1)
  labs_ordered = labs[order(labs$nonzero_prop),]
  labs_ordered[['cell_type_v2']] = factor(labs_ordered[['cell_type_v2']],
                                       levels = labs_ordered[['cell_type_v2']])
  
  p = ggplot(data = labs_ordered, aes(x = cell_type_v2, y = nonzero_prop, fill = nonzero_prop)) +
    geom_bar(stat = 'identity') + coord_flip() + theme_minimal() +
    scale_fill_gradient2(low = '#FFFFFF', high = max_color) +
    ylab(paste0('percentage of ', gene, ' expressing cells')) +
    xlab('cell class')
  
  return(p)
}

# run it
fetal = fetal[, overlap_cells]
fetal@meta.data = df_cells_v2
fetal = fetal[, fetal$cell_type_v2 != 'Outlier']
df_abundance = as.data.frame(table(fetal$cell_type_v2, fetal$development_stage_v2))
colnames(df_abundance) = c('Annotation_v2', 'Timepoint', 'Number of cells')
write.table(df_abundance, 'Number_cells_across_time.txt', sep = '\t')
gene = 'NRXN1'

for (stage in unique(df_cells$development_stage_v2))
{
  organoid = fetal[, fetal$development_stage_v2 == stage]
  p = draw_barchart_genePercentage(organoid, gene = gene, max_color = '#006400')
  p = p + ylim(0,100)
  ggsave(paste0('v2_same_scale/Barplot_percentile_', stage, '_', gene, '_AnnotationV2.pdf'), p,
         width = 18, height = 10, units = 'cm')
}


