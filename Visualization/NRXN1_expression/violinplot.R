# load environment
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(plyr)

source('../../functions/_global.R')
source('../../functions/_visualize.R')

# load data
path = '../v2_filtered/trajectory/'
time_points = c('d22', 'd50', 'd101')
donors = list()
for (i in 1:length(time_points)){
  time_point = time_points[i]
  donors[[i]] = readRDS(paste0(path, 'Donor_', time_point, '_reintegrated.rds'))
}

##############################
### donor KO vs. Ctrl ########
##############################
draw_boxplot = function(donor, cell_types, gene = 'NRXN1'){
  
  donor = donor[, donor$cellclass_draft_v2 %in% cell_types]
  df_meta = data.frame(donor@meta.data[, c('cellclass_draft_v2', 'genotype')])
  colnames(df_meta) = c('cell_type', 'genotype')
  df_meta$genotype = mapvalues(df_meta$genotype, from = c('NRXN_del'), to = ('NRXN1_del'))
  df_meta$cell_type = factor(df_meta$cell_type, levels = cell_types)
  df_meta$genotype = factor(df_meta$genotype, levels = c('Control', 'NRXN1_del'))
  df_meta['expr'] = donor@assays$RNA@data[gene, ]
  
  
  # draw box plot
  # ggplot(data = df_meta, aes(x = genotype, y = expr, fill = genotype)) + 
  #   geom_violin() + 
  #   stat_compare_means(method = 'wilcox.test')
  
  P = ggviolin(df_meta, x = 'genotype', y = 'expr', 
             color = 'genotype', palette = c('#008FC4','#F8766D'), 
             add = 'jitter', facet.by = 'cell_type', 
             ylab = 'NRXN1 expression', size = 0.5, font.label = list(size = 18)) +
      theme(legend.position = 'none',
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      stat_compare_means(method = 'wilcox.test', label = 'p.signif', 
                         label.x.npc = 'middle') 
  P2 = facet(P, facet.by = 'cell_type', ncol = 9, nrow = 3)
  return(P2) 
}

draw_boxplot_acrossTime = function(donors, cell_types, gene){
  p1 = draw_boxplot(donors[[1]], cell_types, gene)
  p2 = draw_boxplot(donors[[2]], cell_types, gene)
  p3 = draw_boxplot(donors[[3]], cell_types, gene)
  
  ggsave(paste0('Donor_d22_comparison_', gene, '.pdf'), p1, width = 35, height = 15, units = 'cm')
  ggsave(paste0('Donor_d50_comparison_', gene, '.pdf'), p2, width = 35, height = 15, units = 'cm')
  ggsave(paste0('Donor_d101_comparison_', gene, '.pdf'), p3, width = 35, height = 15, units = 'cm')
}

################################################
genes = c('PTN', 'NNAT', 'BASP1', 'NDFIP1', 'SKP1', 'SUMO1', 'SUMO2', 'UBB', 'GABARAPL2',
          'UBE2M', 'UBE2V2', 'RNF7', 'STUB1', 'UCHL1')
for (gene in genes)
{
  draw_boxplot_acrossTime(donors, cellclass_ordered_v2, gene)
}
