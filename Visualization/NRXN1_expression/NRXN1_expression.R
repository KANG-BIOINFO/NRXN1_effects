# load environment
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(plyr)

source('../functions/_global.R')
source('../functions/_visualize.R')

# load data
path = 'single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/trajectory/'
time_points = c('d22', 'd50', 'd101')
donors = list()
for (i in 1:length(time_points)){
  time_point = time_points[i]
  donors[[i]] = readRDS(paste0(path, 'Donor_', time_point, '_reintegrated.rds'))
}

### check NRXN1 in donors
# d22 = donors[[1]]
# d50 = donors[[2]]
# d101 = donors[[3]]

violin_gene_times = function(donors, labels, gene, scale = 'area', trim = F)
{
  ps = list()
  
  for (i in 1:length(donors)){
    donor = donors[[i]]
    label = labels[[i]]
    donor = donor[, donor$cellclass_draft_v2 %in% cellclass_ordered_v2]
    
    # only select control
    donor = donor[, donor$genotype == 'Control']
    donor$cellclass_draft_v2 = factor(donor$cellclass_draft_v2, levels = cellclass_ordered_v2)
    df_meta = donor@meta.data[,c('cellclass_draft_v2', 'timepoint_v2')]
    df_meta['expr'] = donor@assays$RNA@data[gene, ]
    df_meta['combined'] = paste0(df_meta['timepoint_v2'], '-', df_meta['cellclass_draft_v2'])
    df_meta['nonzero_expr'] = df_meta$expr > 0
    df_meta$count = 1
    labs = aggregate(count~cellclass_draft_v2, df_meta, sum)
    labs_v2 = aggregate(nonzero_expr~cellclass_draft_v2, df_meta, sum)
    labs$nonzero_expr = labs_v2$nonzero_expr
    labs$nonzero_prop = round(labs$nonzero_expr / labs$count * 100, 1)

    # draw plot 
    y_max = max(df_meta$expr)
    p = ggplot(df_meta, aes(x = cellclass_draft_v2, y = expr, fill = cellclass_draft_v2)) + 
      geom_violin(scale = scale, trim = trim, lwd = 0.1) + 
      geom_text(data = labs, aes(x = cellclass_draft_v2, y = y_max, label = nonzero_prop), size = 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = 'none') + labs(y = paste0('expr_', label))
    # add jitter to each plot
    p = p + geom_jitter(shape = 16, size = 0.15, position = position_jitter(0.05)) +
          scale_fill_manual(values = colour_cellclass_v2)
    # append
    ps[[i]] = p
  }
  P = plot_grid(plotlist = rev(ps), nrow = length(donors), align = 'v', rel_heights = c(1/3, 1/3, 1/3))
  return(P)
}

P = violin_gene_times(donors, time_points, 'NRXN1', scale = 'count', trim = T)
ggsave('expr_NRXN1_donor_ctrl_3time_with_proportion.pdf', plot = P, width = 20, height = 28, units = 'cm')

P = violin_gene_times(donors, time_points, 'NRXN1', scale = 'area', trim = T)
ggsave('expr_NRXN1_donor_ctrl_3time_scale_byArea_with_proportion.pdf', plot = P, width = 20, height = 28, units = 'cm')


##################
### engineered ###
##################
# load data
path = '/data//aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered_new/engineered_4_timepoints/'
time_points = c('d23', 'd50', 'd112')
engs = list()
for (i in 1:length(time_points)){
  time_point = time_points[i]
  engs[[i]] = readRDS(paste0(path, 'Engineered_', time_point, '_reintegrated.rds'))
}

P2 = violin_gene_times(engs, time_points, 'NRXN1', scale = 'count', trim = T)
ggsave('expr_NRXN1_engineered_ctrl_3time_with_proportion.pdf', plot = P2, width = 20, height = 28, units = 'cm')

P2 = violin_gene_times(engs, time_points, 'NRXN1', scale = 'area', trim = T)
ggsave('expr_NRXN1_engineered_ctrl_3time_byArea_with_proportion.pdf', plot = P2, width = 20, height = 28, units = 'cm')


##############################
### engineered KO vs. Ctrl ###
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

########
# Engineered organoids NRXN1 per type regulation
celltypes_d23 = c('IPC1', 'IPC2', 'IPC3', 'IN3', 'IN5', 'CN1', 'CN3', 'CN4', 'CN5')
celltypes_d50 = c('oRG', 'IPC1','IPC2','IPC3','IN1','IN2','IN3','IN4','IN5','IN7','CN1','CN2','CN3','CN4','CN5')
celltypes_d112 = cellclass_ordered_v2

P3 = draw_boxplot(engs[[1]], celltypes_d23, 'NRXN1')
P4 = draw_boxplot(engs[[2]], celltypes_d50, 'NRXN1')
P5 = draw_boxplot(engs[[3]], celltypes_d112, 'NRXN1')

ggsave('violinplot_v2/NRXN1_Engineered_Cre_vs_Flp_d23.pdf', P3, width = 35, height = 15, units = 'cm')
ggsave('violinplot_v2/NRXN1_Engineered_Cre_vs_Flp_d50.pdf', P4, width = 35, height = 15, units = 'cm')
ggsave('violinplot_v2/NRXN1_Engineered_Cre_vs_Flp_d112.pdf', P5, width = 35, height = 15, units = 'cm')

# d160 engineered organoids
eng_d160 = readRDS('../../engineered_4_timepoints/Engineered_d160_integrated.rds')
df_meta = eng_d160@meta.data
df_meta$expr = eng_d160@assays$RNA@data['NRXN1',]
df_meta$sample = factor(df_meta$sample, levels = c('Flp-d160', 'Cre-d160'))
df_meta$cluster = df_meta$integrated_snn_res.0.5
P = ggviolin(df_meta, x = 'sample', y = 'expr', 
             color = 'sample', palette = c('#008FC4','#F8766D'), 
             add = 'jitter', facet.by = 'cluster', 
             ylab = 'NRXN1 expression') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 20)) + 
  stat_compare_means(method = 'wilcox.test', label = 'p.signif', 
                     label.x.npc = 'middle') 
ggsave('NRXN1_Engineered_Cre_vs_Flp_d160_clusters(resolution=0.5).pdf', P, width = 30, height = 30, units = 'cm')

# d50 = engs[[2]]
# d50_IN2 = d50[, d50$cellclass_draft_v2 == 'IN2']
# Idents(d50_IN2) = d50_IN2$sample
# DefaultAssay(d50_IN2) = 'RNA'
# VlnPlot(d50_IN2, features = c('NRXN1'))

################################################
# donor organoids NRXN1 per type regulation
P6 = draw_boxplot(donors[[1]], cellclass_ordered_v2, 'NRXN1')
P7 = draw_boxplot(donors[[2]], cellclass_ordered_v2, 'NRXN1')
P8 = draw_boxplot(donors[[3]], cellclass_ordered_v2, 'NRXN1')

ggsave('violinplot_v2/NRXN1_Donor_NRXN1del_vs_Ctrl_d22.pdf', P6, width = 35, height = 15, units = 'cm')
ggsave('violinplot_v2/NRXN1_Donor_NRXN1del_vs_Ctrl_d50.pdf', P7, width = 35, height = 15, units = 'cm')
ggsave('violinplot_v2/NRXN1_Donor_NRXN1del_vs_Ctrl_d101.pdf', P8, width = 35, height = 15, units = 'cm')

type = 'CN2'
d50_type = d50[, d50$cellclass_draft_v2 == type]
Idents(d50_type) = d50_type$sample
DefaultAssay(d50_type) = 'RNA'
VlnPlot(d50_type, features = c('NRXN1'))

d50_engineered = engs[[2]]
Idents(d50_engineered) = d50_engineered$sample
VlnPlot(d50_engineered, features = c('nFeature_RNA'))


###############################################
## create barcharts for percentage of NRXN1 ###
###############################################
draw_barchart_genePercentage = function(donor, top_n = 10, gene = 'NRXN1', max_color = '#6600CC', control_only = T){
  donor = donor[, donor$cellclass_draft_v2 %in% cellclass_ordered_v2]
  if (control_only){
    donor = donor[, donor$genotype == 'Control']
  }
  df_meta = donor@meta.data[, c('cellclass_draft_v2', 'genotype')]
  df_meta$expr = donor@assays$RNA@data[gene, ]
  df_meta$nonzero_expr = df_meta$expr > 0
  df_meta$count = 1
  labs = aggregate(count~cellclass_draft_v2, df_meta, sum)
  labs_v2 = aggregate(nonzero_expr~cellclass_draft_v2, df_meta, sum)
  labs$nonzero_expr = labs_v2$nonzero_expr
  labs$nonzero_prop = round(labs$nonzero_expr / labs$count * 100, 1)
  labs_ordered = labs[order(labs$nonzero_prop),]
  labs_ordered[['cellclass_draft_v2']] = factor(labs_ordered[['cellclass_draft_v2']],
                                                levels = labs_ordered[['cellclass_draft_v2']])
  
  p = ggplot(data = labs_ordered, aes(x = cellclass_draft_v2, y = nonzero_prop, fill = nonzero_prop)) + 
      geom_bar(stat = 'identity') + coord_flip() + theme_minimal() +
    scale_fill_gradient2(low = '#FFFFFF', high = max_color) +
    ylab(paste0('percentage of ', gene, ' expressing cells')) + 
    xlab('cell class') 
  return(p)
}

time_points = c('d22', 'd50', 'd101')
gene = 'NRXN1'
for (i in 1:length(donors))
{
  donor = donors[[i]]
  time_point = time_points[[i]]
  # p = draw_barchart_genePercentage(donor, gene = gene)
  # ggsave(paste0('Barplot_percentile_DonorControl_',time_point, '_', gene, '.pdf'), p,
  #        width = 18, height = 10, units = 'cm')
  p = draw_barchart_genePercentage(donor, gene = gene)
  p = p + ylim(0,100)
  ggsave(paste0('Barplot_percentile_DonorControl_',time_point, '_', gene, '_sampeScale.pdf'), p,
         width = 18, height = 10, units = 'cm')
}




