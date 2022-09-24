library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(ggridges)

source('../functions/_global.R')
source('../functions/_visualize.R')

# load data for integrated organoid samples
donor = readRDS('../../donor_reintegrate/donor_reintegrated.rds')
donor = donor[, donor$cellclass_draft_v2 %in% cellclass_ordered_v2]
donor$cellclass_draft_v2 = factor(donor$cellclass_draft_v2, levels = cellclass_ordered_v2)
donor$timepoint_v2 = factor(donor$timepoint_v2, levels = c('d022','d050','d101'))
donor$genotype = factor(donor$genotype, levels = c('Control', 'NRXN_del'))

# transfer to monocle3
cds = as.cell_data_set(donor)
cds = cluster_cells(cds)
cds = learn_graph(cds)
cds = order_cells(cds)
saveRDS(cds, 'pseudotime_donors.rds')

# remove outliers on UMAP
umaps = as.data.frame(donor@reductions$UMAP@cell.embeddings)
selected = rownames(umaps[umaps$UMAP_1 < 7.5,])
cds = cds[, selected]
p = plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = F,
               label_roots = F, label_leaves = F, label_branch_points = F) + 
              theme(axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.line = element_blank())
ggsave('UMAP_donor_pseudotime_v2.pdf', p, width = 7.5, height = 7.5, units = 'cm')
ggsave('UMAP_donor_pseudotime_v2.png', p, width = 18, height = 20, units = 'cm')

# show cell abundance density
df = as.data.frame(colData(cds))
df$pseudotime = pseudotime(cds)
df = df[is.finite(df$pseudotime),]


# pseudotime distribution of multiple levels
p2 = ggplot(df, aes(x = pseudotime, y = timepoint_v2)) +
      geom_density_ridges(aes(fill = paste(timepoint_v2, genotype)), alpha = 0.8) + 
      theme(legend.position = c(0.1,0.9),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
ggsave('Ridge_plot_donor_pseudotime.pdf', p2, width = 20, height = 18, units = 'cm')


# show pseudotime for a gene
gene = 'NRXN1'
df$expr = donor[,rownames(df)]@assays$RNA@data[gene,]
df$expr_scale = donor[, rownames(df)]@assays$integrated@scale.data[gene, ]
ggplot(df, aes(x = pseudotime, y = expr_scale, colour = genotype)) + 
  geom_point(size = 0.1) + 
  geom_smooth(method = lm, linetype = 'dashed', color = 'darkred', fill = 'blue') + 
  scale_fill_manual(values = c('red', 'blue')) + 
  theme(legend.position = c(0.1,0.8)) +
  theme_bw()

write.table(df, 'donor_NRXN1_pseudotime_expr.txt', sep = '\t')

gene = 'PTN'
df$expr = donor[,rownames(df)]@assays$RNA@data[gene,]
df$expr_scale = donor[, rownames(df)]@assays$integrated@scale.data[gene, ]
write.table(df, 'donor_PTN_pseudotime_expr.txt', sep = '\t')

# plot
gene = 'NRXN1'
df$expr = donor[,rownames(df)]@assays$RNA@data[gene,]
df$expr_scale = donor[, rownames(df)]@assays$integrated@scale.data[gene, ]
df_ctrl = df[df$genotype == 'Control', ]
p = ggplot(df_ctrl, aes(x = pseudotime, y = expr_scale, colour = cellclass_draft_v2)) + 
        geom_point(size = 0.6) + geom_jitter(size = 0.6) + theme(legend.position = c(0.1,0.8)) + theme_bw()
ggsave('pseudotime_donor_NRXN1_cell_type_control.png', p, width = 17, height = 10, units = 'cm')
ggsave('pseudotime_donor_NRXN1_cell_type_control.pdf', p, width = 17, height = 10, units = 'cm')