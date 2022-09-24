library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)

source('../../v2_filtered_new/Figures/functions/_global.R')
source('../../v2_filtered_new/Figures/functions/_visualize.R')

combined_organoid = readRDS("../combined_all_organoids_filtered.rds")
cds = as.cell_data_set(combined_organoid, assay = "integrated")
cds = cluster_cells(cds)
p1 = plot_cells(cds)
p1

cds = learn_graph(cds)
plot_cells(cds)

cds = order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime")
saveRDS(cds, "cds_monocle3.rds")


# plots
cds = readRDS('cds_monocle3_v2.rds')

# remove cell types
df_meta = read.table('../ctype_v6.txt', sep = '\t', header = 1, row.names = 1)
selected = rownames(df_meta[df_meta$cellclass_draft_v2 %in% cellclass_ordered_v2, ])
cds = cds[, selected]

p = plot_cells(cds[, !is.infinite(pseudotime(cds))], 
           color_cells_by = 'pseudotime',
           show_trajectory_graph = F) + 
          theme(axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line = element_blank())
ggsave('trajectory_monocle3_pseudotime_v3.png', p, width = 14, height = 12, units = 'cm')

# highlight specific pseudotime regions
df_meta = df_meta[selected, ]
df_meta$pseudotime = pseudotime(cds)

# highlight one region
highlight_pseudotime = function(df_meta, cds, min_pseudo, max_pseudo)
{
  highlighted = rownames(df_meta[(df_meta$pseudotime > min_pseudo) & (df_meta$pseudotime < max_pseudo),])
  cds@colData$highlight = 'others'
  cds[, highlighted]@colData$highlight = 'highlight'
  
  ps = plot_cells(cds[, !is.infinite(pseudotime(cds))], 
                 color_cells_by = 'highlight',
                 show_trajectory_graph = F,
                 label_roots = F, 
                 label_leaves = F, 
                 label_branch_points = F,
                 label_groups_by_cluster = F,
                 trajectory_graph_color = c('grey', 'red')) + 
        theme(axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank())
  ggsave(paste0('trajectory_monocle3_pseudotime_highlight_', min_pseudo, '_', max_pseudo, '.png'), ps, 
         width = 14, height = 12, units = 'cm')
}

highlight_pseudotime(df_meta, cds, 17, 20)
highlight_pseudotime(df_meta, cds, 32, 36)

# highlight two regions
cds@colData$highlight = 'others'
# highlighted1 = rownames(df_meta[(df_meta$pseudotime > 17) & (df_meta$pseudotime < 20),])
# highlighted2 = rownames(df_meta[(df_meta$pseudotime > 32) & (df_meta$pseudotime < 36),])

# cds[, highlighted1]@colData$highlight = 'pseudotime_17-20'
# cds[, highlighted2]@colData$highlight = 'pseudotime_32-36'

highlighted1 = rownames(df_meta[(df_meta$pseudotime > 19) & (df_meta$pseudotime < 23),])
highlighted2 = rownames(df_meta[(df_meta$pseudotime > 28) & (df_meta$pseudotime < 32),])

cds[, highlighted1]@colData$highlight = 'pseudotime_19-23'
cds[, highlighted2]@colData$highlight = 'pseudotime_28-32'

ps = plot_cells(cds[, !is.infinite(pseudotime(cds))], 
                color_cells_by = 'highlight',
                show_trajectory_graph = F,
                group_label_size = 0,
                label_roots = F, 
                label_leaves = F, 
                label_branch_points = F,
                label_groups_by_cluster = F) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())

ps = ps + scale_color_manual(values = c('#F5F5F5', '#99CCFF', '#006666'))
ggsave(paste0('trajectory_monocle3_pseudotime_highlight_two_regions_eng.png'), ps, 
       width = 12, height = 12, units = 'cm')
