library(ggplot2)
source('../functions/_global.R')
source('../functions/_visualize.R')

# load data
df_cells = read.table('../../../v2_filtered/ctype_v6.txt',
                      sep = '\t', header = 1, row.names = 1)

df_cells = df_cells[df_cells$cellclass_draft_v2 %in% cellclass_ordered_v2,]
df_cells$cellclass_draft_v2 = factor(df_cells$cellclass_draft_v2, levels = cellclass_ordered_v2)

# create sample group v2 annotations
df_cells$sample_group_v2 = paste0(df_cells$sample_group, '-',df_cells$timepoint_v2)
df_sample_info = df_cells[, c('sample', 'sample_group_v2')]
df_sample_info = distinct(df_sample_info)
rownames(df_sample_info) = df_sample_info$sample
o = 1:12
o2 = 1:8
names(o) = sample_group_ordered
names(o2) = sample_group_ordered_v2
df_sample_info$order = o[df_sample_info$sample_group_v2]
df_sample_info$order_v2 = o2[df_sample_info$sample_group_v2]
df_sample_info = df_sample_info[order(df_sample_info$order),]

# cell type proportion across all samples
df_counts = data.frame(table(df_cells$sample, df_cells$cellclass_draft_v2))
colnames(df_counts) = c('sample', 'cell_class', 'n_cells')
df_counts$sample = factor(df_counts$sample, levels = df_sample_info$sample)
p1 = ggplot(data = df_counts, aes(x = sample, y = n_cells, fill = cell_class)) +
      geom_bar(stat = 'identity', position = 'fill') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
            panel.background = element_rect(fill = 'white', colour = 'black',
                                            size = 1, linetype = 'solid')) +
      ylab('Cell Proportion') + 
      scale_y_continuous(expand = c(0,0)) +    
      scale_fill_manual(values = colour_cellclass_v2)
ggsave('stacked_barplot/all_samples.pdf', p1, height = 15, width = 18, units = 'cm')


# cell type proportion across donor + eng_d112
df_sample_info_v2 = df_sample_info[!is.na(df_sample_info$order_v2),]
df_sample_info_v2 = df_sample_info_v2[order(df_sample_info_v2$order_v2),]
df_counts_v2 = df_counts[df_counts$sample %in% df_sample_info_v2$sample,]
df_counts_v2$sample = factor(df_counts_v2$sample, levels = df_sample_info_v2$sample)
p2 = ggplot(data = df_counts_v2, aes(x = sample, y = n_cells, fill = cell_class)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = 'white', colour = 'black',
                                        size = 1, linetype = 'solid')) +
  ylab('Cell Proportion') + 
  scale_y_continuous(expand = c(0,0)) +    
  scale_fill_manual(values = colour_cellclass_v2)
ggsave('stacked_barplot/selected_samples.pdf', p2, height = 15, width = 18, units = 'cm')
