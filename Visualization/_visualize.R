library(Seurat)
library(ggplot2)
library(cowplot)

### visualize markers (stacked violin) ###
visualize_markers_violin = function(seu_obj, cell_type_key, genes)
{
  DefaultAssay(seu_obj) = 'RNA'
  Idents(seu_obj) = seu_obj@meta.data[[cell_type_key]]
  VlnPlot(seu_obj, features = genes, pt.size = 0, stack = T, flip = T) + 
    theme(legend.position = 'none')
}

### visualize markers (dot plot) ###
visualize_markers_dot = function(seu_obj, cell_type_key, genes)
{
  # color indicates levels; size indicates proportion
  DefaultAssay(seu_obj) = 'RNA'
  Idents(seu_obj) = seu_obj@meta.data[[cell_type_key]]
  DotPlot(seu_obj, features = genes, scale = F) + 
    theme(axis.text.x = element_text(angle = 90))
}

### visualize metadata on umap ####
visualize_dimplot = function(seu_obj, key, colors, output_name, pt.size = 0.25, alpha = 1, width = 10, height = 10)
{
  g = DimPlot(seu_obj, group.by = key, cols = colors, pt.size = pt.size) + 
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      ggtitle("")
  g[[1]]$layers[[1]]$aes_params$alpha = alpha
  ggsave(paste0(output_name, '.png'), width = width, height = height)
}


### barplot of cell abundance ###
plot_2way_barplot = function(df_cells, cell_orders, x = "cellclass_draft_v2", y = "genotype", 
                             y_scale = T, output_name = 'abundance', height = 6, width = 8)
{
  df_freq = data.frame(table(df_cells[[x]], df_cells[[y]]))
  colnames(df_freq) = c(x, y, "n_cells")
  df_freq[["log2_n_cells"]] = log2(df_freq[["n_cells"]]+1)
  df_freq[[x]] = factor(df_freq[[x]], levels = cell_orders)
  count_name = "n_cells"
  if (y_scale){count_name = "log2_n_cells"}
  # pdf(paste0(output_name, '.pdf'), height = height, width = width)
  # png(paste0(output_name, '.png'), height = height, width = width)
  p1 = ggplot(df_freq, aes(x = df_freq[[x]], y = df_freq[[count_name]], fill = df_freq[[y]])) + 
    geom_bar(stat = "identity", position = "fill") + theme_bw()  + ylab("proportion") +
    theme(axis.text.y = element_text(family = "Helvetica", size = 12),
          axis.text.x = element_text(family = "Helvetica", angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_blank(), 
          legend.position = "none",
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) + 
    guides(fill = guide_legend( title = y )) + 
    scale_fill_manual("legend", values = c("Control" = "#99CCFF", "NRXN_del" = "#FF6666")) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "black", size = 0.6) 
  
  p2 = ggplot(df_freq, aes(x = df_freq[[x]], y = df_freq[[count_name]]), color = "#C0C0C0") +
    geom_bar(stat = "identity") + theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(family = "Helvetica", size = 12),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  
  P = plot_grid(p2, p1, nrow = 2, align = 'v', rel_heights = c(1/4, 3/4))
  ggsave(paste0(output_name, '.pdf'), P, height = height, width = width)
  # dev.off()
}


### visualize cell proportion changes (milo) ### 
plot_milo_umap = function(milo, colour_by = 'genotype', text_by = 'cellclass_draft_v2')
{
  milo = buildNhoodGraph(milo)
  umap_pl <- plotReducedDim(milo, dimred = "UMAP", colour_by = colour_by, text_by = text_by, 
                            text_size = 3, point_size=0.5) +guides(fill="none")
  nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1) 
  umap_pl + nh_graph_pl + plot_layout(guides="collect")
  da_results = annotateNhoods(milo, da_results, coldata_col = 'cellclass_draft_v2')
  ggplot(da_results, aes(cellclass_draft_v2_fraction)) + geom_histogram(bins = 50)
}

plot_milo_beeswarm = function(milo, da_results)
{
  da_results$celltype = ifelse(da_results$cellclass_draft_v2_fraction < 0.7, 'Mixed', da_results$cellclass_draft_v2)
  abundant_celltypes = rownames(table(da_results$celltype))[table(da_results$celltype) > 20]
  da_results_filtered = da_results[da_results$celltype %in% abundant_celltypes, ]
  plotDAbeeswarm(da_results_filtered, group.by = 'celltype', alpha = 1)
}





