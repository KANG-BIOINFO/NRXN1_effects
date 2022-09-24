library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggridges)
library(cowplot)

# find trajectory of d101 only organoids
organoids = readRDS("../../v2_filtered/combined_all_organoids_filtered_v3.rds")

organoids_eng = organoids[, organoids$source == "Engineered"]
organoids_d112 = organoids_eng[, organoids_eng$timepoint_v2 == "d112"]
organoids_d50 = organoids_eng[, organoids_eng$timepoint_v2 == "d050"]
organoids_d23 = organoids_eng[, organoids_eng$timepoint_v2 == "d023"]

# find lineage on d101
reintegration = function(organoids_dx, split.by = "sample")
{
  DefaultAssay(organoids_dx) = "RNA"
  
  organoids_dx <- SplitObject(organoids_dx, split.by = split.by)
  for (i in seq_along(organoids_dx)) {
    organoids_dx[[i]] <- NormalizeData(organoids_dx[[i]]) %>% FindVariableFeatures()
  }
  features <- SelectIntegrationFeatures(organoids_dx)
  for (i in seq_along(along.with = organoids_dx)) {
    organoids_dx[[i]] <- ScaleData(organoids_dx[[i]], features = features) %>% RunPCA(features = features)
  }
  anchors <- FindIntegrationAnchors(organoids_dx, anchor.features = features, reduction = "rpca", dims = 1:30)
  integrated <- IntegrateData(anchors, dims = 1:30)
  
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated)
  integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
  integrated <- FindNeighbors(integrated, dims = 1:30)
  # integrated <- FindClusters(integrated)
  return (integrated)
}

new_d112 = reintegration(organoids_d112)
new_d50 = reintegration(organoids_d50)
new_d23 = reintegration(organoids_d23)
saveRDS(new_d23, "Engineered_d23_reintegrated.rds")
saveRDS(new_d50, "Engineered_d50_reintegrated.rds")
saveRDS(new_d112, "Engineered_d112_reintegrated.rds")

# DimPlot(new_d50, group.by = "cellclass_draft_v2", label = T) + theme(legend.position = "none")
# DimPlot(new_d50, group.by = "genotype")


# monocle3
run_monocle3 = function(seu)
{
  cds = as.cell_data_set(seu)
  cds = cluster_cells(cds)
  cds = learn_graph(cds = cds)
  cds = order_cells(cds)
  return(cds)
}

cds_d101 = run_monocle3(new_d101)
cds_d50 = run_monocle3(new_d50)
cds_d22 = run_monocle3(new_d22)
saveRDS(cds_d22, "pseudotime_d22.rds")
saveRDS(cds_d50, "pseudotime_d50.rds")
saveRDS(cds_d101, "pseudotime_d101.rds")

plot_cells(cds_d101, color_cells_by = "pseudotime")

# 
theme_set(theme_minimal())
plot_ridge_pseudotime = function(cds, min_time=10, max_time=15, group.by = "genotype")
{
  df = as.data.frame(colData(cds))
  df$pseudotime = pseudotime(cds)
  
  pseudo_t = pseudotime(cds)
  selected_cells = names(pseudo_t)[(pseudo_t > min_time) & (pseudo_t < max_time)]
  selected_cells = setdiff(colnames(cds), selected_cells)
  cds@colData$pseudotime_v2 = pseudotime(cds)
  cds@colData[selected_cells, "pseudotime_v2"] = Inf
  p1 = plot_cells(cds, color_cells_by = "pseudotime_v2")+
    guides(fill = guide_legend( title = "" ))
  
  df = df[is.finite(df$pseudotime),]
  p2 = ggplot(df, aes(x = pseudotime, y = df[[group.by]], fill = pseudotime)) +
    stat_density_ridges(geom = "density_ridges_gradient",
                        quantile_lines = TRUE,
                        quantiles = 10) +
    theme(legend.position = c(0.1,0.9))+
    guides(fill = guide_legend( title = group.by ))
  
  # p2 = ggplot(df, aes(x = pseudotime, y = df[[group.by]])) +
  #     geom_density_ridges(aes(fill = df[[group.by]]), alpha = 0.8) + 
  #     theme(legend.position = c(0.1,0.9),
  #           panel.background = element_blank(),
  #           axis.line.x = 
  #           axis.ticks = element_blank(),
  #           axis.text = element_blank(),
  #           axis.title = element_blank())+
  #   guides(fill = guide_legend( title = group.by )) 
  
  plot_grid(p1, p2, nrow = 1, ncol = 2, align = 'h', rel_widths = c(1/2, 1/2))
}

plot_ridge_pseudotime(cds_d101)




