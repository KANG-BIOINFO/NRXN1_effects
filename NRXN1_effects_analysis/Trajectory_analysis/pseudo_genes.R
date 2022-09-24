library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(ggridges)
library(dplyr)

cds = readRDS("./cds_monocle3.rds")
df_cells = read.table("../ctype_v6.txt", sep = "\t", header = 1, row.names = 1)[colnames(cds),]
for (col in colnames(df_cells)){
  if (!col %in% colnames(cds@colData)){cds@colData[[col]] = df_cells[[col]]}
}

plot_cells(cds, color_cells_by = "pseudotime", 
           label_roots = F, label_leaves = F, label_branch_points = F) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())

# preparation
combined_organoid = readRDS("../combined_all_organoids_filtered.rds")
cds = estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = rownames(combined_organoid[["RNA"]])
saveRDS(cds, "cds_monocle3_v2.rds")

# pseudotime-associated GENES
cds_pseudoGenes = graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
saveRDS(cds_pseudoGenes, "cds_pseudoGenes.rds")

# pseudotime distribution on different groups
cds_eng = cds[, colData(cds)[["source"]] == "Engineered"]
cds_donor = cds[, colData(cds)[["source"]] == "Donor"]

df = as.data.frame(colData(cds))
df$pseudotime = pseudotime(cds)
df = df[is.finite(df$pseudotime),]

ggplot(df, aes(x = pseudotime, y = source)) +
  geom_density_ridges(aes(fill = source)) + 
  theme(legend.position = c(0.1,0.9),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

# pseudotime distribution of multiple levels
ggplot(df[df$source == "Donor",], aes(x = pseudotime, y = timepoint_v2)) +
  geom_density_ridges(aes(fill = paste(timepoint_v2, genotype)), alpha = 0.8) + 
  theme(legend.position = c(0.1,0.9),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

ggplot(df[df$source == "Engineered",], aes(x = pseudotime, y = timepoint_v2)) +
  geom_density_ridges(aes(fill = paste(timepoint_v2, genotype)), alpha = 0.8) + 
  theme(legend.position = c(0.1,0.9),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
  

# reformat raw counts
organoids = readRDS("../combined_all_organoids_filtered_v2.rds")
genes_2k = rownames(organoids)
raw_matrix = organoids@assays$RNA@counts[genes_2k,]
cds@assays@data$counts = raw_matrix

# pseudotime-associated GENES
subset_cells = sample(colnames(cds), 1000)
cds_subset = cds[, subset_cells]
cds_subset@colData$pseudotime = pseudotime(cds_subset)
cds_subset = estimate_size_factors(cds_subset)
cds_subset = cds_subset[, is.finite(cds_subset@colData[["pseudotime"]])]
gene_fits = fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs = coefficient_table(gene_fits)
fit_coefs_pseudo = fit_coefs %>% filter(term == "pseudotime")
fit_coefs_pseudo$gene = genes_2k
fit_coefs_pseudo = fit_coefs_pseudo %>% filter (q_value < 0.05) %>% select(gene, term, q_value, estimate)
fit_coefs_pseudo = fit_coefs_pseudo[order(fit_coefs_pseudo$q_value),]

plot_cells(cds, genes = as.list(fit_coefs_pseudo$gene[1:3]))
FeaturePlot(organoids, features = fit_coefs_pseudo$gene[1:3], ncol = 3)
cds_subset_genes = cds_subset[fit_coefs_pseudo$gene[1:3],]
rowData(cds_subset_genes)[["gene_short_name"]] = rownames(rowData(cds_subset_genes))
plot_genes_in_pseudotime(cds_subset_genes, color_cells_by = "timepoint_v2", ncol = 3)