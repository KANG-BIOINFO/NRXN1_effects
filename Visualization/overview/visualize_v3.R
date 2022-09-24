library(Seurat)
library(ggplot2)
library(miloR)
library(readxl)
library(gridExtra)
library(cowplot)

setwd("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/")

organoids = readRDS("./combined_all_organoids_filtered_v2.rds")
df_cells = read.table("ctype_v6.txt", sep = "\t", header = 1, row.names = 1)
df_cells = df_cells[colnames(organoids), ]
for (col in colnames(df_cells))
{
  if (! col %in% colnames(organoids@meta.data)){organoids@meta.data[[col]] = df_cells[[col]]}
}
df_order = read_excel("annotation/celltype_cellgroup_match.xlsx", sheet = 1)

organoids_donor = organoids[, organoids$source == "Donor"]
organoids_d22 = organoids_donor[, organoids_donor$timepoint_v2 == "d022"]
organoids_d50 = organoids_donor[, organoids_donor$timepoint_v2 == "d050"]
organoids_d101 = organoids_donor[, organoids_donor$timepoint_v2 == "d101"]

DimPlot(organoids_d22, group.by = "genotype")

# build miloR
# milo_d22 = Milo(organoids_d22)


# barplot
plot_barplot = function(df_cells, x = "cellclass_draft_v2", y = "genotype", y_scale = T)
{
  df_freq = data.frame(table(df_cells[[x]], df_cells[[y]]))
  colnames(df_freq) = c(x, y, "n_cells")
  df_freq[["log2_n_cells"]] = log2(df_freq[["n_cells"]]+1)
  df_freq[[x]] = factor(df_freq[[x]], levels = df_order[["cellclass"]])
  df_freq = df_freq[!df_freq[[x]] %in% c("low sequencing depth", "Unknown"),]
  
  count_name = "n_cells"
  if (y_scale){count_name = "log2_n_cells"}
    
  p1 = ggplot(df_freq, aes(x = df_freq[[x]], y = df_freq[[count_name]], fill = df_freq[[y]])) + 
        geom_bar(stat = "identity", position = "fill") + theme_bw()  + ylab("proportion") +
        theme(axis.text.y = element_text(family = "Arial", size = 12),
              axis.text.x = element_text(family = "Arial", angle = 90, vjust = 0.5, hjust = 1),
              axis.title = element_blank(), 
              legend.position = "none",
              plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) + 
        guides(fill = guide_legend( title = y )) + 
        scale_fill_manual("legend", values = c("Control" = "#99CCFF", "NRXN_del" = "#FF6666")) +
        geom_hline(yintercept = 0.5, linetype = "dashed",
                   color = "black", size = 1.2) 
  
  p2 = ggplot(df_freq, aes(x = df_freq[[x]], y = df_freq[[count_name]]), color = "#C0C0C0") +
        geom_bar(stat = "identity") + theme_bw() + 
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_text(family = "Arial", size = 12),
              axis.title.y = element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  
  plot_grid(p2, p1, nrow = 2, align = 'v', rel_heights = c(1/4, 3/4))
}

plot_barplot(organoids_d50, y_scale = F)
DimPlot(organoids_d101, group.by = "genotype", cols = c("#99CCFF","#FF6666")) 




