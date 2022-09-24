library(ggplot2)
library(Seurat)

# load data
setwd("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/visualize/")
organoids = readRDS("../combined_all_organoids_filtered_v2.rds")
df_cells = read.table("../ctype_v6.txt", sep = "\t", header = 1, row.names = 1)
df_cells = df_cells[colnames(organoids),]
organoids@meta.data = df_cells

organoids$cellclass_draft_v2 = factor(organoids$cellclass_draft_v2,
                                     levels = c("Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                                                "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                                                "Mesenchymal cell", "Microglia",
                                                "OPC", "IPC1","IPC2","IPC3",
                                                "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                                                "CN1","CN2","CN3","CN4","CN5",
                                                "low sequencing depth", "Unknown"))

organoids_eng = organoids[, organoids$source == "Engineered"]
organoids_donor = organoids[, organoids$source == "Donor"]

# vis
DimPlot(organoids, group.by = "cellclass_draft_v2", raster = F) + theme(legend.position = "none")

DimPlot(organoids_donor, group.by = "cellclass_draft_v2", raster = F) + theme(legend.position = "none")
DimPlot(organoids_eng, group.by = "cellclass_draft_v2", raster = F) + theme(legend.position = "none")

# vis each organoids 
colour_cellclass = c("#00008B","#0000CD","#4169E1","#4682B4",
                     "#6495ED","#6A5ACD","#9400D3","#C71585","#C71585",
                     "#9ACD32","#6B8E23",
                     "#8B4513", "#D2B48C","#FFA500","#DAA520",
                     "#FF1493","#FF1493","#FF1493","#FF1493","#FF1493","#FF1493","#FF1493",
                     "#DC143C","#DC143C","#DC143C","#DC143C","#DC143C",
                     "#C0C0C0","#808080")
DimPlot(organoids_donor, group.by = "cellclass_draft_v2",  label.size = 10, cols = colour_cellclass) + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle("")

DimPlot(organoids_eng, group.by = "cellclass_draft_v2",  label.size = 10, cols = colour_cellclass) + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle("")

# proportions between engineered and donor
df_eng = organoids_eng@meta.data
df_donor = organoids_donor@meta.data
df_match = read.table("../annotation/celltype_cellgroup_match.txt", sep = "\t", header = 1, row.names = 1)

df_table1 = data.frame(table(df_cells$cellclass_draft_v2, df_cells$source))
colnames(df_table1) = c("cell_class", "Source", "n_cells")
df_table1$cell_group = df_match[as.character(df_table1$cell_class), "cell.group"]
df_table1$order = df_match[as.character(df_table1$cell_class), "order"]
df_table1 = df_table1[order(df_table1$order),]
df_table1$cell_class = factor(df_table1$cell_class, levels = unique(as.character(df_table1$cell_class)))

ggplot(df_table1, aes(x = cell_class, y = n_cells, fill = Source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme(axis.text.y = element_text(family = "Arial", size = 12),
        axis.text.x = element_text(family = "Arial", angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9)) +
  xlab("cell_type") + ylab("n cells")




