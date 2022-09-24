library(ggplot2)
library(Seurat)

setwd("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/")

# load data
df_cells = read.table("ctype_v4.txt", sep = "\t", header = 1, row.names = 1)
df_cells$source = plyr::mapvalues(df_cells$sample_group, from = c("Engineered-Cre", "Engineered-Flp", "Donor-Ctrl", "Donor-NRXN1del"),
                                  to = c("Engineered", "Engineered", "Patient-derived", "Patient-derived"))
df_cells$timepoint_v2 = plyr::mapvalues(df_cells$timepoint, from = c("d22", "d23", "d50", "d101", "d112"),
                                        to = c("d22/23","d22/23", "d50", "d101/112","d101/112"))
organoid = readRDS("combined_all_organoids_filtered_v2.rds")
df_cells = df_cells[colnames(organoid),]

organoid@meta.data = df_cells
organoid$cellclass_draft_v2 = factor(organoid$cellclass_draft_v2,
                                     levels = c("Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                                                "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                                                "Mesenchymal cell", "Microglia",
                                                "OPC", "IPC1","IPC2","IPC3",
                                                "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                                                "CN1","CN2","CN3","CN4","CN5",
                                                "low sequencing depth", "Unknown"))


# UMAPs (Figure 1)
# umap for cell type / sample group / time point
colour_cellclass = c("#00008B","#0000CD","#4169E1","#4682B4",
                     "#6495ED","#6A5ACD","#9400D3","#C71585","#C71585",
                     "#9ACD32","#6B8E23",
                     "#8B4513", "#D2B48C","#FFA500","#DAA520",
                     "#FF1493","#FF1493","#FF1493","#FF1493","#FF1493","#FF1493","#FF1493",
                     "#DC143C","#DC143C","#DC143C","#DC143C","#DC143C",
                     "#C0C0C0","#808080")
DimPlot(organoid, group.by = "cellclass_draft_v2",  label.size = 10, cols = colour_cellclass) + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle("")

# UMAP FOR timepoint
colour_timepoint = c()
DimPlot(organoid, group.by = "timepoint_v2",  label.size = 10, cols = colour_timepoint) + 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle("")

# cell abundance (INs)
INs = c("IN1","IN2","IN3","IN4","IN5","IN6","IN7")
time_lvls = c("d22/23", "d50", "d101/112")
sample_lvls = c("Donor-Ctrl", "Donor-NRXN1del", "Engineered-Flp", "Engineered-Cre")
df_INs = df_cells[df_cells$cellclass_draft_v2 %in% INs, ]

for (IN in INs)
{
  df_IN = df_INs[df_INs$cellclass_draft_v2 == IN, ]
  
}
test = df_INs[df_INs$cellclass_draft_v2 == "IN1",]
df_abun = data.frame(table(test$sample_group, test$timepoint_v2))
colnames(df_abun) = c("sample_group", "timepoint", "n_cells")
df_abun$timepoint = factor(df_abun$timepoint, levels = rev(time_lvls))
df_abun$sample_group = factor(df_abun$sample_group, levels = sample_lvls)
ggplot(df_abun, aes(x = sample_group, y = n_cells, fill = timepoint)) + 
  geom_bar(stat = "identity", position = "fill", binwi) + 
  theme(axis.text.y = element_text(family = "Arial", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank()) + 
  scale_fill_manual("legend", values = c("d22/23"="#FFE5CC", "d50"="#FFB266", "d101/112"="#990000"))
  







