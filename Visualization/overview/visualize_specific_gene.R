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

DefaultAssay(organoids_eng) = "RNA"
organoids_eng = organoids[, organoids$source == "Engineered"]
organoids_eng$combined = paste0(organoids_eng$cellclass_draft_v2, " - ", organoids_eng$genotype)
orders = unique(organoids_eng$combined)[order(unique(organoids_eng$combined))]
organoids_eng$combined = factor(organoids_eng$combined, levels = orders)
Idents(organoids_eng) = organoids_eng$combined

organoids_eng_d23 = organoids_eng[, organoids_eng$timepoint_v2 == "d023"]
organoids_eng_d50 = organoids_eng[, organoids_eng$timepoint_v2 == "d050"]
organoids_eng_d112 = organoids_eng[, organoids_eng$timepoint_v2 == "d112"]
VlnPlot(organoids_eng, features = c("NRXN1"), pt.size = 0.2) + theme(legend.position = "none")
VlnPlot(organoids_eng_d23, features = c("NRXN1"), pt.size = 0.2) + theme(legend.position = "none")
VlnPlot(organoids_eng_d50, features = c("NRXN1"), pt.size = 0.2) + theme(legend.position = "none")
VlnPlot(organoids_eng_d112, features = c("NRXN1"), pt.size = 0.2) + theme(legend.position = "none")




