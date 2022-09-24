library(Seurat)
library(ggplot2)

# combined = readRDS('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_filtered.rds')
# DimPlot(combined, group.by = 'integrated_snn_res.2', raster = F, label = T) + theme(legend.position = 'none')

# QC of donors
donors = readRDS('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient-derived_integration/patient_ipcs_all_organoids.rds')
donors = FindClusters(donors, resolution = 3.0)
DimPlot(donors, group.by = "integrated_snn_res.3", label = T)
p = VlnPlot(donors, features = c("nCount_RNA", "nFeature_RNA"),
        group.by = "integrated_snn_res.3", ncol = 1, pt.size = 0)
ggsave('quality_controls_donorCombined_res3.pdf', p, height = 20, width = 30, units = 'cm')

# QC of engineered
engs = readRDS('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/organoid_integrated_C2_RPCA.RDS')
p2 = VlnPlot(engs, features = c("nCount_RNA", "nFeature_RNA"),
            group.by = "integrated_snn_res.2", ncol = 1, pt.size = 0)
p2
ggsave('quality_controls_EngineeredCombined_res2.pdf', p2, height = 20, width = 27, units = 'cm')
