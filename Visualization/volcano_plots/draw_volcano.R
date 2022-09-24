library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(plyr)
library(EnhancedVolcano)

source('../functions/_global.R')
source('../functions/_visualize.R')

# try enhancedvolcano on a specific plot
path = 'single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/DEGs/time_samplegroup_celltype_KO/'

deg_files = list.files(path, pattern = '^d')
names = gsub('.txt', '', deg_files)

draw_violinplot = function(deg_file, label, pthred = 1e-50, xmax = 5, percentile = 1/3000)
{
  # load
  df_deg = read.table(paste0(path, deg_file), sep = '\t', header = 1, row.names = 1)
  
  # remove unwanted genes
  RPS_genes = rownames(df_deg)[grepl(pattern = '^RPS', x = rownames(df_deg))]
  RPL_genes = rownames(df_deg)[grepl(pattern = '^RPL', x = rownames(df_deg))]
  MT_genes = rownames(df_deg)[grepl(pattern = '^MT-', x = rownames(df_deg))]
  remove_genes = c(RPS_genes, RPL_genes, MT_genes)
  
  df_deg = df_deg[setdiff(rownames(df_deg), remove_genes), ]
  
  # clip
  pthred = max(quantile(df_deg$pvals_adj, percentile)[[1]], pthred)
  xmax = min(xmax, quantile(abs(df_deg$logfoldchanges), 1-percentile)[[1]])
    
  df_deg[df_deg$pvals_adj < pthred, 'pvals_adj'] = pthred
  df_deg[df_deg$logfoldchanges < -xmax, 'logfoldchanges'] = -xmax
  df_deg[df_deg$logfoldchanges > xmax, 'logfoldchanges'] = xmax
  
  print(dim(df_deg))
  
  p = EnhancedVolcano(df_deg, 
                      lab = rownames(df_deg), 
                      x = 'logfoldchanges',
                      y = 'pvals_adj',
                      title = label,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      xlim = c(-xmax-0.1, xmax+0.1),
                      ylim = c(0, -log10(pthred) * 1.1))
  return(p)
}

for (i in 1:length(deg_files))
{
  deg_file = deg_files[[i]]
  label = names[[i]]
  p = draw_violinplot(deg_file, label, pthred = 1e-50, xmax = 5)
  ggsave(paste0('volcano_figures_v2/volcano_', label, '.eps'), device = cairo_ps, p,
         height = 20, width = 20, units = 'cm')
}




