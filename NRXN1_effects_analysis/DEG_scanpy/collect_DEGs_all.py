import numpy as np
import pandas as pd 
import os

deg_files = [i for i in os.listdir('time_samplegroup_celltype_KO') if (i.startswith('d') and i.endswith('.txt'))]

df_combined = pd.DataFrame()
for deg_file in deg_files:
    timepoint = deg_file.split('_')[0]
    cell_type = deg_file.split('_')[1]
    source = deg_file.split('_')[2]

    df = pd.read_csv('time_samplegroup_celltype_KO/' + deg_file, sep = '\t', header = 0, index_col = 0)
    df['timepoint'] = timepoint
    df['source'] = source
    df['cell_type'] = cell_type
    df = df.loc[df['pvals_adj'] < 0.05, :].copy()
    df_combined = pd.concat([df_combined, df], axis = 0)

df_combined.to_excel('DEG_source_time_celltype_AllDEGs.xlsx')
    