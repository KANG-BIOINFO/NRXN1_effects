import pandas as pd
import numpy as np
import os

file_names = [i for i in os.listdir('results/') if i.endswith('_cellclass_2.txt')]
timepoints = [i.split('_')[1] for i in file_names]
fdr_vals = [i.split('_')[3] for i in file_names]
sources = ['Engineered' if 'Engineered' in i else 'Donor' for i in file_names]

df_combined = pd.DataFrame()
for timepoint, fdr_val, source, file_name in zip(timepoints, fdr_vals, sources, file_names):
    df = pd.read_csv('results/'+file_name, sep = '\t', header = 0)
    df['timepoint'] = timepoint
    df['source'] = source 
    df['fdr_used'] = fdr_val
    cols = np.setdiff1d(df.columns, ['Covariate'])
    df = df.loc[:, cols].copy()

    df_combined = pd.concat([df_combined, df], axis = 0)

df_combined.to_excel('scCODA_DA_analysis_combined.xlsx')
